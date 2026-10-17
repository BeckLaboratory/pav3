"""Entrypoints that drive the lift and subseq operations over data sources.

Each entrypoint takes a prepared source table and an open-able output handler, iterates over data sources (sample/hap
pairs in PAV mode, or a single alignment file in Align mode), lifts coordinates, optionally extracts sequences, and
feeds rows to the output handler.
"""

__all__ = [
    'lift_pav',
    'lift_align',
    'subseq_pav',
    'subseq_align',
]

import logging
from pathlib import Path
from typing import Any, Optional

import polars as pl

from .. import pipeline
from .errors import LiftConfigError, LiftInputError
from .output import LiftOutput, SubseqOutput
from .source import LiftSource, LiftSourceAlign, LiftSourcePav

logger = logging.getLogger(__name__)


def _load_pav_context(pav_dir: Optional[str | Path]) -> tuple[Path, dict, pl.DataFrame]:
    """Load the PAV configuration and assembly table for a run directory.

    :param pav_dir: PAV run directory (defaults to the current directory).

    :returns: A tuple of (pav_dir, config, assembly table).
    """
    pav_dir = Path(pav_dir) if pav_dir is not None else Path('.')

    try:
        config = pipeline.read_pav_config(pipeline.get_pav_config_path(root_path=pav_dir))

        asm_table_path = pipeline.get_pav_asm_table_path(pav_config=config, root_path=pav_dir)
        asm_table = pipeline.read_assembly_table(asm_table_path, config)
    except (FileNotFoundError, ValueError) as e:
        raise LiftConfigError(f'Failed to load PAV configuration from "{pav_dir}": {e}') from e

    return pav_dir, config, asm_table


def _validate_regions(source_table: pl.DataFrame, fai_dict: Optional[dict[str, int]], reverse: bool) -> None:
    """Check that all source regions reference known sequences and fall within bounds.

    :param source_table: Source coordinate table.
    :param fai_dict: Sequence lengths for the source coordinate space.
    :param reverse: If `True`, source coordinates are query; otherwise reference.

    :raises ValueError: If sequence lengths are unavailable or any region is invalid.
    """
    if fai_dict is None:
        raise ValueError('Source sequence lengths are unavailable for region validation')

    space = 'query' if reverse else 'reference'
    problems = []

    for row in source_table.iter_rows(named=True):
        chrom = row['chrom']

        if chrom not in fai_dict:
            problems.append(f'{chrom}:{row["pos"] + 1}-{row["end"]} (unknown {space} sequence)')
            continue

        length = fai_dict[chrom]

        if row['pos'] < 0 or row['end'] < row['pos'] or row['end'] > length:
            problems.append(f'{chrom}:{row["pos"] + 1}-{row["end"]} (out of bounds, length={length})')

    if problems:
        shown = '; '.join(problems[:10])
        more = '' if len(problems) <= 10 else f' (and {len(problems) - 10} more)'
        raise ValueError(f'Invalid {space} regions: {shown}{more}')


def _unlifted_row(row: dict[str, Any], src: LiftSource, reverse: bool, do_seq: bool) -> dict[str, Any]:
    """Build an output row for a region that did not lift (source fields filled, lift fields null).

    :param row: Source-table row.
    :param src: Lift source (for sample/hap).
    :param reverse: If `True`, source coordinates are query; otherwise reference.
    :param do_seq: Include a null `seq` field if `True`.

    :returns: Partial output-row dict.
    """
    source_struct = {'chrom': row['chrom'], 'pos': row['pos'], 'end': row['end']}

    out_row = {
        'sample': src.sample,
        'hap': src.hap,
        'item': row['item'],
        'ref': source_struct if not reverse else None,
        'qry': source_struct if reverse else None,
        'strand': None,
        'align_index': None,
    }

    if do_seq:
        out_row['seq'] = None

    return out_row


def _process(
        source_table: pl.DataFrame,
        src: LiftSource,
        output,
        reverse: bool,
        allow_split: bool,
        do_seq: bool,
) -> None:
    """Lift every source region through one data source and write results.

    :param source_table: Source coordinate table.
    :param src: Lift source.
    :param output: Output handler (already opened).
    :param reverse: Lift query -> reference if `True`, else reference -> query.
    :param allow_split: Allow start/end to lift to different alignment records.
    :param do_seq: Extract sequences (subseq) if `True`.
    """
    _validate_regions(source_table, src.source_fai(reverse), reverse)

    for row in source_table.iter_rows(named=True):
        lifted = src.lift(row, reverse=reverse, allow_split=allow_split)

        if not lifted:
            logger.debug(
                'No lift: item=%s %s:%s-%s (sample=%s hap=%s)',
                row['item'], row['chrom'], row['pos'] + 1, row['end'], src.sample, src.hap,
            )

            if output.keep_unlifted:
                output.add(_unlifted_row(row, src, reverse, do_seq))

            continue

        for lift_row in lifted:
            if do_seq:
                lift_row['seq'] = src.seq(lift_row, reverse=reverse)

            output.add(lift_row)


def _pav_sources(pav_dir, sample, hap) -> tuple[Path, list[tuple[str, str]]]:
    """Resolve the (sample, hap) list for a PAV run, validating explicit selections.

    :param pav_dir: PAV run directory.
    :param sample: Sample name or `None` to expand all samples.
    :param hap: Haplotype name or `None` to expand all haplotypes.

    :returns: A tuple of (pav_dir, list of (sample, hap)).

    :raises ValueError: If an explicit sample/hap pair is not present.
    """
    pav_dir, _config, asm_table = _load_pav_context(pav_dir)

    try:
        source_list = list(pipeline.sample_hap_iter(asm_table, sample, hap))
    except ValueError as e:
        raise LiftConfigError(str(e)) from e

    if sample is not None and hap is not None and (sample, hap) not in source_list:
        raise LiftConfigError(f'Sample/haplotype not found or haplotype not defined: {sample}/{hap}')

    if not source_list:
        logger.warning('No samples/haplotypes to process')

    return pav_dir, source_list


def lift_pav(
        source_table: pl.DataFrame,
        output: LiftOutput,
        *,
        pav_dir: Optional[str | Path] = None,
        sample: Optional[str] = None,
        hap: Optional[str] = None,
        trim: str = 'qry',
        reverse: bool = False,
        allow_split: bool = False,
        ref_file: Optional[str | Path] = None,
) -> None:
    """Lift coordinates in PAV mode.

    :param source_table: Source coordinate table.
    :param output: Lift output handler.
    :param pav_dir: PAV run directory.
    :param sample: Sample name or `None` to expand all samples.
    :param hap: Haplotype name or `None` to expand all haplotypes.
    :param trim: Alignment trim level ("none", "qry", or "qryref").
    :param reverse: Lift query -> reference if `True`.
    :param allow_split: Allow start/end to lift to different alignment records.
    :param ref_file: Ignored in PAV mode (warning emitted if set).
    """
    pav_dir, source_list = _pav_sources(pav_dir, sample, hap)

    if source_table.height == 0:
        logger.warning('No regions to lift')

    with output:
        for src_sample, src_hap in source_list:
            logger.debug('Processing %s/%s', src_sample, src_hap)

            with LiftSourcePav(
                pav_dir=pav_dir, sample=src_sample, hap=src_hap, trim=trim, ref_file=ref_file,
            ) as src:
                _process(source_table, src, output, reverse, allow_split, do_seq=False)


def subseq_pav(
        source_table: pl.DataFrame,
        output: SubseqOutput,
        *,
        pav_dir: Optional[str | Path] = None,
        sample: Optional[str] = None,
        hap: Optional[str] = None,
        trim: str = 'qry',
        reverse: bool = False,
        allow_split: bool = False,
        uppercase: bool = False,
        ref_file: Optional[str | Path] = None,
) -> None:
    """Extract sequences in PAV mode.

    :param source_table: Source coordinate table.
    :param output: Subseq output handler.
    :param pav_dir: PAV run directory.
    :param sample: Sample name or `None` to expand all samples.
    :param hap: Haplotype name or `None` to expand all haplotypes.
    :param trim: Alignment trim level ("none", "qry", or "qryref").
    :param reverse: Lift query -> reference (extract reference sequence) if `True`.
    :param allow_split: Allow start/end to lift to different alignment records.
    :param uppercase: Uppercase extracted sequences if `True`.
    :param ref_file: Ignored in PAV mode (warning emitted if set).
    """
    pav_dir, source_list = _pav_sources(pav_dir, sample, hap)

    if source_table.height == 0:
        logger.warning('No regions to lift')

    with output:
        for src_sample, src_hap in source_list:
            logger.debug('Processing %s/%s', src_sample, src_hap)

            with LiftSourcePav(
                pav_dir=pav_dir, sample=src_sample, hap=src_hap,
                uppercase=uppercase, trim=trim, ref_file=ref_file,
            ) as src:
                _process(source_table, src, output, reverse, allow_split, do_seq=True)


def lift_align(
        source_table: pl.DataFrame,
        output: LiftOutput,
        *,
        align_file: str | Path,
        query_file: Optional[str | Path] = None,
        ref_file: Optional[str | Path] = None,
        sample: Optional[str] = None,
        hap: Optional[str] = None,
        reverse: bool = False,
        allow_split: bool = False,
) -> None:
    """Lift coordinates in Align mode.

    :param source_table: Source coordinate table.
    :param output: Lift output handler.
    :param align_file: Alignment file (SAM/BAM/CRAM).
    :param query_file: Indexed query FASTA/FASTQ (optional).
    :param ref_file: Reference FASTA (required for CRAM).
    :param sample: Sample name written verbatim to output.
    :param hap: Haplotype name written verbatim to output.
    :param reverse: Lift query -> reference if `True`.
    :param allow_split: Allow start/end to lift to different alignment records.
    """
    if source_table.height == 0:
        logger.warning('No regions to lift')

    with output:
        with LiftSourceAlign(
            align_file=align_file, query_file=query_file, ref_file=ref_file,
            sample=sample, hap=hap, no_seq=True,
        ) as src:
            _process(source_table, src, output, reverse, allow_split, do_seq=False)


def subseq_align(
        source_table: pl.DataFrame,
        output: SubseqOutput,
        *,
        align_file: str | Path,
        query_file: Optional[str | Path] = None,
        ref_file: Optional[str | Path] = None,
        sample: Optional[str] = None,
        hap: Optional[str] = None,
        reverse: bool = False,
        allow_split: bool = False,
        uppercase: bool = False,
) -> None:
    """Extract sequences in Align mode.

    :param source_table: Source coordinate table.
    :param output: Subseq output handler.
    :param align_file: Alignment file (SAM/BAM/CRAM).
    :param query_file: Indexed query FASTA/FASTQ (optional).
    :param ref_file: Reference FASTA (required for CRAM and for reverse subseq).
    :param sample: Sample name written verbatim to output.
    :param hap: Haplotype name written verbatim to output.
    :param reverse: Lift query -> reference (extract reference sequence) if `True`.
    :param allow_split: Allow start/end to lift to different alignment records.
    :param uppercase: Uppercase extracted sequences if `True`.

    :raises LiftInputError: If reverse extraction is requested without a reference file.
    """
    if reverse and ref_file is None:
        raise LiftInputError('A reference file (-f/--reference) is required for reverse subseq in Align mode')

    if source_table.height == 0:
        logger.warning('No regions to lift')

    with output:
        with LiftSourceAlign(
            align_file=align_file, query_file=query_file, ref_file=ref_file,
            sample=sample, hap=hap, uppercase=uppercase, no_seq=False,
        ) as src:
            _process(source_table, src, output, reverse, allow_split, do_seq=True)
