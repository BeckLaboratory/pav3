"""Subseq subcommand.

Extract sequences corresponding to a set of regions lifted through assembly alignments. By default, reference regions
yield the corresponding query sequence in reference orientation; with -r/--reverse, query regions yield the
corresponding reference sequence. Runs in PAV mode or Align mode (-a/--align-file).
"""

import argparse
import logging
from typing import Optional

from ..lift import subseq_align, subseq_pav
from ..util import init_logger
from ._common_opt import (
    _add_opt_lift_common,
    _add_opt_reference,
    _add_opt_subseq_extra,
    _add_opt_version,
)
from ._lift_cli import build_source_table, resolve_output
from ._returncode import ExitCode, run_guarded

logger = logging.getLogger(__name__)


def subcommand_subseq(
        coordinates: Optional[list[str]] = None,
        sample: Optional[str] = None,
        hap: Optional[str] = None,
        input: Optional[str] = None,
        input_format: Optional[str] = None,
        output: Optional[str] = None,
        output_format: Optional[str] = None,
        trim: str = 'qry',
        reverse: bool = False,
        pav_directory: Optional[str] = None,
        allow_split: bool = False,
        align_file: Optional[str] = None,
        query_file: Optional[str] = None,
        query_format: Optional[str] = None,
        reference: Optional[str] = None,
        upper: bool = False,
        verbose: bool = False,
) -> int:
    """Run the subseq subcommand.

    :param coordinates: Positional region strings.
    :param sample: Sample name (PAV mode expansion if `None`).
    :param hap: Haplotype name (PAV mode expansion if `None`).
    :param input: Input coordinate file (standard input if `None` and no coordinates).
    :param input_format: Input format override.
    :param output: Output file (standard output if `None`).
    :param output_format: Output format override.
    :param trim: Alignment trim level (PAV mode).
    :param reverse: Lift query -> reference (extract reference sequence) if `True`.
    :param pav_directory: PAV run directory (PAV mode).
    :param allow_split: Allow start/end to lift to different alignment records.
    :param align_file: Alignment file (enables Align mode).
    :param query_file: Query FASTA/FASTQ (Align mode).
    :param query_format: Query file format override (informational).
    :param reference: Reference FASTA (Align mode reverse subseq).
    :param upper: Uppercase extracted sequences if `True`.
    :param verbose: Enable DEBUG logging.

    :returns: Exit code.
    """
    init_logger(logging.DEBUG if verbose else logging.INFO, force=True)

    if align_file is not None and pav_directory is not None:
        logger.error('Cannot combine -a/--align-file with -d/--pav-directory (Align mode uses no PAV directory)')
        return int(ExitCode.USAGE)

    if align_file is None and query_file is not None:
        logger.error('-q/--query-file requires -a/--align-file')
        return int(ExitCode.USAGE)

    def action() -> None:
        source_table = build_source_table(coordinates, input, input_format)
        out_handler = resolve_output(True, output, output_format, reverse)

        if align_file is not None:
            subseq_align(
                source_table, out_handler,
                align_file=align_file, query_file=query_file, ref_file=reference,
                sample=sample, hap=hap, reverse=reverse, allow_split=allow_split, uppercase=upper,
            )
        else:
            subseq_pav(
                source_table, out_handler,
                pav_dir=pav_directory, sample=sample, hap=hap, trim=trim,
                reverse=reverse, allow_split=allow_split, uppercase=upper, ref_file=reference,
            )

    return run_guarded(action)


def _add_subparser_subseq(
        subparsers,
) -> argparse.ArgumentParser:
    """Add the subseq subcommand to a subparser.

    :param subparsers: Subparser object.

    :returns: Configured subparser.
    """
    parser_subseq = subparsers.add_parser(
        'subseq',
        description='Extract sequences for regions lifted through assembly alignments.',
        help='Extract sequences for lifted regions.',
    )

    _add_opt_version(parser_subseq)
    _add_opt_lift_common(parser_subseq)
    _add_opt_reference(parser_subseq)
    _add_opt_subseq_extra(parser_subseq)

    return parser_subseq
