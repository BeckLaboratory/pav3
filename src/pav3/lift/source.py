"""Coordinate-lifting sources for the lift and subseq tools.

`LiftSource` is an abstract interface over a single (sample, haplotype) alignment set. `LiftSourcePav` reads a PAV run
directory; `LiftSourceAlign` reads a single SAM/BAM/CRAM alignment file. Both lift coordinates (reference <-> query) and
extract sequences.
"""

__all__ = [
    'LiftSource',
    'LiftSourcePav',
    'LiftSourceAlign',
]

import logging
import os
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Optional

import Bio.Seq
import polars as pl
import pysam

import agglovar

from .. import pipeline
from ..align.lift import AlignLift, get_lift_pairs
from ..align.tables import sam_to_align_table
from .errors import LiftInputError

logger = logging.getLogger(__name__)


def _strand(is_rev: bool) -> str:
    """Get a strand character for an orientation.

    :param is_rev: `True` if reverse-complemented.

    :returns: "-" if `is_rev` else "+".
    """
    return '-' if is_rev else '+'


def _revcomp(seq: str) -> str:
    """Reverse-complement a sequence.

    :param seq: Sequence.

    :returns: Reverse-complemented sequence.
    """
    return str(Bio.Seq.Seq(seq).reverse_complement())


class LiftSource(ABC):
    """Abstract base class for lifting coordinates and extracting sequences for one (sample, haplotype)."""

    def __init__(
            self,
            sample: Optional[str] = None,
            hap: Optional[str] = None,
            uppercase: bool = False,
    ) -> None:
        """Initialize a lift source.

        :param sample: Sample name (written verbatim to output).
        :param hap: Haplotype name (written verbatim to output).
        :param uppercase: Uppercase extracted sequences if `True`.
        """
        self.sample = sample
        self.hap = hap
        self.uppercase = uppercase
        self._align_lift: Optional[AlignLift] = None

    @property
    @abstractmethod
    def ref_fai(self) -> Optional[dict[str, int]]:
        """Reference sequence lengths keyed by chromosome, or `None` if unavailable."""

    @property
    @abstractmethod
    def qry_fai(self) -> dict[str, int]:
        """Query sequence lengths keyed by query ID."""

    @abstractmethod
    def _extract_seq(self, kind: str, chrom: str, pos: int, end: int, is_rev: bool) -> str:
        """Extract a subsequence from a source.

        :param kind: "ref" for the reference source or "qry" for the query source.
        :param chrom: Sequence name.
        :param pos: Start position (0-based).
        :param end: End position (0-based, exclusive).
        :param is_rev: Reverse-complement the sequence if `True`.

        :returns: Extracted sequence.
        """

    def source_fai(self, reverse: bool) -> Optional[dict[str, int]]:
        """Get the sequence-length dict for the source (input) coordinate space.

        :param reverse: If `True`, source coordinates are query; otherwise reference.

        :returns: Sequence-length dict for the source coordinate space.
        """
        return self.qry_fai if reverse else self.ref_fai

    def lift(
            self,
            row: dict[str, Any],
            reverse: bool = False,
            allow_split: bool = False,
    ) -> list[dict[str, Any]]:
        """Lift a source region through the alignment.

        By default (`reverse=False`) `row` holds reference coordinates and the result gives query coordinates. When
        `reverse=True`, `row` holds query coordinates and the result gives reference coordinates. A region that lifts to
        multiple alignment records yields one output row per record, ordered best-first.

        :param row: Source-table row with `chrom`, `pos`, `end`, `item`.
        :param reverse: Lift query -> reference if `True`, else reference -> query.
        :param allow_split: Allow the start and end to lift to different alignment records (same query ID and strand).

        :returns: List of output-row dicts (empty if the region does not lift).
        """
        align_lift = self._align_lift
        chrom = row['chrom']
        pos = row['pos']
        end = row['end']
        end_coord = end - 1 if end > pos else pos

        if not reverse:
            lift_a = align_lift.to_qry(chrom, pos)
            lift_b = align_lift.to_qry(chrom, end_coord)
        else:
            lift_a = align_lift.to_ref(chrom, pos)
            lift_b = align_lift.to_ref(chrom, end_coord)

        pairs = get_lift_pairs(
            a=lift_a,
            b=lift_b,
            same_align=False,
            single_align=not allow_split,
            multi='best',
            len_=end - pos,
        ) or []

        out_rows = []

        for pair_a, pair_b in pairs:
            # get_lift_pairs matches source chrom/strand, not target ID; enforce a shared target here
            if pair_a['lift_id'] != pair_b['lift_id'] or pair_a['is_rev'] != pair_b['is_rev']:
                continue

            low, high = sorted((pair_a, pair_b), key=lambda lift: lift['lift_pos'])
            is_rev = bool(low['is_rev'])

            # A zero-length source region (pos == end) lifts to a zero-length target, not a 1-bp span.
            target_end = high['lift_pos'] + 1 if end > pos else low['lift_pos']

            target = {'chrom': low['lift_id'], 'pos': low['lift_pos'], 'end': target_end}
            source = {'chrom': chrom, 'pos': pos, 'end': end}

            ref_struct, qry_struct = (target, source) if reverse else (source, target)

            out_rows.append({
                'sample': self.sample,
                'hap': self.hap,
                'item': row['item'],
                'ref': ref_struct,
                'qry': qry_struct,
                'strand': _strand(is_rev),
                'align_index': {'pos': low['align_index'], 'end': high['align_index']},
            })

        return out_rows

    def seq(
            self,
            row: dict[str, Any],
            reverse: bool = False,
    ) -> str:
        """Extract the sequence for a lifted output row.

        The extracted sequence is the target coordinate (query by default, reference if `reverse`) in the orientation
        matching the anchor (reverse-complemented if the alignment is reverse).

        :param row: Output-row dict from :meth:`lift`.
        :param reverse: Matches the `reverse` used for the lift.

        :returns: Extracted sequence (empty string if the target is undefined).
        """
        kind = 'ref' if reverse else 'qry'
        struct = row.get(kind)

        if struct is None:
            return ''

        sequence = self._extract_seq(kind, struct['chrom'], struct['pos'], struct['end'], row.get('strand') == '-')

        return sequence.upper() if self.uppercase else sequence

    def close(self) -> None:  # noqa: B027
        """Release any open resources. Subclasses override as needed."""

    def __enter__(self) -> 'LiftSource':
        """Enter context guard."""
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        """Exit context guard, releasing resources."""
        self.close()


class LiftSourcePav(LiftSource):
    """Lift source backed by a PAV run directory."""

    def __init__(
            self,
            pav_dir: Optional[str | Path] = None,
            sample: Optional[str] = None,
            hap: Optional[str] = None,
            uppercase: bool = False,
            trim: str = 'qry',
            ref_file: Optional[str | Path] = None,
    ) -> None:
        """Initialize a PAV-mode lift source.

        :param pav_dir: PAV run directory (defaults to the current directory).
        :param sample: Sample name.
        :param hap: Haplotype name.
        :param uppercase: Uppercase extracted sequences if `True`.
        :param trim: Alignment trim level ("none", "qry", or "qryref").
        :param ref_file: Ignored in PAV mode (a warning is emitted if given).

        :raises FileNotFoundError: If required PAV files are missing.
        """
        super().__init__(sample=sample, hap=hap, uppercase=uppercase)

        if ref_file is not None:
            logger.warning('Reference file is ignored in PAV mode (using the PAV run reference): %s', ref_file)

        self.pav_dir = Path(pav_dir) if pav_dir is not None else Path('.')

        align_path = self.pav_dir / 'results' / sample / 'align' / hap / f'align_trim-{trim}.parquet'

        if not align_path.is_file():
            raise FileNotFoundError(f'Alignment table not found: {align_path}')

        df_align = pl.read_parquet(align_path)

        qry_files = pipeline.expand_fofn(
            self.pav_dir / 'data' / 'query' / sample / f'query_{hap}.fofn', rel=self.pav_dir,
        )
        self._qry_fa = str(qry_files[0])
        df_qry_fai = agglovar.fa.read_fai(str(qry_files[1]), name='qry_id')

        ref_files = pipeline.expand_fofn(self.pav_dir / 'data' / 'ref' / 'ref.fofn', rel=self.pav_dir)
        self._ref_fa = str(ref_files[0])
        df_ref_fai = agglovar.fa.read_fai(str(ref_files[1]))

        self._align_lift = AlignLift(df_align, df_qry_fai)
        self._ref_fai = dict(zip(df_ref_fai['chrom'], df_ref_fai['len']))
        self._qry_fai = dict(zip(df_qry_fai['qry_id'], df_qry_fai['len']))

        self._qry_handle: Optional[pysam.FastaFile] = None
        self._ref_handle: Optional[pysam.FastaFile] = None

    @property
    def ref_fai(self) -> dict[str, int]:
        """Reference sequence lengths keyed by chromosome."""
        return self._ref_fai

    @property
    def qry_fai(self) -> dict[str, int]:
        """Query sequence lengths keyed by query ID."""
        return self._qry_fai

    def _extract_seq(self, kind: str, chrom: str, pos: int, end: int, is_rev: bool) -> str:
        """Extract a subsequence from the PAV query or reference FASTA."""
        if kind == 'qry':
            if self._qry_handle is None:
                self._qry_handle = pysam.FastaFile(self._qry_fa)
            handle = self._qry_handle
        else:
            if self._ref_handle is None:
                self._ref_handle = pysam.FastaFile(self._ref_fa)
            handle = self._ref_handle

        sequence = handle.fetch(chrom, pos, end)

        return _revcomp(sequence) if is_rev else sequence

    def close(self) -> None:
        """Close open FASTA handles."""
        if self._qry_handle is not None:
            self._qry_handle.close()
            self._qry_handle = None

        if self._ref_handle is not None:
            self._ref_handle.close()
            self._ref_handle = None


class LiftSourceAlign(LiftSource):
    """Lift source backed by a single SAM/BAM/CRAM alignment file."""

    def __init__(
            self,
            align_file: str | Path,
            query_file: Optional[str | Path] = None,
            ref_file: Optional[str | Path] = None,
            sample: Optional[str] = None,
            hap: Optional[str] = None,
            uppercase: bool = False,
            no_seq: bool = False,
    ) -> None:
        """Initialize an Align-mode lift source.

        :param align_file: Alignment file (SAM/BAM/CRAM).
        :param query_file: Indexed query FASTA/FASTQ. If omitted, sequences and lengths are recovered from the
            alignment records.
        :param ref_file: Reference FASTA. Required to decode CRAM and to extract reference sequences (reverse subseq).
        :param sample: Sample name (written verbatim to output).
        :param hap: Haplotype name (written verbatim to output).
        :param uppercase: Uppercase extracted sequences if `True`.
        :param no_seq: Do not retain query sequences in memory (lift only).

        :raises FileNotFoundError: If the query FASTA index is missing.
        :raises LiftInputError: If CRAM is used without a reference.
        """
        super().__init__(sample=sample, hap=hap, uppercase=uppercase)

        self.align_file = str(align_file)
        self.query_file = str(query_file) if query_file is not None else None
        self.ref_file = str(ref_file) if ref_file is not None else None
        self.no_seq = no_seq

        if self.align_file.lower().endswith('.cram') and self.ref_file is None:
            raise LiftInputError('A reference file (-f/--reference) is required to read a CRAM alignment file')

        self._qry_seq: Optional[dict[str, str]] = None
        self._qry_handle: Optional[pysam.FastaFile] = None
        self._ref_handle: Optional[pysam.FastaFile] = None

        # Reference lengths from the alignment header (a header read, not a record scan)
        with pysam.AlignmentFile(self.align_file, reference_filename=self.ref_file) as align_handle:
            self._ref_fai = dict(zip(align_handle.header.references, align_handle.header.lengths))

        if self.query_file is not None:
            fai_path = self.query_file + '.fai'

            if not os.path.isfile(fai_path):
                raise FileNotFoundError(
                    f'Query FASTA/FASTQ index not found (expected "{fai_path}"); index with samtools faidx'
                )

            df_qry_fai = agglovar.fa.read_fai(fai_path, name='qry_id')

            df_align = sam_to_align_table(
                Path(self.align_file),
                df_qry_fai,
                ref_fa_filename=self.ref_file,
                qry_fa_filename=self.query_file,
            )
        else:
            # Single scan: build the alignment table and derive query lengths (and, for subseq, sequences) at once
            df_align, df_qry_fai, self._qry_seq = sam_to_align_table(
                Path(self.align_file),
                None,
                ref_fa_filename=self.ref_file,
                qry_fa_filename=self.query_file,
                collect_qry_fai=True,
                collect_qry_seq=not no_seq,
            )

        self._align_lift = AlignLift(df_align, df_qry_fai)
        self._qry_fai = dict(zip(df_qry_fai['qry_id'], df_qry_fai['len']))

    @property
    def ref_fai(self) -> dict[str, int]:
        """Reference sequence lengths keyed b3y chromosome (from the alignment header)."""
        return self._ref_fai

    @property
    def qry_fai(self) -> dict[str, int]:
        """Query sequence lengths keyed by query ID."""
        return self._qry_fai

    def _extract_seq(self, kind: str, chrom: str, pos: int, end: int, is_rev: bool) -> str:
        """Extract a subsequence from the query (in-memory or FASTA) or reference FASTA."""
        if kind == 'qry':
            if self._qry_seq is not None:
                if chrom not in self._qry_seq:
                    raise LiftInputError(
                        f'Query sequence for "{chrom}" is unavailable (SEQ missing or hard-clipped in the '
                        f'alignment); provide a query file with -q'
                    )

                sequence = self._qry_seq[chrom][pos:end]
            else:
                if self._qry_handle is None:
                    self._qry_handle = pysam.FastaFile(self.query_file)
                sequence = self._qry_handle.fetch(chrom, pos, end)
        else:
            if self.ref_file is None:
                raise LiftInputError('A reference file (-f/--reference) is required to extract reference sequences')

            if self._ref_handle is None:
                self._ref_handle = pysam.FastaFile(self.ref_file)

            sequence = self._ref_handle.fetch(chrom, pos, end)

        return _revcomp(sequence) if is_rev else sequence

    def close(self) -> None:
        """Close open FASTA handles."""
        if self._qry_handle is not None:
            self._qry_handle.close()
            self._qry_handle = None

        if self._ref_handle is not None:
            self._ref_handle.close()
            self._ref_handle = None
