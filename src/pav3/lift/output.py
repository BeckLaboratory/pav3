"""Output handlers for the lift and subseq tools.

Each handler is a context manager. Text handlers (CSV/TSV/BED) stream rows in source order as they arrive. Parquet
handlers accumulate rows, sort by `item` (stable) at close, and (for subseq) offload large sequence buffers to
temporary files. All handlers always create an output (a header line, a valid empty Parquet, or an empty FASTA) even
when no rows are written.
"""

__all__ = [
    'LIFT_SCHEMA',
    'SUBSEQ_SCHEMA',
    'LiftOutput',
    'LiftOutputCsv',
    'LiftOutputBed',
    'LiftOutputParquet',
    'SubseqOutput',
    'SubseqOutputCsv',
    'SubseqOutputBed',
    'SubseqOutputFasta',
    'SubseqOutputParquet',
]

import gzip
import logging
import shutil
import sys
import tempfile
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Optional, TextIO

import polars as pl

from ..io import BGZFWriterIO

logger = logging.getLogger(__name__)

_REGION_STRUCT = pl.Struct({'chrom': pl.String, 'pos': pl.Int64, 'end': pl.Int64})
_ALIGN_INDEX_STRUCT = pl.Struct({'pos': pl.Int32, 'end': pl.Int32})

LIFT_SCHEMA: dict[str, pl.DataType] = {
    'sample': pl.String,
    'hap': pl.String,
    'item': pl.Int32,
    'ref': _REGION_STRUCT,
    'qry': _REGION_STRUCT,
    'strand': pl.String,
    'align_index': _ALIGN_INDEX_STRUCT,
}
"""Schema of the lift output table."""

SUBSEQ_SCHEMA: dict[str, pl.DataType] = {**LIFT_SCHEMA, 'seq': pl.String}
"""Schema of the subseq output table (lift schema plus `seq`)."""

_LIFT_FIELDS = ['sample', 'hap', 'item', 'ref', 'qry', 'strand', 'align_index']
_SUBSEQ_FIELDS = _LIFT_FIELDS + ['seq']

_DEFAULT_CHUNK_BP = 1_000_000


def _region_str(struct: Optional[dict], empty: str = '') -> str:
    """Render a coordinate struct as a 1-based, comma-grouped region string.

    :param struct: Region struct (`chrom`, `pos`, `end`) or `None`.
    :param empty: Placeholder returned when `struct` is `None`.

    :returns: Region string (e.g. "chr1:1,001-2,000") or `empty` if `struct` is `None`.
    """
    if struct is None:
        return empty

    return f'{struct["chrom"]}:{struct["pos"] + 1:,}-{struct["end"]:,}'


def _align_index_str(struct: Optional[dict], empty: str = '') -> str:
    """Render an align_index struct as "pos,end" (no thousands separators).

    :param struct: Align-index struct (`pos`, `end`) or `None`.
    :param empty: Placeholder returned when `struct` is `None`.

    :returns: Formatted string or `empty` if `struct` is `None`.
    """
    if struct is None:
        return empty

    return f'{struct["pos"]},{struct["end"]}'


def _cell(value: Any, empty: str = '') -> str:
    """Render a scalar cell value, mapping `None` to a placeholder.

    :param value: Cell value.
    :param empty: Placeholder returned when `value` is `None`.

    :returns: String representation, or `empty` for `None`.
    """
    return empty if value is None else str(value)


def _delimited_line(cells: list[str], separator: str) -> str:
    """Join cells into a delimited line, quoting any cell that contains the separator, a quote, or a newline.

    Region strings carry thousands-separator commas and align-index cells are rendered "pos,end", so CSV output
    (comma separator) must quote them (RFC 4180) to avoid splitting a value across columns. Tab-delimited output never
    triggers quoting because no rendered cell contains a tab.

    :param cells: Cell strings.
    :param separator: Field separator.

    :returns: The joined line (without a trailing newline).
    """
    quoted = []

    for cell in cells:
        if separator in cell or '"' in cell or '\n' in cell or '\r' in cell:
            quoted.append('"' + cell.replace('"', '""') + '"')
        else:
            quoted.append(cell)

    return separator.join(quoted)


def _open_text(out_file: Optional[str], is_gz: bool) -> tuple[TextIO, bool]:
    """Open a text output handle (file or stdout).

    :param out_file: Output file path, or `None` for stdout.
    :param is_gz: Open with gzip compression if `True`.

    :returns: A tuple of (handle, close_handle) where `close_handle` is `False` for stdout.
    """
    if out_file is None:
        return sys.stdout, False

    if is_gz:
        return gzip.open(out_file, 'wt', encoding='utf-8'), True

    return open(out_file, 'wt', encoding='utf-8'), True


def _target_struct(row: dict, reverse: bool) -> Optional[dict]:
    """Get the target coordinate struct for a row (query by default, reference if `reverse`).

    :param row: Output-row dict.
    :param reverse: If `True`, the target is the reference coordinate.

    :returns: Target coordinate struct.
    """
    return row.get('ref') if reverse else row.get('qry')


class _OutputBase(ABC):
    """Common base for output handlers."""

    keep_unlifted: bool = False

    def __init__(self, out_file: Optional[str] = None, reverse: bool = False) -> None:
        """Initialize an output handler.

        :param out_file: Output file path, or `None` for stdout.
        :param reverse: Matches the `reverse` used for the lift (selects the BED/FASTA target coordinate).
        """
        self.out_file = out_file
        self.reverse = reverse

    def __enter__(self) -> '_OutputBase':
        """Enter context guard and open the output."""
        self._open()
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        """Exit context guard and finalize the output."""
        self._close()

    @abstractmethod
    def _open(self) -> None:
        """Open the output and write any header."""

    @abstractmethod
    def _close(self) -> None:
        """Finalize and close the output."""

    @abstractmethod
    def add(self, row: dict[str, Any]) -> None:
        """Write or buffer one output row.

        :param row: Output-row dict.
        """


class _TextTableOutput(_OutputBase):
    """Base for streaming delimited text output (CSV/TSV/BED)."""

    keep_unlifted = True
    fields: list[str] = _LIFT_FIELDS

    def __init__(
            self,
            out_file: Optional[str] = None,
            reverse: bool = False,
            is_gz: bool = False,
            separator: str = ',',
    ) -> None:
        """Initialize a delimited text output handler.

        :param out_file: Output file path, or `None` for stdout.
        :param reverse: Matches the `reverse` used for the lift.
        :param is_gz: Gzip-compress the output if `True`.
        :param separator: Field separator.
        """
        super().__init__(out_file, reverse)
        self.is_gz = is_gz
        self.separator = separator
        self._handle: Optional[TextIO] = None
        self._close_handle = False

    def _header(self) -> list[str]:
        """Get header column names."""
        return self.fields

    def _open(self) -> None:
        """Open the output and write the header line."""
        self._handle, self._close_handle = _open_text(self.out_file, self.is_gz)
        self._handle.write(_delimited_line(self._header(), self.separator) + '\n')

    def _close(self) -> None:
        """Close the output handle (never closes stdout)."""
        if self._handle is not None and self._close_handle:
            self._handle.close()

        self._handle = None

    def _row_cells(self, row: dict[str, Any]) -> list[str]:
        """Render a row as a list of cell strings.

        :param row: Output-row dict.

        :returns: List of formatted cells matching the handler's fields.
        """
        cells = [
            _cell(row.get('sample')),
            _cell(row.get('hap')),
            _cell(row.get('item')),
            _region_str(row.get('ref')),
            _region_str(row.get('qry')),
            _cell(row.get('strand')),
            _align_index_str(row.get('align_index')),
        ]

        if 'seq' in self.fields:  # subseq output includes the sequence column
            cells.append(_cell(row.get('seq')))

        return cells

    def add(self, row: dict[str, Any]) -> None:
        """Write one row."""
        self._handle.write(_delimited_line(self._row_cells(row), self.separator) + '\n')


class _BedOutput(_TextTableOutput):
    """Base for streaming BED output (target coordinate first, then all fields)."""

    keep_unlifted = False

    def __init__(self, out_file: Optional[str] = None, reverse: bool = False, is_gz: bool = False) -> None:
        """Initialize a BED output handler.

        :param out_file: Output file path, or `None` for stdout.
        :param reverse: Matches the `reverse` used for the lift (selects the target coordinate).
        :param is_gz: Gzip-compress the output if `True`.
        """
        super().__init__(out_file, reverse, is_gz=is_gz, separator='\t')

    def _header(self) -> list[str]:
        """Get the BED header (commented), starting with target chrom/pos/end."""
        return ['#chrom', 'pos', 'end'] + self.fields

    def add(self, row: dict[str, Any]) -> None:
        """Write one BED row (target coordinate first; empty target columns if the region did not lift)."""
        target = _target_struct(row, self.reverse)

        if target is None:
            coord_cells = ['', '', '']
        else:
            coord_cells = [target['chrom'], str(target['pos']), str(target['end'])]

        self._handle.write(_delimited_line(coord_cells + self._row_cells(row), self.separator) + '\n')


class _ParquetOutput(_OutputBase):
    """Base for buffered Parquet output with optional temp-file offloading."""

    keep_unlifted = True
    schema: dict[str, pl.DataType] = LIFT_SCHEMA
    accumulate_seq = False

    def __init__(
            self,
            out_file: Optional[str] = None,
            reverse: bool = False,
            tmp_dir: Optional[str | bool] = None,
            chunk_bp: int = _DEFAULT_CHUNK_BP,
    ) -> None:
        """Initialize a Parquet output handler.

        :param out_file: Output file path (required; Parquet cannot be written to stdout).
        :param reverse: Matches the `reverse` used for the lift.
        :param tmp_dir: Temporary directory for chunk offloading: a path, `None`/`True` for the system temp directory,
            or `False` to disable (accumulate in memory).
        :param chunk_bp: Offload the buffer to a temp file once accumulated sequence length reaches this many bases.

        :raises ValueError: If `out_file` is `None`.
        """
        super().__init__(out_file, reverse)

        if out_file is None:
            raise ValueError('Parquet output requires an output file')

        self.tmp_dir = tmp_dir
        self.chunk_bp = chunk_bp
        self._rows: list[dict[str, Any]] = []
        self._buffer_bp = 0
        self._chunk_paths: list[Path] = []
        self._tmp_root: Optional[Path] = None

    def _open(self) -> None:
        """No header to write for Parquet."""

    def _normalize(self, row: dict[str, Any]) -> dict[str, Any]:
        """Project a row onto the output schema keys.

        :param row: Output-row dict.

        :returns: Row dict limited to schema keys (missing keys filled with `None`).
        """
        return {key: row.get(key) for key in self.schema}

    def add(self, row: dict[str, Any]) -> None:
        """Buffer one row, offloading to a temp file when the sequence buffer is full."""
        normalized = self._normalize(row)
        self._rows.append(normalized)

        if self.accumulate_seq and self.tmp_dir is not False:
            self._buffer_bp += len(normalized.get('seq') or '')

            if self._buffer_bp >= self.chunk_bp:
                self._flush_chunk()

    def _ensure_tmp(self) -> None:
        """Create the temporary directory on first use."""
        if self._tmp_root is not None:
            return

        base = None if self.tmp_dir in (None, True) else str(self.tmp_dir)
        self._tmp_root = Path(tempfile.mkdtemp(prefix='pav_lift_', dir=base))

    def _flush_chunk(self) -> None:
        """Write the current buffer to a temp Parquet file and clear it."""
        if not self._rows:
            return

        self._ensure_tmp()
        chunk_path = self._tmp_root / f'chunk_{len(self._chunk_paths):06d}.parquet'
        pl.DataFrame(self._rows, schema=self.schema).write_parquet(chunk_path)
        self._chunk_paths.append(chunk_path)
        self._rows = []
        self._buffer_bp = 0

    def _close(self) -> None:
        """Concatenate any chunks with the in-memory buffer, sort by `item`, and write the output."""
        try:
            if self._chunk_paths:
                frames = [pl.scan_parquet(path) for path in self._chunk_paths]

                if self._rows:
                    frames.append(pl.LazyFrame(self._rows, schema=self.schema))

                (
                    pl.concat(frames)
                    .sort('item', maintain_order=True)
                    .sink_parquet(self.out_file)
                )
            else:
                frame = (
                    pl.DataFrame(self._rows, schema=self.schema)
                    if self._rows else pl.DataFrame(schema=self.schema)
                )
                frame.sort('item', maintain_order=True).write_parquet(self.out_file)
        finally:
            if self._tmp_root is not None:
                shutil.rmtree(self._tmp_root, ignore_errors=True)
                self._tmp_root = None


#
# Lift output handlers
#

class LiftOutput(_OutputBase):
    """Abstract base for lift output handlers."""


class LiftOutputCsv(_TextTableOutput, LiftOutput):
    """CSV/TSV lift output (separator selects CSV vs TSV)."""

    fields = _LIFT_FIELDS


class LiftOutputBed(_BedOutput, LiftOutput):
    """BED lift output."""

    fields = _LIFT_FIELDS


class LiftOutputParquet(_ParquetOutput, LiftOutput):
    """Parquet lift output."""

    schema = LIFT_SCHEMA
    accumulate_seq = False


#
# Subseq output handlers
#

class SubseqOutput(_OutputBase):
    """Abstract base for subseq output handlers."""


class SubseqOutputCsv(_TextTableOutput, SubseqOutput):
    """CSV/TSV subseq output (includes the `seq` column)."""

    fields = _SUBSEQ_FIELDS


class SubseqOutputBed(_BedOutput, SubseqOutput):
    """BED subseq output (includes the `seq` column). Unlifted regions are written with empty lift fields."""

    keep_unlifted = True
    fields = _SUBSEQ_FIELDS


class SubseqOutputParquet(_ParquetOutput, SubseqOutput):
    """Parquet subseq output with sequence-buffer temp offloading."""

    schema = SUBSEQ_SCHEMA
    accumulate_seq = True


class SubseqOutputFasta(SubseqOutput):
    """FASTA subseq output (record ID is the 1-based target coordinate).

    A region that did not lift is still written: the record is named "UNLIFTED_N" (N counting from 1), its sequence is
    a single "N", and every lift-derived tag is "NA". The output is never indexed, because record names are not
    guaranteed to be unique.
    """

    keep_unlifted = True

    def __init__(self, out_file: Optional[str] = None, reverse: bool = False, is_gz: bool = False) -> None:
        """Initialize a FASTA output handler.

        :param out_file: Output file path, or `None` for stdout.
        :param reverse: Matches the `reverse` used for the lift (selects the target coordinate).
        :param is_gz: bgzip-compress the output if `True`.

        :raises ValueError: If compressed FASTA is requested for stdout.
        """
        super().__init__(out_file, reverse)

        if is_gz and out_file is None:
            raise ValueError('Compressed FASTA output requires an output file')

        self.is_gz = is_gz
        self._handle = None
        self._close_handle = False
        self._unlifted_count = 0

    def _open(self) -> None:
        """Open the FASTA output."""
        if self.is_gz:
            self._handle = BGZFWriterIO(self.out_file)
            self._handle.__enter__()
            self._close_handle = True
        else:
            self._handle, self._close_handle = _open_text(self.out_file, False)

    def add(self, row: dict[str, Any]) -> None:
        """Write one FASTA record (an "N"-sequence record named "UNLIFTED_N" if the region did not lift)."""
        target = _target_struct(row, self.reverse)

        if target is None:
            self._unlifted_count += 1
            name = f'UNLIFTED_{self._unlifted_count}'
            seq = 'N'
        else:
            name = f'{target["chrom"]}:{target["pos"] + 1}-{target["end"]}'
            seq = row.get('seq') or ''

        tags = ' '.join([
            f'ref={_region_str(row.get("ref"), empty="NA")}',
            f'qry={_region_str(row.get("qry"), empty="NA")}',
            f'sample={_cell(row.get("sample"), empty="NA")}',
            f'hap={_cell(row.get("hap"), empty="NA")}',
            f'item={_cell(row.get("item"), empty="NA")}',
            f'strand={_cell(row.get("strand"), empty="NA")}',
            f'align_index={_align_index_str(row.get("align_index"), empty="NA")}',
        ])

        self._handle.write(f'>{name} {tags}\n{seq}\n')

    def _close(self) -> None:
        """Close the FASTA output (never indexed; record names are not guaranteed unique)."""
        if self._handle is not None and self._close_handle:
            self._handle.close()

        self._handle = None
