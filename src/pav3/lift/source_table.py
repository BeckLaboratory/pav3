"""Build source coordinate tables for the lift and subseq tools.

A source table is a Polars DataFrame with columns `chrom`, `pos`, `end` (0-based half-open) and `item` (Int32, the
0-based input order). Coordinates may come from a list of region strings/tuples, a file, or a string (stdin).
"""

__all__ = [
    'SOURCE_SCHEMA',
    'source_table_from_list',
    'source_table_from_str',
    'source_table_from_file',
]

import io as _io
import logging
from pathlib import Path
from typing import Iterable, Optional

import polars as pl

from ..io import PlainOrGzFile
from ..region import region_from_string
from .errors import LiftInputError
from ..format import Format, detect_format, parse_format

logger = logging.getLogger(__name__)

SOURCE_SCHEMA: dict[str, pl.DataType] = {
    'chrom': pl.String,
    'pos': pl.Int64,
    'end': pl.Int64,
    'item': pl.Int32,
}
"""Schema of a source coordinate table."""

_COORD_SCHEMA = {'chrom': pl.String, 'pos': pl.Int64, 'end': pl.Int64}


def _empty_source_table() -> pl.DataFrame:
    """Get an empty source table with the standard schema.

    :returns: Empty source DataFrame.
    """
    return pl.DataFrame(schema=SOURCE_SCHEMA)


def _finalize_source_table(df: pl.DataFrame) -> pl.DataFrame:
    """Validate columns, drop/overwrite `item`, and add the input-order `item` column.

    :param df: Input DataFrame with at least `chrom`, `pos`, and `end` columns.

    :returns: Source table with `chrom`, `pos`, `end`, `item`.

    :raises ValueError: If required columns are missing.
    """
    if 'item' in df.columns:
        logger.warning('Input contains an "item" column; overwriting with input order')
        df = df.drop('item')

    if missing := {'chrom', 'pos', 'end'} - set(df.columns):
        raise ValueError(f'Source table missing required columns: {", ".join(sorted(missing))}')

    return (
        df
        .select(
            pl.col('chrom').cast(pl.String),
            pl.col('pos').cast(pl.Int64),
            pl.col('end').cast(pl.Int64),
        )
        .with_row_index('item')
        .with_columns(pl.col('item').cast(pl.Int32))
    )


def _txt_to_table(lines: Iterable[str]) -> pl.DataFrame:
    """Build a source table from a list of 1-based coordinate strings.

    Blank lines and lines starting with "#" are ignored.

    :param lines: Iterable of coordinate strings (e.g. "chr1:1,000-2,000").

    :returns: Source table.
    """
    rows = []

    for line in lines:
        line = line.strip()

        if not line or line.startswith('#'):
            continue

        region = region_from_string(line)
        rows.append({'chrom': region.chrom, 'pos': region.pos, 'end': region.end})

    if not rows:
        return _empty_source_table()

    return _finalize_source_table(pl.DataFrame(rows, schema=_COORD_SCHEMA))


def _bed_to_table(lines: Iterable[str]) -> pl.DataFrame:
    """Build a source table from BED lines (0-based half-open, first three columns).

    Comment lines ("#") and UCSC "track"/"browser" declaration lines (where "track" or "browser" is the first
    whitespace-delimited token) are ignored. Every other line must have at least three non-empty tab-separated fields
    (a line with fewer, or with an empty column among the first three, is an error; this catches space-delimited
    input). A single leading line whose second and third columns are not integers is treated as a header and skipped;
    a non-integer coordinate on any later line is an error (rather than being silently dropped).

    :param lines: Iterable of BED lines.

    :returns: Source table.

    :raises LiftInputError: If a data line does not have three non-empty tab-separated fields, or a non-header line has
        non-integer coordinates.
    """
    rows = []
    seen_data = False

    for line_num, line in enumerate(lines, start=1):
        line = line.rstrip('\n').rstrip('\r')

        stripped = line.strip()

        if not stripped or stripped.startswith('#') or stripped.split(None, 1)[0] in ('track', 'browser'):
            continue

        fields = line.split('\t')

        if len(fields) < 3 or not (fields[0].strip() and fields[1].strip() and fields[2].strip()):
            raise LiftInputError(
                f'BED line {line_num} does not have three non-empty tab-separated fields '
                f'(is the file tab-delimited?): {line!r}'
            )

        try:
            pos = int(fields[1])
            end = int(fields[2])
        except ValueError:
            if not seen_data:
                seen_data = True
                continue  # Leading header line

            raise LiftInputError(
                f'BED line {line_num} has non-integer coordinates: {line!r}'
            )

        seen_data = True
        rows.append({'chrom': fields[0], 'pos': pos, 'end': end})

    if not rows:
        return _empty_source_table()

    return _finalize_source_table(pl.DataFrame(rows, schema=_COORD_SCHEMA))


def source_table_from_list(regions: Iterable[str | tuple]) -> pl.DataFrame:
    """Build a source table from a list of region strings or (chrom, pos, end) tuples.

    Region strings are 1-based inclusive (UCSC/IGV). Tuples are 0-based half-open (BED), for convenient API use.

    :param regions: Iterable of region strings or (chrom, pos, end) tuples.

    :returns: Source table.

    :raises ValueError: If a region cannot be interpreted.
    """
    rows = []

    for region in regions:
        if isinstance(region, str):
            parsed = region_from_string(region)
            rows.append({'chrom': parsed.chrom, 'pos': parsed.pos, 'end': parsed.end})

        elif isinstance(region, (tuple, list)) and len(region) == 3:
            chrom, pos, end = region
            rows.append({'chrom': str(chrom), 'pos': int(pos), 'end': int(end)})

        else:
            raise ValueError(f'Unrecognized region: {region!r}')

    if not rows:
        return _empty_source_table()

    return _finalize_source_table(pl.DataFrame(rows, schema=_COORD_SCHEMA))


def source_table_from_str(text: str, file_type: Optional[str] = None) -> pl.DataFrame:
    """Build a source table from a string of TXT or BED data (used for stdin).

    When `file_type` is `None`, the format is auto-detected: if the first non-empty line contains a tab character, the
    input is treated as BED, otherwise it is treated as a list of coordinate strings (TXT). A `file_type` may specify
    any format except Parquet, and gzip-compressed input is not allowed.

    :param text: Input text.
    :param file_type: Format string (from -I) or `None` to auto-detect.

    :returns: Source table.

    :raises ValueError: If the format is unsupported for string input.
    """
    if file_type is not None:
        fmt, is_gz = parse_format(file_type)

        if is_gz:
            raise ValueError('Gzip-compressed input is not supported on standard input')

        if fmt == Format.PARQUET:
            raise ValueError('Parquet input is not supported on standard input')

        if fmt == Format.TXT:
            return _txt_to_table(text.splitlines())

        if fmt == Format.BED:
            return _bed_to_table(text.splitlines())

        if fmt in (Format.CSV, Format.TSV):
            sep = ',' if fmt == Format.CSV else '\t'
            return _finalize_source_table(pl.read_csv(_io.StringIO(text), separator=sep))

        raise ValueError(f'Unsupported input format for standard input: {fmt.value}')

    first_line = next((line for line in text.splitlines() if line.strip()), '')

    if '\t' in first_line:
        return _bed_to_table(text.splitlines())

    return _txt_to_table(text.splitlines())


def source_table_from_file(path: str | Path, file_type: Optional[str] = None) -> pl.DataFrame:
    """Build a source table from a file (TXT, BED, CSV, TSV, or Parquet).

    Text formats may be gzip-compressed. When `file_type` is `None`, the format is auto-detected from the filename.

    :param path: Input file path.
    :param file_type: Format string (from -I) or `None` to auto-detect from the filename.

    :returns: Source table.

    :raises ValueError: If the format cannot be determined or is unsupported.
    """
    path = Path(path)

    if file_type is not None:
        fmt, is_gz = parse_format(file_type)
    else:
        try:
            detected = detect_format(path)
        except ValueError as e:
            raise LiftInputError(str(e)) from e

        if detected is None:
            raise ValueError(f'Cannot determine input format from filename: {path}')

        fmt, is_gz = detected

    if fmt == Format.PARQUET:
        return _finalize_source_table(pl.read_parquet(path))

    if fmt in (Format.TXT, Format.BED):
        with PlainOrGzFile(str(path)) as in_file:
            lines = in_file.read().splitlines()

        return _txt_to_table(lines) if fmt == Format.TXT else _bed_to_table(lines)

    if fmt in (Format.CSV, Format.TSV):
        sep = ',' if fmt == Format.CSV else '\t'

        if is_gz:
            with PlainOrGzFile(str(path)) as in_file:
                df = pl.read_csv(_io.StringIO(in_file.read()), separator=sep)
        else:
            df = pl.read_csv(path, separator=sep)

        return _finalize_source_table(df)

    raise ValueError(f'Unsupported input format: {fmt.value}')
