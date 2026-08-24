"""Convert Parquet tables to BED, CSV, TSV, and Excel.

This module backs the ``pav3 maketable`` subcommand. It is designed for PAV variant and alignment tables, but works
with any Parquet file. Text output is written with the Polars streaming engine and sunk directly to the output file
(bgzip-compressed when the output is gzipped), so tables larger than memory can be converted.

Non-scalar columns (list, array, and struct) have no text representation in delimited output, so they are converted
to strings before writing: a list of scalars is joined with ",", and anything deeper is JSON-encoded (see
:func:`nested_to_string`).
"""

__all__ = [
    'BED_COORD_FIELDS',
    'EXCEL_MAX_ROWS',
    'LIST_SEPARATOR',
    'bed_coord_check',
    'is_simple_list',
    'is_sorted',
    'maketable',
    'nested_to_string',
    'sort_check_expr',
]

import logging
from pathlib import Path
from typing import Optional

import polars as pl

from .errors import PavInputError, PavUsageError
from .format import Format, detect_format, parse_format
from .io import BGZFWriterIO

logger = logging.getLogger(__name__)


BED_COORD_FIELDS = ('chrom', 'pos', 'end')
"""Roles of the first three columns of a BED-compatible table."""

EXCEL_MAX_ROWS = 1_048_576
"""Maximum number of rows in an Excel worksheet (including the header row)."""

LIST_SEPARATOR = ','
"""Separator joining the elements of a list or array of scalars."""

_NESTED_TYPES = (pl.List, pl.Array, pl.Struct)

_JSON_WRAP_FIELD = 'v'

_JSON_PREFIX_LEN = len(f'{{"{_JSON_WRAP_FIELD}":')

_OUTPUT_FORMATS = {Format.BED, Format.CSV, Format.TSV, Format.XLSX}
"""Output formats supported by maketable."""

_SEPARATOR = {
    Format.BED: '\t',
    Format.TSV: '\t',
    Format.CSV: ',',
}
"""Field separator by output format."""


def bed_coord_check(
        schema: pl.Schema,
) -> None:
    """Check that the first three columns of a table are BED coordinates (string, integer, integer).

    Integer columns may be of any integer type (signed or unsigned, any width). Columns after the first three are not
    checked.

    :param schema: Table schema.

    :raises ValueError: If the table has fewer than three columns or the first three columns are not a string
        followed by two integers.
    """
    names = schema.names()

    if len(names) < 3:
        raise ValueError(
            f'Table must have at least 3 columns for BED coordinates '
            f'(chrom, pos, end): Found {len(names)}'
        )

    dtypes = [schema[name] for name in names[:3]]

    if dtypes[0] != pl.String:
        raise ValueError(
            f'Column 1 ("{names[0]}") must be a string for BED coordinates (chrom): Found {dtypes[0]}'
        )

    for index in (1, 2):
        if not dtypes[index].is_integer():
            raise ValueError(
                f'Column {index + 1} ("{names[index]}") must be an integer for BED coordinates '
                f'({BED_COORD_FIELDS[index]}): Found {dtypes[index]}'
            )


def sort_check_expr(
        columns: list[str],
) -> pl.Expr:
    """Get an expression testing whether a table is sorted by `columns` in ascending order.

    Polars cannot compare struct values, so the lexicographic comparison against the preceding row is built up from
    per-column comparisons. Null values sort first (matching `sort()` with `nulls_last=False`).

    :param columns: Sort columns in order of precedence.

    :returns: A boolean expression that is `True` if the table is sorted.
    """
    def _le(index: int) -> pl.Expr:
        """Test that row n-1 is ordered at or before row n starting at sort column `index`."""
        col = pl.col(columns[index])
        prev = col.shift(1)

        # Comparisons against a null are null, which would be swallowed by all(), so each is guarded to yield a
        # boolean. Nulls sort first: a null in the preceding row is ordered at or before any value, and a null in
        # the current row is out of order unless the preceding row is also null.
        both = prev.is_not_null() & col.is_not_null()

        lt = (prev.is_null() & col.is_not_null()) | (both & (prev < col))
        eq = (prev.is_null() & col.is_null()) | (both & (prev == col))

        if index + 1 == len(columns):
            return lt | eq

        return lt | (eq & _le(index + 1))

    # The first row has no preceding row: shift() makes every column null, which the null clauses order first
    return _le(0).all()


def is_sorted(
        lf: pl.LazyFrame,
        columns: list[str],
) -> bool:
    """Determine whether a table is already sorted by `columns`.

    Only the sort columns are read, and Parquet projection pushdown keeps this pass inexpensive relative to a sort
    (which must materialize the whole table).

    :param lf: Input table.
    :param columns: Sort columns in order of precedence.

    :returns: `True` if the table is sorted.
    """
    return bool(
        lf
        .select(columns)
        .select(sort_check_expr(columns).alias('is_sorted'))
        .collect(engine='streaming')
        .item()
    )


def is_simple_list(
        dtype: pl.DataType,
) -> bool:
    """Determine whether a data type is a list or array of non-nested elements.

    :param dtype: Column data type.

    :returns: `True` for a list or array whose elements are scalars (e.g. a list of strings or integers), and `False`
        for a struct or for a list of lists or structs.
    """
    return isinstance(dtype, (pl.List, pl.Array)) and not dtype.inner.is_nested()


def nested_to_string(
        lf: pl.LazyFrame,
) -> pl.LazyFrame:
    """Convert all non-scalar (list, array, and struct) columns of a table to strings.

    Delimited text and Excel output have no representation for nested values, and Polars cannot cast a list to a
    string, so each nested column is rendered one of two ways:

    - A list or array of scalars is joined with ",", which is how PAV renders `filter`, `discord`, and `inner`
      elsewhere (e.g. in tracks). A null element is dropped rather than nulling the whole cell.
    - Anything else (a struct, or a list of lists or structs) is JSON-encoded. The column is wrapped in a
      single-field struct, JSON-encoded, and the wrapper is stripped, which yields a JSON representation of the value
      itself for any nesting depth. JSON escapes separators and newlines, so an encoded value never breaks the
      delimited output.

    Null values stay null in both cases rather than rendering as the string "null". A joined empty list renders as an
    empty string, so it is indistinguishable from a null list; JSON-encoded values keep that distinction (`[]` vs an
    empty cell).

    :param lf: Input table.

    :returns: Table with nested columns replaced by string columns. The original column order is preserved.
    """
    schema = lf.collect_schema()

    nested = [
        name for name in schema.names()
        if isinstance(schema[name], _NESTED_TYPES)
    ]

    if not nested:
        return lf

    joined = [name for name in nested if is_simple_list(schema[name])]

    logger.debug(
        'Converting %d non-scalar column(s) to strings: %d joined (%s), %d JSON-encoded (%s)',
        len(nested),
        len(joined), ', '.join(joined) if joined else '-',
        len(nested) - len(joined), ', '.join(name for name in nested if name not in set(joined)) or '-',
    )

    return lf.with_columns([
        (
            pl.when(pl.col(name).is_null())
            .then(None)  # A null value renders as "null" or ""; write it as a null cell like any other column
            .otherwise(_nested_expr(name, schema[name]))
            .alias(name)
        )
        for name in nested
    ])


def _nested_expr(
        name: str,
        dtype: pl.DataType,
) -> pl.Expr:
    """Get the expression converting one non-scalar column to a string.

    :param name: Column name.
    :param dtype: Column data type.

    :returns: A string expression for the column (see :func:`nested_to_string`).
    """
    if is_simple_list(dtype):
        return (
            pl.col(name)
            .cast(pl.List(pl.String))
            .list.join(LIST_SEPARATOR, ignore_nulls=True)
        )

    return (
        pl.struct(pl.col(name).alias(_JSON_WRAP_FIELD))
        .struct.json_encode()
        .str.slice(_JSON_PREFIX_LEN)  # Strip the wrapper prefix: {"v":
        .str.head(-1)                 # Strip the wrapper suffix: }
    )


def _resolve_format(
        output: str,
        output_type: Optional[str],
) -> tuple[Format, bool]:
    """Resolve the output format and compression from the format option or the output filename.

    :param output: Output file name.
    :param output_type: Output format string (from -O), or `None` to detect from the filename.

    :returns: A tuple of (format, is_gz).

    :raises PavUsageError: If the format is not recognized or is not a supported output format.
    """
    is_gz_name = str(output).strip().lower().endswith('.gz')

    if output_type is not None:
        try:
            fmt, is_gz = parse_format(output_type)
        except ValueError as e:
            raise PavUsageError(str(e)) from e

        if is_gz and not is_gz_name:
            logger.warning(
                'Output type "%s" is compressed, but the output filename does not end with ".gz": %s',
                output_type, output
            )

        if not is_gz and is_gz_name:
            logger.warning(
                'Output type "%s" is not compressed, but the output filename ends with ".gz": %s',
                output_type, output
            )

    else:
        try:
            detected = detect_format(output)
        except ValueError as e:
            raise PavUsageError(str(e)) from e

        if detected is None:
            raise PavUsageError(
                f'Cannot determine the output type from the filename (use -O/--output-type): {output}'
            )

        fmt, is_gz = detected

    if fmt not in _OUTPUT_FORMATS:
        raise PavUsageError(
            f'Unsupported output type for maketable: {fmt.value} '
            f'(expected one of: {", ".join(sorted(f.value for f in _OUTPUT_FORMATS))})'
        )

    return fmt, is_gz


def _sink_text(
        lf: pl.LazyFrame,
        output: str,
        fmt: Format,
        is_gz: bool,
) -> None:
    """Sink a table to a delimited text file, optionally bgzip-compressed.

    BED output is written with a commented header line ("#chrom") so the file is both self-describing and readable by
    tools that skip "#" lines. Tab-delimited output is written unquoted, because RFC 4180 quoting is not understood by
    BED readers: nested values are JSON-encoded, which escapes tabs and newlines, so an encoded cell never contains
    the separator (a string column holding a raw tab is written as-is and would not be valid BED or TSV either way).
    CSV output uses RFC 4180 quoting.

    :param lf: Table to write.
    :param output: Output file name.
    :param fmt: Output format.
    :param is_gz: Write bgzip-compressed output if `True`.
    """
    separator = _SEPARATOR[fmt]

    if fmt == Format.BED:
        name = lf.collect_schema().names()[0]

        if not name.startswith('#'):
            lf = lf.rename({name: f'#{name}'})

    sink_args = {
        'separator': separator,
        'include_header': True,
        'null_value': '',
        'quote_style': 'necessary' if fmt == Format.CSV else 'never',
        'check_extension': False,  # Compression is resolved from the output type, not the filename
        'engine': 'streaming',
    }

    if is_gz:
        with BGZFWriterIO(output) as out_file:
            lf.sink_csv(out_file, **sink_args)
    else:
        lf.sink_csv(output, **sink_args)


def _write_excel(
        lf: pl.LazyFrame,
        output: str,
) -> None:
    """Write a table to an Excel workbook.

    Excel output cannot be streamed; the table is collected into memory before writing.

    :param lf: Table to write.
    :param output: Output file name.
    """
    df = lf.collect(engine='streaming')

    if df.height + 1 > EXCEL_MAX_ROWS:
        logger.warning(
            'Table has %d rows, which exceeds the Excel worksheet limit of %d rows (including the header): %s',
            df.height, EXCEL_MAX_ROWS, output
        )

    df.write_excel(output)


def maketable(
        input: str | Path,
        output: str | Path,
        output_type: Optional[str] = None,
        sort: Optional[bool] = None,
) -> None:
    """Convert a Parquet table to a BED, CSV, TSV, or Excel file.

    The first three columns of a BED file are the coordinates (chrom, pos, end), and PAV tables are already in this
    order with 0-based half-open coordinates, so no coordinate transformation is applied. The same three columns are
    the sort key.

    Sorting is on by default for BED output and off by default for all other output types. The coordinate columns are
    checked whenever sorting is enabled or the output is BED (a table that cannot be sorted by coordinates is not a
    valid BED file). If the table is already sorted, which PAV tables always are, the sort is skipped.

    :param input: Input Parquet file. May be a glob matching multiple files with the same schema.
    :param output: Output file name.
    :param output_type: Output type (BED, BEDGZ, CSV, CSVGZ, TSV, TSVGZ, XLSX/EXCEL; not case sensitive). If `None`,
        the type is detected from the output filename.
    :param sort: Sort by the first three columns. If `None`, sort BED output and do not sort other output types.

    :raises PavInputError: If the input file does not exist.
    :raises PavUsageError: If the output type cannot be determined or is unsupported, or if sorting or BED output was
        requested for a table without BED-compatible coordinate columns.
    """
    input = str(input)
    output = str(output)

    if not _is_glob(input) and not Path(input).is_file():
        raise PavInputError(f'Input Parquet file missing or is not a regular file: {input}')

    fmt, is_gz = _resolve_format(output, output_type)

    if sort is None:
        sort = fmt == Format.BED

    logger.info('Reading Parquet table: %s', input)

    lf = pl.scan_parquet(input)
    schema = lf.collect_schema()

    # Check coordinate columns if the table is sorted by them or written as BED
    if sort or fmt == Format.BED:
        try:
            bed_coord_check(schema)
        except ValueError as e:
            raise PavUsageError(
                f'Table is not BED-compatible ({"sorting" if sort else "BED output"} requires '
                f'chrom, pos, and end in the first three columns): {e}'
            ) from e

    if sort:
        sort_columns = schema.names()[:3]

        if is_sorted(lf, sort_columns):
            logger.info('Table is already sorted by (%s)', ', '.join(sort_columns))
        else:
            logger.info('Sorting by (%s)', ', '.join(sort_columns))
            lf = lf.sort(sort_columns)

    lf = nested_to_string(lf)

    logger.info('Writing %s%s: %s', fmt.value, ' (bgzip)' if is_gz else '', output)

    if fmt == Format.XLSX:
        _write_excel(lf, output)
    else:
        _sink_text(lf, output, fmt, is_gz)

    logger.info('Done')


def _is_glob(
        path: str,
) -> bool:
    """Determine whether a path contains glob wildcards.

    :param path: Input path.

    :returns: `True` if the path contains glob wildcards.
    """
    return any(char in path for char in '*?[')
