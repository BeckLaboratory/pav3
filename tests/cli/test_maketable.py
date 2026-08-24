"""Tests for the maketable tool and the ``pav3 maketable`` subcommand."""
from __future__ import annotations

import gzip

import polars as pl
import pytest

from pav3.cli import main
from pav3.cli._returncode import ExitCode
from pav3.errors import PavInputError, PavUsageError
from pav3.maketable import bed_coord_check, is_simple_list, is_sorted, maketable, nested_to_string


@pytest.fixture
def var_parquet(tmp_path):
    """Write an unsorted variant-like table with nested columns and return its path."""
    path = tmp_path / 'var.parquet'

    pl.DataFrame({
        'chrom': ['chr2', 'chr1', 'chr1', 'chr10'],
        'pos': [100, 500, 100, 7],
        'end': [101, 600, 200, 8],
        'id': ['v0', 'v1', 'v2', 'v3'],
    }).with_columns(
        align_source=pl.Series([[1], [2, 3], [], None], dtype=pl.List(pl.Int32)),
        hom_ref=pl.Series(
            [{'up': 1, 'dn': 2}, {'up': 3, 'dn': 4}, {'up': 5, 'dn': 6}, None],
            dtype=pl.Struct({'up': pl.Int32, 'dn': pl.Int32}),
        ),
    ).write_parquet(path)

    return path


@pytest.fixture
def sorted_parquet(tmp_path):
    """Write a sorted BED-compatible table and return its path."""
    path = tmp_path / 'sorted.parquet'

    pl.DataFrame({
        'chrom': ['chr1', 'chr1', 'chr2'],
        'pos': [1, 2, 3],
        'end': [2, 3, 4],
        'x': [1.5, 2.5, 3.5],
    }).write_parquet(path)

    return path


@pytest.fixture
def nonbed_parquet(tmp_path):
    """Write a table without BED-compatible coordinate columns and return its path."""
    path = tmp_path / 'nonbed.parquet'

    pl.DataFrame({'a': [1, 2], 'b': ['x', 'y'], 'c': [3, 4]}).write_parquet(path)

    return path


def read_lines(path, is_gz=False):
    """Read an output file as a list of lines with trailing newlines stripped."""
    if is_gz:
        with gzip.open(path, 'rt') as in_file:  # BGZF is gzip-compatible
            text = in_file.read()
    else:
        text = path.read_text()

    return text.strip('\n').split('\n')


#
# Schema check
#

def test_bed_coord_check_accepts_any_integer_width():
    """Coordinate columns may be any integer type, signed or unsigned."""
    schema = pl.Schema({'chrom': pl.String, 'pos': pl.Int32, 'end': pl.UInt64, 'x': pl.Float64})

    bed_coord_check(schema)


def test_bed_coord_check_rejects_non_string_chrom():
    """A non-string first column is not a BED chrom."""
    with pytest.raises(ValueError, match='must be a string'):
        bed_coord_check(pl.Schema({'chrom': pl.Int64, 'pos': pl.Int64, 'end': pl.Int64}))


def test_bed_coord_check_rejects_non_integer_coordinate():
    """A non-integer coordinate column is not a BED coordinate."""
    with pytest.raises(ValueError, match='must be an integer'):
        bed_coord_check(pl.Schema({'chrom': pl.String, 'pos': pl.Int64, 'end': pl.Float64}))


def test_bed_coord_check_rejects_narrow_table():
    """A table with fewer than three columns has no BED coordinates."""
    with pytest.raises(ValueError, match='at least 3 columns'):
        bed_coord_check(pl.Schema({'chrom': pl.String, 'pos': pl.Int64}))


#
# Sort check
#

@pytest.mark.parametrize(
    'chrom,pos,end,expected',
    [
        (['chr1', 'chr1', 'chr2'], [1, 2, 1], [2, 3, 2], True),
        (['chr1', 'chr1'], [5, 5], [10, 20], True),      # Ties on chrom and pos break on end
        (['chr1', 'chr1'], [5, 5], [20, 10], False),     # End out of order
        (['chr2', 'chr1'], [1, 1], [2, 2], False),       # Chrom out of order
        (['chr1', 'chr1'], [5, 1], [10, 2], False),      # Pos out of order
        ([None, 'chr1'], [1, 1], [2, 2], True),          # Nulls sort first
        (['chr1', None], [1, 1], [2, 2], False),
    ],
)
def test_is_sorted(tmp_path, chrom, pos, end, expected):
    """Sortedness is detected for the first three columns, including null and tie cases."""
    path = tmp_path / 'x.parquet'
    pl.DataFrame({'chrom': chrom, 'pos': pos, 'end': end}).write_parquet(path)

    assert is_sorted(pl.scan_parquet(path), ['chrom', 'pos', 'end']) is expected


#
# Nested columns
#

@pytest.mark.parametrize(
    'dtype,expected',
    [
        (pl.List(pl.String), True),
        (pl.List(pl.Int32), True),
        (pl.Array(pl.Int32, 2), True),
        (pl.List(pl.Struct({'a': pl.Int32})), False),
        (pl.List(pl.List(pl.Int32)), False),
        (pl.Struct({'a': pl.Int32}), False),
    ],
)
def test_is_simple_list(dtype, expected):
    """Lists and arrays of scalars are simple; structs and lists of nested values are not."""
    assert is_simple_list(dtype) is expected


def test_nested_to_string_joins_simple_lists(tmp_path):
    """Lists and arrays of scalars are joined with a comma."""
    path = tmp_path / 'x.parquet'

    pl.DataFrame({'i': [0, 1, 2]}).with_columns(
        s=pl.Series([['A', 'B'], [], None], dtype=pl.List(pl.String)),
        n=pl.Series([[1, 2], [3], None], dtype=pl.List(pl.Int32)),
        arr=pl.Series([[1, 2], [3, 4], [5, 6]], dtype=pl.Array(pl.Int32, 2)),
    ).write_parquet(path)

    df = nested_to_string(pl.scan_parquet(path)).collect()

    assert df['s'].to_list() == ['A,B', '', None]
    assert df['n'].to_list() == ['1,2', '3', None]
    assert df['arr'].to_list() == ['1,2', '3,4', '5,6']


def test_nested_to_string_drops_null_elements(tmp_path):
    """A null element inside a list is dropped rather than nulling the whole cell."""
    path = tmp_path / 'x.parquet'

    pl.DataFrame({'i': [0]}).with_columns(
        s=pl.Series([['A', None, 'B']], dtype=pl.List(pl.String)),
    ).write_parquet(path)

    assert nested_to_string(pl.scan_parquet(path)).collect()['s'].to_list() == ['A,B']


def test_nested_to_string_json_encodes_deep_types(tmp_path):
    """Structs and lists of structs are JSON-encoded; nulls stay null."""
    path = tmp_path / 'x.parquet'

    pl.DataFrame({'i': [0, 1]}).with_columns(
        st=pl.Series([{'a': 1}, None], dtype=pl.Struct({'a': pl.Int32})),
        lst_st=pl.Series([[{'a': 1}], []], dtype=pl.List(pl.Struct({'a': pl.Int32}))),
        lst_lst=pl.Series([[[1, 2]], None], dtype=pl.List(pl.List(pl.Int32))),
    ).write_parquet(path)

    df = nested_to_string(pl.scan_parquet(path)).collect()

    assert df.columns == ['i', 'st', 'lst_st', 'lst_lst']
    assert df['st'].to_list() == ['{"a":1}', None]
    assert df['lst_st'].to_list() == ['[{"a":1}]', '[]']
    assert df['lst_lst'].to_list() == ['[[1,2]]', None]


def test_nested_to_string_leaves_scalar_tables_untouched(sorted_parquet):
    """A table with no nested columns is returned unchanged."""
    lf = pl.scan_parquet(sorted_parquet)

    assert nested_to_string(lf).collect_schema() == lf.collect_schema()


#
# Output formats
#

def test_bed_output(tmp_path, var_parquet):
    """BED output is tab-delimited, sorted, and carries a commented header."""
    out_file = tmp_path / 'out.bed'
    maketable(var_parquet, out_file)

    lines = read_lines(out_file)

    assert lines[0] == '#chrom\tpos\tend\tid\talign_source\thom_ref'
    assert [line.split('\t')[3] for line in lines[1:]] == ['v2', 'v1', 'v3', 'v0']
    assert all(len(line.split('\t')) == 6 for line in lines)


def test_bed_output_nested_and_null_cells(tmp_path, var_parquet):
    """Nested cells are JSON and null cells are empty, with no quoting in tab-delimited output."""
    out_file = tmp_path / 'out.bed'
    maketable(var_parquet, out_file)

    cells = {line.split('\t')[3]: line.split('\t') for line in read_lines(out_file)[1:]}

    assert cells['v1'][4] == '2,3'                    # List of scalars is joined
    assert cells['v1'][5] == '{"up":3,"dn":4}'        # Struct is JSON
    assert cells['v2'][4] == ''                       # Empty list joins to an empty cell
    assert cells['v3'][4] == ''                       # Null list
    assert cells['v3'][5] == ''                       # Null struct
    assert '"' not in cells['v1'][4]


def test_csv_output_is_quoted(tmp_path, var_parquet):
    """CSV output quotes cells containing the separator (RFC 4180) and has a plain header."""
    out_file = tmp_path / 'out.csv'
    maketable(var_parquet, out_file)

    lines = read_lines(out_file)

    assert lines[0] == 'chrom,pos,end,id,align_source,hom_ref'
    assert '"2,3"' in '\n'.join(lines)          # Joined list quoted because it holds the CSV separator
    assert '"{""up"":3,""dn"":4}"' in '\n'.join(lines)


def test_tsv_output(tmp_path, sorted_parquet):
    """TSV output is tab-delimited with a plain (uncommented) header."""
    out_file = tmp_path / 'out.tsv'
    maketable(sorted_parquet, out_file)

    assert read_lines(out_file)[0] == 'chrom\tpos\tend\tx'


@pytest.mark.parametrize('suffix', ['bed.gz', 'csv.gz', 'tsv.gz'])
def test_gz_output_is_bgzf(tmp_path, sorted_parquet, suffix):
    """Compressed output is bgzip-compressed and readable as BGZF."""
    out_file = tmp_path / f'out.{suffix}'
    maketable(sorted_parquet, out_file)

    assert out_file.read_bytes()[:4] == b'\x1f\x8b\x08\x04'  # BGZF magic (gzip with extra field)
    assert len(read_lines(out_file, is_gz=True)) == 4


def test_xlsx_output(tmp_path, var_parquet):
    """Excel output round-trips through a reader with nested columns as JSON strings."""
    out_file = tmp_path / 'out.xlsx'
    maketable(var_parquet, out_file)

    df = pl.read_excel(out_file)

    assert df.height == 4
    assert df.columns == ['chrom', 'pos', 'end', 'id', 'align_source', 'hom_ref']


#
# Sorting
#

def test_bed_sorts_by_default(tmp_path, var_parquet):
    """BED output is sorted by coordinate without an explicit option."""
    out_file = tmp_path / 'out.bed'
    maketable(var_parquet, out_file)

    coords = [line.split('\t')[:3] for line in read_lines(out_file)[1:]]

    assert coords == [['chr1', '100', '200'], ['chr1', '500', '600'], ['chr10', '7', '8'], ['chr2', '100', '101']]


def test_csv_does_not_sort_by_default(tmp_path, var_parquet):
    """Output types other than BED keep the input order unless sorting is requested."""
    out_file = tmp_path / 'out.csv'
    maketable(var_parquet, out_file)

    assert [line.split(',')[0] for line in read_lines(out_file)[1:]] == ['chr2', 'chr1', 'chr1', 'chr10']


def test_csv_sorts_when_requested(tmp_path, var_parquet):
    """An explicit sort applies to output types that do not sort by default."""
    out_file = tmp_path / 'out.csv'
    maketable(var_parquet, out_file, sort=True)

    assert [line.split(',')[0] for line in read_lines(out_file)[1:]] == ['chr1', 'chr1', 'chr10', 'chr2']


def test_bed_no_sort_keeps_input_order(tmp_path, var_parquet):
    """Sorting can be turned off for BED output."""
    out_file = tmp_path / 'out.bed'
    maketable(var_parquet, out_file, sort=False)

    assert [line.split('\t')[0] for line in read_lines(out_file)[1:]] == ['chr2', 'chr1', 'chr1', 'chr10']


def test_sorted_input_skips_sort(tmp_path, sorted_parquet, caplog):
    """A table that is already sorted is not sorted again."""
    out_file = tmp_path / 'out.bed'

    with caplog.at_level('INFO', logger='pav3.maketable'):
        maketable(sorted_parquet, out_file)

    assert any('already sorted' in record.message for record in caplog.records)


#
# Coordinate check by output type and sort option
#

def test_bed_requires_coordinates(tmp_path, nonbed_parquet):
    """BED output of a table without coordinate columns is an error."""
    with pytest.raises(PavUsageError, match='not BED-compatible'):
        maketable(nonbed_parquet, tmp_path / 'out.bed')


def test_bed_requires_coordinates_even_without_sort(tmp_path, nonbed_parquet):
    """The coordinate check still applies to BED output when sorting is off."""
    with pytest.raises(PavUsageError, match='not BED-compatible'):
        maketable(nonbed_parquet, tmp_path / 'out.bed', sort=False)


def test_sort_requires_coordinates(tmp_path, nonbed_parquet):
    """Sorting a table without coordinate columns is an error for any output type."""
    with pytest.raises(PavUsageError, match='not BED-compatible'):
        maketable(nonbed_parquet, tmp_path / 'out.csv', sort=True)


def test_csv_without_sort_skips_coordinate_check(tmp_path, nonbed_parquet):
    """A table without coordinate columns is written when no sort is requested."""
    out_file = tmp_path / 'out.csv'
    maketable(nonbed_parquet, out_file)

    assert read_lines(out_file)[0] == 'a,b,c'


#
# Output type resolution
#

@pytest.mark.parametrize(
    'output_type,is_gz',
    [('bed', False), ('BED', False), ('bedgz', True), ('BEDGZ', True)],
)
def test_output_type_is_case_insensitive(tmp_path, sorted_parquet, output_type, is_gz):
    """Output types are not case sensitive and accept a GZ suffix."""
    out_file = tmp_path / 'out.dat'
    maketable(sorted_parquet, out_file, output_type=output_type)

    assert read_lines(out_file, is_gz=is_gz)[0].startswith('#chrom')


@pytest.mark.parametrize('output_type', ['xlsx', 'excel', 'EXCEL'])
def test_excel_alias(tmp_path, sorted_parquet, output_type):
    """Excel output accepts both the XLSX and EXCEL type names."""
    out_file = tmp_path / 'out.dat'
    maketable(sorted_parquet, out_file, output_type=output_type)

    assert pl.read_excel(out_file).height == 3


def test_unknown_extension(tmp_path, sorted_parquet):
    """An output type that cannot be detected from the filename is a usage error."""
    with pytest.raises(PavUsageError, match='Cannot determine the output type'):
        maketable(sorted_parquet, tmp_path / 'out.dat')


def test_unrecognized_output_type(tmp_path, sorted_parquet):
    """An unrecognized output type is a usage error."""
    with pytest.raises(PavUsageError, match='Unrecognized format'):
        maketable(sorted_parquet, tmp_path / 'out.bed', output_type='nope')


def test_unsupported_output_type(tmp_path, sorted_parquet):
    """Parquet is a valid format but not a maketable output type."""
    with pytest.raises(PavUsageError, match='Unsupported output type'):
        maketable(sorted_parquet, tmp_path / 'out.parquet')


def test_gzipped_excel_rejected(tmp_path, sorted_parquet):
    """Excel output cannot be gzip-compressed."""
    with pytest.raises(PavUsageError, match='cannot be gzip-compressed'):
        maketable(sorted_parquet, tmp_path / 'out.xlsx.gz')


def test_missing_input(tmp_path):
    """A missing input file is an input error."""
    with pytest.raises(PavInputError, match='missing or is not a regular file'):
        maketable(tmp_path / 'nope.parquet', tmp_path / 'out.bed')


def test_uncompressed_type_with_gz_name_warns(tmp_path, sorted_parquet, caplog):
    """An uncompressed output type written to a ".gz" filename warns but still writes."""
    out_file = tmp_path / 'out.bed.gz'

    with caplog.at_level('WARNING', logger='pav3.maketable'):
        maketable(sorted_parquet, out_file, output_type='BED')

    assert any('ends with ".gz"' in record.message for record in caplog.records)
    assert read_lines(out_file)[0].startswith('#chrom')  # Plain text despite the name


def test_compressed_type_with_plain_name_warns(tmp_path, sorted_parquet, caplog):
    """A compressed output type written to a filename without ".gz" warns but still writes."""
    out_file = tmp_path / 'out.bed'

    with caplog.at_level('WARNING', logger='pav3.maketable'):
        maketable(sorted_parquet, out_file, output_type='BEDGZ')

    assert any('does not end with ".gz"' in record.message for record in caplog.records)
    assert read_lines(out_file, is_gz=True)[0].startswith('#chrom')


#
# Multiple inputs
#

def test_glob_input(tmp_path):
    """A glob input reads and sorts across all matching files."""
    pl.DataFrame({'chrom': ['chr1', 'chr3'], 'pos': [5, 1], 'end': [6, 2]}).write_parquet(tmp_path / 'a.parquet')
    pl.DataFrame({'chrom': ['chr2'], 'pos': [9], 'end': [10]}).write_parquet(tmp_path / 'b.parquet')

    out_file = tmp_path / 'out.bed'
    maketable(str(tmp_path / '*.parquet'), out_file)

    assert [line.split('\t')[0] for line in read_lines(out_file)[1:]] == ['chr1', 'chr2', 'chr3']


#
# CLI
#

def test_cli_writes_bed(tmp_path, var_parquet):
    """The subcommand writes a BED file and returns OK."""
    out_file = tmp_path / 'out.bed'

    assert main(['maketable', '-i', str(var_parquet), '-o', str(out_file)]) == int(ExitCode.OK)
    assert read_lines(out_file)[0].startswith('#chrom')


def test_cli_no_sort_flag(tmp_path, var_parquet):
    """-S turns off sorting for BED output."""
    out_file = tmp_path / 'out.bed'

    assert main(['maketable', '-i', str(var_parquet), '-o', str(out_file), '-S']) == int(ExitCode.OK)
    assert [line.split('\t')[0] for line in read_lines(out_file)[1:]] == ['chr2', 'chr1', 'chr1', 'chr10']


def test_cli_sort_flag(tmp_path, var_parquet):
    """-s turns on sorting for output types that do not sort by default."""
    out_file = tmp_path / 'out.csv'

    assert main(['maketable', '-i', str(var_parquet), '-o', str(out_file), '-s']) == int(ExitCode.OK)
    assert [line.split(',')[0] for line in read_lines(out_file)[1:]] == ['chr1', 'chr1', 'chr10', 'chr2']


def test_cli_output_type_flag(tmp_path, sorted_parquet):
    """-O overrides the type detected from the output filename."""
    out_file = tmp_path / 'out.dat'

    assert main(['maketable', '-i', str(sorted_parquet), '-o', str(out_file), '-O', 'tsv']) == int(ExitCode.OK)
    assert read_lines(out_file)[0] == 'chrom\tpos\tend\tx'


def test_cli_sort_flags_are_mutually_exclusive(tmp_path, sorted_parquet):
    """-s and -S cannot be combined."""
    with pytest.raises(SystemExit) as exc_info:
        main(['maketable', '-i', str(sorted_parquet), '-o', str(tmp_path / 'out.bed'), '-s', '-S'])

    assert exc_info.value.code == int(ExitCode.USAGE)


def test_cli_missing_input_exit_code(tmp_path):
    """A missing input file returns the INPUT exit code."""
    args = ['maketable', '-i', str(tmp_path / 'nope.parquet'), '-o', str(tmp_path / 'out.bed')]

    assert main(args) == int(ExitCode.INPUT)


def test_cli_bad_output_type_exit_code(tmp_path, sorted_parquet):
    """An unrecognized output type returns the USAGE exit code."""
    args = ['maketable', '-i', str(sorted_parquet), '-o', str(tmp_path / 'out.bed'), '-O', 'nope']

    assert main(args) == int(ExitCode.USAGE)


def test_cli_incompatible_schema_exit_code(tmp_path, nonbed_parquet):
    """BED output of a table without coordinate columns returns the USAGE exit code."""
    args = ['maketable', '-i', str(nonbed_parquet), '-o', str(tmp_path / 'out.bed')]

    assert main(args) == int(ExitCode.USAGE)
