"""Tests for source-table building, format detection, and Align-mode lifting."""
from __future__ import annotations

import polars as pl
import pytest

from pav3.lift import (
    Format,
    LiftSourceAlign,
    detect_format,
    parse_format,
    source_table_from_file,
    source_table_from_list,
    source_table_from_str,
)


def test_source_table_from_list_string():
    """Region strings are 1-based and comma-tolerant."""
    df = source_table_from_list(['chr1:1,000-2,000'])
    assert df.to_dicts() == [{'chrom': 'chr1', 'pos': 999, 'end': 2000, 'item': 0}]


def test_source_table_from_list_tuple():
    """Tuples are 0-based half-open."""
    df = source_table_from_list([('chr1', 10, 20), ('chr2', 5, 8)])
    assert df['pos'].to_list() == [10, 5]
    assert df['item'].to_list() == [0, 1]


def test_source_table_stdin_sniff_bed():
    """A tab in the first line selects BED (0-based)."""
    df = source_table_from_str('chr1\t10\t20\nchr2\t0\t5\n')
    assert df.to_dicts()[0] == {'chrom': 'chr1', 'pos': 10, 'end': 20, 'item': 0}


def test_source_table_stdin_sniff_txt():
    """No tab in the first line selects a coordinate list (1-based)."""
    df = source_table_from_str('chr1:11-20\n')
    assert df.to_dicts() == [{'chrom': 'chr1', 'pos': 10, 'end': 20, 'item': 0}]


def test_source_table_item_overwrite_warning(tmp_path, caplog):
    """An existing item column is overwritten with a warning."""
    path = tmp_path / 'src.parquet'
    pl.DataFrame({'chrom': ['chr1'], 'pos': [10], 'end': [20], 'item': [99]}).write_parquet(path)

    with caplog.at_level('WARNING'):
        df = source_table_from_file(path)

    assert df['item'].to_list() == [0]
    assert any('item' in rec.message for rec in caplog.records)


def test_detect_format_rejects_gzipped_parquet():
    """A gzipped-Parquet filename is rejected (Parquet is not gzip-compressed)."""
    with pytest.raises(ValueError):
        detect_format('x.parquet.gz')


def test_source_table_bed_track_word_boundary(tmp_path):
    """A contig whose name merely starts with 'track'/'browser' is not mistaken for a declaration line."""
    path = tmp_path / 'src.bed'
    path.write_text('track name=x\ntrackpad\t10\t20\nbrowserCtg\t30\t40\n')

    df = source_table_from_file(path)
    assert df['chrom'].to_list() == ['trackpad', 'browserCtg']
    assert df['pos'].to_list() == [10, 30]


def test_source_table_bed_bad_coordinate_after_data_error(tmp_path):
    """A non-integer coordinate on a line after real data is an error, not a silent skip."""
    from pav3.lift.errors import LiftInputError

    path = tmp_path / 'src.bed'
    path.write_text('chr1\t10\t20\nchr2\tBAD\t40\n')

    with pytest.raises(LiftInputError):
        source_table_from_file(path)


def test_source_table_bed_spaces_error(tmp_path):
    """A space-delimited BED file fails immediately with an input error."""
    from pav3.lift.errors import LiftInputError

    path = tmp_path / 'src.bed'
    path.write_text('chr1 10 20\nchr2 30 40\n')

    with pytest.raises(LiftInputError):
        source_table_from_file(path)


def test_source_table_bed_empty_field_error(tmp_path):
    """A BED line with an empty column among the first three fails."""
    from pav3.lift.errors import LiftInputError

    path = tmp_path / 'src.bed'
    path.write_text('chr1\t\t20\n')

    with pytest.raises(LiftInputError):
        source_table_from_file(path)


def test_source_table_bed_header_skipped(tmp_path):
    """A BED header row (non-integer coordinates) is skipped, not an error."""
    path = tmp_path / 'src.bed'
    path.write_text('chrom\tstart\tend\nchr1\t10\t20\n')

    df = source_table_from_file(path)
    assert df.to_dicts() == [{'chrom': 'chr1', 'pos': 10, 'end': 20, 'item': 0}]


@pytest.mark.parametrize('suffix', ['bed', 'csv', 'tsv', 'parquet'])
def test_source_table_from_file_roundtrip(tmp_path, suffix):
    """Files in each tabular format round-trip to the same coordinates."""
    path = tmp_path / f'src.{suffix}'

    if suffix == 'bed':
        path.write_text('chr1\t10\t20\nchr2\t30\t40\n')
    elif suffix == 'parquet':
        pl.DataFrame({'chrom': ['chr1', 'chr2'], 'pos': [10, 30], 'end': [20, 40]}).write_parquet(path)
    else:
        sep = ',' if suffix == 'csv' else '\t'
        path.write_text(f'chrom{sep}pos{sep}end\nchr1{sep}10{sep}20\nchr2{sep}30{sep}40\n')

    df = source_table_from_file(path)
    assert df['chrom'].to_list() == ['chr1', 'chr2']
    assert df['pos'].to_list() == [10, 30]
    assert df['end'].to_list() == [20, 40]
    assert df['item'].to_list() == [0, 1]


def test_parse_format():
    """Format strings parse case-insensitively with an optional GZ suffix."""
    assert parse_format('txtgz') == (Format.TXT, True)
    assert parse_format('BED') == (Format.BED, False)

    with pytest.raises(ValueError):
        parse_format('PARQUETGZ')

    with pytest.raises(ValueError):
        parse_format('nope')


def test_detect_format():
    """Formats are detected from filename extensions."""
    assert detect_format('x.bed.gz') == (Format.BED, True)
    assert detect_format('x.parquet') == (Format.PARQUET, False)
    assert detect_format('x.unknown') is None


def test_lift_forward(align_fixture):
    """A reference region lifts to query coordinates on a forward-aligned contig."""
    with LiftSourceAlign(align_file=str(align_fixture.sam), no_seq=True) as src:
        rows = src.lift({'chrom': 'chr1', 'pos': 20, 'end': 40, 'item': 0})

    assert len(rows) == 1
    assert rows[0]['ref'] == {'chrom': 'chr1', 'pos': 20, 'end': 40}
    assert rows[0]['qry'] == {'chrom': 'ctg1', 'pos': 0, 'end': 20}
    assert rows[0]['strand'] == '+'
    assert rows[0]['align_index'] == {'pos': 0, 'end': 0}


def test_lift_reverse(align_fixture):
    """A query region lifts back to reference coordinates with -r/--reverse."""
    with LiftSourceAlign(align_file=str(align_fixture.sam), no_seq=True) as src:
        rows = src.lift({'chrom': 'ctg1', 'pos': 0, 'end': 20, 'item': 0}, reverse=True)

    assert len(rows) == 1
    assert rows[0]['ref'] == {'chrom': 'chr1', 'pos': 20, 'end': 40}


def test_lift_reverse_aligned_contig(align_fixture):
    """A region over a reverse-aligned contig reports the minus strand."""
    with LiftSourceAlign(align_file=str(align_fixture.sam), no_seq=True) as src:
        rows = src.lift({'chrom': 'chr1', 'pos': 70, 'end': 80, 'item': 0})

    assert len(rows) == 1
    assert rows[0]['strand'] == '-'
    assert rows[0]['qry']['chrom'] == 'ctg2'


def test_lift_unlifted_region(align_fixture):
    """A reference region outside any alignment yields no lift rows."""
    with LiftSourceAlign(align_file=str(align_fixture.sam), no_seq=True) as src:
        rows = src.lift({'chrom': 'chr1', 'pos': 0, 'end': 5, 'item': 0})

    assert rows == []


def test_lift_zero_length_region(align_fixture):
    """A zero-length reference region (pos == end) lifts to a zero-length query target, not a 1-bp span."""
    with LiftSourceAlign(align_file=str(align_fixture.sam), no_seq=True) as src:
        rows = src.lift({'chrom': 'chr1', 'pos': 30, 'end': 30, 'item': 0})

    assert len(rows) == 1
    assert rows[0]['qry']['pos'] == rows[0]['qry']['end']  # zero-length target
