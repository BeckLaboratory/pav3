"""Tests for Align-mode sequence extraction and output handlers."""
from __future__ import annotations

import polars as pl
import pysam

from pav3.lift import (
    LiftSourceAlign,
    SubseqOutputParquet,
    subseq_align,
)


def test_subseq_forward_matches_reference(align_fixture):
    """Forward extraction returns the query segment in reference orientation.

    Because both contigs align perfectly, the extracted sequence must equal the reference bases for both the
    forward-aligned and reverse-aligned contigs.
    """
    ref = align_fixture.ref

    with LiftSourceAlign(align_file=str(align_fixture.sam)) as src:
        fwd = src.lift({'chrom': 'chr1', 'pos': 20, 'end': 40, 'item': 0})
        rev = src.lift({'chrom': 'chr1', 'pos': 70, 'end': 80, 'item': 0})

        assert src.seq(fwd[0]) == ref[20:40]
        assert src.seq(rev[0]) == ref[70:80]


def test_subseq_reverse_needs_reference(align_fixture):
    """Reverse extraction without a reference raises an error."""
    from pav3.lift.errors import LiftInputError

    with LiftSourceAlign(align_file=str(align_fixture.sam)) as src:
        rows = src.lift({'chrom': 'ctg1', 'pos': 0, 'end': 20, 'item': 0}, reverse=True)

        try:
            src.seq(rows[0], reverse=True)
        except LiftInputError:
            pass
        else:
            raise AssertionError('Expected LiftInputError without a reference file')


def test_subseq_reverse_with_reference(align_fixture):
    """Reverse extraction returns the reference bases for the lifted query region."""
    ref = align_fixture.ref

    with LiftSourceAlign(align_file=str(align_fixture.sam), ref_file=str(align_fixture.ref_fa)) as src:
        rows = src.lift({'chrom': 'ctg1', 'pos': 0, 'end': 20, 'item': 0}, reverse=True)
        assert src.seq(rows[0], reverse=True) == ref[20:40]


def test_subseq_uppercase(align_fixture):
    """The uppercase option upper-cases extracted sequences."""
    with LiftSourceAlign(align_file=str(align_fixture.sam), uppercase=True) as src:
        rows = src.lift({'chrom': 'chr1', 'pos': 20, 'end': 40, 'item': 0})
        assert src.seq(rows[0]).isupper()


def test_subseq_align_fasta_output(align_fixture, tmp_path):
    """The subseq FASTA output names records by target coordinate and carries tags."""
    from pav3.lift import SubseqOutputFasta

    out_fa = tmp_path / 'out.fa'
    src_table = pl.DataFrame(
        {'chrom': ['chr1'], 'pos': [20], 'end': [40], 'item': pl.Series([0], dtype=pl.Int32)}
    )

    subseq_align(
        src_table,
        SubseqOutputFasta(out_file=str(out_fa)),
        align_file=str(align_fixture.sam),
    )

    with pysam.FastaFile(str(out_fa)) as fa_file:
        assert fa_file.references == ['ctg1:1-20']
        assert fa_file.fetch('ctg1:1-20') == align_fixture.ref[20:40]


def test_csv_output_quotes_comma_cells(tmp_path):
    """CSV output quotes region/align_index cells so embedded commas do not split columns."""
    from pav3.lift import SubseqOutputCsv

    out_csv = tmp_path / 'out.csv'
    row = {
        'sample': 'HG1', 'hap': 'h1', 'item': 0,
        'ref': {'chrom': 'chr1', 'pos': 1000, 'end': 2000},
        'qry': {'chrom': 'ctg1', 'pos': 5000, 'end': 6000},
        'strand': '+',
        'align_index': {'pos': 3, 'end': 4},
        'seq': 'ACGT',
    }

    handler = SubseqOutputCsv(out_file=str(out_csv), separator=',')
    with handler:
        handler.add(row)

    # Re-parse as CSV: the region cell (with thousands commas) must survive as a single field
    df = pl.read_csv(out_csv)
    assert df.height == 1
    assert df['ref'].to_list() == ['chr1:1,001-2,000']
    assert df['align_index'].to_list() == ['3,4']
    assert df['seq'].to_list() == ['ACGT']


def test_subseq_unlifted_fasta(align_fixture, tmp_path):
    """Unliftable regions are written to FASTA as single-"N" records named "UNLIFTED_N" (counting from 1)."""
    from pav3.lift import SubseqOutputFasta

    out_fa = tmp_path / 'out.fa'
    src_table = pl.DataFrame(
        {
            'chrom': ['chr1', 'chr1'],
            'pos': [0, 10],
            'end': [10, 20],
            'item': pl.Series([0, 1], dtype=pl.Int32),
        }
    )

    subseq_align(
        src_table,
        SubseqOutputFasta(out_file=str(out_fa)),
        align_file=str(align_fixture.sam),
    )

    lines = out_fa.read_text().splitlines()
    headers = [line for line in lines if line.startswith('>')]

    assert headers[0].startswith('>UNLIFTED_1 ')
    assert headers[1].startswith('>UNLIFTED_2 ')
    assert 'strand=NA' in headers[0]
    assert 'align_index=NA' in headers[0]
    assert 'qry=NA' in headers[0]          # forward mode: query (target) is undefined
    assert 'ref=chr1:1-10' in headers[0]   # forward mode: reference (source) is present
    assert lines[1] == 'N'

    # No FASTA index is written (record names are not guaranteed unique)
    assert not (tmp_path / 'out.fa.fai').exists()


def test_subseq_unlifted_bed(align_fixture, tmp_path):
    """An unliftable region is written to BED with empty target columns and empty lift fields."""
    from pav3.lift import SubseqOutputBed

    out_bed = tmp_path / 'out.bed'
    src_table = pl.DataFrame(
        {'chrom': ['chr1'], 'pos': [0], 'end': [10], 'item': pl.Series([0], dtype=pl.Int32)}
    )

    subseq_align(
        src_table,
        SubseqOutputBed(out_file=str(out_bed)),
        align_file=str(align_fixture.sam),
    )

    rows = [line for line in out_bed.read_text().splitlines() if not line.startswith('#')]
    assert len(rows) == 1
    cells = rows[0].split('\t')
    assert cells[0:3] == ['', '', '']        # empty target coordinate columns
    assert cells[3 + 2] == '0'               # item still present


def test_subseq_unlifted_parquet(align_fixture, tmp_path):
    """An unliftable region is written to Parquet with the item present and null lift fields."""
    out_pq = tmp_path / 'out.parquet'
    src_table = pl.DataFrame(
        {'chrom': ['chr1'], 'pos': [0], 'end': [10], 'item': pl.Series([0], dtype=pl.Int32)}
    )

    subseq_align(src_table, SubseqOutputParquet(out_file=str(out_pq)), align_file=str(align_fixture.sam))

    df = pl.read_parquet(out_pq)
    assert df.height == 1
    assert df['item'].to_list() == [0]
    assert df['qry'].to_list() == [None]     # forward mode: query (target) is null
    assert df['strand'].to_list() == [None]
    assert df['seq'].to_list() == [None]


def test_subseq_align_parquet_temp_chunk(align_fixture, tmp_path):
    """Parquet output with a tiny chunk threshold offloads to temp and still produces the full table."""
    out_pq = tmp_path / 'out.parquet'
    src_table = pl.DataFrame(
        {
            'chrom': ['chr1', 'chr1'],
            'pos': [20, 70],
            'end': [40, 80],
            'item': pl.Series([0, 1], dtype=pl.Int32),
        }
    )

    handler = SubseqOutputParquet(out_file=str(out_pq), tmp_dir=str(tmp_path), chunk_bp=1)

    subseq_align(src_table, handler, align_file=str(align_fixture.sam))

    df = pl.read_parquet(out_pq)
    assert df.height == 2
    assert set(df.columns) == {'sample', 'hap', 'item', 'ref', 'qry', 'strand', 'align_index', 'seq'}
    assert df['item'].to_list() == [0, 1]
    assert df['seq'].to_list()[0] == align_fixture.ref[20:40]
