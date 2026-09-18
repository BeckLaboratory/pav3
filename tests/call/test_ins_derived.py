"""Tests for duplications derived from insertions."""
from __future__ import annotations

import polars as pl
import pysam
import pytest

from pav3.call.integrate import add_ins_derived
from pav3.vcf import standard_info_fields


REF_SEQ = 'ACGT' * 250

_DUP_TYPE = pl.List(pl.Struct({'chrom': pl.String, 'pos': pl.Int64, 'end': pl.Int64, 'is_rev': pl.Boolean}))


@pytest.fixture
def ref_fa(tmp_path):
    """Write a small indexed reference FASTA."""
    path = tmp_path / 'ref.fa'
    path.write_text(f'>chr1\n{REF_SEQ}\n')
    pysam.faidx(str(path))

    return path


def _insdel_frame(vartype: str, varsubtype: str | None, dup: list | None) -> pl.LazyFrame:
    """Build a one-row insertion table."""
    return pl.LazyFrame({
        'id': ['src-1'],
        'vartype': [vartype],
        'varsubtype': [varsubtype],
        'qry_id': ['q1'],
        'qry_rev': [False],
        'var_score': [1.0],
        'align_source': [[1, 2]],
        'filter': [[]],
        'dup': [dup],
    }, schema_overrides={
        'filter': pl.List(pl.String), 'dup': _DUP_TYPE,
        'align_source': pl.List(pl.Int32), 'varsubtype': pl.String,
    })


def _derive(df: pl.LazyFrame, ref_fa) -> pl.DataFrame:
    """Run the derivation and collect the duplication table."""
    collect_list = {'dup': []}
    add_ins_derived(df=df, collect_list=collect_list, ref_fa=ref_fa)

    if not collect_list['dup']:
        return pl.DataFrame()

    return pl.concat(collect_list['dup']).collect()


def test_tandem_insertion_derives_a_tandem_duplication(ref_fa):
    """A tandem-duplication insertion carries its subtype to the derived duplication."""
    out = _derive(_insdel_frame('INS', 'TD', [{'chrom': 'chr1', 'pos': 100, 'end': 140, 'is_rev': None}]), ref_fa)

    assert out.height == 1

    row = out.row(0, named=True)

    assert row['vartype'] == 'DUP'
    assert row['varsubtype'] == 'TD'
    assert (row['chrom'], row['pos'], row['end']) == ('chr1', 100, 140)
    assert row['varlen'] == 40
    assert row['derived'] == 'src-1'
    assert row['call_source'] == 'DERIVED'
    assert 'DERIVED' in row['filter']


def test_plain_insertion_derives_a_plain_duplication(ref_fa):
    """An insertion with no subtype derives a duplication with no subtype."""
    out = _derive(_insdel_frame('INS', None, [{'chrom': 'chr1', 'pos': 100, 'end': 140, 'is_rev': None}]), ref_fa)

    assert out.row(0, named=True)['varsubtype'] is None


def test_derived_duplication_carries_the_reference_sequence(ref_fa):
    """The sequence is the reference over the duplicated region."""
    out = _derive(_insdel_frame('INS', 'TD', [{'chrom': 'chr1', 'pos': 100, 'end': 140, 'is_rev': None}]), ref_fa)

    assert out.row(0, named=True)['seq'] == REF_SEQ[100:140]


def test_query_coordinates_are_left_null(ref_fa):
    """The duplicated region has no query interval to report."""
    row = _derive(_insdel_frame('INS', 'TD', [{'chrom': 'chr1', 'pos': 100, 'end': 140, 'is_rev': None}]),
                  ref_fa).row(0, named=True)

    assert row['qry_pos'] is None
    assert row['qry_end'] is None


def test_multiple_duplicated_regions_derive_one_each(ref_fa):
    """Each duplicated region becomes its own call."""
    out = _derive(_insdel_frame('INS', None, [
        {'chrom': 'chr1', 'pos': 100, 'end': 140, 'is_rev': None},
        {'chrom': 'chr1', 'pos': 300, 'end': 360, 'is_rev': None},
    ]), ref_fa)

    assert out.height == 2
    assert sorted(out['varlen'].to_list()) == [40, 60]


@pytest.mark.parametrize(
    ('vartype', 'dup'),
    [
        ('INS', None),
        ('INS', []),
        ('DEL', [{'chrom': 'chr1', 'pos': 100, 'end': 140, 'is_rev': None}]),
    ],
)
def test_nothing_derived_without_an_insertion_and_regions(ref_fa, vartype, dup):
    """Only insertions annotated with duplicated regions derive a duplication."""
    assert _derive(_insdel_frame(vartype, None, dup), ref_fa).height == 0


def test_derived_duplication_writes_no_info_seq():
    """The sequence stays in the variant tables but is kept out of the VCF."""
    df = pl.LazyFrame({
        'vartype': ['DUP', 'DUP'],
        'call_source': ['DERIVED', 'INTER'],
        'varsubtype': ['TD', 'TD'],
        'pos': [100, 100], 'end': [140, 140], 'varlen': [40, 40],
        'seq': ['acgt', 'acgt'],
    }, schema_overrides={'varlen': pl.Int64}).with_columns(pl.lit(True).alias('_is_sym_alt'))

    info = [';'.join(rec) for rec in standard_info_fields(df).collect()['_vcf_info'].to_list()]

    assert 'SEQ=' not in info[0]
    assert 'SEQ=ACGT' in info[1]
