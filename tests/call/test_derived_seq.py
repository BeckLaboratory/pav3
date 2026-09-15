"""Tests for variant sequences that produced structurally invalid VCF records.

Covers two defects and the guard that catches them:

* CPX-derived deletions carried no "seq", so the VCF REF field was empty.
* Insertions at a tandem duplication carried an empty "seq", so REF matched ALT.
"""
from __future__ import annotations

from types import SimpleNamespace

import polars as pl
import pysam
import pytest

from pav3.call.integrate import add_cpx_derived
from pav3.lgsv.variant import Variant
from pav3.vcf import vcf_fields_seq_insdel


REF_SEQ = 'ACGT' * 250


@pytest.fixture
def ref_fa(tmp_path):
    """Write a small indexed reference FASTA."""
    path = tmp_path / 'ref.fa'
    path.write_text(f'>chr1\n{REF_SEQ}\n')
    pysam.faidx(str(path))

    return path


# --- Variant.seq: tandem duplication ---

def _seq_stub(len_ref: int, qry_seq: str = '', dup_pos: int = 100, dup_end: int = 140) -> SimpleNamespace:
    """Build the attributes "Variant.seq" reads for a non-DEL variant.

    :param len_ref: Reference distance between anchors. Negative for a tandem duplication.
    :param qry_seq: Query sequence between the anchors.
    :param dup_pos: Start of the duplicated reference region.
    :param dup_end: End of the duplicated reference region.
    """
    return SimpleNamespace(
        is_null=False,
        is_patch=False,
        _seq=None,
        vartype='INS',
        qry_rev=False,
        caller_resources=SimpleNamespace(cache_ref={'chr1': REF_SEQ}, cache_qry={'qry1': qry_seq}),
        region_qry=SimpleNamespace(chrom='qry1', pos=0, end=len(qry_seq)),
        region_ref_qrytrim=SimpleNamespace(chrom='chr1', pos=dup_pos, end=dup_end),
        interval=SimpleNamespace(len_ref=len_ref),
    )


def test_tandem_dup_insertion_seq_is_not_empty():
    """A pure tandem duplication has no query bases between anchors but must still carry sequence."""
    seq = Variant.seq.fget(_seq_stub(len_ref=-40))

    assert seq == REF_SEQ[100:140]
    assert len(seq) == 40


def test_tandem_dup_insertion_seq_includes_query_bases():
    """A duplication with additional inserted bases carries both."""
    seq = Variant.seq.fget(_seq_stub(len_ref=-40, qry_seq='TTTT'))

    assert seq == 'TTTT' + REF_SEQ[100:140]


def test_insertion_without_duplication_is_unchanged():
    """An insertion with no reference overlap reads only the query."""
    assert Variant.seq.fget(_seq_stub(len_ref=0, qry_seq='TTTT')) == 'TTTT'


# --- add_cpx_derived: derived deletions ---

def test_add_cpx_derived_fills_seq(ref_fa):
    """Derived variants at the outer breakpoint carry the reference sequence over their span."""
    df = pl.LazyFrame({
        'var_index': [0],
        'chrom': ['chr1'],
        'qry_id': ['qry1'],
        'qry_rev': [False],
        'var_score': [1.0],
        'qry_pos': [0],
        'qry_end': [0],
        'id': ['chr1-101-CPX-50'],
        'filter': [[]],
    }, schema_overrides={'filter': pl.List(pl.String), 'qry_pos': pl.Int64, 'qry_end': pl.Int64})

    # Two anchors leaving a 40 bp reference gap -> a derived deletion of chr1:100-140.
    df_segment = pl.LazyFrame({
        'var_index': [0, 0],
        'chrom': ['chr1', 'chr1'],
        'pos': [50, 140],
        'end': [100, 190],
        'is_anchor': [True, True],
        'is_aligned': [True, True],
        'is_rev': [False, False],
        'qry_id': ['qry1', 'qry1'],
        'qry_pos': [50, 100],
        'qry_end': [100, 150],
        'align_index': [0, 1],
    })

    collect_list = {'insdel': [], 'dup': [], 'inv': []}
    add_cpx_derived(df=df, df_segment=df_segment, collect_list=collect_list, ref_fa=ref_fa)

    df_insdel = pl.concat(collect_list['insdel'], how='diagonal').collect()

    assert df_insdel.height == 1

    row = df_insdel.row(0, named=True)

    assert row['vartype'] == 'DEL'
    assert row['varlen'] == 40
    assert row['seq'] == REF_SEQ[100:140]


# --- vcf_fields_seq_insdel: the guard ---

def _insdel_frame(vartype: str, seq: str | None) -> pl.LazyFrame:
    """Build a one-row insdel table."""
    return pl.LazyFrame(
        {'id': ['v1'], 'chrom': ['chr1'], 'pos': [20], 'end': [24], 'vartype': [vartype], 'seq': [seq]},
        schema_overrides={'seq': pl.String},
    ).with_row_index('_index')


def test_valid_insdel_builds_ref_and_alt(ref_fa):
    """A well-formed insertion and deletion gain the anchor base."""
    out = vcf_fields_seq_insdel(_insdel_frame('DEL', REF_SEQ[20:24]), ref_fa).collect()

    assert out.row(0, named=True)['_vcf_ref'] == REF_SEQ[19:24]
    assert out.row(0, named=True)['_vcf_alt'] == REF_SEQ[19]


@pytest.mark.parametrize(
    ('vartype', 'seq'),
    [
        ('DEL', None),  # Derived deletion with no sequence -> empty REF
        ('DEL', ''),
        ('INS', None),
        ('INS', ''),  # Tandem-duplication insertion with no sequence -> REF matches ALT
    ],
)
def test_missing_seq_is_rejected(ref_fa, vartype, seq):
    """A missing or empty sequence must fail rather than write an invalid record."""
    with pytest.raises(ValueError, match='Invalid REF or ALT'):
        vcf_fields_seq_insdel(_insdel_frame(vartype, seq), ref_fa).collect()
