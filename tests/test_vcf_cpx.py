"""Tests for sequence-resolved complex variants and symbolic REF bases in :mod:`pav3.vcf`.

A sequence-resolved complex variant must reconstruct the haplotype when applied to the reference: the haplotype is
the reference up to the end of the left anchor, the query between the anchors, then the reference from the start of
the right anchor.
"""
from __future__ import annotations

import random
from types import SimpleNamespace

import polars as pl
import pysam
import pytest

from pav3.lgsv.variant import Variant, seq_ref_region
from pav3.vcf import init_vcf_fields, vcf_fields_seq_cpx, vcf_fields_sym

REF_SEQ = ''.join(random.Random(1).choice('ACGT') for _ in range(1000))

QRY_SEQ = 'TTTTGGGGCCCCAAAA'


@pytest.fixture
def ref_fa(tmp_path):
    """Write a small indexed reference FASTA."""
    path = tmp_path / 'ref.fa'
    path.write_text(f'>chr1\n{REF_SEQ}\n')
    pysam.faidx(str(path))

    return path


def _apply_record(pos: int, ref: str, alt: str) -> str:
    """Apply a VCF record (1-based POS) to the reference."""
    assert REF_SEQ[pos - 1:pos - 1 + len(ref)].upper() == ref.upper(), 'REF does not match the reference'
    return REF_SEQ[:pos - 1] + alt + REF_SEQ[pos - 1 + len(ref):]


def _cpx_frame(seq_ref_pos: int, seq_ref_end: int, seq: str) -> pl.LazyFrame:
    return pl.LazyFrame({
        'chrom': ['chr1'],
        'pos': [seq_ref_pos],
        'end': [max(seq_ref_end, seq_ref_pos + 1)],
        'id': ['cpx1'],
        'vartype': ['CPX'],
        'filter': [[]],
        'seq_ref_pos': [seq_ref_pos],
        'seq_ref_end': [seq_ref_end],
        'seq': [seq],
    }, schema_overrides={'filter': pl.List(pl.String)}).with_row_index('_index')


@pytest.mark.parametrize(
    ('seq_ref_pos', 'seq_ref_end'),
    [
        (100, 150),  # Deleted reference bases replaced by the sequence
        (100, 100),  # Insertion
        (0, 40),     # Starts at the first base: padding base on the right
        (0, 0),      # Insertion before the first base
    ],
)
def test_cpx_record_reconstructs_haplotype(ref_fa, seq_ref_pos, seq_ref_end):
    """Applying the record replaces [seq_ref_pos, seq_ref_end) with the variant sequence."""
    row = vcf_fields_seq_cpx(_cpx_frame(seq_ref_pos, seq_ref_end, QRY_SEQ), ref_fa).collect().row(0, named=True)

    assert not row['_is_sym_alt']
    assert _apply_record(row['_vcf_pos'], row['_vcf_ref'], row['_vcf_alt']) == (
        REF_SEQ[:seq_ref_pos] + QRY_SEQ + REF_SEQ[seq_ref_end:]
    )


def test_cpx_sequence_resolved_by_default(ref_fa):
    """Complex variants are sequence-resolved unless symbolic records are requested."""
    df = _cpx_frame(100, 150, QRY_SEQ)

    row = init_vcf_fields(df, ref_fa=ref_fa, vartype='cpx').collect().row(0, named=True)
    assert row['_vcf_alt'] == REF_SEQ[99] + QRY_SEQ
    assert not row['_is_sym_alt']

    row = init_vcf_fields(df, ref_fa=ref_fa, vartype='cpx', cpx_symbolic=True).collect().row(0, named=True)
    assert row['_vcf_alt'] == '<CPX>'
    assert row['_is_sym_alt']


def test_cpx_missing_columns_raises(ref_fa):
    """A complex variant table without the replaced region (earlier PAV versions) is rejected."""
    with pytest.raises(ValueError, match='seq_ref_pos'):
        vcf_fields_seq_cpx(_cpx_frame(100, 150, QRY_SEQ).drop('seq_ref_pos', 'seq_ref_end'), ref_fa)


def _interval(len_ref, region_ref, region_ref_qryref=None, len_ref_qryref=None):
    return SimpleNamespace(
        len_ref=len_ref,
        region_ref=SimpleNamespace(chrom='chr1', pos=region_ref[0], end=region_ref[1]),
        region_ref_qryref=(
            SimpleNamespace(chrom='chr1', pos=region_ref_qryref[0], end=region_ref_qryref[1])
            if region_ref_qryref is not None else None
        ),
        len_ref_qryref=len_ref_qryref,
    )


def _variant_seq(len_ref, region_qry_seq, dup_region, qry_spans_outer_dup) -> str:
    """Variant sequence from "Variant.seq" (query between anchors, with reference duplication when needed)."""
    return Variant.seq.fget(SimpleNamespace(
        is_null=False,
        is_patch=False,
        _seq=None,
        vartype='CPX',
        qry_rev=False,
        qry_spans_outer_dup=qry_spans_outer_dup,
        caller_resources=SimpleNamespace(cache_ref={'chr1': REF_SEQ}, cache_qry={'qry1': region_qry_seq}),
        region_qry=SimpleNamespace(chrom='qry1', pos=0, end=len(region_qry_seq)),
        region_qry_seq=None,
        region_ref_qrytrim=SimpleNamespace(chrom='chr1', pos=dup_region[0], end=dup_region[1]),
        interval=SimpleNamespace(len_ref=len_ref),
    ))


def test_seq_ref_region_anchor_gap():
    """Anchors leave a reference gap (left anchor ends at 100, right anchor starts at 150)."""
    interval = _interval(len_ref=50, region_ref=(100, 150))
    seq = _variant_seq(50, QRY_SEQ, (100, 150), False)

    pos, end = seq_ref_region(interval, False)

    assert (pos, end) == (100, 150)
    assert REF_SEQ[:pos] + seq + REF_SEQ[end:] == REF_SEQ[:100] + QRY_SEQ + REF_SEQ[150:]


def test_seq_ref_region_insertion_site():
    """Anchors meet on the reference (insertion site, region_ref end is pos + 1)."""
    interval = _interval(len_ref=0, region_ref=(100, 101))
    seq = _variant_seq(0, QRY_SEQ, (100, 101), False)

    pos, end = seq_ref_region(interval, False)

    assert (pos, end) == (100, 100)
    assert REF_SEQ[:pos] + seq + REF_SEQ[end:] == REF_SEQ[:100] + QRY_SEQ + REF_SEQ[100:]


def test_seq_ref_region_overlap_reference_dup():
    """Anchors overlap (left anchor ends at 140, right starts at 100) and the duplication comes from the reference.

    The haplotype is ref[:140] + Q + ref[100:]: the reference over [100, 140) appears twice.
    """
    interval = _interval(len_ref=-40, region_ref=(100, 140))
    seq = _variant_seq(-40, QRY_SEQ, (100, 140), False)

    pos, end = seq_ref_region(interval, False)

    assert REF_SEQ[:pos] + seq + REF_SEQ[end:] == REF_SEQ[:140] + QRY_SEQ + REF_SEQ[100:]


def test_seq_ref_region_overlap_qryref():
    """Anchors overlap, and reference trimming cut them to meet at 120.

    The query between the trimmed anchors holds the left anchor's bases over [120, 140), the inserted sequence, and
    the right anchor's bases over [100, 120).
    """
    qryref_seq = REF_SEQ[120:140] + QRY_SEQ + REF_SEQ[100:120]

    interval = _interval(len_ref=-40, region_ref=(100, 140), region_ref_qryref=(120, 121), len_ref_qryref=0)
    seq = _variant_seq(-40, qryref_seq, (100, 140), True)

    pos, end = seq_ref_region(interval, True)

    assert (pos, end) == (120, 120)
    assert REF_SEQ[:pos] + seq + REF_SEQ[end:] == REF_SEQ[:140] + QRY_SEQ + REF_SEQ[100:]


@pytest.mark.parametrize('pos', [10, 0])
def test_symbolic_ref_is_padding_base(ref_fa, pos):
    """Symbolic records use the base before the variant as REF, at POS (the first base at the chromosome start)."""
    df = pl.LazyFrame({
        'chrom': ['chr1'],
        'vartype': ['INV'],
        'varsubtype': [None],
        'pos': [pos],
        'end': [pos + 20],
    }, schema_overrides={'varsubtype': pl.String}).with_row_index('_index')

    row = vcf_fields_sym(df, ref_fa).collect().row(0, named=True)

    assert row['_vcf_pos'] == max(pos, 1)
    assert row['_vcf_ref'] == REF_SEQ[max(pos, 1) - 1]


def test_cpx_matching_reference_is_dropped(ref_fa, caplog):
    """A complex variant whose sequence matches the bases it replaces (in any case) is left out with a warning."""
    df = _cpx_frame(100, 150, REF_SEQ[100:150].lower())

    with caplog.at_level('WARNING', logger='pav3.vcf'):
        assert vcf_fields_seq_cpx(df, ref_fa).collect().height == 0
        stats = {}
        assert init_vcf_fields(df, ref_fa=ref_fa, vartype='cpx', log_stats=stats).collect().height == 0

    assert 'cpx1' in caplog.text
    assert stats == {'n_noop_cpx': 1, 'noop_cpx_label': ['cpx1']}


def test_cpx_without_sequence_matching_reference_raises(ref_fa):
    """REF matching ALT without a sequence is a callset bug, not a variant to skip."""
    with pytest.raises(ValueError, match='Invalid REF or ALT'):
        vcf_fields_seq_cpx(_cpx_frame(100, 100, ''), ref_fa).collect()
