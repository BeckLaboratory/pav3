"""Tests for paired insertion and deletion calls over one alignment interval."""
from __future__ import annotations

from types import SimpleNamespace

import polars as pl
import pytest

from pav3.lgsv.variant import DeletionVariant, InsertionVariant
from pav3.region import Region


REF_SEQ = 'ACGT' * 5000
QRY_SEQ = 'TGCA' * 5000


class _ScoreModel:
    """Linear score model, enough to exercise penalty arithmetic."""

    @staticmethod
    def gap(length: int) -> float:
        """Gap penalty."""
        return -float(abs(length))

    @staticmethod
    def off_variant() -> float:
        """Off-variant weight."""
        return 0.5

    @staticmethod
    def template_switch(n: int = 1) -> float:
        """Template switch penalty."""
        return -100.0


@pytest.fixture
def caller_resources():
    """Caller resources with in-memory sequence caches."""
    return SimpleNamespace(
        score_model=_ScoreModel(),
        pav_params=SimpleNamespace(lg_cpx_min_aligned_prop=0.5),
        verbose=False,
        cache_ref={'chr1': REF_SEQ},
        cache_qry={'q1': QRY_SEQ},
    )


def _interval(len_ins: int, len_del: int) -> SimpleNamespace:
    """Build an interval with inserted query bases and deleted reference bases.

    :param len_ins: Inserted query bases between the anchors.
    :param len_del: Deleted reference bases between the anchors.
    """
    df_segment = pl.DataFrame({
        'is_anchor': [True, True],
        'filter': [[], []],
        'chrom': ['chr1'] * 2,
        'pos': [0, 0],
        'end': [1, 1],
        'is_aligned': [True] * 2,
        'is_rev': [False] * 2,
        'qry_id': ['q1'] * 2,
        'qry_pos': [0, 0],
        'qry_end': [1, 1],
        'align_index': [0, 1],
    }, schema_overrides={'filter': pl.List(pl.String)})

    region_ref = Region('chr1', 1000, 1000 + max(len_del, 1))
    region_qry = Region('q1', 500, 500 + len_ins)

    return SimpleNamespace(
        start_index=0, end_index=1, df_segment=df_segment, is_rev=False, chrom='chr1',
        seg_len_aligned=[], seg_len_unaligned=[len_ins] if len_ins else [],
        seg_len_or_dup=[], seg_len=[len_ins] if len_ins else [],
        len_outer_dup=0, len_outer_del=len_del, len_qry=len_ins, len_ref=len_del,
        region_ref=region_ref, region_qry=region_qry,
        region_ref_qryref=region_ref, region_qry_qryref=region_qry, qryref_spans_outer_dup=True,
        qry_aligned_pass_prop=1.0,
    )


def _call(cls, caller_resources, len_ins: int, len_del: int):
    """Build a completed variant call."""
    variant = cls(_interval(len_ins, len_del), caller_resources, SimpleNamespace())
    variant.complete_anno()

    return variant


def test_insertion_emits_off_target_deletion(caller_resources):
    """An insertion over deleted reference bases emits the deletion as its own record."""
    variant = _call(InsertionVariant, caller_resources, len_ins=500, len_del=32)

    assert variant.off_target_del is not None
    assert variant.off_target_del.vartype == 'DEL'
    assert variant.off_target_del.varlen == 32


def test_deletion_emits_off_target_insertion(caller_resources):
    """A deletion with inserted bases emits the insertion as its own record."""
    variant = _call(DeletionVariant, caller_resources, len_ins=32, len_del=6494)

    assert variant.off_target_ins is not None
    assert variant.off_target_ins.vartype == 'INS'
    assert variant.off_target_ins.varlen == 32


def test_insertion_varlen_is_not_netted(caller_resources):
    """The insertion length does not subtract the deleted bases, which are their own record."""
    variant = _call(InsertionVariant, caller_resources, len_ins=500, len_del=32)

    assert variant.varlen == 500
    assert variant.varlen == len(variant.seq)


def test_var_score_keeps_the_off_target_penalty(caller_resources):
    """The score used to resolve the optimal path is never modified."""
    variant = _call(InsertionVariant, caller_resources, len_ins=500, len_del=32)

    # gap(500) + gap(32) * off_variant()
    assert variant.var_score == pytest.approx(-516.0)


def test_emitted_score_drops_the_off_target_penalty(caller_resources):
    """The emitted record carries the score without the penalty, leaving var_score intact."""
    variant = _call(InsertionVariant, caller_resources, len_ins=500, len_del=32)

    assert variant.var_score_primary == pytest.approx(-500.0)
    assert variant.rows()[0]['var_score'] == pytest.approx(-500.0)
    assert variant.var_score == pytest.approx(-516.0)


def test_off_target_score_has_no_penalty(caller_resources):
    """An off-target call is scored on its own length alone."""
    variant = _call(InsertionVariant, caller_resources, len_ins=500, len_del=32)

    assert variant.off_target_del.var_score == pytest.approx(-32.0)
    assert variant.off_target_del.var_score_primary is None


@pytest.mark.parametrize(
    ('cls', 'len_ins', 'len_del'),
    [
        (InsertionVariant, 500, 0),
        (DeletionVariant, 0, 6494),
    ],
)
def test_no_off_target_emits_one_record(caller_resources, cls, len_ins, len_del):
    """An interval with a single component emits one record and sets no primary score."""
    variant = _call(cls, caller_resources, len_ins, len_del)

    assert len(variant.rows()) == 1
    assert variant.var_score_primary is None
    assert variant.off_target_del is None
    assert variant.off_target_ins is None


def test_rows_emits_both_records(caller_resources):
    """Both components appear, each carrying its own length and sequence."""
    rows = _call(InsertionVariant, caller_resources, len_ins=500, len_del=32).rows()

    assert [row['vartype'] for row in rows] == ['INS', 'DEL']
    assert [row['varlen'] for row in rows] == [500, 32]
    assert [len(row['seq']) for row in rows] == [500, 32]


def test_off_target_does_not_recurse(caller_resources):
    """A variant built as an off-target never looks for one of its own."""
    variant = _call(InsertionVariant, caller_resources, len_ins=500, len_del=32)
    off_target = variant.off_target_del

    off_target.complete_anno()

    assert off_target.off_target_ins is None
    assert off_target.off_target_del is None
