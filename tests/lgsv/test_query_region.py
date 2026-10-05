"""Tests for the query regions reported by large SV calls."""
from __future__ import annotations

from types import SimpleNamespace

import polars as pl
import pytest

from pav3.anno import perfect_homology
from pav3.lgsv.variant import ComplexVariant, DeletionVariant, InsertionVariant
from pav3.region import Region


REF_SEQ = 'ACGTTGCA' * 4000
QRY_SEQ = 'TTGACCAG' * 4000


class _ScoreModel:
    """Linear score model, enough to exercise the call paths."""

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
        df_ref_fai=pl.DataFrame({'chrom': ['chr1'], 'len': [len(REF_SEQ)]}),
    )


def _interval(len_ins: int, len_del: int) -> SimpleNamespace:
    """Build an insertion/deletion interval.

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


#
# Deletion query region
#

@pytest.mark.parametrize('len_ins', [0, 32])
def test_deletion_reports_a_one_base_query_region(caller_resources, len_ins):
    """A deletion reports a 1 bp query region whether or not bases were inserted with it."""
    variant = DeletionVariant(_interval(len_ins, 6494), caller_resources, SimpleNamespace())

    assert variant.qry_pos == 500
    assert variant.qry_end == 501
    assert len(variant.region_qry) == 1


def test_off_target_deletion_reports_a_one_base_query_region(caller_resources):
    """A deletion emitted as an insertion's off-target reports the same 1 bp query region."""
    variant = InsertionVariant(_interval(500, 32), caller_resources, SimpleNamespace())
    variant.complete_anno()

    assert variant.off_target_del.qry_end - variant.off_target_del.qry_pos == 1


def test_deletion_homology_uses_the_region_between_the_anchors(caller_resources):
    """Widening the reported query region does not shift the flanks homology is measured against."""
    interval = _interval(0, 6494)
    variant = DeletionVariant(interval, caller_resources, SimpleNamespace())

    assert variant.region_qry_seq == interval.region_qry

    hom = perfect_homology(
        seq_var=variant.seq,
        seq_ref=QRY_SEQ,
        pos=interval.region_qry.pos,
        end=interval.region_qry.end,
        is_rev=False,
    )

    assert variant.hom_qry == {'up': hom[0], 'dn': hom[1]}


def test_homology_region_defaults_to_the_query_region(caller_resources):
    """A call that does not set a homology region measures over the region it reports."""
    variant = InsertionVariant(_interval(500, 0), caller_resources, SimpleNamespace())

    hom = perfect_homology(
        seq_var=variant.seq,
        seq_ref=QRY_SEQ,
        pos=variant.qry_pos,
        end=variant.qry_end,
        is_rev=False,
    )

    assert variant.region_qry_seq is None
    assert variant.hom_qry == {'up': hom[0], 'dn': hom[1]}


#
# Complex query region
#

def _interval_cpx(len_qry: int, len_dup: int) -> SimpleNamespace:
    """Build a complex interval with an outer tandem duplication.

    :param len_qry: Query bases between the untrimmed anchors.
    :param len_dup: Reference bases covered by both anchors (outer duplication).
    """
    df_segment = pl.DataFrame({
        'seg_index': [0, 1, 2],
        'is_anchor': [True, False, True],
        'is_aligned': [True] * 3,
        'filter': [[], [], []],
        'filter_pass': [True] * 3,
        'chrom': ['chr1'] * 3,
        'pos': [1000, 4000, 2400 - len_dup],
        'end': [2400, 4200, 3800],
        'is_rev': [False] * 3,
        'qry_id': ['q1'] * 3,
        'qry_pos': [0, 2000, 2000 + len_qry],
        'qry_end': [2000, 2200, 2400 + len_qry],
        'align_index': [0, 1, 2],
    }, schema_overrides={'filter': pl.List(pl.String)})

    return SimpleNamespace(
        start_index=0, end_index=2, df_segment=df_segment, is_rev=False, chrom='chr1',
        n_anchor=2, all_anchor_pass=True,
        seg_len_aligned=[200], seg_len_unaligned=[len_qry - 200],
        seg_len_unaligned_or_outer=[len_qry - 200, len_dup], seg_n_aligned=1,
        len_outer_dup=len_dup, len_outer_del=0, len_qry=len_qry, len_ref=-len_dup,
        region_ref=Region('chr1', 2400 - len_dup, 2400),
        region_qry=Region('q1', 2000, 2000 + len_qry),
        region_ref_qryref=Region('chr1', 2400 - len_dup, 2400),
        region_qry_qryref=Region('q1', 2000 - len_dup, 2000 + len_qry),
        qryref_spans_outer_dup=True,
        qry_aligned_pass_prop=1.0,
    )


#
# Tandem duplications reference trimming did not resolve
#

def _interval_td(len_dup: int, len_ins: int = 0, resolved: bool = True) -> SimpleNamespace:
    """Build a tandem-duplication interval.

    :param len_dup: Reference bases covered by both anchors (outer duplication).
    :param len_ins: Query bases between the anchors beyond the duplication.
    :param resolved: Whether reference trimming removed the overlap between the anchors, exposing the duplicated
        bases in the query region between them.
    """
    df_segment = pl.DataFrame({
        'is_anchor': [True, True],
        'filter': [[], []],
        'chrom': ['chr1'] * 2,
        'pos': [1000, 2400 - len_dup],
        'end': [2400, 3800],
        'is_aligned': [True] * 2,
        'is_rev': [False] * 2,
        'qry_id': ['q1'] * 2,
        'qry_pos': [0, 2000 + len_ins],
        'qry_end': [2000, 3400 + len_ins],
        'align_index': [0, 1],
    }, schema_overrides={'filter': pl.List(pl.String)})

    region_ref = Region('chr1', 2400 - len_dup, 2400)
    region_qry = Region('q1', 2000, 2000 + len_ins)

    return SimpleNamespace(
        start_index=0, end_index=1, df_segment=df_segment, is_rev=False, chrom='chr1',
        seg_len_aligned=[], seg_len_unaligned=[len_ins] if len_ins else [],
        seg_len_or_dup=[len_dup], seg_len=[len_ins] if len_ins else [],
        len_outer_dup=len_dup, len_outer_del=0, len_qry=len_ins, len_ref=-len_dup,
        region_ref=region_ref, region_qry=region_qry, region_ref_qryref=region_ref,
        region_qry_qryref=(
            Region('q1', 2000 - len_dup, 2000 + len_ins) if resolved else region_qry
        ),
        qryref_spans_outer_dup=resolved,
        qry_aligned_pass_prop=1.0,
    )


@pytest.mark.parametrize('resolved', [True, False])
def test_tandem_duplication_sequence_is_complete(caller_resources, resolved):
    """A tandem duplication carries its full sequence whether or not trimming exposed the duplicated bases."""
    variant = InsertionVariant(_interval_td(len_dup=1400, resolved=resolved), caller_resources, SimpleNamespace())

    assert variant.varsubtype == 'TD'
    assert variant.varlen == 1400
    assert len(variant.seq) == variant.varlen


def test_unresolved_tandem_duplication_reports_a_one_base_query_region(caller_resources):
    """A duplication whose bases appear in no query region reports 1 bp rather than an empty range."""
    variant = InsertionVariant(_interval_td(len_dup=1400, resolved=False), caller_resources, SimpleNamespace())

    assert variant.qry_end - variant.qry_pos == 1
    assert variant.region_qry_seq.pos == variant.region_qry_seq.end


def test_unresolved_tandem_duplication_takes_bases_from_the_reference(caller_resources):
    """The duplicated bases come from the reference when no query region covers them."""
    interval = _interval_td(len_dup=1400, resolved=False)
    variant = InsertionVariant(interval, caller_resources, SimpleNamespace())

    assert variant.seq == REF_SEQ[interval.region_ref.pos:interval.region_ref.end]


def test_resolved_tandem_duplication_takes_bases_from_the_query(caller_resources):
    """The duplicated bases come from the query when reference trimming exposed them."""
    interval = _interval_td(len_dup=1400, resolved=True)
    variant = InsertionVariant(interval, caller_resources, SimpleNamespace())

    assert variant.region_qry == interval.region_qry_qryref
    assert variant.seq == QRY_SEQ[interval.region_qry_qryref.pos:interval.region_qry_qryref.end]


def test_unresolved_tandem_duplication_with_inserted_bases(caller_resources):
    """Bases inserted alongside an unresolved duplication are kept with it."""
    variant = InsertionVariant(_interval_td(len_dup=1400, len_ins=60, resolved=False), caller_resources,
                               SimpleNamespace())

    assert variant.varlen == 1460
    assert len(variant.seq) == variant.varlen
    assert variant.qry_end - variant.qry_pos == 60


def test_complex_takes_the_query_region_from_trimmed_alignments(caller_resources):
    """A complex call reports the query region that spans duplicated bases at the variant site."""
    interval = _interval_cpx(len_qry=1010, len_dup=602)
    variant = ComplexVariant(interval, caller_resources, SimpleNamespace())

    assert variant.is_variant
    assert variant.region_qry == interval.region_qry_qryref
    assert variant.qry_spans_outer_dup


def test_complex_length_matches_the_query_region(caller_resources):
    """A complex call length is the length of the query region it reports."""
    variant = ComplexVariant(_interval_cpx(len_qry=1010, len_dup=602), caller_resources, SimpleNamespace())

    assert variant.varlen == 1612
    assert variant.varlen == variant.qry_end - variant.qry_pos


def test_complex_sequence_counts_duplicated_bases_once(caller_resources):
    """Duplicated bases already covered by the query region are not added a second time."""
    variant = ComplexVariant(_interval_cpx(len_qry=1010, len_dup=602), caller_resources, SimpleNamespace())

    assert len(variant.seq) == variant.varlen
