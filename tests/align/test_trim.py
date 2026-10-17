"""Tests for alignment trimming preferences (pav3.align.trim).

Alignment tables are built from small synthetic SAM records with "=" CIGAR operations, so no sequence files are
needed. Every overlap is between perfectly matching alignments, so trimming either record removes the same number of
events and the outcome depends only on the trimming preferences:

* A filtered alignment is truncated in favor of a passing (PASS) alignment.
* For two alignments on the same chromosome and in the same orientation, the left-most record is trimmed so that the
  breakpoint is left-aligned in reference coordinates.
* Otherwise, the lower-scoring record is trimmed.
* If one record is contained within another, the filtered record is removed if only one is filtered, otherwise the
  lower-scoring record is removed.
"""

from pathlib import Path

import polars as pl
import pytest

from pav3.align.tables import sam_to_align_table
from pav3.align.trim import trim_alignments_qry, trim_alignments_ref

QRY_LEN = 2000

MIN_MAPQ = 20


def _align_table(tmp_path: Path, records: list[tuple], qry_len: int = QRY_LEN) -> pl.DataFrame:
    """Build an alignment table from synthetic SAM records.

    :param tmp_path: Directory for the SAM file.
    :param records: Tuples of (qry_id, flag, chrom, pos (1-based), mapq, cigar).
    :param qry_len: Length of every query sequence.

    :returns: Alignment table (trim-none) with records below MIN_MAPQ tagged with the MAPQ filter.
    """
    qry_ids = sorted({rec[0] for rec in records})

    sam_path = tmp_path / 'align.sam'
    with open(sam_path, 'wt') as out_file:
        out_file.write('@HD\tVN:1.6\n')
        out_file.write('@SQ\tSN:chr1\tLN:100000\n')
        for qry_id, flag, chrom, pos, mapq, cigar in records:
            out_file.write(f'{qry_id}\t{flag}\t{chrom}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t*\t*\n')

    df_qry_fai = pl.DataFrame({'qry_id': qry_ids, 'len': [qry_len] * len(qry_ids)})

    return sam_to_align_table(sam_path=sam_path, df_qry_fai=df_qry_fai, min_mapq=MIN_MAPQ)


def _qry_fai(df: pl.DataFrame, qry_len: int = QRY_LEN) -> pl.DataFrame:
    qry_ids = df['qry_id'].unique().sort().to_list()
    return pl.DataFrame({'qry_id': qry_ids, 'len': [qry_len] * len(qry_ids)})


def _by_index(df: pl.DataFrame) -> dict[int, dict]:
    return {row['align_index']: row for row in df.iter_rows(named=True)}


def test_min_mapq_tags_filter(tmp_path):
    """Records below the minimum MAPQ get the MAPQ filter, others pass."""
    df = _align_table(tmp_path, [
        ('q1', 0, 'chr1', 1001, 60, '1100=900S'),
        ('q1', 2048, 'chr1', 5001, 5, '1000H1000='),
    ])

    rows = _by_index(df)
    assert rows[0]['filter'] == []
    assert rows[1]['filter'] == ['MAPQ']


def test_qry_trim_truncates_filtered(tmp_path):
    """A partial query overlap between a PASS and a filtered alignment truncates the filtered alignment.

    Both orders are tested so the result does not depend on which record is on the left.
    """
    for mapq_a, mapq_b in [(60, 5), (5, 60)]:
        df = _align_table(tmp_path, [
            ('q1', 0, 'chr1', 1001, mapq_a, '1100=900S'),     # q[0, 1100) -> chr1:1000-2100
            ('q1', 2048, 'chr1', 5001, mapq_b, '1000H1000='),  # q[1000, 2000) -> chr1:5000-6000
        ])

        rows = _by_index(trim_alignments_qry(df, _qry_fai(df)))
        assert rows[0]['qry_end'] == rows[1]['qry_pos'] == (1000 if mapq_a < MIN_MAPQ else 1100), (mapq_a, mapq_b)

        if mapq_a < MIN_MAPQ:
            # Filtered record a truncated
            assert (rows[0]['pos'], rows[0]['end']) == (1000, 2000)
            assert (rows[1]['pos'], rows[1]['end']) == (5000, 6000)
        else:
            # Filtered record b truncated
            assert (rows[0]['pos'], rows[0]['end']) == (1000, 2100)
            assert (rows[1]['pos'], rows[1]['end']) == (5100, 6000)


def test_qry_trim_same_orientation_left_aligns(tmp_path):
    """A tie between same-orientation alignments trims the left-most record to left-align the breakpoint.

    The 100 query bases at q[1000, 1100) align to both chr1:2000-2100 and chr1:5000-5100. Keeping them in the right
    record places the deletion at chr1:2000-5000 (left-aligned). The left record has the higher score (1100 vs 1000
    aligned bases), so a score-based choice would keep them in the left record instead.
    """
    df = _align_table(tmp_path, [
        ('q1', 0, 'chr1', 1001, 60, '1100=900S'),     # q[0, 1100) -> chr1:1000-2100
        ('q1', 2048, 'chr1', 5001, 60, '1000H1000='),  # q[1000, 2000) -> chr1:5000-6000
    ])

    rows = _by_index(trim_alignments_qry(df, _qry_fai(df)))

    assert (rows[0]['pos'], rows[0]['end'], rows[0]['qry_end']) == (1000, 2000, 1000)
    assert (rows[1]['pos'], rows[1]['end'], rows[1]['qry_pos']) == (5000, 6000, 1000)


def test_qry_trim_opposite_orientation_uses_score(tmp_path):
    """A tie between opposite-orientation alignments trims the lower-scoring record.

    The right record is reverse-complemented and has the lower score (1000 vs 1100 aligned bases), so it is trimmed.
    """
    df = _align_table(tmp_path, [
        ('q1', 0, 'chr1', 1001, 60, '1100=900S'),     # q[0, 1100) -> chr1:1000-2100 (+)
        ('q1', 2064, 'chr1', 5001, 60, '1000=1000H'),  # q[1000, 2000) -> chr1:5000-6000 (-)
    ])

    rows = _by_index(df)
    assert not rows[0]['is_rev'] and rows[1]['is_rev']
    assert (rows[1]['qry_pos'], rows[1]['qry_end']) == (1000, 2000)

    rows = _by_index(trim_alignments_qry(df, _qry_fai(df)))

    assert (rows[0]['pos'], rows[0]['end'], rows[0]['qry_end']) == (1000, 2100, 1100)
    assert (rows[1]['qry_pos'], rows[1]['qry_end']) == (1100, 2000)
    assert rows[1]['end'] - rows[1]['pos'] == 900


@pytest.mark.parametrize('mapq_large,expected_kept', [(5, 0), (60, 1)])
def test_qry_trim_contained(tmp_path, mapq_large, expected_kept):
    """A record contained within another in query coordinates.

    The small record (q[100, 600), score of 500 aligned bases) is contained in the large record (q[0, 2000), score
    of 2000 aligned bases). When only the large record is filtered, it is removed. When both pass, the lower-scoring
    small record is removed.
    """
    df = _align_table(tmp_path, [
        ('q1', 2048, 'chr1', 1001, 60, '100H500=1400H'),  # q[100, 600) -> chr1:1000-1500
        ('q1', 0, 'chr1', 5001, mapq_large, '2000='),     # q[0, 2000) -> chr1:5000-7000
    ])

    df_trim = trim_alignments_qry(df, _qry_fai(df))

    assert df_trim['align_index'].to_list() == [expected_kept]


@pytest.mark.parametrize('mapq_large,expected_kept', [(5, 0), (60, 1)])
def test_ref_trim_contained(tmp_path, mapq_large, expected_kept):
    """A record contained within another in reference coordinates (two different query sequences).

    The small record (chr1:1100-1600, score of 500 aligned bases) is contained in the large record (chr1:1000-3000,
    score of 2000 aligned bases).
    """
    df = _align_table(tmp_path, [
        ('q1', 0, 'chr1', 1101, 60, '500=1500S'),     # q1[0, 500) -> chr1:1100-1600
        ('q2', 0, 'chr1', 1001, mapq_large, '2000='),  # q2[0, 2000) -> chr1:1000-3000
    ])

    df_trim = trim_alignments_ref(df, _qry_fai(df))

    assert df_trim['align_index'].to_list() == [expected_kept]
