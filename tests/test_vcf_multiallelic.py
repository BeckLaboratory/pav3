"""Tests for combining records at the same site into multi-allelic records (:func:`pav3.vcf.combine_multiallelic`)."""
from __future__ import annotations

import polars as pl

from pav3.vcf import combine_multiallelic, reformat_vcf_table, standard_info_fields, vcf_record_table

SAMPLE_COL = '_vcf_sample_S'


def _rec(pos, ref, alt, gt, rec_id, info=None, filter_='PASS', sym=False, chrom='chr1'):
    return {
        '_vcf_chrom': chrom,
        '_vcf_pos': pos,
        '_vcf_id': rec_id,
        '_vcf_ref': ref,
        '_vcf_alt': alt,
        '_vcf_filter': filter_,
        '_vcf_info': info if info is not None else [f'SVLEN={abs(len(alt) - len(ref)) or 1}', f'QUERY=q:{rec_id}'],
        '_vcf_format': ['GT'],
        '_is_sym_alt': sym,
        SAMPLE_COL: [gt],
    }


def _combine(records):
    stats = {}
    df = combine_multiallelic(pl.DataFrame(records), log_stats=stats)
    return reformat_vcf_table(df).collect().rows(named=True), stats


def test_het_alt_deletions():
    """Two deletions of different lengths at one position (DEL12 on h1, DEL10 on h2) become one "1|2" record."""
    rows, stats = _combine([
        _rec(100, 'A' + 'C' * 12, 'A', '1|0', 'del12'),
        _rec(100, 'A' + 'C' * 10, 'A', '0|1', 'del10'),
    ])

    assert len(rows) == 1
    row = rows[0]

    assert row['REF'] == 'A' + 'C' * 12
    assert row['ALT'] == 'A,ACC'
    assert row['S'] == '1|2'
    assert row['ID'] == 'del12;del10'
    assert row['INFO'] == 'SVLEN=12,10;QUERY=q:del12,q:del10'
    assert stats['n_combined_site'] == 1
    assert stats['n_combined_record'] == 2


def test_het_alt_insertions():
    """Insertions of 8 T and 7 T at one position."""
    rows, _ = _combine([
        _rec(200, 'G', 'G' + 'T' * 7, '0|1', 'ins7'),
        _rec(200, 'G', 'G' + 'T' * 8, '1|0', 'ins8'),
    ])

    assert len(rows) == 1
    assert rows[0]['ALT'] == 'G' + 'T' * 8 + ',G' + 'T' * 7
    assert rows[0]['S'] == '1|2'
    assert rows[0]['ID'] == 'ins8;ins7'


def test_snv_different_alt():
    rows, _ = _combine([
        _rec(300, 'A', 'C', '1|0', 'snv_c'),
        _rec(300, 'A', 'G', '0|1', 'snv_g'),
    ])

    assert len(rows) == 1
    assert (rows[0]['REF'], rows[0]['ALT'], rows[0]['S']) == ('A', 'C,G', '1|2')


def test_snv_on_insertion_pad_base():
    """An SNV on h1 and an insertion after the same base on h2 share POS and combine."""
    rows, _ = _combine([
        _rec(400, 'A', 'C', '1|0', 'snv'),
        _rec(400, 'A', 'ATTT', '0|1', 'ins'),
    ])

    assert len(rows) == 1
    assert (rows[0]['REF'], rows[0]['ALT'], rows[0]['S']) == ('A', 'C,ATTT', '1|2')


def test_snv_and_deletion_extends_alt():
    """A shorter REF is extended with the reference bases of the longer REF."""
    rows, _ = _combine([
        _rec(500, 'A', 'G', '1|0', 'snv'),
        _rec(500, 'ACGT', 'A', '0|1', 'del'),
    ])

    assert len(rows) == 1
    assert (rows[0]['REF'], rows[0]['ALT'], rows[0]['S']) == ('ACGT', 'GCGT,A', '1|2')


def test_same_haplotype_not_combined():
    """A haplotype carrying two records at one position (SNV plus insertion) leaves the records as they are."""
    rows, stats = _combine([
        _rec(600, 'A', 'C', '1|0', 'snv'),
        _rec(600, 'A', 'ATTT', '1|1', 'ins'),
    ])

    assert len(rows) == 2
    assert stats['n_conflict_site'] == 1
    assert stats['n_combined_site'] == 0


def test_different_filter_not_combined():
    """A passing allele is not combined with a filtered allele."""
    rows, stats = _combine([
        _rec(700, 'A', 'C', '1|0', 'snv_pass'),
        _rec(700, 'A', 'G', '0|1', 'snv_lc', filter_='LCALIGN'),
    ])

    assert len(rows) == 2
    assert {row['FILTER'] for row in rows} == {'PASS', 'LCALIGN'}
    assert stats['n_combined_site'] == 0


def test_overlap_at_different_position_not_combined():
    """Overlapping records starting at different positions are left as they are and counted."""
    rows, stats = _combine([
        _rec(800, 'A' + 'C' * 10, 'A', '1|0', 'del'),
        _rec(805, 'C', 'T', '0|1', 'snv'),
    ])

    assert len(rows) == 2
    assert stats['n_overlap_record'] == 1


def test_uncallable_haplotype():
    """A haplotype carrying no allele is "." only if every record is uncallable there."""
    rows, _ = _combine([
        _rec(900, 'A', 'C', '1|.', 'snv_c'),
        _rec(900, 'A', 'G', '0|.', 'snv_g'),
    ])

    assert rows[0]['S'] == '1|.'

    rows, _ = _combine([
        _rec(1000, 'ACG', 'A', '1|.', 'del'),
        _rec(1000, 'A', 'T', '0|0', 'snv'),
    ])

    assert rows[0]['S'] == '1|0'
    assert rows[0]['ALT'] == 'A,TCG'


def test_identical_alleles_deduplicated():
    """Identical alleles on both haplotypes (left separate by the merge) become one ALT."""
    rows, _ = _combine([
        _rec(1100, 'A', 'AT', '1|0', 'ins_h1'),
        _rec(1100, 'A', 'AT', '0|1', 'ins_h2'),
    ])

    assert len(rows) == 1
    assert (rows[0]['ALT'], rows[0]['S'], rows[0]['ID']) == ('AT', '1|1', 'ins_h1;ins_h2')
    assert rows[0]['INFO'] == 'SVLEN=1;QUERY=q:ins_h1'


def test_symbolic_and_single_records_unchanged():
    """Symbolic records are never combined, and single records pass through."""
    records = [
        _rec(1200, 'A', '<INV>', '1|0', 'inv', sym=True, info=['END=1300']),
        _rec(1200, 'A', 'C', '0|1', 'snv'),
        _rec(1250, 'A', 'G', '1|1', 'snv2'),
    ]

    rows, stats = _combine(records)

    assert [row['ID'] for row in rows] == ['inv', 'snv', 'snv2']
    assert rows[0]['INFO'] == 'END=1300'
    assert stats['n_combined_site'] == 0
    assert stats['n_record_out'] == 3


def test_format_separator():
    """FORMAT keys and sample values are joined with colons."""
    df = pl.DataFrame([_rec(1300, 'A', 'C', '1|0', 'snv')]).with_columns(
        pl.lit(['GT', 'XX']).alias('_vcf_format'),
        pl.lit(['1|0', '5']).alias(SAMPLE_COL),
    )

    row = reformat_vcf_table(df.lazy()).collect().row(0, named=True)

    assert (row['FORMAT'], row['S']) == ('GT:XX', '1|0:5')


def test_extended_allele_overlapping_own_variant_not_combined():
    """Extending the h1 SNV to the h2 deletion's REF would overlap the next h1 SNV inside the deletion."""
    rows, stats = _combine([
        _rec(1400, 'A', 'G', '1|0', 'snv'),
        _rec(1400, 'ACGTA', 'A', '0|1', 'del'),
        _rec(1402, 'G', 'T', '1|0', 'snv_inner'),
    ])

    assert [row['ID'] for row in rows] == ['snv', 'del', 'snv_inner']
    assert stats['n_span_conflict_site'] == 1
    assert stats['n_combined_site'] == 0


def test_variant_inside_longest_allele_does_not_block():
    """A variant on the haplotype carrying the longest REF is not affected by the combine (REF is not extended)."""
    rows, stats = _combine([
        _rec(1500, 'A', 'G', '1|0', 'snv'),
        _rec(1500, 'ACGTA', 'A', '0|1', 'del'),
        _rec(1502, 'G', 'T', '0|1', 'snv_inner'),
    ])

    assert [row['ID'] for row in rows] == ['snv;del', 'snv_inner']
    assert stats['n_combined_site'] == 1


def test_records_without_query_region():
    """Deletions derived from complex events (no query region, so no QUERY field) at one site are combined."""
    df = pl.DataFrame([
        _rec(100, 'A' + 'C' * 12, 'A', '1|0', 'del12', filter_='DERIVED'),
        _rec(100, 'A' + 'C' * 10, 'A', '0|1', 'del10', filter_='DERIVED'),
    ]).drop('_vcf_info').with_columns(
        pl.Series('vartype', ['DEL', 'DEL']),
        pl.Series('varlen', [12, 10]),
        pl.Series('qry_id', ['q1', 'q2']),
        pl.Series('qry_pos', [None, None], dtype=pl.Int64),
        pl.Series('qry_end', [None, None], dtype=pl.Int64),
    )

    rows, stats = _combine(vcf_record_table(standard_info_fields(df.lazy())).collect().to_dicts())

    assert len(rows) == 1
    assert (rows[0]['ALT'], rows[0]['S'], rows[0]['INFO']) == ('A,ACC', '1|2', 'SVLEN=12,10')
    assert stats['n_combined_site'] == 1
