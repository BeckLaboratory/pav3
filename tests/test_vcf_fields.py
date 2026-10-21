"""Tests for VCF field construction in :mod:`pav3.vcf`."""
from __future__ import annotations

import polars as pl
import pysam
import pytest

from pav3.vcf import (
    get_headers, reformat_vcf_table, sort_vcf_table, standard_info_fields, vcf_fields_seq_insdel, vcf_fields_snv,
    vcf_fields_sym,
)


REF_SEQ = 'ACGT' * 250


@pytest.fixture
def ref_fa(tmp_path):
    """Write a small indexed reference FASTA."""
    path = tmp_path / 'ref.fa'
    path.write_text(f'>chr1\n{REF_SEQ}\n')
    pysam.faidx(str(path))

    return path


def _frame(vartype: str, varsubtype: str | None) -> pl.LazyFrame:
    """Build a one-row symbolic variant table."""
    return pl.LazyFrame({
        'chrom': ['chr1'],
        'vartype': [vartype],
        'varsubtype': [varsubtype],
        'pos': [10],
        'end': [20],
        'varlen': [10],
        'seq': ['acgt'],
    }, schema_overrides={'varlen': pl.Int64, 'varsubtype': pl.String}).with_row_index('_index')


def _info(df: pl.LazyFrame) -> str:
    """Render the INFO column for a one-row table."""
    out = standard_info_fields(df.with_columns(pl.lit(True).alias('_is_sym_alt'))).collect()

    return ';'.join(out['_vcf_info'].to_list()[0])


@pytest.mark.parametrize(
    ('vartype', 'varsubtype', 'alt'),
    [
        ('DUP', 'TD', '<DUP:TANDEM>'),
        ('DUP', 'INVDUP', '<DUP:INVDUP>'),
        ('DUP', None, '<DUP>'),
        ('CPX', 'DEL', '<CPX>'),
        ('INV', None, '<INV>'),
    ],
)
def test_duplication_subtype_moves_to_alt(ref_fa, vartype, varsubtype, alt):
    """Duplication subtypes are written into the ALT allele, other types keep a plain allele."""
    assert vcf_fields_sym(_frame(vartype, varsubtype), ref_fa).collect()['_vcf_alt'].to_list() == [alt]


def test_subtype_not_repeated_in_info():
    """A subtype carried by the ALT allele is not also written to INFO/SUBTYPE."""
    assert 'SUBTYPE=' not in _info(_frame('DUP', 'TD'))
    assert 'SUBTYPE=' not in _info(_frame('DUP', 'INVDUP'))


def test_subtype_kept_in_info_for_other_types():
    """Subtypes the ALT allele does not carry stay in INFO/SUBTYPE."""
    assert 'SUBTYPE=DEL' in _info(_frame('CPX', 'DEL'))
    assert 'SUBTYPE=TD' in _info(_frame('INS', 'TD'))


def test_null_subtype_is_not_written():
    """A null subtype produces no INFO/SUBTYPE."""
    assert 'SUBTYPE=' not in _info(_frame('DUP', None))
    assert 'SUBTYPE=' not in _info(_frame('INV', None))


def test_alt_headers_declare_duplication_subtypes():
    """Subtyped duplication alleles are declared in the header."""
    alt_set = {line for line in get_headers() if line.startswith('##ALT')}

    assert any('ID=DUP:TANDEM' in line for line in alt_set)
    assert any('ID=DUP:INVDUP' in line for line in alt_set)


def test_pass_filter_is_not_declared():
    """PASS is implicitly defined by the specification and is not written to the header."""
    assert not any(line.startswith('##FILTER=<ID=PASS,') for line in get_headers())


def _info_for(vartype: str, is_sym: bool) -> str:
    """Render the INFO column for one variant of a given type and ALT style."""
    df = pl.LazyFrame({
        'chrom': ['chr1'], 'vartype': [vartype], 'varsubtype': [None],
        'pos': [10], 'end': [20], 'varlen': [10], 'seq': ['acgt'],
    }, schema_overrides={'varlen': pl.Int64, 'varsubtype': pl.String}).with_columns(
        pl.lit(is_sym).alias('_is_sym_alt')
    )

    return ';'.join(standard_info_fields(df).collect()['_vcf_info'].to_list()[0])


@pytest.mark.parametrize('vartype', ['DEL', 'DUP'])
def test_svclaim_on_symbolic_del_and_dup(vartype):
    """VCF requires a claim for symbolic DEL and DUP. PAV calls from adjacency, never abundance."""
    assert 'SVCLAIM=J' in _info_for(vartype, is_sym=True)


@pytest.mark.parametrize('vartype', ['INS', 'INV', 'CPX'])
def test_svclaim_omitted_for_other_symbolic_types(vartype):
    """INS and INV carry an implicit adjacency claim and CPX is not a specification allele type."""
    assert 'SVCLAIM' not in _info_for(vartype, is_sym=True)


@pytest.mark.parametrize('vartype', ['DEL', 'DUP', 'INS'])
def test_svclaim_omitted_for_sequence_resolved(vartype):
    """A sequence-resolved REF and ALT pair leaves nothing to disambiguate."""
    assert 'SVCLAIM' not in _info_for(vartype, is_sym=False)


def _info_list(**cols) -> list[str | None]:
    """Get the INFO list for one sequence-resolved deletion with extra columns."""
    df = pl.LazyFrame({
        'chrom': ['chr1'], 'vartype': ['DEL'], 'varsubtype': [None],
        'pos': [10], 'end': [20], 'varlen': [10], 'seq': ['acgt'], '_is_sym_alt': [False],
        **{col: [val] for col, val in cols.items()},
    }, schema_overrides={
        'varlen': pl.Int64, 'varsubtype': pl.String, 'qry_pos': pl.Int64, 'qry_end': pl.Int64,
        'hom_ref': pl.Struct({'up': pl.Int32, 'dn': pl.Int32}), 'hom_qry': pl.Struct({'up': pl.Int32, 'dn': pl.Int32}),
    })

    return standard_info_fields(df).collect()['_vcf_info'].to_list()[0]


def test_query_written_with_query_coordinates():
    """QUERY is the 1-based query region."""
    assert 'QUERY=q1:6-15' in _info_list(qry_id='q1', qry_pos=5, qry_end=15)


def test_query_omitted_without_query_coordinates():
    """Variants derived from complex events have no query region and write no QUERY (previously a null field)."""
    info = _info_list(qry_id='q1', qry_pos=None, qry_end=None)

    assert None not in info
    assert not any(field.startswith('QUERY=') for field in info)


def test_null_info_field_dropped():
    """A field built from a missing value is left out rather than written as a null field."""
    info = _info_list(hom_ref={'up': 3, 'dn': None}, hom_qry={'up': 0, 'dn': 0})

    assert None not in info
    assert not any(field.startswith('HOMLEN=') for field in info)


@pytest.mark.parametrize('vartype', ['DEL', 'DUP', 'INV'])
def test_svlen_on_symbolic(vartype):
    """SVLEN gives the length of symbolic alleles."""
    assert 'SVLEN=10' in _info_for(vartype, is_sym=True).split(';')


@pytest.mark.parametrize('vartype', ['DEL', 'INS', 'CPX'])
def test_svlen_omitted_for_sequence_resolved(vartype):
    """REF and ALT give the length of sequence-resolved alleles, for which VCF directs writers to omit SVLEN."""
    assert 'SVLEN' not in _info_for(vartype, is_sym=False)


def test_empty_subtype_and_traces_not_written():
    """Empty strings are left out rather than written as a key with no value."""
    info = _info_list(varsubtype='', trace_ref='', trace_qry='')

    assert not any(field.split('=')[0] in {'SUBTYPE', 'CPX_REF_TRACE', 'CPX_QRY_TRACE'} for field in info)


def test_info_values_percent_encoded():
    """Characters with a special meaning in INFO are percent-encoded in names embedded in values."""
    assert 'QUERY=q%3B1%3Da%2Cb%25:6-15' in _info_list(qry_id='q;1=a,b%', qry_pos=5, qry_end=15)


def test_header_lead_lines():
    """The file date uses the conventional "fileDate" key, and the reference is a file name, not a malformed URI."""
    headers = get_headers(ref_filename='/path/to/hs1.fa.gz')

    assert any(line.startswith('##fileDate=') for line in headers)
    assert '##reference=hs1.fa.gz' in headers


def _ref_info() -> pl.DataFrame:
    """Get reference information with sequences in FASTA order chr1, chr2, chr10 (stored out of order)."""
    return pl.DataFrame({
        'chrom': ['chr10', 'chr1', 'chr2'],
        'len': [300, 100, 200],
        'md5': ['c', 'a', 'b'],
        'order': [2, 0, 1],
    })


def test_contig_headers_in_reference_order():
    """Contig lines follow the reference FASTA order."""
    contigs = [line for line in get_headers(df_ref_info=_ref_info()) if line.startswith('##contig')]

    assert [line.split(',')[0] for line in contigs] == ['##contig=<ID=chr1', '##contig=<ID=chr2', '##contig=<ID=chr10']


def test_records_sorted_in_reference_order():
    """Records follow the contig header order, not string order (chr10 before chr2)."""
    df = pl.LazyFrame({
        '#CHROM': ['chr10', 'chr2', 'chr1', 'chr2'],
        'POS': [5, 9, 7, 3],
        'ID': ['d', 'c', 'a', 'b'],
    })

    out = sort_vcf_table(df, _ref_info()).collect()

    assert out['ID'].to_list() == ['a', 'b', 'c', 'd']
    assert out.columns == ['#CHROM', 'POS', 'ID']


def test_empty_info_is_missing_value():
    """A record with no INFO fields writes the missing value."""
    df = pl.LazyFrame({
        '_vcf_chrom': ['chr1'], '_vcf_pos': [10], '_vcf_id': ['v1'], '_vcf_ref': ['A'], '_vcf_alt': ['C'],
        '_vcf_filter': ['PASS'], '_vcf_info': [[]], '_vcf_format': [['GT']], '_vcf_sample_S': [['0|1']],
    }, schema_overrides={'_vcf_info': pl.List(pl.String)})

    assert reformat_vcf_table(df).collect()['INFO'].to_list() == ['.']


def test_snv_iupac_reference_reduced():
    """An ambiguous reference base is written as the base first alphabetically, and an SNV it matches is dropped."""
    df = pl.LazyFrame({
        'pos': [10, 20, 30], 'ref': ['R', 'r', 'y'], 'alt': ['G', 'A', 'g'],
    }).with_row_index('_index')

    out = vcf_fields_snv(df).collect()

    assert out.select('_vcf_ref', '_vcf_alt').rows() == [('A', 'G'), ('C', 'G')]


def test_insdel_iupac_and_soft_masked_reference(tmp_path):
    """Soft-masked and ambiguous reference bases are upper-cased and reduced in REF and in the padding base of ALT."""
    path = tmp_path / 'ref_iupac.fa'
    path.write_text('>chr1\nacgtRYacgtacgtacgt\n')
    pysam.faidx(str(path))

    df = pl.LazyFrame({
        'chrom': ['chr1'], 'pos': [5], 'end': [6], 'vartype': ['DEL'], 'seq': ['y'], 'id': ['del1'],
    }).with_row_index('_index')

    row = vcf_fields_seq_insdel(df, path).collect().row(0, named=True)

    assert (row['_vcf_pos'], row['_vcf_ref'], row['_vcf_alt']) == (5, 'AC', 'A')
