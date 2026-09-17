"""Tests for VCF field construction in :mod:`pav3.vcf`."""
from __future__ import annotations

import polars as pl
import pysam
import pytest

from pav3.vcf import get_headers, standard_info_fields, vcf_fields_sym


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
