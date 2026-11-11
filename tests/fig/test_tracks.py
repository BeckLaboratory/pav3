"""Tests for variant track AutoSQL definitions (pav3.fig.tracks).

"bedToBigBed" reads the AutoSQL table name as a single identifier. Variant track table names are built from the
assembly and haplotype names, so an assembly name with a dash (e.g. "SAMPLE-01") produced a table line such as
"table PAVCallHapSAMPLE-01h1SvInvdupcpx", and "bedToBigBed" failed with "Expecting ( got 01h1SvInvdupcpx". These
tests check that the table name is a valid identifier while the long description and the "assembly" field keep the
original name.
"""

from pathlib import Path
import re

import polars as pl
import pytest

from pav3.fig.tracks import (
    call_hap_insdelsnv,
    call_hap_invdupcpx,
    write_track_bed,
)

AS_IDENTIFIER = re.compile(r'[A-Za-z_][A-Za-z0-9_]*')

ASM_NAMES = [
    ('sample1', 'sample1'),
    ('SAMPLE-01', 'SAMPLE01'),
    ('asm-1-T', 'asm1T'),
    ('asm.v1 x', 'asmv1x'),
]

DF_FAI = pl.DataFrame({'chrom': ['chr1'], 'len': [1_000_000]})

HOM = pl.Struct({'up': pl.Int32, 'dn': pl.Int32})

REGION = pl.Struct({'chrom': pl.String, 'pos': pl.Int64, 'end': pl.Int64})

DUP = pl.List(pl.Struct({'chrom': pl.String, 'pos': pl.Int64, 'end': pl.Int64, 'is_rev': pl.Boolean}))

SEG = pl.List(pl.Struct({
    'chrom': pl.String, 'pos': pl.Int64, 'end': pl.Int64,
    'qry_id': pl.String, 'qry_pos': pl.Int64, 'qry_end': pl.Int64, 'is_rev': pl.Boolean,
}))

VAR_BASE = {
    'chrom': pl.String,
    'pos': pl.Int64,
    'end': pl.Int64,
    'id': pl.String,
    'vartype': pl.String,
    'varsubtype': pl.String,
    'varlen': pl.Int64,
    'filter': pl.List(pl.String),
    'qry_id': pl.String,
    'qry_pos': pl.Int64,
    'qry_end': pl.Int64,
    'qry_rev': pl.Boolean,
    'call_source': pl.String,
    'var_score': pl.Float32,
    'align_source': pl.List(pl.Int32),
}

# Column layouts of the PAV 3 per-haplotype call tables read by the INV/DUP/CPX track
INVDUPCPX_SCHEMA = {
    'inv': {
        **{k: v for k, v in VAR_BASE.items() if k != 'varsubtype'},
        'outer_ref': REGION, 'outer_qry': REGION, 'align_gap': pl.Int32, 'var_index': pl.UInt32,
        'hom_ref': HOM, 'hom_qry': HOM, 'derived': pl.String, 'discord': pl.List(pl.String),
        'inner': pl.List(pl.String), 'seq': pl.String,
    },
    'dup': {
        **VAR_BASE, 'derived': pl.String, 'seq': pl.String,
    },
    'cpx': {
        **VAR_BASE,
        'seg': SEG, 'dup': DUP, 'seg_n': pl.Int32, 'trace_qry': pl.String, 'trace_ref': pl.String,
        'var_index': pl.UInt32, 'hom_ref': HOM, 'hom_qry': HOM, 'seq': pl.String,
    },
    'segment': {
        'chrom': pl.String, 'pos': pl.Int64, 'end': pl.Int64, 'filter': pl.List(pl.String),
        'filter_pass': pl.Boolean, 'is_anchor': pl.Boolean, 'is_aligned': pl.Boolean, 'qry_id': pl.String,
        'qry_pos': pl.Int64, 'qry_end': pl.Int64, 'is_rev': pl.Boolean, 'len_ref': pl.Int64, 'len_qry': pl.Int64,
        'gap_ref': pl.Int64, 'align_index': pl.Int32, 'qry_order': pl.Int32, 'seg_index': pl.UInt32,
        'score': pl.Float32, 'var_index': pl.UInt32, 'id': pl.String,
    },
}


def _insdel_callset(tmp_path: Path) -> Path:
    """Write a one-record INS callset.

    :param tmp_path: Output directory.

    :returns: Path to the callset Parquet file.
    """
    callset_path = tmp_path / 'call_insdel_h1.parquet'

    pl.DataFrame(
        [{
            'chrom': 'chr1', 'pos': 1000, 'end': 1001, 'id': 'chr1-1001-INS-100', 'vartype': 'INS',
            'varsubtype': None, 'varlen': 100, 'filter': [], 'qry_id': 'h1tg000001l', 'qry_pos': 5000,
            'qry_end': 5100, 'qry_rev': False, 'call_source': 'INTRA', 'var_score': 0.0, 'align_source': [0],
        }],
        schema=VAR_BASE,
        orient='row',
    ).write_parquet(callset_path)

    return callset_path


def _assert_as_header(as_lines: list[str], asm_name: str, table_name: str) -> None:
    """Check the AutoSQL table name and long description.

    :param as_lines: AutoSQL lines returned by a track builder.
    :param asm_name: Assembly name the track was built for.
    :param table_name: Expected table name.
    """
    table_kw, as_table_name = as_lines[0].split(' ', 1)

    assert table_kw == 'table'
    assert AS_IDENTIFIER.fullmatch(as_table_name), f'Not an AutoSQL identifier: {as_table_name!r}'
    assert as_table_name == table_name
    assert f'({asm_name}-h1:' in as_lines[1]
    assert as_lines[2] == '('


@pytest.mark.parametrize('asm_name, asm_table_name', ASM_NAMES)
def test_insdel_table_name(tmp_path: Path, asm_name: str, asm_table_name: str) -> None:
    """INS/DEL track table names drop non-identifier characters, records keep the assembly name."""
    chunk_iter, as_lines = call_hap_insdelsnv(
        {
            'asm_name': asm_name,
            'hap': 'h1',
            'callset_path': _insdel_callset(tmp_path),
            'varclass': 'sv',
            'vartype': 'insdel',
        },
        df_fai=DF_FAI,
    )

    _assert_as_header(as_lines, asm_name, f'PAVCallHap{asm_table_name}h1SvInsdel')

    bed_path = tmp_path / 'track.bed'
    assert write_track_bed(chunk_iter, bed_path) == 1

    as_fields = [line.split(' ')[1].rstrip(';') for line in as_lines[3:-1]]
    bed_record = dict(zip(as_fields, bed_path.read_text().rstrip('\n').split('\t'), strict=True))

    assert bed_record['assembly'] == asm_name


@pytest.mark.parametrize('asm_name, asm_table_name', ASM_NAMES)
def test_invdupcpx_table_name(tmp_path: Path, asm_name: str, asm_table_name: str) -> None:
    """INV/DUP/CPX track table names drop non-identifier characters."""
    callset_def = {
        'asm_name': asm_name,
        'hap': 'h1',
        'varclass': 'sv',
        'vartype': 'invdupcpx',
    }

    for key, schema in INVDUPCPX_SCHEMA.items():
        callset_path = tmp_path / f'{key}.parquet'
        pl.DataFrame(schema=schema).write_parquet(callset_path)
        callset_def[f'callset_path_{key}'] = callset_path

    _, as_lines = call_hap_invdupcpx(callset_def, df_fai=DF_FAI)

    _assert_as_header(as_lines, asm_name, f'PAVCallHap{asm_table_name}h1SvInvdupcpx')
