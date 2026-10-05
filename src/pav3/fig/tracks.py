"""Routines for creating tracks."""

__all__ = [
    'TRACK_RESOURCE',
    'call_hap_insdelsnv',
    'call_hap_invdupcpx',
    'align',
    'set_default_as_types',
    'write_track_bed',

]

from collections.abc import (
    Iterable,
    Iterator,
)
import itertools
import importlib.resources
from pathlib import Path
import re
from typing import Optional

import polars as pl
import polars.selectors as cs

from . import const as fig_const
from . import util as fig_util
from ..call import expr

from .. import schema
from .. import util

TRACK_RESOURCE = 'pav3.data.tracks'

MG_STAT_FIELDS: tuple[str, ...] = tuple(
    struct_field.name for struct_field in schema.VARIANT['mg_stat'].inner.fields
)
"""Merge statistic fields in the order they are formatted into the "mg_stat" track column."""

CONTROL_CHAR_PATTERN = r'[\x00-\x1f\x7f]'
"""Control characters (including tab, newline, and carriage return) stripped from BED track fields."""

NUMERIC_AS_TYPES: frozenset[str] = frozenset({'int', 'uint', 'short', 'ushort', 'byte', 'ubyte', 'float'})
"""AutoSQL types that must be backed by a numeric column (array suffixes stripped before matching)."""

ITEM_RGB_FIELDS: frozenset[str] = frozenset({'color'})
"""Fields exempt from the numeric check: UCSC itemRgb fields are declared "uint" but hold "R,G,B" strings."""

ALIGN_TRACK_TABLE_FILENAME = 'alignment_track_fields.tsv'

VAR_TRACK_TABLE_FILENAME = 'variant_track_fields.tsv'

TRIM_DESC = {
    'none': 'None',
    'qry': 'Qry',
    'qryref': 'Qry-Ref',
}

VARCLASS_FMT = {
    'sv': 'SV',
    'indel': 'Indel',
    'svindel': 'SV/Indel',
    'snv': 'SNV',
}

VARTYPE_FMT = {
    'ins': 'INS',
    'del': 'DEL',
    'insdel': 'INS/DEL',
    'snv': 'SNV',
    'cpx': 'CPX',
    'inv': 'INV',
    'dup': 'DUP',
    'invdupcpx': 'INV/DUP/CPX',
}


#
# Track output
#

def write_track_bed(
        chunks: Iterable[pl.LazyFrame],
        bed_path: Path | str,
) -> int:
    """Write track chunks to a BED file for conversion to BigBed.

    Each chunk is collected and appended in turn, so peak memory is bounded by the largest chunk
    instead of the whole track. The track builders emit one chunk per reference chromosome in the
    order "bedToBigBed" expects (chromosomes sorted lexicographically, records sorted by position
    within each chromosome), so appending chunks in order yields a correctly sorted BED.

    :param chunks: Track chunks in output order. Every chunk must have the same columns, already
        cast to strings by `set_default_as_types`.
    :param bed_path: Path of the BED file to write.

    :return: Number of records written.
    """
    n_record = 0

    with open(bed_path, 'wb') as out_file:
        for chunk in chunks:
            df_chunk = chunk.collect()

            if df_chunk.height == 0:
                continue

            df_chunk.write_csv(out_file, separator='\t', include_header=False)

            n_record += df_chunk.height

    return n_record


#
# Variant tracks
#

def call_hap_insdelsnv(
        callset_def: dict[str, str | Path],
        df_fai: pl.LazyFrame | pl.DataFrame,
        field_path: Optional[Path | str] = None,
        thick_lines: bool = False,
        center_ins: bool = True,
) -> tuple[Iterator[pl.LazyFrame], list[str]]:
    """Create a track of variant calls.

    The track is returned as one chunk per reference chromosome in output order (see
    `write_track_bed`) rather than as a single frame. Sorting is a blocking operation, so a
    whole-track frame must be fully materialized before it can be written, and every column is
    converted to a string for BED output, which costs far more than the source table (each string
    column carries an 8-byte offset per row). Chunking by chromosome bounds peak memory to the
    largest chromosome.

    :param callset_def: Callset definition.
    :param df_fai: Reference table with "chrom" and "len" columns.
    :param field_path: Path to field table describing track fields.
    :param thick_lines: If True, draw variants with thick lines. Otherwise, use thick for insertion site.
    :param center_ins: If True, center insertions on the variant position. Otherwise, use the start position.

    :return: Tuple of an iterator over per-chromosome variant LazyFrames and AS file lines.
    """
    if isinstance(df_fai, pl.DataFrame):
        df_fai = df_fai.lazy()

    asm_name = asm_name if (asm_name := callset_def.get('asm_name', '').strip()) else 'Unknown'
    hap = hap if (hap := callset_def.get('hap', '').strip()) else 'Unknown'
    vartype = (vartype if (vartype := callset_def.get('vartype', '').strip()) else 'unknown').lower()
    varclass = (varclass if (varclass := callset_def.get('varclass', '').strip()) else 'svindel').lower()

    track_desc_short = f'PAVCallHap{asm_name}{hap}{varclass.capitalize()}{vartype.capitalize()}'
    track_desc_long = f'PAV Call Hap ({asm_name}-{hap}: {VARCLASS_FMT[varclass]} {VARTYPE_FMT[vartype]})'

    df_as = _read_df_as(VAR_TRACK_TABLE_FILENAME, field_path)
    color_table = _read_var_color_table()

    # Read callset
    callset_path = Path(p) if (p := callset_def.get('callset_path', None)) is not None else None

    if callset_path is None or not callset_path.is_file():
        raise FileNotFoundError(f'Callset file not found: {callset_path}')

    # Format table
    filter_list = []

    if varclass == 'sv':
        filter_list.append(pl.col('varlen') >= 50)
    elif varclass == 'indel':
        filter_list.append(pl.col('varlen') < 50)
    elif varclass not in {'svindel', 'snv'}:
        raise ValueError(f'Unknown variant class: {varclass}')

    df = (
        pl.scan_parquet(callset_path)
        .filter(*filter_list)
        .with_columns((pl.col('filter').list.len() == 0).alias('_is_pass'))
        .join(color_table, on=['vartype', '_is_pass'], how='left')
        .join(
            df_fai.select(
                'chrom', pl.col('len').alias('_chrom_len')
            ),
            on='chrom', how='left',
        )
    )

    # Add standard fields
    df = df.with_columns(
        *_var_fields(df, asm_name, hap)
    )

    # Set POS and END (including thick lines)
    df = _set_lines(
        df,
        thick_lines,
        center_ins,
        vartype == 'snv',
    )

    df = df.drop(
        '_is_pass', '_chrom_len',
        'qry_id', 'qry_pos', 'qry_end', 'qry_rev',
        strict=False,
    )

    # Check columns
    col_set = set(df.collect_schema().names())

    if (missing_cols := set(df_as.filter('required').select('field').to_series()) - col_set):
        n = len(missing_cols)
        missing_cols = ', '.join(sorted(missing_cols)[:3]) + ('...' if n > 0 else '')

        raise ValueError(f'Missing {n} columns in alignment track definition: {missing_cols}')

    df_as = df_as.filter(
        pl.col('field').is_in(col_set)
    )

    chunk_iter = _chrom_chunks(df, df_as, _chrom_list(callset_path))

    # AS file lines
    as_lines = [
        f'table {track_desc_short}',
        f'"{track_desc_long}"',
        '('
    ] + [
        '{type} {name}; "{desc}"'.format(**row)
        for row in df_as.iter_rows(named=True)
    ] + [
        ')'
    ]

    return chunk_iter, as_lines


def call_hap_invdupcpx(
        callset_def: dict[str, str | Path],
        df_fai: pl.LazyFrame | pl.DataFrame,
        field_path: Optional[Path | str] = None,
) -> tuple[Iterator[pl.LazyFrame], list[str]]:
    """Create a track of inversion, duplication, and complex variant calls.

    See `call_hap_insdelsnv` for why the track is chunked by chromosome.

    :param callset_def: Callset definition.
    :param df_fai: Reference table with "chrom" and "len" columns.
    :param field_path: Path to field table describing track fields.

    :return: Tuple of an iterator over per-chromosome variant LazyFrames and AS file lines.
    """
    if isinstance(df_fai, pl.DataFrame):
        df_fai = df_fai.lazy()

    asm_name = asm_name if (asm_name := callset_def.get('asm_name', '').strip()) else 'Unknown'
    hap = hap if (hap := callset_def.get('hap', '').strip()) else 'Unknown'
    vartype = (vartype if (vartype := callset_def.get('vartype', '').strip()) else 'unknown').lower()
    varclass = (varclass if (varclass := callset_def.get('varclass', '').strip()) else 'invdupcpx').lower()

    track_desc_short = f'PAVCallHap{asm_name}{hap}{varclass.capitalize()}{vartype.capitalize()}'
    track_desc_long = f'PAV Call Hap ({asm_name}-{hap}: {VARCLASS_FMT[varclass]} {VARTYPE_FMT[vartype]})'

    df_as = _read_df_as(VAR_TRACK_TABLE_FILENAME, field_path)
    color_table = _read_var_color_table()

    # Read callset
    callset_path = {}

    for key in ('inv', 'dup', 'cpx', 'segment'):
        callset_path[key] = Path(p) if (p := callset_def.get(f'callset_path_{key}', None)) is not None else None

        if callset_path[key] is None or not callset_path[key].is_file():
            raise FileNotFoundError(f'Callset file not found for {key}: {callset_path[key]}')

    # Read tables, initial formatting
    df_inv = (
        pl.scan_parquet(callset_path['inv'])
        .with_columns(
            pl.concat_str(
                pl.lit('REF='),
                pl.col('chrom'), pl.lit(':'), pl.col('pos'), pl.lit('-'), pl.col('end'),
                pl.lit(' - QRY='),
                expr.qry_region_expr()
            ).alias('inv_inner'),
            pl.concat_str(
                pl.lit('REF='),
                pl.col('outer_ref').struct.field('chrom'),
                pl.lit(':'),
                pl.col('outer_ref').struct.field('pos'),
                pl.lit('-'),
                pl.col('outer_ref').struct.field('end'),
                pl.lit(' - QRY='),
                pl.col('outer_qry').struct.field('chrom'),
                pl.lit(':'),
                pl.col('outer_qry').struct.field('pos'),
                pl.lit('-'),
                pl.col('outer_qry').struct.field('end')
            ).alias('inv_outer'),
            pl.col('pos').alias('pos_thick'),
            pl.col('end').alias('end_thick'),
            pl.col('outer_ref').struct.field('pos').alias('pos'),
            pl.col('outer_ref').struct.field('end').alias('end'),
            pl.col('align_gap').cast(pl.String),
        )
    )

    df_cpx = (
        pl.scan_parquet(callset_path['cpx'])
        .join(df_fai.select('chrom', pl.col('len').alias('_chrom_len')), on='chrom', how='left',)
        .with_columns(
            pl.col('pos').alias('pos_thick'),
            pl.col('end').alias('end_thick'),
            (
                (pl.col('pos') + pl.col('end')) // 2 + pl.col('varlen') // 2
            ).clip(
                pl.col('end'), pl.col('_chrom_len')
            ).alias('end')
        )
        .with_columns(
            (pl.col('end') - pl.col('varlen')).clip(0, pl.col('pos_thick')).alias('pos')
        )
        .drop('_chrom_len')
    )

    df_seg = pl.scan_parquet(callset_path['segment'])

    df_seg = (
        df_seg
        .drop('filter')
        .join(
            df_cpx.select(
                'var_index',
                pl.col('qry_rev').alias('_cpx_is_rev'),
                pl.col('filter')
            ),
            on='var_index', how='left'
        )
        .join(
            (
                df_seg
                .group_by('var_index')
                .agg(pl.col('seg_index').max().alias('_n_seg'))
            ),
            on='var_index', how='left'
        )
        .filter(pl.col('_cpx_is_rev').is_not_null())  # Only CPX segments, drop INV segments
        .sort(['var_index', 'seg_index'])
        .with_columns(
            pl.col('pos').fill_null(strategy='forward').alias('_pos_fwd'),
            pl.col('pos').fill_null(strategy='backward').alias('_pos_bak'),
            pl.col('chrom').fill_null(strategy='forward').alias('_chrom_fwd'),
            pl.col('chrom').fill_null(strategy='backward').alias('_chrom_bak'),
        )
        .with_columns(
            pl.col('pos').fill_null(
                pl.when(pl.col('_cpx_is_rev'))
                .then('_pos_bak')
                .otherwise('_pos_fwd')
            ),
            pl.col('chrom').fill_null(
                pl.when(pl.col('_cpx_is_rev'))
                .then('_chrom_bak')
                .otherwise('_chrom_fwd')
            ),
            (pl.col('is_rev') ^ pl.col('_cpx_is_rev')).alias('is_rev'),
        )
        .filter(~ pl.col('is_anchor'))
        .select(
            'chrom', 'pos',
            pl.col('end').fill_null(pl.col('pos') + 1).alias('end'),
            (
                pl.when(~ pl.col('is_aligned'))
                .then(pl.lit('CPX_INS'))
                .when(pl.col('is_rev'))
                .then(pl.lit('CPX_INVDUP'))
                .otherwise(pl.lit('CPX_DUP'))
                .alias('vartype')
            ),
            pl.concat_str('id', pl.lit(':'), 'seg_index').alias('id'),
            'qry_id', 'qry_pos', 'qry_end',
            'filter',
            (
                pl.when(pl.col('align_index').is_not_null())
                .then(pl.concat_list('align_index'))
                .otherwise(pl.lit([]))
            ).alias('align_source'),
            pl.concat_str(
                pl.lit('Segment '), pl.col('seg_index') + 1, pl.lit('/'), '_n_seg',
                pl.lit(', len_ref='), 'len_ref', pl.lit(', len_qry='), 'len_qry',
            ).alias('cpx_info'),
            'var_index'
        )
        .with_columns(
            pl.col('pos').alias('pos_thick'),
            pl.col('end').alias('end_thick'),
        )
    )

    df_dup = (
        pl.scan_parquet(callset_path['dup']).with_row_index('_index')
        .join(df_cpx, left_on='derived', right_on='id', how='anti')
    )

    # Concat tables
    df = (
        pl.concat(
            [
                df_inv,
                df_cpx,
                df_dup,
                (
                    df_seg
                    .join(
                        df_cpx, on='var_index', how='left'
                    )
                    .drop(
                        cs.ends_with('_right'),
                        'seg', 'dup', 'trace_ref', 'trace_qry', 'hom_ref', 'hom_qry', 'seq'
                    )
                ),
            ],
            how='diagonal',
        )
        .with_columns(
            (pl.col('filter').list.len() == 0).alias('_is_pass'),
            pl.col('pos_thick').fill_null(pl.col('pos')),
            pl.col('end_thick').fill_null(pl.col('end')),
        )
        .join(color_table, on=['vartype', '_is_pass'], how='left')
    )

    df = (
        df.with_columns(
            *_var_fields(df, asm_name, hap)
        )
    )

    # Check columns
    col_set = set(df.collect_schema().names())

    if (missing_cols := set(df_as.filter('required').select('field').to_series()) - col_set):
        n = len(missing_cols)
        missing_cols = ', '.join(sorted(missing_cols)[:3]) + ('...' if n > 0 else '')

        raise ValueError(f'Missing {n} columns in alignment track definition: {missing_cols}')

    df_as = df_as.filter(
        pl.col('field').is_in(col_set)
    )

    chunk_iter = _chrom_chunks(
        df, df_as,
        _chrom_list(
            *(callset_path[key] for key in ('inv', 'dup', 'cpx'))
        )
    )

    # AS file lines
    as_lines = [
        f'table {track_desc_short}',
        f'"{track_desc_long}"',
        '('
    ] + [
        '{type} {name}; "{desc}"'.format(**row)
        for row in df_as.iter_rows(named=True)
    ] + [
        ')'
    ]

    return chunk_iter, as_lines


#
# Alignment tracks
#

def align(
        align_defs: list[dict[str, str | Path]],
        field_path: Optional[Path | str] = None,
) -> tuple[pl.LazyFrame, list[str]]:
    """Create alignment track.

    :param trim: Trim type (none, qry, qryref).
    :param field_path: Path to field table describing track fields.

    :return: Tuple of alignment DataFrame and AS file lines.
    """
    mpl = util.require_optional_module('matplotlib')

    # Get track description
    track_desc_short, track_desc_long = get_track_desc_align(
        'PAV Align',
        {
            re.sub(r'[^a-zA-Z0-9_]', '', str(d.get('trim', None)).strip().lower())
            for d in align_defs if d.get('trim', None) is not None
        }
    )

    # Read field table
    df_as = _read_df_as(ALIGN_TRACK_TABLE_FILENAME, field_path)

    # Initialise colors
    colormap = mpl.colormaps[fig_const.ALIGN_COLORMAP]
    rotate_iter = rotate_map()

    # Read alignments
    df_list = list()

    for i, align_def in enumerate(align_defs):
        align_path = Path(p) if (p := align_def.get('align_path', None)) is not None else None
        hap = align_def.get('hap', None)
        asm_name = align_def.get('asm_name', '')
        trim = TRIM_DESC.get(align_def.get('trim', None), 'Unknown')

        color = colormap(next(rotate_iter))

        color_pass = pl.lit(
            fig_util.color_to_ucsc_string(color)
        ).cast(pl.String)

        color_fail = pl.lit(
            fig_util.color_to_ucsc_string(
                fig_util.lighten_color(color, 0.25)
            )
        ).cast(pl.String)

        if hap is None or not (hap := str(hap).strip()):
            raise ValueError(f'Haplotype not found for alignment definition index {i}')

        if align_path is None or not align_path.is_file():
            raise FileNotFoundError(f'Alignment file not found for alignment definition index {i}: {str(align_path)}')

        df_list.append(
            pl.scan_parquet(align_path)
            .with_row_index('_index')
            .drop('asm_name', 'hap', 'trim', 'align_ops', strict=False)
            .with_columns(
                pl.lit(asm_name).alias('asm_name'),
                pl.lit(hap).alias('hap'),
                pl.lit(trim).alias('trim'),
                (
                    pl.when(pl.col('filter').list.len() == 0)
                    .then(color_pass)
                    .otherwise(color_fail)
                ).alias('color'),
                pl.lit(i).alias('_align_def_index'),
                expr.qry_region_expr().alias('qry_region'),
            )
        )

    df = (
        pl.concat(df_list)
        .sort('chrom', 'pos', 'end', '_align_def_index', '_index')

        .with_columns(
            pl.col('pos').alias('pos_thick'),
            pl.col('end').alias('end_thick'),
            pl.col('filter').list.join(',').cast(pl.String).alias('filter'),
            pl.concat_str(
                'qry_id',
                pl.lit(', hap='),
                'hap',
                pl.lit(', idx='),
                'align_index',
                pl.lit(', qry_ord='),
                'qry_order',
            ).alias('id'),
            (pl.col('match_prop') * 1000).clip(0.0, 1000.0).cast(pl.Int32).alias('track_score'),
            pl.col('rg').fill_null('.').cast(pl.String).alias('rg'),
            (pl.when('is_rev').then(pl.lit('-')).otherwise(pl.lit('+'))).cast(pl.String).alias('strand'),
        )
        .drop('_align_def_index', '_index', 'qry_id', 'qry_pos', 'qry_end', 'is_rev')
        .collect()
    )

    # Check columns
    col_set = set(df.columns)

    if (missing_cols := set(df_as.filter('required').select('field').to_series()) - col_set):
        n = len(missing_cols)
        missing_cols = ', '.join(sorted(missing_cols)[:3]) + ('...' if n > 0 else '')

        raise ValueError(f'Missing {n} columns in variant track table: {missing_cols}')

    df_as = df_as.filter(
        pl.col('field').is_in(df.columns)
    )

    df = df.select(df_as.select('field').to_series().to_list())
    df = set_default_as_types(df, df_as, strict=True)

    # AS file lines
    as_lines = [
        f'table {track_desc_short}',
        f'"{track_desc_long}"',
        '('
    ] + [
        '{type} {name}; "{desc}"'.format(**row)
        for row in df_as.iter_rows(named=True)
    ] + [
        ')'
    ]

    return df.lazy(), as_lines


def get_track_desc_align(
        track_desc: str = 'PAV Align',
        trim: Iterable[str] | str = (),
) -> tuple[str, str]:
    """Get track description with trim levels.

    :param track_desc: Track description.
    :param trim: Trim levels.

    :return: Tuple of short and long track descriptions.
    """

    trim_list = list(dict.fromkeys(
        [
            t2 for t2 in {
                str(t).strip() for t in
                ((trim,) if isinstance(trim, str) else trim)
                if t is not None
            } if t2
        ]
    ).keys())

    trim_list = [  # Each element is a tuple of the identifier and a description
        (t, TRIM_DESC.get(t.lower(), 'Unknown'))
        for t in trim_list
    ]

    track_desc_list = [
        t for t in re.split(r'[^A-Za-z_]+', track_desc) if t
    ]

    # Get track description
    if len(trim_list) == 1:
        trim, trim_desc = trim_list[0]
        trim_desc_short = trim.capitalize()

    elif len(trim) > 1:
        trim_desc_short = 'Multi'
        trim_desc = '/'.join(t[1] for t in trim_list)

    else:
        trim_desc_short = 'Unknown'
        trim_desc = 'Unknown'

    track_desc_short = f'{"".join([t.capitalize() for t in track_desc_list])}{trim_desc_short}'
    track_desc_long = f'{" ".join(track_desc_list)} (Trim {trim_desc})'

    return track_desc_short, track_desc_long


#
# Track utilities
#

def rotate_map() -> Iterator[float]:
    """Generate intervals in the range 0.0-1.0.

    This function aims to accomplish two goals:
    1) Iterate through a colormap in a way such that haplotype "n" always has the same color
        independent of the number of haplotypes in the callset.
    2) Separate colors as much as possible.

    The naive solution is to iterate through the colormap linearly in steps of 1/n, but this does
    not accomplish the second goal possibly coloring haplotypes in unexpected ways. For example,
    if a sample has a third "unphased" haplotype, then all colors for that sample would be
    different from samples with two haplotypes making it more difficult to visualize.

    :return: Infinite iterator of intervals in the range 0.0-1.0.
    """
    next_offset_iter = itertools.cycle([4, 8, 16, 32, 64])

    width = next(next_offset_iter)
    step = 1 / width

    i = 0.0

    while True:
        if i >= 1.0:
            width = next(next_offset_iter)
            step = 1 / width
            i = step / 2

        yield i
        i += step


def set_default_as_types(
        df: pl.LazyFrame | pl.DataFrame,
        df_as: pl.DataFrame | pl.LazyFrame,
        strict: bool = True,
) -> pl.LazyFrame:
    """Set default values for BED table columns based on AutoSQL types.

    Each column is coerced to the type specified in the AutoSQL table and missing values are
    filled with an appropriate default value for the type ("." for string/char, 0 for integer,
    and 0.0 for float types).

    Columns typed as numbers in the AutoSQL table must hold numeric values, which is checked before
    the cast to string. A column name passed as a string literal to a Polars expression promotes the
    whole column to strings and reaches "bedToBigBed" as a field it cannot parse. UCSC itemRgb fields
    (`ITEM_RGB_FIELDS`) are exempt: they are declared "uint" but hold "R,G,B" strings.

    Control characters are stripped from string-typed columns. A tab or newline reaching a BED
    field breaks "bedToBigBed -tab", which splits records on raw tabs and newlines and knows
    nothing about the CSV quoting Polars applies when writing such a value. An embedded newline
    splits one record across two lines and surfaces as "Expecting N words line L ... got M",
    where M is the 1-based index of the offending field.

    :param df: Table of the BED track table to be converted to BigBed format.
    :param df_as: Table of the data types including columns "field" (column name in `df` and
        "type" (AutoSQL type as a string).
    :param strict: Whether to raise an error if a column is missing in the AutoSQL table. If False,
        unknown columns are skipped.
    """
    if isinstance(df, pl.DataFrame):
        df = df.lazy()

    if isinstance(df_as, pl.LazyFrame):
        df_as = df_as.collect()

    schema = df.collect_schema()

    bracket_pattern = re.compile(r'\[[1-9]\d*]$')

    for col_name, col_type in schema.items():
        try:
            as_select = df_as.filter(pl.col('field') == col_name).select('type')
        except pl.exceptions.ColumnNotFoundError:
            raise ValueError('AS table is missing "field" or "type" columns')

        if as_select.is_empty():
            if strict:
                raise ValueError(f'Missing column "{col_name}" (dtype {col_type}) in AS table')
            continue

        try:
            as_type = as_select.item()
        except ValueError:
            raise ValueError(f'Multiple values for column "{col_name}" in AS table')

        as_base_type = bracket_pattern.sub('', as_type)

        default_value = {
            'string': '.',
            'lstring': '.',
            'char': '.',
            'int': '0',
            'uint': '0',
            'short': '0',
            'ushort': '0',
            'byte': '0',
            'ubyte': '0',
            'float': '0',
        }.get(
            as_base_type,
            None
        )

        if default_value is None:
            raise ValueError(f'Unsupported column type in AS table: {as_type} (column "{col_name}")')

        if (
            as_base_type in NUMERIC_AS_TYPES
            and col_name not in ITEM_RGB_FIELDS
            and not col_type.is_numeric()
        ):
            raise ValueError(
                f'Column "{col_name}" is typed "{as_type}" in the AS table but holds {col_type} values. A '
                f'non-numeric value reaching a numeric BED field surfaces as "Trailing characters parsing '
                f'integer in field N" from "bedToBigBed"'
            )

        col_expr = pl.col(col_name).cast(pl.String)

        if as_base_type in {'string', 'lstring', 'char'}:
            col_expr = col_expr.str.replace_all(CONTROL_CHAR_PATTERN, ' ')

        df = df.with_columns(
            col_expr.fill_null(default_value)
        )

    return df


def _chrom_list(
        *callset_paths: Path | str,
) -> list[str]:
    """Get the sorted list of reference chromosomes present in one or more callsets.

    The list is taken from the callsets rather than from the reference index so that no records are
    silently dropped when a track is chunked by chromosome. Chromosomes are sorted lexicographically
    to match the order "bedToBigBed" expects.

    :param callset_paths: Paths of callset Parquet files. Each must have a "chrom" column.

    :return: Sorted list of chromosome names.
    """
    chrom_set: set[str] = set()

    for callset_path in callset_paths:
        chrom_set |= set(
            pl.scan_parquet(callset_path)
            .select(pl.col('chrom').unique().drop_nulls())
            .collect()
            .to_series()
            .to_list()
        )

    return sorted(chrom_set)


def _chrom_chunks(
        df: pl.LazyFrame,
        df_as: pl.DataFrame,
        chrom_list: Iterable[str],
) -> Iterator[pl.LazyFrame]:
    """Split a track into one sorted chunk per chromosome.

    Records are sorted by position within each chromosome and chunks are emitted in chromosome
    order, which together give the record order "bedToBigBed" expects. The sort is stable, so
    records tied on position keep their order in the source callset.

    :param df: Track table before the AS field selection.
    :param df_as: AS field table filtered to the fields present in `df`.
    :param chrom_list: Chromosomes in output order.

    :return: Iterator over per-chromosome track chunks.
    """
    field_list = df_as.select('field').to_series().to_list()

    for chrom in chrom_list:
        yield set_default_as_types(
            (
                df
                .filter(pl.col('chrom') == chrom)
                .sort('pos', 'end', maintain_order=True)
                .select(field_list)
            ),
            df_as,
            strict=True,
        )


def _read_df_as(
        track_table_filename: str,
        field_path: Optional[Path | str] = None,
) -> pl.DataFrame:
    """Get the table of fields for variant tracks.

    :param track_table_filename: Name of the default track table file in the resource.
    :param field_path: Path to the field table file. If None, the default track table file is used.

    :return: DataFrame of the field table.
    """
    # Read field table
    if field_path is None:
        if not importlib.resources.is_resource(TRACK_RESOURCE, track_table_filename):
            raise FileNotFoundError(f'Track field table not found in resource: {TRACK_RESOURCE}.{track_table_filename}')

        field_path = importlib.resources.files(TRACK_RESOURCE).joinpath(track_table_filename)
    else:
        field_path = Path(field_path)

    if not field_path.is_file():
        raise FileNotFoundError(f'Field table not found: {field_path}')

    return pl.read_csv(field_path, separator='\t')


def _read_var_color_table():
    return pl.DataFrame(
        [
            {
                'vartype': vartype,
                '_is_pass': is_pass,
                'color': (
                    fig_util.color_to_ucsc_string(
                        color if is_pass else fig_util.lighten_color(color, 0.25)
                    )
                )
            }
            for vartype, color in fig_const.COLOR_VARTYPE.items()
            for is_pass in (True, False)
        ],
        orient='row',
        schema={
            'vartype': pl.String,
            '_is_pass': pl.Boolean,
            'color': pl.String,
        },
    ).lazy()


def _var_fields(
        df: pl.LazyFrame,
        asm_name: str,
        hap: str,
) -> list[pl.Expr]:
    col_set = set()

    if df is not None:
        col_set |= set(df.collect_schema().names())

    field_list = [
        pl.col('color').fill_null(fig_util.color_to_ucsc_string((0, 0, 0))),
        pl.col('filter').list.join(',').cast(pl.String).alias('filter'),
        pl.lit(asm_name).alias('asm_name'),
        pl.lit(hap).alias('hap'),
        expr.qry_region_expr().alias('qry_region'),
        pl.col('align_source').cast(pl.List(pl.String)).list.join(',').cast(pl.String).alias('align_source'),
        pl.lit(1000).alias('track_score'),
    ]

    # Variants have no strand: an INV or a reverse-oriented CPX segment is a property of the call, not an
    # orientation of the reference feature, so every variant track record is unstranded.
    field_list.append(
        pl.lit('.').alias('strand')
    )

    if 'call_source' in col_set:
        field_list.append(
            pl.col('call_source').fill_null('.')
        )

    if 'varsubtype' in col_set:
        field_list.append(
            pl.col('varsubtype').fill_null('.')
        )

    if 'hom_ref' in col_set:
        field_list.append(
            pl.when(pl.col('hom_ref').is_not_null())
            .then(pl.concat_str(
                pl.col('hom_ref').struct.field('up'),
                pl.lit(','),
                pl.col('hom_ref').struct.field('dn'),
            ))
            .otherwise(pl.lit('.'))
            .alias('hom_ref')
        )

    if 'hom_qry' in col_set:
        field_list.append(
            pl.when(pl.col('hom_qry').is_not_null())
            .then(pl.concat_str(
                pl.col('hom_qry').struct.field('up'),
                pl.lit(','),
                pl.col('hom_qry').struct.field('dn'),
            ))
            .otherwise(pl.lit('.'))
            .alias('hom_qry')
        )

    if 'discord' in col_set:
        field_list.append(pl.col('discord').fill_null([]).list.join(','))

    if 'inner' in col_set:
        field_list.append(pl.col('inner').fill_null([]).list.join(','))

    if 'dup' in col_set:
        field_list.append(
            pl.col('dup').fill_null([]).list.eval(
                pl.concat_str(
                    pl.element().struct.field('chrom'),
                    pl.lit(':'),
                    pl.element().struct.field('pos'),
                    pl.lit('-'),
                    pl.element().struct.field('end'),
                    pl.lit('('),
                    (
                        pl.when(pl.element().struct.field('is_rev'))
                        .then(pl.lit('-'))
                        .otherwise(pl.lit('+'))
                    ),
                    pl.lit(')'),
                )
            ).list.join(',')
        )

    if 'seg' in col_set:
        field_list.append(
            pl.col('seg').fill_null([]).list.eval(
                pl.concat_str(
                    pl.element().struct.field('chrom'),
                    pl.lit(':'),
                    pl.element().struct.field('pos'),
                    pl.lit('-'),
                    pl.element().struct.field('end'),
                    pl.lit('('),
                    (
                        pl.when(pl.element().struct.field('is_rev'))
                        .then(pl.lit('-'))
                        .otherwise(pl.lit('+'))
                    ),
                    pl.lit(')>'),
                    pl.element().struct.field('qry_id'),
                    pl.lit(':'),
                    pl.element().struct.field('qry_pos'),
                    pl.lit('-'),
                    pl.element().struct.field('qry_end'),
                )
            ).list.join(',')
        )

    # Merge fields (merged callsets only, see pav3.schema.VARIANT)
    if 'mg_src' in col_set:
        field_list.append(
            pl.col('mg_src').fill_null([]).list.eval(
                pl.concat_str(
                    pl.element().struct.field('sample'),
                    pl.lit(':'),
                    pl.element().struct.field('hap'),
                    pl.lit(':'),
                    pl.element().struct.field('var_id'),
                )
            ).list.join(',')
        )

    if 'mg_stat' in col_set:
        stat_expr_list: list[pl.Expr] = []

        for stat_field in MG_STAT_FIELDS:
            if stat_expr_list:
                stat_expr_list.append(pl.lit(':'))

            stat_expr_list.append(pl.element().struct.field(stat_field))

        field_list.append(
            pl.col('mg_stat').fill_null([]).list.eval(
                # Statistics a merge does not compute are null (e.g. "match_prop" for variant types merged without
                # sequence matching). "ignore_nulls" leaves an empty slot between separators instead of nulling the
                # whole entry, which would drop it from the joined list and misalign the remaining fields.
                pl.concat_str(*stat_expr_list, ignore_nulls=True)
            ).list.join(',')
        )

    return field_list


def _set_lines(
        df: pl.LazyFrame,
        thick_lines: bool = False,
        center_ins: bool = True,
        is_snv: bool = False,
) -> pl.LazyFrame:
    """
    Set POS, END, POS_THICK, END_THICK based on variant type and thick-line settings.

    :param df: Input DataFrame
    :param thick_lines: Whether to use thick lines
    :param center_ins: Whether to center INS variants (only if thick_lines is False)

    :return: Table with positions set for a track.
    """

    if is_snv:
        return df.with_columns(
            pl.col('pos').alias('pos_thick'),
            pl.col('end').alias('end_thick'),
        )

    if thick_lines:
        return (
            df
            .with_columns(
                (
                    pl.when(pl.col('vartype') == 'INS')
                    .then(pl.col('pos') + pl.col('varlen'))
                    .otherwise(pl.col('end'))
                )
                .clip(0, pl.col('_chrom_len'))
                .alias('end'),
            )
            .with_columns(
                pl.min_horizontal('pos', 'end').alias('pos_thick'),
                pl.max_horizontal('pos', 'end').alias('end_thick'),
            )
        )

    # Thin lines
    divider = 2 if center_ins else 1

    return (
        df
        .with_columns(
            pl.col('pos').alias('pos_thick'),
            (
                pl.when(pl.col('vartype') == 'INS')
                .then('end')
                .otherwise('pos')  # No thick part of pos_thick == end_thick
                .alias('end_thick')
            ),
            (
                pl.when(pl.col('vartype') == 'INS')
                .then(
                    (pl.col('end') + pl.col('varlen') // divider)
                    .clip(0, pl.col('_chrom_len'))
                )
                .otherwise('end')
                .alias('end')
            )
        )
        .with_columns(
            pl.when(pl.col('vartype') == 'INS')
            .then(
                # Shift back if centered or runs off chrom
                (pl.col('end') - pl.col('varlen')).clip(0)
            )
            .otherwise('pos')
            .alias('pos')
        )
    )
