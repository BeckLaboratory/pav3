"""Transform variant calls from Parquet files to VCF files."""

__all__ = [
    'InfoField',
    'FilterField',
    'FormatField',
    'AltField',
    'VCF_VERSION',
    'INFO_FIELDS',
    'FILTER_FIELDS',
    'FORMAT_FIELDS',
    'ALT_FIELDS',
    'VARTYPES',
    'VARTYPE_COMPAT',
    'REF_SPAN_VARTYPE',
    'SVCLAIM_VARTYPE',
    'SEQ_VARTYPE',
    'alt_subtype_expr',
    'get_headers',
    'init_vcf_fields',
    'vcf_fields_seq_insdel',
    'vcf_fields_seq_cpx',
    'vcf_fields_sym',
    'vcf_fields_snv',
    'gt_column_asm',
    'standard_info_fields',
    'VCF_RECORD_COLUMNS',
    'vcf_record_table',
    'combine_multiallelic',
    'reformat_vcf_table',
    'sort_vcf_table',
    'resolve_hap_list',
]


import collections
from collections.abc import Iterable
import datetime
import logging
from pathlib import Path
import threading
from typing import Optional

import agglovar
import numpy as np
import polars as pl
import polars.selectors as cs
import pysam

from . import __version__
from .const import FILTER_REASON

logger = logging.getLogger(__name__)

# Note on number fields (PAV does not use most of these)
#
# Number fields:
# * Int: Exact number
# * A: One per ALT allele
# * R: One per allele including REF
# * G: One per each possible genotype
# * .: Varies, unknown, or unbounded
#
# Format fields have additional Number options:
# * LA: A except only ALT alleles in LAA are considered present.
# * LR: R except only ALT alleles in LAA are considered present.
# * LG: G except only ALT alleles in LAA are considered present.
# * P: One value for each allele value in GT.
# * M: One value for each possible base modification.


@agglovar.meta.decorators.immutable
class InfoField:
    name: str = agglovar.meta.descriptors.CheckedString(match=r'[A-Z]([A-Z_]*[A-Z])?')
    number: str = agglovar.meta.descriptors.CheckedString(match=r'([1-9][0-9]*)|[ARG.]')
    type_: str = agglovar.meta.descriptors.CheckedString(match={'Integer', 'Float', 'Flag', 'Character', 'String'})
    description: str = agglovar.meta.descriptors.CheckedString(min_len=1, strip=True)

    def __init__(self, name, number, type_, description):
        self.name = str(name)
        self.number = str(number)
        self.type_ = str(type_)
        self.description = str(description)

        if self.type_ == 'Flag' and self.number != '0':
            raise ValueError(f'Flag type must have number 0: Number={self.number}')

    def __str__(self) -> str:
        return f'##INFO=<ID={self.name},Number={self.number},Type={self.type_},Description="{self.description}">'


@agglovar.meta.decorators.immutable
class FilterField:
    name: str = agglovar.meta.descriptors.CheckedString(match=r'[A-Z]([A-Z_]*[A-Z])?')
    description: str = agglovar.meta.descriptors.CheckedString(min_len=1, strip=True)

    def __init__(self, name, description):
        self.name = str(name)
        self.description = str(description)

    def __str__(self) -> str:
        return f'##FILTER=<ID={self.name},Description="{self.description}">'


@agglovar.meta.decorators.immutable
class FormatField:
    name: str = agglovar.meta.descriptors.CheckedString(match=r'[A-Z]([A-Z_]*[A-Z])?')
    number: str = agglovar.meta.descriptors.CheckedString(match=r'([1-9][0-9]*)|[ARGPM.]|L[ARG]')
    type_: str = agglovar.meta.descriptors.CheckedString(match={'Integer', 'Float', 'Character', 'String'})
    description: str = agglovar.meta.descriptors.CheckedString(min_len=1, strip=True)

    def __init__(self, name, number, type_, description):
        self.name = str(name)
        self.number = str(number)
        self.type_ = str(type_)
        self.description = str(description)

    def __str__(self) -> str:
        return f'##FORMAT=<ID={self.name},Number={self.number},Type={self.type_},Description="{self.description}">'


@agglovar.meta.decorators.immutable
class AltField:
    name: str = agglovar.meta.descriptors.CheckedString(match=r'[A-Z]([A-Z:_]*[A-Z])?')
    description: str = agglovar.meta.descriptors.CheckedString(min_len=1, strip=True)

    def __init__(self, name, description):
        self.name = str(name)
        self.description = str(description)

    def __str__(self) -> str:
        return f'##ALT=<ID={self.name},Description="{self.description}">'


VCF_VERSION = '4.5'

INFO_FIELDS = [
    # InfoField(
    #     'ID', '1', 'String',
    #     'Unique variant ID',
    # ),
    InfoField(
        'SVLEN', 'A', 'Integer',
        'Variant length for symbolic alleles: affected reference span for DEL, INV, and DUP, inserted bases '
        '(query-side) for CPX. Does not equal the POS/END reference span for CPX. Not written for sequence-resolved '
        'alleles, where REF and ALT give the length.',
    ),
    InfoField(
        'END', '1', 'Integer',
        'Variant end position on the reference (1-based, inclusive)',
    ),
    InfoField(
        'SVCLAIM', 'A', 'String',
        'Claim made by the structural variant call. Valid values are D, J, DJ for abundance, adjacency and both '
        'respectively',
    ),
    InfoField(
        'SUBTYPE', 'A', 'String',
        'Variant subtype. TD for an insertion of a tandem duplication. For CPX, the colon-separated segment types of '
        'the reference trace (e.g. DEL, NML, INV, DUP).',
    ),
    InfoField(
        'CPX_REF_TRACE', 'A', 'String',
        'Reference trace for complex variants',
    ),
    InfoField(
        'CPX_QRY_TRACE', 'A', 'String',
        'Query trace for complex variants',
    ),
    InfoField(
        'QUERY', 'A', 'String',
        'Query region of the variant as \\"contig:start-end\\" (1-base, closed coordinates). Deletions report a 1 bp '
        'region at the breakpoint. Taken from one haplotype for variants merged across haplotypes.',
    ),
    InfoField(
        'QUERYREV', 'A', 'Character',
        'Query sequence orientation vs the reference (+ or -, * if unknown)',
    ),
    InfoField(
        'SEQ', 'A', 'String',
        'Variant sequence, written for symbolic alleles only (REF and ALT carry it otherwise). For INV, the inverted '
        'allele from the query. Not written for derived variants.',
    ),
    InfoField(
        'DUP', 'A', 'String',
        'Slash-separated list of regions duplicated by this variant, each region in format '
        '\\"chrom:pos-end\\" (1-base, closed coordinates)'
    ),
    InfoField(
        'HOMLEN', 'A', 'Integer',
        'Length of breakpoint homology, upstream plus downstream. From reference breakpoints for DEL and DUP '
        'and from query breakpoints for INS and CPX. Not reported for SNV or INV.'
    ),
]

# PASS is implicitly defined as the first filter by the VCF specification, so it is not declared here.
FILTER_FIELDS = [
    FilterField(name, description)
    for name, description in FILTER_REASON.items()
]

FORMAT_FIELDS = [
    FormatField(
        'GT', '1', 'String',
        'Genotype',
    ),
]

ALT_FIELDS = [
    AltField('INS', 'Insertion'),
    AltField('DEL', 'Deletion'),
    AltField('INV', 'Inversion'),
    AltField('DUP', 'Duplication'),
    AltField('DUP:TANDEM', 'Tandem duplication'),
    AltField('DUP:INVDUP', 'Inverted duplication'),
    AltField('CPX', 'Complex variant'),
]

VARTYPES = ('insdel', 'inv', 'snv', 'cpx', 'dup')

VARTYPE_COMPAT = {
    'insdel': {'INS', 'DEL'},
    'ins': {'INS'},
    'del': {'DEL'},
    'inv': {'INV'},
    'snv': {'SNV'},
    'cpx': {'CPX'},
    'dup': {'DUP'},
}

_MAX_INVALID_REPORT = 10
"""Maximum number of variant IDs named when REF and ALT fail validation."""

IUPAC_REDUCE = {
    'R': 'A', 'Y': 'C', 'S': 'C', 'W': 'A', 'K': 'G', 'M': 'A', 'B': 'C', 'D': 'A', 'H': 'A', 'V': 'A',
}
"""IUPAC ambiguity codes mapped to the base VCF writes for them.

VCF allows only A, C, G, T, and N in REF and ALT, and an ambiguous reference base is reduced to the base that is first
alphabetically (e.g. R is written as A).
"""

INFO_ENCODE = {
    '%': '%25', ',': '%2C', ';': '%3B', '=': '%3D', '\t': '%09', '\n': '%0A', '\r': '%0D',
}
"""Characters percent-encoded in INFO values.

These have a special meaning in the INFO column, so a literal occurrence in a value (e.g. in a contig name) is encoded.
"""


def _concrete_bases(expr: pl.Expr) -> pl.Expr:
    """Upper-case a sequence and reduce IUPAC ambiguity codes to a single base (see `IUPAC_REDUCE`).

    REF and ALT receive the same reduction so the padding base they share stays identical.

    :param expr: Sequence expression.

    :returns: Expression with the sequence VCF writes in REF or ALT.
    """
    return expr.str.to_uppercase().str.replace_many(IUPAC_REDUCE)


def _encode_info_value(expr: pl.Expr) -> pl.Expr:
    """Percent-encode characters with a special meaning in INFO values (see `INFO_ENCODE`).

    :param expr: String expression.

    :returns: Encoded string expression.
    """
    return expr.str.replace_many(INFO_ENCODE)


def alt_subtype_expr() -> pl.Expr:
    """Get an expression for the subtype carried by a symbolic ALT allele, or null where there is none.

    Duplication subtypes are encoded in the ALT allele ("<DUP:TANDEM>") rather than INFO/SUBTYPE, following the VCF
    convention of writing a symbolic allele as a colon-separated list of types and subtypes. Records resolving to a
    non-null subtype here do not repeat it in INFO/SUBTYPE.

    :returns: Expression resolving to the ALT subtype, or null when the subtype stays in INFO/SUBTYPE.
    """
    return (
        pl.when(pl.col('vartype') != 'DUP')
        .then(None)
        .when(pl.col('varsubtype') == 'TD')
        .then(pl.lit('TANDEM'))
        .when(pl.col('varsubtype') == 'INVDUP')
        .then(pl.lit('INVDUP'))
        .otherwise(None)
        .cast(pl.String)
    )


SEQ_VARTYPE = frozenset({'INS', 'DEL', 'SNV'})
"""Variant types always written with sequence-resolved REF and ALT. Complex variants are also sequence-resolved unless
symbolic complex variants are requested (see `init_vcf_fields`)."""

REF_SPAN_VARTYPE = frozenset({'DEL', 'INV', 'DUP', 'CPX'})
"""Variant types spanning reference bases, which receive an INFO/END annotation.

INS spans no reference bases and SNVs are fully described by REF and ALT, so neither is annotated.
"""

SVCLAIM_VARTYPE = frozenset({'DEL', 'DUP'})
"""Variant types requiring an INFO/SVCLAIM annotation on a symbolic ALT allele.

VCF requires a claim for symbolic DEL and DUP alleles. INS and INV carry an implicit adjacency claim, and a
sequence-resolved REF and ALT pair leaves nothing to disambiguate.
"""


def _ref_info_order(
    df_ref_info: pl.DataFrame,
) -> pl.DataFrame:
    """Order a reference information table by the order of sequences in the reference FASTA.

    :param df_ref_info: Reference information table. Tables without an "order" column are returned as they are.

    :returns: Reference information table in reference order.
    """
    return df_ref_info.sort('order') if 'order' in df_ref_info.columns else df_ref_info


def get_headers(
    ref_filename: Optional[str | Path] = None,
    df_ref_info: Optional[pl.DataFrame] = None,
) -> list[str]:

    # Lead headers
    headers = [
        f'##fileformat=VCFv{VCF_VERSION}',
        f'##fileDate={datetime.datetime.now().strftime("%Y%m%d")}',
        f'##source=PAVv{__version__}',
    ]

    # Reference
    if ref_filename is not None:
        # Only the file name, do not leak filesystem paths for security reasons. Not written as a "file://" URI,
        # which would read the file name as a host name.
        headers.append(
            f'##reference={Path(ref_filename).name}'
        )

    # Reference contigs in reference order (records are sorted in the same order, see `sort_vcf_table`)
    if df_ref_info is not None:
        headers.extend(
            f'##contig=<ID={chrom},length={length},md5={md5}>'
            for chrom, length, md5 in _ref_info_order(df_ref_info).select(['chrom', 'len', 'md5']).iter_rows()
        )

    # ALT headers
    headers.extend(
        str(alt_field) for alt_field in ALT_FIELDS
    )

    # INFO headers
    headers.extend(
        str(info_field) for info_field in INFO_FIELDS
    )

    # FILTER headers
    headers.extend(
        str(filter_field) for filter_field in FILTER_FIELDS
    )

    # FORMAT headers
    headers.extend(
        str(format_field) for format_field in FORMAT_FIELDS
    )

    return headers


def init_vcf_fields(
    df: pl.LazyFrame,
    ref_fa: Optional[str | Path] = None,
    use_sym: Optional[bool] = None,
    vartype: Optional[str] = None,
    cpx_symbolic: bool = False,
    log_stats: Optional[dict] = None,
) -> pl.LazyFrame:
    """Initialize VCF fields.

    :param df: Variant table with an "_index" column.
    :param ref_fa: Reference FASTA filename. Required for non-SNV variants.
    :param use_sym: Whether to use symbolic alternate records. If None, INS, DEL, SNV, and CPX (unless
        `cpx_symbolic` is set) are not symbolic, and all other variant types are.
    :param vartype: Expected variant type or types (i.e. "insdel" allows "INS" and "DEL"). Setting
        this field ensures the variant types is handled consisntently if all variant types
        are not in table `df` or `df` is empty.
    :param cpx_symbolic: Write complex variants as symbolic "<CPX>" records instead of sequence-resolved REF and ALT
        when `use_sym` is None.
    :param log_stats: If not None, add statistics to this dictionary (see `vcf_fields_seq_cpx`).
    """
    seq_vartype_set = SEQ_VARTYPE | ({'CPX'} if not cpx_symbolic else set())

    if df.collect_schema().get('_index') is None:
        df = df.with_row_index('_index')

    vartype_set = set(
        df.select(pl.col('vartype').unique()).collect().to_series()
    )

    if vartype is not None:
        vartype_compat = VARTYPE_COMPAT.get(vartype, None)

        if vartype_compat is None:
            raise ValueError(f'Unknown variant type: {vartype}')

        if no_compat_set := vartype_set - vartype_compat:
            raise ValueError(
                f'Declared variant type "{vartype}" is not compatible with types found '
                f'in the variant table: {", ".join(sorted([str(_).upper() for _ in no_compat_set]))}'
            )

        vartype_set = vartype_compat

    if 'SNV' in vartype_set and len(vartype_set) > 1:
        raise ValueError(
            f'Cannot mix SNV variants with other variant types: '
            f'Found variant types {", ".join(vartype_set)}'
        )

    is_snv = 'SNV' in vartype_set

    if is_snv and use_sym:
        raise ValueError('Variant type "SNV" cannot use symbolic ALTs')

    if use_sym is None:

        no_sym_set = vartype_set & seq_vartype_set
        is_sym_set = vartype_set - seq_vartype_set

        if not (no_sym_set or is_sym_set):
            use_sym = False

        elif no_sym_set and is_sym_set:
            raise ValueError(
                f'Unable to automatically determine whether to use symbolic ALTs or sequence ALTs given variant types '
                f'using sequence ALTs ({", ".join(sorted(seq_vartype_set))}) and symbolic ALTs (all others): '
                f'Found variant types {", ".join(vartype_set)}'
            )
        else:
            use_sym = bool(is_sym_set)

    if use_sym and is_snv:
        raise ValueError('Variant type "SNV" cannot use symbolic ALTs')

    if not use_sym and (vartype_set - SEQ_VARTYPE - {'CPX'}):
        raise ValueError(
            f'Can only use sequence ALTs for INS, DEL, SNV, and CPX: '
            f'Found variant types {", ".join(vartype_set)}'
        )

    # Get POS, REF, and ALT
    if is_snv:
        df_pos_ref_alt = vcf_fields_snv(df)

    elif use_sym:
        if ref_fa is None:
            raise ValueError('Missing reference FASTA filename: Required for symbolic ALTs')

        df_pos_ref_alt = vcf_fields_sym(df, ref_fa)

    else:
        if ref_fa is None:
            raise ValueError('Missing reference FASTA filename: Required for non-SNV sequence ALTs')

        df_pos_ref_alt_list = []

        if vartype_set & {'INS', 'DEL'}:
            df_pos_ref_alt_list.append(
                vcf_fields_seq_insdel(df.filter(pl.col('vartype').is_in(['INS', 'DEL'])), ref_fa)
                .select('_index', '_vcf_pos', '_vcf_ref', '_vcf_alt', '_is_sym_alt')
            )

        if 'CPX' in vartype_set:
            df_pos_ref_alt_list.append(
                vcf_fields_seq_cpx(df.filter(pl.col('vartype') == 'CPX'), ref_fa, log_stats=log_stats)
                .select('_index', '_vcf_pos', '_vcf_ref', '_vcf_alt', '_is_sym_alt')
            )

        df_pos_ref_alt = (
            pl.concat(df_pos_ref_alt_list)
            if df_pos_ref_alt_list else
            vcf_fields_seq_insdel(df, ref_fa).select('_index', '_vcf_pos', '_vcf_ref', '_vcf_alt', '_is_sym_alt')
        )

    # Initialize fields
    df_vcf = (
        df
        .with_columns(
            '_index',
            pl.col('chrom').alias('_vcf_chrom'),
            pl.col('id').alias('_vcf_id'),
            pl.lit('.').alias('_vcf_qual'),
            (
                pl.when(pl.col('filter').list.len() > 0)
                .then(pl.col('filter').list.join(separator=';'))
                .otherwise(pl.lit('PASS'))
                .cast(pl.String)
                .alias('_vcf_filter')
            ),
            (
                pl.lit([])
                .cast(pl.List(pl.String))
                .alias('_vcf_info')
            ),
            (
                pl.lit([])
                .cast(pl.List(pl.String))
                .alias('_vcf_format')
            )
        )
    )

    # Add POS, REF, and ALT. Records describing no change are dropped by the field functions above (e.g. a complex
    # variant with a sequence matching the reference), so the join keeps only records with fields.
    return (
        df_vcf
        .drop(['_vcf_pos', '_vcf_ref', '_vcf_alt', '_is_sym_alt'], strict=False)
        .join(
            (
                df_pos_ref_alt
                .select('_index', '_vcf_pos', '_vcf_ref', '_vcf_alt', '_is_sym_alt')
            ), on='_index', how='inner', maintain_order='left'
        )
    )


def vcf_fields_seq_insdel(
    df: pl.LazyFrame,
    ref_fa: Path | str,
) -> pl.LazyFrame:
    """Set POS, REF, and ALT VCF fields for insertions and deletions (non-symbolic).

    VCF records where "REF" and "ALT" contain sequences (i.e. not symbolic ALTs like "<INS>")
    require a reference base immediately upstream or downstream of the variant to set the context.
    This function sets these fields and extracts the extra reference base.

    The reference base in field "_ref_base" is the base immediately before the variant (if variant
    position > 0) or immediately after (variant position == 0). Field "_ref_base_left" is a boolean
    flag indicating whether "_ref_base" is the base before the variant or after. In both cases,
    "_ref_base_pos" is the position of the variant in the reference base in reference
    coordinates (0-based).

    Returned table has fields:

    * _index: Index of variant.
    * _vcf_pos: Variant position for VCF field "POS".
    * _vcf_ref: Reference sequence for VCF field "REF".
    * _vcf_alt: Alternate sequence for VCF field "ALT".
    * _ref_base: Reference base.
    * _ref_base_pos: Position of variant in reference base in reference coordinates (0-based).
    * _ref_base_left: Whether the reference base is the base before the variant (True) or
      after (False). This is used to determine whether to put the reference base before or
      after the sequence when constructing the REF and ALT fields.
    * _is_sym_alt: Whether the alternate sequence is a symbolic alternate (always False).

    :param df: Variant table with an "_index" column.
    :param ref_fa: Reference FASTA filename.

    :returns: Variant table with reference base and position fields added.
    """

    if ref_fa is None:
        raise ValueError('Missing reference FASTA filename')

    ref_fa = str(ref_fa)

    # Get position of the reference base
    df_fields = (
        df.with_columns(
            (pl.col('pos') > 0).alias('_ref_base_left')
        )
        .with_columns(
            (
                pl.when('_ref_base_left')
                .then(pl.col('pos') - 1)
                .otherwise('end')
            ).alias('_ref_base_pos')
        )
    )

    # Carry an identifying column through for error reporting (dropped by the final select).
    label_col_list = ['id'] if 'id' in df.collect_schema().names() else []

    # Get fields
    with _LockingPysam(ref_fa) as ref_file:
        df_fields = (
            df_fields
            .select(
                '_index',
                (
                    pl.struct([pl.col('chrom'), pl.col('_ref_base_pos')])
                    .map_elements(
                        lambda vals: ref_file.fetch(
                            str(vals['chrom']), int(vals['_ref_base_pos']), int(vals['_ref_base_pos'] + 1)
                        ),
                        return_dtype=pl.String
                    )
                    .cast(pl.String)
                    .alias('_ref_base')
                ),
                '_ref_base_pos',
                '_ref_base_left',
                'pos',
                'vartype',
                'seq',
                *label_col_list,
            )
            .collect()  # Must collect while ref_file is open
        )

    # Set REF and ALT fields
    df_fields = (
        df_fields
        .with_columns(
            _concrete_bases(pl.col('_ref_base')),
            _concrete_bases(pl.col('seq')),
        )
        .with_columns(  # Set REF and ALT fields
            (
                pl.when(pl.col('vartype') == 'INS')
                .then(pl.lit(''))
                .when(pl.col('vartype') == 'DEL')
                .then(pl.col('seq'))
                .cast(pl.String)
                .alias('_vcf_ref')
            ),
            (
                pl.when(pl.col('vartype') == 'INS')
                .then(pl.col('seq'))
                .when(pl.col('vartype') == 'DEL')
                .then(pl.lit(''))
                .cast(pl.String)
                .alias('_vcf_alt')
            )
        )
        .select(
            '_index',
            (
                pl.when(pl.col('_ref_base_left'))
                .then(pl.col('pos'))
                .otherwise(pl.col('pos') + 1)
                .alias('_vcf_pos')
            ),
            (
                pl.when(pl.col('_ref_base_left'))
                .then(pl.concat_str('_ref_base', '_vcf_ref'))
                .otherwise(pl.concat_str('_vcf_ref', '_ref_base'))
                .cast(pl.String)
                .alias('_vcf_ref')
            ),
            (
                pl.when(pl.col('_ref_base_left'))
                .then(pl.concat_str('_ref_base', '_vcf_alt'))
                .otherwise(pl.concat_str('_vcf_alt', '_ref_base'))
                .cast(pl.String)
                .alias('_vcf_alt')
            ),
            pl.lit(False).alias('_is_sym_alt'),
            '_ref_base',
            '_ref_base_pos',
            '_ref_base_left',
            *label_col_list,
        )
    )

    _check_ref_alt(df_fields, label_col_list[0] if label_col_list else '_index')

    return df_fields.drop(label_col_list).lazy()


def vcf_fields_seq_cpx(
    df: pl.LazyFrame,
    ref_fa: Path | str,
    log_stats: Optional[dict] = None,
) -> pl.LazyFrame:
    """Set POS, REF, and ALT VCF fields for sequence-resolved complex variants.

    A complex variant replaces reference bases [seq_ref_pos, seq_ref_end) with its sequence ("seq"), and the
    sequence-resolved record describes exactly that, so applying it to the reference (e.g. "bcftools consensus")
    reconstructs the haplotype. The reference region is empty for a complex variant deleting no reference bases.

    As with insertions and deletions, REF and ALT start with the reference base before the replaced region, or end
    with the base after it when the region starts at the first base of the chromosome.

    Returned table has fields:

    * _index: Index of variant.
    * _vcf_pos: Variant position for VCF field "POS".
    * _vcf_ref: Reference sequence for VCF field "REF".
    * _vcf_alt: Alternate sequence for VCF field "ALT".
    * _is_sym_alt: Whether the alternate sequence is a symbolic alternate (always False).

    :param df: Variant table with an "_index" column and "seq_ref_pos", "seq_ref_end", and "seq" columns.
    :param ref_fa: Reference FASTA filename.
    :param log_stats: If not None, add statistics to this dictionary: "n_noop_cpx" (complex variants left out because
        their sequence matches the reference bases they replace) and "noop_cpx_label" (labels of up to 10 of them).

    :returns: Variant table with VCF fields.

    :raises ValueError: If columns are missing or the REF and ALT fields are invalid for any record.
    """
    if ref_fa is None:
        raise ValueError('Missing reference FASTA filename')

    if missing_cols := {'seq_ref_pos', 'seq_ref_end', 'seq'} - set(df.collect_schema().names()):
        raise ValueError(
            f'Missing columns for sequence-resolved complex variants: {", ".join(sorted(missing_cols))} '
            f'(complex variant tables from an earlier PAV version must be regenerated, or write symbolic complex '
            f'variants with config parameter "vcf_cpx_symbolic")'
        )

    label_col_list = ['id'] if 'id' in df.collect_schema().names() else []

    df_fields = df.with_columns(
        (pl.col('seq_ref_pos') > 0).alias('_ref_base_left')
    )

    # Fetch the replaced reference bases with the padding base (left if possible, right otherwise)
    with _LockingPysam(str(ref_fa)) as ref_file:
        df_fields = (
            df_fields
            .select(
                '_index',
                '_ref_base_left',
                'seq_ref_pos',
                'seq',
                (
                    pl.struct(['chrom', 'seq_ref_pos', 'seq_ref_end', '_ref_base_left'])
                    .map_elements(
                        lambda vals: ref_file.fetch(
                            str(vals['chrom']),
                            int(vals['seq_ref_pos'] - (1 if vals['_ref_base_left'] else 0)),
                            int(vals['seq_ref_end'] + (0 if vals['_ref_base_left'] else 1)),
                        ),
                        return_dtype=pl.String
                    )
                    .cast(pl.String)
                    .alias('_vcf_ref')
                ),
                *label_col_list,
            )
            .collect()  # Must collect while ref_file is open
        )

    df_fields = df_fields.with_columns(
        _concrete_bases(pl.col('_vcf_ref')),
        _concrete_bases(pl.col('seq')),
    ).select(
        '_index',
        (
            pl.when(pl.col('_ref_base_left'))
            .then(pl.col('seq_ref_pos'))
            .otherwise(pl.col('seq_ref_pos') + 1)
            .alias('_vcf_pos')
        ),
        '_vcf_ref',
        (
            pl.when(pl.col('_ref_base_left'))
            .then(pl.concat_str(pl.col('_vcf_ref').str.head(1), 'seq'))
            .otherwise(pl.concat_str('seq', pl.col('_vcf_ref').str.tail(1)))
            .cast(pl.String)
            .alias('_vcf_alt')
        ),
        pl.lit(False).alias('_is_sym_alt'),
        (pl.col('seq').str.len_chars() > 0).alias('_has_seq'),
        *label_col_list,
    )

    # A complex variant whose sequence matches the reference bases it replaces describes no change. It is left out of
    # the VCF with a warning rather than failing it. REF matching ALT without a sequence is still an error (see
    # `_check_ref_alt`).
    is_noop = pl.col('_has_seq').fill_null(False) & (pl.col('_vcf_ref') == pl.col('_vcf_alt')).fill_null(False)

    if (n_noop := df_fields.select(is_noop.sum()).item()) > 0:
        label_col = label_col_list[0] if label_col_list else '_index'
        label_list = [
            str(label) for label in df_fields.filter(is_noop)[label_col].head(_MAX_INVALID_REPORT).to_list()
        ]

        logger.warning(
            'Skipping %d complex variant(s) with a sequence matching the reference bases it replaces: %s%s',
            n_noop, ', '.join(label_list), ', ...' if n_noop > len(label_list) else ''
        )

        if log_stats is not None:
            log_stats['n_noop_cpx'] = n_noop
            log_stats['noop_cpx_label'] = label_list

        df_fields = df_fields.filter(~is_noop)

    df_fields = df_fields.drop('_has_seq')

    _check_ref_alt(df_fields, label_col_list[0] if label_col_list else '_index')

    return df_fields.drop(label_col_list).lazy()


def _check_ref_alt(
    df_fields: pl.DataFrame,
    label_col: str,
) -> None:
    """Reject structurally invalid sequence-resolved records.

    REF and ALT must both be non-empty, and a record with REF matching ALT describes no variant. Both occur when "seq"
    is null or empty for a variant type that needs it, which is a callset bug worth failing on rather than writing out
    for a downstream tool to trip over. REF and ALT are compared as written (upper case, see `_concrete_bases`), so a
    difference in case does not hide a match.

    :param df_fields: Table with "_vcf_ref" and "_vcf_alt" columns.
    :param label_col: Column identifying records in error messages.

    :raises ValueError: If any record is invalid.
    """
    df_invalid = df_fields.filter(
        pl.col('_vcf_ref').is_null()
        | pl.col('_vcf_alt').is_null()
        | (pl.col('_vcf_ref').str.len_chars() == 0)
        | (pl.col('_vcf_alt').str.len_chars() == 0)
        | (pl.col('_vcf_ref') == pl.col('_vcf_alt'))
    )

    if df_invalid.height > 0:
        label_list = [
            str(label) for label in df_invalid[label_col].head(_MAX_INVALID_REPORT).to_list()
        ]

        raise ValueError(
            f'Invalid REF or ALT in {df_invalid.height} VCF record(s) '
            f'(empty REF, empty ALT, or REF matching ALT): '
            f'{", ".join(label_list)}'
            f'{", ..." if df_invalid.height > len(label_list) else ""}'
        )


def vcf_fields_sym(
    df: pl.LazyFrame,
    ref_fa: Path | str,
) -> pl.LazyFrame:
    """Get POS, REF, and ALT fields for symbolic alternate records.

    Symbolic alts such as "<INS>" and "<DEL>" can be set for the ALT VCF field instead of exact
    sequences. In these cases, the SEQ info field should also be included (not handled by this
    function).

    * _index: Index of variant.
    * _vcf_pos: Variant position for VCF field "POS".
    * _vcf_ref: Reference sequence for VCF field "REF".
    * _vcf_alt: Symbolic alternate representation.
    * _is_sym_alt: Whether the alternate sequence is a symbolic alternate (always True).

    :param df: Variant table with an "_index" column.
    :param ref_fa: Reference FASTA filename.

    :returns: LazyFrame with VCF fields set.
    """
    subtype_expr = (
        alt_subtype_expr()
        if 'varsubtype' in df.collect_schema().names() else
        pl.lit(None).cast(pl.String)
    )

    # POS is the padding base before the variant (0-based pos - 1, or 1-based pos), and REF is that base. A variant
    # starting at the first base of the chromosome has no base before it and takes the first base.
    with _LockingPysam(ref_fa) as ref_file:
        return (
            df
            .with_columns(
                pl.max_horizontal(pl.col('pos'), pl.lit(1)).alias('_vcf_pos')
            )
            .select(
                '_index',
                '_vcf_pos',
                _concrete_bases(
                    pl.struct([pl.col('chrom'), pl.col('_vcf_pos')])
                    .map_elements(
                        lambda vals: ref_file.fetch(
                            str(vals['chrom']), int(vals['_vcf_pos'] - 1), int(vals['_vcf_pos'])
                        ),
                        return_dtype=pl.String
                    )
                    .cast(pl.String)
                ).alias('_vcf_ref'),
                (
                    pl.concat_str(
                        pl.lit('<'),
                        pl.col('vartype').str.to_uppercase(),
                        pl.when(subtype_expr.is_not_null())
                        .then(pl.concat_str(pl.lit(':'), subtype_expr))
                        .otherwise(pl.lit('')),
                        pl.lit('>')
                    )
                    .alias('_vcf_alt')
                ),
                pl.lit(True).alias('_is_sym_alt'),
            )
            .collect().lazy()  # Must collect while ref_file is open
        )


def vcf_fields_snv(
        df: pl.LazyFrame,
) -> pl.LazyFrame:
    """Set REF and ALT for VCF fields for SNVs.

    VCF records where "REF" and "ALT" contain sequences (i.e. not symbolic ALTs like "<INS>")
    require a reference base immediately upstream or downstream of the variant to set the context.
    This function sets these fields and extracts the extra reference base.

    The reference base in field "_ref_base" is the base immediately before the variant (if variant
    position > 0) or immediately after (variant position == 0). Field "_ref_base_left" is a boolean
    flag indicating whether "_ref_base" is the base before the variant or after. In both cases,
    "_ref_base_pos" is the position of the variant in the reference base in reference
    coordinates (0-based).

    Returned table has fields:

    * _index: Index of variant.
    * _vcf_pos: Variant position for VCF field "POS".
    * _vcf_ref: Reference sequence for VCF field "REF".
    * _vcf_alt: Alternate sequence for VCF field "ALT".
    * _is_sym_alt: Whether the alternate sequence is a symbolic alternate (always False).

    :param df: Variant table with an "_index" column.

    :returns: Variant table with reference base and position fields added.
    """

    return (
        df
        .select(
            '_index',
            (pl.col('pos') + 1).alias('_vcf_pos'),
            _concrete_bases(pl.col('ref')).alias('_vcf_ref'),
            _concrete_bases(pl.col('alt')).alias('_vcf_alt'),
            pl.lit(False).alias('_is_sym_alt'),
        )
        # An SNV on an ambiguous reference base (e.g. R) matching the base it reduces to (A) describes no change
        .filter(pl.col('_vcf_ref') != pl.col('_vcf_alt'))
    )


def gt_column_asm(
    df: pl.LazyFrame,
    asm_name: str,
    callable_dict: dict[tuple[str, str], pl.LazyFrame],
    col_name: str = 'gt',
    separator: str = '|',
) -> pl.LazyFrame:
    """Get a table of genotypes for a single assembly.

    The genotype string contains "1" (alternate allele), "0" (reference allele), and "." (no call).
    The order of alleles is defined by keys in `callable_dict` for all records matching the
    assembly name.

    :param df: Variant table with an "_index" column.
    :param asm_name: Assembly name.
    :param callable_dict: A dictionary keys tuple(assembly name, haplotype name) and values
        containing a polars LazyFrame of callable regions.
    :param col_name: Name of the genotype column in the output table.
    :param separator: Genotype separator character ("|" for phased, "/" for unphased).

    :return: A polars LazyFrame with two columns, "_index" and the genotype string (column name is
        set by `col_name`).
    """
    if df.collect_schema().get('_index') is None:
        df = df.with_row_index('_index')

    # Start with an index, add genotype columns
    df_gt = df.select('_index')

    for (callable_asm_name, hap), df_callable in callable_dict.items():
        if callable_asm_name != asm_name:
            continue

        # Determine if each value matches this haplotype.
        callable_hap = (
            agglovar.bed.intersect.as_proportion(
                df,
                df_callable,
                'callable_hap'
            )
            .with_columns(
                pl.col('callable_hap') >= 0.5
            )
        )

        df_hap = (
            df.select(
                '_index', 'chrom', 'pos', 'end',
                (
                    # Match the source's sample and haplotype directly from the "mg_src" struct
                    # fields (see pav3.schema.VARIANT). These are resolved from merge metadata when
                    # the callset is built, which avoids splitting a "sample-hap" name on "-".
                    pl.col('mg_src')
                    .list.filter(
                        (pl.element().struct.field('sample') == asm_name)
                        & (pl.element().struct.field('hap') == hap)
                    ).list.len() > 0
                ).alias('in_hap'),
            )
            .join(
                callable_hap, on='_index', how='left'
            )
            .with_columns(
                pl.col('callable_hap').fill_null(False)
            )
            .select(
                '_index',
                pl.when(pl.col('in_hap'))
                .then(pl.lit('1'))
                .when(pl.col('callable_hap'))
                .then(pl.lit('0'))
                .otherwise(pl.lit('.'))
                .cast(pl.String)
                .alias(hap)
            )
        )

        df_gt = df_gt.join(df_hap, on='_index', how='left')

    # Join tables and get haplotypes
    return (
        df_gt
        .select(
            '_index',
            pl.concat_str(cs.exclude('_index').fill_null(pl.lit('.')), separator=separator).alias(col_name)
        )
    )


def standard_info_fields(
    df: pl.LazyFrame,
) -> pl.LazyFrame:
    """Set standard INFO columns.

    Appends to the "_vcf_info" list for known INFO annotations found in variant tables.

    :param df: Variant table.

    :returns: Variant table with INFO columns added.
    """

    cols = set(df.collect_schema())

    if '_vcf_info' not in cols:
        df = df.with_columns(
            pl.lit([])
            .cast(pl.List(pl.String))
            .alias('_vcf_info')
        )

    # if 'id' in cols:
    #     df = df.with_columns(
    #         pl.col('_vcf_info').list.concat(
    #             pl.concat_str(pl.lit('ID='), 'id')
    #         )
    #     )

    # SVTYPE is not emitted. VCF 4.5 deprecates it as redundant with ALT, and the type is recoverable from the
    # symbolic ALT (<INV>, <DUP>, <CPX>) or from the REF and ALT lengths for sequence-resolved INS, DEL, and SNV.

    if len(cols & {'varlen', '_is_sym_alt'}) == 2:
        # Symbolic records only. VCF defines SVLEN for symbolic structural variant alleles and directs writers to use
        # the missing value for all others, where REF and ALT give the length. Writing it for sequence-resolved
        # records puts "." in multi-allelic records combining an SNV with an indel, which some tools fail on.
        #
        # SVLEN is always positive. VCF 4.5 takes the absolute value of SVLEN and treats a negative value as
        # positive, so deletions are no longer negated.
        df = df.with_columns(
            pl.when(pl.col('_is_sym_alt'))
            .then(
                pl.col('_vcf_info').list.concat(
                    pl.concat_str(pl.lit('SVLEN='), 'varlen')
                )
            )
            .otherwise('_vcf_info')
            .alias('_vcf_info')
        )

    if len(cols & {'end', 'vartype', '_is_sym_alt'}) == 3:
        # END spans the affected reference bases, which htslib uses to compute "rlen" for indexing. Without it,
        # symbolic records index as 1 bp and region queries can miss them.
        #
        # Only symbolic records need it. VCF 4.4 and later make END a computed field, and for a sequence-resolved
        # record REF already gives the span, so writing it there is redundant.
        #
        # Taken from "end" rather than POS plus SVLEN: "pos" and "end" are 0-based half-open, so "end" is already
        # the 1-based inclusive VCF END, and for CPX "varlen" is a query-side length, not a reference span.
        df = df.with_columns(
            pl.when(pl.col('_is_sym_alt') & pl.col('vartype').is_in(REF_SPAN_VARTYPE))
            .then(
                pl.col('_vcf_info').list.concat(
                    pl.concat_str(pl.lit('END='), 'end')
                )
            )
            .otherwise('_vcf_info')
            .alias('_vcf_info')
        )

    if len(cols & {'vartype', '_is_sym_alt'}) == 2:
        # PAV calls deletions and duplications from non-reference adjacencies between aligned segments, never from
        # a measure of DNA abundance, so the claim is always adjacency ("J"). A deleted segment relocated
        # elsewhere in the genome would leave abundance unchanged, and PAV has no way to observe that.
        df = df.with_columns(
            pl.when(pl.col('_is_sym_alt') & pl.col('vartype').is_in(SVCLAIM_VARTYPE))
            .then(
                pl.col('_vcf_info').list.concat(pl.lit('SVCLAIM=J'))
            )
            .otherwise('_vcf_info')
            .alias('_vcf_info')
        )

    if len(cols & {'varsubtype', 'vartype'}) == 2:
        # Subtypes the ALT allele already carries (duplications) are not repeated here.
        df = df.with_columns(
            pl.when(
                (pl.col('varsubtype').str.len_chars() > 0).fill_null(False) & alt_subtype_expr().is_null()
            )
            .then(
                pl.col('_vcf_info').list.concat(
                    pl.concat_str(pl.lit('SUBTYPE='), 'varsubtype')
                )
            )
            .otherwise('_vcf_info')
            .alias('_vcf_info')
        )

    if 'trace_ref' in cols:
        df = df.with_columns(
            pl.when((pl.col('trace_ref').str.len_chars() > 0).fill_null(False))
            .then(
                pl.col('_vcf_info').list.concat(
                    pl.concat_str(
                        pl.lit('CPX_REF_TRACE='),
                        _encode_info_value(pl.col('trace_ref').str.replace_all(',', ''))
                    )
                )
            )
            .otherwise('_vcf_info')
            .alias('_vcf_info')
        )

    if 'trace_qry' in cols:
        df = df.with_columns(
            pl.when((pl.col('trace_qry').str.len_chars() > 0).fill_null(False))
            .then(
                pl.col('_vcf_info').list.concat(
                    pl.concat_str(
                        pl.lit('CPX_QRY_TRACE='),
                        _encode_info_value(pl.col('trace_qry').str.replace_all(',', ''))
                    )
                )
            )
            .otherwise('_vcf_info')
            .alias('_vcf_info')
        )

    if len(cols & {'qry_id', 'qry_pos', 'qry_end'}) == 3:
        # Variants derived from complex events at the outer breakpoint have no query region (see
        # `pav3.call.integrate.add_cpx_derived`).
        df = df.with_columns(
            pl.when(pl.all_horizontal(pl.col('qry_id', 'qry_pos', 'qry_end').is_not_null()))
            .then(
                pl.col('_vcf_info').list.concat(
                    pl.concat_str(
                        pl.lit('QUERY='),
                        _encode_info_value(pl.col('qry_id')),
                        pl.lit(':'),
                        pl.col('qry_pos') + 1,
                        pl.lit('-'),
                        pl.col('qry_end')
                    )
                )
            )
            .otherwise('_vcf_info')
            .alias('_vcf_info')
        )

    if 'qry_rev' in cols:
        df = df.with_columns(
            pl.col('_vcf_info').list.concat(
                pl.concat_str(
                    pl.lit('QUERYREV='),
                    pl.col('qry_rev').replace_strict({True: '-', False: '+'}, default='*')
                )
            )
        )

    if len(cols & {'_is_sym_alt', 'seq'}) == 2:
        # Symbolic records only. REF and ALT already carry the sequence for sequence-resolved records, where SEQ
        # duplicates it exactly. Upper-cased to match REF and ALT, which are not soft-masked.
        #
        # Derived records are left out. A derived duplication carries the reference over the duplicated region,
        # which does not describe an insertion placed inside the region it duplicates. POS, END, and SVLEN
        # describe the variant, and the sequence remains available in the variant tables.
        not_derived = (
            (pl.col('call_source') != 'DERIVED')
            if 'call_source' in cols else
            pl.lit(True)
        )

        df = df.with_columns(
            pl.when(pl.col('_is_sym_alt') & pl.col('seq').is_not_null() & not_derived)
            .then(
                pl.col('_vcf_info').list.concat(
                    pl.concat_str(pl.lit('SEQ='), pl.col('seq').str.to_uppercase())
                )
            )
            .otherwise('_vcf_info')
            .alias('_vcf_info')
        )

    if 'dup' in cols:
        df = df.with_columns(
            pl.when(pl.col('dup').is_not_null() & (pl.col('dup').list.len() > 0))
            .then(
                pl.col('_vcf_info').list.concat(
                    pl.concat_str(
                        pl.lit('DUP='),
                        pl.col('dup').list.eval(
                            pl.concat_str(
                                _encode_info_value(pl.element().struct.field('chrom')),
                                pl.lit(':'),
                                (pl.element().struct.field('pos') + 1).cast(pl.String),
                                pl.lit('-'),
                                pl.element().struct.field('end').cast(pl.String),
                            )
                        ).list.join('/')
                    )
                )
            )
            .otherwise('_vcf_info')
            .alias('_vcf_info')
        )

    if len(cols & {'hom_ref', 'hom_qry', 'vartype'}) == 3:
        # Breakpoint homology as a single length, summing the upstream and downstream sides.
        #
        # Reference breakpoints describe a deletion or duplication best and query breakpoints an insertion or
        # complex variant. SNVs have no homology, and PAV does not refine inversion breakpoints finely enough to
        # report it.
        hom_expr = (
            pl.when(pl.col('vartype').is_in({'DEL', 'DUP'}))
            .then(pl.col('hom_ref'))
            .when(pl.col('vartype').is_in({'INS', 'CPX'}))
            .then(pl.col('hom_qry'))
            .otherwise(None)
        )

        df = df.with_columns(
            pl.when(hom_expr.is_not_null())
            .then(
                pl.col('_vcf_info').list.concat(
                    pl.concat_str(
                        pl.lit('HOMLEN='),
                        hom_expr.struct.field('up') + hom_expr.struct.field('dn'),
                    )
                )
            )
            .otherwise('_vcf_info')
            .alias('_vcf_info')
        )

    # TODO: Adjust inner and derived annotations for merged variants.
    #
    # Merging can get variant IDs out of sync with annotated derived and inner IDs. For example, if Variant A is
    # inside variant B and the representations for variants A and B come from different haplotypes, then this annotation
    # may not make sense. For example, if Variant A is derived from haplotype "h2" and variant B derived from haplotype
    # "h1" (e.g. variant B is homozygous and A is not), then the variant IDs will be out of sync. Unless merging
    # corrects this, do not report in the VCF.
    #
    # if 'derived' in cols:
    #     df = df.with_columns(
    #         pl.when(pl.col('derived').is_not_null())
    #         .then(
    #             pl.col('_vcf_info').list.concat(
    #                 pl.concat_str(pl.lit('DERIVED='), 'derived')
    #             )
    #         )
    #         .otherwise('_vcf_info')
    #         .alias('_vcf_info')
    #     )
    #
    # if 'inner' in cols:
    #     df = df.with_columns(
    #         pl.when(pl.col('inner').list.len() > 0)
    #         .then(
    #             pl.col('_vcf_info').list.concat(
    #                 pl.concat_str(
    #                     pl.lit('INNER='),
    #                     pl.col('inner').list.join(','),
    #                 )
    #             )
    #         )
    #         .otherwise('_vcf_info')
    #         .alias('_vcf_info')
    #     )
    #
    # if 'discord' in cols:
    #     df = df.with_columns(
    #         pl.when(pl.col('discord').list.len() > 0)
    #         .then(
    #             pl.col('_vcf_info').list.concat(
    #                 pl.concat_str(
    #                     pl.lit('DISCORD='),
    #                     pl.col('inner').list.join(','),
    #                 )
    #             )
    #         )
    #         .otherwise('_vcf_info')
    #         .alias('_vcf_info')
    #     )

    # A field built from a missing value is null (e.g. HOMLEN with a missing breakpoint side). Leave it out of INFO,
    # which the VCF writer did implicitly, so records combined by `combine_multiallelic` never see null fields.
    return df.with_columns(
        pl.col('_vcf_info').list.drop_nulls()
    )


VCF_RECORD_COLUMNS = (
    '_vcf_chrom', '_vcf_pos', '_vcf_id', '_vcf_ref', '_vcf_alt', '_vcf_filter', '_vcf_info', '_vcf_format',
    '_is_sym_alt',
)
"""Columns of a VCF record table before it is formatted (see `vcf_record_table`), excluding sample columns."""


def vcf_record_table(
    df: pl.LazyFrame,
    sample_col_prefix: str = '_vcf_sample_',
) -> pl.LazyFrame:
    """Select VCF record fields from a variant table with VCF fields set.

    The table has one row per VCF record with fields that are not yet formatted (INFO, FORMAT, and sample fields are
    lists), so records from different variant types may be concatenated and combined (see `combine_multiallelic`)
    before they are formatted by `reformat_vcf_table`.

    :param df: Variant table with VCF fields ("_vcf_" columns, see `init_vcf_fields`).
    :param sample_col_prefix: Prefix of sample columns.

    :returns: Table of VCF records.
    """
    return df.select(
        *VCF_RECORD_COLUMNS,
        cs.starts_with(sample_col_prefix),
    )


def combine_multiallelic(
    df: pl.LazyFrame | pl.DataFrame,
    sample_col_prefix: str = '_vcf_sample_',
    log_stats: Optional[dict[str, int]] = None,
) -> pl.LazyFrame:
    """Combine sequence-resolved records at the same position into multi-allelic records.

    Haplotypes carrying different alleles at one site are written as separate records by default, one per allele
    (e.g. "1|0" and "0|1"). Comparison tools (e.g. vcfeval) and "bcftools consensus" do not reliably resolve
    overlapping records, so records at the same position are combined into one record with an ALT per distinct allele
    (e.g. "1|2").

    Records are combined if they are sequence-resolved, have the same CHROM, POS, and FILTER, and if no haplotype
    carries more than one of them. Records with different filters are not combined so that a passing allele is never
    filtered by an allele on another haplotype. If a haplotype carries more than one record at the same position
    (e.g. an SNV and an insertion after it), the records are left as they are. Records are also left as they are if
    extending a haplotype's allele to the combined REF would overlap another variant on that haplotype (e.g. an SNV
    combined with a deletion on the other haplotype, and a second SNV inside the deletion on the SNV's haplotype),
    which would make the two variants on one haplotype overlap. Overlapping records starting at different positions
    are not combined.

    In a combined record:

    * REF: The longest REF. All records start at the same position, so shorter REF alleles are a prefix of it.
    * ALT: One per distinct allele, each extended by the reference bases the longer REF adds, in order of haplotypes
      carrying them. Alleles no haplotype carries (e.g. haplotypes excluded from the VCF) follow.
    * ID: IDs separated by semicolons in ALT order.
    * INFO: Per-allele (Number=A) values in ALT order. Values from the first record are used when two records have
      the same ALT allele. A value missing from a record is ".". Fields that are not per-allele are dropped.
    * GT: The ALT allele index for haplotypes carrying an allele. Haplotypes carrying no allele are "0" unless no
      record is callable on that haplotype ("."), for which every record has ".".

    :param df: Table of VCF records (see `vcf_record_table`) with genotypes as the first element of each sample
        column (alleles separated by "|" or "/").
    :param sample_col_prefix: Prefix of sample columns.
    :param log_stats: If not None, add statistics to this dictionary: "n_record_in", "n_record_out",
        "n_combined_record" (records combined into multi-allelic records), "n_combined_site" (multi-allelic records
        written), "n_conflict_site" (sites with records that could not be combined because a haplotype carries more
        than one), "n_span_conflict_site" (sites not combined because an extended allele would overlap another variant
        on its haplotype), and "n_overlap_record" (records overlapping the REF of a preceding record at a different
        position).

    :returns: Table of VCF records (LazyFrame evaluating `df` again when collected).
    """
    if isinstance(df, pl.DataFrame):
        df = df.lazy()

    sample_cols = [col for col in df.collect_schema().names() if col.startswith(sample_col_prefix)]

    info_number = {info_field.name: info_field.number for info_field in INFO_FIELDS}

    group_cols = ['_vcf_chrom', '_vcf_pos', '_vcf_filter']

    # Flag records sharing a site. Only these are collected, other records pass through.
    df = (
        df.with_row_index('_rec_index')
        .with_columns(
            (
                (~ pl.col('_is_sym_alt'))
                & (pl.col('_vcf_pos').count().over([*group_cols, '_is_sym_alt']) > 1)
            ).alias('_is_multi')
        )
    )

    df_schema = df.collect_schema()

    df_multi = df.filter(pl.col('_is_multi')).collect().sort('_rec_index')

    stats = collections.Counter({
        'n_record_in': df.select(pl.len()).collect().item(),
        'n_combined_record': 0,
        'n_combined_site': 0,
        'n_conflict_site': 0,
        'n_span_conflict_site': 0,
    })

    hap_var_index = _HapVariantIndex(df, sample_cols) if df_multi.height > 0 else None

    out_rows = []

    for _, df_group in df_multi.group_by(group_cols, maintain_order=True):
        rows = df_group.rows(named=True)

        combined = _combine_site(rows, sample_cols, info_number, hap_var_index)

        if combined is None:
            stats['n_conflict_site'] += 1
            out_rows.extend(rows)
        elif combined == 'span':
            stats['n_span_conflict_site'] += 1
            out_rows.extend(rows)
        else:
            stats['n_combined_record'] += len(rows)
            stats['n_combined_site'] += 1
            out_rows.append(combined)

    df_out = (
        pl.concat(
            [
                df.filter(~ pl.col('_is_multi')),
                pl.DataFrame(out_rows, schema=df_schema).lazy(),
            ]
        )
        .sort('_vcf_chrom', '_vcf_pos', '_rec_index')
        .drop('_rec_index', '_is_multi')
    )

    if log_stats is not None:
        stats['n_record_out'] = stats['n_record_in'] - stats['n_combined_record'] + stats['n_combined_site']

        # Count sequence-resolved records overlapping the REF of a preceding record at a different position. These
        # are not combined.
        stats['n_overlap_record'] = (
            df_out
            .filter(~ pl.col('_is_sym_alt'))
            .select(
                '_vcf_chrom', '_vcf_pos',
                (pl.col('_vcf_pos') + pl.col('_vcf_ref').str.len_chars()).alias('_ref_end')
            )
            .with_columns(
                pl.col('_ref_end').cum_max().shift(1).over('_vcf_chrom').alias('_prev_end'),
                pl.col('_vcf_pos').shift(1).over('_vcf_chrom').alias('_prev_pos'),
            )
            .filter(
                (pl.col('_vcf_pos') < pl.col('_prev_end'))
                & (pl.col('_vcf_pos') != pl.col('_prev_pos'))
            )
            .select(pl.len())
            .collect()
            .item()
        )

        log_stats.update(stats)

    return df_out


class _HapVariantIndex:
    """Positions of records carried by each haplotype, for finding variants on one haplotype in a region."""

    def __init__(
        self,
        df: pl.LazyFrame,
        sample_cols: list[str],
    ) -> None:
        """Index positions of all records carried by each haplotype.

        :param df: Table of VCF records with genotypes as the first element of each sample column.
        :param sample_cols: Sample columns.
        """
        df_gt = (
            df.select(
                '_vcf_chrom', '_vcf_pos',
                *[
                    pl.col(sample_col).list.first().str.replace_all('/', '|').str.split('|').alias(sample_col)
                    for sample_col in sample_cols
                ]
            )
            .collect()
        )

        self.pos_dict: dict[tuple[str, int, str], np.ndarray] = {}

        for sample_col in sample_cols:
            n_hap = df_gt.select(pl.col(sample_col).list.len().max()).item() or 0

            for hap_i in range(n_hap):
                for chrom, pos_list in (
                    df_gt
                    .filter(pl.col(sample_col).list.get(hap_i, null_on_oob=True) == '1')
                    .group_by('_vcf_chrom')
                    .agg(pl.col('_vcf_pos').sort())
                    .iter_rows()
                ):
                    self.pos_dict[(sample_col, hap_i, chrom)] = np.asarray(pos_list)

    def has_variant(
        self,
        sample_col: str,
        hap_i: int,
        chrom: str,
        start: int,
        end: int,
    ) -> bool:
        """Determine if a haplotype carries a record with POS in [start, end)."""
        pos_arr = self.pos_dict.get((sample_col, hap_i, chrom), None)

        if pos_arr is None:
            return False

        i = np.searchsorted(pos_arr, start, side='left')

        return bool(i < len(pos_arr) and pos_arr[i] < end)


def _combine_site(
    rows: list[dict],
    sample_cols: list[str],
    info_number: dict[str, str],
    hap_var_index: Optional[_HapVariantIndex] = None,
) -> Optional[dict | str]:
    """Combine records at one site into a multi-allelic record.

    :param rows: Records at the same site (same CHROM, POS, and FILTER), in output order.
    :param sample_cols: Sample columns.
    :param info_number: INFO "Number" for each INFO field.
    :param hap_var_index: Positions of records carried by each haplotype. If None, extended alleles are not checked
        for overlapping variants on the same haplotype.

    :returns: A combined record, "span" if an extended allele would overlap another variant on its haplotype, or None
        if the records cannot be combined for another reason.
    """
    # All records must have the same FORMAT with GT first, and REF alleles must share a prefix.
    if any(row['_vcf_format'] != rows[0]['_vcf_format'] for row in rows) or rows[0]['_vcf_format'] != ['GT']:
        return None

    ref = max((row['_vcf_ref'].upper() for row in rows), key=len)

    if not all(ref.startswith(row['_vcf_ref'].upper()) for row in rows):
        return None

    # Genotypes per record and sample
    rec_gt = []  # [record][sample] -> list of allele strings
    sample_sep = {}

    for row in rows:
        gt_list = []

        for sample_col in sample_cols:
            gt = row[sample_col][0]
            sep = '/' if '/' in gt else '|'
            sample_sep.setdefault(sample_col, sep)
            gt_list.append(gt.split(sep))

        rec_gt.append(gt_list)

    # Each haplotype may carry at most one record
    n_hap = {sample_col: len(rec_gt[0][i]) for i, sample_col in enumerate(sample_cols)}

    for i, sample_col in enumerate(sample_cols):
        if any(len(gt[i]) != n_hap[sample_col] for gt in rec_gt):
            return None

        for hap_i in range(n_hap[sample_col]):
            if sum(gt[i][hap_i] == '1' for gt in rec_gt) > 1:
                return None

    # Extending an allele to the combined REF must not overlap another variant on the haplotype carrying it
    if hap_var_index is not None:
        site_pos = rows[0]['_vcf_pos']
        site_chrom = rows[0]['_vcf_chrom']

        for i, sample_col in enumerate(sample_cols):
            for hap_i in range(n_hap[sample_col]):
                for rec_i, gt in enumerate(rec_gt):
                    ref_len = len(rows[rec_i]['_vcf_ref'])

                    if gt[i][hap_i] == '1' and ref_len < len(ref) and hap_var_index.has_variant(
                        sample_col, hap_i, site_chrom, site_pos + ref_len, site_pos + len(ref)
                    ):
                        return 'span'

    # Extended ALT allele for each record
    rec_alt = [row['_vcf_alt'].upper() + ref[len(row['_vcf_ref']):] for row in rows]

    # Order alleles by haplotypes carrying them, then alleles not carried by any haplotype
    rec_order = []

    for i, sample_col in enumerate(sample_cols):
        for hap_i in range(n_hap[sample_col]):
            for rec_i, gt in enumerate(rec_gt):
                if gt[i][hap_i] == '1' and rec_i not in rec_order:
                    rec_order.append(rec_i)

    rec_order.extend(rec_i for rec_i in range(len(rows)) if rec_i not in rec_order)

    alt_list = []        # Distinct ALT alleles in order
    alt_rec = []         # Record providing fields for each ALT allele
    rec_allele = {}      # Record index -> 1-based allele index

    for rec_i in rec_order:
        if rec_alt[rec_i] in alt_list:
            rec_allele[rec_i] = alt_list.index(rec_alt[rec_i]) + 1
        else:
            alt_list.append(rec_alt[rec_i])
            alt_rec.append(rec_i)
            rec_allele[rec_i] = len(alt_list)

    # ID
    id_list = []

    for rec_i in rec_order:
        rec_id = rows[rec_i]['_vcf_id']

        if rec_id is not None and rec_id != '.' and rec_id not in id_list:
            id_list.append(rec_id)

    # INFO: Per-allele fields in ALT order
    rec_info = []

    for rec_i in alt_rec:
        info_dict = {}

        for info in rows[rec_i]['_vcf_info']:
            key, val = info.split('=', 1) if '=' in info else (info, None)
            info_dict[key] = val

        rec_info.append(info_dict)

    info_keys = list(dict.fromkeys(key for info_dict in rec_info for key in info_dict.keys()))

    info_list = [
        f'{key}=' + ','.join(
            info_dict.get(key) if info_dict.get(key) is not None else '.' for info_dict in rec_info
        )
        for key in info_keys if info_number.get(key) == 'A'
    ]

    # Genotypes
    combined = dict(rows[0])

    for i, sample_col in enumerate(sample_cols):
        allele_list = []

        for hap_i in range(n_hap[sample_col]):
            carrier = [rec_i for rec_i, gt in enumerate(rec_gt) if gt[i][hap_i] == '1']

            if carrier:
                allele_list.append(str(rec_allele[carrier[0]]))
            elif all(gt[i][hap_i] == '.' for gt in rec_gt):
                allele_list.append('.')
            else:
                allele_list.append('0')

        combined[sample_col] = [sample_sep[sample_col].join(allele_list)] + list(rows[0][sample_col][1:])

    combined['_vcf_ref'] = ref
    combined['_vcf_alt'] = ','.join(alt_list)
    combined['_vcf_id'] = ';'.join(id_list) if id_list else None
    combined['_vcf_info'] = info_list

    return combined


def reformat_vcf_table(
    df: pl.LazyFrame,
    sample_col_prefix: str = '_vcf_sample_'
) -> pl.LazyFrame:
    """Reformat VCF table to VCF format."""
    schema = df.collect_schema()

    strip_len = len(sample_col_prefix)

    sample_col_exprs = [
        pl.col(sample_col).list.join(':').alias(sample_col[strip_len:])
        for sample_col in schema.keys() if sample_col.startswith(sample_col_prefix)
    ]

    return (
        df
        .select(
            pl.col('_vcf_chrom').alias('#CHROM'),
            pl.col('_vcf_pos').alias('POS'),
            pl.col('_vcf_id').fill_null('.').alias('ID'),
            pl.col('_vcf_ref').alias('REF').str.to_uppercase(),
            pl.col('_vcf_alt').alias('ALT').str.to_uppercase(),
            pl.lit('.').alias('QUAL'),
            pl.col('_vcf_filter').alias('FILTER'),
            (
                pl.when(pl.col('_vcf_info').list.len() > 0)
                .then(
                    pl.col('_vcf_info')
                    .list.eval(  # Values are encoded where they are built, but don't break VCFs if one is missed
                        pl.element().str.replace_all(';', '%3B', literal=True)
                    )
                    .list.join(';')
                )
                .otherwise(pl.lit('.'))
            ).alias('INFO'),
            pl.col('_vcf_format').list.join(':').alias('FORMAT'),
            *sample_col_exprs,
        )
    )


def sort_vcf_table(
    df: pl.LazyFrame,
    df_ref_info: pl.DataFrame,
) -> pl.LazyFrame:
    """Sort formatted VCF records in reference order.

    Chromosomes are ordered as they are in the reference FASTA, which is the order of "##contig" header lines (see
    `get_headers`). VCF only requires the records for each chromosome to be contiguous, but tools validating records
    against the sequence dictionary (e.g. Picard) reject other orders. Records are sorted by position within each
    chromosome and by ID at the same position.

    :param df: Formatted VCF records (see `reformat_vcf_table`).
    :param df_ref_info: Reference information table with a "chrom" column and optionally an "order" column.

    :returns: Sorted VCF records.
    """
    chrom_order = (
        _ref_info_order(df_ref_info)
        .select(
            pl.col('chrom').alias('#CHROM'),
            pl.int_range(pl.len(), dtype=pl.Int64).alias('_chrom_order'),
        )
        .lazy()
    )

    return (
        df
        .join(chrom_order, on='#CHROM', how='left', maintain_order='left')
        .sort('_chrom_order', '#CHROM', 'POS', 'ID', nulls_last=True)
        .drop('_chrom_order')
    )


def resolve_hap_list(
    hap_list: Iterable[str],
    vcf_haplotypes: Optional[str],
) -> list[str]:
    """Get a list of haplotypes in the order they should appear for a sample.

    :param hap_list: List of haplotypes defined for this assembly.
    :param vcf_haplotypes: Optional, a string of haplotypes separated by commas to include in the
        VCF haplotypes. If defined, all haplotypes must appear in `hap_list` and the returned list
        is subset and ordered by this list.

    :returns: List of haplotypes.
    """
    hap_list = list(hap_list)

    if vcf_haplotypes is not None:
        vcf_hap_list = [hap.strip() for hap in vcf_haplotypes.split(',') if hap.strip()]

        if missing_hap := [hap for hap in vcf_hap_list if hap not in set(hap_list)]:
            missing_n = len(missing_hap)

            missing_hap = ', '.join([
                f'"{hap}"' for hap in missing_hap
            ]) + ('...' if missing_n > 3 else '')

            raise ValueError(
                f'Configuration item "vcf_haplotypes" defines {missing_n} haplotypes: '
                f'{missing_hap}'
            )

        return vcf_hap_list
    return hap_list


class _LockingPysam():
    """Serves pysam.FastaFile.fetch() under a lock

    Pysam is not thread-safe, and will fail when called through Polars. This class explicitly locks it.
    """
    filename: Path

    def __init__(
            self,
            filename: Path | str,
    ):
        """Create a new locking object.

        :param filename: FASTA file name.
        """
        self.filename = Path(filename)
        self.ref_file = None
        self.lock = threading.Lock()

    def __enter__(self):
        if self.ref_file is not None:
            raise RuntimeError('Ref file already open')

        self.ref_file = pysam.FastaFile(str(self.filename))

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.ref_file is None:
            raise RuntimeError('Ref file not open')

        self.ref_file.close()
        self.ref_file = None

    def fetch(
            self,
            chrom: str,
            start: int,
            end: int,
    ):
        if self.ref_file is None:
            raise RuntimeError('Ref file not open')

        with self.lock:
            return self.ref_file.fetch(chrom, start, end)
