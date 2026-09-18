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
    'alt_subtype_expr',
    'get_headers',
    'init_vcf_fields',
    'vcf_fields_seq_insdel',
    'vcf_fields_sym',
    'vcf_fields_snv',
    'gt_column_iterator',
    'standard_info_fields',
    'reformat_vcf_table',
]


import collections
from collections.abc import Iterable
import datetime
from pathlib import Path
import threading
from typing import Optional

import agglovar
import polars as pl
import polars.selectors as cs
import pysam

from . import __version__
from .const import FILTER_REASON

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
        'Variant length: inserted bases for INS and CPX (query-side), affected reference span for DEL, INV, '
        'and DUP. Does not equal the POS/END reference span for INS or CPX.',
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
        'Variant subtype. TD for a tandem duplication, or for CPX the type of the reference breakpoint (DEL or DUP)',
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
        'Query sequence position',
    ),
    InfoField(
        'QUERYREV', 'A', 'Character',
        'Query sequence orientation vs the reference (+ or -, * if unknown)',
    ),
    InfoField(
        'SEQ', 'A', 'String',
        'Variant sequence.',
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


REF_SPAN_VARTYPE = frozenset({'DEL', 'INV', 'DUP', 'CPX'})
"""Variant types spanning reference bases, which receive an INFO/END annotation.

INS spans no reference bases and SNVs are fully described by REF and ALT, so neither is annotated.
"""

SVCLAIM_VARTYPE = frozenset({'DEL', 'DUP'})
"""Variant types requiring an INFO/SVCLAIM annotation on a symbolic ALT allele.

VCF requires a claim for symbolic DEL and DUP alleles. INS and INV carry an implicit adjacency claim, and a
sequence-resolved REF and ALT pair leaves nothing to disambiguate.
"""


def get_headers(
    ref_filename: Optional[str | Path] = None,
    df_ref_info: Optional[pl.DataFrame] = None,
) -> list[str]:

    # Lead headers
    headers = [
        f'##fileformat=VCFv{VCF_VERSION}',
        f'##filedate={datetime.datetime.now().strftime("%Y%m%d")}',
        f'##source=PAVv{__version__}',
    ]

    # Reference
    if ref_filename is not None:
        # Only filename, do not leak filesystem paths for security reasons
        headers.append(
            f'##reference=file://{Path(ref_filename).name}'
        )

    # Reference contigs
    if df_ref_info is not None:
        headers.extend(
            f'##contig=<ID={chrom},length={length},md5={md5}>'
            for chrom, length, md5 in df_ref_info.select(['chrom', 'len', 'md5']).iter_rows()
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
) -> pl.LazyFrame:
    """Initialize VCF fields.

    :param df: Variant table with an "_index" column.
    :param ref_fa: Reference FASTA filename. Required for non-SNV variants.
    :param use_sym: Whether to use symbolic alternate records. If None, INS, DEL, and SNV ore not
        symbolic, and all other variant types are.
    :param vartype: Expected variant type or types (i.e. "insdel" allows "INS" and "DEL"). Setting
        this field ensures the variant types is handled consisntently if all variant types
        are not in table `df` or `df` is empty.
    """
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

        no_sym_set = vartype_set & {'INS', 'DEL', 'SNV'}
        is_sym_set = vartype_set - {'INS', 'DEL', 'SNV'}

        if not (no_sym_set or is_sym_set):
            use_sym = False

        elif no_sym_set and is_sym_set:
            raise ValueError(
                f'Unable to automatically determine whether to use symbolic ALTs or sequence ALTs given variant types '
                f'using sequence ALTs (INS, DEL, SNV) and symbolic ALTs (all others): '
                f'Found variant types {", ".join(vartype_set)}'
            )
        else:
            use_sym = bool(is_sym_set)

    if use_sym and is_snv:
        raise ValueError('Variant type "SNV" cannot use symbolic ALTs')

    if not use_sym and (vartype_set - {'INS', 'DEL', 'SNV'}):
        raise ValueError(
            f'Can only use sequence ALTs for INS, DEL, and SNV: '
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

        df_pos_ref_alt = vcf_fields_seq_insdel(df, ref_fa)

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

    # Add POS, REF, and ALT
    return (
        df_vcf
        .drop(['_vcf_pos', '_vcf_ref', '_vcf_alt', '_is_sym_alt'], strict=False)
        .join(
            (
                df_pos_ref_alt
                .select('_index', '_vcf_pos', '_vcf_ref', '_vcf_alt', '_is_sym_alt')
            ), on='_index', how='left'
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

    # Reject structurally invalid records.
    #
    # REF and ALT must both be non-empty, and a record with REF matching ALT describes no variant. Both occur when
    # "seq" is null or empty for a variant type that needs it, which is a callset bug worth failing on rather than
    # writing out for a downstream tool to trip over.
    df_invalid = df_fields.filter(
        pl.col('_vcf_ref').is_null()
        | pl.col('_vcf_alt').is_null()
        | (pl.col('_vcf_ref').str.len_chars() == 0)
        | (pl.col('_vcf_alt').str.len_chars() == 0)
        | (pl.col('_vcf_ref') == pl.col('_vcf_alt'))
    )

    if df_invalid.height > 0:
        label_col = label_col_list[0] if label_col_list else '_index'

        label_list = [
            str(label) for label in df_invalid[label_col].head(_MAX_INVALID_REPORT).to_list()
        ]

        raise ValueError(
            f'Invalid REF or ALT in {df_invalid.height} VCF record(s) '
            f'(empty REF, empty ALT, or REF matching ALT): '
            f'{", ".join(label_list)}'
            f'{", ..." if df_invalid.height > len(label_list) else ""}'
        )

    return df_fields.drop(label_col_list).lazy()


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

    with _LockingPysam(ref_fa) as ref_file:
        return (
            df
            .select(
                '_index',
                pl.col('pos').alias('_vcf_pos'),
                (
                    pl.struct([pl.col('chrom'), pl.col('pos')])
                    .map_elements(
                        lambda vals: ref_file.fetch(str(vals['chrom']), int(vals['pos']), int(vals['pos'] + 1)),
                        return_dtype=pl.String
                    )
                    .cast(pl.String)
                    .alias('_vcf_ref')
                ),
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
            pl.col('ref').alias('_vcf_ref'),
            pl.col('alt').alias('_vcf_alt'),
            pl.lit(False).alias('_is_sym_alt'),
        )
    )


def gt_column_iterator(
    df: pl.LazyFrame | pl.DataFrame,
    hap_order: Iterable[tuple[str, str]],
    callable_dict: dict[tuple[str, str], pl.LazyFrame | pl.DataFrame],
    gt_col_name: str = 'gt',
) -> Iterable[tuple[str, pl.LazyFrame]]:
    """Get genotype columns for all assemblies/samples.

    Argument `hap_order` is a list of assembly/haplotype pairs, and it determines the order of
    assemblies (or "samples") and the order of haplotypes for each sample. Generally, it should
    list all the haplotypes for a sample sequentially (e.g.
    `[("a", "h1"), ("a", "h2"), ("b", "h1"), ("b", "h2")]`). If assemblies are not listed
    sequentially, then the first accurance of each assembly name determines the order, but
    haplotype order for each assembly is still preserved.

    Returns an iterator of assembly name and genotype table. The genotype table contains two
    columns: Column "_index" is the index of `df`, and the second column is string genotypes for
    each record with a column name defined by argument `gt_col_name`.

    :param df: Table of merged variants.
    :param hap_order: An iterable of assembly/haplotype tuples in the order they will appear.
    :param callable_dict: A dict keyed by assembly name and haplotype (tuple) pointing to tables
        of callable reference regions.
    :param gt_col_name: Name of the genotype column in the returned genotype tables. May not be
        "_index".

    :returns: An iterator of tuples (assembly name, and genotype table).
    """
    if gt_col_name is None or not (gt_col_name := gt_col_name.strip()):
        gt_col_name = 'gt'

    if gt_col_name == '_index':
        raise ValueError(
            'Genotype column argument "gt_col_name" cannot use reserved '
            'column name for row indexes, "_index"'
        )

    # Set and check assembly and haplotype orders
    hap_order = list(hap_order)

    if dup_asm_hap := [
        asm_hap for asm_hap, count in collections.Counter(f'{asm_name}-{hap}' for asm_name, hap in hap_order).items()
        if count > 1
    ]:
        n = len(dup_asm_hap)
        dup_asm_hap = ', '.join(dup_asm_hap[:3]) + ('...' if n > 3 else '')
        raise ValueError(
            f'Assembly/haplotype pairs must be unique when separated by a dash: '
            f'Found {n} duplicate pairs: {dup_asm_hap}'
        )

    asm_order = list(dict.fromkeys((asm_name for asm_name, hap in hap_order)).keys())

    # List of ordered haplotypes per assembly
    hap_dict = {
        asm_name_ordered: list(dict.fromkeys((
            hap for asm_name, hap in hap_order if asm_name == asm_name_ordered
        )).keys())
        for asm_name_ordered in asm_order
    }

    # Callable dict should have LazyFrame
    callable_dict = dict(callable_dict)

    for (asm_name, hap), df_callable in callable_dict.values():
        if isinstance(df_callable, pl.DataFrame):
            df_callable[(asm_name, hap)] = df_callable.lazy()

    for asm_name in hap_dict.keys():
        yield asm_name, gt_column_asm(
            df,
            asm_name=asm_name,
            callable_dict=callable_dict,
            col_name='gt',
            separator='|',
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

    if len(cols & {'varlen', 'vartype'}) == 2:
        # SVLEN is always positive. VCF 4.5 takes the absolute value of SVLEN and treats a negative value as
        # positive, so deletions are no longer negated.
        df = df.with_columns(
            pl.when(pl.col('vartype') != 'SNV')
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
            pl.when(pl.col('varsubtype').is_not_null() & alt_subtype_expr().is_null())
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
            pl.when(pl.col('trace_ref').is_not_null())
            .then(
                pl.col('_vcf_info').list.concat(
                    pl.concat_str(
                        pl.lit('CPX_REF_TRACE='),
                        pl.col('trace_ref').str.replace_all(',', '')
                    )
                )
            )
            .otherwise('_vcf_info')
            .alias('_vcf_info')
        )

    if 'trace_qry' in cols:
        df = df.with_columns(
            pl.when(pl.col('trace_qry').is_not_null())
            .then(
                pl.col('_vcf_info').list.concat(
                    pl.concat_str(
                        pl.lit('CPX_QRY_TRACE='),
                        pl.col('trace_qry').str.replace_all(',', '')
                    )
                )
            )
            .otherwise('_vcf_info')
            .alias('_vcf_info')
        )

    if len(cols & {'qry_id', 'qry_pos', 'qry_end'}) == 3:
        df = df.with_columns(
            pl.col('_vcf_info').list.concat(
                pl.concat_str(
                    pl.lit('QUERY='),
                    pl.col('qry_id'),
                    pl.lit(':'),
                    pl.col('qry_pos') + 1,
                    pl.lit('-'),
                    pl.col('qry_end')
                )
            )
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
                                pl.element().struct.field('chrom'),
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

    return df


def reformat_vcf_table(
    df: pl.LazyFrame,
    sample_col_prefix: str = '_vcf_sample_'
) -> pl.LazyFrame:
    """Reformat VCF table to VCF format."""
    schema = df.collect_schema()

    strip_len = len(sample_col_prefix)

    sample_col_exprs = [
        (
            pl.col(sample_col)
            .list.eval(  # Replacing semicolons should not be needed, but don't break VCFs if they are present
                pl.element().str.replace_all(';', ':')
            )
        ).list.join(';').alias(sample_col[strip_len:])
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
                pl.col('_vcf_info')
                .list.eval(  # Replacing semicolons should not be needed, but don't break VCFs if they are present
                    pl.element().str.replace_all(';', ':')
                )
            ).list.join(';').alias('INFO'),
            pl.col('_vcf_format').list.join(';').alias('FORMAT'),
            *sample_col_exprs,
        )
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
