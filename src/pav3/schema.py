"""Table schemas."""

from typing import (
    overload,
    Union,
)

import polars as pl

import agglovar

PolarsDataType = Union[pl.DataType, type[pl.DataType]]

__all__ = [
    'ALIGN',
    'VARIANT',
    'cast',
]


ALIGN: dict[str, PolarsDataType] = {
    'chrom': agglovar.schema.VARIANT['chrom'],
    'pos': agglovar.schema.VARIANT['pos'],
    'end': agglovar.schema.VARIANT['end'],
    'align_index': pl.Int32,
    'filter': agglovar.schema.VARIANT['filter'],
    'qry_id': pl.String,
    'qry_pos': pl.Int64,
    'qry_end': pl.Int64,
    'qry_order': pl.Int32,
    'rg': pl.String,
    'mapq': pl.UInt8,
    'is_rev': pl.Boolean,
    'flags': pl.UInt16,
    'align_ops': pl.Struct({
        'op_code': pl.List(pl.UInt8),
        'op_len': pl.List(pl.Int32)
    })
}
"""Schema of alignment tables excluding features (added by align.features)."""

VARIANT: dict[str, PolarsDataType] = {
    'chrom': agglovar.schema.VARIANT['chrom'],
    'pos': agglovar.schema.VARIANT['pos'],
    'end': agglovar.schema.VARIANT['end'],
    'id': agglovar.schema.VARIANT['id'],
    'vartype': agglovar.schema.VARIANT['vartype'],
    'varsubtype': agglovar.schema.VARIANT['vartype'],  # Same as vartype, not defined in agglovar
    'varlen': agglovar.schema.VARIANT['varlen'],
    'ref': agglovar.schema.VARIANT['ref'],
    'alt': agglovar.schema.VARIANT['alt'],
    'filter': ALIGN['filter'],
    'qry_id': ALIGN['qry_id'],
    'qry_pos': ALIGN['qry_pos'],
    'qry_end': ALIGN['qry_end'],
    'qry_rev': pl.Boolean,
    'call_source': pl.String,
    'var_score': pl.Float32,
    'outer_ref': pl.Struct({
        'chrom': pl.String,
        'pos': pl.Int64,
        'end': pl.Int64,
    }),
    'outer_qry': pl.Struct({
        'chrom': pl.String,
        'pos': pl.Int64,
        'end': pl.Int64,
    }),
    'templ_res': pl.Boolean,
    'align_source': pl.List(ALIGN['align_index']),
    'seg': pl.List(
        pl.Struct({
            'chrom': pl.String,
            'pos': pl.Int64,
            'end': pl.Int64,
            'qry_id': pl.String,
            'qry_pos': pl.Int64,
            'qry_end': pl.Int64,
            'is_rev': pl.Boolean,
        })
    ),
    'dup': pl.List(
        pl.Struct({
            'chrom': pl.String,
            'pos': pl.Int64,
            'end': pl.Int64,
            'is_rev': pl.Boolean,
        })
    ),
    'seg_n': pl.Int32,
    'trace_qry': pl.String,
    'trace_ref': pl.String,
    'align_gap': pl.Int32,
    'var_index': pl.UInt32,
    'hom_ref': pl.Struct({'up': pl.Int32, 'dn': pl.Int32}),
    'hom_qry': pl.Struct({'up': pl.Int32, 'dn': pl.Int32}),
    'derived': pl.String,
    'discord': pl.List(pl.String),
    'inner': pl.List(pl.String),
    'seq_ref_pos': pl.Int64,
    'seq_ref_end': pl.Int64,
    'seq': agglovar.schema.VARIANT['seq'],
    'mg_src': pl.List(
        pl.Struct({
            'src_index': pl.Int32,
            'sample': pl.String,
            'hap': pl.String,
            'var_index': pl.UInt32,
            'var_id': pl.String,
        })
    ),
    'mg_stat': pl.List(
        pl.Struct({
            'ro': pl.Float32,
            'offset_dist': pl.Int32,
            'offset_prop': pl.Float32,
            'size_ro': pl.Float32,
            'match_prop': pl.Float32,
            'weight': pl.Float32,
        })
    ),
    'mg_src_lead': pl.UInt32,
}
"""
Schema for variant tables.

Coordinates:
All positions (pos, end, qry_pos, and qry_end) are 0-based, half-open (BED-like coordinates).
Query coordinates are in the original query sequence coordinates even if it was reverse-complemented during alignment
(differs from SAM/BAM). By preserving the original query coordinates, it is possible to extract a variant locus
directly from the assembly sequence without first translating the coordinates or extracting and transforming a whole
assembly sequence (e.g. use pysam or samtools faidx to extract the variant locus, then reverse-complement the result if
necessary). This also allows variant calls to be directly used as annotations on the assembly sequence without
transforming the coordinates.

Outer inversions:
Inversions have inner- and outer- coordinates. Inner coordinates identify the uniquely-inverted region of an inversion,
and these coordinates are found in pos, end, qry_pos, and qry_end. Outer coordinates identify the inversion including
inverted repeats flanking the uniquely inverted region, and these coordinates are found in outer_ref and outer_qry.
The outer coordinates, if defined, should encapsulated the inner coordinates. Some variant callers use the outer
coordinates to define the inversion locus. The true breakpoint is typically somewhere between the inner- and
outer-coordinates.

Fields:

    * chrom: Reference chromosome.
    * pos: Reference position (see notes above - Coordinates).
    * end: End position (see notes above - Coordinates).
    * id: Variant ID. Must be unique in the final variant table, intermediate tables may have duplicate IDs before
        variant integration resolves them.
    * vartype: Variant type, uppercase. Should be one of "SNV", "INS", "DEL", "INV", "CPX".
    * varsubtype: Variant subtype, uppercase. "TD" indicates a tandem duplication.
    * varlen: Length of the variant. Must be defined for all but SNV variants, which are 1 bp by definition.
    * ref: Reference base for SNVs.
    * alt: Alternate base for SNVs.
    * filter: List of filters applied to this variant. An empty list indicates that the variant passed all filters (i.e.
        a "PASS" variant). Filters should be uppercase. Null values should not be generated by PAV, but can be assumed
        to mean that filter status is unknown and should not be interpreted as a "PASS" variant.
    * qry_id: Query ID. Name of the assembly sequence variant was derived from.
    * qry_pos: Query position (see notes above - Coordinates).
    * qry_end: Query end position (see notes above - Coordinates).
    * qry_rev: If the sequence was reverse-complemented during alignment (see notes above - Coordinates).
    * call_source: A short string in all caps indicating how the variant was identified. "INTRA" is intra-alignment
        variant calling, "INTER" is inter-alignment variant calling from split alignment records.
    * var_score: Variant score based on a scoring model.
    * outer_ref: Outer reference region for inversions (see notes above - Outer inversions).
    * outer_qry: Outer query region for inversions (see notes above - Outer inversions).
    * templ_res: Template resolution flag. True if the variant is aligned across its segments, False if variant
        segments are unaligned. Common in complex variants where smaller template switches are left unaligned.
    * align_source: List of "align_index" values from the alignment table this variant was derived from.
        Intra-alignment variants will have a single value. Inter-alignment variants will have 0 or more values and will
        exlude anchors. Indexes will be in order of the alignment chain from left to right in reference orientation.
    * seg: A list of segments belonging to this variant. Each element links a reference region and query region. Field
        "is_rev" is true when the segment is inverted relative to the reference (assembly sequence orientation does not
        affect this value).
    * dup: List of duplicated reference regions where a precise query coordinate is unknown.
    * seg_n: Number of distinct segments in a complex variant.
    * trace_qry: Formatted string tracing the query structure of a complex variant.
    * trace_ref: Formatted string tracing the reference structure of a complex variant.
    * align_gap: Difference in the number of bases between the query and reference sequences for inter-alignment
        inversions. Large numbers suggest a more complex arrangement than a balanced inversion.
    * var_index: A unique index for variants within a callset.
    * hom_ref: Perfect homology flanking the variant breakpoints in reference coordinates. Field "up" is the number of
        homologous bases upstream of the variant and "dn" is the number downstream.
    * hom_qry: Perfect homology flanking the variant breakpoints in query coordinates. Fields are the same as
        "hom_ref".
    * derived: ID of a variant this variant was derived from (where applicable). If non-null, then "DERIVED" must also
        be present in "filter".
    * discord: IDs of variants this variant was discordant with. For example, a small variant inside a larger deletion
        will have the deletion variant ID in this list). If not empty, then "DISCORD" must also be present in "filter".
    * inner: IDs of a variant this variant was inside of. These are typically small variants inside aligned segments
        that are part of a larger variant (i.e. duplicated region of an insertion or complex event). If not empty, then
        "INNER" must also be present in "filter".
    * seq_ref_pos: Complex variants only. Start of the reference region "seq" replaces (see seq_ref_end).
    * seq_ref_end: Complex variants only. End of the reference region "seq" replaces. Replacing reference bases
        [seq_ref_pos, seq_ref_end) with "seq" reconstructs the haplotype across the variant. The region is empty
        (seq_ref_pos == seq_ref_end) when nothing is deleted. It differs from [pos, end) when the anchoring alignments
        overlap on the reference (duplication at the outer breakpoint), where [pos, end) is the overlap.
    * seq: Variant sequence in reference orientation. Null if not defined (never "*" or other placeholders).

Merge fields:

The last three fields ("mg_src", "mg_stat", and "mg_src_lead") are only present in merged callsets (i.e. haplotypes
merged into a sample callset by `pav3.workflow.call.merge_haplotypes`). They are absent from per-haplotype callsets.

    * mg_src: List of source entries for a merged variant, one per haplotype that contributed to the merged call.
        This is a PAV-specific reshaping of the merge source struct emitted by agglovar: the agglovar source-name and
        metadata fields are resolved into explicit "sample" and "hap" fields (the metadata is dropped once consumed).
        Each struct has fields:

            * src_index: Source index (0-based input order of the merged callset).
            * sample: Sample (assembly) name the source variant came from.
            * hap: Haplotype the source variant came from.
            * var_index: Index of the source variant within its callset ("var_index" of the source table).
            * var_id: ID of the source variant within its callset.

    * mg_stat: List of merge statistics emitted by agglovar, one per source variant merged into the lead variant
        (i.e. one fewer entry than "mg_src"). The struct fields are the join statistics computed by the pairwise
        overlap engine. The struct does not vary by variant type: agglovar always emits the full set of default join
        columns, and a statistic that a merge does not compute is emitted as a null rather than being dropped from
        the struct. The fields are:

            * ro: Reciprocal overlap with the lead variant.
            * offset_dist: Breakpoint offset distance (max of the start and end offsets) in bp.
            * offset_prop: Offset distance over the smaller of the two variant lengths.
            * size_ro: Reciprocal size overlap (size of the smaller over the size of the larger).
            * match_prop: Proportion of matched bases in an alignment of the two variant sequences. Only computed
                when the merge parameters request sequence matching ("match_prop_min"), which for the defaults in
                `pav3.const.DEFAULT_MERGE_PARAMS` is INS/DEL and CPX. Null for SNV, INV, and DUP.
            * weight: Weight assigned to the match by the merge parameters, which chooses among candidate matches.

        The struct is fixed by the merge parameters, not by PAV, so non-default merge parameters can change it. The
        one field the defaults use but never report is "seg_ro" (segment reciprocal overlap, used by the CPX merge):
        agglovar computes it to filter candidate matches but does not carry it into the output statistics.

    * mg_src_lead: Index into "mg_src" of the lead (representative) source entry. The lead source is the variant whose
        columns were carried into the merged record, so
        `pl.col('mg_src').list.get(pl.col('mg_src_lead'))` identifies the callset and variant the merged record's
        coordinates, sequence, and annotations came from.
"""


@overload
def cast(
    df: pl.DataFrame,
    schema: dict[str, PolarsDataType],
    do_sort: bool = True,
) -> pl.DataFrame:
    ...


@overload
def cast(
    df: pl.LazyFrame,
    schema: dict[str, PolarsDataType],
    do_sort: bool = True,
) -> pl.LazyFrame:
    ...


def cast(
    df: pl.DataFrame | pl.LazyFrame,
    schema: dict[str, PolarsDataType],
    do_sort: bool = True,
) -> pl.DataFrame | pl.LazyFrame:
    """
    Cast a DataFrame or LazyFrame to a schema.

    :param df: DataFrame or LazyFrame to cast.
    :param schema: Schema to cast to.
    :param do_sort: Whether to sort the DataFrame or LazyFrame by the order defined in the schema.
        If the table contains columns not defined by the schema, they are retained in order after
        columns defined by the schema.
    """

    if isinstance(df, pl.DataFrame):
        existing_schema = df.schema
    else:
        existing_schema = df.collect_schema()

    df = df.cast(
        {
            col: schema.get(col, type_)
            for col, type_ in existing_schema.items()
        }
    )

    if do_sort:
        df = (
            df
            .select(
                *[  # Known columns first
                    col for col in schema.keys() if col in existing_schema.keys()
                ],
                *[  # Unknown columns last
                    col for col in existing_schema.keys() if col not in schema.keys()
                ]
            )
        )

    return df
