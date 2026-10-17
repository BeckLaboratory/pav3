"""Variant calling workflow tasks"""

from collections import defaultdict
from collections.abc import Iterable
import itertools
import logging
import os
from pathlib import Path
import polars as pl
from typing import (
    Any,
    Optional,
)

import pav3
import agglovar

from .. import const

from ..io import TempDirContainer

from ..call.expr import sort_expr

logger = logging.getLogger(__name__)


def _materialize(
        df: pl.LazyFrame,
        temp_file_container: TempDirContainer,
        label: str,
) -> pl.LazyFrame:
    """Write a lazy frame to a temporary Parquet file and return a streaming scan of it.

    Used as a materialization barrier between lazy operations. It serves two purposes:

    1. Memory: Variant tables may be very large. Writing intermediate results to a temporary file and re-reading
       them as streaming scans keeps them off the Python heap instead of accumulating collected frames in memory
       across the integration loop.

    2. Correctness: `apply_trim_filter` and `apply_discord_and_inner_filter` each add a row-index column
       (`with_row_index`). When both are chained into a single lazy plan and that plan is sunk, the Polars
       new-streaming engine mis-optimizes the chained row-index operations and raises a spurious
       "duplicate name when building schema" / "column not found" error. The same plan collects correctly with
       the in-memory engine, confirming the logical plan is valid and the fault is in the streaming optimizer.
       Breaking the plan with a Parquet round-trip isolates each helper's row index and avoids the bug.

    :param df: Lazy frame to materialize.
    :param temp_file_container: Active temporary file container (must be inside its context).
    :param label: Label used in the temporary file name.

    :return: Streaming scan over the materialized Parquet file.
    """
    temp_file = temp_file_container.next(prefix=f'integrate_{label}_', suffix='.parquet')
    df.sink_parquet(temp_file)
    return pl.scan_parquet(temp_file)


def integrate_sources(
        pav_params,
        input: dict[str, str],
        output: dict[str, str],
        temp_dir_name: Optional[str] = None,
) -> None:
    """Integrate variant calls from all sources into the final per-haplotype call tables.

    :param pav_params: PAV parameters.
    :param input: Input file dictionary (Snakemake rule input).
    :param output: Output file dictionary (Snakemake rule output).
    :param temp_dir_name: Parent directory for temporary files. Defaults to the system temporary directory.
    """
    # Intermediate variant tables are materialized to a temporary directory and re-read as streaming scans (see
    # `_materialize`). This bounds memory for large callsets and works around a Polars new-streaming engine bug.
    with TempDirContainer(temp_dir_name, prefix='pav_integrate_') as temp_file_container:
        _integrate_sources(pav_params, input, output, temp_file_container)


def _integrate_sources(
        pav_params,
        input: dict[str, str],
        output: dict[str, str],
        temp_file_container: TempDirContainer,
) -> None:
    inv_min = pav_params.inv_min
    inv_max = pav_params.inv_max if pav_params.inv_max > 0 else float('inf')

    # Read alignments
    df_align_none = pl.scan_parquet(input['align_none'])
    df_align_qry = pl.scan_parquet(input['align_qry'])
    df_align_qryref = pl.scan_parquet(input['align_qryref'])

    # Read trimmed regions (regions are tuples of alignment records and coordinates within the record).
    # IntervalTree where coordinates are tuples - (index, pos):(index, end)
    df_trim = pav3.call.integrate.read_trim_table(
        df_align_none,
        df_align_qry,
        df_align_qryref,
    ).lazy()

    del df_align_none

    # Create a table of regions added to the DISCORD filter
    df_discord_schema = {
        key: val for key, val in pav3.schema.VARIANT.items() if key in {'chrom', 'pos', 'end', 'id'}
    }

    df_discord = pl.LazyFrame([], schema=df_discord_schema)

    # Save alignment records for INNER variants
    df_inner_schema = {
        'align_index': pav3.schema.ALIGN['align_index'],
        'id': pav3.schema.VARIANT['id'],
    }

    df_inner = pl.LazyFrame([], schema=df_inner_schema)

    # Table connecting var_index to ids
    df_var_index_id_schema = {
        'var_index': pav3.schema.VARIANT['var_index'],
        'id': pav3.schema.VARIANT['id'],
    }

    df_var_index_id = pl.LazyFrame([], schema=df_var_index_id_schema)

    # Read segment table
    df_segment = pl.scan_parquet(input['inter_segment'])

    # List of variant tables to collect from multiple sources
    collect_list = defaultdict(list)

    # Parameters controlling how variants are integrated and in what order
    param_dict = {  # do_write, add_discord, filter_discord, filter_inner
        'inter_cpx':    (True,     True,        False,          False),
        'inter_insdel': (False,    True,        False,          False),
        'inter_inv':    (False,    True,        False,          False),
        'intra_inv':    (True,     True,        True,           True),
        'intra_insdel': (True,     False,       True,           True),
        'intra_snv':    (True,     False,       True,           True),
    }

    # do_write: Write variant call table. Set to False for INS/DEL or INV variants until they are all collected
    # add_discord: Add variant regions to the DISCORD regions.
    # filter_discord: Apply DISCORD filter.
    # filter_inner: Apply INNER filter.
    #
    # Note: add_inner is implied by vartype == 'cpx'

    for sourcetype_vartype in param_dict.keys():
        # if sourcetype_vartype == 'intra_snv':
        #     raise RuntimeError(f'Stopping at {sourcetype_vartype}')

        if pav_params.debug:
            print(f'Processing {sourcetype_vartype}')

        do_write, add_discord, filter_discord, filter_inner = param_dict[sourcetype_vartype]

        sourcetype, vartype = sourcetype_vartype.rsplit('_', 1)

        is_lg = sourcetype == 'inter'
        is_insdel, is_inv, is_snv, is_cpx = (vartype == val for val in ('insdel', 'inv', 'snv', 'cpx'))

        # Read variant table
        df = pav3.call.integrate.set_base_cols(
            pl.scan_parquet(input[sourcetype_vartype]),
            sourcetype,
            vartype,
        )

        # Apply variant length filters
        if is_inv:
            df = (
                df
                .with_columns(
                    pl.when((pl.col('varlen') < inv_min) | (pl.col('varlen') > inv_max))
                    .then(pl.col('filter').list.concat([pl.lit('VARLEN')]))
                    .otherwise(pl.col('filter'))
                    .alias('filter')
                )
            )

        # Filter TRIMREF & TRIMQRY
        if not (is_lg or is_inv):
            df = pav3.call.integrate.apply_trim_filter(df, df_trim)

            # Materialization barrier between the two row-index helpers. `apply_trim_filter` and
            # `apply_discord_and_inner_filter` both call `with_row_index`; chaining them into one sunk plan trips a
            # Polars streaming-engine bug (see `_materialize`). The round-trip isolates the two row indices.
            df = _materialize(df, temp_file_container, f'{sourcetype_vartype}_trim')

        # Filter DISCORD
        if filter_discord or filter_inner:
            df = pav3.call.integrate.apply_discord_and_inner_filter(
                df,
                df_discord if filter_discord else None,
                df_inner if filter_inner else None,
            )

        # Version variant IDs prioritizing PASS over non-PASS.
        df = pav3.call.integrate.id_and_version(
            df=df,
            is_snv=is_snv,
            existing_ids=collect_list[vartype]
        )

        # Materialize the fully-processed variant table and continue from a streaming scan. Every downstream use
        # (DISCORD regions, CPX-derived variants, writes, and cross-source aggregation in `collect_list`) reads from
        # this scan. This keeps large tables off the heap and replaces the previous in-memory `df.collect().lazy()`.
        df = _materialize(df, temp_file_container, sourcetype_vartype)

        # Read CPX segment table
        if sourcetype == 'inter':
            df_segment_var = df_segment.join(df.select(['var_index', 'id']), on='var_index', how='left')
        else:
            df_segment_var = None

        # Get discord expr
        if add_discord and not pav_params.redundant_callset:
            update_discord_frame = (
                pl.concat(
                    [
                        df_discord,
                        (
                            df
                            .filter(
                                pl.col('filter').list.len() == 0,
                                pl.col('end') > (pl.col('pos') + 1)
                            )
                            .select(['chrom', 'pos', 'end', 'id'])
                        )
                    ]
                )
                .sort(['chrom', 'pos', 'end'])
            )
        else:
            update_discord_frame = df_discord

        # Add variants derived from partial complex events
        if is_cpx:
            pav3.call.integrate.add_cpx_derived(
                df=df,
                df_segment=df_segment_var,
                collect_list=collect_list,
                ref_fa=pav3.pipeline.expand_fofn(input['ref_fofn'])[0],
            )

        # Add duplications derived from insertions placing a copy of a reference region
        pav3.call.integrate.add_ins_derived(
            df=df,
            collect_list=collect_list,
            ref_fa=pav3.pipeline.expand_fofn(input['ref_fofn'])[0],
        )

        # Aggregate if type is split over multiple inputs. `df` is already a disk-backed scan (materialized above),
        # so this keeps memory bounded even when several sources accumulate in `collect_list`.
        if not do_write or len(collect_list[vartype]) > 0:
            collect_list[vartype].append(df)

        # Write
        write_list = [update_discord_frame]

        if sourcetype == 'inter':
            write_list.append(df.select(['var_index', 'id']))
            collect_index_var_index = len(write_list) - 1
        else:
            collect_index_var_index = None

        if df_segment_var is not None:
            write_list.append(
                df_segment_var
                .filter(
                    ~ pl.col('is_anchor')
                    & pl.col('is_aligned')
                    & pl.col('align_index').is_not_null()
                    & pl.col('id').is_not_null()
                )
                .join(  # Passing variants only
                    (
                        df
                        .filter(pl.col('filter').list.len() == 0)
                        .select('var_index')
                    ),
                    on='var_index',
                    how='inner'
                )
                .select(
                    'align_index', 'id'
                )
            )
            collect_index_inner = len(write_list) - 1
        else:
            collect_index_inner = None

        if do_write:
            if collect_list[vartype]:
                df = pl.concat(collect_list[vartype], how='diagonal')

            # Sort and order columns
            col_names = df.collect_schema().names()

            df = (
                df
                .with_row_index('_concat_index')
                .with_columns(
                    pl.col('filter').list.unique().list.sort()
                )
                .sort(
                    pl.col('filter').list.len(),
                    '_concat_index',
                    descending=(True, False),
                )
            )

            df = pav3.call.integrate.id_and_version(
                df=df, is_snv=is_snv,
            )

            df = (
                df
                .sort(pav3.call.expr.sort_expr())
                .select([col for col in pav3.schema.VARIANT.keys() if col in col_names])
            )

            if 'inner' in col_names:
                df = df.with_columns(pl.col('inner').fill_null([]))

            if 'discord' in col_names:
                df = df.with_columns(pl.col('discord').fill_null([]))

            # Create output objects
            write_list.append(df.sink_parquet(output[vartype], lazy=True))

        # Collect and write
        # Always run even if not do_write to update discord regions
        collect_all_list = pl.collect_all(write_list)

        df_discord = collect_all_list[0].lazy()

        if collect_index_var_index is not None:
            df_var_index_id = pl.concat(
                [
                    df_var_index_id.collect(),
                    collect_all_list[collect_index_var_index],
                ]
            ).lazy()

        if collect_index_inner is not None:
            df_inner = pl.concat(
                [
                    df_inner.collect(),
                    collect_all_list[collect_index_inner],
                ]
            ).lazy()

    # Write duplications
    if pav_params.debug:
        print('Writing dup')

    df = pl.concat(collect_list['dup'], how='diagonal')

    (
        df
        .with_columns(pl.col('filter').list.unique().sort())
        .sort(pav3.call.expr.sort_expr())
        .select([col for col in pav3.schema.VARIANT.keys() if col in df.collect_schema().names()])
        .sink_parquet(output['dup'])
    )

    # Write segment and ref_trace tables
    if pav_params.debug:
        print('Writing segment & trace tables')

    (
        df_segment
        .join(df_var_index_id, on='var_index', how='left')
        .sink_parquet(output['inter_segment'])
    )

    (
        pl.scan_parquet(input['inter_ref_trace'])
        .join(df_var_index_id, on='var_index', how='left')
        .sink_parquet(output['inter_ref_trace'])
    )


def _mg_src_sample_hap_expr() -> pl.Expr:
    """Reshape the agglovar merge source list into PAV's "mg_src" schema.

    agglovar emits each "mg_src" entry as a struct with "src_name" and "src_meta" (a JSON string
    holding the callset metadata PAV attaches per source, i.e. "sample" and "hap"). This resolves
    "sample" and "hap" into explicit struct fields from the metadata and drops both "src_name" and
    "src_meta", producing the "mg_src" struct defined by `pav3.schema.VARIANT`. Reading the
    haplotype from metadata avoids splitting the source name on "-", which breaks when a sample or
    haplotype name itself contains "-".

    :return: An expression replacing the "mg_src" column with the PAV-specific struct schema.
    """
    meta = pl.element().struct.field('src_meta')

    return (
        pl.col('mg_src').list.eval(
            pl.struct(
                pl.element().struct.field('src_index').cast(pl.Int32).alias('src_index'),
                meta.str.json_path_match('$.sample').cast(pl.String).alias('sample'),
                meta.str.json_path_match('$.hap').cast(pl.String).alias('hap'),
                pl.element().struct.field('var_index').cast(pl.UInt32).alias('var_index'),
                pl.element().struct.field('var_id').cast(pl.String).alias('var_id'),
            )
        ).alias('mg_src')
    )


def get_merge_params(
        vartype: str,
        merge_params_override: Optional[dict[str, Any]] = None,
) -> list[dict[str, Any]]:
    """Get haplotype merge parameters for a variant type.

    :param vartype: Variant type ("insdel", "snv", "inv", "cpx", "dup"; "ins" and "del" are aliases for "insdel").
    :param merge_params_override: Parameters overriding `pav3.const.DEFAULT_MERGE_PARAMS` keyed by variant type
        (config parameter "merge_params"). Each value is a list of pairwise overlap stages or a single stage. Variant
        types missing from this dictionary keep the defaults.

    :returns: A list of pairwise overlap stage definitions.

    :raises ValueError: If `vartype` or a key in `merge_params_override` is not a known variant type, or if a value is
        not a stage or a list of stages.
    """
    if vartype is not None and (vartype := vartype.lower().strip()) in {'ins', 'del'}:
        vartype = 'insdel'

    if vartype not in const.DEFAULT_MERGE_PARAMS:
        raise ValueError(f'Unknown variant type for merge parameters: "{vartype}"')

    if merge_params_override is not None:
        unknown_keys = set(merge_params_override.keys()) - set(const.DEFAULT_MERGE_PARAMS.keys())

        if unknown_keys:
            raise ValueError(
                f'Unknown variant type(s) in merge parameters: {", ".join(sorted(unknown_keys))} '
                f'(expected {", ".join(const.DEFAULT_MERGE_PARAMS.keys())})'
            )

        if vartype in merge_params_override:
            stage_list = merge_params_override[vartype]

            if isinstance(stage_list, dict):
                stage_list = [stage_list]

            if not isinstance(stage_list, list) or not stage_list or not all(
                isinstance(stage, dict) for stage in stage_list
            ):
                raise ValueError(
                    f'Merge parameters for "{vartype}" must be a stage (object) or a non-empty list of stages'
                )

            return [dict(stage) for stage in stage_list]

    return [dict(stage) for stage in const.DEFAULT_MERGE_PARAMS[vartype]]


def merge_haplotypes(
        vartype: str,
        callsets: Iterable[agglovar.merge.base.CallsetDefInputType],
        ref_path: str | Path,
        out_path: str | Path,
        merge_params: Optional[Iterable[dict[str, Any]]] = None,
        merge_name: Optional[str] = None,
        temp_file_container: Optional[TempDirContainer] = None,
) -> None:
    """Create a single callset for a phased assembly from all haplotypes.

    :param vartype: Variant type. Used to make some decisions about merging beyond the merge parameters.
    :param callsets: An iterable of callsets to merge. Each element is a two-element tuple with a
        LazyFrame (first element) and the name of the source  callset (second element). The source
        callset is typically "NAME_HAP" (assembly name and haplotype), but can be anything. A default
        source name is set if missing.
    :param ref_path: Path to reference file. Must have a ".fai" index.
    :param out_path: Path to output file.
    :param merge_params: Parameters controlling this merge. If not set, parameters are retrieved from
        `pav3.const.DEFAULT_MERGE_PARAMS` with key `vartype` ("ins" and "del" are aliases for "insdel").
    :param merge_name: Name of this merge for tracking in log messages. If not set, omitted from log messages.
    :param temp_file_container: Manages temporary files written by the merge. If not set, a default container
        is used.
    """
    known_vartype = ['ins', 'del', 'insdel', 'snv', 'inv', 'cpx', 'dup']

    merge_sort_expr = sort_expr(has_id=False)

    if vartype is None or (vartype := vartype.lower().strip()) not in known_vartype:
        raise ValueError(f'Unknown variant type "{vartype}" (expected one of {", ".join(known_vartype)})')

    if merge_params is None:
        merge_params = get_merge_params(vartype)
    else:
        merge_params = list(merge_params)

    ref_path = Path(ref_path)

    if not ref_path.is_file():
        raise FileNotFoundError(f'Reference Path does not exist or is not a regular file: {ref_path}')

    ref_fai_path = Path(str(ref_path) + '.fai')

    if not ref_fai_path.is_file():
        raise FileNotFoundError(f'Reference Path to FAI index does not exist or is not a regular file: {ref_fai_path}')

    callsets = tuple(callsets)

    if merge_name is not None and (merge_name := str(merge_name).strip()):
        log_prefix = f'{merge_name}: '
    else:
        merge_name = None
        log_prefix = ''

    if out_path is None:
        raise ValueError('Output file is missing')

    out_path = Path(out_path)

    os.makedirs(out_path.parent, exist_ok=True)

    # Setup join object
    pairwise_join = agglovar.pairwise.overlap.PairwiseOverlap.from_definition(
        merge_params
    )

    merge_runner = agglovar.merge.cumulative.MergeCumulative(
        pairwise_join,
        lead_strategy=agglovar.merge.cumulative.LeadStrategy.LEFT,
    )

    # Get chromosome list
    with open(ref_fai_path, 'rt') as in_file:
        chrom_list = sorted([
            chrom for chrom in (
                line.split('\t')[0].strip() for line in in_file
            ) if chrom
        ])

    # Get filters for variant type
    if vartype == 'insdel':
        vartype_list = ['INS', 'DEL']
    else:
        vartype_list = [vartype.upper(),]

    if temp_file_container is None:
        prefix_name = 'pav3_merge_hap'
        if merge_name is not None:
            prefix_name += f'_{merge_name}'

        temp_file_container = TempDirContainer(prefix=prefix_name)

    with temp_file_container:
        concat_chrom_list = []

        for chrom in chrom_list:
            df_chrom_list = []

            for vartype, filter_pass in itertools.product(vartype_list, (True, False)):
                logger.debug(
                    '%sMerging chromosome %s (vartype=%s, pass=%s): ' % (log_prefix, chrom, vartype, filter_pass)
                )

                # Remove duplicate calls within a haplotype. The key includes the variant type and sequence so
                # that distinct alleles at the same site (e.g. a 1 bp INS and DEL at one position, or two INS of the
                # same length with different sequences) are not dropped.
                if vartype == 'SNV':
                    dedup_cols = ['chrom', 'pos', 'alt']
                elif vartype in {'INS', 'DEL', 'CPX'}:
                    dedup_cols = ['chrom', 'pos', 'end', 'varlen', 'vartype', 'seq']
                else:
                    dedup_cols = ['chrom', 'pos', 'end', 'varlen', 'vartype']

                pre_filter = [
                    pl.col('chrom') == chrom,
                    pl.col('vartype').str.to_uppercase() == vartype,
                    pl.struct(dedup_cols).is_first_distinct(),
                    (
                        pl.col('filter').list.len() == 0
                        if filter_pass else
                        pl.col('filter').list.len() > 0
                    )
                ]

                next_filename = temp_file_container.next(
                    prefix=f'split_{chrom}_{vartype if vartype else vartype}_{filter_pass}_'
                )

                df_chrom_list.append(next_filename)

                (
                    merge_runner(
                        callsets,
                        retain_index=True,
                        pre_filter=pre_filter,
                    )
                    .with_columns(_mg_src_sample_hap_expr())
                    .sort(merge_sort_expr)
                    .sink_parquet(next_filename)
                )

            # Merge and sort this chromosome
            chrom_next_filename = temp_file_container.next(
                prefix=f'concat_{chrom}_'
            )

            concat_chrom_list.append(chrom_next_filename)

            # Version IDs across all passes for this chromosome, not within each one.
            #
            # Passing and non-passing variants are merged separately, so a base ID present in both would take the
            # unversioned ID in each pass and collide. Versioning after the concatenation ranks them together, and
            # "id_version_expr" ranks on the filter count, so passing variants keep the unversioned ID. IDs carry
            # the chromosome and variant type, so they cannot collide across chromosomes or merge passes.
            (
                pl.concat([
                    pl.scan_parquet(filename)
                    for filename in df_chrom_list
                ])
                .sort('chrom', 'pos', 'end', 'id')
                .with_columns(agglovar.util.var.id_version_expr())
                .sink_parquet(chrom_next_filename)
            )

        # Merge all
        logger.debug('%sConcatenating chromosomes' % log_prefix)

        (
            pl.concat([
                pl.scan_parquet(filename)
                for filename in concat_chrom_list
            ])
            .sink_parquet(out_path)
        )

        logger.debug('%sDone' % log_prefix)
