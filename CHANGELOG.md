# Changelog

<!--
Style guide:
https://common-changelog.org/
-->

<!-- TEMPLATE
## [x.y.z] - YYYY-MM-DD

### Changed
- 

### Added
- 

### Removed
-

### Fixed
-

-->

## 3.0.0.dev27

### Changed
- Low-confidence alignment filtering (`LCALIGN`) is off by default (config parameter `lc_model` defaults to
  `none`). It removed correct calls in highly diverged regions such as HLA class II. Set `lc_model` to `default` for
  the previous behavior.
- Alignment records with MAPQ below 20 are filtered by default (new config parameter `min_mapq`). Their variants
  are filtered with `MAPQ`, and their regions are not callable. Low-MAPQ fragments in satellites and repeats
  produced false calls that broke haplotype reconstruction from the VCF.
- Haplotype merge defaults (`pav3.const.DEFAULT_MERGE_PARAMS`) for `insdel` and `cpx` merge only matching alleles:
  same position and end, size and sequence identity of at least 0.998 for insertions and deletions, and 0.98 for
  complex variants (with segment overlap 0.8). The previous population-style tolerances merged different alleles
  from the two haplotypes into one homozygous record and moved variants up to 200 bp to the other haplotype's
  position. Override with the new `merge_params` config parameter.
- VCF records at the same position are combined into multi-allelic records when they have the same filter and no
  haplotype carries more than one of them, so different alleles on two haplotypes are written as one record (e.g.
  `1|2`) instead of two overlapping records (`1|0` and `0|1`).
- Complex variants are written as sequence-resolved REF and ALT alleles by default. The record replaces the
  reference bases skipped by the variant with the variant sequence, so applying it to the reference (e.g.
  `bcftools consensus`) reconstructs the haplotype. Set `vcf_cpx_symbolic` to write `<CPX>` alleles.
- Complex variant tables have new columns `seq_ref_pos` and `seq_ref_end`, the reference region the variant
  sequence replaces.
- Tighter compliance with VCF 4.5.
- VCF INFO headers declare an accurate `Number`. `END` is `Number=1` and `SVLEN` is `Number=A`, matching the
  reserved keys; all others are `Number=A`. All were previously `Number=.`.
- VCF `INFO/SVLEN` is always positive. Deletions are no longer negated.
- VCF INFO values no longer contain literal commas, which VCF reserves as the value delimiter. `DUP` separates
  regions with `/`, and `CPX_REF_TRACE` and `CPX_QRY_TRACE` drop thousands separators. The `trace_ref` and
  `trace_qry` variant table columns are unchanged.
- VCF `INFO/END` and `INFO/SEQ` are written only for symbolic ALT alleles, where REF and ALT do not carry the
  span or the sequence. `SEQ` is upper-cased to match REF and ALT, and is not written for derived variants, where
  it does not describe the alternate allele.
- Tandem and inverted duplications are written as `<DUP:TANDEM>` and `<DUP:INVDUP>` rather than `INFO/SUBTYPE`.
- Replaced VCF `INFO/HOMREF` and `INFO/HOMQRY` with the reserved `INFO/HOMLEN`, summing upstream and downstream
  homology from reference breakpoints for `DEL` and `DUP` and query breakpoints for `INS` and `CPX`.
- Variant tracks report `.` for strand. The strand column was never set, and a variant carries no orientation of
  the reference feature it is drawn on.

### Added
- YAML configuration file `pav.yaml` (or `pav.yml`), accepting the same parameters as `pav.json`. PAV reads
  `pav.json` only if no YAML configuration is found, and warns and uses the YAML configuration if both are present.
- Config parameters `min_mapq`, `merge_params`, and `vcf_cpx_symbolic`.
- VCF `INFO/SVCLAIM` is written as `J` for symbolic `DEL` and `DUP` alleles, which VCF requires. PAV calls these
  from non-reference adjacencies rather than from a measure of DNA abundance.
- Duplications derived from insertions placing a copy of a reference region. A tandem-duplication insertion
  derives a tandem duplication.
- Insertions over deleted reference bases and deletions with inserted bases now emit both variants instead of the
  larger one alone. The smaller component was previously discarded, leaving an allele that reported no inserted
  sequence.
- Documentation page describing the VCF output.
- Config parameter `perf_log` logs the CPU and memory usage of each rule from input assemblies to the VCF to
  `log/perf/{asm_name}/{rule}_{wildcards}.parquet` (e.g. `log/perf/HG00171/call_intra_snv_insdel_h1.parquet`).
  `true` takes a snapshot every 30 seconds, a positive integer sets the seconds between snapshots, and `false` or `0`
  (default) disables it. The table is rewritten after every snapshot, so a job killed for exceeding its memory keeps
  its samples. Adjunct rules (tracks, alignment depth) are not logged.

### Removed
- Snakemake `benchmark` directives (`log/benchmark/`), replaced by `perf_log`.
- VCF `INFO/SVTYPE`. VCF 4.5 deprecates it as redundant with ALT, which carries the type in the symbolic allele
  or in the REF and ALT lengths.
- The `##FILTER=<ID=PASS>` header. VCF defines PASS implicitly.

### Fixed
- Alignment trimming truncated the passing alignment instead of the filtered one where a passing and a filtered
  alignment overlapped (the filtered flag was inverted).
- Alignment trimming applied the left-aligning tie-break to alignments in opposite orientations instead of
  alignments in the same orientation.
- Alignment trimming removed a passing alignment contained in a higher-scoring filtered alignment. The filtered
  alignment is now removed.
- The haplotype merge dropped one of two variants with the same position, end, and length on one haplotype (e.g. a
  1 bp insertion and deletion at one site). Duplicates are now identified with the variant type and sequence.
- Symbolic VCF records (INV, DUP) had the first base of the variant as `REF` instead of the base at `POS`.
- VCF FORMAT keys and sample values were separated by `;` instead of `:`.
- An unknown `lc_model` name is an error instead of silently disabling LC filtering (path-like names) or raising an
  unrelated import error (plain names).
- Variant IDs were duplicated when a variant merged into both the passing and non-passing haplotype merges. IDs
  are now versioned across both merges rather than within each.
- Deletions derived from complex variants were written with an empty VCF `REF`.
- Insertions at a tandem duplication were written with `REF` matching `ALT`, encoding no variant.
- Writing a VCF now fails, naming the variants, if a record has an empty `REF`, an empty `ALT`, or a `REF`
  matching its `ALT`.
- Query homology was read from the reference homology field.
- VCF `INFO/DUP` was dropped for variants duplicating an even number of regions.
- `DUP` and `CPX` records carried no `INFO/END`, so they indexed as 1 bp and could be missed by region queries.
- Large deletions reported an empty query region, which wrote `INFO/QUERY` as an inverted range. Deletions now
  report a 1 bp query region at the breakpoint, matching intra-alignment deletions. Query homology is still
  measured over the query region between the anchors.
- Complex variants reported the query region between the untrimmed anchors, which excludes bases duplicated at
  the variant site, leaving `INFO/QUERY` shorter than `INFO/SVLEN`. They now report the region from the
  query-reference-trimmed alignments, as insertions do, when that region spans the duplicated bases.
- Insertions and complex variants took the duplicated bases at the variant site from the query region of the
  query-reference-trimmed alignments without checking that reference trimming resolved the overlap between the
  anchors. Where it did not, a tandem duplication was written with no sequence at all. The duplicated bases are
  again read from the reference in that case, and the variant reports a 1 bp query region rather than an empty
  one.
- Perfect homology was reported as zero when it ran to the end of the compared flank, understating `INFO/HOMLEN`
  for short variants in low-complexity sequence.
- The assertion comparing the query-trimmed and query-reference-trimmed reference regions of an alignment
  interval joined its conditions with `or`, so it passed whenever any one of them held.
- INV/DUP/CPX tracks failed to build with `Trailing characters parsing integer in field 7` once the duplication
  callset carried duplications derived from insertions. Duplications are the only records without their own
  thick-line coordinates, and the missing values were filled with the literal text `pos` and `end` instead of the
  variant position, which `bedToBigBed` rejected.
- Track columns typed as numbers in the AutoSQL table are now checked before the BED is written, naming the
  column rather than failing inside `bedToBigBed`. UCSC itemRgb columns, which are declared `uint` but hold
  `R,G,B` strings, are exempt.

## 3.0.0.dev26 - 2026-08-28

### Changed
- Version bump to mitigate version conflicts with built containers with a stale version number.

## 3.0.0.dev25 - 2026-08-28

### Added
- `lift`, `subseq`, and `maketable` subcommands.

### Changed
- Repository moved to BeckLaboratory/pav3 and images deployed to GHCR.
- Updated Agglovar version and merge fields.
- Documentation fixes for Read the Docs and a callout about container permissions.

### Fixed
- `pav3.workflow.link_fasta()` returned an empty list for uncompressed FASTA files.

## 3.0.0.dev23 - 2026-06-18

### Added
- Docker and Apptainer build processes.

### Changed
- Routines locating PAV resources accept a path other than the working directory.
- Updated Polars and other dependency versions.

### Fixed
- Workaround for a Polars 1.41 streaming bug.

## 3.0.0.dev22 - 2026-05-06

### Added
- Functions converting alignment operations to sequence operations (`M` to `=`/`X`).

### Changed
- VCF writing supports multiple samples.
- Schema harmonized with Agglovar.

### Fixed
- Reading `M` alignment operations from SAM files.
- Intra-alignment variant performance.
- `from_definition` typo.

## 3.0.0.dev21 - 2026-03-30

### Fixed
- Off-by-one error lifting regions.
- Float32/Float64 mismatch with Polars 1.38 and newer.
- Snakemake `namedlist` import location varies by version.

## 3.0.0.dev20 - 2026-03-23

### Changed
- Genotype procedure supports multi-sample merges.
- Removed Snakemake from `rules/align.smk`.

### Fixed
- VCF irregularities.

## 3.0.0.dev19 - 2026-03-19

### Changed
- Tuned inter-alignment variant models for robustness.
- Score model accepts n template switches.
- `integrate_sources` moved from Snakemake into a library routine.

### Fixed
- Inter-alignment variant discovery and scoring model bugs.
- Interval reported the absolute value for `ref_len`, which is negative for a tandem duplication at a breakpoint.
- Float32/Float64 type mismatch between NumPy and Polars.

## 3.0.0.dev18 - 2026-02-18

### Fixed
- Complex variant structures, including duplications at a breakpoint.
- VCF failures with empty variant tables.
- Dependabot vulnerability alert.

## 3.0.0.dev17 - 2026-02-10

### Fixed
- Missing inversion schema for intra-alignment inversions with no variants.

## 3.0.0.dev16 - 2026-02-09

### Changed
- Schema cast sorts columns by schema order.

### Fixed
- Intra-alignment inversion discovery.
- Dependabot alerts (urllib3, filelock, virtualenv).

## 3.0.0.dev15 - 2026-02-08

### Added
- Variant and alignment tracks.
- Tools for interactive development and troubleshooting.
- Configuration file flexibility.

### Changed
- Reworked the workflow framework.
- Alignment lift framework rewritten for efficiency.

### Fixed
- Inversion detection, including false negatives around alignment breaks with no aligned inverted segment.
- Intra-alignment variant calling and coordinate lifts.
- Variant sort order at intermediate stages.
- Variant table schema mismatches.

## 3.0.0.dev14 - 2025-12-08

### Changed
- Semicolons are converted when writing VCFs.

## 3.0.0.dev13 - 2025-12-05

### Added
- VCF writer.

### Changed
- `call` subcommand renamed to `batch`.
- Configuration file renamed from `config.json` to `pav.json`, with a warning for the old name.

## 3.0.0.dev12 - 2025-11-27

### Fixed
- Inversion flag bugs.
- Errors from PAV and Agglovar versions getting out of sync.

## 3.0.0.dev11 - 2025-11-26

### Changed
- Updated Agglovar version.

## 3.0.0.dev10 - 2025-11-06

### Changed
- Updated Agglovar version.

## 3.0.0.dev8 - 2025-11-06

### Changed
- Variant ID expressions moved to Agglovar.

### Fixed
- QC routines causing inversion loss.
- Build not finding PAV files.

## 3.0.0.dev6, 3.0.0.dev7 - 2025-10-21

### Fixed
- Removed a line raising `NotImplementedError`.

## 3.0.0.dev5 - 2025-10-21

### Changed
- Haplotype variant tables relocated.

### Fixed
- Variant table handling.

## 3.0.0.dev4 - 2025-10-08

### Added
- Command-line interface.
- `CONTRIBUTING.md`.

### Changed
- Package data relocated and packaging updated for upload.

## 3.0.0.dev3 - 2025-10-06

### Added
- PAV 3 complete through haplotype call tables (rule `call_tables_hap`).

## 3.0.0.dev2 - 2025-06-23

### Added
- Low-confidence alignment prediction model for filtering alignment records.
- QC for large SV calling and tracks.

### Changed
- Alignment trimming performance for deep overlapping alignment records.
- Reworked feature handling and alignment operation handling.
- CIGAR and non-large inversion calls integrated with PAV 3.
- I/O avoids pysam, which breaks on large operation lengths.
- Renamed `.snakemake` files to `.smk`.

### Removed
- `keep-going` from the default Snakemake profiles.

### Fixed
- Large variant calls and alignment filter aliases.

## 3.0.0.dev1 - 2024-12-22

### Changed
- Altered PAV project structure including library.
- Replaced Pandas with Polars.
- Increased performance with Polars and vectorization.
- Repository transitioned to the Beck laboratory.

### Added
- Variant intersects are handled by agglovar.

### Removed
- SVPop was replaced by agglovar and some functions were moved into PAV.
