<h1 align="center"><img width="300px" src="img/logo/PAVLogo_Full.png"/></h1>
<p align="center">Phased Assembly Variant Caller</p>

***
<!-- Templated header from the pbsv github page: https://github.com/PacificBiosciences/pbsv -->

Variant caller for assembled genomes.


## Install

PAV runs from a container or from a Python environment.

### Containers

Public container images include PAV and all of its dependencies.

Docker:

```bash
docker pull ghcr.io/becklaboratory/pav3:latest

docker run --rm -v "$PWD:$PWD" -w "$PWD" --user "$(id -u):$(id -g)" \
    ghcr.io/becklaboratory/pav3:latest batch
```

Apptainer (or Singularity):

```bash
apptainer pull pav3.sif oras://ghcr.io/becklaboratory/pav3-sif:latest

apptainer run pav3.sif batch
```

The container runs `pav3`, so arguments after the image are the same as on the command line (`batch`, `--help`, and
so on).

Both examples give PAV access to the current directory only. If the reference or the assemblies are in another
directory, add `-v /path:/path` (Docker) or `--bind /path` (Apptainer), replacing `/path` with that directory.

Replace `latest` with a version (e.g. `3.0.0`) to pin a release.

### Install locally

```bash
pip install pav3
```

PAV requires Python 3.12 or newer and `minimap2` in the `PATH`. The `tracks` targets also require `bedToBigBed`
([UCSC executable downloads](https://hgdownload.soe.ucsc.edu/admin/exe/)). The containers include both. All other
dependencies are Python packages installed with PAV.


## Configuration

PAV runs in a directory with a configuration file (`pav.yaml`) and an assembly table (`assemblies.csv`), and writes
its output to the same directory.

### Configuration file

`pav.yaml` sets PAV options. The only required option is `reference`, the path to the reference FASTA file. All other
options have defaults.

```yaml
reference: /path/to/hg38.no_alt.fa.gz
```

PAV also reads `pav.yml`, or `pav.json` if no YAML file exists. If both a YAML and a JSON file exist, PAV prints a
warning and uses the YAML file. `pav3 batch --config FILE` reads the configuration from another path.

### Reference genome

Use a primary assembly of chromosome scaffolds, including unplaced and unlocalized contigs if it has them (GRCh/hg
references). Do not use references with ALT, PATCH, or decoy scaffolds. The redundant sequence in these contigs may
increase callset errors.

Suitable references:
* GRCh38 no-ALT reference shared by the HGSVC:
  https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/technical/reference/20200513_hg38_NoALT/
* T2T-CHM13v2.0 (UCSC hs1), unaltered.

### Assembly table

The assembly table lists the input assemblies, one row per assembly. PAV uses the first of these files it finds in
the run directory:

1. `assemblies.csv`
2. `assemblies.xlsx`
3. `assemblies.tsv`

To use another file, set its path with the `assembly_table` option in `pav.yaml`. It may be CSV, Excel, TSV, or
Parquet, and CSV and TSV files may be gzipped.

Columns:
* `name`: Assembly or sample name.
* `hap_*`: One column per haplotype, named for the haplotype (e.g. `hap_h1` or `hap_mat`).
* `config` (optional): Configuration options for this assembly (see
  [Per-assembly configuration](#per-assembly-configuration)).

Example with two haplotypes:

| name    | hap_h1                   | hap_h2                   |
|---------|--------------------------|--------------------------|
| sample1 | /path/to/sample1_hap1.fa | /path/to/sample1_hap2.fa |
| sample2 | /path/to/sample2_hap1.fa | /path/to/sample2_hap2.fa |

Haplotype names are labels and do not change how input is processed. An assembly may have any number of haplotypes,
and a blank entry means the assembly does not have that haplotype. Each VCF genotype has one allele per haplotype with
an input path, in table column order. With the table above, `0|1` means the variant is in `h2` but not `h1`.

### Haplotype input formats

Each `hap_*` entry is a path to a FASTA, FASTQ, GFA, or FOFN (file of file names) file. Separate multiple files with
semicolons (e.g. `file1.fa.gz;file2.fa.gz`); file types may be mixed. A FOFN may list files of any of these types,
including other FOFNs (nesting is not recommended, but cycles are detected and broken).

PAV is fastest with a single indexed FASTA file, uncompressed or bgzipped. For any other input, PAV first writes and
indexes a FASTA file for the haplotype, which adds to the run time.

### Per-assembly configuration

The optional `config` column sets options for one assembly, overriding `pav.yaml` and PAV defaults. It holds
semicolon-separated `key=value` pairs (e.g. `min_mapq=30;lc_model=default`). `reference` can only be set in
`pav.yaml`.


## Run

From the directory with `pav.yaml` and the assembly table:

```bash
pav3 batch
```

If the `pav3` command is not on the `PATH`, run `python -m pav3 batch`.

Use `-c N` to limit PAV to N cores, `-n` for a dry run, and `--profile DIR` to pass a Snakemake profile (e.g. for
cluster execution). `pav3 batch --help` lists all options.


## Output

For each assembly `NAME`:
* `NAME.vcf.gz`: Variant calls, indexed (`NAME.vcf.gz.csi`).
* `results/NAME/call/`: Variant tables merged across haplotypes.
* `results/NAME/call_hap/`: Variant tables for each haplotype, and callable regions in reference (`callable_ref`)
  and assembly (`callable_qry`) coordinates.
* `results/NAME/align/`: Alignment tables.

Tables are in Parquet format (`pav3 maketable` converts them to BED, CSV, TSV, or Excel).


## Targets

`pav3 batch` builds the VCF for every assembly. To build other outputs or stop at an earlier step, name a target:

```bash
pav3 batch tracks
```

* `vcf` (default): VCF for each assembly (`NAME.vcf.gz`).
* `call_tables`: Variant tables merged across haplotypes (`results/NAME/call/`). Stops before writing the VCF.
* `call_tables_hap`: Variant tables for each haplotype (`results/NAME/call_hap/`). Stops before merging haplotypes.
* `callable`: Callable regions for each haplotype in reference and assembly coordinates
  (`results/NAME/call_hap/callable_*`).
* `align`: Alignment tables (`results/NAME/align/`). Stops before calling variants.
* `data_init`: Reference files shared by all assemblies (`data/ref/`).
* `tracks`: All UCSC tracks (`tracks_var` and `tracks_align`). Requires `bedToBigBed`.
* `tracks_var`: Variant tracks for each haplotype (`tracks/call_hap/`).
* `tracks_align`: Alignment tracks (`tracks/align/`).

A target can also be an output file. To run each assembly as an independent job, run `data_init` once, then run one
job per assembly in the same directory. These jobs can run at the same time without locking each other:

```bash
pav3 batch data_init

pav3 batch sample1.vcf.gz
pav3 batch sample2.vcf.gz
```

Without `data_init` first, each job tries to build the shared reference files, and Snakemake stops the later jobs
with a "Directory cannot be locked" error.


## Other commands

* `pav3 lift`: Lift regions from reference to assembly coordinates, or from assembly to reference coordinates with
  `-r`. Uses the alignments in a PAV run directory (`-d`, default: the current directory; requires the `align`
  target) or a SAM, BAM, or CRAM file (`-a`).
* `pav3 subseq`: Lift regions as `pav3 lift` does and output their sequences: assembly sequences by default, or
  reference sequences with `-r`.
* `pav3 maketable`: Convert a PAV Parquet table to BED, CSV, TSV, or Excel.
* `pav3 license`: Print the PAV license (`brief`, `full`, or `spdx`).
* `pav3 call`: Deprecated. Prints a warning and runs `pav3 batch`.

Run `pav3 COMMAND --help` for the options of each command.


## PAV versions

PAV follows Python's packaging version scheme (https://packaging.python.org/en/latest/discussions/versioning/) with
three fields:

- Major: Large changes or new features.
- Minor: Small changes that may affect the PAV API or command-line interface.
- Patch: Bug fixes and small features. Patch releases do not break API or command-line compatibility, but may change
  results.

Pre-release versions have a suffix: development (".devN"), alpha ("aN"), beta ("bN"), or release candidate ("rcN"),
where "N" is an integer greater than 0. For example, "3.0.0.dev1" is a development version before 3.0.0.


## Reporting bugs

Please report problems on the [issues page](https://github.com/BeckLaboratory/pav3/issues).

Thank you for your bug reports and feedback. They really do help improve PAV!


## Cite PAV

Cite PAV 3 with the bioRxiv preprint:

> Peter Audano, Carolyn Paisie, The Human Genome Structural Variation Consortium (HGSVC), and Christine Beck. 2024. “Large Complex Structural Rearrangements in Human Genomes Harbor Cryptic Structures.” In bioRxiv. December 22. https://doi.org/10.1101/2024.12.19.629504.


### PAV 1 & 2

Cite earlier versions of PAV as:

> Ebert et al., “Haplotype-Resolved Diverse Human Genomes and Integrated Analysis of Structural Variation,”
> Science, February 25, 2021, eabf7117, https://doi.org/10.1126/science.abf7117 (PMID: 33632895).


## License

PAV is released under the MIT License. See [LICENSE](LICENSE).
