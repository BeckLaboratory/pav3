<h1 align="center"><img width="300px" src="img/logo/PAVLogo_Full.png"/></h1>
<p align="center">Phased Assembly Variant Caller</p>

***
<!-- Templated header from the pbsv github page: https://github.com/PacificBiosciences/pbsv -->

Variant caller for assembled genomes.


## Install

PAV can be run from a container or integrated into a Python or conda environment.

### Containers

Public container images include PAV and all dependencies.

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

The container runs `pav3` itself, so arguments follow the image exactly as they would on the command line (`batch`, `--help`, and so on).

Both examples give PAV access to the current directory only. If the reference genome or the assemblies are somewhere
else, add `-v /path:/path` for Docker or `--bind /path` for Apptainer
and change "/path" to the actual path on the host machine.

Replace `latest` with a version (e.g. `3.0.0.dev24`) to pin a specific release.


### Install locally

Alternatively, PAV can be installed directly from PyPI into existing Python or conda environments without requiring a container.

```bash
pip install pav3
```

Python 3.12 or newer is required. See [Dependencies](#dependencies) below.

`minimap2` will need to be available in the environment where PAV is run.


## Configuration

PAV currenly supports batch processing of assemblies (`pav batch` command). A configuration file and a table of assemlbies are needed.

### Configuration file

A YAML configuration file, `pav.yaml`, configures PAV. Default options are built-in, and the only required option is
`reference` pointing to a reference FASTA file.

Example:

```
reference: /path/to/hg38.no_alt.fa.gz
```

Older versions read `pav.json`, which PAV will read if no `pav.yaml` file is found. If both are present, a warning is
printed and `pav.yaml` is used.

### Assembly table

The assembly table points to input assemblies. It may be in TSV, CSV, Excel, or parquet formats (TSV and CSV may
optionally be gzipped). Each assembly has one row in the table.

Columns:
* name: Assembly or sample name.
* hap_\*: One column for each assembled haplotype.

For example, a table with two haplotypes might look like:

| name    | hap_h1                   | hap_h2                   |
|---------|--------------------------|--------------------------|
| sample1 | /path/to/sample1_hap1.fa | /path/to/sample1_hap2.fa |
| sample2 | /path/to/sample2_hap1.fa | /path/to/sample2_hap2.fa |

Haplotype names can be anything, such as "pat" and "mat" (columns `hap_pat` and `hap_mat`). Haplotype names do not change how input is handled.

The output VCF will have one allele per defined haplotype in the order they are defined in the table. Using the above table, "0|1" means the the variant is in "h2", but not "h1". Missing haplotypes are acceptable and blank entries are ignored. Expected ploidy per sample is the number of haplotypes with defined paths. Any number of haplotypes may be input for each sample.


## Run

To run PAV:

```bash
pav3 batch
```

Some Python environments may require you to run `pav3` through the `python` command:

```bash
python -m pav3 batch
```

## Output
PAV will output a VCF file `NAME.vcf.gz` to the run directory for each sample.

* `results/NAME/call_hap`: Unmerged variant calls.
  * Includes tables of callable regions in reference ("callable_ref") and query ("callable_qry") coordinates.
* `results/NAME/call`: Variant calls merged across haplotypes.

All tables are in parquet format.

## Additional information

### Dependencies

Currently, PAV needs `minimap2` in the environment where it is run (the containers above already include it). This may
change in future releases.

Optionally, PAV can build UCSC tracks, which requires `bedToBigBed` in the environment (see [UCSC executable downloads](http://hgdownload.soe.ucsc.edu/admin/exe/)).

All other dependencies are PyPI packages and are installed automatically with `pav3`.

unique value.


#### Haplotype input formats

Each "hap_" column contains paths to input files in FASTA, FASTQ, GFA, or FOFN format. FOFN may contain paths to these
same file types including other FOFNs (recursive FOFNs are not recommended, although cycles are detected and broken). Multiple files
can be input by separating them by semi-colons (i.e. "path/to/file1.fasta.gz;path/to/file2.fasta.gz") and a mix of
file types is possible.

PAV will be fastest if the input is a single indexed FASTA file, which may be plain text or bgzipped (Samtools bgzip).
For all other input including multiple FASTA files or any GFA files for one haplotype, PAV will build a single FASTA
file from the input and index it, which adds to run time.


#### Configuration column for global overrides

An optional "config" column in the assembly table can override global configuration parameters per assembly. Global configuration parameters
are defined in `pav.yaml` or are PAV default values if not defined. Values in this column are semicolon-separated lists
of key-value pairs (i.e. "key1=val1;key2=val2"). The "reference" parameter can only be set in `pav.yaml`.


## Suitable genome references

Align against primary assemblies composed of chromosome scaffolds including unplaced contigs and unlocalized contigs
if they exist (GRCh/hg references). Do not use references with ALT, PATCH, or DECOY scaffolds; reference redundancy
from these non-primary contigs may increase callset errors.

The GRCh38 HGSVC no-ALT reference for long reads was shared by the HGSVC consortium:
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/technical/reference/20200513_hg38_NoALT/

The T2T-CHM13v2.0 (hs1 on UCSC) is suitable without alteration.


## PAV versions

PAV uses Python package versioning with three fields:

- Major: Major changes or new features.
- Minor: Small changes, but may affect PAV's API or command-line interfaces.
- Patch: Small changes and minor new features. Patch should do not break API or command-line compatibility, but may
  add new features or fix bugs that alter results.

PAV follows Python's packaging versioning scheme (https://packaging.python.org/en/latest/discussions/versioning/).

PAV may use pre-release versions with a suffix for development releases (".devN"), alpha ("aN"), beta ("bN"), or
release-candidate ("rcN") where "N" is an integer greater than 0. For example, "3.0.0.dev1" is a pre-3.0.0 development version.


## Reporting bugs

Please report problems you encounter on the
[issues page](https://github.com/BeckLaboratory/pav3/issues).

Thank you for your bug reports and feedback, it really does help improve PAV!


## Cite PAV

PAV 3 can be cited using the bioRxiv preprint:

[https://www.biorxiv.org/content/10.1101/2024.12.19.629504v1](https://www.biorxiv.org/content/10.1101/2024.12.19.629504v1).

> Peter Audano, Carolyn Paisie, The Human Genome Structural Variation Consortium (HGSVC), and Christine Beck. 2024. “Large Complex Structural Rearrangements in Human Genomes Harbor Cryptic Structures.” In bioRxiv. December 22. https://doi.org/10.1101/2024.12.19.629504.


### PAV 1 & 2

Earlier versions of PAV can be cited as:

> Ebert et al., “Haplotype-Resolved Diverse Human Genomes and Integrated Analysis of Structural Variation,”
> Science, February 25, 2021, eabf7117, https://doi.org/10.1126/science.abf7117 (PMID: 33632895).


## License

PAV is released under the MIT License. See [LICENSE](LICENSE).
