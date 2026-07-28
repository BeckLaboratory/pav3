"""Common options.

Common options are defined here. This allows each subcommand to explicitly add common options it uses without relying
on global parent parsers. Some subcommands will use common options that others will not.
"""

import argparse

from .. import __version__


def _add_opt_version(parser: argparse.ArgumentParser) -> None:
    """Add version option to parser."""
    parser.add_argument(
        '--version',
        action='version',
        version=f'PAV {__version__}',
        help='Show PAV version and exit.'
    )


def _add_opt_verbose(
        parser: argparse.ArgumentParser
) -> None:
    """Add verbose option to parser."""
    parser.add_argument(
        '--verbose', '-v',
        default=False, action='store_true',
        help='Generate verbose output',
    )


def _add_opt_debug(
        parser: argparse.ArgumentParser
) -> None:
    """Add debug option to parser."""
    parser.add_argument(
        '--debug',
        default=False, action='store_true',
        help='Generate very verbose debugging output'
    )


def _add_opt_lift_common(
        parser: argparse.ArgumentParser
) -> None:
    """Add options common to the "lift" and "subseq" subcommands."""
    parser.add_argument(
        'coordinates',
        nargs='*',
        help='Regions as "chrom:pos-end" (1-based, commas allowed). Mutually exclusive with -i/--input.',
    )

    parser.add_argument(
        '--sample', '-s',
        default=None,
        help='Sample name. If omitted in PAV mode, all samples are processed.',
    )

    parser.add_argument(
        '--hap', '-p',
        default=None,
        help='Haplotype name. If omitted in PAV mode, all haplotypes are processed.',
    )

    parser.add_argument(
        '--input', '-i',
        default=None,
        help='Input coordinate file. Defaults to standard input if no coordinates are given.',
    )

    parser.add_argument(
        '--input-format', '-I',
        default=None, dest='input_format',
        help='Input format (TXT, BED, CSV, TSV, PARQUET; may append "GZ"). Overrides filename detection.',
    )

    parser.add_argument(
        '--output', '-o',
        default=None,
        help='Output file. Defaults to standard output.',
    )

    parser.add_argument(
        '--output-format', '-O',
        default=None, dest='output_format',
        help='Output format. Overrides filename detection.',
    )

    parser.add_argument(
        '--trim', '-m',
        default='qry', choices=['none', 'qry', 'qryref'],
        help='Alignment trim level (PAV mode only). Ignored in Align mode. Default: qry.',
    )

    parser.add_argument(
        '--reverse', '-r',
        default=False, action='store_true',
        help='Reverse: input is query coordinates, lift/extract to reference.',
    )

    parser.add_argument(
        '--pav-directory', '-d',
        default=None, dest='pav_directory',
        help='PAV run directory (default: current directory). Error if given in Align mode.',
    )

    parser.add_argument(
        '--allow-split', '-l',
        default=False, action='store_true', dest='allow_split',
        help='Allow pos and end to lift to different alignment records (same query ID and strand).',
    )

    parser.add_argument(
        '--align-file', '-a',
        default=None, dest='align_file',
        help='Alignment file (SAM/BAM/CRAM). Enables Align mode. Error if combined with -d/--pav-directory.',
    )

    parser.add_argument(
        '--verbose', '-v',
        default=False, action='store_true',
        help='Generate verbose (DEBUG) output.',
    )


def _add_opt_reference(
        parser: argparse.ArgumentParser
) -> None:
    """Add the reference option, common to "lift" and "subseq"."""
    parser.add_argument(
        '--reference', '-f',
        default=None,
        help=(
            'Reference FASTA. Required to read a CRAM alignment file (lift and subseq), and to extract reference '
            'sequences with -r/--reverse in subseq. Ignored in PAV mode.'
        ),
    )


def _add_opt_subseq_extra(
        parser: argparse.ArgumentParser
) -> None:
    """Add options specific to the "subseq" subcommand (excluding -f/--reference)."""
    parser.add_argument(
        '--query-file', '-q',
        default=None, dest='query_file',
        help='Query FASTA/FASTQ file (Align mode). Requires -a/--align-file.',
    )

    parser.add_argument(
        '--query-format', '-Q',
        default=None, dest='query_format',
        help='Query file format (FASTA or FASTQ; may append "GZ"). Overrides filename detection.',
    )

    parser.add_argument(
        '--upper', '-u',
        default=False, action='store_true',
        help='Uppercase all extracted sequences.',
    )
