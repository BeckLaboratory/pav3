"""Maketable subcommand.

Convert a Parquet table to a BED, CSV, TSV, or Excel file. Designed for PAV variant and alignment tables, but works
with any Parquet file.
"""

import argparse
import logging
from typing import Optional

from ..maketable import maketable
from ..util import init_logger
from ._common_opt import _add_opt_verbose, _add_opt_version
from ._returncode import run_guarded

logger = logging.getLogger(__name__)


def subcommand_maketable(
        input: str,
        output: str,
        output_type: Optional[str] = None,
        sort: Optional[bool] = None,
        verbose: bool = False,
) -> int:
    """Run the maketable subcommand.

    :param input: Input Parquet file.
    :param output: Output file.
    :param output_type: Output type override (BED, CSV, TSV, XLSX; may append "GZ"). Detected from the output
        filename if `None`.
    :param sort: Sort by the first three columns (chrom, pos, end). Defaults to `True` for BED output and `False`
        for all other output types if `None`.
    :param verbose: Enable DEBUG logging.

    :returns: Exit code.
    """
    init_logger(logging.DEBUG if verbose else logging.INFO, force=True)

    def action() -> None:
        maketable(input, output, output_type=output_type, sort=sort)

    return run_guarded(action)


def _add_subparser_maketable(
        subparsers,
) -> argparse.ArgumentParser:
    """Add the maketable subcommand to a subparser.

    :param subparsers: Subparser object.

    :returns: Configured subparser.
    """
    parser_maketable = subparsers.add_parser(
        'maketable',
        description=(
            'Convert a Parquet table to a BED, CSV, TSV, or Excel file. The output type is detected from the output '
            'filename unless -O/--output-type is given.'
        ),
        help='Convert a Parquet table to BED, CSV, TSV, or Excel.',
    )

    _add_opt_version(parser_maketable)

    parser_maketable.add_argument(
        '--input', '-i',
        required=True,
        help='Input Parquet file.',
    )

    parser_maketable.add_argument(
        '--output', '-o',
        required=True,
        help='Output file.',
    )

    parser_maketable.add_argument(
        '--output-type', '-O',
        default=None, dest='output_type',
        help=(
            'Output type (BED, BEDGZ, CSV, CSVGZ, TSV, TSVGZ, XLSX or EXCEL; not case sensitive). Overrides the type '
            'detected from the output filename. Compressed types ("GZ") are bgzip-compressed.'
        ),
    )

    sort_group = parser_maketable.add_mutually_exclusive_group()

    sort_group.add_argument(
        '--sort', '-s',
        default=None, dest='sort', action='store_true',
        help=(
            'Sort the output by the first three columns (chrom, pos, end), which must be a string followed by two '
            'integers. Default for BED output.'
        ),
    )

    sort_group.add_argument(
        '--no-sort', '-S',
        default=None, dest='sort', action='store_false',
        help='Do not sort the output. Default for output types other than BED.',
    )

    _add_opt_verbose(parser_maketable)

    return parser_maketable
