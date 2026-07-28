"""Lift subcommand.

Lift a set of genomic regions through assembly alignments (reference -> query by default, query -> reference with
-r/--reverse). Runs in PAV mode (a PAV run directory) or Align mode (-a/--align-file).
"""

import argparse
import logging
from typing import Optional

from ..lift import lift_align, lift_pav
from ..util import init_logger
from ._common_opt import _add_opt_lift_common, _add_opt_reference, _add_opt_version
from ._lift_cli import build_source_table, resolve_output, run_guarded
from ._returncode import ExitCode

logger = logging.getLogger(__name__)


def subcommand_lift(
        coordinates: Optional[list[str]] = None,
        sample: Optional[str] = None,
        hap: Optional[str] = None,
        input: Optional[str] = None,
        input_format: Optional[str] = None,
        output: Optional[str] = None,
        output_format: Optional[str] = None,
        trim: str = 'qry',
        reverse: bool = False,
        pav_directory: Optional[str] = None,
        allow_split: bool = False,
        align_file: Optional[str] = None,
        reference: Optional[str] = None,
        verbose: bool = False,
) -> int:
    """Run the lift subcommand.

    :param coordinates: Positional region strings.
    :param sample: Sample name (PAV mode expansion if `None`).
    :param hap: Haplotype name (PAV mode expansion if `None`).
    :param input: Input coordinate file (standard input if `None` and no coordinates).
    :param input_format: Input format override.
    :param output: Output file (standard output if `None`).
    :param output_format: Output format override.
    :param trim: Alignment trim level (PAV mode).
    :param reverse: Lift query -> reference if `True`.
    :param pav_directory: PAV run directory (PAV mode).
    :param allow_split: Allow start/end to lift to different alignment records.
    :param align_file: Alignment file (enables Align mode).
    :param reference: Reference FASTA (required to read a CRAM in Align mode; ignored in PAV mode).
    :param verbose: Enable DEBUG logging.

    :returns: Exit code.
    """
    init_logger(logging.DEBUG if verbose else logging.INFO, force=True)

    if align_file is not None and pav_directory is not None:
        logger.error('Cannot combine -a/--align-file with -d/--pav-directory (Align mode uses no PAV directory)')
        return int(ExitCode.USAGE)

    def action() -> None:
        source_table = build_source_table(coordinates, input, input_format)
        out_handler = resolve_output(False, output, output_format, reverse)

        if align_file is not None:
            lift_align(
                source_table, out_handler,
                align_file=align_file, ref_file=reference, sample=sample, hap=hap,
                reverse=reverse, allow_split=allow_split,
            )
        else:
            lift_pav(
                source_table, out_handler,
                pav_dir=pav_directory, sample=sample, hap=hap, trim=trim,
                reverse=reverse, allow_split=allow_split, ref_file=reference,
            )

    return run_guarded(action)


def _add_subparser_lift(
        subparsers,
) -> argparse.ArgumentParser:
    """Add the lift subcommand to a subparser.

    :param subparsers: Subparser object.

    :returns: Configured subparser.
    """
    parser_lift = subparsers.add_parser(
        'lift',
        description='Lift genomic regions through assembly alignments.',
        help='Lift regions between reference and query coordinates.',
    )

    _add_opt_version(parser_lift)
    _add_opt_lift_common(parser_lift)
    _add_opt_reference(parser_lift)

    return parser_lift
