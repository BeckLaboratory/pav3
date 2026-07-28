"""Shared CLI helpers for the ``lift`` and ``subseq`` subcommands.

These helpers build the source coordinate table, resolve the output handler, and map exceptions to exit codes, keeping
the subcommand modules thin wrappers over the :mod:`pav3.lift` API.
"""

import logging
import sys
from typing import Callable, Optional

import polars as pl

from ..lift import (
    Format,
    LiftOutputBed,
    LiftOutputCsv,
    LiftOutputParquet,
    SubseqOutputBed,
    SubseqOutputCsv,
    SubseqOutputFasta,
    SubseqOutputParquet,
    detect_format,
    parse_format,
    source_table_from_file,
    source_table_from_list,
    source_table_from_str,
)
from ..lift.errors import LiftConfigError, LiftError, LiftInputError, LiftOutputError, LiftUsageError
from ._returncode import ExitCode

logger = logging.getLogger(__name__)

_LIFT_OUTPUT_FORMATS = {Format.TSV, Format.CSV, Format.BED, Format.PARQUET}
_SUBSEQ_OUTPUT_FORMATS = _LIFT_OUTPUT_FORMATS | {Format.FASTA}


def build_source_table(
        coordinates: Optional[list[str]],
        input_file: Optional[str],
        input_format: Optional[str],
) -> pl.DataFrame:
    """Build a source coordinate table from the command line, a file, or standard input.

    :param coordinates: Positional region strings (mutually exclusive with `input_file`).
    :param input_file: Input coordinate file, or `None`.
    :param input_format: Input format string (from -I), or `None` to auto-detect.

    :returns: Source coordinate table.

    :raises LiftUsageError: If both coordinates and an input file are given.
    """
    if coordinates and input_file:
        raise LiftUsageError('Cannot provide both coordinates on the command line and an input file (-i/--input)')

    if coordinates:
        return source_table_from_list(coordinates)

    if input_file:
        return source_table_from_file(input_file, input_format)

    logger.info('Reading regions from standard input')

    return source_table_from_str(sys.stdin.read(), input_format)


def _resolve_format(out_file: Optional[str], out_format: Optional[str]) -> tuple[Format, bool]:
    """Resolve the output format and compression from the format option or filename.

    :param out_file: Output file path, or `None` for stdout.
    :param out_format: Output format string (from -O), or `None`.

    :returns: A tuple of (format, is_gz). Defaults to BED for stdout when neither is given.

    :raises LiftUsageError: If the format cannot be determined.
    """
    if out_format is not None:
        try:
            return parse_format(out_format)
        except ValueError as e:
            raise LiftUsageError(str(e)) from e

    if out_file is not None:
        try:
            detected = detect_format(out_file)
        except ValueError as e:
            raise LiftUsageError(str(e)) from e

        if detected is None:
            raise LiftUsageError(f'Cannot determine output format from filename: {out_file}')

        return detected

    return Format.BED, False


def resolve_output(
        is_subseq: bool,
        out_file: Optional[str],
        out_format: Optional[str],
        reverse: bool,
):
    """Build an output handler for the lift or subseq operation.

    :param is_subseq: `True` for a subseq handler, `False` for a lift handler.
    :param out_file: Output file path, or `None` for stdout.
    :param out_format: Output format string (from -O), or `None` to auto-detect.
    :param reverse: Matches the `reverse` used for the lift.

    :returns: An output handler (context manager).

    :raises LiftUsageError: If the format is invalid for the operation or the destination is unsupported.
    """
    fmt, is_gz = _resolve_format(out_file, out_format)

    valid = _SUBSEQ_OUTPUT_FORMATS if is_subseq else _LIFT_OUTPUT_FORMATS

    if fmt not in valid:
        raise LiftUsageError(
            f'Invalid output format for {"subseq" if is_subseq else "lift"}: {fmt.value}'
        )

    try:
        if fmt == Format.PARQUET:
            handler_cls = SubseqOutputParquet if is_subseq else LiftOutputParquet
            return handler_cls(out_file=out_file, reverse=reverse)

        if fmt == Format.BED:
            handler_cls = SubseqOutputBed if is_subseq else LiftOutputBed
            return handler_cls(out_file=out_file, reverse=reverse, is_gz=is_gz)

        if fmt in (Format.CSV, Format.TSV):
            separator = ',' if fmt == Format.CSV else '\t'
            handler_cls = SubseqOutputCsv if is_subseq else LiftOutputCsv
            return handler_cls(out_file=out_file, reverse=reverse, is_gz=is_gz, separator=separator)

        # FASTA (subseq only)
        return SubseqOutputFasta(out_file=out_file, reverse=reverse, is_gz=is_gz)

    except ValueError as e:
        raise LiftUsageError(str(e)) from e


def run_guarded(action: Callable[[], None]) -> int:
    """Run an action, mapping exceptions to exit codes.

    :param action: Zero-argument callable that performs the operation.

    :returns: An integer exit code.
    """
    try:
        action()
    except LiftUsageError as e:
        logger.error('%s', e)
        return int(ExitCode.USAGE)
    except LiftConfigError as e:
        logger.error('%s', e)
        return int(ExitCode.CONFIG)
    except LiftOutputError as e:
        logger.error('%s', e)
        return int(ExitCode.OUTPUT)
    except (LiftInputError, FileNotFoundError) as e:
        logger.error('%s', e)
        return int(ExitCode.INPUT)
    except LiftError as e:
        logger.error('%s', e)
        return int(ExitCode.ERROR)
    except ValueError as e:
        logger.error('%s', e)
        return int(ExitCode.ERROR)
    except Exception:
        logger.exception('Unexpected error')
        return int(ExitCode.ERROR)

    return int(ExitCode.OK)
