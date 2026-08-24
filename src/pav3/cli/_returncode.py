"""CLI return codes.

Shared exit codes for PAV subcommands and the guard that maps exceptions onto them. ``USAGE`` matches argparse's exit
code for command-line usage errors.
"""

import enum
import logging
from typing import Callable

from ..errors import (
    PavConfigError,
    PavError,
    PavInputError,
    PavOutputError,
    PavUsageError,
)

__all__ = [
    'ExitCode',
    'run_guarded',
]

logger = logging.getLogger(__name__)


class ExitCode(enum.IntEnum):
    """PAV CLI exit codes."""

    OK = 0
    """Success."""

    ERROR = 1
    """General or uncategorized error."""

    USAGE = 2
    """Command-line usage error (matches argparse)."""

    CONFIG = 3
    """Bad, missing, or ambiguous assembly table or PAV configuration."""

    INPUT = 4
    """Missing or unreadable input."""

    OUTPUT = 5
    """Failure writing an output."""


def run_guarded(
        action: Callable[[], None]
) -> int:
    """Run an action, mapping exceptions to exit codes.

    :param action: Zero-argument callable that performs the operation.

    :returns: An integer exit code.
    """
    try:
        action()
    except PavUsageError as e:
        logger.error('%s', e)
        return int(ExitCode.USAGE)
    except PavConfigError as e:
        logger.error('%s', e)
        return int(ExitCode.CONFIG)
    except PavOutputError as e:
        logger.error('%s', e)
        return int(ExitCode.OUTPUT)
    except (PavInputError, FileNotFoundError) as e:
        logger.error('%s', e)
        return int(ExitCode.INPUT)
    except PavError as e:
        logger.error('%s', e)
        return int(ExitCode.ERROR)
    except ValueError as e:
        logger.error('%s', e)
        return int(ExitCode.ERROR)
    except Exception:
        logger.exception('Unexpected error')
        return int(ExitCode.ERROR)

    return int(ExitCode.OK)
