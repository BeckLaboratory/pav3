"""CLI return codes.

Shared exit codes for PAV subcommands. ``USAGE`` matches argparse's exit code for command-line usage errors.
"""

import enum

__all__ = [
    'ExitCode',
]


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
