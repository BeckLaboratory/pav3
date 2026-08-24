"""Exceptions for the lift and subseq tools, mapped by the CLI to exit codes.

These names are aliases for the PAV exceptions in :mod:`pav3.errors`, which the lift and subseq tools predate.
"""

__all__ = [
    'LiftError',
    'LiftUsageError',
    'LiftConfigError',
    'LiftInputError',
    'LiftOutputError',
]

from ..errors import (
    PavConfigError,
    PavError,
    PavInputError,
    PavOutputError,
    PavUsageError,
)

LiftError = PavError
"""Base class for lift/subseq errors."""

LiftUsageError = PavUsageError
"""A command-line usage error (maps to exit code USAGE)."""

LiftConfigError = PavConfigError
"""A bad, missing, or ambiguous PAV configuration or assembly table (maps to exit code CONFIG)."""

LiftInputError = PavInputError
"""A missing or unreadable input (maps to exit code INPUT)."""

LiftOutputError = PavOutputError
"""A failure writing an output (maps to exit code OUTPUT)."""
