"""Exceptions for the lift and subseq tools, mapped by the CLI to exit codes."""

__all__ = [
    'LiftError',
    'LiftUsageError',
    'LiftConfigError',
    'LiftInputError',
    'LiftOutputError',
]


class LiftError(Exception):
    """Base class for lift/subseq errors."""


class LiftUsageError(LiftError):
    """A command-line usage error (maps to exit code USAGE)."""


class LiftConfigError(LiftError):
    """A bad, missing, or ambiguous PAV configuration or assembly table (maps to exit code CONFIG)."""


class LiftInputError(LiftError):
    """A missing or unreadable input (maps to exit code INPUT)."""


class LiftOutputError(LiftError):
    """A failure writing an output (maps to exit code OUTPUT)."""
