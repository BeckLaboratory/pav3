"""PAV exceptions mapped by the CLI to exit codes.

Library routines raise these to signal how a failure should be reported. :func:`pav3.cli.run_guarded` maps each to an
:class:`~pav3.cli.ExitCode`; any other exception is reported as an unexpected error.
"""

__all__ = [
    'PavError',
    'PavUsageError',
    'PavConfigError',
    'PavInputError',
    'PavOutputError',
]


class PavError(Exception):
    """Base class for PAV errors."""


class PavUsageError(PavError):
    """A command-line usage error (maps to exit code USAGE)."""


class PavConfigError(PavError):
    """A bad, missing, or ambiguous PAV configuration or assembly table (maps to exit code CONFIG)."""


class PavInputError(PavError):
    """A missing or unreadable input (maps to exit code INPUT)."""


class PavOutputError(PavError):
    """A failure writing an output (maps to exit code OUTPUT)."""
