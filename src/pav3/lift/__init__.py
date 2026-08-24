"""Coordinate lifting and sequence extraction for the ``pav3 lift`` and ``pav3 subseq`` tools.

This package exposes the API used by the CLI: source-table builders, lift sources (PAV and Align modes), output
handlers, and the run entrypoints. It is distinct from :mod:`pav3.align.lift`, which provides the low-level alignment
lift engine (:class:`~pav3.align.lift.AlignLift`) that this package builds on.
"""

__all__ = [
    # Formats
    'Format',
    'parse_format',
    'detect_format',
    # Source tables
    'SOURCE_SCHEMA',
    'source_table_from_list',
    'source_table_from_str',
    'source_table_from_file',
    # Sources
    'LiftSource',
    'LiftSourcePav',
    'LiftSourceAlign',
    # Output schemas and handlers
    'LIFT_SCHEMA',
    'SUBSEQ_SCHEMA',
    'LiftOutput',
    'LiftOutputCsv',
    'LiftOutputBed',
    'LiftOutputParquet',
    'SubseqOutput',
    'SubseqOutputCsv',
    'SubseqOutputBed',
    'SubseqOutputFasta',
    'SubseqOutputParquet',
    # Entrypoints
    'lift_pav',
    'lift_align',
    'subseq_pav',
    'subseq_align',
]

from ..format import Format, detect_format, parse_format

from .source_table import (
    SOURCE_SCHEMA,
    source_table_from_file,
    source_table_from_list,
    source_table_from_str,
)

from .source import (
    LiftSource,
    LiftSourceAlign,
    LiftSourcePav,
)

from .output import (
    LIFT_SCHEMA,
    SUBSEQ_SCHEMA,
    LiftOutput,
    LiftOutputBed,
    LiftOutputCsv,
    LiftOutputParquet,
    SubseqOutput,
    SubseqOutputBed,
    SubseqOutputCsv,
    SubseqOutputFasta,
    SubseqOutputParquet,
)

from .run import (
    lift_align,
    lift_pav,
    subseq_align,
    subseq_pav,
)
