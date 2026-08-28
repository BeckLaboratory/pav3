"""PAV command-line interface (CLI), argument parsing, and subcommand routines."""

__all__ = [
    'build_parser',
    'parse_arguments',
    'main',
    'subcommand_batch',
    'subcommand_call',
    'subcommand_lift',
    'subcommand_subseq',
    'subcommand_maketable',
    'subcommand_license',
]

from ._cli import (
    build_parser, parse_arguments, main,
)

from ._subcommand_call import subcommand_call

from ._subcommand_batch import subcommand_batch

from ._subcommand_lift import subcommand_lift

from ._subcommand_subseq import subcommand_subseq

from ._subcommand_maketable import subcommand_maketable

from ._subcommand_license import subcommand_license
