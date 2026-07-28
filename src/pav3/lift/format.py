"""File format detection and parsing for the lift and subseq tools."""

__all__ = [
    'Format',
    'parse_format',
    'detect_format',
    'TEXT_TABLE_FORMATS',
    'SEQ_FORMATS',
]

import enum
from pathlib import Path
from typing import Optional


class Format(enum.Enum):
    """Recognized input, output, and sequence file formats."""

    TXT = 'TXT'
    BED = 'BED'
    TSV = 'TSV'
    CSV = 'CSV'
    PARQUET = 'PARQUET'
    FASTA = 'FASTA'
    FASTQ = 'FASTQ'


TEXT_TABLE_FORMATS = {Format.TXT, Format.BED, Format.TSV, Format.CSV}
"""Text formats that may be gzip-compressed."""

SEQ_FORMATS = {Format.FASTA, Format.FASTQ}
"""Sequence formats (may be bgzip-compressed)."""

_NAME_MAP = {fmt.value: fmt for fmt in Format}

_EXT_MAP = {
    'txt': Format.TXT,
    'bed': Format.BED,
    'tsv': Format.TSV,
    'csv': Format.CSV,
    'parquet': Format.PARQUET,
    'fasta': Format.FASTA,
    'fa': Format.FASTA,
    'fna': Format.FASTA,
    'fastq': Format.FASTQ,
    'fq': Format.FASTQ,
}


def parse_format(fmt_str: str) -> tuple[Format, bool]:
    """Parse a format option string (from -I, -O, or -Q).

    The value is case-insensitive and may end with "GZ" to indicate gzip/bgzip compression (e.g. "TXTGZ", "bedgz").

    :param fmt_str: Format string.

    :returns: A tuple of (format, is_gz).

    :raises ValueError: If the format is unrecognized or an unsupported gzip suffix is given.
    """
    if fmt_str is None or (value := fmt_str.strip().upper()) == '':
        raise ValueError('Empty format string')

    is_gz = False

    if value not in _NAME_MAP and value.endswith('GZ'):
        is_gz = True
        value = value[:-2]

    if value not in _NAME_MAP:
        raise ValueError(f'Unrecognized format: {fmt_str}')

    fmt = _NAME_MAP[value]

    if is_gz and fmt == Format.PARQUET:
        raise ValueError('Parquet cannot be gzip-compressed')

    return fmt, is_gz


def detect_format(filename: str | Path) -> Optional[tuple[Format, bool]]:
    """Detect a file format from a filename by its extension (case-insensitive).

    :param filename: File name.

    :returns: A tuple of (format, is_gz), or `None` if the format cannot be determined.

    :raises ValueError: If a gzip-compressed Parquet file is indicated (e.g. "x.parquet.gz"); Parquet is not
        gzip-compressed.
    """
    name = str(filename).strip().lower()
    is_gz = False

    if name.endswith('.gz'):
        is_gz = True
        name = name[:-3]

    if '.' not in name:
        return None

    ext = name.rsplit('.', 1)[1]
    fmt = _EXT_MAP.get(ext)

    if fmt is None:
        return None

    if is_gz and fmt == Format.PARQUET:
        raise ValueError(f'Parquet files cannot be gzip-compressed: {filename}')

    return fmt, is_gz
