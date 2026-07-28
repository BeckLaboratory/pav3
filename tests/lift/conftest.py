"""Fixtures for the lift/subseq tool tests.

Builds a small self-contained Align-mode fixture (reference FASTA + SAM) with one forward-aligned contig and one
reverse-aligned contig, so lifts and extracted sequences can be checked against known reference bases.
"""
from __future__ import annotations

import random
from dataclasses import dataclass
from pathlib import Path

import pysam
import pytest


@dataclass
class AlignFixture:
    """Paths and reference sequence for the Align-mode test fixture."""

    sam: Path
    ref_fa: Path
    ref: str


@pytest.fixture
def align_fixture(tmp_path: Path) -> AlignFixture:
    """Create a reference FASTA and a SAM with a forward and a reverse alignment.

    - ``chr1`` is a 100 bp reference.
    - ``ctg1`` (40 bp) aligns forward to chr1:21-60 (1-based), CIGAR 40=.
    - ``ctg2`` (30 bp) aligns reverse to chr1:71-100 (1-based), CIGAR 30=.
    """
    rng = random.Random(1234)
    ref = ''.join(rng.choice('ACGT') for _ in range(100))

    ref_fa = tmp_path / 'ref.fa'
    with open(ref_fa, 'w') as out_file:
        out_file.write('>chr1\n')
        for i in range(0, len(ref), 60):
            out_file.write(ref[i:i + 60] + '\n')
    pysam.faidx(str(ref_fa))

    seq_fwd = ref[20:60]        # ctg1 SEQ (reference-forward)
    seq_rev = ref[70:100]       # ctg2 SEQ (reference-forward; contig is the reverse complement)

    sam = tmp_path / 'aln.sam'
    with open(sam, 'w') as out_file:
        out_file.write('@HD\tVN:1.6\tSO:coordinate\n')
        out_file.write('@SQ\tSN:chr1\tLN:100\n')
        out_file.write(f'ctg1\t0\tchr1\t21\t60\t40=\t*\t0\t0\t{seq_fwd}\t*\n')
        out_file.write(f'ctg2\t16\tchr1\t71\t60\t30=\t*\t0\t0\t{seq_rev}\t*\n')

    return AlignFixture(sam=sam, ref_fa=ref_fa, ref=ref)
