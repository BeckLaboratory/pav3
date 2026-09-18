"""Tests for variant annotation utilities."""
from __future__ import annotations

import numpy as np
import pytest

from pav3.anno import perfect_homology, _true_len


#
# _true_len
#

@pytest.mark.parametrize(
    ('values', 'expected'),
    [
        ([], 0),
        ([False, True, True], 0),
        ([True, False, True], 1),
        ([True, True, False], 2),
        ([True, True, True], 3),
    ],
)
def test_true_len_counts_the_leading_run(values, expected):
    """The leading run of True values is counted, up to the end of the array."""
    assert _true_len(np.array(values, dtype=np.bool_)) == expected


#
# perfect_homology
#

def test_homology_runs_to_the_end_of_the_flank():
    """A flank matching the variant over its full length is fully homologous."""
    # Deleting one "AC" from a dinucleotide repeat: both flanks repeat the deleted bases.
    assert perfect_homology(seq_var='AC', seq_ref='AC' * 10, pos=8, end=8, is_rev=False) == (2, 2)


def test_homology_stops_at_a_mismatch():
    """Homology ends at the first base that does not match."""
    assert perfect_homology(seq_var='AAAA', seq_ref='CCCC' + 'AAATTTTT', pos=4, end=4, is_rev=False) == (0, 3)


def test_homology_is_zero_without_a_match():
    """A flank sharing no leading base with the variant has no homology."""
    assert perfect_homology(seq_var='AAAA', seq_ref='CCCC' + 'GGGG', pos=4, end=4, is_rev=False) == (0, 0)


def test_homology_at_a_sequence_boundary():
    """A variant at the start of the sequence has no upstream flank to match."""
    assert perfect_homology(seq_var='AC', seq_ref='AC' * 10, pos=0, end=0, is_rev=False) == (0, 2)


def test_homology_swaps_flanks_on_the_reverse_strand():
    """On the reverse strand the flanks are reverse-complemented and exchanged."""
    fwd = perfect_homology(seq_var='AAAA', seq_ref='CCCC' + 'AAATTTTT', pos=4, end=4, is_rev=False)
    rev = perfect_homology(seq_var='TTTT', seq_ref='CCCC' + 'AAATTTTT', pos=4, end=4, is_rev=True)

    assert fwd == (0, 3)
    assert rev == (3, 0)
