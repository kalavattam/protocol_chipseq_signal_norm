#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_no_negative_signal.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude (Opus 5) was used in design, development, and documentation,
# with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


import numpy as np
import pytest

from protocol_chipseq_signal_norm.cli.compute_signal import (
    calc_sig_chrom_array,
    calc_sig_chrom_direct_sparse_np,
)

# Every quantity 'compute_signal' accumulates is non-negative: a fragment
# covers base pairs in a bin, weighted by '1/L' under 'norm' and 'frag', or by
# 1 under 'unadj' (where 'L' is fragment length). Nothing subtracts, so a
# negative output bin cannot be real data; it can only be arithmetic residue.
#
# Earlier implementations produced exactly that. Interiors are filled with a
# difference array ('+value' at the interior start, '-value' one past the end)
# reconstructed by a left-to-right 'np.cumsum'. In exact arithmetic the running
# total returns to zero between fragments; in 'float64' it does not, because
# many rounded additions separate a fragment's two marks. The residue is of
# order machine epsilon times the peak running sum, its sign is arbitrary, and
# the exact 'sig != 0.0' emit filter wrote it out: on real data, up to ~3% of
# bins near 1e-20, about half of them negative.
#
# The fix masks the 'float64' cumulative sum with the 'int64' interior-coverage
# count. Interior coverage is an integer question, so that indicator is exact
# and no tolerance has to be tuned.
RNG_SEED = 20260816


def _fragments(n_frag, chrom_size, min_len, max_len, seed=RNG_SEED):
    """
    Fragments scattered with gaps between them.

    Gaps are what matter here: each one is a place the running total has to
    return to zero, and therefore a place residue can appear. A densely tiled
    chromosome would hide the defect.
    """
    rng = np.random.default_rng(seed)
    starts = np.sort(rng.integers(0, chrom_size - max_len - 1, size=n_frag))
    lengths = rng.integers(min_len, max_len, size=n_frag)

    return (
        starts.astype(np.int64),
        (starts + lengths).astype(np.int64),
        lengths.astype(np.float64),
    )


# 'long fragments, sparse' is the case that reproduced the defect: long
# interiors make each fragment's '+value' and '-value' far apart in the
# cumulative sum, and wide gaps give the residue somewhere to land.
CASES = [
    pytest.param(20_000, 2_000_000, 40, 120, id="short-fragments-wide-gaps"),
    pytest.param(60_000, 1_000_000, 140, 200, id="nucleosome-scale-dense"),
    pytest.param(5_000, 3_000_000, 400, 900, id="long-fragments-sparse"),
    pytest.param(1, 100_000, 300, 301, id="single-fragment"),
]


@pytest.mark.parametrize(("n_frag", "size", "lo", "hi"), CASES)
@pytest.mark.parametrize("is_len", [False, True], ids=["unweighted", "length"])
@pytest.mark.parametrize("is_norm", [False, True], ids=["raw", "norm"])
def test_calc_sig_chrom_array_is_never_negative(
    n_frag, size, lo, hi, is_len, is_norm
):
    """Coverage from the array accumulator is non-negative in every bin."""
    starts, ends, lengths = _fragments(n_frag, size, lo, hi)

    out = calc_sig_chrom_array(
        chrom="I",
        starts=starts,
        ends=ends,
        lengths=lengths,
        chrom_size=size,
        siz_bin=10,
        is_len=is_len,
        is_norm=is_norm,
        fragment_count=float(n_frag),
        scl_fct=None,
    )

    values = np.fromiter(out.values(), dtype=np.float64, count=len(out))
    negative = values[values < 0]

    assert negative.size == 0, (
        f"{negative.size} bin(s) carry negative coverage, worst "
        f"{negative.min():.6g}; coverage is a sum of non-negative "
        f"contributions and cannot be below zero"
    )


@pytest.mark.parametrize(("n_frag", "size", "lo", "hi"), CASES)
@pytest.mark.parametrize("use_bincount", [False, True], ids=["B1", "B3"])
def test_direct_sparse_is_never_negative(n_frag, size, lo, hi, use_bincount):
    """
    Both sparse result formats are non-negative in every bin.

    'B1' and 'B3' build the difference array differently (in-place accumulation
    versus a subtraction of two 'np.bincount' sums), and the second is the more
    cancellation-prone of the two, so both are covered.
    """
    starts, ends, lengths = _fragments(n_frag, size, lo, hi)

    tag, parts = calc_sig_chrom_direct_sparse_np(
        chrom="I",
        starts=starts,
        ends=ends,
        lengths=lengths,
        chrom_size=size,
        siz_bin=10,
        is_len=True,
        use_bincount=use_bincount,
    )

    if not parts:
        pytest.skip("no covered bins for this configuration")

    values = np.concatenate(
        [np.asarray(part[2], dtype=np.float64) for part in parts]
    )
    negative = values[values < 0]

    assert negative.size == 0, (
        f"{negative.size} bin(s) carry negative coverage, worst "
        f"{negative.min():.6g} (use_bincount={use_bincount})"
    )


def test_zero_coverage_gaps_are_absent_not_tiny():
    """
    A gap between fragments emits no row, rather than a near-zero one.

    This is the other half of the same defect. Suppressing the sign is not
    enough: residue of '+1e-20' is exactly as spurious as '-1e-20', and a bin
    with no fragment over it should not appear in the output at all.
    """
    # Two well-separated fragments, so the bins between them are unambiguously
    # uncovered and any row there is residue.
    starts = np.array([1_000, 500_000], dtype=np.int64)
    ends = np.array([1_600, 500_600], dtype=np.int64)
    lengths = np.array([600.0, 600.0], dtype=np.float64)

    out = calc_sig_chrom_array(
        chrom="I",
        starts=starts,
        ends=ends,
        lengths=lengths,
        chrom_size=1_000_000,
        siz_bin=10,
        is_len=False,
        is_norm=True,
        fragment_count=2.0,
        scl_fct=None,
    )

    covered = set()
    for start, end in zip(starts.tolist(), ends.tolist()):
        covered.update(range(start // 10, ((end - 1) // 10) + 1))

    stray = {
        position // 10
        for (_, position) in out
        if position // 10 not in covered
    }

    assert not stray, (
        f"{len(stray)} bin(s) emitted outside any fragment; nearest strays "
        f"{sorted(stray)[:5]}"
    )
