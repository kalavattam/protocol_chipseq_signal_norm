#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: algotest_variants_test.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.6);
# - Anthropic Claude Code (Opus 5, Fable 5).
#
# Distributed under the MIT license.


"""
Exercise the quarantined 'compute_signal' algorithm-testing variants.

These tests moved out of 'tests/unit/compute_signal/' together with the
experimental machinery they prove. They are not part of the maintained safe
suite; run them directly with 'pytest dev/algorithm_testing/tests'. The
maintained suite must stay green without them.
"""

import numpy as np
import pytest
from dev.algorithm_testing.compute_signal_algotest import (
    calc_sig_chrom,
    calc_sig_chrom_array,
    calc_sig_chrom_direct_dense_np,
    calc_sig_chrom_direct_sparse_np,
    calc_sig_chrom_event_np,
    calc_sig_task,
    materialize_event_signal_parts,
)

from protocol_chipseq_signal_norm.cli.compute_signal import (
    calc_sig_chrom_direct_sparse_np as calc_sig_chrom_direct_sparse_b1,
)

RNG_SEED = 20260816

# Shared fragment geometries reused from the production non-negativity
# suite: gaps are where cumsum residue can land, and long interiors are what
# reproduced the pre-fix defect.
CASES = [
    pytest.param(20_000, 2_000_000, 40, 120, id="short-fragments-wide-gaps"),
    pytest.param(60_000, 1_000_000, 140, 200, id="nucleosome-scale-dense"),
    pytest.param(5_000, 3_000_000, 400, 900, id="long-fragments-sparse"),
    pytest.param(1, 100_000, 300, 301, id="single-fragment"),
]

# Variant configurations of the full-featured sparse accumulator, keyed by
# the retired '--prototype_result_format' spellings they realized.
SPARSE_VARIANTS = [
    pytest.param({"local_span": True}, id="direct_sparse_local_np"),
    pytest.param({"use_bincount": True}, id="direct_sparse_bincount_np"),
    pytest.param(
        {"local_span": True, "use_bincount": True},
        id="direct_sparse_local_bincount_np",
    ),
    pytest.param({"emit_touched_only": True}, id="direct_sparse_touched_np"),
    pytest.param(
        {"emit_touched_only": True, "use_bincount": True},
        id="direct_sparse_touched_bincount_np",
    ),
]


def _fragments(
    n_frag: int,
    chrom_size: int,
    min_len: int,
    max_len: int,
    seed: int = RNG_SEED,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Fragments scattered with gaps between them.
    """

    rng = np.random.default_rng(seed)
    starts = np.sort(rng.integers(0, chrom_size - max_len - 1, size=n_frag))
    lengths = rng.integers(min_len, max_len, size=n_frag)

    return (
        starts.astype(np.int64),
        (starts + lengths).astype(np.int64),
        lengths.astype(np.float64),
    )


def _baseline_sparse(
    starts: np.ndarray,
    ends: np.ndarray,
    lengths: np.ndarray,
    chrom_size: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute the production B1 result as concatenated starts and values.
    """

    _tag, parts = calc_sig_chrom_direct_sparse_b1(
        chrom="I",
        starts=starts,
        ends=ends,
        lengths=lengths,
        chrom_size=chrom_size,
        siz_bin=10,
        is_len=True,
    )

    if not parts:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.float64)

    return parts[0][1], parts[0][2]


def test_calc_sig_chrom_accumulates_fragment_overlap_by_bin() -> None:
    observed = calc_sig_chrom(
        "chrI",
        [(0, 10, 10), (5, 15, 10)],
        fragment_count=2,
        siz_bin=10,
        is_len=False,
        is_norm=False,
    )

    assert observed[("chrI", 0)] == 15.0
    assert observed[("chrI", 10)] == 5.0


def test_calc_sig_chrom_length_and_count_normalization() -> None:
    observed = calc_sig_chrom(
        "chrI",
        [(0, 10, 10), (5, 15, 10)],
        fragment_count=2,
        siz_bin=10,
        is_len=True,
        is_norm=True,
        scl_fct=4.0,
    )

    assert observed[("chrI", 0)] == pytest.approx(3.0)
    assert observed[("chrI", 10)] == pytest.approx(1.0)


def test_calc_sig_task_dispatches_to_calc_sig_chrom() -> None:
    task = ("chrI", [(0, 5, 5)], 1, 10, False, False, None)

    assert calc_sig_task(task)[("chrI", 0)] == 5.0


def test_calc_sig_chrom_rejects_invalid_bin_size() -> None:
    with pytest.raises(ValueError, match="siz_bin"):
        calc_sig_chrom("chrI", [], 0, 0, False, False)


def test_calc_sig_chrom_array_matches_python_kernel() -> None:
    fragments = [(0, 10, 10), (5, 25, 20), (12, 15, 3), (79, 80, 1)]
    expected = calc_sig_chrom(
        "chrI",
        fragments,
        fragment_count=len(fragments),
        siz_bin=10,
        is_len=True,
        is_norm=True,
        scl_fct=4.0,
    )

    observed = calc_sig_chrom_array(
        "chrI",
        np.asarray([frag[0] for frag in fragments]),
        np.asarray([frag[1] for frag in fragments]),
        np.asarray([frag[2] for frag in fragments], dtype=float),
        chrom_size=80,
        fragment_count=len(fragments),
        siz_bin=10,
        is_len=True,
        is_norm=True,
        scl_fct=4.0,
    )

    assert observed.keys() == expected.keys()

    for key in expected:
        assert observed[key] == pytest.approx(expected[key])


def test_calc_sig_chrom_array_emits_only_touched_bins() -> None:
    fragments = [
        (0, 26040, 3519),
        (26080, 40000, 3519),
    ]
    expected = calc_sig_chrom(
        "chrI",
        fragments,
        fragment_count=len(fragments),
        siz_bin=10,
        is_len=True,
        is_norm=True,
    )

    observed = calc_sig_chrom_array(
        "chrI",
        np.asarray([frag[0] for frag in fragments]),
        np.asarray([frag[1] for frag in fragments]),
        np.asarray([frag[2] for frag in fragments], dtype=float),
        chrom_size=50000,
        fragment_count=len(fragments),
        siz_bin=10,
        is_len=True,
        is_norm=True,
    )

    assert observed.keys() == expected.keys()
    assert ("chrI", 26040) not in observed
    assert ("chrI", 26050) not in observed
    assert ("chrI", 26060) not in observed
    assert ("chrI", 26070) not in observed

    for key in expected:
        assert observed[key] == pytest.approx(expected[key])


@pytest.mark.parametrize("kwargs", SPARSE_VARIANTS)
def test_sparse_variants_match_the_public_b1_baseline(
    kwargs: dict[str, bool],
) -> None:
    # The retired hidden CLI formats asserted rendered-output equality with
    # the public baseline; relocated here, the same claim is made at the
    # function level. Emitted coordinates must agree exactly (post-fix, the
    # emission gate is the same integer indicator in every branch), and
    # values must agree to float64 rounding, the documented order-of-
    # summation difference between the scatter and bincount branches.
    starts, ends, lengths = _fragments(5_000, 3_000_000, 400, 900)
    base_starts, base_values = _baseline_sparse(
        starts,
        ends,
        lengths,
        3_000_000,
    )

    _tag, parts = calc_sig_chrom_direct_sparse_np(
        chrom="I",
        starts=starts,
        ends=ends,
        lengths=lengths,
        chrom_size=3_000_000,
        siz_bin=10,
        is_len=True,
        **kwargs,
    )

    assert parts

    np.testing.assert_array_equal(parts[0][1], base_starts)
    np.testing.assert_allclose(parts[0][2], base_values, rtol=1e-12)


def test_idx_variant_carries_bin_indices_for_the_same_bins() -> None:
    # The 'direct_sparse_idx_np' format is the same arithmetic as B1 with
    # coordinates reported as bin indices instead of base-pair starts.
    starts, ends, lengths = _fragments(5_000, 3_000_000, 400, 900)
    base_starts, base_values = _baseline_sparse(
        starts,
        ends,
        lengths,
        3_000_000,
    )

    tag, parts = calc_sig_chrom_direct_sparse_np(
        chrom="I",
        starts=starts,
        ends=ends,
        lengths=lengths,
        chrom_size=3_000_000,
        siz_bin=10,
        is_len=True,
        return_bin_indices=True,
    )

    assert tag == "direct_sparse_idx_np"
    assert parts

    np.testing.assert_array_equal(
        parts[0][1].astype(np.int64) * 10,
        base_starts,
    )
    np.testing.assert_allclose(parts[0][2], base_values, rtol=1e-12)


def test_dense_variant_matches_the_public_b1_baseline() -> None:
    starts, ends, lengths = _fragments(5_000, 3_000_000, 400, 900)
    base_starts, base_values = _baseline_sparse(
        starts,
        ends,
        lengths,
        3_000_000,
    )

    _tag, parts = calc_sig_chrom_direct_dense_np(
        chrom="I",
        starts=starts,
        ends=ends,
        lengths=lengths,
        chrom_size=3_000_000,
        siz_bin=10,
        is_len=True,
    )

    assert parts

    _chrom, values, touched = parts[0]
    idx = np.flatnonzero(touched & (values != 0.0))

    np.testing.assert_array_equal(idx.astype(np.int64) * 10, base_starts)
    np.testing.assert_allclose(values[idx], base_values, rtol=1e-12)


def test_event_variant_materializes_to_the_public_b1_baseline() -> None:
    starts, ends, lengths = _fragments(5_000, 3_000_000, 400, 900)
    base_starts, base_values = _baseline_sparse(
        starts,
        ends,
        lengths,
        3_000_000,
    )

    _tag, parts = calc_sig_chrom_event_np(
        chrom="I",
        starts=starts,
        ends=ends,
        lengths=lengths,
        chrom_size=3_000_000,
        siz_bin=10,
        is_len=True,
    )

    assert parts

    _tag_merged, merged = materialize_event_signal_parts(
        {"I": [tuple(parts[0][1:])]},
        10,
    )

    assert merged

    np.testing.assert_array_equal(merged[0][1], base_starts)
    np.testing.assert_allclose(merged[0][2], base_values, rtol=1e-12)


@pytest.mark.parametrize(("n_frag", "size", "lo", "hi"), CASES)
@pytest.mark.parametrize("is_len", [False, True], ids=["unweighted", "length"])
@pytest.mark.parametrize("is_norm", [False, True], ids=["raw", "norm"])
def test_calc_sig_chrom_array_is_never_negative(
    n_frag: int,
    size: int,
    lo: int,
    hi: int,
    is_len: bool,
    is_norm: bool,
) -> None:
    """
    Coverage from the array accumulator is non-negative in every bin.
    """

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
def test_b3_direct_sparse_is_never_negative(
    n_frag: int,
    size: int,
    lo: int,
    hi: int,
) -> None:
    """
    The 'B3' bincount branch is non-negative in every bin.

    'B1' and 'B3' build the difference array differently (in-place accumulation
    versus a subtraction of two 'np.bincount' sums), and the second is the more
    cancellation-prone of the two. The 'B1' arm of this assertion lives with
    the production accumulator in
    'tests/unit/compute_signal/test_no_negative_signal.py'.
    """

    starts, ends, lengths = _fragments(n_frag, size, lo, hi)

    _tag, parts = calc_sig_chrom_direct_sparse_np(
        chrom="I",
        starts=starts,
        ends=ends,
        lengths=lengths,
        chrom_size=size,
        siz_bin=10,
        is_len=True,
        use_bincount=True,
    )

    if not parts:
        pytest.skip("no covered bins for this configuration")

    values = np.concatenate(
        [np.asarray(part[2], dtype=np.float64) for part in parts],
    )
    negative = values[values < 0]

    assert negative.size == 0, (
        f"{negative.size} bin(s) carry negative coverage, worst "
        f"{negative.min():.6g} (use_bincount=True)"
    )
