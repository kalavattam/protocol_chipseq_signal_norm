#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: compute_signal_algotest.py
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
#   GPT-5.6);
# - Anthropic Claude Code (Opus 5, Fable 5).
#
# Distributed under the MIT license.


"""
Quarantined algorithm-testing variants of the 'compute_signal' accumulators.

This module preserves, unchanged, the experimental signal machinery removed
from 'protocol_chipseq_signal_norm.cli.compute_signal' when the production file
was reduced to its single public arithmetic ('direct_sparse_np', the
private-roadmap shorthand 'B1'). It is the full comparison harness: the
alternative result formats (including 'B3', 'direct_sparse_bincount_np'), the
dense and event accumulators, the non-default merge strategies, and the
digest/profile-only write helpers all live here, next to the reference dict
kernel they were originally benchmarked against.

Nothing here is production code. The module exists so that the July 2026
benchmark and correctness comparisons stay runnable from source, and it is
slated for adoption by a future archive repository. Production behavior must
never import from this module; the dependency direction is one-way, with this
module importing shared helpers from the installed package.

The S3-MIG-009 cumsum-residue masking fix (commit '20c337f') is retained at
every accumulation site preserved here.
"""

from __future__ import annotations

import hashlib
import math
import os
import time
from collections import defaultdict
from collections.abc import Callable, Iterable, Iterator
from concurrent.futures import ProcessPoolExecutor
from typing import Any

import numpy as np

from protocol_chipseq_signal_norm.cli.compute_signal import (
    coalesce_parts,
    collect_frg_arr,
    format_bdg_rows,
    format_bdg_rows_task,
    format_bdg_value,
    is_sparse_sig,
    iter_aln_frg,
    iter_idx_frg,
)
from protocol_chipseq_signal_norm.utilities.utils_bdg import (
    key_bin,
    write_bdg,
)
from protocol_chipseq_signal_norm.utilities.utils_check import (
    validate_comparison,
)
from protocol_chipseq_signal_norm.utilities.utils_chrom import sort_chrom
from protocol_chipseq_signal_norm.utilities.utils_io import open_out

# Every prototype result format the benchmarking harness understood. Only
# 'direct_sparse_np' (B1) survives in production; the rest are exercised
# through this module.
PROTOTYPE_RESULT_FORMAT_CHOICES = (
    "dict",
    "sparse_np",
    "dense_np",
    "direct_sparse_np",
    "direct_sparse_idx_np",
    "direct_sparse_local_np",
    "direct_sparse_touched_np",
    "direct_sparse_bincount_np",
    "direct_sparse_touched_bincount_np",
    "direct_sparse_local_bincount_np",
    "direct_dense_np",
    "event_np",
)
PROTOTYPE_MERGE_STRATEGY_CHOICES = (
    "dict_merge",
    "chrom_array_merge",
    "vectorized_merge",
    "array_sparse_merge_legacy",
    "array_sparse_merge",
    "array_dense_merge",
    "event_diff_merge",
)
PROTOTYPE_WRITE_MODE_CHOICES = ("full", "digest", "profile_only")
DIRECT_SPARSE_RESULT_FORMATS = (
    "direct_sparse_np",
    "direct_sparse_idx_np",
    "direct_sparse_local_np",
    "direct_sparse_touched_np",
    "direct_sparse_bincount_np",
    "direct_sparse_touched_bincount_np",
    "direct_sparse_local_bincount_np",
)


def signal_dict_to_result(
    sig: dict[tuple[str, int], float],
    result_format: str,
    siz_bin: int,
) -> object:
    """
    Convert a signal dict to one private benchmark result representation.

    Parameters
    ----------
    sig : dict[tuple[str, int], float]
        Sparse signal keyed by chromosome and bin start.
    result_format : str
        Requested private result representation.
    siz_bin : int
        Bin width used to derive array coordinates.

    Returns
    -------
    result : object
        Signal in the requested private benchmark representation.

    Raises
    ------
    ValueError
        If 'result_format' is not recognized.
    """

    if result_format == "dict":
        return sig

    by_chrom: dict[str, list[tuple[int, float]]] = defaultdict(list)

    for (chrom, start), value in sig.items():
        by_chrom[chrom].append((start, value))

    if result_format in ("sparse_np", "direct_sparse_np"):
        return (
            result_format,
            [
                (
                    chrom,
                    np.asarray([start for start, _ in rows], dtype=np.int64),
                    np.asarray([value for _, value in rows], dtype=np.float64),
                )
                for chrom, rows in by_chrom.items()
            ],
        )

    if result_format == "dense_np":
        dense_parts = []

        for chrom, rows in by_chrom.items():
            starts = [start for start, _ in rows]
            first = min(starts)
            last = max(starts)
            values = np.zeros(
                ((last - first) // siz_bin) + 1,
                dtype=np.float64,
            )

            for start, value in rows:
                values[(start - first) // siz_bin] += value

            dense_parts.append((chrom, first, values))

        return "dense_np", dense_parts

    raise ValueError(f"Unknown prototype result format: {result_format!r}.")


def count_result_bins(result: Any) -> int:
    """
    Count signal bins for standard dict results and prototype tuple results.

    Parameters
    ----------
    result : Any
        Supported mapping- or tuple-backed signal result.

    Returns
    -------
    count : int
        Number of represented signal bins.
    """

    if isinstance(result, tuple) and result:
        if len(result) == 2 and isinstance(result[1], int):
            return count_result_bins(result[0])

        tag = result[0]
        if tag in ("sparse_np", "direct_sparse_np", "direct_sparse_idx_np"):
            return sum(int(starts.size) for _, starts, _ in result[1])

        if tag == "sparse_merged":
            return sum(int(starts.size) for _, starts, _ in result[1])

        if tag == "direct_dense_np":
            return sum(
                int(np.count_nonzero(touched & (values != 0.0)))
                for _, values, touched in result[1]
            )

        if tag == "array_dense_merged":
            return sum(
                int(np.count_nonzero(touched & (values != 0.0)))
                for _, values, touched in result[2]
            )

        if tag == "event_np":
            return sum(
                int(edge_bins.size + diff_bins.size + touch_bins.size)
                for (
                    _,
                    _,
                    edge_bins,
                    _,
                    diff_bins,
                    _,
                    touch_bins,
                    _,
                ) in result[1]
            )

        if tag == "dense_np":
            return sum(
                int(np.count_nonzero(values)) for _, _, values in result[1]
            )

        if isinstance(tag, dict):
            return len(tag)

    return len(result)


def est_result_bytes(result: Any) -> int:
    """
    Estimate NumPy-array payload bytes for private benchmark result formats.

    Parameters
    ----------
    result : Any
        Supported mapping- or array-backed signal result.

    Returns
    -------
    byte_count : int
        Sum of owned NumPy array payload bytes.
    """

    if isinstance(result, tuple) and result:
        tag = result[0]
        if tag in ("sparse_np", "direct_sparse_np", "direct_sparse_idx_np"):
            return sum(
                int(starts.nbytes + values.nbytes)
                for _, starts, values in result[1]
            )

        if tag == "sparse_merged":
            return sum(
                int(starts.nbytes + values.nbytes)
                for _, starts, values in result[1]
            )

        if tag == "direct_dense_np":
            return sum(
                int(values.nbytes + touched.nbytes)
                for _, values, touched in result[1]
            )

        if tag == "array_dense_merged":
            return sum(
                int(values.nbytes + touched.nbytes)
                for _, values, touched in result[2]
            )

        if tag == "event_np":
            return sum(
                int(
                    edge_bins.nbytes
                    + edge_values.nbytes
                    + diff_bins.nbytes
                    + diff_values.nbytes
                    + touch_bins.nbytes
                    + touch_values.nbytes,
                )
                for (
                    _,
                    _,
                    edge_bins,
                    edge_values,
                    diff_bins,
                    diff_values,
                    touch_bins,
                    touch_values,
                ) in result[1]
            )

        if tag == "dense_np":
            return sum(int(values.nbytes) for _, _, values in result[1])

        if isinstance(tag, dict):
            return 0

    return 0


def collect_fragment_arrays(
    fil_aln: str,
    siz_chr: dict[str, int],
    usr_frg: int | None = None,
    ref_fa: str | None = None,
    allow_sec: bool = True,
    allow_supp: bool = False,
    allow_dup: bool = True,
) -> tuple[dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]], int]:
    """
    Parse a BAM or CRAM once into per-chromosome NumPy fragment arrays.
    """

    return collect_frg_arr(
        iter_aln_frg(
            fil_aln=fil_aln,
            siz_chr=siz_chr,
            usr_frg=usr_frg,
            ref_fa=ref_fa,
            allow_sec=allow_sec,
            allow_supp=allow_supp,
            allow_dup=allow_dup,
        ),
    )


def calc_sig_chrom(
    chrom: str,
    fragment_records: list[tuple[int, int, int]],
    fragment_count: int,
    siz_bin: int,
    is_len: bool,
    is_norm: bool,
    scl_fct: float | None = None,
) -> dict[tuple[str, int], float]:
    """
    Compute one chromosome of exact fragment-bin overlap signal.

    Function respects half-open '[start, end)' fragments and avoids per-base
    loops. Overlap is measured in bases, and output is per-bin sums of per-base
    contributions.

    If provided, 'scl_fct' is applied after any optional normalization.
    ('scl_fct > 0' is required to avoid silent zeroing.)

    Parameters
    ----------
    chrom : str
        Chromosome name.
    fragment_records : list[tuple[int, int, int]]
        Start, end, and fragment-length tuples.
    fragment_count : int
        Total number of fragments (used when 'is_norm=True').
    siz_bin : int
        Bin size in base pairs.
    is_len : bool
        If 'True', normalize by fragment length.
    is_norm : bool
        If 'True', normalize by both fragment length and total fragment count
        so that the genome-wide summed signal is approximately 1.
    scl_fct : float | None
        Scaling factor applied to signal.

    Returns
    -------
    sig_bin : dict[tuple[str, int], float]
        A dictionary of binned signal data, where keys are
        '(chrom, bin_start)' and values are per-bin signal scores.

    Raises
    ------
    ValueError
        - If 'siz_bin' <= 0.
        - If 'is_len' or 'is_norm' is True and any fragment length is
          nonpositive.
        - If 'is_norm' is True and 'fragment_count' <= 0.
        - If a provided 'scl_fct' <= 0.

    Notes
    -----
    - If 'is_len=True' or 'is_norm=True', each fragment contributes one divided
      by the fragment length per covered base.
    - If 'is_norm=True', signal is additionally divided by 'fragment_count' so
      that genome-wide summed signal is approximately 1.
    - If 'scl_fct' is provided, signal values are scaled accordingly.
    - Signal is accumulated per bin according to the number of overlapping
      bases contributed by each fragment.
    """

    validate_comparison(siz_bin, "gt", 0, "siz_bin", allow_none=False)

    sig_bin = defaultdict(float)

    for fragment_start, fragment_end, fragment_length in fragment_records:
        if fragment_end <= fragment_start:
            continue

        if (is_len or is_norm) and fragment_length <= 0:
            raise ValueError("'frg_len' must be > 0 when using normalization.")

        per_bas = (1.0 / fragment_length) if (is_len or is_norm) else 1.0

        first_bin_start = (fragment_start // siz_bin) * siz_bin
        last_bin_start = ((fragment_end - 1) // siz_bin) * siz_bin

        for bin_start in range(first_bin_start, last_bin_start + 1, siz_bin):
            bin_end = bin_start + siz_bin
            overlap_start = max(fragment_start, bin_start)
            overlap_end = min(fragment_end, bin_end)
            overlap = overlap_end - overlap_start

            if overlap > 0:
                sig_bin[(chrom, bin_start)] += per_bas * overlap

    # Normalized signal sums to approximately 1 across the genome.
    if is_norm:
        if fragment_count <= 0:
            raise ValueError(
                "Normalization requires non-zero total fragments.",
            )

        for k in sig_bin:
            sig_bin[k] /= fragment_count

    if scl_fct is not None:
        validate_comparison(scl_fct, "gt", 0, "scl_fct")

        for k in sig_bin:
            sig_bin[k] *= scl_fct

    return sig_bin


def calc_sig_chrom_array(
    chrom: str,
    starts: np.ndarray,
    ends: np.ndarray,
    lengths: np.ndarray,
    chrom_size: int,
    fragment_count: int,
    siz_bin: int,
    is_len: bool,
    is_norm: bool,
    scl_fct: float | None = None,
) -> dict[tuple[str, int], float]:
    """
    Compute binned signal for one chromosome with vectorized range additions.

    Parameters
    ----------
    chrom : str
        Chromosome represented by the fragment arrays.
    starts, ends, lengths : np.ndarray
        Fragment starts, ends, and lengths in matching order.
    chrom_size : int
        Chromosome length used to clip bins.
    fragment_count : int
        Total accepted fragments used for optional normalization.
    siz_bin : int
        Signal-bin width.
    is_len : bool
        Whether fragments contribute inverse-length weights.
    is_norm : bool
        Whether to apply depth normalization.
    scl_fct : float | None
        Optional multiplicative scale factor.

    Returns
    -------
    signal : dict[tuple[str, int], float]
        Sparse nonzero signal keyed by chromosome and bin start.

    Raises
    ------
    ValueError
        If chromosome size, fragment lengths, or normalization inputs are
        invalid.
    """

    validate_comparison(siz_bin, "gt", 0, "siz_bin", allow_none=False)

    if chrom_size <= 0:
        raise ValueError(f"Chromosome size must be > 0 for {chrom!r}.")

    n_bins = math.ceil(chrom_size / siz_bin)
    if starts.size == 0:
        return {}

    valid = ends > starts

    if not np.all(valid):
        starts = starts[valid]
        ends = ends[valid]
        lengths = lengths[valid]

    if starts.size == 0:
        return {}

    if is_len or is_norm:
        if np.any(lengths <= 0):
            raise ValueError("'frg_len' must be > 0 when using normalization.")

        weights = 1.0 / lengths
    else:
        weights = np.ones(starts.shape, dtype=np.float64)

    sig = np.zeros(n_bins, dtype=np.float64)
    touched = np.zeros(n_bins, dtype=bool)
    start_bins = starts // siz_bin
    end_bins = (ends - 1) // siz_bin

    same = start_bins == end_bins

    if np.any(same):
        np.add.at(
            sig,
            start_bins[same],
            (ends[same] - starts[same]) * weights[same],
        )
        touched[start_bins[same]] = True

    multi = ~same

    if np.any(multi):
        starts_subset = starts[multi]
        ends_subset = ends[multi]
        start_bins_subset = start_bins[multi]
        end_bins_subset = end_bins[multi]
        weights_subset = weights[multi]

        left_overlap = (
            (start_bins_subset + 1) * siz_bin - starts_subset
        ) * weights_subset
        right_overlap = (
            ends_subset - end_bins_subset * siz_bin
        ) * weights_subset
        np.add.at(sig, start_bins_subset, left_overlap)
        np.add.at(sig, end_bins_subset, right_overlap)
        touched[start_bins_subset] = True
        touched[end_bins_subset] = True

        has_interior = end_bins_subset > start_bins_subset + 1

        if np.any(has_interior):
            diff = np.zeros(n_bins + 1, dtype=np.float64)
            touched_diff = np.zeros(n_bins + 1, dtype=np.int64)
            first = start_bins_subset[has_interior] + 1
            last_exclusive = end_bins_subset[has_interior]
            value = siz_bin * weights_subset[has_interior]

            np.add.at(diff, first, value)
            np.add.at(diff, last_exclusive, -value)

            np.add.at(touched_diff, first, 1)
            np.add.at(touched_diff, last_exclusive, -1)

            # Interior coverage is exact in 'int64'; the 'float64' value sum is
            # not. Floating-point cancellation can leave ~1e-20 residue in bins
            # whose true value is zero, which 'sig != 0.0' would emit as data.
            # Masking by the exact indicator restores zero.
            interior = np.cumsum(touched_diff[:-1]) > 0

            sig += np.where(interior, np.cumsum(diff[:-1]), 0.0)
            touched |= interior

    if is_norm:
        if fragment_count <= 0:
            raise ValueError(
                "Normalization requires non-zero total fragments.",
            )

        sig /= fragment_count

    if scl_fct is not None:
        validate_comparison(scl_fct, "gt", 0, "scl_fct")

        sig *= scl_fct

    nonzero_indices = np.flatnonzero(touched & (sig != 0.0))

    return {
        (chrom, int(index) * siz_bin): float(sig[index])
        for index in nonzero_indices
    }


def calc_sig_chrom_direct_sparse_np(
    chrom: str,
    starts: np.ndarray,
    ends: np.ndarray,
    lengths: np.ndarray,
    chrom_size: int,
    siz_bin: int,
    is_len: bool,
    return_bin_indices: bool = False,
    local_span: bool = False,
    use_bincount: bool = False,
    emit_touched_only: bool = False,
) -> object:
    """
    Compute binned signal and return sparse NumPy arrays without a dict.

    Parameters
    ----------
    chrom : str
        Chromosome identifier for the selected signal calculation.
    starts : np.ndarray
        Zero-based interval start coordinates.
    ends : np.ndarray
        Zero-based interval end coordinates.
    lengths : np.ndarray
        Fragment lengths aligned with the interval coordinates.
    chrom_size : int
        Chromosome length in base pairs.
    siz_bin : int
        Positive signal-bin width in base pairs.
    is_len : bool
        Whether signal weights are normalized by fragment length.
    return_bin_indices : bool
        Whether sparse results use bin indices instead of base coordinates.
    local_span : bool
        Whether allocation is bounded to the observed local span.
    use_bincount : bool
        Whether aggregation uses NumPy bincount.
    emit_touched_only : bool
        Whether sparse output omits untouched or zero-valued bins.

    Returns
    -------
    result : object
        Strategy tag and sparse chromosome result payload.

    Raises
    ------
    ValueError
        If a supplied value violates the validated contract.
    """

    validate_comparison(siz_bin, "gt", 0, "siz_bin", allow_none=False)

    if chrom_size <= 0:
        raise ValueError(f"Chromosome size must be > 0 for {chrom!r}.")

    genome_n_bins = math.ceil(chrom_size / siz_bin)
    if starts.size == 0:
        return (
            "direct_sparse_idx_np"
            if return_bin_indices
            else "direct_sparse_np",
            [],
        )

    valid = ends > starts

    if not np.all(valid):
        starts = starts[valid]
        ends = ends[valid]
        lengths = lengths[valid]

    if starts.size == 0:
        return (
            "direct_sparse_idx_np"
            if return_bin_indices
            else "direct_sparse_np",
            [],
        )

    if is_len:
        if np.any(lengths <= 0):
            raise ValueError("'frg_len' must be > 0 when using normalization.")

        weights = 1.0 / lengths
    else:
        weights = np.ones(starts.shape, dtype=np.float64)

    start_bins_global = starts // siz_bin
    end_bins_global = (ends - 1) // siz_bin

    if local_span:
        first_bin = int(np.min(start_bins_global))
        last_bin = int(np.max(end_bins_global))
    else:
        first_bin = 0
        last_bin = genome_n_bins - 1

    first_bin = max(0, first_bin)
    last_bin = min(genome_n_bins - 1, last_bin)
    n_bins = last_bin - first_bin + 1
    if n_bins <= 0:
        return (
            "direct_sparse_idx_np"
            if return_bin_indices
            else "direct_sparse_np",
            [],
        )

    sig = np.zeros(n_bins, dtype=np.float64)
    touched = np.zeros(n_bins, dtype=bool) if emit_touched_only else None
    start_bins = start_bins_global - first_bin
    end_bins = end_bins_global - first_bin

    same = start_bins == end_bins

    if np.any(same):
        same_values = (ends[same] - starts[same]) * weights[same]

        if use_bincount:
            sig += np.bincount(
                start_bins[same],
                weights=same_values,
                minlength=n_bins,
            )[:n_bins]
        else:
            np.add.at(sig, start_bins[same], same_values)

        if touched is not None:
            touched[start_bins[same]] = True

    multi = ~same

    if np.any(multi):
        starts_subset = starts[multi]
        ends_subset = ends[multi]
        start_bins_subset = start_bins[multi]
        end_bins_subset = end_bins[multi]
        global_start_bins_subset = start_bins_global[multi]
        global_end_bins_subset = end_bins_global[multi]
        weights_subset = weights[multi]

        left_overlap = (
            (global_start_bins_subset + 1) * siz_bin - starts_subset
        ) * weights_subset
        right_overlap = (
            ends_subset - global_end_bins_subset * siz_bin
        ) * weights_subset

        if use_bincount:
            sig += np.bincount(
                start_bins_subset,
                weights=left_overlap,
                minlength=n_bins,
            )[:n_bins]
            sig += np.bincount(
                end_bins_subset,
                weights=right_overlap,
                minlength=n_bins,
            )[:n_bins]
        else:
            np.add.at(sig, start_bins_subset, left_overlap)
            np.add.at(sig, end_bins_subset, right_overlap)

        if touched is not None:
            touched[start_bins_subset] = True
            touched[end_bins_subset] = True

        has_interior = end_bins_subset > start_bins_subset + 1

        if np.any(has_interior):
            first = start_bins_subset[has_interior] + 1
            last_exclusive = end_bins_subset[has_interior]
            value = siz_bin * weights_subset[has_interior]

            if use_bincount:
                diff = np.bincount(
                    first,
                    weights=value,
                    minlength=n_bins + 1,
                ).astype(np.float64, copy=False)
                diff -= np.bincount(
                    last_exclusive,
                    weights=value,
                    minlength=n_bins + 1,
                ).astype(np.float64, copy=False)
            else:
                diff = np.zeros(n_bins + 1, dtype=np.float64)
                np.add.at(diff, first, value)
                np.add.at(diff, last_exclusive, -value)

            # The 'int64' interior sum is exact; 'float64' rounding can leave
            # ~1e-20 residue in true-zero bins that 'sig != 0.0' would emit.
            # Masking by the integer sum restores exact zero, tolerance-free.
            interior_diff = np.zeros(n_bins + 1, dtype=np.int64)
            np.add.at(interior_diff, first, 1)
            np.add.at(interior_diff, last_exclusive, -1)
            interior = np.cumsum(interior_diff[:-1]) > 0

            sig += np.where(interior, np.cumsum(diff[:-1]), 0.0)

            if touched is not None:
                touched |= interior

    if touched is None:
        idx = np.flatnonzero(sig != 0.0)
    else:
        idx = np.flatnonzero(touched & (sig != 0.0))

    if idx.size == 0:
        return (
            "direct_sparse_idx_np"
            if return_bin_indices
            else "direct_sparse_np",
            [],
        )

    if return_bin_indices:
        index_dtype = (
            np.int32 if n_bins <= np.iinfo(np.int32).max else np.int64
        )

        return (
            "direct_sparse_idx_np",
            [
                (
                    chrom,
                    (idx + first_bin).astype(index_dtype, copy=False),
                    sig[idx].astype(np.float64, copy=False),
                ),
            ],
        )

    return (
        "direct_sparse_np",
        [
            (
                chrom,
                ((idx + first_bin) * siz_bin).astype(np.int64, copy=False),
                sig[idx].astype(np.float64, copy=False),
            ),
        ],
    )


def calc_sig_chrom_direct_dense_np(
    chrom: str,
    starts: np.ndarray,
    ends: np.ndarray,
    lengths: np.ndarray,
    chrom_size: int,
    siz_bin: int,
    is_len: bool,
) -> object:
    """
    Compute binned signal and return one dense chromosome array.

    Parameters
    ----------
    chrom : str
        Chromosome represented by the fragment arrays.
    starts, ends, lengths : np.ndarray
        Fragment starts, ends, and lengths in matching order.
    chrom_size : int
        Chromosome length used to clip bins.
    siz_bin : int
        Signal-bin width.
    is_len : bool
        Whether fragments contribute inverse-length weights.

    Returns
    -------
    result : object
        Tagged chromosome result containing dense signal values.

    Raises
    ------
    ValueError
        If chromosome size or fragment lengths are invalid.
    """

    validate_comparison(siz_bin, "gt", 0, "siz_bin", allow_none=False)

    if chrom_size <= 0:
        raise ValueError(f"Chromosome size must be > 0 for {chrom!r}.")

    n_bins = math.ceil(chrom_size / siz_bin)
    if starts.size == 0:
        return "direct_dense_np", []

    valid = ends > starts

    if not np.all(valid):
        starts = starts[valid]
        ends = ends[valid]
        lengths = lengths[valid]

    if starts.size == 0:
        return "direct_dense_np", []

    if is_len:
        if np.any(lengths <= 0):
            raise ValueError("'frg_len' must be > 0 when using normalization.")

        weights = 1.0 / lengths
    else:
        weights = np.ones(starts.shape, dtype=np.float64)

    values = np.zeros(n_bins, dtype=np.float64)
    touched = np.zeros(n_bins, dtype=bool)
    start_bins = starts // siz_bin
    end_bins = (ends - 1) // siz_bin

    same = start_bins == end_bins

    if np.any(same):
        np.add.at(
            values,
            start_bins[same],
            (ends[same] - starts[same]) * weights[same],
        )
        touched[start_bins[same]] = True

    multi = ~same

    if np.any(multi):
        starts_subset = starts[multi]
        ends_subset = ends[multi]
        start_bins_subset = start_bins[multi]
        end_bins_subset = end_bins[multi]
        weights_subset = weights[multi]

        left_overlap = (
            (start_bins_subset + 1) * siz_bin - starts_subset
        ) * weights_subset
        right_overlap = (
            ends_subset - end_bins_subset * siz_bin
        ) * weights_subset
        np.add.at(values, start_bins_subset, left_overlap)
        np.add.at(values, end_bins_subset, right_overlap)
        touched[start_bins_subset] = True
        touched[end_bins_subset] = True

        has_interior = end_bins_subset > start_bins_subset + 1

        if np.any(has_interior):
            diff = np.zeros(n_bins + 1, dtype=np.float64)
            touched_diff = np.zeros(n_bins + 1, dtype=np.int64)
            first = start_bins_subset[has_interior] + 1
            last_exclusive = end_bins_subset[has_interior]
            value = siz_bin * weights_subset[has_interior]

            np.add.at(diff, first, value)
            np.add.at(diff, last_exclusive, -value)

            np.add.at(touched_diff, first, 1)
            np.add.at(touched_diff, last_exclusive, -1)

            # Interior coverage is exact in 'int64'; the 'float64' value sum is
            # not. Floating-point cancellation can leave ~1e-20 residue in bins
            # whose true value is zero, which 'values != 0.0' would emit as
            # data. Masking by the exact indicator restores zero.
            interior = np.cumsum(touched_diff[:-1]) > 0

            values += np.where(interior, np.cumsum(diff[:-1]), 0.0)
            touched |= interior

    if not np.any(touched & (values != 0.0)):
        return "direct_dense_np", []

    return "direct_dense_np", [(chrom, values, touched)]


def calc_sig_chrom_event_np(
    chrom: str,
    starts: np.ndarray,
    ends: np.ndarray,
    lengths: np.ndarray,
    chrom_size: int,
    siz_bin: int,
    is_len: bool,
) -> object:
    """
    Return bin edge and range-add events for parent-side materialization.

    Parameters
    ----------
    chrom : str
        Chromosome represented by the fragment arrays.
    starts, ends, lengths : np.ndarray
        Fragment starts, ends, and lengths in matching order.
    chrom_size : int
        Chromosome length used to clip events.
    siz_bin : int
        Signal-bin width.
    is_len : bool
        Whether fragments contribute inverse-length weights.

    Returns
    -------
    result : object
        Tagged chromosome result containing range-add event arrays.

    Raises
    ------
    ValueError
        If chromosome size or fragment lengths are invalid.
    """

    validate_comparison(siz_bin, "gt", 0, "siz_bin", allow_none=False)

    if chrom_size <= 0:
        raise ValueError(f"Chromosome size must be > 0 for {chrom!r}.")

    n_bins = math.ceil(chrom_size / siz_bin)
    if starts.size == 0:
        return "event_np", []

    valid = ends > starts

    if not np.all(valid):
        starts = starts[valid]
        ends = ends[valid]
        lengths = lengths[valid]

    if starts.size == 0:
        return "event_np", []

    if is_len:
        if np.any(lengths <= 0):
            raise ValueError("'frg_len' must be > 0 when using normalization.")

        weights = 1.0 / lengths
    else:
        weights = np.ones(starts.shape, dtype=np.float64)

    start_bins = starts // siz_bin
    end_bins = (ends - 1) // siz_bin

    edge_bin_parts = []
    edge_value_parts = []

    diff_bin_parts = []
    diff_value_parts = []

    touch_bin_parts = []
    touch_value_parts = []

    same = start_bins == end_bins

    if np.any(same):
        same_bins = start_bins[same].astype(np.int64, copy=False)
        edge_bin_parts.append(same_bins)
        edge_value_parts.append(
            ((ends[same] - starts[same]) * weights[same]).astype(
                np.float64,
                copy=False,
            ),
        )

    multi = ~same

    if np.any(multi):
        starts_subset = starts[multi]
        ends_subset = ends[multi]
        start_bins_subset = start_bins[multi]
        end_bins_subset = end_bins[multi]
        weights_subset = weights[multi]

        edge_bin_parts.extend(
            [
                start_bins_subset.astype(np.int64, copy=False),
                end_bins_subset.astype(np.int64, copy=False),
            ],
        )
        edge_value_parts.extend(
            [
                (
                    ((start_bins_subset + 1) * siz_bin - starts_subset)
                    * weights_subset
                ).astype(np.float64, copy=False),
                (
                    (ends_subset - end_bins_subset * siz_bin) * weights_subset
                ).astype(np.float64, copy=False),
            ],
        )
        has_interior = end_bins_subset > start_bins_subset + 1

        if np.any(has_interior):
            first = (start_bins_subset[has_interior] + 1).astype(
                np.int64,
                copy=False,
            )
            last_exclusive = end_bins_subset[has_interior].astype(
                np.int64,
                copy=False,
            )
            value = (siz_bin * weights_subset[has_interior]).astype(
                np.float64,
                copy=False,
            )
            diff_bin_parts.extend([first, last_exclusive])
            diff_value_parts.extend([value, -value])
            touch_bin_parts.extend([first, last_exclusive])
            touch_value_parts.extend(
                [
                    np.ones(first.shape, dtype=np.int64),
                    -np.ones(last_exclusive.shape, dtype=np.int64),
                ],
            )

    edge_bins = (
        np.concatenate(edge_bin_parts)
        if edge_bin_parts
        else np.empty(0, dtype=np.int64)
    )
    edge_values = (
        np.concatenate(edge_value_parts)
        if edge_value_parts
        else np.empty(0, dtype=np.float64)
    )
    diff_bins = (
        np.concatenate(diff_bin_parts)
        if diff_bin_parts
        else np.empty(0, dtype=np.int64)
    )
    diff_values = (
        np.concatenate(diff_value_parts)
        if diff_value_parts
        else np.empty(0, dtype=np.float64)
    )
    touch_bins = (
        np.concatenate(touch_bin_parts)
        if touch_bin_parts
        else np.empty(0, dtype=np.int64)
    )
    touch_values = (
        np.concatenate(touch_value_parts)
        if touch_value_parts
        else np.empty(0, dtype=np.int64)
    )

    if edge_bins.size == 0 and diff_bins.size == 0:
        return "event_np", []

    return (
        "event_np",
        [
            (
                chrom,
                n_bins,
                edge_bins,
                edge_values,
                diff_bins,
                diff_values,
                touch_bins,
                touch_values,
            ),
        ],
    )


def calc_sig_task(data: tuple[object, ...]) -> object:
    """
    Unpack one worker-task tuple and dispatch to 'calc_sig_chrom()'.
    """

    return calc_sig_chrom(*data)


def calc_sig_array_task(data: tuple[object, ...]) -> object:
    """
    Unpack one worker-task tuple and dispatch to 'calc_sig_chrom_array()'.
    """

    return calc_sig_chrom_array(*data)


def calc_sig_indexed_fetch_task(data: tuple[object, ...]) -> object:
    """
    Fetch one indexed region, parse fragments, and compute unscaled signal.

    Parameters
    ----------
    data : tuple[object, ...]
        Serialized indexed-fetch task arguments.

    Returns
    -------
    result : object
        Tagged task result for parent-side signal assembly.
    """

    (
        fil_aln,
        chrom,
        siz_chr,
        usr_frg,
        ref_fa,
        start,
        end,
        siz_bin,
        is_len,
        result_format,
    ) = data

    fragment_arrays, fragment_count = collect_frg_arr(
        iter_idx_frg(
            fil_aln=fil_aln,
            chrom=chrom,
            siz_chr=siz_chr,
            usr_frg=usr_frg,
            ref_fa=ref_fa,
            start=start,
            end=end,
        ),
    )

    arrays = fragment_arrays.get(chrom)

    if arrays is None:
        if result_format in DIRECT_SPARSE_RESULT_FORMATS:
            tag = (
                "direct_sparse_idx_np"
                if result_format == "direct_sparse_idx_np"
                else "direct_sparse_np"
            )
            sig = (tag, [])
        elif result_format in ("direct_dense_np", "event_np"):
            sig = (result_format, [])
        else:
            sig = {}
    elif result_format in DIRECT_SPARSE_RESULT_FORMATS:
        local_span = result_format in (
            "direct_sparse_local_np",
            "direct_sparse_local_bincount_np",
        )
        use_bincount = result_format in (
            "direct_sparse_bincount_np",
            "direct_sparse_touched_bincount_np",
            "direct_sparse_local_bincount_np",
        )
        emit_touched_only = result_format in (
            "direct_sparse_touched_np",
            "direct_sparse_touched_bincount_np",
        )
        sig = calc_sig_chrom_direct_sparse_np(
            chrom=chrom,
            starts=arrays[0],
            ends=arrays[1],
            lengths=arrays[2],
            chrom_size=siz_chr[chrom],
            siz_bin=siz_bin,
            is_len=is_len,
            return_bin_indices=result_format == "direct_sparse_idx_np",
            local_span=local_span,
            use_bincount=use_bincount,
            emit_touched_only=emit_touched_only,
        )
    elif result_format == "direct_dense_np":
        sig = calc_sig_chrom_direct_dense_np(
            chrom=chrom,
            starts=arrays[0],
            ends=arrays[1],
            lengths=arrays[2],
            chrom_size=siz_chr[chrom],
            siz_bin=siz_bin,
            is_len=is_len,
        )
    elif result_format == "event_np":
        sig = calc_sig_chrom_event_np(
            chrom=chrom,
            starts=arrays[0],
            ends=arrays[1],
            lengths=arrays[2],
            chrom_size=siz_chr[chrom],
            siz_bin=siz_bin,
            is_len=is_len,
        )
    else:
        sig = calc_sig_chrom_array(
            chrom=chrom,
            starts=arrays[0],
            ends=arrays[1],
            lengths=arrays[2],
            chrom_size=siz_chr[chrom],
            fragment_count=1,
            siz_bin=siz_bin,
            is_len=is_len,
            is_norm=False,
            scl_fct=None,
        )

    if result_format in (
        *DIRECT_SPARSE_RESULT_FORMATS,
        "direct_dense_np",
        "event_np",
    ):
        return sig, fragment_count

    return signal_dict_to_result(sig, result_format, siz_bin), fragment_count


def apply_sig_adj(
    signal_bins: dict[tuple[str, int], float],
    fragment_count: int,
    is_norm: bool,
    scl_fct: float | None = None,
) -> None:
    """
    Apply global depth normalization and scaling to already merged signal.

    Parameters
    ----------
    signal_bins : dict[tuple[str, int], float]
        Mutable sparse signal values.
    fragment_count : int
        Total accepted fragments used for normalization.
    is_norm : bool
        Whether to apply depth normalization.
    scl_fct : float | None
        Optional multiplicative scale factor.

    Raises
    ------
    ValueError
        If normalization is requested with no accepted fragments.
    """

    if is_array_signal(signal_bins):
        if is_norm:
            if fragment_count <= 0:
                raise ValueError(
                    "Normalization requires non-zero total fragments.",
                )

            for values in iter_sig_arr(signal_bins):
                values /= fragment_count

        if scl_fct is not None:
            validate_comparison(scl_fct, "gt", 0, "scl_fct")

            for values in iter_sig_arr(signal_bins):
                values *= scl_fct

        return

    if is_norm:
        if fragment_count <= 0:
            raise ValueError(
                "Normalization requires non-zero total fragments.",
            )

        for key in signal_bins:
            signal_bins[key] /= fragment_count

    if scl_fct is not None:
        validate_comparison(scl_fct, "gt", 0, "scl_fct")

        for key in signal_bins:
            signal_bins[key] *= scl_fct


def is_array_dense_signal(signal_bins: object) -> bool:
    """
    Return True for the merged dense-array signal container.
    """

    return (
        isinstance(signal_bins, tuple)
        and len(signal_bins) == 3
        and signal_bins[0] == "array_dense_merged"
    )


def is_array_signal(signal_bins: object) -> bool:
    """
    Return True for merged signal containers that avoid the final dict.
    """

    return is_sparse_sig(signal_bins) or is_array_dense_signal(
        signal_bins,
    )


def iter_sig_arr(signal_bins: Any) -> Iterator[np.ndarray]:
    """
    Yield mutable value arrays from merged array-backed signal containers.
    """

    if is_sparse_sig(signal_bins):
        for _, _, values in signal_bins[1]:
            yield values
    elif is_array_dense_signal(signal_bins):
        for _, values, _ in signal_bins[2]:
            yield values


def n_sig_rows(signal_bins: Any) -> int:
    """
    Count output bedGraph rows without assuming a dict-backed signal.
    """

    if is_sparse_sig(signal_bins):
        return sum(int(starts.size) for _, starts, _ in signal_bins[1])

    if is_array_dense_signal(signal_bins):
        return sum(
            int(np.count_nonzero(touched & (values != 0.0)))
            for _, values, touched in signal_bins[2]
        )

    return len(signal_bins)


def sparse_payload_to_starts(
    tag: str,
    coords: np.ndarray,
    siz_bin: int,
) -> np.ndarray:
    """
    Convert sparse payload coordinates to bedGraph bin starts.
    """

    if tag == "direct_sparse_idx_np":
        return coords.astype(np.int64, copy=False) * siz_bin

    return coords.astype(np.int64, copy=False)


def materialize_event_signal_parts(
    parts_by_chrom: dict[str, list[tuple]],
    siz_bin: int,
) -> object:
    """
    Combine event arrays and materialize final sparse signal arrays.

    Parameters
    ----------
    parts_by_chrom : dict[str, list[tuple]]
        Range-add event-array parts grouped by chromosome.
    siz_bin : int
        Signal-bin width used to derive output coordinates.

    Returns
    -------
    signal : object
        Ordered sparse signal arrays grouped by chromosome.
    """

    merged = []

    for chrom in sorted(parts_by_chrom, key=sort_chrom):
        n_bins = max(part[0] for part in parts_by_chrom[chrom])
        values = np.zeros(n_bins, dtype=np.float64)
        touched = np.zeros(n_bins, dtype=bool)

        edge_bins_parts = [part[1] for part in parts_by_chrom[chrom]]
        edge_values_parts = [part[2] for part in parts_by_chrom[chrom]]
        edge_bins = (
            np.concatenate(edge_bins_parts)
            if edge_bins_parts
            else (np.empty(0, dtype=np.int64))
        )
        edge_values = (
            np.concatenate(edge_values_parts)
            if edge_values_parts
            else np.empty(0, dtype=np.float64)
        )

        if edge_bins.size > 0:
            np.add.at(values, edge_bins, edge_values)
            touched[edge_bins] = True

        diff_bins_parts = [part[3] for part in parts_by_chrom[chrom]]
        diff_values_parts = [part[4] for part in parts_by_chrom[chrom]]
        diff_bins = (
            np.concatenate(diff_bins_parts)
            if diff_bins_parts
            else (np.empty(0, dtype=np.int64))
        )
        diff_values = (
            np.concatenate(diff_values_parts)
            if diff_values_parts
            else np.empty(0, dtype=np.float64)
        )

        touch_bins_parts = [part[5] for part in parts_by_chrom[chrom]]
        touch_values_parts = [part[6] for part in parts_by_chrom[chrom]]
        touch_bins = (
            np.concatenate(touch_bins_parts)
            if touch_bins_parts
            else (np.empty(0, dtype=np.int64))
        )
        touch_values = (
            np.concatenate(touch_values_parts)
            if touch_values_parts
            else np.empty(0, dtype=np.int64)
        )

        # Resolve the exact 'int64' interior indicator before the 'float64'
        # value sum that uses it. Floating-point cancellation can leave ~1e-20
        # residue in truly zero bins, which 'values != 0.0' would emit as data.
        # Masking by the indicator restores exact zero.
        interior = None

        if touch_bins.size > 0:
            touch_diff = np.zeros(n_bins + 1, dtype=np.int64)
            np.add.at(touch_diff, touch_bins, touch_values)
            interior = np.cumsum(touch_diff[:-1]) > 0
            touched |= interior

        if diff_bins.size > 0:
            diff = np.zeros(n_bins + 1, dtype=np.float64)
            np.add.at(diff, diff_bins, diff_values)
            summed = np.cumsum(diff[:-1])
            values += (
                summed if interior is None else np.where(interior, summed, 0.0)
            )

        idx = np.flatnonzero(touched & (values != 0.0))

        if idx.size > 0:
            merged.append(
                (
                    chrom,
                    (idx * siz_bin).astype(np.int64, copy=False),
                    values[idx].astype(np.float64, copy=False),
                ),
            )

    return "sparse_merged", merged


def _append_mapping_sparse_parts(
    signal_result: object,
    sparse_parts: dict[str, list[tuple[np.ndarray, np.ndarray]]],
) -> None:
    """
    Convert one mapping result into chromosome-keyed sparse arrays.

    Parameters
    ----------
    signal_result : object
        Mapping from chromosome/start keys to signal values.
    sparse_parts : dict[str, list[tuple[np.ndarray, np.ndarray]]]
        Mutable chromosome inventory receiving start and value arrays.
    """

    by_chrom: dict[str, list[tuple[int, float]]] = defaultdict(list)

    for (chrom, start), value in signal_result.items():
        by_chrom[chrom].append((start, value))

    for chrom, rows in by_chrom.items():
        sparse_parts[chrom].append(
            (
                np.asarray(
                    [start for start, _ in rows],
                    dtype=np.int64,
                ),
                np.asarray(
                    [value for _, value in rows],
                    dtype=np.float64,
                ),
            ),
        )


def _accumulate_sparse_merge_result(
    signal_result: object,
    sparse_parts: dict[str, list[tuple[np.ndarray, np.ndarray]]],
    siz_bin: int,
) -> None:
    """
    Append one worker result to sparse-merge chromosome parts.

    Mapping, sparse-array, and dense-array worker representations are
    normalized to chromosome-keyed start and value arrays.

    Parameters
    ----------
    signal_result : object
        Mapping- or tagged-tuple worker result.
    sparse_parts : dict[str, list[tuple[np.ndarray, np.ndarray]]]
        Mutable chromosome inventory receiving sparse arrays.
    siz_bin : int
        Signal-bin width used to convert bin indexes to starts.
    """

    if not isinstance(signal_result, tuple) or not signal_result:
        _append_mapping_sparse_parts(signal_result, sparse_parts)

        return

    tag = signal_result[0]

    if tag in (
        "sparse_np",
        "direct_sparse_np",
        "direct_sparse_idx_np",
    ):
        for chrom, coords, values in signal_result[1]:
            sparse_parts[chrom].append(
                (
                    sparse_payload_to_starts(tag, coords, siz_bin),
                    values,
                ),
            )

        return

    if tag == "dense_np":
        for chrom, first, values in signal_result[1]:
            idx = np.flatnonzero(values != 0.0)
            sparse_parts[chrom].append(
                (
                    first + idx.astype(np.int64) * siz_bin,
                    values[idx],
                ),
            )

        return

    _append_mapping_sparse_parts(signal_result, sparse_parts)


def _accumulate_dense_merge_result(
    signal_result: object,
    ensure_array: Callable[[str], np.ndarray],
    ensure_touched: Callable[[str], np.ndarray],
    siz_bin: int,
) -> None:
    """
    Accumulate one worker result into dense arrays and touch masks.

    Parameters
    ----------
    signal_result : object
        Mapping- or tagged-tuple worker result.
    ensure_array : Callable[[str], np.ndarray]
        Resolver for a chromosome's mutable dense signal array.
    ensure_touched : Callable[[str], np.ndarray]
        Resolver for a chromosome's mutable boolean touch mask.
    siz_bin : int
        Signal-bin width used to convert starts to array indexes.
    """

    if not isinstance(signal_result, tuple) or not signal_result:
        for (chrom, start), value in signal_result.items():
            ensure_array(chrom)[start // siz_bin] += value
            ensure_touched(chrom)[start // siz_bin] = True

        return

    tag = signal_result[0]

    if tag == "direct_dense_np":
        for chrom, values, touched in signal_result[1]:
            arr = ensure_array(chrom)
            arr[: values.size] += values
            ensure_touched(chrom)[: touched.size] |= touched

        return

    if tag == "dense_np":
        for chrom, first, values in signal_result[1]:
            arr = ensure_array(chrom)
            first_idx = first // siz_bin
            arr[first_idx : first_idx + values.size] += values
            ensure_touched(chrom)[first_idx : first_idx + values.size] |= (
                values != 0.0
            )

        return

    if tag in (
        "sparse_np",
        "direct_sparse_np",
        "direct_sparse_idx_np",
    ):
        for chrom, coords, values in signal_result[1]:
            if tag == "direct_sparse_idx_np":
                idx = coords.astype(np.int64, copy=False)
            else:
                idx = coords // siz_bin

            np.add.at(ensure_array(chrom), idx, values)
            ensure_touched(chrom)[idx] = True

        return

    for (chrom, start), value in signal_result.items():
        ensure_array(chrom)[start // siz_bin] += value
        ensure_touched(chrom)[start // siz_bin] = True


def _accumulate_event_merge_result(
    signal_result: object,
    event_parts: dict[str, list[tuple]],
    sparse_parts: dict[str, list[tuple[np.ndarray, np.ndarray]]],
    siz_bin: int,
) -> None:
    """
    Append one worker result to event or sparse fallback parts.

    Parameters
    ----------
    signal_result : object
        Mapping- or tagged-tuple worker result.
    event_parts : dict[str, list[tuple]]
        Mutable chromosome inventory receiving event-array payloads.
    sparse_parts : dict[str, list[tuple[np.ndarray, np.ndarray]]]
        Mutable chromosome inventory receiving fallback sparse arrays.
    siz_bin : int
        Signal-bin width used to convert bin indexes to starts.
    """

    if not isinstance(signal_result, tuple) or not signal_result:
        for (chrom, start), value in signal_result.items():
            sparse_parts[chrom].append(
                (
                    np.asarray([start], dtype=np.int64),
                    np.asarray([value], dtype=np.float64),
                ),
            )

        return

    tag = signal_result[0]

    if tag == "event_np":
        for (
            chrom,
            n_bins,
            edge_bins,
            edge_values,
            diff_bins,
            diff_values,
            touch_bins,
            touch_values,
        ) in signal_result[1]:
            event_parts[chrom].append(
                (
                    n_bins,
                    edge_bins,
                    edge_values,
                    diff_bins,
                    diff_values,
                    touch_bins,
                    touch_values,
                ),
            )

        return

    if tag in (
        "sparse_np",
        "direct_sparse_np",
        "direct_sparse_idx_np",
    ):
        for chrom, coords, values in signal_result[1]:
            sparse_parts[chrom].append(
                (
                    sparse_payload_to_starts(tag, coords, siz_bin),
                    values,
                ),
            )

        return

    for (chrom, start), value in signal_result.items():
        sparse_parts[chrom].append(
            (
                np.asarray([start], dtype=np.int64),
                np.asarray([value], dtype=np.float64),
            ),
        )


def _accumulate_dict_merge_result(
    signal_result: object,
    combined: defaultdict[tuple[str, int], float],
    siz_bin: int,
) -> None:
    """
    Convert and accumulate one worker result into a signal dictionary.

    Parameters
    ----------
    signal_result : object
        Mapping- or tagged-tuple worker result.
    combined : defaultdict[tuple[str, int], float]
        Mutable chromosome/start signal mapping.
    siz_bin : int
        Signal-bin width used to convert array indexes to starts.
    """

    if not isinstance(signal_result, tuple) or not signal_result:
        for key, value in signal_result.items():
            combined[key] += value

        return

    tag = signal_result[0]

    if tag in (
        "sparse_np",
        "direct_sparse_np",
        "direct_sparse_idx_np",
    ):
        for chrom, coords, values in signal_result[1]:
            starts = sparse_payload_to_starts(tag, coords, siz_bin)

            for start, value in zip(starts, values, strict=True):
                combined[(chrom, int(start))] += float(value)

        return

    if tag == "direct_dense_np":
        for chrom, values, touched in signal_result[1]:
            idx = np.flatnonzero(touched & (values != 0.0))

            for index in idx:
                combined[(chrom, int(index) * siz_bin)] += float(values[index])

        return

    if tag == "dense_np":
        for chrom, first, values in signal_result[1]:
            for offset, value in enumerate(values):
                if value != 0.0:
                    combined[(chrom, first + offset * siz_bin)] += float(value)

        return

    for key, value in signal_result.items():
        combined[key] += value


def _accumulate_array_merge_result(
    signal_result: object,
    ensure_array: Callable[[str], np.ndarray],
    siz_bin: int,
) -> None:
    """
    Accumulate one worker result into chromosome arrays.

    Parameters
    ----------
    signal_result : object
        Mapping- or tagged-tuple worker result.
    ensure_array : Callable[[str], np.ndarray]
        Resolver for a chromosome's mutable dense signal array.
    siz_bin : int
        Signal-bin width used to convert starts to array indexes.
    """

    if not isinstance(signal_result, tuple) or not signal_result:
        for (chrom, start), value in signal_result.items():
            ensure_array(chrom)[start // siz_bin] += value

        return

    tag = signal_result[0]

    if tag in (
        "sparse_np",
        "direct_sparse_np",
        "direct_sparse_idx_np",
    ):
        for chrom, coords, values in signal_result[1]:
            arr = ensure_array(chrom)

            if tag == "direct_sparse_idx_np":
                idx = coords.astype(np.int64, copy=False)
            else:
                idx = coords // siz_bin

            np.add.at(arr, idx, values)

        return

    if tag == "dense_np":
        for chrom, first, values in signal_result[1]:
            arr = ensure_array(chrom)
            first_idx = first // siz_bin
            arr[first_idx : first_idx + values.size] += values

        return

    for (chrom, start), value in signal_result.items():
        ensure_array(chrom)[start // siz_bin] += value


def merge_sig(
    results: Iterable[Any],
    merge_strategy: str,
    siz_chr: dict[str, int],
    siz_bin: int,
    is_prototype_strategy: bool,
    profile: dict[str, Any] | None = None,
) -> tuple[Any, int, int, int]:
    """
    Merge worker signal results and return signal, fragments, bins, and bytes.

    Parameters
    ----------
    results : Iterable[Any]
        Per-chromosome signal results to validate and merge.
    merge_strategy : str
        Named implementation used to combine worker results.
    siz_chr : dict[str, int]
        Chromosome lengths keyed by chromosome identifier.
    siz_bin : int
        Positive signal-bin width in base pairs.
    is_prototype_strategy : bool
        Whether the selected merge strategy is explicitly experimental.
    profile : dict[str, Any] | None
        Optional mutable timing and merge-statistics payload.

    Returns
    -------
    combined_signal, fragment_count, bin_count, payload_bytes : tuple[
        Any, int, int, int
    ]
        Mapping- or array-backed merged signal and execution statistics.

    Notes
    -----
    The first return value intentionally remains 'Any' because benchmark
    strategies use several tagged tuple representations in addition to the
    production signal mapping.
    """

    combined_signal: Any = defaultdict(float)
    combined_arrays: dict[str, np.ndarray] = {}
    array_touched: dict[str, np.ndarray] = {}
    sparse_parts: dict[str, list[tuple[np.ndarray, np.ndarray]]] = defaultdict(
        list,
    )
    event_parts: dict[str, list[tuple]] = defaultdict(list)

    fragment_count = 0
    bin_count = 0
    payload_bytes = 0

    meta_s = 0.0
    accum_s = 0.0
    array_to_dict_seconds = 0.0
    array_coalesce_seconds = 0.0
    dense_finalize_seconds = 0.0
    event_materialize_seconds = 0.0

    sparse_coalesce_stats: dict[str, int] = {}

    def ensure_array(chrom: str) -> np.ndarray:
        """
        Ensure a dense chromosome array exists for one chromosomal key.

        Parameters
        ----------
        chrom : str
            Chromosome name used as the lookup key.

        Returns
        -------
        array : np.ndarray
            Mutable float64 array sized to the chromosome span under the
            current bin width.
        """

        if chrom not in combined_arrays:
            combined_arrays[chrom] = np.zeros(
                math.ceil(siz_chr[chrom] / siz_bin),
                dtype=np.float64,
            )

        return combined_arrays[chrom]

    def ensure_touched(chrom: str) -> np.ndarray:
        """
        Ensure a boolean touch mask exists for one chromosomal key.

        Parameters
        ----------
        chrom : str
            Chromosome name used as the lookup key.

        Returns
        -------
        mask : np.ndarray
            Mutable boolean mask sized to the chromosome span under the current
            bin width.
        """

        if chrom not in array_touched:
            array_touched[chrom] = np.zeros(
                math.ceil(siz_chr[chrom] / siz_bin),
                dtype=bool,
            )

        return array_touched[chrom]

    for result in results:
        time_meta = time.perf_counter()

        if is_prototype_strategy:
            signal_result, result_fragment_count = result
            fragment_count += result_fragment_count
        else:
            signal_result = result

        bin_count += count_result_bins(signal_result)
        payload_bytes += est_result_bytes(signal_result)

        if profile is not None:
            meta_s += time.perf_counter() - time_meta

        time_accumulate = time.perf_counter()

        if merge_strategy in (
            "array_sparse_merge",
            "array_sparse_merge_legacy",
        ):
            _accumulate_sparse_merge_result(
                signal_result,
                sparse_parts,
                siz_bin,
            )
        elif merge_strategy == "array_dense_merge":
            _accumulate_dense_merge_result(
                signal_result,
                ensure_array,
                ensure_touched,
                siz_bin,
            )
        elif merge_strategy == "event_diff_merge":
            _accumulate_event_merge_result(
                signal_result,
                event_parts,
                sparse_parts,
                siz_bin,
            )
        elif merge_strategy == "dict_merge":
            _accumulate_dict_merge_result(
                signal_result,
                combined_signal,
                siz_bin,
            )
        else:
            _accumulate_array_merge_result(
                signal_result,
                ensure_array,
                siz_bin,
            )

        if profile is not None:
            accum_s += time.perf_counter() - time_accumulate

    dict_converting_merges = ("chrom_array_merge", "vectorized_merge")

    if merge_strategy in dict_converting_merges:
        time_convert = time.perf_counter()

        for chrom, array in combined_arrays.items():
            nonzero_indexes = np.flatnonzero(array != 0.0)

            for index in nonzero_indexes:
                combined_signal[(chrom, int(index) * siz_bin)] = float(
                    array[index],
                )

        if profile is not None:
            array_to_dict_seconds = time.perf_counter() - time_convert

    if merge_strategy in ("array_sparse_merge", "array_sparse_merge_legacy"):
        time_coalesce = time.perf_counter()
        combined_signal = coalesce_parts(
            sparse_parts,
            optimize=merge_strategy == "array_sparse_merge",
            stats=sparse_coalesce_stats,
        )

        if profile is not None:
            array_coalesce_seconds = time.perf_counter() - time_coalesce

    if merge_strategy == "array_dense_merge":
        time_dense = time.perf_counter()
        combined_signal = (
            "array_dense_merged",
            siz_bin,
            [
                (
                    chrom,
                    combined_arrays[chrom],
                    array_touched[chrom],
                )
                for chrom in sorted(combined_arrays, key=sort_chrom)
                if np.any(
                    array_touched[chrom] & (combined_arrays[chrom] != 0.0),
                )
            ],
        )

        if profile is not None:
            dense_finalize_seconds = time.perf_counter() - time_dense

    if merge_strategy == "event_diff_merge":
        time_event = time.perf_counter()
        combined_signal = materialize_event_signal_parts(
            event_parts,
            siz_bin,
        )

        if sparse_parts:
            _, sparse_merged = coalesce_parts(sparse_parts)
            combined_parts: dict[str, list[tuple[np.ndarray, np.ndarray]]] = (
                defaultdict(list)
            )

            for chrom, starts, values in combined_signal[1] + sparse_merged:
                combined_parts[chrom].append((starts, values))

            combined_signal = coalesce_parts(combined_parts)

        if profile is not None:
            event_materialize_seconds = time.perf_counter() - time_event

    if profile is not None:
        phases = profile["phases_s"]
        phases["merge_result_metadata"] = meta_s
        phases["merge_accumulate_results"] = accum_s
        phases["merge_array_to_dict"] = array_to_dict_seconds
        phases["merge_array_coalesce"] = array_coalesce_seconds
        phases["merge_array_dense_finalize"] = dense_finalize_seconds
        phases["merge_event_materialize"] = event_materialize_seconds
        profile["sparse_coalesce_stats"] = sparse_coalesce_stats

    return (
        combined_signal,
        fragment_count,
        bin_count,
        payload_bytes,
    )


def format_bdg_dict_rows(
    rows: list[tuple[tuple[str, int], float]],
    siz_bin: int,
    decimal_places: int,
    siz_chr: dict[str, int] | None,
) -> str:
    """
    Format a bedGraph row shard with the same row semantics as write_bdg().
    """

    lines = []

    for (chrom, bin_start), value in rows:
        bin_end = bin_start + siz_bin

        if siz_chr is not None:
            chrom_size = siz_chr.get(chrom)
            if chrom_size is None:
                raise ValueError(
                    f"Missing chromosome size for bedGraph row: {chrom!r}.",
                )

            if bin_start < 0 or bin_start >= chrom_size:
                raise ValueError(
                    "bedGraph bin start is outside chromosome bounds: "
                    f"{chrom}:{bin_start} (chromosome size {chrom_size}).",
                )

            bin_end = min(bin_end, chrom_size)

        lines.append(
            f"{chrom}\t{bin_start}\t{bin_end}\t"
            f"{format_bdg_value(value, decimal_places)}\n",
        )

    return "".join(lines)


def format_bdg_dict_rows_task(data: tuple[object, ...]) -> str:
    """
    ProcessPool-friendly wrapper for formatting one bedGraph shard.
    """

    return format_bdg_dict_rows(*data)


def iter_sig_chunks(
    signal_bins: Any,
    target_chunks: int,
) -> Iterator[tuple[str, np.ndarray, np.ndarray]]:
    """
    Yield ordered array-backed signal chunks for parallel formatting.

    Parameters
    ----------
    signal_bins : Any
        Supported array-backed signal representation.
    target_chunks : int
        Approximate number of ordered output chunks.

    Yields
    ------
    chunk : tuple[str, np.ndarray, np.ndarray]
        Chromosome with matching bin-index and value arrays.
    """

    total_rows = max(n_sig_rows(signal_bins), 1)
    chunk_size = max(1, math.ceil(total_rows / max(target_chunks, 1)))

    if is_sparse_sig(signal_bins):
        for chrom, starts, values in signal_bins[1]:
            for i in range(0, starts.size, chunk_size):
                yield (
                    chrom,
                    starts[i : i + chunk_size],
                    values[i : i + chunk_size],
                )
    elif is_array_dense_signal(signal_bins):
        dense_siz_bin = signal_bins[1]

        for chrom, values, touched in signal_bins[2]:
            idx = np.flatnonzero(touched & (values != 0.0))
            starts = (idx * dense_siz_bin).astype(np.int64, copy=False)

            for i in range(0, idx.size, chunk_size):
                chunk_starts = starts[i : i + chunk_size]
                chunk_values = values[idx[i : i + chunk_size]]

                yield (chrom, chunk_starts, chunk_values)


def digest_bdg_rows(
    coverage: Any,
    siz_bin: int,
    decimal_places: int,
    siz_chr: dict[str, int] | None,
) -> tuple[int, str]:
    """
    Return row count and SHA-256 for deterministic bedGraph-formatted rows.

    Parameters
    ----------
    coverage : Any
        Supported signal representation.
    siz_bin : int
        Signal-bin width.
    decimal_places : int
        Decimal precision for rendered values.
    siz_chr : dict[str, int] | None
        Optional chromosome sizes used to clip final intervals.

    Returns
    -------
    row_count, digest : tuple[int, str]
        Rendered row count and hexadecimal SHA-256 digest.
    """

    if is_array_signal(coverage):
        digest = hashlib.sha256()
        row_count = 0

        for chrom, starts, values in iter_sig_chunks(coverage, 1):
            text = format_bdg_rows(
                chrom,
                starts,
                values,
                siz_bin,
                decimal_places,
                siz_chr.get(chrom) if siz_chr is not None else None,
            )
            row_count += text.count("\n")
            digest.update(text.encode("utf-8"))

        return row_count, digest.hexdigest()

    digest = hashlib.sha256()
    row_count = 0
    items = sorted(
        coverage.items(),
        key=lambda item: key_bin(item[0][0], item[0][1]),
    )
    chunk_size = 100000

    for index in range(0, len(items), chunk_size):
        text = format_bdg_dict_rows(
            items[index : index + chunk_size],
            siz_bin,
            decimal_places,
            siz_chr,
        )
        row_count += text.count("\n")
        digest.update(text.encode("utf-8"))

    return row_count, digest.hexdigest()


def write_bdg_sparse(
    signal_bins: Any,
    fil_out: str,
    siz_bin: int,
    decimal_places: int,
    siz_chr: dict[str, int] | None,
) -> None:
    """
    Write a sparse-array signal directly to bedGraph.
    """

    with open_out(fil_out) as handle:
        for chrom, starts, values in iter_sig_chunks(signal_bins, 1):
            handle.write(
                format_bdg_rows(
                    chrom,
                    starts,
                    values,
                    siz_bin,
                    decimal_places,
                    siz_chr.get(chrom) if siz_chr is not None else None,
                ),
            )


def write_bdg_par(
    coverage: Any,
    fil_out: str,
    siz_bin: int,
    decimal_places: int,
    siz_chr: dict[str, int] | None,
    workers: int,
) -> None:
    """
    Format sorted bedGraph shards in worker processes and write in order.

    Parameters
    ----------
    coverage : Any
        Supported signal representation.
    fil_out : str
        Output bedGraph path.
    siz_bin : int
        Signal-bin width.
    decimal_places : int
        Decimal precision for rendered values.
    siz_chr : dict[str, int] | None
        Optional chromosome sizes used to clip final intervals.
    workers : int
        Number of formatting worker processes.
    """

    validate_comparison(
        workers,
        "ge",
        1,
        "wrk_writer",
        allow_none=False,
    )

    if is_array_signal(coverage):
        if workers == 1:
            write_bdg_sparse(
                coverage,
                fil_out,
                siz_bin,
                decimal_places,
                siz_chr,
            )

            return

        task_data = [
            (
                chrom,
                starts,
                values,
                siz_bin,
                decimal_places,
                siz_chr.get(chrom) if siz_chr is not None else None,
            )
            for chrom, starts, values in iter_sig_chunks(
                coverage,
                workers,
            )
        ]

        try:
            with ProcessPoolExecutor(
                max_workers=min(workers, os.cpu_count() or 1),
            ) as ex:
                formatted = list(ex.map(format_bdg_rows_task, task_data))
        except OSError:
            write_bdg_sparse(
                coverage,
                fil_out,
                siz_bin,
                decimal_places,
                siz_chr,
            )

            return

        with open_out(fil_out) as handle:
            for shard in formatted:
                handle.write(shard)

        return

    if workers == 1:
        write_bdg(
            coverage,
            fil_out,
            siz_bin,
            decimal_places,
            siz_chr=siz_chr,
        )

        return

    items = sorted(
        coverage.items(),
        key=lambda item: key_bin(item[0][0], item[0][1]),
    )

    if not items:
        with open_out(fil_out):
            return

    chunk_size = math.ceil(len(items) / workers)
    chunks = [
        items[index : index + chunk_size]
        for index in range(0, len(items), chunk_size)
    ]
    task_data = [(chunk, siz_bin, decimal_places, siz_chr) for chunk in chunks]

    try:
        with ProcessPoolExecutor(
            max_workers=min(workers, os.cpu_count() or 1),
        ) as ex:
            formatted = list(ex.map(format_bdg_dict_rows_task, task_data))
    except OSError:
        write_bdg(
            coverage,
            fil_out,
            siz_bin,
            decimal_places,
            siz_chr=siz_chr,
        )

        return

    with open_out(fil_out) as handle:
        for shard in formatted:
            handle.write(shard)
