#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: utils_stabilizer.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5-series models; most recent: GPT-5.6) were
# used in development and documentation.
#
# Distributed under the MIT license.


"""
Stabilizer-statistic helpers for ratio workflows.

Notes
-----
Helpers support pseudocount and denominator-floor computations used by
'compute_pseudo.py' and related scripts.
"""

from __future__ import annotations

import math
import sys
from collections.abc import Iterable, Iterator

from protocol_chipseq_signal_norm.utilities.utils_bdg import iter_rows_bdg
from protocol_chipseq_signal_norm.utilities.utils_io import (
    DEF_SKP_PFX,
    is_header,
    open_in,
)

assert sys.version_info >= (3, 11), "Python >= 3.11 required."


def iter_vals_bdg(
    path: str,
    eps: float = 0.0,
    mode_nz: str = "closed",
    skp_pfx: tuple[str, ...] | None = None,
    nz_policy: str = "abs",
) -> Iterator[float]:
    """
    Yield bedGraph values (col 4) from a file, with configurable “nonzero”
    filtering.

    Parameters
    ----------
        path : str
            bedGraph-like file path ('.gz' handled), or '-' for stdin.
        eps : float
            Epsilon threshold.
        mode_nz : str
            {'closed', 'open', 'off'} meaning:
                - closed: treat “zero” as <= eps (or |v| <= eps under abs)
                - open:   treat “zero” as <  eps (or |v| <  eps under abs)
                - off:    no epsilon-based filtering
        skp_pfx : tuple[str, ...] | None
            Header prefixes to skip; defaults to DEF_SKP_PFX.
        nz_policy : str
            How “nonzero” is defined:
                - "abs": use |v| vs eps (good for pseudocounts)
                - "pos": use v vs eps (good for denominator floors)
                - "off": do not drop based on sign/eps (but drop nonfinite)

    Yields:
        v_num : float
            float values (finite only).

    Raises
    ------
        ValueError on unknown mode_nz / nz_policy.
    """
    if skp_pfx is None:
        skp_pfx = DEF_SKP_PFX

    def skp_prd(line: str) -> bool:
        """
        Predicate used by 'iter_vals_bdg()' to decide whether a raw input line
        should be skipped (e.g., header, metadata, or blank line).
        """
        return is_header(line, skp_pfx)

    def is_zero(v: float) -> bool:
        """
        Return True when a numeric value is treated as zero under the current
        epsilon/nonzero policy; otherwise False.
        """
        if mode_nz == "off" or nz_policy == "off":
            return False

        if nz_policy == "abs":
            x = abs(v)
        elif nz_policy == "pos":
            x = v
        else:
            raise ValueError(f"Unknown nz_policy: {nz_policy!r}")

        if mode_nz == "closed":
            return x <= eps
        if mode_nz == "open":
            return x < eps

        raise ValueError(f"Unknown mode_nz: {mode_nz!r}")

    with open_in(path) as fh:
        for _chrom, _s, _e, _tok, v_num in iter_rows_bdg(fh, skp_prd):
            if v_num is None:
                continue
            if not math.isfinite(v_num):
                continue
            if is_zero(v_num):
                continue
            yield v_num


def compute_stats_robust(vals: Iterable[float]) -> dict[str, float | int]:
    """
    Return simple robust summary statistics for finite values.

    Parameters
    ----------
        vals : Iterable[float]
            Input values. Non-finite values are ignored.

    Returns
    -------
        dict[str, float | int]
            Dictionary with keys:
                - 'n'
                - 'median'
                - 'mean'

            If no finite values remain, returns:
                '{"n": 0, "median": nan, "mean": nan}'.

    Raises
    ------
        None.

    Notes
    -----
        This function materializes 'vals' into a list.
    """
    xs = [v for v in vals if math.isfinite(v)]
    if not xs:
        return {"n": 0, "median": float("nan"), "mean": float("nan")}

    xs.sort()
    n = len(xs)

    if n % 2:
        median = xs[n // 2]
    else:
        median = 0.5 * (xs[n // 2 - 1] + xs[n // 2])

    mean = sum(xs) / n
    return {"n": n, "median": median, "mean": mean}


def determine_coef_eff(method: str, coef: float | None) -> float | None:
    """
    Resolve effective coefficient for methods that use one.

    If coef is None:
        - frc_mdn_nz, frc_avg_nz -> 0.01
        - min_nz                 -> 1.0
        - qntl_nz                -> None (no coefficient)

    Raises
    ------
        ValueError
            If 'method' is unrecognized.
    """
    if coef is not None:
        return coef

    if method in ("frc_mdn_nz", "frc_avg_nz"):
        return 0.01

    if method == "min_nz":
        return 1.0

    if method == "qntl_nz":
        return None

    raise ValueError(f"Unknown stabilizer-selection method: {method!r}")


def median_sorted(xs: list[float]) -> float:
    """
    Return the median of an already sorted numeric list.

    Parameters
    ----------
        xs : list[float]
            Sorted numeric values. Must be non-empty.

    Returns
    -------
        float
            Median of 'xs'.

    Raises
    ------
        IndexError
            If 'xs' is empty.

    Notes
    -----
        This helper assumes the caller has already sorted the values.
    """
    n = len(xs)
    if n % 2:
        return xs[n // 2]
    return 0.5 * (xs[n // 2 - 1] + xs[n // 2])


def pick_stabilizer(
    vals: Iterable[float],
    method: str,
    coef: float | None = None,
    qntl_pct: float = 1.0,
    floor: float = 0.0,
    qntl_rule: str = "round",
) -> float:
    """
    Choose a stabilizer value from a collection of values.

    Parameters
    ----------
        vals : Iterable[float]
            Iterable of values; non-finite values are ignored.
        method : str
            One of {'frc_mdn_nz', 'qntl_nz', 'frc_avg_nz', 'min_nz'}.
        coef : float | None
            Coefficient for the frc_* and min_* methods. If None, resolved by
            determine_coef_eff().
        qntl_pct : float
            Quantile in percent for qntl_nz (0..100). Decimals allowed.
        floor : float
            Lower bound applied to the result: max(value, floor).
        qntl_rule : str
            {'round', 'floor'} selection rule on sorted values:
                k = round(p*(n-1))  or  k = floor(p*(n-1))

    Returns
    -------
        float (may be nan if no finite values exist after filtering).

    Raises
    ------
        ValueError for invalid method, invalid qntl_pct, invalid qntl_rule.
    """
    xs = [v for v in vals if math.isfinite(v)]
    if not xs:
        return float("nan")
    xs.sort()

    if coef is None:
        coef = determine_coef_eff(method, None)

    if method == "frc_mdn_nz":
        return max((coef if coef is not None else 0.01) * median_sorted(xs),
                   floor)

    if method == "frc_avg_nz":
        mean = sum(xs) / len(xs)
        return max((coef if coef is not None else 0.01) * mean, floor)

    if method == "min_nz":
        return max((coef if coef is not None else 1.0) * xs[0], floor)

    if method == "qntl_nz":
        if not math.isfinite(qntl_pct) or not (0.0 <= qntl_pct <= 100.0):
            raise ValueError("Error: qntl_pct must be finite and in [0, 100].")

        p = qntl_pct / 100.0

        if qntl_rule == "round":
            k = round(p * (len(xs) - 1))
        elif qntl_rule == "floor":
            k = math.floor(p * (len(xs) - 1))
        else:
            raise ValueError(f"Error: Unknown qntl_rule: {qntl_rule!r}")

        k = max(0, min(len(xs) - 1, k))
        return max(xs[k], floor)

    raise ValueError(f"Error: Unknown method: {method!r}")
