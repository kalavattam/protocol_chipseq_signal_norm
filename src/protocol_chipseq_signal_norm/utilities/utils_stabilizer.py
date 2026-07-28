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
    Yield filtered bedGraph values from column four.

    Epsilon handling and sign policy determine which finite values remain.

    Parameters
    ----------
    path : str
        BedGraph-like path, or '-' for standard input.
    eps : float
        Epsilon threshold.
    mode_nz : str
        Boundary mode: 'closed', 'open', or 'off'.
    skp_pfx : tuple[str, ...] | None
        Header prefixes to skip. The default is 'DEF_SKP_PFX'.
    nz_policy : str
        Sign policy: 'abs', 'pos', or 'off'.

    Yields
    ------
    value : float
        Finite values retained by the selected policies.

    Raises
    ------
    ValueError
        If 'mode_nz' or 'nz_policy' is unknown.
    """

    if skp_pfx is None:
        skp_pfx = DEF_SKP_PFX

    def skip_predicate(line: str) -> bool:
        """
        Return whether a raw input row is non-data content.
        """

        return is_header(line, skp_pfx)

    def is_zero(v: float) -> bool:
        """
        Return whether a value is zero under the active policy.
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

    with open_in(path) as handle:
        for (
            _chromosome,
            _start,
            _end,
            _token,
            numeric_value,
        ) in iter_rows_bdg(handle, skip_predicate):
            if numeric_value is None:
                continue

            if not math.isfinite(numeric_value):
                continue

            if is_zero(numeric_value):
                continue

            yield numeric_value


def compute_stats_robust(
    values: Iterable[float],
) -> dict[str, float | int]:
    """
    Return simple robust summary statistics for finite values.

    Parameters
    ----------
    values : Iterable[float]
        Input values. Non-finite values are ignored.

    Returns
    -------
    statistics : dict[str, float | int]
        Dictionary with keys:
            - 'n'
            - 'median'
            - 'mean'

        If no finite values remain, returns:
            '{"n": 0, "median": nan, "mean": nan}'.

    Notes
    -----
    This function materializes 'values' into a list.
    """

    finite_values = [value for value in values if math.isfinite(value)]

    if not finite_values:
        return {"n": 0, "median": float("nan"), "mean": float("nan")}

    finite_values.sort()
    count = len(finite_values)

    if count % 2:
        median = finite_values[count // 2]
    else:
        median = 0.5 * (
            finite_values[count // 2 - 1] + finite_values[count // 2]
        )

    mean = sum(finite_values) / count

    return {"n": count, "median": median, "mean": mean}


def determine_coef_eff(method: str, coef: float | None) -> float | None:
    """
    Resolve effective coefficient for methods that use one.

    If coef is None:
        - frc_mdn_nz, frc_avg_nz -> 0.01
        - min_nz                 -> 1.0
        - qntl_nz                -> None (no coefficient)

    Parameters
    ----------
    method : str
        Selected stabilization method.
    coef : float | None
        Optional coefficient override for the selected method.

    Returns
    -------
    coefficient : float | None
        The effective method-specific coefficient.

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


def median_sorted(values: list[float]) -> float:
    """
    Return the median of an already sorted numeric list.

    Parameters
    ----------
    values : list[float]
        Sorted numeric values. Must be non-empty.

    Returns
    -------
    median : float
        Median of 'values'.

    Raises
    ------
    IndexError
        If 'values' is empty.

    Notes
    -----
    This helper assumes the caller has already sorted the values.
    """

    count = len(values)
    if count % 2:
        return values[count // 2]

    return 0.5 * (values[count // 2 - 1] + values[count // 2])


def pick_stabilizer(
    values: Iterable[float],
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
    values : Iterable[float]
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
    value : float
        Selected stabilizer, which may be NaN if filtering removes all finite
        values.

    Raises
    ------
    ValueError for invalid method, invalid qntl_pct, invalid qntl_rule.
    """

    finite_values = [value for value in values if math.isfinite(value)]
    if not finite_values:
        return float("nan")

    finite_values.sort()

    if coef is None:
        coef = determine_coef_eff(method, None)

    if method == "frc_mdn_nz":
        return max(
            (coef if coef is not None else 0.01)
            * median_sorted(finite_values),
            floor,
        )

    if method == "frc_avg_nz":
        mean = sum(finite_values) / len(finite_values)
        return max((coef if coef is not None else 0.01) * mean, floor)

    if method == "min_nz":
        return max(
            (coef if coef is not None else 1.0) * finite_values[0],
            floor,
        )

    if method == "qntl_nz":
        if not math.isfinite(qntl_pct) or not (0.0 <= qntl_pct <= 100.0):
            raise ValueError("Error: qntl_pct must be finite and in [0, 100].")

        p = qntl_pct / 100.0

        if qntl_rule == "round":
            index = round(p * (len(finite_values) - 1))
        elif qntl_rule == "floor":
            index = math.floor(p * (len(finite_values) - 1))
        else:
            raise ValueError(f"Error: Unknown qntl_rule: {qntl_rule!r}")

        index = max(0, min(len(finite_values) - 1, index))

        return max(finite_values[index], floor)

    raise ValueError(f"Error: Unknown method: {method!r}")
