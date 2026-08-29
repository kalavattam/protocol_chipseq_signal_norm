#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: utils_stabilizer.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5-series models; most recent: GPT-5.6);
# - Anthropic Claude Code (Opus 5).
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

# Aliases for normalized coverage, matching 'compute_signal.METHOD_CANON' so
# one substrate has one vocabulary across the codebase. 'norm' is canonical;
# 'nc' is shorthand used in prose and is accepted rather than privileged.
NORM_CANON = {
    "CPM": "CPM",
    "BPM": "BPM",
    "RPKM": "RPKM",
    "None": "None",
    "RPGC": "RPGC",
    "n": "norm",
    "nc": "norm",
    "nrm": "norm",
    "norm": "norm",
    "normalized": "norm",
}
NORM_CHOICES = tuple(NORM_CANON.keys())


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


def canonicalize_norm(norm: str) -> str:
    """
    Map a normalization alias onto its canonical name.

    Parameters
    ----------
    norm : str
        Any key of 'NORM_CANON'.

    Returns
    -------
    canonical : str
        The canonical name.

    Raises
    ------
    ValueError
        If 'norm' is not a recognized alias.
    """

    if norm not in NORM_CANON:
        raise ValueError(f"Error: Unknown --normalization: {norm!r}")

    return NORM_CANON[norm]


def compute_pseudo_edger(
    lib_a: float,
    lib_b: float,
    prior_count: float = 2.0,
    norm: str = "CPM",
    siz_bin: int = 10,
    scale_a: float | None = None,
    scale_b: float | None = None,
    frg_a: float | None = None,
    frg_b: float | None = None,
) -> dict[str, object]:
    """
    Derive edgeR-equivalent scale factors and pseudocounts for deepTools.

    Parameters
    ----------
    lib_a : float
        Library size for track A, as edgeR's 'lib.size'.
    lib_b : float
        Library size for track B.
    prior_count : float
        edgeR's 'prior.count' before library-size scaling.
    norm : str
        Target deepTools normalization: 'CPM', 'BPM', 'RPKM', 'None', or
        'RPGC'.
    siz_bin : int
        Bin width in base pairs; used by 'RPKM'.
    scale_a, scale_b : float | None
        Externally supplied deepTools scale factors, required for 'RPGC'.
    frg_a, frg_b : float | None
        Fragment counts, required for 'norm' (normalized coverage). This is the
        same denominator 'compute_signal' divides by, not the number of
        alignment records; the two coincide only when one record per fragment
        survives filtering. Not derivable from a normalized-coverage track,
        which sums to 1 by construction.

    Returns
    -------
    result : dict[str, object]
        Keys 'scale_A', 'scale_B', 'pseudo_A', 'pseudo_B', 'prior_scaled_A',
        'prior_scaled_B', 'is_edger', and 'note'. Normalized coverage adds
        'k_A' and 'k_B'; every other mode omits them, so a consumer tests for
        the key rather than assuming it.

    Raises
    ------
    ValueError
        For a nonpositive library size, a negative 'prior.count', an unknown
        'norm', 'RPGC' without both scale factors, or normalized coverage
        without both fragment counts.

    Notes
    -----
    edgeR's 'cpm(log=TRUE)' computes '(y_i + y0_i) / (L_i + 2 * y0_i) * 1e6'
    with 'y0_i = prior.count * L_i / mean(L)'. That is linear in 'y_i', so it
    splits into a slope and an intercept that deepTools can express as a scale
    factor and a pseudocount:

        s_i = 1e6 / (L_i + 2 * y0_i)
        p_i = s_i * y0_i = 1e6 * prior.count / (mean(L) + 2 * prior.count)

    'p_i' carries no per-sample index, so edgeR's rule is symmetric in
    normalized units. The scale factor is not deepTools' CPM factor '1e6 / N',
    so '--normalizeUsing' cannot reach it; the pair must be passed as
    '--scaleFactors A:B --pseudocount P P'.

    Symmetry holds in real arithmetic, not in 'float64': 'p_A' and 'p_B' reach
    the same value by different routes, so they can differ in the last bits.
    Compare them with a tolerance, never with '=='.

    'prior_scaled_A' and 'prior_scaled_B' are 'y0_i' itself: the per-sample
    prior in count space, before 's_i' converts it to the output substrate.
    Nothing downstream consumes them; they exist to be read. Their ratio
    'prior_scaled_A / prior_scaled_B' is 'L_A / L_B', the depth imbalance in
    prior units; their sum is '2 * prior.count' always, so a violation means a
    library size is wrong; and 'pseudo_i' factors as 's_i * prior_scaled_i',
    the decomposition the pseudocount alone hides: a tiny value may come from
    the depth correction or from the scale factor, and only the pair tells them
    apart.

    Under 'norm' they scale by 'frg_i / mean(frg)' rather than by
    'L_i / mean(L)', so there they report fragment imbalance, not bin-sum
    imbalance. That is also the one mode where they do not feed 'pseudo_i',
    which is symmetric by construction, so they carry the only per-sample
    information it returns.

    A single track is expressed by passing its library size as both 'lib_a' and
    'lib_b'. That is not an approximation: edgeR averages library sizes over
    all columns and scales each prior by 'L_i / mean(L)'
    ('add_prior_count.c:88'), so with one column the ratio is exactly 1, the
    prior stays nominal, and 'mean(L)' is that column's own library size.
    Passing 'L' twice reproduces both. Because 'mean(L)' differs between the
    two framings, a track's one-track pseudocount is not its two-track
    pseudocount; that is edgeR's behavior too, not an artifact here.

    'norm' (aliases 'nc', 'n', 'nrm', 'normalized') is normalized coverage:
    each fragment deposits exactly 1.0 across its footprint and the track is
    divided by the fragment count, so it sums to 1, not to a library size. Its
    own total is therefore useless as a denominator; the library size that
    matters is 'N', the fragment count that 'compute_signal' divided by, which
    must be supplied. That is not the alignment-record count unless exactly one
    record per fragment survives filtering, which is what '--samFlagInclude 64'
    arranges for paired-end data and what its absence undoes.

    Normalized coverage also needs a correction edgeR does not make. A fragment
    spanning 'k' bins deposits '1/k' into each, so an 'nc' bin is a sum of
    fractional shares rather than a count of events, and is under-dispersed by
    about 'k' (measured 'Var/E' 0.109 against 1.89 for counts). A prior
    calibrated on Poisson counts is therefore about 'k'-fold too strong, giving

        p_nc = prior.count / (k_bar * N_bar), k = L / N

    with 'k_bar' and 'N_bar' averaged over the pair. The 'prior.count' and the
    '1/N_bar' come from edgeR; the '1/k_bar' does not, which is why 'is_edger'
    is False for this mode.

    'None' and 'RPGC' also do not reproduce the estimator and are likewise
    reported with 'is_edger' False. edgeR adjusts the denominator to
    'L_i + 2 * y0_i'; RPGC's denominator is 'N * F / G', whose meaning is
    one-fold genome coverage, and substituting a bin-matrix column sum for an
    alignment count makes that meaning false. Both apply the prior's magnitude
    only, in the proportional form 'p_i = s_i * y0_i'.
    """

    for label, lib in (("lib_a", lib_a), ("lib_b", lib_b)):
        if not math.isfinite(lib) or lib <= 0.0:
            raise ValueError(f"{label!r} must be finite and positive.")

    if not math.isfinite(prior_count) or prior_count < 0.0:
        raise ValueError("'prior_count' must be finite and nonnegative.")

    norm = canonicalize_norm(norm)

    if norm == "norm":
        for label, frg in (("frg_a", frg_a), ("frg_b", frg_b)):
            if frg is None or not math.isfinite(frg) or frg <= 0.0:
                raise ValueError(
                    f"{label!r} must be finite and positive for 'norm'; it is "
                    "the fragment count, which a normalized-coverage track "
                    "cannot supply because it sums to 1.",
                )

        k_a = lib_a / frg_a
        k_b = lib_b / frg_b
        k_mean = 0.5 * (k_a + k_b)
        frg_mean = 0.5 * (frg_a + frg_b)
        pseudo = prior_count / (k_mean * frg_mean)

        return {
            "scale_A": 1.0,
            "scale_B": 1.0,
            "pseudo_A": pseudo,
            "pseudo_B": pseudo,
            "prior_scaled_A": prior_count * frg_a / frg_mean,
            "prior_scaled_B": prior_count * frg_b / frg_mean,
            "k_A": k_a,
            "k_B": k_b,
            "is_edger": False,
            "note": (
                "edgeR's prior divided by k_bar for the under-dispersion of "
                "normalized coverage, a correction edgeR does not make"
            ),
        }

    lib_mean = 0.5 * (lib_a + lib_b)
    prior_a = prior_count * lib_a / lib_mean
    prior_b = prior_count * lib_b / lib_mean

    if norm in ("CPM", "BPM"):
        scale_a = 1e6 / (lib_a + 2.0 * prior_a)
        scale_b = 1e6 / (lib_b + 2.0 * prior_b)
        is_edger = True
        note = "exact; BPM is CPM in deepTools"
    elif norm == "RPKM":
        if siz_bin <= 0:
            raise ValueError("'siz_bin' must be positive for 'RPKM'.")

        scale_a = 1e9 / ((lib_a + 2.0 * prior_a) * siz_bin)
        scale_b = 1e9 / ((lib_b + 2.0 * prior_b) * siz_bin)
        is_edger = True
        note = "exact"
    elif norm == "None":
        scale_a = 1.0
        scale_b = 1.0
        is_edger = False
        note = (
            "reproduces edgeR's ratio up to a constant log offset, since the "
            "denominator adjustment is absent"
        )
    else:
        if scale_a is None or scale_b is None:
            raise ValueError(
                "'RPGC' requires both 'scale_a' and 'scale_b'; read them from "
                "'bamCoverage --verbose' as the final scaling factor.",
            )

        is_edger = False
        note = (
            "prior magnitude only, in the proportional form; RPGC's one-fold "
            "denominator has no edgeR analog"
        )

    return {
        "scale_A": scale_a,
        "scale_B": scale_b,
        "pseudo_A": scale_a * prior_a,
        "pseudo_B": scale_b * prior_b,
        "prior_scaled_A": prior_a,
        "prior_scaled_B": prior_b,
        "is_edger": is_edger,
        "note": note,
    }


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
