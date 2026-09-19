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
Helpers support the pseudocount and denominator-floor computations used by
'compute_pseudo.py', 'compute_pseudo_deeptools.py' and
'compute_input_floor.py'.
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

# The union of both tools' substrates (i.e., types of signal, adjusted or not);
# each CLI restricts its own 'choices' to the subset it serves. Spellings
# mirror 'compute_signal.METHOD_CANON'.
# fmt: off
SUBSTRATE_CANON = {
    # Fractional deposition: a partly covered bin takes a share proportional to
    # the overlap.
    "unadj": "unadj",
    "frag": "frag",
    "norm": "norm",
    "nc": "norm",

    # Whole-count deposition: one count per touched bin.
    "count": "count",
    "cpm": "cpm",

    # Written by deepTools '--normalizeUsing'.
    "CPM": "CPM",
    "BPM": "BPM",
    "RPKM": "RPKM",
    "None": "None",
    "RPGC": "RPGC",
}
# fmt: on
SUBSTRATE_CHOICES = tuple(SUBSTRATE_CANON.keys())
# TODO: Here and elsewhere, in code, docs, API, comments, docstrings, etc.,
# rename 'substrate' to something clearer, e.g., 'typ_sig'; ergo, first poll
# if/how 'typ_sig' is already being used across the codebase.


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


def canonicalize_substrate(substrate: str) -> str:
    """
    Map a substrate alias onto its canonical name.

    Parameters
    ----------
    substrate : str
        Any key of 'SUBSTRATE_CANON'.

    Returns
    -------
    canonical : str
        The canonical name.

    Raises
    ------
    ValueError
        If 'substrate' is not a recognized alias.
    """

    if substrate not in SUBSTRATE_CANON:
        raise ValueError(f"Error: Unknown substrate: {substrate!r}")

    return SUBSTRATE_CANON[substrate]


# TODO: Here and elsewhere across the codebase (code, comments, docs, etc.),
# 'n_bin_(a|b)' becomes 'n_ovlp_(a|b)', and 'substrate' potentially becomes
# 'typ_sig' or, if/when appropriate in prose, "type of signal", "signal type",
# and the like.
def compute_pseudo_edger(
    substrate: str = "norm",
    prior_count: float = 2.0,
    siz_bin: int | None = None,
    n_bin_a: float | None = None,
    n_bin_b: float | None = None,
    n_frg_a: float | None = None,
    n_frg_b: float | None = None,
    total_a: float | None = None,
    total_b: float | None = None,
    scale_a: float | None = None,
    scale_b: float | None = None,
) -> dict[str, object]:
    """
    Derive edgeR-equivalent scale factors and pseudocounts for a substrate.

    Parameters
    ----------
    substrate : str
        Target substrate. This project's are the fractional 'unadj', 'frag' and
        'norm', and the whole-count 'count' and 'cpm'; deepTools' are 'CPM',
        'BPM', 'RPKM', 'RPGC' and 'None'. Each CLI restricts its own accepted
        subset, so this function spans both. Defaults to 'norm', but every call
        in this project passes it explicitly.
    prior_count : float
        edgeR's 'prior.count', before the per-sample scaling each substrate
        applies.
    siz_bin : int | None
        Bin width in base pairs, required for 'RPKM' and unused otherwise. No
        default: a wrong width rescales an 'RPKM' pair in silence.
    n_bin_a, n_bin_b : float | None
        Required for every substrate. Fragment-bin overlap count 'L' for each
        track, the bin-matrix column sum and edgeR's 'lib.size': how many bins
        each fragment spans, added up over fragments. Not the fragment count,
        which is 'n_frg_a'.
    n_frg_a, n_frg_b : float | None
        Fragment count 'N' for each track, required by the fractional
        substrates. It counts the fragments 'compute_signal' admitted, which is
        what it divides a 'norm' track by. It is not the alignment-record
        count: the two agree only where filtering leaves one record per
        fragment. A fractional track cannot supply it, its column total being a
        deposition total rather than a count.
    total_a, total_b : float | None
        The substrate's own column total, required for 'unadj', where it is the
        total fragment base pairs. Only the track reports it: no report flag
        writes it, and it is not the fragment-bin overlap count 'L' that 'k'
        needs, which 'n_bin_a' carries.
    scale_a, scale_b : float | None
        Externally supplied deepTools scale factors, required for 'RPGC'.

    Returns
    -------
    result : dict[str, object]
        Keys 'scale_A', 'scale_B', 'pseudo_A', 'pseudo_B', 'prior_scaled_A',
        'prior_scaled_B', 'is_edger', and 'note'. The fractional substrates add
        'k_A' and 'k_B'; the whole-count and deepTools ones omit them, so a
        consumer tests for the key rather than assuming it.

    Raises
    ------
    ValueError
        For a nonpositive fragment-bin overlap count; a negative 'prior.count';
        an unknown substrate; a fractional substrate without both fragment
        counts; 'unadj' whose column total is absent, nonpositive, or below the
        fragment count; 'RPKM' without a positive 'siz_bin'; or 'RPGC' without
        both scale factors.

    Notes
    -----
    Every substrate gets a scale factor and a pseudocount: add the pseudocount
    to a track, then multiply by the scale. The pair differs by substrate
    because a prior must be denominated in the units it is added to.

    edgeR's 'cpm(log=TRUE)' computes '(y + y0_i) / (L_i + 2 y0_i) * 1e6' with
    'y0_i = pc * L_i / L_bar', writing
    - 'y' for one bin's count in track 'i',
    - 'L_i' for that track's fragment-bin overlap count,
    - 'y0_i' for its scaled prior,
    - 'pc' for 'prior.count', and
    - 'L_bar' for the mean of 'L_i' over the tracks.

    edgeR's 'addPriorCount' computes the same 'y0_i' without the '1e6'. The
    expression is linear in 'y', so it splits into a slope and an intercept,
    which is what a scale factor and a pseudocount are. 'is_edger' reports
    whether the returned pair reproduces it.

    Below, 'B' is the bin width, 'T' the substrate's own column total, 's' the
    scale factor returned alongside the pseudocount, and 'p_nc' the fractional
    prior derived below. This project's substrates return the following:
    - 'unadj', 'frag', 'norm': scale 1, pseudocount 'p_nc * T'.
    - 'count': scale 1, pseudocount 'y0_i'.
    - 'cpm': scale '1 / (1 + 2 pc / L_bar)', pseudocount 'pc * 1e6 / L_bar'.

    The deepTools substrates, which 'compute_pseudo_deeptools' serves:
    - 'CPM', 'BPM': scale '1e6 / (L_i + 2 y0_i)', pseudocount 's * y0_i'.
    - 'RPKM': scale '1e9 / ((L_i + 2 y0_i) B)', pseudocount 's * y0_i'.
    - 'None': scale 1, pseudocount 'y0_i'; this is 'count' reached elsewhere.
    - 'RPGC': scale supplied, pseudocount 's * y0_i'.

    So 'is_edger' is True for 'cpm', 'CPM', 'BPM' and 'RPKM'. Elsewhere the
    pair departs deliberately in three ways:
    1. The fractional substrates divide the prior by 'k_bar' to answer the
       under-dispersion of fractional deposition, giving
       'p_nc = pc / (k_bar * N_bar)' with 'k_i = L_i / N_i', averaged over the
       pair. Its 'pc' and '1 / N_bar' come from edgeR, but its '1 / k_bar' does
       not.
    2. 'count' and 'None' carry edgeR's prior without its denominator
       adjustment, so every bin's log2 ratio is offset by exactly
       'log2(L_B / L_A)'. That is the depth imbalance, which is what
       normalizing removes: relative comparison between bins is unaffected, but
       absolute fold change is not.
    3. 'RPGC' has no edgeR analogue for its 'N * F / G' denominator, so only
       the prior's magnitude crosses over.

    A closed substrate's pseudocount is symmetric in real arithmetic but not
    always in 'float64': 'pseudo_A' and 'pseudo_B' can land a bit apart, and
    which substrates do depends on the counts, so check symmetry with a
    tolerance. 'cpm' is the exception: it's symmetric bit for bit because its
    scale is computed once from the closed form rather than per sample.

    Single-track mode passes one track's overlap count as both 'n_bin_a' and
    'n_bin_b'. That is exact, not approximate: edgeR scales each prior by
    'L_i / L_bar' ('add_prior_count.c:88'), which is 1 with one column. But
    'L_bar' is that track's own count in single-track mode and the mean of both
    in two-track mode, so a track's single-track pseudocount is not its
    two-track one.

    Below, 'prior_scaled_A' and 'prior_scaled_B' expose 'y0_i' for reading;
    nothing consumes them. Their sum is always '2 * prior.count'; their ratio
    is 'L_A / L_B', or 'N_A / N_B' under the fractional substrates, which scale
    by 'N_i / N_bar' instead.
    """

    for label, lib in (("n_bin_a", n_bin_a), ("n_bin_b", n_bin_b)):
        if lib is None or not math.isfinite(lib) or lib <= 0.0:
            raise ValueError(
                f"{label!r} must be finite and positive; every substrate "
                "needs both fragment-bin overlap counts.",
            )

    if not math.isfinite(prior_count) or prior_count < 0.0:
        raise ValueError("'prior_count' must be finite and nonnegative.")

    substrate = canonicalize_substrate(substrate)

    if substrate in ("norm", "frag", "unadj"):
        for label, frg in (("n_frg_a", n_frg_a), ("n_frg_b", n_frg_b)):
            if frg is None or not math.isfinite(frg) or frg <= 0.0:
                raise ValueError(
                    f"{label!r} must be finite and positive for {substrate!r};"
                    " it is the fragment count, which a fractional track "
                    "cannot supply because its own total is not a count.",
                )

        k_a = n_bin_a / n_frg_a
        k_b = n_bin_b / n_frg_b
        k_mean = 0.5 * (k_a + k_b)
        n_frg_mean = 0.5 * (n_frg_a + n_frg_b)
        pseudo = prior_count / (k_mean * n_frg_mean)

        # Family members differ only by their column total 'T', so the prior
        # scales with it: 'p_i = p_nc * T_i'. 'T' is 1 for 'norm', the fragment
        # count for 'frag', and the total fragment base pairs for 'unadj'.
        if substrate == "frag":
            pseudo_a, pseudo_b = pseudo * n_frg_a, pseudo * n_frg_b
        elif substrate == "unadj":
            totals = (
                ("total_a", total_a, n_frg_a),
                ("total_b", total_b, n_frg_b),
            )

            for label, total, frg in totals:
                if total is None or not math.isfinite(total) or total <= 0.0:
                    raise ValueError(
                        f"{label!r} must be finite and positive for 'unadj'; "
                        "it is the track's own column sum, the total fragment "
                        "base pairs, which no report flag writes.",
                    )

                # Every fragment spans at least one base pair, so 'T >= N'
                # holds. A total below 'N' came off another substrate and would
                # rescale the prior in silence. It fires even when 'n_bin' is
                # supplied, which the inferred-'L' check in 'compute_pseudo'
                # cannot.
                if total < frg:
                    raise ValueError(
                        f"{label!r} is {total}, below the fragment count "
                        f"{frg}; an 'unadj' track sums to the total fragment "
                        "base pairs and every fragment spans at least one, so "
                        "this total cannot have come from one.",
                    )

            pseudo_a, pseudo_b = pseudo * total_a, pseudo * total_b
        else:
            pseudo_a = pseudo_b = pseudo

        return {
            "scale_A": 1.0,
            "scale_B": 1.0,
            "pseudo_A": pseudo_a,
            "pseudo_B": pseudo_b,
            "prior_scaled_A": prior_count * n_frg_a / n_frg_mean,
            "prior_scaled_B": prior_count * n_frg_b / n_frg_mean,
            "k_A": k_a,
            "k_B": k_b,
            "is_edger": False,
            "note": (
                "edgeR's prior divided by k_bar for the under-dispersion of "
                "fractional deposition, a correction edgeR does not make, "
                f"then denominated in {substrate!r}'s own column total"
            ),
        }

    n_bin_mean = 0.5 * (n_bin_a + n_bin_b)
    prior_a = prior_count * n_bin_a / n_bin_mean
    prior_b = prior_count * n_bin_b / n_bin_mean

    if substrate in ("CPM", "BPM"):
        scale_a = 1e6 / (n_bin_a + 2.0 * prior_a)
        scale_b = 1e6 / (n_bin_b + 2.0 * prior_b)
        is_edger = True
        note = "exact; BPM reduces to CPM with fixed bin width"
    elif substrate == "RPKM":
        if siz_bin is None or siz_bin <= 0:
            raise ValueError(
                "'siz_bin' must be given and positive for 'RPKM'; it is the "
                "bin width the per-kilobase denominator divides by.",
            )

        scale_a = 1e9 / ((n_bin_a + 2.0 * prior_a) * siz_bin)
        scale_b = 1e9 / ((n_bin_b + 2.0 * prior_b) * siz_bin)
        is_edger = True
        note = "exact"
    elif substrate in ("None", "count"):
        scale_a = 1.0
        scale_b = 1.0
        is_edger = False
        note = (
            "reproduces edgeR's ratio up to a constant log offset, since the "
            "denominator adjustment is absent"
        )
    elif substrate == "cpm":
        # Our 'cpm' track divides by plain 'L_i' where edgeR divides by
        # 'L_i + 2 y0_i', so the scale is their ratio. Since 'y0_i / L_i' is
        # 'pc / L_bar' for every track, that ratio is '1 / (1 + 2 pc / L_bar)',
        # with no sample index. Computing it once from that closed form keeps
        # the pair bit-identical; the per-sample quotient agrees only in real
        # arithmetic, but not always in 'float64'.
        scale_a = 1.0 / (1.0 + 2.0 * prior_count / n_bin_mean)
        scale_b = scale_a
        is_edger = True
        note = (
            "exact; the scale factor restores edgeR's adjusted denominator, "
            "and being symmetric it cancels in any A-over-B ratio"
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

    if substrate == "cpm":
        # The pseudocount 'y0_i * 1e6 / L_i' reduces to 'pc * 1e6 / L_bar', so
        # it too is symmetric.
        pseudo_a = prior_count * 1e6 / n_bin_mean
        pseudo_b = pseudo_a
    else:
        pseudo_a = scale_a * prior_a
        pseudo_b = scale_b * prior_b

    return {
        "scale_A": scale_a,
        "scale_B": scale_b,
        "pseudo_A": pseudo_a,
        "pseudo_B": pseudo_b,
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
    ValueError
        For an unrecognized 'method' or 'qntl_rule', or a 'qntl_pct' that is
        nonfinite or outside [0, 100].
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
