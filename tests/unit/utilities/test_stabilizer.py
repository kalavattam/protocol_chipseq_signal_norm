#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_stabilizer.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.6);
# - Anthropic Claude Code (Opus 5).
#
# Distributed under the MIT license.


import math
from pathlib import Path

import pytest

from protocol_chipseq_signal_norm.utilities.utils_stabilizer import (
    NORM_CANON,
    canonicalize_norm,
    compute_pseudo_edger,
    compute_stats_robust,
    determine_coef_eff,
    iter_vals_bdg,
    median_sorted,
    pick_stabilizer,
)

# Verified against edgeR 4.4.0 to 1.589e-12 over 4,000 bins, 4.4e-16 in
# process. The same pair makes the 'float64' asymmetry visible.
LIB_A = 24495.0
LIB_B = 13605.0

# Counts 'c(5, 3, 0, 8)' sum to 16 and 'c(4, 2, 1, 6)' to 13, per edgeR 4.4.0.
# The two framings differ because the mean library size does.
LIB_ONE_TRACK = 16.0
LIB_PARTNER = 13.0

# Fragment counts for 'norm'. The 3:2 ratio is not the tracks' ratio, so a test
# cannot pass by reading library sizes.
FRG_A = 1200.0
FRG_B = 800.0

NORM_SYMMETRIC = ("CPM", "BPM", "RPKM", "norm")
NORM_ASYMMETRIC = ("None", "RPGC")
NORM_ALL = NORM_SYMMETRIC + NORM_ASYMMETRIC


def _edger(norm: str, **kwargs: float) -> dict[str, object]:
    """
    Call the estimator with the extra inputs each normalization requires.

    'RPKM' needs a bin width, 'RPGC' needs both deepTools scale factors, and
    'norm' needs both fragment counts, so a table-driven test cannot pass one
    argument set to every branch.
    """

    extra: dict[str, float] = {"lib_a": LIB_A, "lib_b": LIB_B, "norm": norm}

    if norm == "RPKM":
        extra["siz_bin"] = 10
    elif norm == "RPGC":
        extra.update(scale_a=0.7, scale_b=1.3)
    elif norm == "norm":
        extra.update(frg_a=FRG_A, frg_b=FRG_B)

    extra.update(kwargs)

    return compute_pseudo_edger(**extra)


def test_compute_stats_robust_ignores_nonfinite_values() -> None:
    stats = compute_stats_robust([1.0, float("nan"), 3.0, float("inf")])

    assert stats == {"n": 2, "median": 2.0, "mean": 2.0}


def test_determine_coef_eff_and_median_sorted() -> None:
    assert determine_coef_eff("frc_mdn_nz", None) == 0.01
    assert determine_coef_eff("min_nz", None) == 1.0
    assert determine_coef_eff("qntl_nz", None) is None
    assert median_sorted([1.0, 2.0, 10.0, 20.0]) == 6.0


def test_pick_stabilizer_quantile_and_floor() -> None:
    assert pick_stabilizer([1.0, 2.0, 10.0], "qntl_nz", qntl_pct=50) == 2.0
    assert pick_stabilizer([1.0, 2.0], "min_nz", floor=5.0) == 5.0
    assert math.isnan(pick_stabilizer([float("nan")], "qntl_nz"))


def test_pick_stabilizer_rejects_bad_quantile() -> None:
    with pytest.raises(ValueError, match="qntl_pct"):
        pick_stabilizer([1.0], "qntl_nz", qntl_pct=101)


def test_iter_vals_bdg_filters_by_positive_policy(tmp_path: Path) -> None:
    path = tmp_path / "values.bdg"
    path.write_text(
        "chrI 0 10 0\nchrI 10 20 0.1\nchrI 20 30 2\nchrI 30 40 nan\n",
        encoding="utf-8",
    )

    assert list(
        iter_vals_bdg(
            str(path),
            eps=0.1,
            mode_nz="closed",
            nz_policy="pos",
        ),
    ) == [2.0]


def test_compute_pseudo_edger_reproduces_the_edger_verified_pair() -> None:
    """
    Pin the decomposition edgeR itself was checked against.

    These are the values a deepTools run consumes as '--scaleFactors A:B' and
    '--pseudocount P P'. Drift here is drift away from edgeR, which no other
    test in this file would report.
    """

    result = _edger("CPM")

    assert result["scale_A"] == pytest.approx(40.8160877863, abs=1e-10)
    assert result["scale_B"] == pytest.approx(73.4869584951, abs=1e-10)
    assert result["pseudo_A"] == pytest.approx(104.9648367797, abs=1e-10)


def test_compute_pseudo_edger_symmetry_survives_only_a_tolerance() -> None:
    """
    Keep a live witness for why the symmetry tests use 'math.isclose'.

    The equality is exact in real arithmetic (the per-sample 'L_i / L_bar'
    scaling cancels) but not in 'float64', where the two sides take different
    routes, so exact '==' holds only for some library-size pairs.

    If this fails, 'float64' has started agreeing on this pair and the
    tolerance has lost its witness. Pick sizes that reproduce the difference;
    do not restore '=='.
    """

    result = _edger("CPM")

    assert result["pseudo_A"] != result["pseudo_B"]
    assert math.isclose(result["pseudo_A"], result["pseudo_B"], rel_tol=1e-12)


@pytest.mark.parametrize("norm", NORM_SYMMETRIC)
def test_compute_pseudo_edger_is_symmetric_where_edger_is(norm: str) -> None:
    result = _edger(norm)

    assert math.isclose(result["pseudo_A"], result["pseudo_B"], rel_tol=1e-12)


@pytest.mark.parametrize("norm", NORM_ASYMMETRIC)
def test_compute_pseudo_edger_tracks_depth_without_the_adjustment(
    norm: str,
) -> None:
    """
    'None' and 'RPGC' apply the prior's magnitude only, so depth survives.

    Neither adjusts the denominator to 'L_i + 2 * y0_i', which is the term that
    cancels the per-sample scaling everywhere else.
    """

    result = _edger(norm)

    assert result["pseudo_A"] != result["pseudo_B"]


@pytest.mark.parametrize("norm", NORM_ALL)
def test_compute_pseudo_edger_scales_the_prior_by_relative_depth(
    norm: str,
) -> None:
    """
    The per-sample prior is 'y0_i', so its ratio is the depth ratio.

    Under 'norm' the depth that matters is the fragment count rather than the
    bin-sum, which is why the expected ratio switches with the branch.
    """

    result = _edger(norm)
    prior_a = result["prior_scaled_A"]
    prior_b = result["prior_scaled_B"]
    expected = FRG_A / FRG_B if norm == "norm" else LIB_A / LIB_B

    assert math.isclose(prior_a / prior_b, expected, rel_tol=1e-12)
    assert math.isclose(prior_a + prior_b, 2 * 2.0, rel_tol=1e-12)


@pytest.mark.parametrize("norm", ("CPM", "BPM", "RPKM", "None", "RPGC"))
def test_compute_pseudo_edger_factors_pseudo_into_scale_and_prior(
    norm: str,
) -> None:
    result = _edger(norm)

    assert math.isclose(
        result["pseudo_A"],
        result["scale_A"] * result["prior_scaled_A"],
        rel_tol=1e-12,
    )
    assert math.isclose(
        result["pseudo_B"],
        result["scale_B"] * result["prior_scaled_B"],
        rel_tol=1e-12,
    )


def test_compute_pseudo_edger_does_not_factor_pseudo_for_normalized() -> None:
    """
    Pin the one branch where the decomposition does not hold.

    'norm' derives a symmetric pseudocount from the fragment counts instead of
    from 'scale_i * prior_scaled_i', so a consumer cannot recover the
    per-sample prior by dividing. That is why the CLI emits it rather than
    leaving it to be derived.
    """

    result = _edger("norm")

    assert result["scale_A"] == 1.0
    assert result["scale_B"] == 1.0
    assert result["pseudo_A"] != result["prior_scaled_A"]


def test_compute_pseudo_edger_reproduces_edger_for_one_library() -> None:
    """
    One track is expressed by passing its library size twice, not estimated.

    edgeR averages library sizes over all columns and scales each prior by
    'L_i / mean(L)', so with one column that ratio is exactly 1 and the prior
    stays nominal. Both constants were checked against edgeR 4.4.0.
    """

    one = compute_pseudo_edger(LIB_ONE_TRACK, LIB_ONE_TRACK, norm="CPM")
    two = compute_pseudo_edger(LIB_ONE_TRACK, LIB_PARTNER, norm="CPM")

    assert one["pseudo_A"] == 100000.0
    assert one["prior_scaled_A"] == 2.0
    assert two["pseudo_A"] == pytest.approx(108108.1081081081, abs=1e-9)
    assert one["pseudo_A"] != two["pseudo_A"]


@pytest.mark.parametrize(
    ("norm", "is_edger"),
    [
        ("CPM", True),
        ("BPM", True),
        ("RPKM", True),
        ("None", False),
        ("RPGC", False),
        ("norm", False),
    ],
)
def test_compute_pseudo_edger_reports_whether_it_reproduces_edger(
    norm: str,
    is_edger: bool,
) -> None:
    """
    'is_edger' drives the CLI's stderr warning, so it is user-visible.
    """

    result = _edger(norm)

    assert result["is_edger"] is is_edger
    assert result["note"]


def test_compute_pseudo_edger_returns_k_only_for_normalized_coverage() -> None:
    """
    'k_A' and 'k_B' are absent outside 'norm', so a consumer must test first.
    """

    assert "k_A" not in _edger("CPM")
    assert _edger("norm")["k_A"] == pytest.approx(LIB_A / FRG_A, rel=1e-12)


# Each row is a call and the fragment its rejection must name. Short fragments
# keep a wording change from failing the test.
EDGER_REJECTIONS = (
    ({"lib_a": 0.0, "lib_b": LIB_B}, "lib_a"),
    ({"lib_a": -1.0, "lib_b": LIB_B}, "lib_a"),
    ({"lib_a": float("nan"), "lib_b": LIB_B}, "lib_a"),
    ({"lib_a": float("inf"), "lib_b": LIB_B}, "lib_a"),
    ({"lib_a": LIB_A, "lib_b": 0.0}, "lib_b"),
    ({"lib_a": LIB_A, "lib_b": LIB_B, "prior_count": -1.0}, "prior_count"),
    (
        {"lib_a": LIB_A, "lib_b": LIB_B, "norm": "RPKM", "siz_bin": 0},
        "siz_bin",
    ),
    ({"lib_a": LIB_A, "lib_b": LIB_B, "norm": "RPGC"}, "RPGC"),
    (
        {"lib_a": LIB_A, "lib_b": LIB_B, "norm": "RPGC", "scale_a": 0.7},
        "RPGC",
    ),
    ({"lib_a": LIB_A, "lib_b": LIB_B, "norm": "norm"}, "frg_a"),
    (
        {"lib_a": LIB_A, "lib_b": LIB_B, "norm": "norm", "frg_a": FRG_A},
        "frg_b",
    ),
    (
        {
            "lib_a": LIB_A,
            "lib_b": LIB_B,
            "norm": "norm",
            "frg_a": 0.0,
            "frg_b": FRG_B,
        },
        "frg_a",
    ),
    ({"lib_a": LIB_A, "lib_b": LIB_B, "norm": "bogus"}, "Unknown"),
)


@pytest.mark.parametrize(("kwargs", "fragment"), EDGER_REJECTIONS)
def test_compute_pseudo_edger_rejects_unusable_inputs(
    kwargs: dict[str, float],
    fragment: str,
) -> None:
    with pytest.raises(ValueError, match=fragment):
        compute_pseudo_edger(**kwargs)


@pytest.mark.parametrize(("alias", "canonical"), tuple(NORM_CANON.items()))
def test_canonicalize_norm_maps_every_registered_alias(
    alias: str,
    canonical: str,
) -> None:
    """
    The CLI builds '--normalization' choices from this mapping.

    A typo here changes the accepted command-line surface, so every key is
    covered rather than a representative few.
    """

    assert canonicalize_norm(alias) == canonical


def test_canonicalize_norm_rejects_an_unregistered_name() -> None:
    with pytest.raises(ValueError, match="Unknown"):
        canonicalize_norm("cpm")


NC_PUBLISHED = (
    (
        "Hho1_6336",
        318446200.0,
        264987537.0,
        13492934.476769,
        12851814.515126,
        6.867215221172314e-09,
    ),
    (
        "Hho1_6337",
        327429285.0,
        251418845.0,
        13655960.053875,
        12029951.222678,
        6.940275297124279e-09,
    ),
    (
        "Hmo1_7750",
        336446562.0,
        162564489.0,
        14276247.374153,
        7475008.643619,
        8.116474095065805e-09,
    ),
)


@pytest.mark.parametrize(
    ("sample", "lib_a", "lib_b", "frg_a", "frg_b", "expected"),
    NC_PUBLISHED,
    ids=[row[0] for row in NC_PUBLISHED],
)
def test_compute_pseudo_edger_reproduces_the_published_pseudocounts(
    sample: str,
    lib_a: float,
    lib_b: float,
    frg_a: float,
    frg_b: float,
    expected: float,
) -> None:
    """
    Pin the normalized-coverage values a real run produced.

    These are the only assertions made from measured sequencing depths, so they
    catch a change that stays arithmetically self-consistent but no longer
    matches the data. The tolerance is tight because this pins one computation
    against its own recorded output: at 'rel_tol=1e-12' a change in the
    thirteenth significant digit still passed.
    """

    result = compute_pseudo_edger(
        lib_a=lib_a,
        lib_b=lib_b,
        prior_count=2.0,
        norm="norm",
        frg_a=frg_a,
        frg_b=frg_b,
    )

    assert math.isclose(result["pseudo_A"], expected, rel_tol=1e-15)
    assert result["pseudo_A"] == result["pseudo_B"]
