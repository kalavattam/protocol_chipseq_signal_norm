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
    SUBSTRATE_CANON,
    canonicalize_substrate,
    compute_pseudo_edger,
    compute_stats_robust,
    determine_coef_eff,
    iter_vals_bdg,
    median_sorted,
    pick_stabilizer,
)

# Verified against edgeR 4.4.0 to 1.589e-12 over 4,000 bins, 4.4e-16 in
# process. The pair's scale factors and pseudocount were re-derived under
# edgeR 4.8.2 on 2026-09-19 and are unchanged to all ten pinned digits. The
# same pair makes the 'float64' asymmetry visible.
OVLP_A = 24495.0
OVLP_B = 13605.0

# Counts 'c(5, 3, 0, 8)' sum to 16 and 'c(4, 2, 1, 6)' to 13, per edgeR 4.4.0
# and again under 4.8.2. Single-track and two-track values differ because the
# mean fragment-bin overlap count does.
OVLP_ONE_TRACK = 16.0
OVLP_PARTNER = 13.0

# Fragment counts for the fractional substrates. The 3:2 ratio is not the
# tracks' ratio, so a test cannot pass by reading fragment-bin overlap counts.
FRG_A = 1200.0
FRG_B = 800.0

# Column totals for 'unadj', the only substrate that needs them. Both clear
# their fragment count, as the 'T >= N' guard requires, and their 6:5 ratio is
# a third distinct ratio, so a test cannot pass by reading 'FRG_*' or 'OVLP_*'.
TOTAL_A = 240000.0
TOTAL_B = 200000.0

SUBSTRATE_FRACTIONAL = ("unadj", "frag", "norm")
SUBSTRATE_SYMMETRIC = ("CPM", "BPM", "RPKM", "norm", "cpm")
SUBSTRATE_ASYMMETRIC = ("None", "RPGC", "unadj", "frag", "count")
SUBSTRATE_ALL = SUBSTRATE_SYMMETRIC + SUBSTRATE_ASYMMETRIC


def _edger(substrate: str, **kwargs: float) -> dict[str, object]:
    """
    Call the estimator with the extra inputs each substrate requires.

    'RPKM' needs a bin width, 'RPGC' needs both deepTools scale factors, every
    fractional substrate needs both fragment counts, and 'unadj' needs the
    column totals on top of those, so a table-driven test cannot pass one
    argument set to every branch. 'count' and 'cpm' need nothing extra.
    """

    extra: dict[str, float] = {
        "n_ovlp_a": OVLP_A,
        "n_ovlp_b": OVLP_B,
        "substrate": substrate,
    }

    if substrate == "RPKM":
        extra["siz_bin"] = 10
    elif substrate == "RPGC":
        extra.update(scale_a=0.7, scale_b=1.3)
    elif substrate in SUBSTRATE_FRACTIONAL:
        extra.update(n_frg_a=FRG_A, n_frg_b=FRG_B)

        if substrate == "unadj":
            extra.update(total_a=TOTAL_A, total_b=TOTAL_B)

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
    routes, so exact '==' holds only for some overlap-count pairs.

    If this fails, 'float64' has started agreeing on this pair and the
    tolerance has lost its witness. Pick sizes that reproduce the difference;
    do not restore '=='.
    """

    result = _edger("CPM")

    assert result["pseudo_A"] != result["pseudo_B"]
    assert math.isclose(result["pseudo_A"], result["pseudo_B"], rel_tol=1e-12)


@pytest.mark.parametrize("substrate", SUBSTRATE_SYMMETRIC)
def test_compute_pseudo_edger_is_symmetric_where_the_substrate_closes(
    substrate: str,
) -> None:
    """
    Symmetry follows from closure, not from reproducing edgeR.

    A closed substrate divides by its own column total, so the per-sample
    scaling cancels and both tracks take the same pseudocount. 'norm' closes
    without reproducing edgeR and 'cpm' does both, so the two properties are
    tested apart.
    """

    result = _edger(substrate)

    assert math.isclose(result["pseudo_A"], result["pseudo_B"], rel_tol=1e-12)


@pytest.mark.parametrize("substrate", SUBSTRATE_ASYMMETRIC)
def test_compute_pseudo_edger_tracks_depth_without_the_adjustment(
    substrate: str,
) -> None:
    """
    An unclosed substrate keeps the depth difference in its pseudocount.

    It applies the prior's magnitude only: none of these adjusts the
    denominator to 'L_i + 2 * y0_i', the term that cancels the per-sample
    scaling everywhere else. 'None' and 'RPGC' leave it out; 'unadj', 'frag'
    and 'count' have a column total that carries the sample index in the first
    place.
    """

    result = _edger(substrate)

    assert result["pseudo_A"] != result["pseudo_B"]


@pytest.mark.parametrize("substrate", SUBSTRATE_ALL)
def test_compute_pseudo_edger_scales_the_prior_by_relative_depth(
    substrate: str,
) -> None:
    """
    The per-sample prior is 'y0_i', so its ratio is the depth ratio.

    Under a fractional substrate the depth that matters is the fragment count
    rather than the bin-sum, which is why the expected ratio switches with the
    branch.
    """

    result = _edger(substrate)
    prior_a = result["prior_scaled_A"]
    prior_b = result["prior_scaled_B"]

    if substrate in SUBSTRATE_FRACTIONAL:
        expected = FRG_A / FRG_B
    else:
        expected = OVLP_A / OVLP_B

    assert math.isclose(prior_a / prior_b, expected, rel_tol=1e-12)
    assert math.isclose(prior_a + prior_b, 2 * 2.0, rel_tol=1e-12)


@pytest.mark.parametrize(
    "substrate", ("CPM", "BPM", "RPKM", "None", "RPGC", "count")
)
def test_compute_pseudo_edger_factors_pseudo_into_scale_and_prior(
    substrate: str,
) -> None:
    result = _edger(substrate)

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


@pytest.mark.parametrize("substrate", SUBSTRATE_FRACTIONAL)
def test_compute_pseudo_edger_does_not_factor_pseudo_for_fractional(
    substrate: str,
) -> None:
    """
    Pin the branches where the decomposition does not hold.

    A fractional substrate derives its pseudocount from the fragment counts
    instead of from 'scale_i * prior_scaled_i', so a consumer cannot recover
    the per-sample prior by dividing. That is why the CLI emits it rather than
    leaving it to be derived.
    """

    result = _edger(substrate)

    assert result["scale_A"] == 1.0
    assert result["scale_B"] == 1.0
    assert result["pseudo_A"] != result["prior_scaled_A"]


def test_compute_pseudo_edger_does_not_factor_pseudo_for_cpm() -> None:
    """
    'cpm' also refuses the decomposition, but for the opposite reason.

    Its scale is not 1: it carries the correction that restores edgeR's
    adjusted denominator, and its pseudocount is computed from the closed form
    rather than from 'scale_i * prior_scaled_i'.
    """

    result = _edger("cpm")

    assert result["scale_A"] != 1.0
    assert result["pseudo_A"] != result["scale_A"] * result["prior_scaled_A"]


def test_compute_pseudo_edger_reproduces_edger_for_one_track() -> None:
    """
    A single track passes its overlap count twice rather than estimating.

    edgeR averages fragment-bin overlap counts over all columns and scales each
    prior by 'L_i / L_bar', so with one column that ratio is exactly 1 and the
    prior stays nominal. Both constants were checked against edgeR 4.4.0, and
    again against 4.8.2 on 2026-09-19.
    """

    one = compute_pseudo_edger(
        substrate="CPM", n_ovlp_a=OVLP_ONE_TRACK, n_ovlp_b=OVLP_ONE_TRACK
    )
    two = compute_pseudo_edger(
        substrate="CPM", n_ovlp_a=OVLP_ONE_TRACK, n_ovlp_b=OVLP_PARTNER
    )

    assert one["pseudo_A"] == 100000.0
    assert one["prior_scaled_A"] == 2.0
    assert two["pseudo_A"] == pytest.approx(108108.1081081081, abs=1e-9)
    assert one["pseudo_A"] != two["pseudo_A"]


@pytest.mark.parametrize(
    ("substrate", "is_edger"),
    [
        ("CPM", True),
        ("BPM", True),
        ("RPKM", True),
        ("None", False),
        ("RPGC", False),
        ("unadj", False),
        ("frag", False),
        ("norm", False),
        ("count", False),
        ("cpm", True),
    ],
)
def test_compute_pseudo_edger_reports_whether_it_reproduces_edger(
    substrate: str,
    is_edger: bool,
) -> None:
    """
    'is_edger' drives the CLI's stderr warning, so it is user-visible.
    """

    result = _edger(substrate)

    assert result["is_edger"] is is_edger
    assert result["note"]


@pytest.mark.parametrize("substrate", SUBSTRATE_FRACTIONAL)
def test_compute_pseudo_edger_returns_k_only_for_fractional_substrates(
    substrate: str,
) -> None:
    """
    Return 'k' only for the fractional substrates.

    'k_A' and 'k_B' are absent elsewhere, so a consumer must test for the key
    rather than assume it.
    """

    assert "k_A" not in _edger("CPM")
    assert "k_A" not in _edger("count")
    assert "k_A" not in _edger("cpm")
    assert _edger(substrate)["k_A"] == pytest.approx(OVLP_A / FRG_A, rel=1e-12)


# Each row is a call and the fragment its rejection must name. Short fragments
# keep a wording change from failing the test.
EDGER_REJECTIONS = (
    ({"n_ovlp_a": 0.0, "n_ovlp_b": OVLP_B, "substrate": "CPM"}, "n_ovlp_a"),
    ({"n_ovlp_a": -1.0, "n_ovlp_b": OVLP_B, "substrate": "CPM"}, "n_ovlp_a"),
    (
        {"n_ovlp_a": float("nan"), "n_ovlp_b": OVLP_B, "substrate": "CPM"},
        "n_ovlp_a",
    ),
    (
        {"n_ovlp_a": float("inf"), "n_ovlp_b": OVLP_B, "substrate": "CPM"},
        "n_ovlp_a",
    ),
    ({"n_ovlp_a": OVLP_A, "n_ovlp_b": 0.0, "substrate": "CPM"}, "n_ovlp_b"),
    (
        {
            "n_ovlp_a": OVLP_A,
            "n_ovlp_b": OVLP_B,
            "substrate": "CPM",
            "prior_count": -1.0,
        },
        "prior_count",
    ),
    (
        {
            "n_ovlp_a": OVLP_A,
            "n_ovlp_b": OVLP_B,
            "substrate": "RPKM",
            "siz_bin": 0,
        },
        "siz_bin",
    ),
    ({"n_ovlp_a": OVLP_A, "n_ovlp_b": OVLP_B, "substrate": "RPGC"}, "RPGC"),
    (
        {
            "n_ovlp_a": OVLP_A,
            "n_ovlp_b": OVLP_B,
            "substrate": "RPGC",
            "scale_a": 0.7,
        },
        "RPGC",
    ),
    ({"n_ovlp_a": OVLP_A, "n_ovlp_b": OVLP_B, "substrate": "norm"}, "n_frg_a"),
    (
        {
            "n_ovlp_a": OVLP_A,
            "n_ovlp_b": OVLP_B,
            "substrate": "norm",
            "n_frg_a": FRG_A,
        },
        "n_frg_b",
    ),
    (
        {
            "n_ovlp_a": OVLP_A,
            "n_ovlp_b": OVLP_B,
            "substrate": "norm",
            "n_frg_a": 0.0,
            "n_frg_b": FRG_B,
        },
        "n_frg_a",
    ),
    (
        {"n_ovlp_a": OVLP_A, "n_ovlp_b": OVLP_B, "substrate": "bogus"},
        "Unknown",
    ),
)


@pytest.mark.parametrize(("kwargs", "fragment"), EDGER_REJECTIONS)
def test_compute_pseudo_edger_rejects_unusable_inputs(
    kwargs: dict[str, float],
    fragment: str,
) -> None:
    with pytest.raises(ValueError, match=fragment):
        compute_pseudo_edger(**kwargs)


@pytest.mark.parametrize(
    ("alias", "canonical"), tuple(SUBSTRATE_CANON.items())
)
def test_canonicalize_substrate_maps_every_registered_alias(
    alias: str,
    canonical: str,
) -> None:
    """
    Both CLIs canonicalize '--substrate' through this mapping.

    Each restricts its own accepted subset, so a typo here changes what one or
    both command-line surfaces accept. Every key is covered rather than a
    representative few.
    """

    assert canonicalize_substrate(alias) == canonical


def test_frag_prior_is_the_norm_prior_times_the_fragment_count() -> None:
    """
    Check the fractional family's unclosed member against the closed form.

    'frag' is 'norm' before the depth division, so a 'frag' bin is an 'nc' bin
    times 'N_i'. The prior rides the same scalar, which is what keeps the two
    substrates' ratios a constant apart rather than differently shaped.
    """

    n_ovlp_a, n_ovlp_b = 6.0, 18.0
    n_frg_a, n_frg_b = 3.0, 6.0
    prior_count = 2.0

    closed = compute_pseudo_edger(
        n_ovlp_a=n_ovlp_a,
        n_ovlp_b=n_ovlp_b,
        prior_count=prior_count,
        substrate="norm",
        n_frg_a=n_frg_a,
        n_frg_b=n_frg_b,
    )
    frag = compute_pseudo_edger(
        n_ovlp_a=n_ovlp_a,
        n_ovlp_b=n_ovlp_b,
        prior_count=prior_count,
        substrate="frag",
        n_frg_a=n_frg_a,
        n_frg_b=n_frg_b,
    )

    assert frag["pseudo_A"] == closed["pseudo_A"] * n_frg_a
    assert frag["pseudo_B"] == closed["pseudo_B"] * n_frg_b
    assert frag["pseudo_A"] != frag["pseudo_B"]
    assert frag["scale_A"] == 1.0
    assert frag["is_edger"] is False


def test_frag_ratio_differs_from_nc_by_exactly_log2_of_depth_ratio() -> None:
    """
    Check the identity the per-sample prior exists to preserve.

    Scaling a track and its prior by the same per-sample factor multiplies the
    A-over-B ratio by 'N_A / N_B' and changes nothing else, so the two log2
    ratios differ by a constant. A prior that failed to ride the scalar would
    break this on the low-count bins first.
    """

    n_ovlp_a, n_ovlp_b = 6.0, 18.0
    n_frg_a, n_frg_b = 3.0, 6.0
    kwargs = {
        "n_ovlp_a": n_ovlp_a,
        "n_ovlp_b": n_ovlp_b,
        "prior_count": 2.0,
        "n_frg_a": n_frg_a,
        "n_frg_b": n_frg_b,
    }
    closed = compute_pseudo_edger(substrate="norm", **kwargs)
    frag = compute_pseudo_edger(substrate="frag", **kwargs)

    for nc_a, nc_b in ((0.0, 0.5), (0.25, 0.25), (1e-3, 7e-2)):
        ratio_nc = (nc_a + closed["pseudo_A"]) / (nc_b + closed["pseudo_B"])
        ratio_frag = (nc_a * n_frg_a + frag["pseudo_A"]) / (
            nc_b * n_frg_b + frag["pseudo_B"]
        )

        assert ratio_frag == pytest.approx(
            ratio_nc * n_frg_a / n_frg_b, rel=1e-12
        )


def test_unadj_prior_is_the_norm_prior_times_the_substrate_total() -> None:
    """
    Check the base-pair member against its own column total.

    An 'unadj' track sums to the total fragment base pairs, so its prior is the
    closed member's prior times that total. At one fixed fragment length the
    total is 'l * N', which reduces the prior to 'l' times the 'frag' prior;
    the general form carries no such assumption.
    """

    n_ovlp_a, n_ovlp_b = 6.0, 18.0
    n_frg_a, n_frg_b = 3.0, 6.0
    frg_len = 50.0
    kwargs = {
        "n_ovlp_a": n_ovlp_a,
        "n_ovlp_b": n_ovlp_b,
        "prior_count": 2.0,
        "n_frg_a": n_frg_a,
        "n_frg_b": n_frg_b,
    }
    closed = compute_pseudo_edger(substrate="norm", **kwargs)
    frag = compute_pseudo_edger(substrate="frag", **kwargs)
    unadj = compute_pseudo_edger(
        substrate="unadj",
        total_a=frg_len * n_frg_a,
        total_b=frg_len * n_frg_b,
        **kwargs,
    )

    assert unadj["pseudo_A"] == closed["pseudo_A"] * frg_len * n_frg_a
    assert unadj["pseudo_B"] == closed["pseudo_B"] * frg_len * n_frg_b

    # At a fixed fragment length 'unadj' is 'l' times 'frag', bin for bin.
    assert unadj["pseudo_A"] == pytest.approx(
        frag["pseudo_A"] * frg_len, rel=1e-12
    )
    assert unadj["is_edger"] is False


def test_unadj_requires_its_column_total() -> None:
    """
    Check that the one substrate needing a third quantity says so.
    """

    with pytest.raises(ValueError, match="total_a"):
        compute_pseudo_edger(
            n_ovlp_a=6.0,
            n_ovlp_b=18.0,
            prior_count=2.0,
            substrate="unadj",
            n_frg_a=3.0,
            n_frg_b=6.0,
        )


def test_unadj_refuses_a_total_that_cannot_be_base_pairs() -> None:
    """
    Refuse an 'unadj' column total that is impossible, not merely small.

    An 'unadj' track sums to the total fragment base pairs, and every fragment
    spans at least one base pair, so that total is never below 'N'. The same
    theorem guards the inferred fragment-bin overlap count one substrate over;
    here, it catches the case that guard cannot see, where '--n_ovlp' was
    supplied and the total still came off the wrong track. Measured, handing
    'unadj' a normalized track returns a prior 150 times too small in silence.
    """

    kwargs = {
        "n_ovlp_a": 6.0,
        "n_ovlp_b": 18.0,
        "prior_count": 2.0,
        "substrate": "unadj",
        "n_frg_a": 3.0,
        "n_frg_b": 6.0,
    }

    with pytest.raises(ValueError, match="total_a"):
        compute_pseudo_edger(total_a=1.0, total_b=150.0, **kwargs)

    with pytest.raises(ValueError, match="total_b"):
        compute_pseudo_edger(total_a=150.0, total_b=1.0, **kwargs)

    # The boundary a one-base fragment set reaches is 'total == n_frg', so it
    # must pass: the test is impossibility, not implausibility.
    result = compute_pseudo_edger(total_a=3.0, total_b=6.0, **kwargs)

    assert result["pseudo_A"] > 0.0


def test_cpm_reproduces_edger_exactly() -> None:
    """
    Pin the exactness claim the whole-count family exists to support.

    edgeR's 'cpm(log = TRUE)' inner value is
    '(y + y0_i) * 1e6 / (L_i + 2 * y0_i)' with
    'y0_i = prior.count * L_i / L_bar'. Our 'cpm' track is 'y * 1e6 / L_i', so
    recovering edgeR needs both the pseudocount and the scale factor this
    returns.
    """

    for n_ovlp_a, n_ovlp_b, prior_count in (
        (41318705.0, 39204118.0, 2.0),
        (6.0, 18.0, 2.0),
        (1e5, 9e7, 0.5),
    ):
        result = compute_pseudo_edger(
            n_ovlp_a=n_ovlp_a,
            n_ovlp_b=n_ovlp_b,
            prior_count=prior_count,
            substrate="cpm",
        )
        n_ovlp_mean = 0.5 * (n_ovlp_a + n_ovlp_b)

        for label, n_ovlp in (("A", n_ovlp_a), ("B", n_ovlp_b)):
            scale = result[f"scale_{label}"]
            pseudo = result[f"pseudo_{label}"]

            for count in (0.0, 1.0, 37.0, 1e4, 9.5e6):
                prior = prior_count * n_ovlp / n_ovlp_mean
                expected = (count + prior) * 1e6 / (n_ovlp + 2.0 * prior)
                actual = (count * 1e6 / n_ovlp + pseudo) * scale

                assert actual == pytest.approx(expected, rel=1e-12)

        assert result["is_edger"] is True


def test_cpm_scale_and_pseudo_are_bit_identical_across_samples() -> None:
    """
    Check the correction carries no sample index, to the last bit.

    'y0_i / L_i' is 'prior.count / L_bar' for every sample, so the adjustment
    is '1 / (1 + 2 * prior.count / L_bar)' and the closed member stays closed.
    Deriving it per sample instead would differ in the last bits, which is the
    same float64 trap 'compute_pseudo_edger' already documents for the prior.
    """

    result = compute_pseudo_edger(
        n_ovlp_a=41318705.0,
        n_ovlp_b=39204118.0,
        prior_count=2.0,
        substrate="cpm",
    )

    assert result["scale_A"] == result["scale_B"]
    assert result["pseudo_A"] == result["pseudo_B"]


def test_count_takes_the_prior_in_count_units() -> None:
    """
    Check the unclosed whole-count member carries a per-sample prior.

    A 'count' track is edgeR's own 'y', so its pseudocount is edgeR's own
    'y0_i' with no rescaling. It cannot reach edgeR's adjusted denominator in
    count units, which is why the closed sibling is the exact one.
    """

    n_ovlp_a, n_ovlp_b, prior_count = 6.0, 18.0, 2.0

    result = compute_pseudo_edger(
        n_ovlp_a=n_ovlp_a,
        n_ovlp_b=n_ovlp_b,
        prior_count=prior_count,
        substrate="count",
    )
    n_ovlp_mean = 0.5 * (n_ovlp_a + n_ovlp_b)

    assert result["scale_A"] == 1.0
    assert result["scale_B"] == 1.0
    assert result["pseudo_A"] == prior_count * n_ovlp_a / n_ovlp_mean
    assert result["pseudo_B"] == prior_count * n_ovlp_b / n_ovlp_mean
    assert result["pseudo_A"] != result["pseudo_B"]
    assert result["is_edger"] is False


def test_project_substrates_match_compute_signal_exactly() -> None:
    """
    Pin the claim 'SUBSTRATE_CANON' makes in its own comment.

    One substrate has one vocabulary across the codebase, so the project half
    of 'SUBSTRATE_CANON' must equal 'compute_signal.METHOD_CANON' key for key.
    The two maps drifted apart by hand twice during the split, once in each
    direction, so the comment is enforced here rather than trusted.
    """

    from protocol_chipseq_signal_norm.cli.compute_signal import METHOD_CANON

    deeptools = {"CPM", "BPM", "RPKM", "RPGC", "None"}
    project = {
        alias: canonical
        for alias, canonical in SUBSTRATE_CANON.items()
        if canonical not in deeptools
    }

    assert project == METHOD_CANON


def test_canonicalize_substrate_rejects_an_unregistered_name() -> None:
    # Not a spelling of anything: 'cpm' was this case until the whole-count
    # family registered it, and the retired aliases are the live regression.
    for unregistered in ("tpm", "n", "nrm", "normalized", "cnt", "ct", "cp"):
        with pytest.raises(ValueError, match="Unknown"):
            canonicalize_substrate(unregistered)


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
    ("sample", "n_ovlp_a", "n_ovlp_b", "n_frg_a", "n_frg_b", "expected"),
    NC_PUBLISHED,
    ids=[row[0] for row in NC_PUBLISHED],
)
def test_compute_pseudo_edger_reproduces_the_published_pseudocounts(
    sample: str,
    n_ovlp_a: float,
    n_ovlp_b: float,
    n_frg_a: float,
    n_frg_b: float,
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
        n_ovlp_a=n_ovlp_a,
        n_ovlp_b=n_ovlp_b,
        prior_count=2.0,
        substrate="norm",
        n_frg_a=n_frg_a,
        n_frg_b=n_frg_b,
    )

    assert math.isclose(result["pseudo_A"], expected, rel_tol=1e-15)
    assert result["pseudo_A"] == result["pseudo_B"]


# Restated so the reimplementation below never reads the module under test.
PRIOR_DEFAULT = 2.0

# 'BPM' is 'CPM' once the shared bin width cancels; 'RPKM' is 'CPM' per
# kilobase, and '_edger' passes a 10 bp bin.
SUBSTRATE_UNIT = (
    pytest.param("CPM", 1.0, id="CPM"),
    pytest.param("BPM", 1.0, id="BPM"),
    pytest.param("RPKM", 1e3 / 10.0, id="RPKM"),
)


@pytest.mark.parametrize(
    ("substrate", "unit"),
    SUBSTRATE_UNIT,
)
@pytest.mark.parametrize("count", (0.0, 1.0, 7.5, 1234.0))
def test_compute_pseudo_edger_decomposes_the_published_formula(
    substrate: str,
    unit: float,
    count: float,
) -> None:
    """
    's_i * y + p_i' must equal edgeR's '(y + y0_i) / (L_i + 2 * y0_i) * 1e6'.

    edgeR publishes one expression; this module returns a scale factor and a
    pseudocount, because deepTools takes them separately. Assert the identity
    at several counts rather than comparing the decomposition against itself,
    which would prove only that one expression equals itself.
    """

    result = _edger(substrate)
    n_ovlp_mean = 0.5 * (OVLP_A + OVLP_B)
    sides = (
        (OVLP_A, "scale_A", "pseudo_A"),
        (OVLP_B, "scale_B", "pseudo_B"),
    )

    for n_ovlp, scale_key, pseudo_key in sides:
        prior_scaled = PRIOR_DEFAULT * n_ovlp / n_ovlp_mean
        published = (
            (count + prior_scaled) / (n_ovlp + 2.0 * prior_scaled) * 1e6 * unit
        )
        decomposed = result[scale_key] * count + result[pseudo_key]

        assert math.isclose(decomposed, published, rel_tol=1e-12)


def test_compute_pseudo_edger_gives_bpm_the_cpm_result() -> None:
    """
    'BPM' reduces to 'CPM' with fixed bin width, so the fields must agree.

    BPM divides each bin by its width before summing; with fixed bin width that
    cancels, leaving CPM. The reduction follows from deepTools' '--binSize',
    not from a definition.

    Compare the whole result, so a future branch that diverges in the note or
    in 'is_edger' fails here too.
    """

    assert _edger("BPM") == _edger("CPM")


@pytest.mark.parametrize("siz_bin", (1, 10, 50, 200))
def test_compute_pseudo_edger_keeps_cpm_ratios_under_rpkm(
    siz_bin: int,
) -> None:
    """
    'RPKM' rescales both tracks by one bin-width constant, so ratios survive.

    The values legitimately differ, which is why this asserts the ratio. The
    third assertion names the constant, so a change that preserved the ratio
    while corrupting the scale would still fail.
    """

    cpm = _edger("CPM")
    rpkm = _edger("RPKM", siz_bin=siz_bin)
    per_kilobase = 1e3 / siz_bin

    assert math.isclose(
        rpkm["pseudo_A"] / rpkm["pseudo_B"],
        cpm["pseudo_A"] / cpm["pseudo_B"],
        rel_tol=1e-12,
    )
    assert math.isclose(
        rpkm["scale_A"] / rpkm["scale_B"],
        cpm["scale_A"] / cpm["scale_B"],
        rel_tol=1e-12,
    )
    assert math.isclose(
        rpkm["scale_A"], cpm["scale_A"] * per_kilobase, rel_tol=1e-12
    )
