#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_ratio.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.6);
# - Anthropic Claude Code (Opus 5, Opus 5.5).
#
# Distributed under the MIT license.


import argparse
import math
from pathlib import Path

import pytest

from protocol_chipseq_signal_norm.cli.compute_signal_ratio import (
    METHOD_CANON,
    calc_rat_bin,
    comp_sig_rat,
    main,
    parse_args,
    parse_pair,
)
from protocol_chipseq_signal_norm.utilities.utils_cli import (
    CapArgumentParser,
)


def test_parse_pair_accepts_single_or_pair_values() -> None:
    assert parse_pair("2", 1.0) == (2.0, 1.0)
    assert parse_pair("2:3.5", 1.0) == (2.0, 3.5)

    with pytest.raises(Exception, match="Expected"):
        parse_pair("2:3:4", 1.0)


def test_calc_rat_bin_handles_scaling_and_log_reciprocal() -> None:
    scaled = calc_rat_bin(
        4.0,
        2.0,
        2.0,
        1.0,
        1.0,
        0.0,
        None,
        False,
        False,
        None,
    )
    reciprocal = calc_rat_bin(
        4.0,
        2.0,
        None,
        None,
        None,
        None,
        None,
        True,
        True,
        None,
    )

    assert scaled == 5.0
    assert reciprocal == -1.0


def test_calc_rat_bin_zero_handling() -> None:
    zero_result = calc_rat_bin(
        0,
        0,
        None,
        None,
        None,
        None,
        None,
        False,
        False,
        "pre_scale",
    )

    assert zero_result is None
    assert math.isnan(
        calc_rat_bin(1, 0, None, None, None, None, None, False, False, None),
    )


# Measured on the pre-reorder tree at '627d12a'. Scaling moved to after the
# pseudocount add, so every unscaled case has to stay put.
BASELINE_UNSCALED = [
    (
        (4.0, 2.0), (1.0, 1.0), (1.0, 1.0), None, False, False, None, 0.0,
        1.6666666666666667
    ),
    (
        (4.0, 2.0), (None, None), (0.5, 0.4), None, True, False, None, 0.0,
        0.9068905956085185
    ),
    (
        (4.0, 0.0), (1.0, 1.0), (1.0, 0.0), None, True, False, None, 0.0,
        "nan"
    ),
    (
        (4.0, 0.0), (None, None), (None, None), None, True, False, None, 0.0,
        "nan"
    ),
    (
        (0.0, 2.0), (1.0, 1.0), (2.0, 3.0), 2.0, False, False, None, 0.0, 0.4
    ),
    (
        (0.0, 0.0), (1.0, 1.0), (1.0, 1.0), None, True, False, "post_scale",
        0.0, None
    ),
    (
        (1e-09, 3.0), (None, None), (1.0, 1.0), None, True, False, "pre_scale",
        1e-12, -1.999999998557305
    ),
    (
        (7.5, 0.25), (1.0, 1.0), (0.5, 0.4), 0.5, True, True, None, 0.0,
        -3.62148837674627
    ),
]


@pytest.mark.parametrize(
    "sig, scl, psc, dep_min, log2, recip, skip_00, eps, want",
    BASELINE_UNSCALED,
)
def test_unscaled_ratios_match_the_pre_reorder_baseline(
    sig: tuple[float, float],
    scl: tuple[float | None, float | None],
    psc: tuple[float | None, float | None],
    dep_min: float | None,
    log2: bool,
    recip: bool,
    skip_00: str | None,
    eps: float,
    want: float | str | None,
) -> None:
    got = calc_rat_bin(
        sig[0], sig[1], scl[0], scl[1], psc[0], psc[1], dep_min, log2, recip,
        skip_00, eps,
    )

    if want == "nan":
        assert math.isnan(got)
    else:
        assert got == want


def test_scaled_log2_ratio_is_an_additive_offset_on_the_unscaled_one() -> None:
    sig_a, sig_b, psc_a, psc_b = 4.0, 2.0, 1.0, 0.5

    plain = calc_rat_bin(
        sig_a, sig_b, None, None, psc_a, psc_b, None, True, False, None,
    )

    # Powers of two over a ratio that divides exactly, so binary64 loses
    # nothing and the identity pins exactly rather than to a tolerance.
    exact = calc_rat_bin(
        sig_a, sig_b, 2.0, 4.0, psc_a, psc_b, None, True, False, None,
    )

    assert exact == math.log2(2.0 / 4.0) + plain

    # Real spike-in coefficients are not powers of two, so the two groupings
    # differ in the last bits. 1e-15 covers that and nothing coarser.
    rough = calc_rat_bin(
        7.5, 0.25, 1.1413, 1.0, 1.0, 0.4, None, True, False, None,
    )
    rough_plain = calc_rat_bin(
        7.5, 0.25, None, None, 1.0, 0.4, None, True, False, None,
    )

    assert rough == pytest.approx(math.log2(1.1413) + rough_plain, rel=1e-15)


def test_post_scale_zero_test_reads_the_scaled_pair() -> None:
    # Scaling moved after the pseudocount, so 'post_scale' computes the scaled
    # pair on the side. A tiny scale drops this bin; the raw pair would not.
    dropped = calc_rat_bin(
        4.0, 2.0, 1e-4, 1e-4, None, None, None, False, False, "post_scale",
        1e-2,
    )
    kept = calc_rat_bin(
        4.0, 2.0, None, None, None, None, None, False, False, "post_scale",
        1e-2,
    )

    assert dropped is None
    assert kept == 2.0


def test_dep_min_clamps_the_scaled_regularized_denominator() -> None:
    # 'dep_min' still sits after both steps, so under scaling it now sees
    # 'scl_b * (sig_b + psc_b)' rather than 'scl_b * sig_b + psc_b'.
    clamped = calc_rat_bin(
        4.0, 1.0, 2.0, 3.0, 1.0, 0.5, 10.0, False, False, None,
    )
    unclamped = calc_rat_bin(
        4.0, 1.0, 2.0, 3.0, 1.0, 0.5, 4.0, False, False, None,
    )

    assert clamped == (2.0 * (4.0 + 1.0)) / 10.0
    assert unclamped == (2.0 * (4.0 + 1.0)) / (3.0 * (1.0 + 0.5))


def test_each_pseudocount_scales_with_the_side_it_regularizes() -> None:
    got = calc_rat_bin(4.0, 2.0, 2.0, 3.0, 1.0, 0.5, None, False, False, None)

    assert got == (2.0 * (4.0 + 1.0)) / (3.0 * (2.0 + 0.5))


def test_scaling_reaches_a_bin_whose_numerator_is_only_pseudocount() -> None:
    # 'A' is empty, so before the reorder the numerator was the bare
    # pseudocount and a siQ-ChIP coefficient never reached the bin.
    siq = calc_rat_bin(
        0.0, 2.0, 0.0146, 1.0, 1.0, 1.0, None, True, False, None,
    )
    unity = calc_rat_bin(0.0, 2.0, 1.0, 1.0, 1.0, 1.0, None, True, False, None)

    assert siq != unity
    assert siq == pytest.approx(math.log2(0.0146) + unity, rel=1e-15)


def test_comp_sig_rat_writes_small_bedgraph(tmp_path: Path) -> None:
    fil_a = tmp_path / "a.bdg"
    fil_b = tmp_path / "b.bdg"
    fil_out = tmp_path / "ratio.bdg"
    fil_a.write_text("chrI 0 10 4\nchrI 10 20 0\n", encoding="utf-8")
    fil_b.write_text("chrI 0 10 2\nchrI 10 20 0\n", encoding="utf-8")

    comp_sig_rat(
        str(fil_a),
        str(fil_b),
        str(fil_out),
        scl_a=1.0,
        scl_b=1.0,
        psc_a=0.0,
        psc_b=0.0,
        dep_min=None,
        decimal_places=3,
        log2=False,
        recip=False,
        skip_00="pre_scale",
        eps=0.0,
        track=False,
        drp_nan=False,
    )

    assert fil_out.read_text(encoding="utf-8") == "chrI\t0\t10\t2\n"


def test_comp_sig_rat_strict_bins_rejects_end_mismatch(tmp_path: Path) -> None:
    fil_a = tmp_path / "a.bdg"
    fil_b = tmp_path / "b.bdg"
    fil_out = tmp_path / "ratio.bdg"
    fil_a.write_text("chrI 0 10 4\nchrI 10 20 2\n", encoding="utf-8")
    fil_b.write_text("chrI 0 10 2\nchrI 10 21 1\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Mismatched bedGraph grids"):
        comp_sig_rat(
            str(fil_a),
            str(fil_b),
            str(fil_out),
            scl_a=1.0,
            scl_b=1.0,
            psc_a=0.0,
            psc_b=0.0,
            dep_min=None,
            decimal_places=3,
            log2=False,
            recip=False,
            skip_00=None,
            eps=0.0,
            track=False,
            drp_nan=False,
            strict_bins=True,
        )


def test_comp_sig_rat_chr_sizes_rejects_out_of_bounds_input(
    tmp_path: Path,
) -> None:
    fil_a = tmp_path / "a.bdg"
    fil_b = tmp_path / "b.bdg"
    fil_out = tmp_path / "ratio.bdg"
    fil_a.write_text("chrI 70 81 4\n", encoding="utf-8")
    fil_b.write_text("chrI 70 80 2\n", encoding="utf-8")

    with pytest.raises(ValueError, match="extends beyond"):
        comp_sig_rat(
            str(fil_a),
            str(fil_b),
            str(fil_out),
            scl_a=1.0,
            scl_b=1.0,
            psc_a=0.0,
            psc_b=0.0,
            dep_min=None,
            decimal_places=3,
            log2=False,
            recip=False,
            skip_00=None,
            eps=0.0,
            track=False,
            drp_nan=False,
            chrom_sizes={"chrI": 80},
        )


def test_comp_sig_rat_rejects_dash_io(tmp_path: Path) -> None:
    fil_a = tmp_path / "a.bdg"
    fil_b = tmp_path / "b.bdg"
    fil_out = tmp_path / "ratio.bdg"
    fil_a.write_text("chrI 0 10 4\n", encoding="utf-8")
    fil_b.write_text("chrI 0 10 2\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Dash input/output"):
        comp_sig_rat(
            "-",
            str(fil_b),
            str(fil_out),
            scl_a=1.0,
            scl_b=1.0,
            psc_a=0.0,
            psc_b=0.0,
            dep_min=None,
            decimal_places=3,
            log2=False,
            recip=False,
            skip_00=None,
            eps=0.0,
            track=False,
            drp_nan=False,
        )


# Each hidden hyphen alias must restate its primary's 'type', 'choices', and
# action; every pair is tested, and a guard fails on an alias missing here.
HYPHEN_ALIASES = (
    pytest.param("--chr_siz", "--chr-siz", "value", id="chr_siz"),
    pytest.param("--dep_min", "--dep-min", "0.5", id="dep_min"),
    pytest.param("--drp_nan", "--drp-nan", None, id="drp_nan"),
    pytest.param("--fil_A", "--fil-A", "value", id="fil_A"),
    pytest.param("--fil_B", "--fil-B", "value", id="fil_B"),
    pytest.param("--fil_out", "--fil-out", "value", id="fil_out"),
    pytest.param("--scl_fct", "--scl-fct", "value", id="scl_fct"),
    pytest.param("--skip_00", "--skip-00", "pre_scale", id="skip_00"),
    pytest.param("--skp_pfx", "--skp-pfx", "value", id="skp_pfx"),
    pytest.param("--strict_bins", "--strict-bins", None, id="strict_bins"),
)


def _alias_parser() -> argparse.ArgumentParser:
    """
    Return the parser built by the module under test.
    """

    captured: dict[str, argparse.ArgumentParser] = {}
    original = CapArgumentParser.parse_args

    def capture(
        self: CapArgumentParser,
        *args: object,
        **kwargs: object,
    ) -> None:
        captured["parser"] = self

        raise SystemExit(0)

    CapArgumentParser.parse_args = capture

    try:
        parse_args(["-fA", "a.bdg", "-fB", "b.bdg", "-fo", "c.bdg"])
    except SystemExit:
        pass
    finally:
        CapArgumentParser.parse_args = original

    return captured["parser"]


def _registered_hyphen_aliases() -> set[str]:
    """
    Return hidden hyphen spellings that alias a visible option.

    A hidden *option* may also carry a hyphen spelling; those are excluded, as
    they alias nothing the user can see.
    """

    parser = _alias_parser()
    visible = {
        action.dest
        for action in parser._actions
        if action.help is not argparse.SUPPRESS
    }

    return {
        option
        for action in parser._actions
        if action.help is argparse.SUPPRESS and action.dest in visible
        for option in action.option_strings
        if option.startswith("--") and "_" not in option
    }


@pytest.mark.parametrize(("primary", "alias", "value"), HYPHEN_ALIASES)
def test_hidden_hyphen_alias_matches_its_primary(
    primary: str,
    alias: str,
    value: str | None,
) -> None:
    """
    Each hidden hyphen alias parses to the primary's value and type.
    """

    base = list(["-fA", "a.bdg", "-fB", "b.bdg", "-fo", "c.bdg"])
    supplied = [primary] if value is None else [primary, value]
    aliased = [alias] if value is None else [alias, value]
    destination = primary.lstrip("-")

    from_primary = getattr(parse_args(base + supplied), destination)
    from_alias = getattr(parse_args(base + aliased), destination)

    assert from_primary == from_alias
    assert type(from_primary) is type(from_alias)


def test_every_hidden_hyphen_alias_is_covered() -> None:
    """
    No hidden hyphen alias may exist without a row in 'HYPHEN_ALIASES'.
    """

    covered = {row.values[1] for row in HYPHEN_ALIASES}
    registered = _registered_hyphen_aliases()

    assert registered == covered


def test_hidden_hyphen_aliases_stay_out_of_rendered_help(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """
    Hidden aliases remain usable but are never advertised in help.
    """

    with pytest.raises(SystemExit):
        parse_args(["--help"])

    rendered = capsys.readouterr().out
    registered = _registered_hyphen_aliases()

    for alias in registered:
        assert alias not in rendered


def test_hidden_aliases_satisfy_the_required_ratio_options() -> None:
    """
    A run supplying only hyphen spellings is accepted.

    'argparse' enforces 'required' per action, so a hidden alias cannot satisfy
    a required primary. The requirement is enforced on the parsed values in
    'main' instead, and this pins that the alias route works.
    """

    args = parse_args(
        ["--fil-A", "a.bdg", "--fil-B", "b.bdg", "--fil-out", "c.bdg"],
    )

    assert args.fil_A == "a.bdg"
    assert args.fil_B == "b.bdg"
    assert args.fil_out == "c.bdg"


@pytest.mark.parametrize(
    ("omitted", "supplied"),
    [
        pytest.param(
            "--fil_A",
            ["--fil_B", "b.bdg", "--fil_out", "c.bdg"],
            id="fil_A",
        ),
        pytest.param(
            "--fil_B",
            ["--fil_A", "a.bdg", "--fil_out", "c.bdg"],
            id="fil_B",
        ),
        pytest.param(
            "--fil_out",
            ["--fil_A", "a.bdg", "--fil_B", "b.bdg"],
            id="fil_out",
        ),
    ],
)
def test_missing_required_ratio_option_is_rejected(
    omitted: str,
    supplied: list[str],
) -> None:
    """
    Dropping the requirement from argparse must not drop the requirement.
    """

    with pytest.raises(SystemExit) as error:
        main(supplied)

    assert omitted in str(error.value)


@pytest.mark.parametrize("pseudo", ["1", "0.5", "0"])
def test_a_bare_pseudocount_is_refused(tmp_path: Path, pseudo: str) -> None:
    """
    A bare 'A' would regularize file A alone, so it is refused, naming 'A:B'.
    """

    fil_a = tmp_path / "a.bdg"
    fil_b = tmp_path / "b.bdg"
    fil_out = tmp_path / "ratio.bdg"
    fil_a.write_text("chrI 0 10 4\nchrI 10 20 5\n", encoding="utf-8")
    fil_b.write_text("chrI 0 10 2\nchrI 10 20 0\n", encoding="utf-8")

    with pytest.raises(SystemExit) as error:
        main(
            [
                "--fil_A",
                str(fil_a),
                "--fil_B",
                str(fil_b),
                "--fil_out",
                str(fil_out),
                "--pseudo",
                pseudo,
            ],
        )

    assert "'A:B'" in str(error.value)
    assert not fil_out.exists()


def test_a_paired_pseudocount_reaches_both_tracks(tmp_path: Path) -> None:
    """
    'A:A' regularizes the zero denominator that a bare 'A' would leave.
    """

    fil_a = tmp_path / "a.bdg"
    fil_b = tmp_path / "b.bdg"
    fil_out = tmp_path / "ratio.bdg"
    fil_a.write_text("chrI 0 10 5\n", encoding="utf-8")
    fil_b.write_text("chrI 0 10 0\n", encoding="utf-8")

    assert (
        main(
            [
                "--fil_A",
                str(fil_a),
                "--fil_B",
                str(fil_b),
                "--fil_out",
                str(fil_out),
                "--pseudo",
                "1:1",
                "--dp",
                "3",
            ],
        )
        == 0
    )
    assert fil_out.read_text(encoding="utf-8") == "chrI\t0\t10\t6\n"


def test_the_ratio_vocabulary_is_exactly_the_ruled_set() -> None:
    """
    Every accepted spelling is one the vocabulary ruling admits.

    The ratio methods differ from each other in scale, not in adjustment, so
    the plain quotient is named 'linear' rather than 'unadj'. That also keeps
    the token 'unadj' meaning one thing across the tools, where it names a
    deposition rather than a scale.
    """

    # A list comparison pins membership, count and sequence at once, since
    # argparse renders the choices display straight from this mapping.
    assert list(METHOD_CANON) == [
        "linear",
        "log2",
        "l2",
        "linear_r",
        "log2_r",
        "l2_r",
    ]
    assert set(METHOD_CANON.values()) == {
        "linear",
        "log2",
        "linear_r",
        "log2_r",
    }


@pytest.mark.parametrize(
    "retired",
    [
        "unadj",
        "unadj_r",
        "unadjusted",
        "r",
        "raw",
        "u",
        "s",
        "smp",
        "simple",
        "2",
        "lg2",
        "rr",
        "ur",
        "sr",
        "2r",
        "l2r",
        "lg2_r",
    ],
)
def test_retired_ratio_spellings_are_rejected(retired: str) -> None:
    assert retired not in METHOD_CANON

    with pytest.raises(SystemExit):
        parse_args(
            [
                "--fil_A",
                "a.bedGraph",
                "--fil_B",
                "b.bedGraph",
                "--fil_out",
                "o.bedGraph",
                "--method",
                retired,
            ],
        )


@pytest.mark.parametrize(
    ("spelling", "canonical"),
    [
        ("linear", "linear"),
        ("log2", "log2"),
        ("l2", "log2"),
        ("linear_r", "linear_r"),
        ("log2_r", "log2_r"),
        ("l2_r", "log2_r"),
    ],
)
def test_kept_ratio_spellings_standardize(
    spelling: str,
    canonical: str,
) -> None:
    args = parse_args(
        [
            "--fil_A",
            "a.bedGraph",
            "--fil_B",
            "b.bedGraph",
            "--fil_out",
            "o.bedGraph",
            "--method",
            spelling,
        ],
    )

    assert METHOD_CANON[args.method] == canonical


def test_method_help_names_every_ratio_method(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """
    Rendered help names every ratio method and spelling the parser accepts.
    """

    with pytest.raises(SystemExit):
        parse_args(["--help"])

    rendered = capsys.readouterr().out

    for spelling in METHOD_CANON:
        assert f"'{spelling}'" in rendered

    assert "'unadj'" not in rendered
    assert "'unadj_r'" not in rendered
