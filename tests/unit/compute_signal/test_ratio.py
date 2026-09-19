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
# - Anthropic Claude Code (Opus 5).
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

    assert scaled == 4.5
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


# Hidden hyphen aliases are separate 'add_argument' calls that must restate the
# primary's 'type', 'choices', and action; argparse enforces none of that. Each
# pair is exercised below, and a completeness guard fails when an alias reaches
# the parser without a row here.
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
