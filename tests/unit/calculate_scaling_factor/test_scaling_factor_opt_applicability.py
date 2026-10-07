#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_scaling_factor_opt_applicability.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Tests of 'HELP.PARAMETER.APPLICABILITY' in the scaling-factor Python layer.

Each case gives a parameter where it has no effect and checks that it is
refused or warned about, never ignored silently, and that the result is the
same as without it.
"""

from __future__ import annotations

import warnings

import pytest

from protocol_chipseq_signal_norm.cli import (
    calculate_scaling_factor_siqchip as siq,
)
from protocol_chipseq_signal_norm.cli import (
    calculate_scaling_factor_spike as spike,
)
from protocol_chipseq_signal_norm.cli import parse_metadata_siqchip as meta

SIQ_ARGS = [
    "--mass_ip", "2", "--mass_in", "4", "--vol_all", "300", "--vol_in", "20",
    "--len_ip", "200", "--len_in", "180",
]  # fmt: skip
SPIKE_ARGS = ["--spike_ip", "10", "--spike_in", "20"]
MAINS = ["--main_ip", "100", "--main_in", "100"]
META_ARGS = [
    "--cfg", "parser.yml", "--alignment", "sample.bam",
    "--tbl_met", "metadata.tsv",
]  # fmt: skip


# Depths have no effect with an equation that has no depth terms.
@pytest.mark.parametrize("eqn", ["5nd", "6nd"])
@pytest.mark.parametrize("flag", ["--dep_ip", "--dep_in"])
def test_siqchip_refuses_depths_without_depth_terms(
    eqn: str,
    flag: str,
) -> None:
    with pytest.raises(SystemExit) as error:
        siq.main([*SIQ_ARGS, "--eqn", eqn, flag, "3"])

    assert str(error.value) == (
        f"'{flag}' is for '--eqn 5' or '--eqn 6'; it has no effect with "
        f"'--eqn {eqn}'."
    )


@pytest.mark.parametrize("eqn", ["5", "6"])
def test_siqchip_accepts_depths_with_depth_terms(
    eqn: str,
    capsys: pytest.CaptureFixture[str],
) -> None:
    assert (
        siq.main([*SIQ_ARGS, "--eqn", eqn, "--dep_ip", "3", "--dep_in", "5"])
        == 0
    )
    assert capsys.readouterr().err == ""


@pytest.mark.parametrize("eqn", ["5nd", "6nd"])
def test_calculate_alpha_refuses_depths_without_depth_terms(eqn: str) -> None:
    with pytest.raises(ValueError, match="have no effect with 'eqn'"):
        siq.calculate_alpha(eqn, 2, 4, 300, 20, 3, None, 200, 180)


# A main count a coefficient does not use is ignored with a note, and the
# coefficient is the same as without it.
UNUSED_MAINS = (
    ("chiprx_alpha_ratio", ["--main_ip", "--main_in"]),
    ("chiprx_alpha_ip", ["--main_ip", "--main_in"]),
    ("chiprx_alpha_in", ["--main_ip", "--main_in"]),
    ("rxinput_alpha", ["--main_ip"]),
)


@pytest.mark.parametrize(("coef", "flags"), UNUSED_MAINS)
def test_spike_warns_about_unused_main_counts(
    coef: str,
    flags: list[str],
    capsys: pytest.CaptureFixture[str],
) -> None:
    needed = ["--main_in", "100"] if coef == "rxinput_alpha" else []

    assert spike.main(["--coef", coef, *SPIKE_ARGS, *needed]) == 0
    bare = capsys.readouterr()
    assert bare.err == ""

    assert spike.main(["--coef", coef, *SPIKE_ARGS, *MAINS]) == 0
    full = capsys.readouterr()

    assert full.out == bare.out
    assert full.err.splitlines() == [
        f"Note: '{flag}' has no effect with '--coef {coef}' and is ignored."
        for flag in flags
    ]


def test_spike_does_not_check_or_pass_on_an_unused_main_count(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """
    An ignored count is not validated, and the library never sees it.
    """

    argv = ["--coef", "chiprx_alpha_ratio", *SPIKE_ARGS, "--main_ip", "-1"]

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert spike.main(argv) == 0

    assert capsys.readouterr().out == "2\n"


MISSING_MAINS = (
    ("fractional", ["--main_in", "100"], "--main_ip"),
    ("fractional", ["--main_ip", "100"], "--main_in"),
    ("main_per_spike", ["--main_ip", "100"], "--main_in"),
    ("rxinput_alpha", [], "--main_in"),
    ("all", ["--main_in", "100"], "--main_ip"),
)


@pytest.mark.parametrize(("coef", "given", "flag"), MISSING_MAINS)
def test_spike_requires_main_counts_a_coefficient_uses(
    coef: str,
    given: list[str],
    flag: str,
) -> None:
    fmt = ["--format", "tsv"] if coef == "all" else []

    with pytest.raises(SystemExit) as error:
        spike.main(["--coef", coef, *fmt, *SPIKE_ARGS, *given])

    assert str(error.value) == f"'{flag}' is required for '--coef {coef}'."


def test_calculate_scaling_factors_warns_about_unused_main_counts() -> None:
    with pytest.warns(UserWarning, match="'main_ip' has no effect"):
        vals = spike.calculate_scaling_factors(
            100, 10, None, 20, required=("chiprx_alpha_ratio",)
        )

    assert vals == {"chiprx_alpha_ratio": 2.0}

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert spike.calculate_scaling_factors(
            None, 10, None, 20, required=("chiprx_alpha_ratio",)
        ) == {"chiprx_alpha_ratio": 2.0}


def test_calculate_scaling_factors_requires_used_main_counts() -> None:
    with pytest.raises(ValueError, match="'main_in' is required"):
        spike.calculate_scaling_factors(
            None, 10, None, 20, required=("rxinput_alpha",)
        )


# A validation run reads only the configuration.
@pytest.mark.parametrize(
    ("flag", "value"),
    [
        ("--alignment", ["sample.bam"]),
        ("--tbl_met", ["metadata.tsv"]),
        ("--shell", []),
        ("--skp_pfx", ["#"]),
    ],
)
def test_parse_metadata_warns_about_lookup_options_when_validating(
    flag: str,
    value: list[str],
    capsys: pytest.CaptureFixture[str],
) -> None:
    meta.parse_args(["--cfg", "parser.yml", "--validate_cfg", flag, *value])

    assert capsys.readouterr().err == (
        f"Note: '{flag}' has no effect with '--validate_cfg' and is ignored.\n"
    )


def test_parse_metadata_validation_alone_is_silent(
    capsys: pytest.CaptureFixture[str],
) -> None:
    meta.parse_args(["--cfg", "parser.yml", "--validate_cfg"])

    assert capsys.readouterr().err == ""


# Each underscore option also takes a hidden hyphen spelling, which parses to
# the same value and stays out of '--help'.
HIDDEN = (
    (siq, "--mass-ip", "mass_ip", 2.0),
    (siq, "--vol-all", "vol_all", 300.0),
    (siq, "--dep-ip", "dep_ip", 3),
    (siq, "--len-in", "len_in", 180.0),
    (siq, "--lib-vol-ip", "lib_vol_ip", 2.0),
    (spike, "--main-ip", "main_ip", 100),
    (spike, "--spike-in", "spike_in", 20),
    (meta, "--tbl-met", "tbl_met", "metadata.tsv"),
    (meta, "--skp-pfx", "skp_pfx", "#"),
)
BASE = {
    siq: [*SIQ_ARGS, "--eqn", "6", "--dep_ip", "3", "--dep_in", "5"],
    spike: [*SPIKE_ARGS, *MAINS],
    meta: META_ARGS,
}


@pytest.mark.parametrize(("module", "flag", "dest", "value"), HIDDEN)
def test_hidden_hyphen_spellings_parse_and_stay_out_of_help(
    module: object,
    flag: str,
    dest: str,
    value: object,
    capsys: pytest.CaptureFixture[str],
) -> None:
    # Drop the visible spelling and its value, if the base gives them.
    base = list(BASE[module])
    if f"--{dest}" in base:
        i = base.index(f"--{dest}")
        del base[i : i + 2]

    args = module.parse_args([*base, flag, str(value)])

    assert getattr(args, dest) == value
    assert dest in args.supplied

    with pytest.raises(SystemExit):
        module.parse_args(["--help"])

    assert flag not in capsys.readouterr().out


def test_hidden_validate_cfg_spelling_parses() -> None:
    args = meta.parse_args(["--cfg", "parser.yml", "--validate-cfg"])

    assert args.validate_cfg is True


# A required option stays required, whichever spelling would supply it.
@pytest.mark.parametrize(
    ("module", "argv", "flag"),
    [
        (
            siq,
            [a for a in SIQ_ARGS if a not in {"--len_in", "180"}],
            "--len_in",
        ),
        (spike, ["--spike_ip", "10"], "--spike_in"),
    ],
)
def test_required_options_are_still_required(
    module: object,
    argv: list[str],
    flag: str,
) -> None:
    with pytest.raises(SystemExit) as error:
        module.parse_args(argv)

    assert str(error.value) == f"'{flag}' is required."
