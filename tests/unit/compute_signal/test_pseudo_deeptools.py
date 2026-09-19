#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_pseudo_deeptools.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


import argparse
import json
from pathlib import Path

import pytest

import protocol_chipseq_signal_norm.cli.compute_pseudo_deeptools as interop
from protocol_chipseq_signal_norm.cli.compute_pseudo_deeptools import (
    SUB_CHOICES,
    main,
    parse_args,
)

ROOT = Path(__file__).resolve().parents[3]
BEDGRAPH = ROOT / "tests" / "fixtures" / "compute_pseudo" / "bedgraph"

# The same fixture pair 'test_pseudo.py' reads, carrying overlap counts 6 and
# 18, so 'L_bar' is 12 and the two per-sample priors are '2 * 6 / 12' and
# '2 * 18 / 12' exactly. A fixture is consumed by hard failure rather than by a
# skip, so a missing generation step fails loudly instead of turning the suite
# green.
FIL_A = str(BEDGRAPH / "pair_A.bdg")
FIL_B = str(BEDGRAPH / "pair_B.bdg")

# Every substrate this tool serves, with the extra arguments each one needs.
SUB_EXTRA = {
    "CPM": (),
    "BPM": (),
    "RPKM": (),
    "None": (),
    "RPGC": ("--sf_A", "0.7", "--sf_B", "1.3"),
}


def _capture(argv: list[str], capsys: pytest.CaptureFixture[str]) -> str:
    """
    Run the interop CLI and return its standard output.
    """

    assert main(argv) == 0

    return capsys.readouterr().out


def _build_parser() -> argparse.ArgumentParser:
    """
    Return the parser 'parse_args' builds, without parsing anything.
    """

    captured: dict[str, argparse.ArgumentParser] = {}
    parser_type = interop.CapArgumentParser

    def capture(*args: object, **kwargs: object) -> object:
        parser = parser_type(*args, **kwargs)
        captured["parser"] = parser

        return parser

    interop.CapArgumentParser = capture

    try:
        parse_args(["--fil_A", FIL_A])
    finally:
        interop.CapArgumentParser = parser_type

    return captured["parser"]


def test_substrate_choices_exclude_the_core_vocabulary() -> None:
    """
    Check that each CLI accepts only its own substrates.

    Post-split the tool is the namespace, so 'norm' and its aliases belong to
    'compute_pseudo' and reach nothing here.
    """

    assert SUB_CHOICES == ("CPM", "BPM", "RPKM", "RPGC", "None")

    for rejected in ("norm", "nc", "n", "nrm", "normalized"):
        with pytest.raises(SystemExit):
            parse_args(["--fil_A", FIL_A, "--substrate", rejected])


def test_no_normalization_spelling_survives() -> None:
    """
    Check that '--normalization' is excised rather than hidden.
    """

    for spelling in ("--normalization", "--normalisation", "-nm"):
        with pytest.raises(SystemExit):
            parse_args(["--fil_A", FIL_A, spelling, "CPM"])


# Frozen golden values, captured 2026-09-17 from this tool at its commit-1
# state, while the differential test against the pre-split 'compute_pseudo' was
# still green. They are therefore the pre-split arithmetic, and pinning them
# proves the split was a rename rather than a reimplementation.
#
# Never regenerate them: a regenerated value captures the drift it exists to
# catch, turning the guard into a tautology. A failure here means the
# arithmetic moved, which is what must be explained.
GOLDEN_PAIR = {
    "BPM": "125000:125000",
    "CPM": "125000:125000",
    "None": "1:3",
    "RPGC": "0.699999999999999955591079:3.900000000000000355271368",
    "RPKM": "12500000:12500000",
}
GOLDEN_ONE_TRACK_RPKM = "20000000"


@pytest.mark.parametrize("substrate", sorted(GOLDEN_PAIR))
def test_matches_the_pre_split_arithmetic(
    substrate: str,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """
    Pin the split: this tool still computes what 'compute_pseudo' computed.

    The comparison was live against the core tool until the core surface
    dropped its deepTools substrates, at which point the two could no longer be
    run side by side. Freezing the values keeps the guarantee past that window
    rather than losing it with the comparison.

    Equality is exact at full 'float64' precision, because last-bit drift is
    what this exists to catch.
    """

    out = _capture(
        [
            "--fil_A",
            FIL_A,
            "--fil_B",
            FIL_B,
            "--substrate",
            substrate,
            *SUB_EXTRA[substrate],
        ],
        capsys,
    )

    assert out == GOLDEN_PAIR[substrate] + "\n"


def test_matches_the_pre_split_arithmetic_for_one_track(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """
    Pin the single-track value, where 'L_bar' collapses to 'L_A'.
    """

    out = _capture(
        [
            "--fil_A",
            FIL_A,
            "--substrate",
            "RPKM",
        ],
        capsys,
    )

    assert out == GOLDEN_ONE_TRACK_RPKM + "\n"


def test_prt_arg_writes_the_two_track_argument_string(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """
    Check the 'bamCompare' argument string against the closed form.
    """

    out = _capture(
        [
            "--fil_A",
            FIL_A,
            "--fil_B",
            FIL_B,
            "--dp",
            "6",
            "--prt_arg",
        ],
        capsys,
    )

    # 'L_A' is 6 and 'L_B' is 18, so 'L_bar' is 12 and the per-sample priors
    # are 1 and 3. CPM then gives 's_A = 1e6 / (6 + 2)' and
    # 's_B = 1e6 / (18 + 6)', and 'p_i = s_i * prior_i' is 125000 either way,
    # which is the symmetry the estimator guarantees.
    assert out == (
        "--scaleFactors 125000.000000:41666.666667 --pseudocount 125000 "
        "125000\n"
    )


def test_single_track_refuses_prt_arg() -> None:
    """
    Check that the two-track argument string has no single-track spelling.
    """

    with pytest.raises(SystemExit) as excinfo:
        main(["--fil_A", FIL_A, "--prt_arg"])

    assert "no single-track form" in str(excinfo.value)


def test_rpgc_requires_both_scale_factors() -> None:
    """
    Check that 'RPGC' refuses to run without externally supplied factors.
    """

    with pytest.raises(SystemExit) as excinfo:
        main(
            [
                "--fil_A",
                FIL_A,
                "--fil_B",
                FIL_B,
                "--substrate",
                "RPGC",
            ],
        )

    assert "requires both '--sf_A' and '--sf_B'" in str(excinfo.value)

    with pytest.raises(SystemExit) as excinfo:
        main(["--fil_A", FIL_A, "--substrate", "RPGC"])

    assert "requires '--sf_A'" in str(excinfo.value)


def test_rpgc_applies_the_supplied_scale_factors(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """
    Check that 'RPGC' takes the prior's magnitude in proportional form.
    """

    out = _capture(
        [
            "--fil_A",
            FIL_A,
            "--fil_B",
            FIL_B,
            "--substrate",
            "RPGC",
            "--sf_A",
            "0.7",
            "--sf_B",
            "1.3",
            "--dp",
            "6",
        ],
        capsys,
    )

    # Each 'pseudo_i' is 's_i * prior_scaled_i', with priors 1 and 3.
    assert out == "0.7:3.9\n"


def test_rpkm_without_a_track_read_requires_siz_bin() -> None:
    """
    Check that a run reading no track still needs a width for 'RPKM'.
    """

    with pytest.raises(SystemExit) as excinfo:
        main(
            [
                "--fil_A",
                FIL_A,
                "--fil_B",
                FIL_B,
                "--substrate",
                "RPKM",
                "--n_bin_A",
                "6",
                "--n_bin_B",
                "18",
            ],
        )

    assert "'--siz_bin' is required" in str(excinfo.value)


def test_is_one_track_reads_an_unsupplied_fil_b_as_none() -> None:
    """
    Check that either second-track argument leaves single-track mode.
    """

    args = parse_args(["--fil_A", FIL_A])

    assert interop._is_one_track(args) is True

    args = parse_args(["--fil_A", FIL_A, "--fil_B", FIL_B])

    assert interop._is_one_track(args) is False

    # A supplied '--n_bin_B' describes a second track just as '--fil_B' does.
    args = parse_args(["--fil_A", FIL_A, "--n_bin_B", "18"])

    assert interop._is_one_track(args) is False


def test_json_keeps_one_shape_across_both_modes(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """
    Check that the JSON summary does not change shape between modes.
    """

    out = _capture(
        [
            "--fil_A",
            FIL_A,
            "--fil_B",
            FIL_B,
            "--prt_jsn",
        ],
        capsys,
    )
    pair = json.loads(out.splitlines()[-1])

    out = _capture(["--fil_A", FIL_A, "--prt_jsn"], capsys)
    single = json.loads(out.splitlines()[-1])

    assert pair.keys() == single.keys()
    assert pair["params"].keys() == single["params"].keys()
    assert pair["one_track"] is False
    assert single["one_track"] is True

    # The substrate is named for what it is. The fractional substrates own 'k',
    # and this tool serves none of them.
    assert pair["params"]["substrate"] == "CPM"
    assert "normalization" not in pair["params"]
    assert "k" not in pair


def test_hidden_aliases_restate_the_full_action_contract() -> None:
    """
    Check that every hidden alias restates what argparse does not inherit.

    A hidden alias is a separate action, so an omitted 'type' or 'choices'
    silently accepts what the canonical option refuses.
    """

    parser = _build_parser()
    canonical: dict[str, argparse.Action] = {}
    hidden: list[argparse.Action] = []

    for action in parser._actions:
        if action.help is argparse.SUPPRESS:
            hidden.append(action)
        else:
            canonical[action.dest] = action

    assert hidden, "no hidden aliases registered"

    for action in hidden:
        sibling = canonical[action.dest]

        assert type(action) is type(sibling), action.option_strings
        assert action.type is sibling.type, action.option_strings
        assert action.choices == sibling.choices, action.option_strings
        assert action.nargs == sibling.nargs, action.option_strings
        assert action.const == sibling.const, action.option_strings


def test_hidden_hyphen_aliases_reach_the_canonical_dest() -> None:
    """
    Check that each hidden hyphen spelling writes the canonical 'dest'.
    """

    args = parse_args(
        [
            "--fil-A",
            FIL_A,
            "--fil-B",
            FIL_B,
            "--prior-count",
            "3",
            "--siz-bin",
            "10",
            "--n-bin-A",
            "6",
            "--n-bin-B",
            "18",
            "--sf-A",
            "0.7",
            "--sf-B",
            "1.3",
            "--skp-pfx",
            "#",
            "--prt-jsn",
            "--prt-arg",
        ],
    )

    assert args.fil_A == FIL_A
    assert args.fil_B == FIL_B
    assert args.prior_count == 3.0
    assert args.siz_bin == 10
    assert args.n_bin_A == 6.0
    assert args.n_bin_B == 18.0
    assert args.sf_A == 0.7
    assert args.sf_B == 1.3
    assert args.skp_pfx == "#"
    assert args.prt_jsn is True
    assert args.prt_arg is True


def test_fil_a_is_required_through_the_hidden_spelling_too() -> None:
    """
    Check that the parsed value carries the requirement argparse cannot.
    """

    with pytest.raises(SystemExit) as excinfo:
        main(["--prior_count", "2"])

    assert "'--fil_A' is required" in str(excinfo.value)


def test_no_argument_prints_help_and_exits_zero(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """
    Check that a bare invocation prints help to stderr and exits zero.
    """

    with pytest.raises(SystemExit) as excinfo:
        parse_args([])

    assert excinfo.value.code == 0
    assert "compute_pseudo_deeptools" in capsys.readouterr().err


def test_parser_preserves_complete_action_contract() -> None:
    """
    Pin this tool's own public surface, rather than sharing the sibling's.

    A shared pin would drift silently: the two tools deliberately differ, so a
    single assertion could only track whichever one moved last.
    """

    parser = _build_parser()
    actual = {
        action.dest: (
            tuple(action.option_strings),
            type(action).__name__,
            None if action.type is None else action.type.__name__,
            action.default,
            tuple(action.choices or ()),
        )
        for action in parser._actions
        if action.help is not argparse.SUPPRESS
    }

    assert actual == {
        "help": (
            ("-h", "--help"),
            "_HelpAction",
            None,
            "==SUPPRESS==",
            (),
        ),
        "verbose": (
            ("-v", "--verbose"),
            "_StoreTrueAction",
            None,
            False,
            (),
        ),
        "fil_A": (
            ("-fA", "--fil_A"),
            "_StoreAction",
            None,
            "==SUPPRESS==",
            (),
        ),
        "fil_B": (
            ("-fB", "--fil_B"),
            "_StoreAction",
            None,
            None,
            (),
        ),
        "skp_pfx": (
            ("-sp", "--skp_pfx"),
            "_StoreAction",
            "str",
            "#,track,browser",
            (),
        ),
        "substrate": (
            ("-su", "--substrate"),
            "_StoreAction",
            None,
            "CPM",
            ("CPM", "BPM", "RPKM", "RPGC", "None"),
        ),
        "prior_count": (
            ("-pc", "--prior_count"),
            "_StoreAction",
            "float",
            2.0,
            (),
        ),
        "siz_bin": (
            ("-sb", "--siz_bin"),
            "_StoreAction",
            "int",
            None,
            (),
        ),
        "n_bin_A": (
            ("-nbA", "--n_bin_A"),
            "_StoreAction",
            "float",
            None,
            (),
        ),
        "n_bin_B": (
            ("-nbB", "--n_bin_B"),
            "_StoreAction",
            "float",
            None,
            (),
        ),
        "sf_A": (
            ("-sfA", "--sf_A"),
            "_StoreAction",
            "float",
            None,
            (),
        ),
        "sf_B": (
            ("-sfB", "--sf_B"),
            "_StoreAction",
            "float",
            None,
            (),
        ),
        "dp": (
            ("-dp", "--dp"),
            "_StoreAction",
            "int",
            24,
            (),
        ),
        "prt_jsn": (
            ("-pj", "--prt_jsn"),
            "_StoreTrueAction",
            None,
            False,
            (),
        ),
        "prt_arg": (
            ("-pa", "--prt_arg"),
            "_StoreTrueAction",
            None,
            False,
            (),
        ),
    }
