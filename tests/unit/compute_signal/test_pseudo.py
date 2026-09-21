#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_pseudo.py
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
import inspect
import json
import math
from pathlib import Path

import pytest

import protocol_chipseq_signal_norm.cli.compute_pseudo as compute_pseudo
from protocol_chipseq_signal_norm.cli.compute_pseudo import (
    combine_pseudo_sym,
    main,
    parse_args,
)

ROOT = Path(__file__).resolve().parents[3]
BEDGRAPH = ROOT / "tests" / "fixtures" / "compute_pseudo" / "bedgraph"

# Overlap counts the fixture pair carries, so 'L_bar' is 12 and the two
# per-sample priors are '2 * 6 / 12' and '2 * 18 / 12' exactly. A fixture is
# consumed by hard failure rather than by a skip, so a missing generation step
# fails loudly instead of turning the suite green.
FIL_A = str(BEDGRAPH / "pair_A.bdg")
FIL_B = str(BEDGRAPH / "pair_B.bdg")


def test_combine_pseudo_sym_returns_unmodified_for_none_mode() -> None:
    assert combine_pseudo_sym(1.0, 2.0, "none") == (1.0, 2.0)


def test_combine_pseudo_sym_applies_symmetric_modes() -> None:
    assert combine_pseudo_sym(1.0, 3.0, "max") == (3.0, 3.0)
    assert combine_pseudo_sym(1.0, 3.0, "min") == (1.0, 1.0)
    assert combine_pseudo_sym(1.0, 3.0, "arith") == (2.0, 2.0)
    assert combine_pseudo_sym(1.0, 4.0, "geom") == (2.0, 2.0)
    assert combine_pseudo_sym(1.0, 4.0, "harm") == (1.6, 1.6)
    assert combine_pseudo_sym(1.0, 3.0, "use_A") == (1.0, 1.0)
    assert combine_pseudo_sym(1.0, 3.0, "use_B") == (3.0, 3.0)


def test_combine_pseudo_sym_warns_for_defined_mean_fallbacks(
    capsys: pytest.CaptureFixture[str],
) -> None:
    assert combine_pseudo_sym(-1.0, 4.0, "geom") == (-1.0, -1.0)

    geometric = capsys.readouterr()

    assert geometric.out == ""
    assert geometric.err == (
        "Geometric mean undefined for negative values; falling back to "
        "min(pseudo_A, pseudo_B).\n"
    )

    assert combine_pseudo_sym(0.0, 4.0, "harm") == (0.0, 0.0)

    harmonic = capsys.readouterr()

    assert harmonic.out == ""
    assert harmonic.err == (
        "Harmonic mean undefined for nonpositive values; falling back to "
        "min(pseudo_A, pseudo_B).\n"
    )


def test_combine_pseudo_sym_mirrors_single_finite_value(
    capsys: pytest.CaptureFixture[str],
) -> None:
    assert combine_pseudo_sym(math.nan, 2.0, "max") == (2.0, 2.0)
    assert "nonfinite" in capsys.readouterr().err


def test_combine_pseudo_sym_nonfinite_paths_do_not_validate_mode(
    capsys: pytest.CaptureFixture[str],
) -> None:
    assert combine_pseudo_sym(math.nan, 2.0, "bad") == (2.0, 2.0)

    one_finite = capsys.readouterr()

    assert one_finite.out == ""
    assert one_finite.err == (
        "pseudo_A is nonfinite; mirroring pseudo_B in symmetric mode 'bad'.\n"
    )

    result = combine_pseudo_sym(math.nan, math.inf, "bad")
    neither_finite = capsys.readouterr()

    assert math.isnan(result[0])
    assert math.isinf(result[1])
    assert neither_finite.out == ""
    assert neither_finite.err == (
        "Both pseudocounts are nonfinite; returning as-is.\n"
    )


def test_combine_pseudo_sym_rejects_unknown_mode() -> None:
    with pytest.raises(ValueError, match="Unknown"):
        combine_pseudo_sym(1.0, 2.0, "bad")


def test_parser_preserves_complete_action_contract(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    parser_type = compute_pseudo.CapArgumentParser
    captured: dict[str, object] = {}

    def capture_parser(*args: object, **kwargs: object) -> object:
        parser = parser_type(*args, **kwargs)
        captured["parser"] = parser

        return parser

    monkeypatch.setattr(compute_pseudo, "CapArgumentParser", capture_parser)
    parse_args(["--fil_A", "signal_A.bdg"])
    parser = captured["parser"]
    actions = getattr(parser, "_actions")
    actual = {
        action.dest: (
            tuple(action.option_strings),
            action.required,
            type(action).__name__,
            None if action.type is None else action.type.__name__,
            action.default,
            tuple(action.choices or ()),
            action.const,
        )
        for action in actions
        if action.help is not argparse.SUPPRESS
    }

    assert actual == {
        "help": (
            ("-h", "--help"),
            False,
            "_HelpAction",
            None,
            "==SUPPRESS==",
            (),
            None,
        ),
        "verbose": (
            ("-v", "--verbose"),
            False,
            "_StoreTrueAction",
            None,
            False,
            (),
            True,
        ),
        "fil_A": (
            ("-fA", "--fil_A"),
            False,
            "_StoreAction",
            None,
            "==SUPPRESS==",
            (),
            None,
        ),
        "fil_B": (
            ("-fB", "--fil_B"),
            False,
            "_StoreAction",
            None,
            None,
            (),
            None,
        ),
        "skp_pfx": (
            ("-sp", "--skp_pfx"),
            False,
            "_StoreAction",
            "str",
            "#,track,browser",
            (),
            None,
        ),
        "method": (
            ("-m", "--method"),
            False,
            "_StoreAction",
            None,
            "edger",
            ("edger", "frc_mdn_nz", "qntl_nz", "frc_avg_nz", "min_nz"),
            None,
        ),
        "qntl_nz": (
            ("-q", "--qntl_nz"),
            False,
            "_StoreAction",
            "float",
            1.0,
            (),
            None,
        ),
        "coef": (
            ("-c", "--coef"),
            False,
            "_StoreAction",
            "float",
            None,
            (),
            None,
        ),
        "floor": (
            ("-fl", "--floor"),
            False,
            "_StoreAction",
            "float",
            0.0,
            (),
            None,
        ),
        "eps": (
            ("-e", "--eps"),
            False,
            "_StoreAction",
            "float",
            0.0,
            (),
            None,
        ),
        "mode_nz": (
            ("-mz", "--mode_nz"),
            False,
            "_StoreAction",
            None,
            "closed",
            ("closed", "open", "off"),
            None,
        ),
        "sym": (
            ("-s", "--sym"),
            False,
            "_StoreAction",
            None,
            "none",
            ("none", "max", "min", "arith", "geom", "harm", "use_A", "use_B"),
            None,
        ),
        "substrate": (
            ("-su", "--substrate"),
            False,
            "_StoreAction",
            "_check_substrate",
            "norm",
            ("unadj", "frag", "norm", "nc", "count", "cpm"),
            None,
        ),
        "prior_count": (
            ("-pc", "--prior_count"),
            False,
            "_StoreAction",
            "float",
            2.0,
            (),
            None,
        ),
        "siz_bin": (
            ("-sb", "--siz_bin"),
            False,
            "_StoreAction",
            "int",
            None,
            (),
            None,
        ),
        "n_ovlp_A": (
            ("-noA", "--n_ovlp_A"),
            False,
            "_StoreAction",
            "float",
            None,
            (),
            None,
        ),
        "n_ovlp_B": (
            ("-noB", "--n_ovlp_B"),
            False,
            "_StoreAction",
            "float",
            None,
            (),
            None,
        ),
        "n_frg_A": (
            ("-nfA", "--n_frg_A"),
            False,
            "_StoreAction",
            "float",
            None,
            (),
            None,
        ),
        "n_frg_B": (
            ("-nfB", "--n_frg_B"),
            False,
            "_StoreAction",
            "float",
            None,
            (),
            None,
        ),
        "dp": (("-dp", "--dp"), False, "_StoreAction", "int", 24, (), None),
        "prt_jsn": (
            ("-pj", "--prt_jsn"),
            False,
            "_StoreTrueAction",
            None,
            False,
            (),
            True,
        ),
    }


def test_help_channels_examples_and_semantic_order(
    capsys: pytest.CaptureFixture[str],
) -> None:
    with pytest.raises(SystemExit) as no_arguments:
        parse_args([])

    no_argument_capture = capsys.readouterr()

    with pytest.raises(SystemExit) as explicit_help:
        parse_args(["--help"])

    help_text = capsys.readouterr().out
    ordered = (
        "--help",
        "--verbose",
        "--fil_A",
        "--fil_B",
        "--skp_pfx",
        "--method",
        "--qntl_nz",
        "--coef",
        "--floor",
        "--eps",
        "--mode_nz",
        "--sym",
        "--dp",
        "--prt_jsn",
    )

    assert no_arguments.value.code == 0
    assert no_argument_capture.out == ""
    assert no_argument_capture.err.startswith("Usage\n-----\n  compute_pseudo")
    assert explicit_help.value.code == 0
    assert "Examples\n--------" in help_text
    assert "compute_pseudo --fil_A signal_A.bdg" in help_text
    assert "--sym max" in help_text

    positions = [help_text.index(item) for item in ordered]

    assert positions == sorted(positions)


def test_callable_docstring_matches_signature_and_examples() -> None:
    docstring = inspect.getdoc(combine_pseudo_sym)

    assert docstring is not None
    assert list(inspect.signature(combine_pseudo_sym).parameters) == [
        "pseudo_a",
        "pseudo_b",
        "mode",
    ]
    assert "Examples\n--------" in docstring
    assert 'mode="arith"' in docstring
    assert "(2.0, 2.0)" in docstring
    assert 'mode="use_A"' in docstring
    assert "(1.0, 1.0)" in docstring


def test_primary_pair_json_and_verbose_output_are_preserved(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    first = tmp_path / "first.bdg"
    second = tmp_path / "second.bdg"
    first.write_text("chrI\t0\t10\t2\n", encoding="utf-8")
    second.write_text("chrI\t0\t10\t3\n", encoding="utf-8")

    status = main(["--fil_A", str(first), "--method", "min_nz", "--dp", "2"])
    single = capsys.readouterr()

    status_pair = main(
        [
            "--verbose",
            "--fil_A",
            str(first),
            "--fil_B",
            str(second),
            "--method",
            "min_nz",
            "--sym",
            "max",
            "--dp",
            "3",
            "--prt_jsn",
        ],
    )
    pair = capsys.readouterr()
    reported = [
        line.split(maxsplit=1)[0]
        for line in pair.err.splitlines()
        if line.startswith("--")
    ]

    payload = json.loads(pair.out.splitlines()[1])

    assert status == 0
    assert single.out == "2\n"
    assert status_pair == 0
    assert pair.out.splitlines()[0] == "3:3"
    assert payload["pseudocounts"]["pseudo_A_str"] == "3"
    assert payload["pseudocounts"]["pseudo_B_str"] == "3"
    assert reported == [
        "--verbose",
        "--fil_A",
        "--fil_B",
        "--skp_pfx",
        "--method",
        "--coef",
        "--floor",
        "--eps",
        "--mode_nz",
        "--sym",
        "--dp",
        "--prt_jsn",
    ]


def test_main_skips_malformed_rows_and_handles_strict_json_failure(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    mixed = tmp_path / "mixed.bdg"
    empty = tmp_path / "empty.bdg"
    mixed.write_text(
        "chrI\t0\t10\tnot-a-number\nchrI\t10\t20\t2\n", encoding="utf-8"
    )
    empty.write_text("chrI\t0\t10\tnot-a-number\n", encoding="utf-8")

    assert (
        main(["--fil_A", str(mixed), "--method", "min_nz", "--dp", "2"]) == 0
    )

    malformed = capsys.readouterr()

    assert malformed.out == "2\n"
    assert malformed.err == ""

    # Name the method explicitly: the default is 'edger', whose overlap-count
    # path is not what this case exercises.
    status_empty = main(
        ["--fil_A", str(empty), "--method", "frc_mdn_nz", "--prt_jsn"],
    )

    assert status_empty == 0

    strict_json = capsys.readouterr()

    assert strict_json.out == "nan\n"
    assert "No finite values in A after filtering" in strict_json.err
    assert strict_json.err.endswith(
        "Strict JSON disallows nan and inf; adjust '--floor' and '--coef', or "
        "just skip '--prt_jsn'.\n",
    )


def test_main_rejects_invalid_quantile_with_a_stable_error(
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "input.bdg"
    input_path.write_text("chrI\t0\t10\t2\n", encoding="utf-8")

    with pytest.raises(
        SystemExit, match="'--qntl_nz' must be finite and in \\[0, 100\\]."
    ):
        main(
            [
                "--fil_A",
                str(input_path),
                "--method",
                "qntl_nz",
                "--qntl_nz",
                "101",
            ],
        )


def _bdg_pair(tmp_path: Path) -> tuple[str, str]:
    """
    Write a bedGraph pair on a 10 bp grid.
    """

    first = tmp_path / "first.bdg"
    second = tmp_path / "second.bdg"
    first.write_text(
        "chrI\t0\t10\t2\nchrI\t10\t20\t4\n",
        encoding="utf-8",
    )
    second.write_text(
        "chrI\t0\t10\t3\nchrI\t10\t20\t5\n",
        encoding="utf-8",
    )

    return str(first), str(second)


def test_verbose_banner_marks_an_inferred_bin_width(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """
    The banner reports the width used, not the flag that went unsupplied.
    """

    fil_a, fil_b = _bdg_pair(tmp_path)

    status = main(
        [
            "--verbose",
            "--method",
            "edger",
            "--substrate",
            "count",
            "--fil_A",
            fil_a,
            "--fil_B",
            fil_b,
        ],
    )
    banner = capsys.readouterr().err

    assert status == 0
    assert "--siz_bin 10  ## inferred from track ##" in banner
    assert "--siz_bin None" not in banner


def test_verbose_banner_reports_a_supplied_bin_width_unmarked(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    fil_a, fil_b = _bdg_pair(tmp_path)

    status = main(
        [
            "--verbose",
            "--method",
            "edger",
            "--substrate",
            "count",
            "--fil_A",
            fil_a,
            "--fil_B",
            fil_b,
            "--siz_bin",
            "10",
        ],
    )
    banner = capsys.readouterr().err

    assert status == 0
    assert "--siz_bin 10\n" in banner
    assert "inferred from track" not in banner


def test_verbose_banner_reports_an_unset_bin_width(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """
    Supplying both overlap counts reads no track, so no width is resolved.
    """

    fil_a, fil_b = _bdg_pair(tmp_path)

    status = main(
        [
            "--verbose",
            "--method",
            "edger",
            "--substrate",
            "count",
            "--fil_A",
            fil_a,
            "--fil_B",
            fil_b,
            "--n_ovlp_A",
            "100",
            "--n_ovlp_B",
            "200",
        ],
    )
    banner = capsys.readouterr().err

    assert status == 0
    assert "--siz_bin (unset)" in banner


def test_verbose_banner_survives_a_failure_during_resolution(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """
    A verbose run that fails still says what it was asked to do.

    The banner once printed only after the width resolved, so a contradicted
    '--siz_bin' produced the error with no record of the request beside it,
    which is the one case '--verbose' exists for.
    """

    fil_a, fil_b = _bdg_pair(tmp_path)

    with pytest.raises(SystemExit, match="disagrees with"):
        main(
            [
                "--verbose",
                "--method",
                "edger",
                "--substrate",
                "count",
                "--fil_A",
                fil_a,
                "--fil_B",
                fil_b,
                "--siz_bin",
                "20",
            ],
        )

    assert "--siz_bin 20" in capsys.readouterr().err


def test_json_payload_reports_the_per_sample_prior(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """
    The payload carries edgeR's 'y0_i', not only the pseudocount it feeds.

    The fixture pair is imbalanced 1:3, so a regression that dropped the
    per-sample scaling would return the nominal 'prior.count' for both tracks
    rather than the 1.0 and 3.0 asserted here.
    """

    status = main(
        [
            "--method",
            "edger",
            "--substrate",
            "count",
            "--fil_A",
            FIL_A,
            "--fil_B",
            FIL_B,
            "--prt_jsn",
        ],
    )
    payload = json.loads(capsys.readouterr().out.splitlines()[1])

    assert status == 0
    assert payload["n_ovlp"] == {"A": 6.0, "B": 18.0}
    assert payload["prior_scaled"] == {"A": 1.0, "B": 3.0}


def test_json_payload_prior_is_not_derivable_for_normalized_coverage(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """
    Pin why 'prior_scaled' is emitted rather than left to be derived.

    'pseudo_i / scale_i' recovers it under 'count' alone among this tool's
    substrates. Under 'norm' both scale factors are 1.0 and the pseudocount is
    symmetric, so that quotient returns the shared pseudocount instead. The
    fragment counts invert the overlap-count imbalance here (3:1 against the
    tracks' 1:3), so a prior that tracked the tracks could not produce these
    values.
    """

    status = main(
        [
            "--method",
            "edger",
            "--substrate",
            "nc",
            "--fil_A",
            FIL_A,
            "--fil_B",
            FIL_B,
            "--n_frg_A",
            "3",
            "--n_frg_B",
            "1",
            "--prt_jsn",
        ],
    )
    payload = json.loads(capsys.readouterr().out.splitlines()[1])

    assert status == 0
    assert payload["prior_scaled"] == {"A": 3.0, "B": 1.0}
    assert payload["scale_factors"] == {"A": 1.0, "B": 1.0}
    assert payload["pseudocounts"]["pseudo_A"] != payload["prior_scaled"]["A"]


# Every spelling argparse accepts for an option the edgeR path ignores. The
# note existed but read long forms only, so '-c 0.05' ran silently: the case
# its own docstring names as the reason it exists.
IGNORED_SPELLINGS = (
    ("--coef", "0.05"),
    ("--coef=0.05",),
    ("-c", "0.05"),
    ("-c=0.05",),
    ("-c0.05",),
    ("--qntl_nz", "5"),
    ("-q", "5"),
    ("--floor", "1"),
    ("-fl", "1"),
    ("--eps", "1"),
    ("-e", "1"),
    ("--mode_nz", "open"),
    ("-mz", "open"),
    ("--sym", "max"),
    ("-s", "max"),
)


@pytest.mark.parametrize(
    "tokens",
    IGNORED_SPELLINGS,
    ids=[" ".join(row) for row in IGNORED_SPELLINGS],
)
def test_warn_inapplicable_detects_every_spelling(
    tokens: tuple[str, ...],
    capsys: pytest.CaptureFixture[str],
) -> None:
    status = main(
        [
            "--method",
            "edger",
            "--substrate",
            "count",
            "--fil_A",
            FIL_A,
            "--fil_B",
            FIL_B,
            *tokens,
        ],
    )
    note = capsys.readouterr().err

    assert status == 0
    assert "do not apply to '--method edger'" in note


# Options whose short forms begin with '-s', which is '--sym'. Resolving by the
# longest registered prefix is what keeps these from being reported as '--sym'
# carrying an attached value. Each row is the whole argument list, so the test
# does not branch to assemble one.
SYM_PREFIXED = (
    ("-sp", ("-sp", "#,track,browser")),
    ("-sb", ("-sb", "10")),
)


@pytest.mark.parametrize(("short_form", "tokens"), SYM_PREFIXED)
def test_warn_inapplicable_does_not_confuse_sym_with_longer_options(
    short_form: str,
    tokens: tuple[str, ...],
    capsys: pytest.CaptureFixture[str],
) -> None:
    """
    '-s' prefixes '-sp', '-sb', and '-su', which all apply to edgeR.

    Reporting one of them as '--sym' would send a user looking for a flag they
    never passed.
    """

    status = main(
        [
            "--method",
            "edger",
            "--substrate",
            "count",
            "--fil_A",
            FIL_A,
            "--fil_B",
            FIL_B,
            *tokens,
        ],
    )

    assert status == 0
    assert "do not apply" not in capsys.readouterr().err


def test_warn_inapplicable_stays_silent_for_a_distribution_method(
    capsys: pytest.CaptureFixture[str],
) -> None:
    status = main(
        [
            "--method",
            "frc_mdn_nz",
            "--fil_A",
            FIL_A,
            "--fil_B",
            FIL_B,
            "-c",
            "0.05",
        ],
    )

    assert status == 0
    assert "do not apply" not in capsys.readouterr().err


def test_hidden_alias_restates_its_primary(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """
    Every hidden alias restates the action contract of the option it aliases.

    A hidden spelling is a separate 'add_argument' call sharing a 'dest', and
    argparse enforces nothing between the two. An alias that omits 'action'
    silently demands a value where its primary is a flag; one that omits
    'choices' accepts any string the primary would reject. Both shipped, so the
    agreement is asserted rather than assumed.
    """

    parser_type = compute_pseudo.CapArgumentParser
    captured: dict[str, object] = {}

    def capture_parser(*args: object, **kwargs: object) -> object:
        parser = parser_type(*args, **kwargs)
        captured["parser"] = parser

        return parser

    monkeypatch.setattr(compute_pseudo, "CapArgumentParser", capture_parser)
    parse_args(["--fil_A", "signal_A.bdg"])
    actions = getattr(captured["parser"], "_actions")

    def contract(action: argparse.Action) -> tuple[object, ...]:
        return (
            type(action).__name__,
            None if action.type is None else action.type.__name__,
            tuple(action.choices or ()),
            action.nargs,
            action.const,
        )

    primaries = {
        action.dest: action
        for action in reversed(actions)
        if action.help is not argparse.SUPPRESS
    }
    aliased = [
        action
        for action in actions
        if action.help is argparse.SUPPRESS and action.dest in primaries
    ]

    assert aliased, "no hidden aliases found; the guard would pass vacuously"

    mismatched = {
        action.option_strings[0]: (
            contract(primaries[action.dest]),
            contract(action),
        )
        for action in aliased
        if contract(primaries[action.dest]) != contract(action)
    }

    assert mismatched == {}


def test_ignored_option_constants_match_the_parser(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """
    The constants restate the parser, so a test must prove they still agree.

    'compute_pseudo' cannot hand the parser to '_warn_inapplicable' without
    moving the 'add_argument' calls out of 'parse_args', which three alias
    auditors and 'PY.CLI.HELP.LAYOUT' all locate by that function's name.
    Restating the spellings and checking them here keeps both intact.
    """

    parser_type = compute_pseudo.CapArgumentParser
    captured: dict[str, object] = {}

    def capture_parser(*args: object, **kwargs: object) -> object:
        parser = parser_type(*args, **kwargs)
        captured["parser"] = parser

        return parser

    monkeypatch.setattr(compute_pseudo, "CapArgumentParser", capture_parser)
    parse_args(["--fil_A", "signal_A.bdg"])
    actions = getattr(captured["parser"], "_actions")
    registered = {
        action.dest: tuple(action.option_strings)
        for action in actions
        if action.help is not argparse.SUPPRESS
    }
    shorts = tuple(
        option
        for action in actions
        for option in action.option_strings
        if not option.startswith("--")
    )

    assert {
        dest: registered[dest] for dest in compute_pseudo.OPT_IGNORED_EDGER
    } == compute_pseudo.OPT_IGNORED_EDGER
    assert set(compute_pseudo.OPT_SHORT_ALL) == set(shorts)


# Single-track mode is selected from '--fil_B' and '--n_ovlp_B', so an overlap
# count supplied to skip a read also decides the mode. Each row is a whole
# argument list parsed by the real parser, so a flag rename fails here rather
# than quietly reporting two tracks.
ONE_TRACK_ARGV = (
    pytest.param((), True, id="neither_B_option"),
    pytest.param(("--fil_B", FIL_B), False, id="fil_B_long"),
    pytest.param(("-fB", FIL_B), False, id="fil_B_short"),
    pytest.param(("--n_ovlp_B", "18"), False, id="n_ovlp_B_long"),
    pytest.param(("-noB", "18"), False, id="n_ovlp_B_short"),
    pytest.param(("--fil_B", FIL_B, "--n_ovlp_B", "18"), False, id="both"),
)


@pytest.mark.parametrize(("tokens", "expected"), ONE_TRACK_ARGV)
def test_is_one_track_reads_both_b_options(
    tokens: tuple[str, ...],
    expected: bool,
) -> None:
    """
    '--n_ovlp_B' without '--fil_B' is two-track, the asymmetric case.

    Reading only '--fil_B' would mirror A onto B and emit one value for a run
    the user described with two fragment-bin overlap counts.
    """

    args = parse_args(["--fil_A", FIL_A, *tokens])

    assert compute_pseudo._is_one_track(args) is expected


def test_is_one_track_treats_a_zero_overlap_count_as_supplied() -> None:
    """
    '--n_ovlp_B 0' is falsy but supplied, and the check reads 'is None' for it.

    A truthiness test would report single-track and discard the named B track.
    """

    args = parse_args(["--fil_A", FIL_A, "--n_ovlp_B", "0"])

    assert args.n_ovlp_B == 0.0
    assert compute_pseudo._is_one_track(args) is False


# The other B options carry a value for a track rather than asserting one
# exists, so naming them alone leaves the run single-track.
ONE_TRACK_UNRELATED = (pytest.param(("--n_frg_B", "1000"), id="n_frg_B"),)


@pytest.mark.parametrize("tokens", ONE_TRACK_UNRELATED)
def test_is_one_track_ignores_the_other_b_options(
    tokens: tuple[str, ...],
) -> None:
    args = parse_args(["--fil_A", FIL_A, *tokens])

    assert compute_pseudo._is_one_track(args) is True


def test_is_one_track_reads_an_unsupplied_fil_b_as_none() -> None:
    """
    An unsupplied '--fil_B' is bound to None rather than left off.

    'CapArgumentParser' sets 'argument_default=argparse.SUPPRESS', so an option
    that declares no default is absent from the namespace and plain attribute
    access raises. '--fil_B' declares 'default=None' against that, which is
    what lets every reader treat it like '--n_ovlp_B'. Dropping the declaration
    fails this assertion with 'AttributeError'.
    """

    args = parse_args(["--fil_A", FIL_A])

    assert args.fil_B is None
    assert args.n_ovlp_B is None
    assert compute_pseudo._is_one_track(args) is True


# Overlap counts the fixture pair carries. Expected pseudocounts below are
# written from edgeR's published definition rather than by calling the
# implementation under test, so the assertion is evidence and not arithmetic.
OVLP_A = 6.0
OVLP_B = 18.0
PRIOR = 2.0


def _count_prior(
    n_ovlp: float,
    n_ovlp_mean: float,
    prior: float = PRIOR,
) -> float:
    """
    Return edgeR's pseudocount in count units for one overlap count.

    A 'count' track is edgeR's own 'y', so its pseudocount is edgeR's own
    'y0_i = prior * L_i / L_bar' with no rescaling. Unlike the closed
    substrates this carries a per-sample index, which is what makes the
    single-track and two-track values differ.

    Parameters
    ----------
    n_ovlp : float
        This sample's fragment-bin overlap count, edgeR's 'lib.size'.
    n_ovlp_mean : float
        Mean fragment-bin overlap count over the columns being compared, i.e.,
        edgeR's 'L_bar'.
    prior : float
        Nominal 'prior.count'.

    Returns
    -------
    pseudo : float
        Pseudocount in count units.
    """

    return prior * n_ovlp / n_ovlp_mean


def _edger_pseudo(n_ovlp_mean: float, prior: float = PRIOR) -> float:
    """
    Return edgeR's pseudocount on the CPM scale for one mean overlap count.

    'y0_i = prior * L_i / L_bar' and 's_i = 1e6 / (L_i + 2 * y0_i)' give
    'p_i = s_i * y0_i', which reduces to the expression below and carries no
    per-sample index. That is why an edgeR pair is symmetric.

    Parameters
    ----------
    n_ovlp_mean : float
        Mean fragment-bin overlap count over the columns being compared, i.e.,
        edgeR's 'L_bar'.
    prior : float
        Nominal 'prior.count'.

    Returns
    -------
    pseudo : float
        Pseudocount in normalized units.
    """

    return 1e6 * prior / (n_ovlp_mean + 2.0 * prior)


def _norm_like(tmp_path: Path) -> tuple[str, str]:
    """
    Write a fixture pair rescaled to sum to one, as normalized coverage does.
    """

    out = []

    for name in ("pair_A.bdg", "pair_B.bdg"):
        rows = [
            row.split("\t")
            for row in (BEDGRAPH / name).read_text().strip().splitlines()
        ]
        total = sum(float(row[3]) for row in rows)
        path = tmp_path / f"norm_{name}"
        path.write_text(
            "\n".join(
                "\t".join([*row[:3], repr(float(row[3]) / total)])
                for row in rows
            )
            + "\n",
        )
        out.append(str(path))

    return out[0], out[1]


def test_inferred_overlap_count_refuses_a_track_that_cannot_be_counts(
    tmp_path: Path,
) -> None:
    """
    Refuse an inferred 'L' that is impossible rather than merely surprising.

    Every fragment touches at least one bin, so 'L = k * N' is never below 'N'.
    A summed column total under '--n_frg' therefore proves the track is not a
    whole-count one, and inferring 'L' from it would scale the prior silently:
    measured, a normalized track gives a prior ten times the right one. This is
    an impossibility test, not a plausibility heuristic, so a genuine
    whole-count track cannot trip it.
    """

    fil_a, fil_b = _norm_like(tmp_path)

    with pytest.raises(SystemExit) as excinfo:
        main(
            [
                "--method",
                "edger",
                "--substrate",
                "norm",
                "--fil_A",
                fil_a,
                "--fil_B",
                fil_b,
                "--n_frg_A",
                "3",
                "--n_frg_B",
                "6",
            ],
        )

    assert "cannot be a whole-count track" in str(excinfo.value)


def test_inferred_overlap_count_still_works_for_a_whole_count_track(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """
    Keep the inference path the refusal is not meant to touch.

    'pseudo_recheck_2026-08-26/02_recheck_inferred.sh' passes deepTools-raw
    tracks and omits '--n_ovlp', which is exactly this shape. A whole-count
    track satisfies 'colSum >= N' by construction, so it stays green.
    """

    status = main(
        [
            "--method",
            "edger",
            "--substrate",
            "norm",
            "--fil_A",
            FIL_A,
            "--fil_B",
            FIL_B,
            "--n_frg_A",
            "3",
            "--n_frg_B",
            "6",
            "--dp",
            "6",
        ],
    )

    assert status == 0
    assert capsys.readouterr().out.strip() == "0.177778:0.177778"


def test_unadj_requires_the_overlap_counts() -> None:
    """
    Refuse an 'unadj' run that would conflate two different totals.

    'unadj' reads its own track to get the substrate total, the fragment base
    pairs. The fragment-bin overlap count 'L' that 'k' needs is a different
    number off the same file, so letting it default to the summed track
    silently substitutes one for the other: measured on the fixture pair,
    supplying 'L' ten times the column sum moves the prior by exactly ten.
    """

    with pytest.raises(SystemExit) as excinfo:
        main(
            [
                "--method",
                "edger",
                "--substrate",
                "unadj",
                "--fil_A",
                FIL_A,
                "--fil_B",
                FIL_B,
                "--n_frg_A",
                "3",
                "--n_frg_B",
                "6",
            ],
        )

    assert "'--n_ovlp_A' and '--n_ovlp_B'" in str(excinfo.value)


# Fragment counts for the fractional substrates. Their 1:2 ratio is not the
# tracks' 1:3, so a prior that read the overlap counts cannot land on these
# values.
FRG_A = 3.0
FRG_B = 6.0


def _fractional_prior(
    n_frg_a: float = FRG_A,
    n_frg_b: float = FRG_B,
    prior: float = PRIOR,
) -> float:
    """
    Return the closed fractional prior 'p_nc' for the fixture pair.

    The fractional substrates divide edgeR's prior by 'k_bar' to answer the
    under-dispersion of fractional deposition, giving
    'p_nc = prior / (k_bar * N_bar)' with 'k_i = L_i / N_i', averaged over the
    pair. Each member then scales 'p_nc' by its own column total. Written from
    that rule rather than from the tool's output, so the assertion is evidence
    and not arithmetic.

    Parameters
    ----------
    n_frg_a, n_frg_b : float
        Fragment count 'N' for each track.
    prior : float
        Nominal 'prior.count'.

    Returns
    -------
    prior_closed : float
        The prior for a substrate whose column total is one.
    """

    k_bar = 0.5 * (OVLP_A / n_frg_a + OVLP_B / n_frg_b)
    n_frg_bar = 0.5 * (n_frg_a + n_frg_b)

    return prior / (k_bar * n_frg_bar)


# One row per project substrate the other CLI tests never run: the extra
# arguments it needs, and its expected pair. 'unadj' scales 'p_nc' by the
# track's own base-pair total, 'frag' by the fragment count, and 'cpm' closes
# on 'L_bar' with no per-sample index. 'norm' and 'count' are covered
# elsewhere in this file.
#
# The fixture's 'unadj' totals happen to equal its overlap counts, so this row
# cannot separate 'T' from 'L'; 'test_unadj_requires_the_overlap_counts' covers
# that conflation by refusing the run instead.
SUBSTRATE_CLI_CASES = (
    pytest.param(
        "unadj",
        [
            "--n_ovlp_A",
            "6",
            "--n_ovlp_B",
            "18",
            "--n_frg_A",
            "3",
            "--n_frg_B",
            "6",
        ],
        (_fractional_prior() * OVLP_A, _fractional_prior() * OVLP_B),
        id="unadj",
    ),
    pytest.param(
        "frag",
        [
            "--n_frg_A",
            "3",
            "--n_frg_B",
            "6",
        ],
        (_fractional_prior() * FRG_A, _fractional_prior() * FRG_B),
        id="frag",
    ),
    pytest.param(
        "cpm",
        [],
        (
            PRIOR * 1e6 / (0.5 * (OVLP_A + OVLP_B)),
            PRIOR * 1e6 / (0.5 * (OVLP_A + OVLP_B)),
        ),
        id="cpm",
    ),
)


@pytest.mark.parametrize(
    ("substrate", "extra", "expected"), SUBSTRATE_CLI_CASES
)
def test_every_substrate_reaches_stdout_with_its_own_arithmetic(
    substrate: str,
    extra: list[str],
    expected: tuple[float, float],
    capsys: pytest.CaptureFixture[str],
) -> None:
    """
    Run each remaining substrate end to end, not only through the estimator.

    'test_stabilizer.py' covers the arithmetic per substrate; what this covers
    is the plumbing between the CLI and it, where a total can be passed where
    an overlap count belongs. That substitution is silent, so only a pinned
    pair catches it.
    """

    status = main(
        [
            "--method",
            "edger",
            "--substrate",
            substrate,
            "--fil_A",
            FIL_A,
            "--fil_B",
            FIL_B,
            *extra,
        ],
    )

    pseudo_a, pseudo_b = (
        float(part) for part in capsys.readouterr().out.strip().split(":")
    )

    assert status == 0
    assert math.isclose(pseudo_a, expected[0], rel_tol=1e-12)
    assert math.isclose(pseudo_b, expected[1], rel_tol=1e-12)


def test_single_track_emits_one_value_rather_than_a_pair(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """
    Omitting both '--fil_B' and '--n_ovlp_B' emits one value, not 'A:B'.
    """

    status = main(
        [
            "--method",
            "edger",
            "--substrate",
            "count",
            "--fil_A",
            FIL_A,
        ],
    )

    emitted = capsys.readouterr().out.strip()

    assert status == 0
    assert ":" not in emitted
    assert math.isclose(
        float(emitted), _count_prior(OVLP_A, OVLP_A), rel_tol=1e-12
    )


def test_single_track_value_is_not_the_two_track_value(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """
    A track's single-track pseudocount is not its two-track pseudocount.

    'L_bar' is the track's own fragment-bin overlap count in single-track mode
    and the mean of both in two-track mode, so a refactor that quietly made
    them equal would still produce plausible numbers. This is the assertion
    that catches it.
    """

    main(["--method", "edger", "--substrate", "count", "--fil_A", FIL_A])

    alone = float(capsys.readouterr().out.strip())

    main(
        [
            "--method",
            "edger",
            "--substrate",
            "count",
            "--fil_A",
            FIL_A,
            "--fil_B",
            FIL_B,
        ],
    )

    paired = float(capsys.readouterr().out.strip().split(":")[0])

    n_ovlp_mean = (OVLP_A + OVLP_B) / 2.0

    assert math.isclose(alone, _count_prior(OVLP_A, OVLP_A), rel_tol=1e-12)
    assert math.isclose(
        paired, _count_prior(OVLP_A, n_ovlp_mean), rel_tol=1e-12
    )
    assert alone != paired


def _edger_payload(
    argv: list[str],
    capsys: pytest.CaptureFixture[str],
) -> dict:
    """
    Run one edgeR invocation with '--prt_jsn' and return its JSON summary.

    Parameters
    ----------
    argv : list[str]
        Arguments after '--method edger'.
    capsys : pytest.CaptureFixture[str]
        Capture fixture owned by the calling test.

    Returns
    -------
    payload : dict
        Decoded JSON summary, which is the last line written to stdout.
    """

    status = main(
        [
            "--method",
            "edger",
            "--substrate",
            "count",
            "--prt_jsn",
            *argv,
        ],
    )

    assert status == 0

    return json.loads(capsys.readouterr().out.splitlines()[-1])


def test_json_mirrors_b_onto_a_in_single_track_mode(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """
    The B fields mirror A rather than being dropped, and 'one_track' says why.

    A consumer reading 'pseudo_B' must not have to know which mode produced the
    file, so the fields stay populated and equal.
    """

    payload = _edger_payload(["--fil_A", FIL_A], capsys)
    pseudocounts = payload["pseudocounts"]

    assert payload["one_track"] is True
    assert payload["fil_B"] is None
    assert payload["n_ovlp"]["B"] == payload["n_ovlp"]["A"]
    assert payload["scale_factors"]["B"] == payload["scale_factors"]["A"]
    assert pseudocounts["pseudo_B"] == pseudocounts["pseudo_A"]
    assert math.isclose(payload["n_ovlp"]["A"], OVLP_A, rel_tol=1e-12)


def test_json_keeps_one_shape_across_both_modes(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """
    The schema does not change shape between single-track and two-track runs.

    Compare keys recursively rather than spot-checking one field, so a dropped
    or added member fails here rather than in a consumer.
    """

    alone = _edger_payload(["--fil_A", FIL_A], capsys)
    paired = _edger_payload(["--fil_A", FIL_A, "--fil_B", FIL_B], capsys)
    shapes = []

    for payload in (alone, paired):
        shapes.append(
            {
                key: sorted(value) if isinstance(value, dict) else None
                for key, value in sorted(payload.items())
            },
        )

    pseudo_alone = alone["pseudocounts"]["pseudo_A"]
    pseudo_paired = paired["pseudocounts"]["pseudo_A"]

    assert shapes[0] == shapes[1]
    assert alone["one_track"] is True
    assert paired["one_track"] is False
    assert pseudo_paired != pseudo_alone


# Each row is a whole substrate case: the flag that single-track mode still
# requires, a value for it, and the flag a two-track run would also need. The
# message must name only the first.
#
# 'N' must not exceed the fixture's column sum of 6: every fragment touches at
# least one bin, so 'L >= N', and the inference guard refuses an impossible
# pair before this case can reach the message it is checking.
SINGLE_TRACK_REQUIRED = (
    pytest.param("norm", "--n_frg_A", "3", "--n_frg_B", id="norm"),
)


@pytest.mark.parametrize(
    ("substrate", "flag", "value", "paired_flag"), SINGLE_TRACK_REQUIRED
)
def test_single_track_requires_only_the_a_side_flag(
    substrate: str,
    flag: str,
    value: str,
    paired_flag: str,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """
    Single-track mode asks for the A flag alone, and says so when it is absent.

    A message naming both flags would send the user looking for a track they
    deliberately did not supply.
    """

    status = main(
        [
            "--method",
            "edger",
            "--substrate",
            "count",
            "--fil_A",
            FIL_A,
            "--substrate",
            substrate,
            flag,
            value,
        ],
    )

    capsys.readouterr()

    with pytest.raises(SystemExit) as missing:
        main(
            [
                "--method",
                "edger",
                "--substrate",
                "count",
                "--fil_A",
                FIL_A,
                "--substrate",
                substrate,
            ],
        )

    message = str(missing.value)

    assert status == 0
    assert f"'{flag}'" in message
    assert paired_flag not in message
    assert "both" not in message


def test_warn_inapplicable_reports_an_explicitly_passed_default(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """
    A flag whose value happens to be its default still warns when typed.

    Detection reads the supplied tokens rather than the parsed namespace, so it
    can tell a default apart from a choice. A user who typed '--sym none' asked
    the edgeR path for something it does not do, and silence would imply the
    request was honored.
    """

    status = main(
        [
            "--method",
            "edger",
            "--substrate",
            "count",
            "--fil_A",
            FIL_A,
            "--fil_B",
            FIL_B,
            "--sym",
            "none",
        ],
    )

    note = capsys.readouterr().err

    assert status == 0
    assert "--sym" in note
    assert "do not apply to '--method edger'" in note
