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

# Library sizes the fixture pair carries, so 'L_bar' is 12 and the two
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
    assert (
        one_finite.err == (
            "pseudo_A is nonfinite; mirroring pseudo_B in symmetric mode "
            "'bad'.\n"
        )
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
            True,
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
        "normalization": (
            ("-nm", "--normalization"),
            False,
            "_StoreAction",
            None,
            "CPM",
            (
                "CPM",
                "BPM",
                "RPKM",
                "None",
                "RPGC",
                "n",
                "nc",
                "nrm",
                "norm",
                "normalized",
            ),
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
        "lib_A": (
            ("-lA", "--lib_A"),
            False,
            "_StoreAction",
            "float",
            None,
            (),
            None,
        ),
        "lib_B": (
            ("-lB", "--lib_B"),
            False,
            "_StoreAction",
            "float",
            None,
            (),
            None,
        ),
        "sf_A": (
            ("-sfA", "--sf_A"),
            False,
            "_StoreAction",
            "float",
            None,
            (),
            None,
        ),
        "sf_B": (
            ("-sfB", "--sf_B"),
            False,
            "_StoreAction",
            "float",
            None,
            (),
            None,
        ),
        "frg_A": (
            ("-gA", "--frg_A"),
            False,
            "_StoreAction",
            "float",
            None,
            (),
            None,
        ),
        "frg_B": (
            ("-gB", "--frg_B"),
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
        "prt_arg": (
            ("-pa", "--prt_arg"),
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

    # Name the method explicitly: the default is 'edger', whose
    # library-size path is not what this case exercises.
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
    Supplying both library sizes reads no track, so no width is resolved.
    """

    fil_a, fil_b = _bdg_pair(tmp_path)

    status = main(
        [
            "--verbose",
            "--method",
            "edger",
            "--fil_A",
            fil_a,
            "--fil_B",
            fil_b,
            "--lib_A",
            "100",
            "--lib_B",
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
            "--fil_A",
            FIL_A,
            "--fil_B",
            FIL_B,
            "--prt_jsn",
        ],
    )
    payload = json.loads(capsys.readouterr().out.splitlines()[1])

    assert status == 0
    assert payload["lib_sizes"] == {"A": 6.0, "B": 18.0}
    assert payload["prior_scaled"] == {"A": 1.0, "B": 3.0}


def test_json_payload_prior_is_not_derivable_for_normalized_coverage(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """
    Pin why 'prior_scaled' is emitted rather than left to be derived.

    'pseudo_i / scale_i' recovers it in every other mode. Under 'norm' both
    scale factors are 1.0 and the pseudocount is symmetric, so that quotient
    returns the shared pseudocount instead. The fragment counts invert the
    library-size imbalance here (3:1 against the tracks' 1:3), so a prior that
    tracked the tracks could not produce these values.
    """

    status = main(
        [
            "--method",
            "edger",
            "--normalization",
            "nc",
            "--fil_A",
            FIL_A,
            "--fil_B",
            FIL_B,
            "--frg_A",
            "3",
            "--frg_B",
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
        ["--method", "edger", "--fil_A", FIL_A, "--fil_B", FIL_B, *tokens],
    )
    note = capsys.readouterr().err

    assert status == 0
    assert "do not apply to '--method edger'" in note


# Options whose short forms begin with '-s', which is '--sym'. Resolving by
# the longest registered prefix is what keeps these from being reported as
# '--sym' carrying an attached value. Each row is the whole argument list, so
# the test does not branch to assemble one.
SYM_PREFIXED = (
    ("-sp", ("-sp", "#,track,browser")),
    ("-sb", ("-sb", "10")),
    ("-sfA", ("-nm", "RPGC", "-sfA", "0.5", "-sfB", "0.5")),
    ("-sfB", ("-nm", "RPGC", "-sfA", "0.5", "-sfB", "0.5")),
)


@pytest.mark.parametrize(("short_form", "tokens"), SYM_PREFIXED)
def test_warn_inapplicable_does_not_confuse_sym_with_longer_options(
    short_form: str,
    tokens: tuple[str, ...],
    capsys: pytest.CaptureFixture[str],
) -> None:
    """
    '-s' prefixes '-sp', '-sb', '-sfA', and '-sfB', which all apply to edgeR.

    Reporting one of them as '--sym' would send a user looking for a flag they
    never passed.
    """

    status = main(
        ["--method", "edger", "--fil_A", FIL_A, "--fil_B", FIL_B, *tokens],
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
        action.dest: tuple(action.option_strings) for action in actions
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


# Single-track mode is selected from '--fil_B' and '--lib_B', so a library size
# supplied to skip a read also decides the mode. Each row is a whole argument
# list parsed by the real parser, so a flag rename fails here rather than
# quietly reporting two tracks.
ONE_TRACK_ARGV = (
    pytest.param((), True, id="neither_B_option"),
    pytest.param(("--fil_B", FIL_B), False, id="fil_B_long"),
    pytest.param(("-fB", FIL_B), False, id="fil_B_short"),
    pytest.param(("--lib_B", "18"), False, id="lib_B_long"),
    pytest.param(("-lB", "18"), False, id="lib_B_short"),
    pytest.param(("--fil_B", FIL_B, "--lib_B", "18"), False, id="both"),
)


@pytest.mark.parametrize(("tokens", "expected"), ONE_TRACK_ARGV)
def test_is_one_track_reads_both_b_options(
    tokens: tuple[str, ...],
    expected: bool,
) -> None:
    """
    '--lib_B' without '--fil_B' is two-track, the asymmetric case.

    Reading only '--fil_B' would mirror A onto B and emit one value for a run
    the user described with two library sizes.
    """

    args = parse_args(["--fil_A", FIL_A, *tokens])

    assert compute_pseudo._is_one_track(args) is expected


def test_is_one_track_treats_a_zero_library_size_as_supplied() -> None:
    """
    '--lib_B 0' is falsy but supplied, and the check reads 'is None' for it.

    A truthiness test would report single-track and discard the named B track.
    """

    args = parse_args(["--fil_A", FIL_A, "--lib_B", "0"])

    assert args.lib_B == 0.0
    assert compute_pseudo._is_one_track(args) is False


# The other B options carry a value for a track rather than asserting one
# exists, so naming them alone leaves the run single-track.
ONE_TRACK_UNRELATED = (
    pytest.param(("--frg_B", "1000"), id="frg_B"),
    pytest.param(("--sf_B", "0.5"), id="sf_B"),
)


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
    what lets every reader treat it like '--lib_B'. Dropping the declaration
    fails this assertion with 'AttributeError'.
    """

    args = parse_args(["--fil_A", FIL_A])

    assert args.fil_B is None
    assert args.lib_B is None
    assert compute_pseudo._is_one_track(args) is True


# Library sizes the fixture pair carries. Expected pseudocounts below are
# written from edgeR's published definition rather than by calling the
# implementation under test, so the assertion is evidence and not arithmetic.
LIB_A = 6.0
LIB_B = 18.0
PRIOR = 2.0


def _edger_pseudo(mean_lib: float, prior: float = PRIOR) -> float:
    """
    Return edgeR's pseudocount on the CPM scale for one mean library size.

    'y0_i = prior * L_i / mean(L)' and 's_i = 1e6 / (L_i + 2 * y0_i)' give
    'p_i = s_i * y0_i', which reduces to the expression below and carries no
    per-sample index. That is why an edgeR pair is symmetric.

    Parameters
    ----------
    mean_lib : float
        Mean library size over the columns being compared.
    prior : float
        Nominal 'prior.count'.

    Returns
    -------
    pseudo : float
        Pseudocount in normalized units.
    """

    return 1e6 * prior / (mean_lib + 2.0 * prior)


def test_single_track_emits_one_value_rather_than_a_pair(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """
    Omitting both '--fil_B' and '--lib_B' emits one value, not 'A:B'.
    """

    status = main(["--method", "edger", "--fil_A", FIL_A])

    emitted = capsys.readouterr().out.strip()

    assert status == 0
    assert ":" not in emitted
    assert math.isclose(
        float(emitted), _edger_pseudo(LIB_A), rel_tol=1e-12
    )


def test_single_track_value_is_not_the_two_track_value(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """
    A track's one-track pseudocount is not its two-track pseudocount.

    'mean(L)' is the track's own library size in one framing and the mean of
    both in the other, so a refactor that quietly made them equal would still
    produce plausible numbers. This is the assertion that catches it.
    """

    main(["--method", "edger", "--fil_A", FIL_A])

    alone = float(capsys.readouterr().out.strip())

    main(["--method", "edger", "--fil_A", FIL_A, "--fil_B", FIL_B])

    paired = float(capsys.readouterr().out.strip().split(":")[0])

    assert math.isclose(alone, _edger_pseudo(LIB_A), rel_tol=1e-12)
    assert math.isclose(
        paired, _edger_pseudo((LIB_A + LIB_B) / 2.0), rel_tol=1e-12
    )
    assert alone != paired


def test_single_track_refuses_prt_arg(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """
    '--prt_arg' is refused rather than silently truncated.

    'bamCompare' takes a pair for each of '--scaleFactors' and '--pseudocount',
    and 'bamCoverage --scaleFactor' accepts no pseudocount, so emitting a
    one-track argument string would drop the value the user asked for. Assert
    nothing was emitted, not merely that it failed.
    """

    with pytest.raises(SystemExit) as refused:
        main(["--method", "edger", "--fil_A", FIL_A, "--prt_arg"])

    emitted = capsys.readouterr().out
    message = str(refused.value)

    assert "no single-track form" in message
    assert "--scaleFactors" not in emitted
    assert "--pseudocount" not in emitted


def test_prt_arg_writes_the_two_track_argument_string(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """
    The two-track form is a deepTools argument string, not a bare pair.
    """

    status = main(
        ["--method", "edger", "--fil_A", FIL_A, "--fil_B", FIL_B, "--prt_arg"],
    )

    emitted = capsys.readouterr().out.strip()
    expected = _edger_pseudo((LIB_A + LIB_B) / 2.0)
    pseudo = emitted.split("--pseudocount ")[1].split()

    assert status == 0
    assert emitted.startswith("--scaleFactors ")
    assert all(
        math.isclose(float(value), expected, rel_tol=1e-12) for value in pseudo
    )


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

    status = main(["--method", "edger", "--prt_jsn", *argv])

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
    assert payload["lib_sizes"]["B"] == payload["lib_sizes"]["A"]
    assert payload["scale_factors"]["B"] == payload["scale_factors"]["A"]
    assert pseudocounts["pseudo_B"] == pseudocounts["pseudo_A"]
    assert math.isclose(
        payload["lib_sizes"]["A"], LIB_A, rel_tol=1e-12
    )


def test_json_keeps_one_shape_across_both_modes(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """
    The schema does not change shape between one- and two-track runs.

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

    assert shapes[0] == shapes[1]
    assert alone["one_track"] is True
    assert paired["one_track"] is False
    assert paired["pseudocounts"]["pseudo_A"] != (
        alone["pseudocounts"]["pseudo_A"]
    )


# Each row is a whole normalization case: the flag that single-track mode still
# requires, a value for it, and the flag a two-track run would also need. The
# message must name only the first.
SINGLE_TRACK_REQUIRED = (
    pytest.param("norm", "--frg_A", "100", "--frg_B", id="norm"),
    pytest.param("RPGC", "--sf_A", "0.5", "--sf_B", id="RPGC"),
)


@pytest.mark.parametrize(
    ("normalization", "flag", "value", "paired_flag"), SINGLE_TRACK_REQUIRED
)
def test_single_track_requires_only_the_a_side_flag(
    normalization: str,
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
            "--fil_A",
            FIL_A,
            "--normalization",
            normalization,
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
                "--fil_A",
                FIL_A,
                "--normalization",
                normalization,
            ],
        )

    message = str(missing.value)

    assert status == 0
    assert f"'{flag}'" in message
    assert paired_flag not in message
    assert "both" not in message
