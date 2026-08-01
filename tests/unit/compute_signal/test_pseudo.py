#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_pseudo.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
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
        one_finite.err
        == "pseudo_A is nonfinite; mirroring pseudo_B in symmetric mode 'bad'.\n"
    )

    result = combine_pseudo_sym(math.nan, math.inf, "bad")
    neither_finite = capsys.readouterr()
    assert math.isnan(result[0])
    assert math.isinf(result[1])
    assert neither_finite.out == ""
    assert (
        neither_finite.err
        == "Both pseudocounts are nonfinite; returning as-is.\n"
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
            "==SUPPRESS==",
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
            "frc_mdn_nz",
            ("frc_mdn_nz", "qntl_nz", "frc_avg_nz", "min_nz"),
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
    assert single.out == "2.00\n"
    assert status_pair == 0
    assert pair.out.splitlines()[0] == "3.000:3.000"
    assert payload["pseudocounts"]["pseudo_A_str"] == "3.000"
    assert payload["pseudocounts"]["pseudo_B_str"] == "3.000"
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
    assert malformed.out == "2.00\n"
    assert malformed.err == ""

    assert main(["--fil_A", str(empty), "--prt_jsn"]) == 0
    strict_json = capsys.readouterr()
    assert strict_json.out == "nan\n"
    assert "No finite values in A after filtering" in strict_json.err
    assert strict_json.err.endswith(
        "Strict JSON disallows nan and inf; adjust '--floor' and '--coef', or just skip '--prt_jsn'.\n"
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
            ]
        )
