#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_input_floor_opt_applicability.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Test option applicability in 'compute_input_floor'.

Per 'HELP.PARAMETER.APPLICABILITY', the CLI refuses an option that the mode or
method does not use when it would change the result, and ignores a supporting
input with a warning, without checking it or passing it on. The library
function refuses or warns about the parameters whose 'None' default shows they
were given.
"""

import warnings
from pathlib import Path

import pytest

from protocol_chipseq_signal_norm.cli.compute_input_floor import (
    compute_input_floor,
    main,
)

ROOT = Path(__file__).resolve().parents[3]
FIXTURES = ROOT / "tests" / "fixtures" / "compute_signal"
BAM_SE = FIXTURES / "bam" / "se" / "tiny_se.bam"
CRAM_SE = FIXTURES / "cram" / "se" / "tiny_se.cram"
REFERENCE_FASTA = FIXTURES / "reference" / "tiny.fa"
BDG = FIXTURES / "bedgraph" / "ratio_A.bedGraph"

DIST_ONLY = (
    ("--method", "frc_mdn_nz"),
    ("--qntl_nz", "5"),
    ("--coef", "0.5"),
    ("--eps", "0.1"),
    ("--mode_nz", "open"),
    ("--floor", "1"),
)


@pytest.fixture
def bed(tmp_path: Path) -> Path:
    path = tmp_path / "input.bed"
    path.write_text("chrI\t0\t10\nchrI\t5\t15\n", encoding="utf-8")

    return path


def mode_args(mode: str, bed: Path) -> list[str]:
    """Return a minimal valid command line for one mode."""

    if mode == "dist":
        return ["--mode", "dist", "--fil_in", str(BDG)]

    if mode == "frag":
        return ["--mode", "frag", "--fil_in", str(bed)]

    return ["--mode", "norm"]


def refusal(option: str, where: str, now: str) -> str:
    return f"'{option}' is for {where}; it has no effect with {now}."


def note(option: str, now: str) -> str:
    return f"Note: '{option}' has no effect with {now} and is ignored.\n"


# CLI refusals: each stops before any work and names the option.
@pytest.mark.parametrize("mode", ("frag", "norm"))
@pytest.mark.parametrize(("option", "value"), DIST_ONLY)
def test_dist_options_are_refused_in_other_modes(
    bed: Path,
    mode: str,
    option: str,
    value: str,
) -> None:
    with pytest.raises(SystemExit) as error:
        main([*mode_args(mode, bed), option, value])

    assert str(error.value) == refusal(
        option,
        "'--mode dist'",
        f"'--mode {mode}'",
    )


@pytest.mark.parametrize("option", ("--siz_bin", "--siz_gen", "--siz-bin"))
def test_dimensions_are_refused_in_dist(bed: Path, option: str) -> None:
    with pytest.raises(SystemExit) as error:
        main([*mode_args("dist", bed), option, "30"])

    assert str(error.value) == refusal(
        option.replace("-bin", "_bin"),
        "'--mode frag' or '--mode norm'",
        "'--mode dist'",
    )


@pytest.mark.parametrize("mode", ("dist", "norm"))
@pytest.mark.parametrize("option", ("--flags_pe", "--flags_se"))
def test_fragment_flags_are_refused_once_outside_frag(
    bed: Path,
    capsys: pytest.CaptureFixture[str],
    mode: str,
    option: str,
) -> None:
    with pytest.raises(SystemExit) as error:
        main([*mode_args(mode, bed), option, "99"])

    assert str(error.value) == refusal(
        option,
        "'--mode frag'",
        f"'--mode {mode}'",
    )
    assert capsys.readouterr().err == ""


def test_fil_in_is_refused_in_norm(bed: Path) -> None:
    with pytest.raises(SystemExit) as error:
        main(["--mode", "norm", "--fil_in", str(bed)])

    assert str(error.value) == refusal(
        "--fil_in",
        "'--mode dist' or '--mode frag'",
        "'--mode norm'",
    )


@pytest.mark.parametrize("option", ("--flags_pe", "--flags_se"))
def test_fragment_flags_are_refused_with_bed_input(
    bed: Path,
    option: str,
) -> None:
    with pytest.raises(SystemExit) as error:
        main([*mode_args("frag", bed), option, "99"])

    assert str(error.value) == refusal(
        option,
        "BAM or CRAM input",
        "BED input",
    )


@pytest.mark.parametrize(
    ("args", "option", "where", "now"),
    (
        (
            ["--method", "frc_mdn_nz", "--qntl_nz", "5"],
            "--qntl_nz",
            "'--method qntl_nz'",
            "'--method frc_mdn_nz'",
        ),
        (
            ["--coef", "0.5"],
            "--coef",
            "'--method frc_mdn_nz', 'frc_avg_nz', or 'min_nz'",
            "'--method qntl_nz'",
        ),
    ),
)
def test_method_options_are_refused_with_the_other_method(
    bed: Path,
    args: list[str],
    option: str,
    where: str,
    now: str,
) -> None:
    with pytest.raises(SystemExit) as error:
        main([*mode_args("dist", bed), *args])

    assert str(error.value) == refusal(option, where, now)


@pytest.mark.parametrize(
    "args",
    (
        ["--method", "qntl_nz"],
        ["--qntl_nz", "1"],
        ["--eps", "0"],
        ["--mode_nz=closed"],
        ["-f0"],
    ),
)
def test_a_default_given_explicitly_is_still_refused(
    bed: Path,
    args: list[str],
) -> None:
    with pytest.raises(SystemExit) as error:
        main([*mode_args("frag", bed), *args])

    assert "has no effect with '--mode frag'" in str(error.value)


# CLI warnings: a supporting input is ignored, not checked, and the result is
# unchanged.
@pytest.mark.parametrize(
    ("mode", "args", "option", "now"),
    (
        ("dist", ["--ref_fa", "missing.fa"], "--ref_fa", "bedGraph input"),
        ("frag", ["--ref_fa", "missing.fa"], "--ref_fa", "BED input"),
        ("norm", ["--ref_fa", "missing.fa"], "--ref_fa", "'--mode norm'"),
        ("norm", ["--skp_pfx", "#"], "--skp_pfx", "'--mode norm'"),
        ("norm", ["--fmt_in", "bed"], "--fmt_in", "'--mode norm'"),
        (
            "dist",
            ["--fmt_in", "bedGraph"],
            "--fmt_in",
            "a named '--fil_in' path",
        ),
    ),
)
def test_supporting_options_warn_and_leave_the_result(
    bed: Path,
    capsys: pytest.CaptureFixture[str],
    mode: str,
    args: list[str],
    option: str,
    now: str,
) -> None:
    assert main(mode_args(mode, bed)) == 0
    baseline = capsys.readouterr().out

    # The CLI drops an ignored option, so the library raises no warning of its
    # own about it.
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert main([*mode_args(mode, bed), *args]) == 0

    captured = capsys.readouterr()

    assert captured.out == baseline
    assert captured.err == note(option, now)


@pytest.mark.parametrize("fil_in", (BAM_SE, CRAM_SE))
def test_skp_pfx_warns_with_alignment_input(
    capsys: pytest.CaptureFixture[str],
    fil_in: Path,
) -> None:
    base = [
        "--mode",
        "frag",
        "--fil_in",
        str(fil_in),
        "--ref_fa",
        str(REFERENCE_FASTA),
    ]

    if fil_in == BAM_SE:
        base = base[:4]

    assert main([*base, "--skp_pfx", "#"]) == 0

    label = "BAM" if fil_in == BAM_SE else "CRAM"

    assert capsys.readouterr().err == note("--skp_pfx", f"{label} input")


def test_verbose_norm_reports_no_input_options(
    capsys: pytest.CaptureFixture[str],
) -> None:
    assert main(["--verbose", "--mode", "norm", "--ref_fa", "missing.fa"]) == 0

    reported = [
        line.split(maxsplit=1)[0]
        for line in capsys.readouterr().err.splitlines()
        if line.startswith("--")
    ]

    assert reported == [
        "--verbose",
        "--mode",
        "--siz_bin",
        "--siz_gen",
        "--dp",
    ]


# Controls: an option where it acts draws no message.
@pytest.mark.parametrize(
    "args",
    (
        ["--mode", "dist", "--fil_in", str(BDG), "--method", "frc_mdn_nz"],
        [
            "--mode",
            "dist",
            "--fil_in",
            str(BDG),
            "--method",
            "frc_mdn_nz",
            "--coef",
            "0.5",
        ],
        ["--mode", "dist", "--fil_in", str(BDG), "--qntl_nz", "5"],
        ["--mode", "norm", "--siz_bin", "30"],
        ["--mode", "frag", "--fil_in", str(BAM_SE), "--flags_se", "0"],
        [
            "--mode",
            "frag",
            "--fil_in",
            str(CRAM_SE),
            "--ref_fa",
            str(REFERENCE_FASTA),
        ],
    ),
)
def test_options_where_they_act_draw_no_message(
    capsys: pytest.CaptureFixture[str],
    args: list[str],
) -> None:
    assert main(args) == 0
    assert capsys.readouterr().err == ""


def test_fmt_in_with_stdin_draws_no_message(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    monkeypatch.setattr("sys.stdin", BDG.open(encoding="utf-8"))

    assert main(["--mode", "dist", "--fil_in", "-", "--fmt_in", "bedGraph"]) == 0
    assert capsys.readouterr().err == ""


# Library: refuse with 'ValueError', warn with 'warnings.warn'.
def test_library_refuses_flags_with_bed_input(bed: Path) -> None:
    with pytest.raises(ValueError) as error:
        compute_input_floor(str(bed), 10, 100, mode="frag", single_flags={0})

    assert str(error.value) == (
        "'paired_flags' and 'single_flags' are for BAM or CRAM input; they "
        "have no effect with BED input."
    )


@pytest.mark.parametrize(
    ("kwargs", "message"),
    (
        (
            {"mode": "frag", "coef": 0.5},
            "'coef' is for mode 'dist'; it has no effect with 'frag'.",
        ),
        (
            {"mode": "norm", "coef": 0.5},
            "'coef' is for mode 'dist'; it has no effect with 'norm'.",
        ),
        (
            {"mode": "dist", "coef": 0.5},
            "'coef' is for methods 'frc_mdn_nz', 'frc_avg_nz', and 'min_nz'; "
            "it has no effect with 'qntl_nz'.",
        ),
        (
            {"mode": "dist", "paired_flags": {99}},
            "'paired_flags' and 'single_flags' are for mode 'frag'; they have "
            "no effect with 'dist'.",
        ),
        (
            {"mode": "norm", "single_flags": {0}},
            "'paired_flags' and 'single_flags' are for mode 'frag'; they have "
            "no effect with 'norm'.",
        ),
    ),
)
def test_library_refuses_given_parameters_that_cannot_act(
    kwargs: dict[str, object],
    message: str,
) -> None:
    fil_in = str(BDG) if kwargs["mode"] == "dist" else "-"

    with pytest.raises(ValueError) as error:
        compute_input_floor(fil_in, 10, 100, **kwargs)

    assert str(error.value) == message


@pytest.mark.parametrize(
    ("fil_in", "kwargs", "message"),
    (
        (
            BDG,
            {"mode": "dist", "ref_fa": "missing.fa"},
            "'ref_fa' has no effect with bedGraph input and is ignored.",
        ),
        (
            BAM_SE,
            {"mode": "frag", "ref_fa": "missing.fa"},
            "'ref_fa' has no effect with BAM input and is ignored.",
        ),
        (
            BDG,
            {"mode": "dist", "fmt_in": "bedGraph"},
            "'fmt_in' has no effect with a named input path and is ignored.",
        ),
    ),
)
def test_library_warns_about_supporting_parameters(
    bed: Path,
    fil_in: Path | None,
    kwargs: dict[str, object],
    message: str,
) -> None:
    path = str(fil_in or bed)
    expected = compute_input_floor(path, 10, 100, mode=kwargs["mode"])

    with pytest.warns(UserWarning) as record:
        result = compute_input_floor(path, 10, 100, **kwargs)

    assert [str(item.message) for item in record] == [message]
    assert result == expected


@pytest.mark.parametrize(
    ("fil_in", "kwargs"),
    (
        (BDG, {"mode": "dist", "method": "frc_mdn_nz", "coef": 0.5}),
        (BAM_SE, {"mode": "frag", "single_flags": {0}}),
        (
            CRAM_SE,
            {"mode": "frag", "ref_fa": str(REFERENCE_FASTA)},
        ),
    ),
)
def test_library_parameters_where_they_act_raise_no_warning(
    fil_in: Path,
    kwargs: dict[str, object],
) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        compute_input_floor(str(fil_in), 10, 100, **kwargs)
