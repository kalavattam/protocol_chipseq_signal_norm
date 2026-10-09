#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_opt_applicability.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Tests of 'HELP.PARAMETER.APPLICABILITY' in the compute_signal Python layer.

Each case gives a parameter where it has no effect and checks that it is
refused or warned about, never ignored silently.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pytest

from protocol_chipseq_signal_norm.cli import compute_pseudo as pseudo
from protocol_chipseq_signal_norm.cli import compute_pseudo_deeptools as dtools
from protocol_chipseq_signal_norm.cli import compute_signal as signal
from protocol_chipseq_signal_norm.utilities.utils_cli import (
    check_opts_apply,
    find_supplied,
)
from protocol_chipseq_signal_norm.utilities.utils_stabilizer import (
    compute_pseudo_edger,
    determine_coef_eff,
    iter_vals_bdg,
    pick_stabilizer,
)

ROOT = Path(__file__).resolve().parents[3]
FIXTURES = ROOT / "tests" / "fixtures"
FIL_A = str(FIXTURES / "compute_pseudo" / "bedgraph" / "pair_A.bdg")
FIL_B = str(FIXTURES / "compute_pseudo" / "bedgraph" / "pair_B.bdg")
BAM_SE = str(FIXTURES / "compute_signal" / "bam" / "se" / "tiny_se.bam")
REF_FA = str(FIXTURES / "compute_signal" / "reference" / "tiny.fa")


# 'find_supplied' reads every spelling argparse accepts, and stops at '--'.
SPELLINGS = (
    (["--coef", "0.05"], {"coef"}),
    (["--coef=0.05"], {"coef"}),
    (["-c", "0.05"], {"coef"}),
    (["-c0.05"], {"coef"}),
    (["-sp", "#"], {"skp_pfx"}),
    (["-s", "max"], {"sym"}),
    (["--", "--coef"], set()),
)


@pytest.mark.parametrize(("argv", "expected"), SPELLINGS)
def test_find_supplied_reads_every_spelling(
    argv: list[str],
    expected: set[str],
) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--coef", dest="coef")
    parser.add_argument("-sp", "--skp_pfx", dest="skp_pfx")
    parser.add_argument("-s", "--sym", dest="sym")

    assert find_supplied(parser, argv) == expected


def test_check_opts_apply_refuses_warns_and_passes(
    capsys: pytest.CaptureFixture[str],
) -> None:
    opts = {"coef": "--coef"}

    with pytest.raises(SystemExit) as caught:
        check_opts_apply({"coef"}, "refuse", opts, "'A'", "'B'")

    assert str(caught.value.code) == (
        "'--coef' is for 'A'; it has no effect with 'B'."
    )

    check_opts_apply({"coef"}, "warn", opts, "'A'", "'B'")
    check_opts_apply(set(), "refuse", opts, "'A'", "'B'")

    assert capsys.readouterr().err == (
        "Note: '--coef' has no effect with 'B' and is ignored.\n"
    )


# Library parameters given where they cannot act raise.
RAISING = (
    pytest.param(
        lambda: pick_stabilizer([1.0], "frc_mdn_nz", qntl_pct=5.0),
        id="pick_stabilizer qntl_pct",
    ),
    pytest.param(
        lambda: pick_stabilizer([1.0], "min_nz", qntl_rule="floor"),
        id="pick_stabilizer qntl_rule",
    ),
    pytest.param(
        lambda: pick_stabilizer([1.0], "qntl_nz", coef=0.5),
        id="pick_stabilizer coef",
    ),
    pytest.param(
        lambda: determine_coef_eff("bogus", 0.5),
        id="determine_coef_eff method",
    ),
    pytest.param(
        lambda: list(iter_vals_bdg(FIL_A, eps=1.0, mode_nz="off")),
        id="iter_vals_bdg eps",
    ),
    pytest.param(
        lambda: list(iter_vals_bdg(FIL_A, mode_nz="open", nz_policy="off")),
        id="iter_vals_bdg mode_nz",
    ),
    pytest.param(
        lambda: compute_pseudo_edger(
            typ_sig="count",
            n_ovlp_a=6.0,
            n_ovlp_b=18.0,
            n_frg_a=2.0,
        ),
        id="compute_pseudo_edger n_frg",
    ),
    pytest.param(
        lambda: compute_pseudo_edger(
            typ_sig="norm",
            n_ovlp_a=6.0,
            n_ovlp_b=18.0,
            n_frg_a=2.0,
            n_frg_b=6.0,
            total_a=6.0,
        ),
        id="compute_pseudo_edger total",
    ),
    pytest.param(
        lambda: compute_pseudo_edger(
            typ_sig="CPM",
            n_ovlp_a=6.0,
            n_ovlp_b=18.0,
            scale_a=1.0,
        ),
        id="compute_pseudo_edger scale",
    ),
)


@pytest.mark.parametrize("call", RAISING)
def test_library_refuses_a_parameter_that_cannot_act(call) -> None:
    with pytest.raises(ValueError, match="no effect|Unknown"):
        call()


def test_library_warns_about_an_unused_bin_width() -> None:
    with pytest.warns(UserWarning, match="'siz_bin' has no effect with 'CPM'"):
        compute_pseudo_edger(
            typ_sig="CPM",
            n_ovlp_a=6.0,
            n_ovlp_b=18.0,
            siz_bin=10,
        )


# CLI options given where they have no effect: argv and the message expected.
REFUSED_CLI = (
    pytest.param(
        pseudo.main,
        ["--method", "frc_mdn_nz", "--prior_count", "1"],
        "'--prior_count' is for '--method edger'",
        id="compute_pseudo prior_count",
    ),
    pytest.param(
        pseudo.main,
        ["--method", "qntl_nz", "--coef", "0.5"],
        "'--coef' is for '--method frc_mdn_nz'",
        id="compute_pseudo coef",
    ),
    pytest.param(
        pseudo.main,
        ["--method", "min_nz", "--qntl_nz", "5"],
        "'--qntl_nz' is for '--method qntl_nz'",
        id="compute_pseudo qntl_nz",
    ),
    pytest.param(
        pseudo.main,
        ["--method", "min_nz", "--mode_nz", "off", "--eps", "1"],
        "'--eps' is for '--mode_nz closed' or 'open'",
        id="compute_pseudo eps",
    ),
    pytest.param(
        dtools.main,
        ["--typ_sig", "CPM", "--sf_A", "1", "--sf_B", "1"],
        "'--sf_A' is for '--typ_sig RPGC'",
        id="compute_pseudo_deeptools sf",
    ),
)


@pytest.mark.parametrize(("main", "argv", "message"), REFUSED_CLI)
def test_pseudo_clis_refuse(main, argv: list[str], message: str) -> None:
    with pytest.raises(SystemExit) as caught:
        main(["--fil_A", FIL_A, "--fil_B", FIL_B, *argv])

    assert message in str(caught.value.code)


WARNED_CLI = (
    pytest.param(
        pseudo.main,
        ["--method", "edger", "--typ_sig", "count", "--siz_bin", "10"],
        "'--siz_bin' has no effect",
        id="compute_pseudo siz_bin",
    ),
    pytest.param(
        pseudo.main,
        ["--method", "edger", "--typ_sig", "count", "--n_frg_A", "2"],
        "'--n_frg_A' has no effect",
        id="compute_pseudo n_frg",
    ),
    pytest.param(
        dtools.main,
        ["--typ_sig", "CPM", "--siz_bin", "10"],
        "'--siz_bin' has no effect",
        id="compute_pseudo_deeptools siz_bin",
    ),
)


@pytest.mark.parametrize(("main", "argv", "message"), WARNED_CLI)
def test_pseudo_clis_warn_when_no_track_is_read(
    main,
    argv: list[str],
    message: str,
    capsys: pytest.CaptureFixture[str],
) -> None:
    status = main(
        [
            "--fil_A",
            FIL_A,
            "--fil_B",
            FIL_B,
            "--n_ovlp_A",
            "6",
            "--n_ovlp_B",
            "18",
            *argv,
        ],
    )

    assert status == 0
    assert message in capsys.readouterr().err


# 'compute_signal' options that do nothing for the output being written.
SIGNAL_REFUSED = (
    (["--method", "unadj"], "'--method' is for bedGraph output"),
    (["--scl_fct", "2"], "'--scl_fct' is for bedGraph output"),
    (["--siz_bin", "10"], "'--siz_bin' is for bedGraph output"),
)


@pytest.mark.parametrize(("argv", "message"), SIGNAL_REFUSED)
def test_compute_signal_refuses_signal_options_for_bed(
    argv: list[str],
    message: str,
    tmp_path: Path,
) -> None:
    with pytest.raises(SystemExit) as caught:
        signal.main(
            ["--fil_in", BAM_SE, "--fil_out", str(tmp_path / "x.bed"), *argv],
        )

    assert message in str(caught.value.code)


# A report-only run (no '--fil_out') builds no track either.
REPORT_ONLY_REFUSED = (
    (["--method", "unadj"], "'--method' is for bedGraph output"),
    (["--scl_fct", "2"], "'--scl_fct' is for bedGraph output"),
)


@pytest.mark.parametrize(("argv", "message"), REPORT_ONLY_REFUSED)
def test_compute_signal_refuses_signal_options_for_report_only(
    argv: list[str],
    message: str,
    tmp_path: Path,
) -> None:
    with pytest.raises(SystemExit) as caught:
        signal.main(
            [
                "--fil_in",
                BAM_SE,
                "--report_n_frg",
                str(tmp_path / "x.n_frg.txt"),
                *argv,
            ],
        )

    assert message in str(caught.value.code)
    assert "a report-only run" in str(caught.value.code)


# An ignored value is not checked, so an invalid one still only warns.
@pytest.mark.parametrize(
    "argv",
    [
        ["--dp", "3"],
        ["--engine", "window"],
        ["--siz_win", "20"],
        ["--dp", "-1"],
        ["--siz_win", "0"],
    ],
    ids=["dp", "engine", "siz_win", "invalid dp", "invalid siz_win"],
)
def test_compute_signal_warns_for_report_only(
    argv: list[str],
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    status = signal.main(
        [
            "--fil_in",
            BAM_SE,
            "--report_n_frg",
            str(tmp_path / "x.n_frg.txt"),
            *argv,
        ],
    )

    assert status == 0
    assert (tmp_path / "x.n_frg.txt").exists()
    assert (
        f"Note: '{argv[0]}' has no effect with a report-only run"
        in capsys.readouterr().err
    )


# '--siz_bin' sizes only the track and the overlap report.
@pytest.mark.parametrize(
    "fil_out", [None, "x.bed"], ids=["report-only", "bed"]
)
def test_compute_signal_refuses_siz_bin_without_overlap_report(
    fil_out: str | None,
    tmp_path: Path,
) -> None:
    argv = ["--fil_in", BAM_SE, "--report_n_frg", str(tmp_path / "x.txt")]

    if fil_out is not None:
        argv += ["--fil_out", str(tmp_path / fil_out)]

    with pytest.raises(SystemExit) as caught:
        signal.main([*argv, "--siz_bin", "20"])

    assert "'--siz_bin' is for bedGraph output or '--report_n_ovlp'" in str(
        caught.value.code,
    )


def test_compute_signal_takes_siz_bin_for_overlap_report(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    status = signal.main(
        [
            "--fil_in",
            BAM_SE,
            "--report_n_ovlp",
            str(tmp_path / "x.n_ovlp.txt"),
            "--siz_bin",
            "20",
        ],
    )

    assert status == 0
    assert "has no effect" not in capsys.readouterr().err


SIGNAL_WARNED = (
    pytest.param("x.bed", ["--dp", "3"], "'--dp'", id="bed dp"),
    pytest.param(
        "x.bed", ["--engine", "window"], "'--engine'", id="bed engine"
    ),
    pytest.param(
        "x.bdg", ["--siz_win", "20"], "'--siz_win'", id="chrom siz_win"
    ),
    pytest.param("x.bdg", ["--ref_fa", REF_FA], "'--ref_fa'", id="BAM ref_fa"),
    pytest.param(
        "x.bdg",
        ["--ref_fa", "missing.fa"],
        "'--ref_fa'",
        id="BAM missing ref_fa",
    ),
    pytest.param(
        "x.bdg", ["--siz_win", "0"], "'--siz_win'", id="chrom invalid siz_win"
    ),
    pytest.param(
        "x.bed", ["--siz_win", "0"], "'--siz_win'", id="bed invalid siz_win"
    ),
)


@pytest.mark.parametrize(("name", "argv", "flag"), SIGNAL_WARNED)
def test_compute_signal_warns_and_runs(
    name: str,
    argv: list[str],
    flag: str,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    status = signal.main(
        ["--fil_in", BAM_SE, "--fil_out", str(tmp_path / name), *argv],
    )

    assert status == 0
    assert f"Note: {flag} has no effect with" in capsys.readouterr().err


def test_compute_signal_stays_silent_for_options_that_apply(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    status = signal.main(
        [
            "--fil_in",
            BAM_SE,
            "--fil_out",
            str(tmp_path / "x.bdg"),
            "--engine",
            "window",
            "--siz_win",
            "20",
            "--dp",
            "3",
        ],
    )

    assert status == 0
    assert "has no effect" not in capsys.readouterr().err
