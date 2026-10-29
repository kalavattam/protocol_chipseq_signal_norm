#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_merge_bins.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6);
# - Anthropic Claude Code (Opus 5.5).
#
# Distributed under the MIT license.


from pathlib import Path

import pytest

from protocol_chipseq_signal_norm.cli.merge_bins_bdg import (
    format_value,
    merge_bins,
    validate_merge_options,
)
from protocol_chipseq_signal_norm.utilities.utils_io import DEF_SKP_PFX

ROOT = Path(__file__).resolve().parents[3]
MERGE_BINS = ROOT / "tests" / "fixtures" / "bedgraph" / "merge_bins"


def run_merge(
    tmp_path: Path,
    name: str,
    *,
    decimal_places: int | None = None,
    eps: float | None = None,
    skp_pfx: tuple[str, ...] = DEF_SKP_PFX,
) -> str:
    fil_in = MERGE_BINS / name
    fil_out = tmp_path / "output.bedGraph"

    merge_bins(
        str(fil_in),
        str(fil_out),
        decimal_places=decimal_places,
        eps=eps,
        skp_pfx=skp_pfx,
    )

    return fil_out.read_text(encoding="utf-8")


def test_format_value_rounds_and_suppresses_negative_zero() -> None:
    assert format_value("-0.0001", 2, -0.0001) == "0.00"
    assert format_value("1.234", 2, 1.234) == "1.23"


def test_format_value_canonicalizes_nonfinite_tokens() -> None:
    assert format_value("NaN", None, None) == "nan"
    assert format_value("+Inf", None, None) == "inf"
    assert format_value("-INF", None, None) == "-inf"


def test_format_value_preserves_unrounded_tokens() -> None:
    assert format_value("0.500000", None, 0.5) == "0.500000"


def test_validate_merge_options_accepts_none_and_zero_values() -> None:
    validate_merge_options(None, None)
    validate_merge_options(0, 0.0)


def test_merge_bins_exact_text_mode(tmp_path: Path) -> None:
    observed = run_merge(
        tmp_path,
        "exact_text.bedGraph",
    )

    assert observed == ("chrI\t0\t20\t1\nchrI\t20\t40\t1.0\n")


def test_merge_bins_flushes_on_midstream_header(tmp_path: Path) -> None:
    observed = run_merge(
        tmp_path,
        "midstream_header.bedGraph",
    )

    assert observed == (
        "chrI\t0\t20\t2\ntrack type=bedGraph\nchrI\t20\t40\t2\n"
    )


def test_merge_bins_rounding_mode(tmp_path: Path) -> None:
    observed = run_merge(
        tmp_path,
        "rounding.bedGraph",
        decimal_places=2,
    )

    assert observed == ("chrI\t0\t20\t1.00\nchrI\t20\t30\t1.01\n")


def test_merge_bins_epsilon_mode(tmp_path: Path) -> None:
    observed = run_merge(
        tmp_path,
        "epsilon.bedGraph",
        eps=0.001,
    )

    assert observed == ("chrI\t0\t20\t1.0000\nchrI\t20\t30\t1.0020\n")


def test_merge_bins_rounding_plus_epsilon_mode(tmp_path: Path) -> None:
    observed = run_merge(
        tmp_path,
        "rounding_epsilon.bedGraph",
        decimal_places=2,
        eps=0.01,
    )

    assert observed == ("chrI\t0\t20\t1.00\nchrI\t20\t30\t1.02\n")


def test_merge_bins_never_merges_nonfinite_values(tmp_path: Path) -> None:
    observed = run_merge(
        tmp_path,
        "nonfinite.bedGraph",
    )

    assert observed == (
        "chrI\t0\t10\tnan\n"
        "chrI\t10\t20\tnan\n"
        "chrI\t20\t30\tinf\n"
        "chrI\t30\t40\tinf\n"
    )


def test_merge_bins_rejects_malformed_rows(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="Malformed bedGraph line"):
        run_merge(tmp_path, "malformed.bedGraph")


def test_merge_bins_rejects_nonpositive_intervals(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="Non-positive interval width"):
        run_merge(tmp_path, "zero_width.bedGraph")


def test_merge_bins_rejects_negative_rounding_precision(
    tmp_path: Path,
) -> None:
    with pytest.raises(ValueError, match="'--dp' must be >= 0"):
        run_merge(
            tmp_path,
            "one_row.bedGraph",
            decimal_places=-1,
        )


def test_merge_bins_rejects_negative_epsilon(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="'--eps' must be >= 0"):
        run_merge(tmp_path, "one_row.bedGraph", eps=-0.1)
