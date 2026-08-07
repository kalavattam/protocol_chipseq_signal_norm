#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_merge_bins.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in design, development,
# and documentation, with all output reviewed, edited, and approved by the
# author.
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


def run_merge(
    tmp_path: Path,
    text: str,
    *,
    decimal_places: int | None = None,
    eps: float | None = None,
    skp_pfx: tuple[str, ...] = DEF_SKP_PFX,
) -> str:
    fil_in = tmp_path / "input.bdg"
    fil_out = tmp_path / "output.bdg"

    fil_in.write_text(text, encoding="utf-8")
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
        "\n".join(
            [
                "chrI 0 10 1",
                "chrI 10 20 1",
                "chrI 20 30 1.0",
                "chrI 30 40 1.0",
                "",
            ],
        ),
    )

    assert observed == ("chrI\t0\t20\t1\nchrI\t20\t40\t1.0\n")


def test_merge_bins_flushes_on_midstream_header(tmp_path: Path) -> None:
    observed = run_merge(
        tmp_path,
        "\n".join(
            [
                "chrI 0 10 2",
                "chrI 10 20 2",
                "track type=bedGraph",
                "chrI 20 30 2",
                "chrI 30 40 2",
                "",
            ],
        ),
    )

    assert observed == (
        "chrI\t0\t20\t2\ntrack type=bedGraph\nchrI\t20\t40\t2\n"
    )


def test_merge_bins_rounding_mode(tmp_path: Path) -> None:
    observed = run_merge(
        tmp_path,
        "\n".join(
            [
                "chrI 0 10 1.001",
                "chrI 10 20 1.002",
                "chrI 20 30 1.014",
                "",
            ],
        ),
        decimal_places=2,
    )

    assert observed == ("chrI\t0\t20\t1.00\nchrI\t20\t30\t1.01\n")


def test_merge_bins_epsilon_mode(tmp_path: Path) -> None:
    observed = run_merge(
        tmp_path,
        "\n".join(
            [
                "chrI 0 10 1.0000",
                "chrI 10 20 1.0004",
                "chrI 20 30 1.0020",
                "",
            ],
        ),
        eps=0.001,
    )

    assert observed == ("chrI\t0\t20\t1.0000\nchrI\t20\t30\t1.0020\n")


def test_merge_bins_rounding_plus_epsilon_mode(tmp_path: Path) -> None:
    observed = run_merge(
        tmp_path,
        "\n".join(
            [
                "chrI 0 10 1.001",
                "chrI 10 20 1.009",
                "chrI 20 30 1.021",
                "",
            ],
        ),
        decimal_places=2,
        eps=0.01,
    )

    assert observed == ("chrI\t0\t20\t1.00\nchrI\t20\t30\t1.02\n")


def test_merge_bins_never_merges_nonfinite_values(tmp_path: Path) -> None:
    observed = run_merge(
        tmp_path,
        "\n".join(
            [
                "chrI 0 10 NaN",
                "chrI 10 20 nan",
                "chrI 20 30 inf",
                "chrI 30 40 inf",
                "",
            ],
        ),
    )

    assert observed == (
        "chrI\t0\t10\tnan\n"
        "chrI\t10\t20\tnan\n"
        "chrI\t20\t30\tinf\n"
        "chrI\t30\t40\tinf\n"
    )


def test_merge_bins_rejects_malformed_rows(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="Malformed bedGraph line"):
        run_merge(tmp_path, "chrI 0 10\n")


def test_merge_bins_rejects_nonpositive_intervals(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="Non-positive interval width"):
        run_merge(tmp_path, "chrI 10 10 1\n")


def test_merge_bins_rejects_negative_rounding_precision(
    tmp_path: Path,
) -> None:
    with pytest.raises(ValueError, match="'--dp' must be >= 0"):
        run_merge(
            tmp_path,
            "chrI 0 10 1\n",
            decimal_places=-1,
        )


def test_merge_bins_rejects_negative_epsilon(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="'--eps' must be >= 0"):
        run_merge(tmp_path, "chrI 0 10 1\n", eps=-0.1)
