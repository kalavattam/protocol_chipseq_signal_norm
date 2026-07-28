#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_check.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


from pathlib import Path

import pytest

from protocol_chipseq_signal_norm.utilities.utils_check import (
    as_tuple,
    check_exists,
    pair_values_and_thresholds,
    validate_comparison,
    validate_output_path,
)


def test_as_tuple_treats_strings_and_paths_as_scalars() -> None:
    assert as_tuple("abc") == ("abc",)
    assert as_tuple(Path("x")) == (Path("x"),)
    assert as_tuple([1, 2]) == (1, 2)


def test_pair_values_and_thresholds_broadcasts_scalar_threshold() -> None:
    assert pair_values_and_thresholds([1, 2, 3], 0) == [(1, 0), (2, 0), (3, 0)]


def test_validate_comparison_validates_scalars_and_iterables() -> None:
    validate_comparison([1, 2, None], "gt", 0, "value")

    with pytest.raises(ValueError, match="'--value' must be > 0"):
        validate_comparison(0, "gt", 0, "value")


def test_validate_output_path_preserves_extension_and_gzip() -> None:
    assert validate_output_path("signal.bedGraph.gz") == (
        "signal.bedGraph.gz",
        "bedGraph",
        True,
    )

    with pytest.raises(ValueError, match="Invalid extension"):
        validate_output_path("signal.txt")


def test_check_exists_accepts_files_and_rejects_missing(
    tmp_path: Path,
) -> None:
    path = tmp_path / "input.txt"
    path.write_text("x", encoding="utf-8")

    check_exists(path, kind="file")

    with pytest.raises(FileNotFoundError):
        check_exists(tmp_path / "missing.txt", kind="file")
