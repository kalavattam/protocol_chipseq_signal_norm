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
    as_iter,
    check_cmp,
    check_exists,
    check_parse_fil_out,
    pair_val_thresh,
)


def test_as_iter_treats_strings_and_paths_as_scalars():
    assert as_iter("abc") == ("abc",)
    assert as_iter(Path("x")) == (Path("x"),)
    assert as_iter([1, 2]) == (1, 2)


def test_pair_val_thresh_broadcasts_scalar_threshold():
    assert pair_val_thresh([1, 2, 3], 0) == [(1, 0), (2, 0), (3, 0)]


def test_check_cmp_validates_scalars_and_iterables():
    check_cmp([1, 2, None], "gt", 0, "value")

    with pytest.raises(ValueError, match="'--value' must be > 0"):
        check_cmp(0, "gt", 0, "value")


def test_check_parse_fil_out_preserves_extension_and_gzip():
    assert check_parse_fil_out("signal.bedGraph.gz") == (
        "signal.bedGraph.gz",
        "bedGraph",
        True,
    )

    with pytest.raises(ValueError, match="Invalid extension"):
        check_parse_fil_out("signal.txt")


def test_check_exists_accepts_files_and_rejects_missing(tmp_path):
    path = tmp_path / "input.txt"
    path.write_text("x", encoding="utf-8")

    check_exists(path, kind="file")
    with pytest.raises(FileNotFoundError):
        check_exists(tmp_path / "missing.txt", kind="file")
