#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_stabilizer.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


import math
from pathlib import Path

import pytest

from protocol_chipseq_signal_norm.utilities.utils_stabilizer import (
    compute_stats_robust,
    determine_coef_eff,
    iter_vals_bdg,
    median_sorted,
    pick_stabilizer,
)


def test_compute_stats_robust_ignores_nonfinite_values() -> None:
    stats = compute_stats_robust([1.0, float("nan"), 3.0, float("inf")])

    assert stats == {"n": 2, "median": 2.0, "mean": 2.0}


def test_determine_coef_eff_and_median_sorted() -> None:
    assert determine_coef_eff("frc_mdn_nz", None) == 0.01
    assert determine_coef_eff("min_nz", None) == 1.0
    assert determine_coef_eff("qntl_nz", None) is None
    assert median_sorted([1.0, 2.0, 10.0, 20.0]) == 6.0


def test_pick_stabilizer_quantile_and_floor() -> None:
    assert pick_stabilizer([1.0, 2.0, 10.0], "qntl_nz", qntl_pct=50) == 2.0
    assert pick_stabilizer([1.0, 2.0], "min_nz", floor=5.0) == 5.0
    assert math.isnan(pick_stabilizer([float("nan")], "qntl_nz"))


def test_pick_stabilizer_rejects_bad_quantile() -> None:
    with pytest.raises(ValueError, match="qntl_pct"):
        pick_stabilizer([1.0], "qntl_nz", qntl_pct=101)


def test_iter_vals_bdg_filters_by_positive_policy(tmp_path: Path) -> None:
    path = tmp_path / "values.bdg"
    path.write_text(
        "chrI 0 10 0\nchrI 10 20 0.1\nchrI 20 30 2\nchrI 30 40 nan\n",
        encoding="utf-8",
    )

    assert list(
        iter_vals_bdg(
            str(path),
            eps=0.1,
            mode_nz="closed",
            nz_policy="pos",
        ),
    ) == [2.0]
