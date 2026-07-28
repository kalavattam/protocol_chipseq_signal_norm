#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_sum.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


from pathlib import Path

from protocol_chipseq_signal_norm.cli.sum_bdg import sum_bdg
from protocol_chipseq_signal_norm.utilities.utils_io import DEF_SKP_PFX


def test_sum_bdg_sums_finite_values(tmp_path: Path) -> None:
    path = tmp_path / "input.bdg"
    path.write_text(
        "track name=x\nchrI 0 10 1\nchrI 10 20 2\nchrI 20 30 nan\n",
        encoding="utf-8",
    )

    assert sum_bdg(str(path), weight=False, skp_pfx=DEF_SKP_PFX) == 3.0


def test_sum_bdg_can_weight_by_interval_width(tmp_path: Path) -> None:
    path = tmp_path / "input.bdg"
    path.write_text("chrI 0 10 1\nchrI 10 30 2\n", encoding="utf-8")

    assert sum_bdg(str(path), weight=True, skp_pfx=DEF_SKP_PFX) == 50.0
