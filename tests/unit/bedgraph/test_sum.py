#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_sum.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.6);
# - Anthropic Claude Code (Opus 5.5).
#
# Distributed under the MIT license.


from pathlib import Path

from protocol_chipseq_signal_norm.cli.sum_bdg import sum_bdg
from protocol_chipseq_signal_norm.utilities.utils_io import DEF_SKP_PFX

ROOT = Path(__file__).resolve().parents[3]
BEDGRAPH_UNIT = ROOT / "tests" / "fixtures" / "bedgraph"


def test_sum_bdg_sums_finite_values() -> None:
    path = BEDGRAPH_UNIT / "sum" / "header_and_nan.bedGraph"

    assert sum_bdg(str(path), weight=False, skp_pfx=DEF_SKP_PFX) == 3.0


def test_sum_bdg_can_weight_by_interval_width() -> None:
    path = BEDGRAPH_UNIT / "sum" / "two_widths.bedGraph"

    assert sum_bdg(str(path), weight=True, skp_pfx=DEF_SKP_PFX) == 50.0
