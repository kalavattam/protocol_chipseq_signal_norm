#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_relativize.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


import pytest

from protocol_chipseq_signal_norm.cli.relativize_scaling_factors import (
    determine_scaling_column,
    format_scaled,
    relativize,
)


def test_format_scaled_strips_trailing_zeros() -> None:
    assert format_scaled(0.5000001, 3) == "0.5"


def test_determine_scaling_column_prefers_siq() -> None:
    assert determine_scaling_column(["sample", "siq", "spike"]) == "siq"


def test_relativize_scales_ip_and_optionally_input() -> None:
    header = ["sample", "spike"]
    rows = [
        {"sample": "IP_A", "spike": "2"},
        {"sample": "IP_B", "spike": "4"},
        {"sample": "in_A", "spike": "1"},
    ]

    out_header, out_rows = relativize(header, rows, "spike", True, 3)

    assert out_header == ["sample", "spike", "scaled"]
    assert [row["scaled"] for row in out_rows] == ["0.5", "1", "0.25"]


def test_relativize_rejects_missing_sample_column() -> None:
    with pytest.raises(ValueError, match="sample"):
        relativize(["spike"], [{"spike": "1"}], "spike", False, 3)
