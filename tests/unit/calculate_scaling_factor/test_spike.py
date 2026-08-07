#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_spike.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


import pytest

from protocol_chipseq_signal_norm.cli.calculate_scaling_factor_spike import (
    calculate_scaling_factors,
    normalize_coef,
    round_value,
)


def test_normalize_coef_accepts_documented_aliases() -> None:
    assert normalize_coef("chiprx-ratio") == "chiprx_alpha_ratio"
    assert normalize_coef("bio_protocol") == "fractional"
    assert normalize_coef("rxi") == "rxinput_alpha"


def test_calculate_scaling_factors_computes_core_coefficients() -> None:
    vals = calculate_scaling_factors(100, 10, 100, 20)

    assert vals["fractional"] == pytest.approx(11 / 6)
    assert vals["chiprx_alpha_ip"] == pytest.approx(100000.0)
    assert vals["chiprx_alpha_in"] == pytest.approx(50000.0)
    assert vals["chiprx_alpha_ratio"] == pytest.approx(2.0)
    assert vals["rxinput_alpha"] == pytest.approx(50000 / 3)


def test_calculate_scaling_factors_rejects_invalid_inputs() -> None:
    with pytest.raises(ValueError, match="must be >= 0"):
        calculate_scaling_factors(-1, 10, 100, 20)

    with pytest.raises(ZeroDivisionError, match="spike_ip"):
        calculate_scaling_factors(100, 0, 100, 20)


def test_round_value_collapses_negative_zero() -> None:
    assert round_value(-0.0001, 2) == 0.0
