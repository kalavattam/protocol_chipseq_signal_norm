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
    COEF_ORDER,
    calculate_scaling_factors,
    normalize_coef,
    round_value,
)


def test_normalize_coef_accepts_documented_aliases() -> None:
    assert normalize_coef("chiprx-ratio") == "chiprx_alpha_ratio"
    assert normalize_coef("bio_protocol") == "fractional"
    assert normalize_coef("rxi") == "rxinput_alpha"
    assert normalize_coef("main_spike_ratio") == "main_per_spike"
    assert normalize_coef("main-per-spike") == "main_per_spike"
    assert normalize_coef("mps") == "main_per_spike"
    assert normalize_coef("m") == "main_per_spike"


def test_coef_order_pins_the_all_output_order() -> None:
    """
    '--coef all' emits one row per name in 'COEF_ORDER', in that order, so the
    tuple is a consumer-facing contract. Pinning it makes a reordering a
    deliberate edit rather than a silent change for anything reading those rows
    by position.
    """

    assert COEF_ORDER == (
        "fractional",
        "main_per_spike",
        "chiprx_alpha_ip",
        "chiprx_alpha_in",
        "chiprx_alpha_ratio",
        "rxinput_alpha",
    )


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


def test_main_per_spike_is_yeast_material_per_reference() -> None:
    """
    'main_per_spike' divides main by spike counts within each sample, then
    takes the IP-over-input ratio.
    """

    vals = calculate_scaling_factors(100, 10, 100, 20)

    assert vals["main_per_spike"] == pytest.approx((100 / 10) / (100 / 20))


def test_main_per_spike_restores_what_norm_divided_out() -> None:
    """
    A 'norm' track divides by its own main count, so the coefficient that
    cancels that division is 'chiprx_alpha_ratio' times the main-count ratio.
    This identity is why 'main_per_spike' is the exact partner of 'norm'.
    """

    main_ip, spike_ip, main_in, spike_in = 16568363, 284655, 18576743, 364252
    vals = calculate_scaling_factors(main_ip, spike_ip, main_in, spike_in)

    assert vals["main_per_spike"] == pytest.approx(
        vals["chiprx_alpha_ratio"] * (main_ip / main_in), rel=1e-15
    )


def test_fractional_approaches_main_per_spike_as_spike_share_falls() -> None:
    """
    'fractional' divides by the whole library, 'main_per_spike' by the main
    counts alone, so they differ by the ratio of main shares. The gap shrinks
    with the spike-in share and vanishes in the limit.
    """

    main_ip, spike_ip, main_in, spike_in = 16568363, 284655, 18576743, 364252
    vals = calculate_scaling_factors(main_ip, spike_ip, main_in, spike_in)

    share_ip = spike_ip / (main_ip + spike_ip)
    share_in = spike_in / (main_in + spike_in)

    assert vals["fractional"] * (1 - share_ip) / (1 - share_in) == (
        pytest.approx(vals["main_per_spike"], rel=1e-15)
    )

    #  A thousandfold smaller spike-in share collapses the gap.
    lean = calculate_scaling_factors(main_ip, 285, main_in, 364)
    assert lean["fractional"] == pytest.approx(
        lean["main_per_spike"], rel=1e-4
    )


def test_main_per_spike_rejects_a_zero_main_input() -> None:
    """
    'main_per_spike' is the only coefficient that divides by 'main_in', so it
    is the only one that has to guard against it being zero.
    """

    with pytest.raises(ZeroDivisionError, match="main_in"):
        calculate_scaling_factors(100, 10, 0, 20, required=("main_per_spike",))
