#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_ratio.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


import math

import pytest

from protocol_chipseq_signal_norm.cli.compute_signal_ratio import (
    calc_rat_bin,
    comp_sig_rat,
    parse_pair,
)


def test_parse_pair_accepts_single_or_pair_values():
    assert parse_pair("2", 1.0) == (2.0, 1.0)
    assert parse_pair("2:3.5", 1.0) == (2.0, 3.5)

    with pytest.raises(Exception, match="Expected"):
        parse_pair("2:3:4", 1.0)


def test_calc_rat_bin_handles_scaling_pseudocounts_and_log_reciprocal():
    assert calc_rat_bin(
        4.0, 2.0, 2.0, 1.0, 1.0, 0.0, None, False, False, None
    ) == 4.5
    assert calc_rat_bin(
        4.0, 2.0, None, None, None, None, None, True, True, None
    ) == -1.0


def test_calc_rat_bin_zero_handling():
    assert calc_rat_bin(0, 0, None, None, None, None, None, False, False,
                        "pre_scale") is None
    assert math.isnan(calc_rat_bin(
        1, 0, None, None, None, None, None, False, False, None
    ))


def test_comp_sig_rat_writes_small_bedgraph(tmp_path):
    fil_a = tmp_path / "a.bdg"
    fil_b = tmp_path / "b.bdg"
    fil_out = tmp_path / "ratio.bdg"
    fil_a.write_text("chrI 0 10 4\nchrI 10 20 0\n", encoding="utf-8")
    fil_b.write_text("chrI 0 10 2\nchrI 10 20 0\n", encoding="utf-8")

    comp_sig_rat(
        str(fil_a),
        str(fil_b),
        str(fil_out),
        scl_A=1.0,
        scl_B=1.0,
        psc_A=0.0,
        psc_B=0.0,
        dep_min=None,
        dp=3,
        log2=False,
        recip=False,
        skip_00="pre_scale",
        eps=0.0,
        track=False,
        drp_nan=False,
    )

    assert fil_out.read_text(encoding="utf-8") == "chrI\t0\t10\t2\n"


def test_comp_sig_rat_strict_bins_rejects_end_mismatch(tmp_path):
    fil_a = tmp_path / "a.bdg"
    fil_b = tmp_path / "b.bdg"
    fil_out = tmp_path / "ratio.bdg"
    fil_a.write_text("chrI 0 10 4\nchrI 10 20 2\n", encoding="utf-8")
    fil_b.write_text("chrI 0 10 2\nchrI 10 21 1\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Mismatched bedGraph grids"):
        comp_sig_rat(
            str(fil_a),
            str(fil_b),
            str(fil_out),
            scl_A=1.0,
            scl_B=1.0,
            psc_A=0.0,
            psc_B=0.0,
            dep_min=None,
            dp=3,
            log2=False,
            recip=False,
            skip_00=None,
            eps=0.0,
            track=False,
            drp_nan=False,
            strict_bins=True,
        )


def test_comp_sig_rat_chr_sizes_rejects_out_of_bounds_input(tmp_path):
    fil_a = tmp_path / "a.bdg"
    fil_b = tmp_path / "b.bdg"
    fil_out = tmp_path / "ratio.bdg"
    fil_a.write_text("chrI 70 81 4\n", encoding="utf-8")
    fil_b.write_text("chrI 70 80 2\n", encoding="utf-8")

    with pytest.raises(ValueError, match="extends beyond"):
        comp_sig_rat(
            str(fil_a),
            str(fil_b),
            str(fil_out),
            scl_A=1.0,
            scl_B=1.0,
            psc_A=0.0,
            psc_B=0.0,
            dep_min=None,
            dp=3,
            log2=False,
            recip=False,
            skip_00=None,
            eps=0.0,
            track=False,
            drp_nan=False,
            chrom_sizes={"chrI": 80},
        )


def test_comp_sig_rat_rejects_dash_io(tmp_path):
    fil_a = tmp_path / "a.bdg"
    fil_b = tmp_path / "b.bdg"
    fil_out = tmp_path / "ratio.bdg"
    fil_a.write_text("chrI 0 10 4\n", encoding="utf-8")
    fil_b.write_text("chrI 0 10 2\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Dash input/output"):
        comp_sig_rat(
            "-",
            str(fil_b),
            str(fil_out),
            scl_A=1.0,
            scl_B=1.0,
            psc_A=0.0,
            psc_B=0.0,
            dep_min=None,
            dp=3,
            log2=False,
            recip=False,
            skip_00=None,
            eps=0.0,
            track=False,
            drp_nan=False,
        )
