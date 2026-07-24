#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_pseudo.py
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

from protocol_chipseq_signal_norm.cli.compute_pseudo import combine_pseudo_sym


def test_combine_pseudo_sym_returns_unmodified_for_none_mode():
    assert combine_pseudo_sym(1.0, 2.0, "none") == (1.0, 2.0)


def test_combine_pseudo_sym_applies_symmetric_modes():
    assert combine_pseudo_sym(1.0, 3.0, "max") == (3.0, 3.0)
    assert combine_pseudo_sym(1.0, 3.0, "arith") == (2.0, 2.0)
    assert combine_pseudo_sym(1.0, 4.0, "geom") == (2.0, 2.0)
    assert combine_pseudo_sym(1.0, 4.0, "harm") == (1.6, 1.6)


def test_combine_pseudo_sym_mirrors_single_finite_value(capsys):
    assert combine_pseudo_sym(math.nan, 2.0, "max") == (2.0, 2.0)
    assert "nonfinite" in capsys.readouterr().err


def test_combine_pseudo_sym_rejects_unknown_mode():
    with pytest.raises(ValueError, match="Unknown"):
        combine_pseudo_sym(1.0, 2.0, "bad")
