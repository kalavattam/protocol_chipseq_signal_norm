#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_format.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


from protocol_chipseq_signal_norm.utilities.utils_format import format_value


def test_format_value_strips_trailing_zeros() -> None:
    assert format_value(1.2300, 4) == "1.23"
    assert format_value(2.0, 4) == "2"


def test_format_value_collapses_negative_zero() -> None:
    assert format_value(-0.0001, 2) == "0"
