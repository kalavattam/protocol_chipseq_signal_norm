#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: utils_format.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.

"""
Script
------
utils_format.py


Description
-----------
Shared formatting helpers for numeric command-line and table output.


Functions
---------
format_value()
    Round and format a numeric value with at most a requested number of
    decimal places, stripping non-informative trailing zeros and any trailing
    decimal point.


Example
-------
'''python
from scripts.functions.utils_format import format_value

format_value(1.0, 24)      # "1"
format_value(0.125000, 6)  # "0.125"
format_value(-0.0, 6)      # "0"
'''


Notes
-----
- 'format_value()' treats 'rnd' as a maximum printed precision.
- The helper does not validate whether 'rnd' is non-negative; callers should
  validate CLI input before calling it.
"""

from __future__ import annotations
import sys

assert sys.version_info >= (3, 10), "Python >= 3.10 required."


def format_value(x: float, rnd: int) -> str:
    """
    Round and format a numeric value with no non-informative trailing zeros.

    Args:
        x : float
            Numeric value to format.
        rnd : int
            Maximum number of decimal places to retain.

    Returns:
        str
            Rounded value with trailing zeros and any trailing decimal point
            stripped. Negative zero is returned as '0'.

    Raises:
        None.

    Notes:
        This helper expects callers to validate that 'rnd' is non-negative.
    """
    y = round(x, rnd)
    if y == 0.0:
        y = 0.0
    s = f"{y:.{rnd}f}"
    if "." in s:
        s = s.rstrip("0").rstrip(".")
    if s == "-0":
        s = "0"
    return s
