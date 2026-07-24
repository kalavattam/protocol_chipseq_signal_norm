#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: utils_format.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


"""
Numeric formatting helpers.

Notes
-----
Helpers format finite values with bounded decimal precision and strip
non-informative trailing zeros.
"""

from __future__ import annotations

import sys

assert sys.version_info >= (3, 11), "Python >= 3.11 required."


def format_value(x: float, dp: int) -> str:
    """
    Round and format a numeric value with no non-informative trailing zeros.

    Parameters
    ----------
        x : float
            Numeric value to format.
        dp : int
            Maximum number of decimal places to retain.

    Returns
    -------
        str
            Rounded value with trailing zeros and any trailing decimal point
            stripped. Negative zero is returned as '0'.

    Raises
    ------
        None.

    Notes
    -----
        This helper expects callers to validate that 'dp' is non-negative.
    """
    y = round(x, dp)
    if y == 0.0:
        y = 0.0
    s = f"{y:.{dp}f}"
    if "." in s:
        s = s.rstrip("0").rstrip(".")
    if s == "-0":
        s = "0"
    return s
