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


def format_value(value: float, decimal_places: int) -> str:
    """
    Round and format a numeric value with no non-informative trailing zeros.

    Parameters
    ----------
    value : float
        Numeric value to format.
    decimal_places : int
        Maximum number of decimal places to retain.

    Returns
    -------
    text : str
        Rounded value with trailing zeros and any trailing decimal point
        stripped. Negative zero is returned as '0'.

    Notes
    -----
    This helper expects callers to validate that 'decimal_places' is
    non-negative.
    """

    rounded_value = round(value, decimal_places)

    if rounded_value == 0.0:
        rounded_value = 0.0

    rendered = f"{rounded_value:.{decimal_places}f}"

    if "." in rendered:
        rendered = rendered.rstrip("0").rstrip(".")

    if rendered == "-0":
        rendered = "0"

    return rendered
