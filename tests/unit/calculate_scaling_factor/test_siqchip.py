#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_siqchip.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.6);
# - Anthropic Claude Code (Opus 5).
#
# Distributed under the MIT license.


import pytest

from protocol_chipseq_signal_norm.cli.calculate_scaling_factor_siqchip import (
    calculate_alpha,
    parse_args,
)

# The six mass, volume, and length options are required; the depths are not.
REQUIRED_SIQCHIP = (
    "--mass_ip",
    "1",
    "--mass_in",
    "1",
    "--vol_all",
    "1",
    "--vol_in",
    "1",
    "--len_ip",
    "1",
    "--len_in",
    "1",
)


def test_parse_args_binds_unsupplied_depths_to_none() -> None:
    """
    '--dep_ip' and '--dep_out' are read as 'getattr(args, ..., None)'.

    'CapArgumentParser' sets 'argument_default=argparse.SUPPRESS', so without
    an explicit 'default=None' an unsupplied depth is absent from the namespace
    and plain attribute access raises instead of reporting None.
    """

    args = parse_args(list(REQUIRED_SIQCHIP))

    assert args.dep_ip is None
    assert args.dep_in is None


def test_calculate_alpha_equations() -> None:
    assert calculate_alpha("5", 2, 4, 100, 10, 50, 100, 20, 40) == 0.2
    assert calculate_alpha("5nd", 2, 4, 100, 10, None, None, 20, 40) == 0.1
    assert calculate_alpha(
        "6nd", 2, 4, 100, 10, None, None, 20, 40
    ) == pytest.approx(1 / 9)


def test_calculate_alpha_applies_library_volume_ratio() -> None:
    observed_alpha = calculate_alpha(
        "5nd",
        2,
        4,
        100,
        10,
        None,
        None,
        20,
        40,
        lib_vol_ip=5,
        lib_vol_in=10,
    )

    assert observed_alpha == 0.2


def test_calculate_alpha_rejects_bad_equation_and_volume() -> None:
    with pytest.raises(ValueError, match="Unsupported equation"):
        calculate_alpha("bad", 1, 1, 10, 1, 1, 1, 1, 1)

    with pytest.raises(ValueError, match="vol_all"):
        calculate_alpha("6", 1, 1, 10, 10, 1, 1, 1, 1)
