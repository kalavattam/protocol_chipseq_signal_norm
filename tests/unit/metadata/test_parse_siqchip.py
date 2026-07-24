#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_parse_siqchip.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


import pytest

pytest.importorskip("yaml")

from protocol_chipseq_signal_norm.cli.parse_metadata_siqchip import (
    collect_outputs,
    find_row_matching,
    normalize_input_map,
    parse_filename,
    validate_id,
)


def cfg_minimal():
    return {
        "filename": {
            "delimiter": "_",
            "fields": ["sample", "rep"],
            "strip_extensions": [".bam"],
            "strip_suffixes": [".sc"],
        },
        "matching": {"fields": ["sample", "rep"]},
        "field_to_column": {},
        "calculator_inputs": {
            "siqchip": {
                "required": {"mass_ip": "mass_ip", "mass_in": "mass_in"},
                "optional": {
                    "lib_vol_ip": "lib_vol_ip",
                    "lib_vol_in": "lib_vol_in",
                },
            }
        },
    }


def test_parse_filename_strips_extension_and_suffix():
    cfg = cfg_minimal()

    assert parse_filename("IP_1.sc.bam", cfg) == {"sample": "IP", "rep": "1"}


def test_find_row_matching_returns_unique_match():
    cfg = cfg_minimal()
    rows = [
        {"sample": "IP", "rep": "1", "mass_ip": "2", "mass_in": "4"},
        {"sample": "IP", "rep": "2", "mass_ip": "3", "mass_in": "4"},
    ]

    assert find_row_matching(rows, {"sample": "IP", "rep": "2"}, cfg) is rows[1]


def test_collect_outputs_requires_paired_optional_fields():
    cfg = cfg_minimal()

    with pytest.raises(ValueError, match="provided together"):
        collect_outputs(
            {
                "mass_ip": "2",
                "mass_in": "4",
                "lib_vol_ip": "10",
                "lib_vol_in": "",
            },
            cfg,
        )


def test_normalize_input_map_and_validate_id_reject_bad_names():
    assert normalize_input_map(["mass_ip"], "required") == {
        "mass_ip": "mass_ip"
    }

    with pytest.raises(ValueError, match="shell-safe"):
        validate_id("bad-name", "field")
