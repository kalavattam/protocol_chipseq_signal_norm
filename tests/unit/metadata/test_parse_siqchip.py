#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_parse_siqchip.py
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


"""
Test siQ-ChIP metadata parsing helpers.
"""

import pytest

from protocol_chipseq_signal_norm.cli.parse_metadata_siqchip import (
    collect_outputs,
    find_row_matching,
    normalize_input_map,
    parse_args,
    parse_filename,
    validate_id,
)


def minimal_configuration() -> dict[str, object]:
    """
    Return the minimal configuration needed by metadata helper tests.

    Returns
    -------
    configuration : dict[str, object]
        Nested filename, matching, and calculator-input settings.
    """

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
            },
        },
    }


def test_parse_filename_strips_extension_and_suffix() -> None:
    configuration = minimal_configuration()

    assert parse_filename("IP_1.sc.bam", configuration) == {
        "sample": "IP",
        "rep": "1",
    }


def test_find_row_matching_returns_unique_match() -> None:
    configuration = minimal_configuration()
    rows = [
        {"sample": "IP", "rep": "1", "mass_ip": "2", "mass_in": "4"},
        {"sample": "IP", "rep": "2", "mass_ip": "3", "mass_in": "4"},
    ]

    matched = find_row_matching(
        rows,
        {"sample": "IP", "rep": "2"},
        configuration,
    )

    assert matched is rows[1]


def test_collect_outputs_requires_paired_optional_fields() -> None:
    configuration = minimal_configuration()

    with pytest.raises(ValueError, match="provided together"):
        collect_outputs(
            {
                "mass_ip": "2",
                "mass_in": "4",
                "lib_vol_ip": "10",
                "lib_vol_in": "",
            },
            configuration,
        )


def test_normalize_input_map_and_validate_id_reject_bad_names() -> None:
    assert normalize_input_map(["mass_ip"], "required") == {
        "mass_ip": "mass_ip",
    }

    with pytest.raises(ValueError, match="shell-safe"):
        validate_id("bad-name", "field")


# Both options are required unless '--validate_cfg' is given, which
# 'parse_args' enforces itself rather than through 'required=True'.
CONDITIONAL_REQUIRED = (
    pytest.param("--alignment", ("--tbl_met", "metadata.tsv"), id="alignment"),
    pytest.param("--tbl_met", ("--alignment", "sample.bam"), id="tbl_met"),
)


@pytest.mark.parametrize(("omitted", "supplied"), CONDITIONAL_REQUIRED)
def test_parse_args_reports_a_conditionally_required_option(
    omitted: str,
    supplied: tuple[str, ...],
    capsys: pytest.CaptureFixture[str],
) -> None:
    """
    Omitting one reports a usage error rather than raising 'AttributeError'.

    'CapArgumentParser' sets 'argument_default=argparse.SUPPRESS', so until
    both options declared 'default=None' the guard reading them raised before
    it could call 'parser.error', and this message was unreachable.
    """

    with pytest.raises(SystemExit) as error:
        parse_args(["--cfg", "parser.yml", *supplied])

    note = capsys.readouterr().err
    expected = f"'{omitted}' is required unless '--validate_cfg' is supplied."

    assert error.value.code == 2
    assert expected in " ".join(note.split())
