#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_help_option_order.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.

"""
Test registered semantic option-order realizations.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest
from dev.audit.help_option_order import (
    python_reporting_order,
    render_help,
    validate_order,
)

ROOT = Path(__file__).resolve().parents[3]
CONFIG = ROOT / "dev/config/help_contracts.json"


def data() -> dict:
    return json.loads(CONFIG.read_text(encoding="utf-8"))


def test_accepted_pilot_realizations_match_source_projections() -> None:
    assert validate_order(ROOT, data()) == []


def test_unreviewed_roles_remain_visible() -> None:
    payload = copy.deepcopy(data())
    payload["option_realizations"][0]["roles_reviewed"] = False

    assert "HELP.OPTION.ORDER.ROLE_UNREVIEWED" in {
        item["rule_id"] for item in validate_order(ROOT, payload)
    }


def emitted(payload: dict) -> set[str]:
    """
    Return exact diagnostics for one rejected registered realization.
    """

    return {item["rule_id"] for item in validate_order(ROOT, payload)}


def test_category_precedence_uses_explicit_roles() -> None:
    payload = copy.deepcopy(data())
    record = payload["option_realizations"][0]
    record["roles"]["dp"] = "execution_environment_resources"

    assert "HELP.OPTION.ORDER.CATEGORY" in emitted(payload)


def test_immediate_adjacency_is_stronger_than_same_group() -> None:
    payload = copy.deepcopy(data())
    record = payload["option_realizations"][0]
    record["relationships"][0]["members"] = ["fil_in", "ref_fa"]

    assert "HELP.OPTION.ORDER.ADJACENCY" in emitted(payload)

    record["relationships"][0]["kind"] = "same_group"

    assert "HELP.OPTION.ORDER.ADJACENCY" not in emitted(payload)
    assert "HELP.OPTION.ORDER.GROUP" not in emitted(payload)


def test_contiguous_group_rejects_an_intervening_role() -> None:
    payload = copy.deepcopy(data())
    record = payload["option_realizations"][0]
    record["relationships"][1]["members"] = ["method", "floor", "dp"]

    assert "HELP.OPTION.ORDER.CONTIGUITY" in emitted(payload)


def test_same_group_requires_one_reviewed_role() -> None:
    payload = copy.deepcopy(data())
    record = payload["option_realizations"][0]
    record["relationships"][2]["members"] = ["siz_bin", "dp"]

    assert "HELP.OPTION.ORDER.GROUP" in emitted(payload)


def test_four_registered_surface_orders_require_exact_parity() -> None:
    payload = copy.deepcopy(data())
    record = payload["option_realizations"][0]
    record["surface_orders"]["usage"][3:5] = ["infmt", "fil_in"]

    assert "HELP.OPTION.ORDER.SURFACE_PARITY" in emitted(payload)


def test_usage_rows_preserve_registered_group_boundaries_and_order() -> None:
    payload = data()
    command = payload["option_realizations"][0]["adapters"]["help_command"]
    rendered = render_help(ROOT, command)

    merged = rendered.replace(
        "[--help] [--verbose]\n    [--mode",
        "[--help] [--verbose] [--mode",
        1,
    )
    split = rendered.replace(
        "[--help] [--verbose]",
        "[--help]\n    [--verbose]",
        1,
    )
    misordered = rendered.replace(
        "[--help] [--verbose]\n    [--mode",
        "[--mode]\n    [--help] [--verbose",
        1,
    )

    for fault in (merged, split, misordered):
        findings = validate_order(
            ROOT,
            payload,
            {"compute_input_floor_cli": {"rendered_help": fault}},
        )
        assert "HELP.OPTION.ORDER.GROUP" in {
            item["rule_id"] for item in findings
        }


def test_actual_python_parser_divergence_is_not_hidden_by_config_arrays() -> (
    None
):
    payload = data()
    source = (
        ROOT / "src/protocol_chipseq_signal_norm/cli/compute_input_floor.py"
    ).read_text(encoding="utf-8")
    source = source.replace('"--fil_in",', '"--fil_input_fault",', 1)
    findings = validate_order(
        ROOT,
        payload,
        {
            "compute_input_floor_cli": {
                "parser_text": source,
            },
        },
    )
    assert "HELP.OPTION.ORDER.SURFACE_PARITY" in {
        item["rule_id"] for item in findings
    }


def test_actual_usage_and_parameter_sections_fail_independently() -> None:
    payload = data()
    command = payload["option_realizations"][0]["adapters"]["help_command"]
    rendered = render_help(ROOT, command)

    usage_fault = rendered.replace("[--fil_in FIL_IN] ", "", 1)
    usage_findings = validate_order(
        ROOT,
        payload,
        {
            "compute_input_floor_cli": {
                "rendered_help": usage_fault,
            },
        },
    )
    assert any("usage differs" in item["message"] for item in usage_findings)

    option_fault = rendered.replace("  -fi, --fil_in FIL_IN", "  FIL_IN", 1)
    option_findings = validate_order(
        ROOT,
        payload,
        {
            "compute_input_floor_cli": {
                "rendered_help": option_fault,
            },
        },
    )
    assert any(
        "options_parameters differs" in item["message"]
        for item in option_findings
    )


def test_actual_reporting_order_fault_is_detected() -> None:
    payload = data()
    reporting = "\n".join(
        [
            "--mode dist",
            "--verbose",
            "--fil_in signal.bdg",
        ],
    )
    findings = validate_order(
        ROOT,
        payload,
        {
            "compute_input_floor_cli": {
                "reporting_texts": [reporting],
            },
        },
    )
    assert any(
        "reporting_order differs" in item["message"] for item in findings
    )


def test_reporter_bindings_are_closed_and_checked_before_invocation() -> None:
    record = next(
        item
        for item in data()["option_realizations"]
        if item["surface_id"] == "compute_pseudo_cli"
    )
    adapters = record["adapters"]

    with pytest.raises(RuntimeError, match="closed data binding"):
        python_reporting_order(
            ROOT,
            Path(adapters["parser_path"]),
            adapters["reporting_symbol"],
            adapters["reporting_probes"],
            [{"source": "environment_variable"}],
        )
    with pytest.raises(RuntimeError, match="do not match"):
        python_reporting_order(
            ROOT,
            Path(adapters["parser_path"]),
            adapters["reporting_symbol"],
            adapters["reporting_probes"],
            [{"source": "parsed_args"}],
        )


def test_actual_shell_parser_missing_option_is_detected() -> None:
    payload = data()
    source = (ROOT / "bin/combine_parts_scaling_factor.sh").read_text(
        encoding="utf-8",
    )
    source = source.replace("-f|--force)", "-f)", 1)
    findings = validate_order(
        ROOT,
        payload,
        {
            "combine_parts_scaling_factor_help": {
                "parser_text": source,
            },
        },
    )
    assert "HELP.OPTION.ORDER.SURFACE_PARITY" in {
        item["rule_id"] for item in findings
    }
