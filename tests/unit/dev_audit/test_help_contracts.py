#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_help_contracts.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Test help applicability, permitted emitters, and checker non-overlap.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path

from dev.audit.help_aliases import check_document
from dev.audit.help_contracts import validate_contract
from dev.audit.help_style import check_rendered_help

ROOT = Path(__file__).resolve().parents[3]
CONFIG = ROOT / "dev/config/help_contracts.json"
FIXTURES = ROOT / "tests/fixtures/help_contracts"


def data() -> dict:
    """
    Return an independent config copy.
    """

    return json.loads(CONFIG.read_text(encoding="utf-8"))


def test_example_forms_separate_applicability_from_representation() -> None:
    payload = data()
    callable_record = next(
        item
        for item in payload["examples"]
        if item["surface_id"] == "compute_input_floor_callable"
    )

    assert payload["schema_version"] == 3
    assert callable_record == {
        "surface_id": "compute_input_floor_callable",
        "language": "python",
        "owner": "compute_input_floor",
        "example_form": "callable_source_language",
        "representation_owner": "PY.DOCSTRING.NUMPY",
        "lifecycle": "deferred_migration",
        "deferred_record": "S3-MIG-001",
        "disposition": "required_two",
        "minimum_count": 2,
        "existing_example_count": 0,
        "example_fingerprints": [],
        "preservation": "exact",
        "authority": "semantic_reaudit_correction",
    }
    assert validate_contract(ROOT, payload) == []


def test_fresh_pilot_records_require_closed_schema_three_data() -> None:
    payload = data()
    callable_record = next(
        item
        for item in payload["examples"]
        if item["surface_id"] == "compute_pseudo_callable"
    )
    realization = next(
        item
        for item in payload["option_realizations"]
        if item["surface_id"] == "compute_pseudo_cli"
    )

    assert callable_record["authority"] == "authorized_fresh_pilot"
    assert callable_record["lifecycle"] == "active_enforced"
    assert realization["adapters"]["reporting_call_arguments"] == [
        {"source": "parsed_args"},
        {"literal": None},
        {"literal": ["track"]},
    ]
    assert validate_contract(ROOT, payload) == []


def test_schema_three_rejects_unknown_authority_and_reporter_binding() -> None:
    payload = data()
    record = next(
        item
        for item in payload["examples"]
        if item["surface_id"] == "compute_pseudo_callable"
    )
    record["authority"] = "unreviewed"

    assert "HELP.CONTRACT.SCHEMA" in {
        item["rule_id"] for item in validate_contract(ROOT, payload)
    }

    payload = data()
    record = next(
        item
        for item in payload["option_realizations"]
        if item["surface_id"] == "compute_pseudo_cli"
    )
    record["adapters"]["reporting_call_arguments"][0] = {
        "source": "environment_variable"
    }

    assert "HELP.CONTRACT.SCHEMA" in {
        item["rule_id"] for item in validate_contract(ROOT, payload)
    }


def test_deferred_example_record_is_bounded_to_approved_destination() -> None:
    payload = data()
    callable_record = next(
        item
        for item in payload["examples"]
        if item["surface_id"] == "compute_input_floor_callable"
    )
    callable_record["deferred_record"] = "S3-MIG-999"

    assert {item["rule_id"] for item in validate_contract(ROOT, payload)} == {
        "HELP.CONTRACT.APPLICABILITY"
    }


def test_checker_assignment_rejects_two_permitted_emitters() -> None:
    payload = data()
    duplicate = copy.deepcopy(payload["checker_assignments"][0])
    duplicate["permitted_emitter"] = "dev.audit.help_contracts"
    payload["checker_assignments"].append(duplicate)

    assert {item["rule_id"] for item in validate_contract(ROOT, payload)} == {
        "HELP.CONTRACT.APPLICABILITY"
    }

    payload = data()
    payload["checker_assignments"].append(
        copy.deepcopy(payload["checker_assignments"][0]),
    )

    assert {item["rule_id"] for item in validate_contract(ROOT, payload)} == {
        "HELP.CONTRACT.APPLICABILITY"
    }


def test_incomplete_nested_record_is_a_schema_finding() -> None:
    payload = data()
    del payload["format_vocabulary"][0]["recognized_suffixes"]

    assert "HELP.CONTRACT.SCHEMA" in {
        item["rule_id"] for item in validate_contract(ROOT, payload)
    }

    payload = data()
    del payload["option_realizations"][0]["usage_rows"]

    assert "HELP.CONTRACT.SCHEMA" in {
        item["rule_id"] for item in validate_contract(ROOT, payload)
    }


def test_usage_rows_require_complete_unique_membership_and_rationales() -> (
    None
):
    payload = data()
    payload["option_realizations"][0]["usage_rows"][0]["members"].append(
        "mode",
    )

    assert "HELP.CONTRACT.APPLICABILITY" in {
        item["rule_id"] for item in validate_contract(ROOT, payload)
    }

    payload = data()
    record = next(
        item
        for item in payload["option_realizations"]
        if item["surface_id"] == "combine_parts_scaling_factor_help"
    )
    record["usage_rows"][0]["rationale"] = " "

    assert "HELP.CONTRACT.APPLICABILITY" in {
        item["rule_id"] for item in validate_contract(ROOT, payload)
    }


def test_schema_enum_and_extra_property_faults_are_not_ignored() -> None:
    payload = data()
    payload["surfaces"][0]["audience"] = "invented"
    assert {item["rule_id"] for item in validate_contract(ROOT, payload)} == {
        "HELP.CONTRACT.SCHEMA"
    }

    payload = data()
    payload["surfaces"][0]["extra"] = True
    assert {item["rule_id"] for item in validate_contract(ROOT, payload)} == {
        "HELP.CONTRACT.SCHEMA"
    }


def test_reference_audience_and_realization_faults_are_distinct() -> None:
    payload = data()
    payload["examples"][0]["surface_id"] = "missing"
    assert "HELP.CONTRACT.REFERENCE" in {
        item["rule_id"] for item in validate_contract(ROOT, payload)
    }

    payload = data()
    payload["surfaces"][0]["availability"] = "public_repository"
    assert "HELP.CONTRACT.APPLICABILITY" in {
        item["rule_id"] for item in validate_contract(ROOT, payload)
    }

    payload = data()
    payload["option_realizations"][0]["language"] = "shell"
    assert "HELP.CONTRACT.APPLICABILITY" in {
        item["rule_id"] for item in validate_contract(ROOT, payload)
    }


def test_consumer_matrix_and_registry_owner_are_enforced_separately() -> None:
    payload = data()
    payload["checker_consumers"].pop()
    assert "HELP.CONTRACT.APPLICABILITY" in {
        item["rule_id"] for item in validate_contract(ROOT, payload)
    }

    payload = data()
    payload["checker_assignments"][0]["owner"] = "HELP.EXAMPLES"
    assert "HELP.CONTRACT.APPLICABILITY" in {
        item["rule_id"] for item in validate_contract(ROOT, payload)
    }


def test_shared_data_does_not_grant_parallel_alias_emission() -> None:
    payload = data()
    alias = next(
        item
        for item in payload["checker_assignments"]
        if item["diagnostic_id"] == "HELP.PARAMETER.ALIAS_DUPLICATE"
    )
    alias["permitted_emitter"] = "dev.audit.help_style"

    assert "HELP.CONTRACT.APPLICABILITY" in {
        item["rule_id"] for item in validate_contract(ROOT, payload)
    }


def test_shell_alias_duplicate_has_only_authoritative_diagnostic() -> None:
    text = (FIXTURES / "shell.sh.fixture").read_text(encoding="utf-8")
    alias_findings, _ = check_document(
        ROOT,
        "lib/bash/core/fixture.sh",
        text,
    )
    style_ids = {
        finding.rule_id
        for finding in check_rendered_help(
            text,
            registry_path=ROOT / "dev/config/command_names.json",
        )
    }

    assert [finding.message for finding in alias_findings] == [
        "duplicate alias",
    ]
    assert [finding.diagnostic_id for finding in alias_findings] == [
        "HELP.PARAMETER.ALIAS_DUPLICATE",
    ]
    assert [
        finding.format().split(":", 1)[0] for finding in alias_findings
    ] == ["HELP.PARAMETER.ALIAS_DUPLICATE"]
    assert "HELP.PARAMETER.ALIAS_DUPLICATE" not in style_ids


def test_parameter_disposition_identity_and_enrollment_faults_are_owned() -> (
    None
):
    payload = data()
    duplicate = copy.deepcopy(payload["parameter_occurrence_dispositions"][0])
    duplicate["id"] = "different_id_same_identity"
    payload["parameter_occurrence_dispositions"].append(duplicate)
    assert "HELP.CONTRACT.APPLICABILITY" in {
        item["rule_id"] for item in validate_contract(ROOT, payload)
    }


def test_parameter_elsewhere_in_file_does_not_satisfy_recorded_symbol() -> (
    None
):
    payload = data()
    record = next(
        item
        for item in payload["parameter_occurrence_dispositions"]
        if item["id"] == "source_helpers_dir_fnc_88"
    )
    record["symbol"] = "_source_helper_err"
    messages = [item["message"] for item in validate_contract(ROOT, payload)]
    assert any(
        "not exactly once in recorded symbol" in message
        for message in messages
    )

    payload = data()
    payload["parameter_occurrence_dispositions"][0]["symbol"] = (
        "missing_symbol"
    )
    assert "HELP.CONTRACT.APPLICABILITY" in {
        item["rule_id"] for item in validate_contract(ROOT, payload)
    }

    payload = data()
    payload["parameter_occurrence_dispositions"][0]["parameter"] = (
        "missing_parameter"
    )
    assert "HELP.CONTRACT.APPLICABILITY" in {
        item["rule_id"] for item in validate_contract(ROOT, payload)
    }

    payload = data()
    payload["parameter_occurrence_dispositions"][0]["disposition"] = (
        "non_applicable_same_name"
    )
    assert "HELP.CONTRACT.APPLICABILITY" in {
        item["rule_id"] for item in validate_contract(ROOT, payload)
    }
