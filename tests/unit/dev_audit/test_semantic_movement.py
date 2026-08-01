#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_semantic_movement.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.

"""
Test semantic-movement completeness and human authority.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path

from dev.audit.semantic_movement import validate_record

ROOT = Path(__file__).resolve().parents[3]
CASES = ROOT / "tests/fixtures/semantic_movement/cases.json"
SCHEMA = ROOT / "dev/schemas/semantic_movement.schema.json"


def schema() -> dict:
    return json.loads(SCHEMA.read_text(encoding="utf-8"))


def test_accepted_movement_is_complete() -> None:
    cases = json.loads(CASES.read_text(encoding="utf-8"))
    assert validate_record(cases["accepted"], schema()) == []


def test_consequential_delta_requires_human_approval() -> None:
    cases = json.loads(CASES.read_text(encoding="utf-8"))
    record = copy.deepcopy(cases["accepted"])
    record["consequential_delta"] = True
    record["human_authorization"]["required"] = True
    record["human_authorization"]["status"] = "pending"

    assert validate_record(record, schema()) == [
        "consequential delta lacks explicit human approval",
    ]


def test_llm_analysis_cannot_be_the_authorizing_role() -> None:
    cases = json.loads(CASES.read_text(encoding="utf-8"))
    record = copy.deepcopy(cases["accepted"])
    record["llm_role"] = "authorizing_change"

    errors = validate_record(record, schema())
    assert any(
        "supporting_comparison" in error or "supporting comparison" in error
        for error in errors
    )


def test_supplied_schema_rejects_missing_and_extra_fields() -> None:
    cases = json.loads(CASES.read_text(encoding="utf-8"))
    record = copy.deepcopy(cases["accepted"])
    del record["deltas"]
    assert validate_record(record, schema())[0].startswith("schema ")

    record = copy.deepcopy(cases["accepted"])
    record["invented"] = True
    assert validate_record(record, schema())[0].startswith("schema ")
