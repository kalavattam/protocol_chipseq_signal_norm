#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: semantic_movement.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Validate lossless semantic-movement evidence without authorizing changes.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from jsonschema import Draft202012Validator


def schema_errors(
    record: dict[str, Any],
    schema: dict[str, Any],
) -> list[str]:
    """
    Return deterministic Draft 2020-12 validation errors.
    """

    validator = Draft202012Validator(schema)
    return [
        (
            f"schema {'/'.join(str(part) for part in error.absolute_path) or '/'}"
            f": {error.message}"
        )
        for error in sorted(
            validator.iter_errors(record),
            key=lambda item: (
                tuple(str(part) for part in item.absolute_path),
                item.message,
            ),
        )
    ]


def validate_record(
    record: dict[str, Any],
    schema: dict[str, Any] | None = None,
) -> list[str]:
    """
    Return schema, completeness, preservation, and authority failures.
    """

    if schema is not None:
        errors = schema_errors(record, schema)
        if errors:
            return errors

    errors: list[str] = []
    if record.get("old_to_new_complete") is not True:
        errors.append("old-to-new dispositions are incomplete")
    if not record.get("old_to_new_dispositions"):
        errors.append("old-to-new disposition evidence is empty")
    if not record.get("new_to_old_sources"):
        errors.append("new-to-old provenance is empty")
    if not record.get("new_to_old_provenance"):
        errors.append("new-to-old fact provenance is empty")
    if record.get("availability_reduced") is not False:
        errors.append("audience availability was reduced")
    if record.get("behavior_changed") is not False:
        errors.append("protected behavior changed")
    if record.get("llm_role") != "supporting_comparison":
        errors.append("LLM role may only be supporting comparison")

    roles = {
        item.get("role")
        for item in record.get("evidence_roles", [])
        if isinstance(item, dict)
    }
    if "supporting_comparison" not in roles:
        errors.append("supporting comparison evidence role is absent")

    authorization = record.get("human_authorization", {})
    if record.get("consequential_delta") and (
        not authorization.get("required")
        or authorization.get("status") != "approved"
        or not authorization.get("reference")
    ):
        errors.append("consequential delta lacks explicit human approval")
    if authorization.get("status") == "approved" and not authorization.get(
        "reference",
    ):
        errors.append("approved human authorization lacks an exact source")

    return errors


def main(argv: list[str] | None = None) -> int:
    """
    Validate one record or an array of records against a supplied schema.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("--schema", type=Path, required=True)
    parser.add_argument("--validate", type=Path, required=True)
    args = parser.parse_args(argv)
    schema = json.loads(args.schema.read_text(encoding="utf-8"))
    payload = json.loads(args.validate.read_text(encoding="utf-8"))
    records = payload if isinstance(payload, list) else [payload]
    errors = [
        f"record[{index}]: {error}"
        for index, record in enumerate(records)
        for error in validate_record(record, schema)
    ]
    for error in errors:
        print(error)
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
