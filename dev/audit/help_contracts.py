#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: help_contracts.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Validate shared help applicability and permitted diagnostic emitters.
"""

from __future__ import annotations

import argparse
import ast
import json
import re
import tomllib
from collections import Counter
from pathlib import Path
from typing import Any

from jsonschema import Draft202012Validator

RULE_SCHEMA = "HELP.CONTRACT.SCHEMA"
RULE_REFERENCE = "HELP.CONTRACT.REFERENCE"
RULE_APPLICABILITY = "HELP.CONTRACT.APPLICABILITY"
ALLOWED_DEFERRED_EXAMPLE_RECORDS = {"S3-MIG-001"}

EXPECTED_CONSUMERS = {
    "dev.audit.help_aliases",
    "dev.audit.help_contracts",
    "dev.audit.help_examples",
    "dev.audit.help_option_order",
    "dev.audit.help_style",
    "tests.contract.repository.test_parameter_docs_consistency",
}
AUDIENCE_COMBINATIONS = {
    "installed_help": ("command_user", "installed", "complete"),
    "callable_documentation": (
        "python_caller",
        "installed",
        "direct_call_contract",
    ),
    "public_design": (
        "design_reviewer",
        "public_repository",
        "extended_derivation",
    ),
}


def _module_path(root: Path, module: str) -> Path:
    """
    Resolve one configured Python-module or Shell-contract identity.
    """

    suffix = ".sh" if module.startswith("tests.contract.") else ".py"
    return root / (module.replace(".", "/") + suffix)


def _symbol_exists(path: Path, language: str, symbol: str) -> bool:
    """
    Return whether one registered source symbol exists.
    """

    if symbol == "document":
        return language == "markdown"

    text = path.read_text(encoding="utf-8")
    if language == "python":
        try:
            tree = ast.parse(text)
        except SyntaxError:
            return False
        return any(
            isinstance(
                node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)
            )
            and node.name == symbol
            for node in ast.walk(tree)
        )

    if language == "shell":
        return (
            re.search(
                rf"(?m)^(?:function\s+)?{re.escape(symbol)}\s*(?:\(\))?\s*\{{",
                text,
            )
            is not None
        )

    return False


def _symbol_text(path: Path, language: str, symbol: str) -> str | None:
    """Return one bounded configured symbol body, or 'None' when absent."""

    text = path.read_text(encoding="utf-8")
    if symbol == "<file>":
        return text
    if language == "python":
        try:
            tree = ast.parse(text)
        except SyntaxError:
            return None
        for node in ast.walk(tree):
            if (
                isinstance(
                    node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)
                )
                and node.name == symbol
            ):
                return ast.get_source_segment(text, node)
        return None
    match = re.search(
        rf"(?m)^(?:function\s+)?{re.escape(symbol)}\s*(?:\(\))?\s*\{{",
        text,
    )
    if match is None:
        return None
    start = match.start()
    depth = 0
    for index in range(match.start(), len(text)):
        if text[index] == "{":
            depth += 1
        elif text[index] == "}":
            depth -= 1
            if depth == 0:
                return text[start : index + 1]
    return None


def _parameter_occurrences(text: str, language: str, parameter: str) -> int:
    """Count declared help/argument occurrences in one bounded source surface."""

    if language == "python":
        try:
            tree = ast.parse(text)
        except SyntaxError:
            return 0
        return sum(
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "add_argument"
            and any(
                isinstance(argument, ast.Constant)
                and argument.value == f"--{parameter}"
                for argument in node.args
            )
            for node in ast.walk(tree)
        )
    return len(
        re.findall(
            rf"(?m)^[^\n]*--{re.escape(parameter)}\s*:\s*|"
            rf"^\s*\d+\s+{re.escape(parameter)}\s*:\s*",
            text,
        )
    )


def _registry_diagnostics(root: Path) -> dict[str, str]:
    """
    Map every registered diagnostic ID to its unique authoritative owner.
    """

    registry = root / "dev/config/rules.toml"
    data = tomllib.loads(registry.read_text(encoding="utf-8"))
    owners: dict[str, str] = {}
    duplicates: set[str] = set()
    for rule in data.get("rule", []):
        owner = str(rule.get("rule_id", ""))
        for diagnostic in rule.get("diagnostic_ids", []):
            diagnostic = str(diagnostic)
            if diagnostic in owners:
                duplicates.add(diagnostic)
            owners[diagnostic] = owner

    if duplicates:
        raise ValueError(
            "standards registry has duplicate diagnostic IDs: "
            + ", ".join(sorted(duplicates)),
        )
    return owners


def _schema_findings(
    data: dict[str, Any],
    schema: dict[str, Any],
) -> list[dict[str, Any]]:
    """
    Return deterministic Draft 2020-12 schema findings.
    """

    validator = Draft202012Validator(schema)
    findings: list[dict[str, Any]] = []
    for error in sorted(
        validator.iter_errors(data),
        key=lambda item: (
            tuple(str(part) for part in item.absolute_path),
            item.message,
        ),
    ):
        locator = "/".join(str(part) for part in error.absolute_path) or "/"
        findings.append(
            {
                "rule_id": RULE_SCHEMA,
                "message": f"{locator}: {error.message}",
            },
        )
    return findings


def validate_contract(
    root: Path,
    data: dict[str, Any],
    schema: dict[str, Any] | None = None,
) -> list[dict[str, Any]]:
    """
    Return exact schema, reference, audience, and assignment findings.

    This checker validates shared data. It never diagnoses delegated alias,
    description, example, vocabulary, token-quoting, or option-order defects.
    """

    root = root.resolve()
    if schema is None:
        schema = json.loads(
            (root / "dev/schemas/help_contracts.schema.json").read_text(
                encoding="utf-8"
            ),
        )

    findings = _schema_findings(data, schema)
    if findings:
        return findings

    def add(rule_id: str, message: str) -> None:
        findings.append({"rule_id": rule_id, "message": message})

    surfaces = data["surfaces"]
    dispositions = data["parameter_occurrence_dispositions"]
    occurrence_counts = Counter(item["id"] for item in dispositions)
    for occurrence, count in sorted(occurrence_counts.items()):
        if count != 1:
            add(
                RULE_APPLICABILITY,
                f"parameter occurrence must have exactly one disposition: {occurrence}",
            )
    identity_counts = Counter(
        (item["path"], item["symbol"], item["parameter"])
        for item in dispositions
    )
    for identity, count in sorted(identity_counts.items()):
        if count != 1:
            add(
                RULE_APPLICABILITY,
                "parameter occurrence identity must have exactly one disposition: "
                f"{identity[0]}::{identity[1]}::{identity[2]}",
            )
    family_members = {
        (member["path"], member["symbol"], member["parameter"])
        for family in data["parameter_families"]
        if family["status"] == "registered_shared"
        for member in family["members"]
    }
    for record in dispositions:
        path = root / record["path"]
        identity = (record["path"], record["symbol"], record["parameter"])
        if not path.is_file():
            add(
                RULE_APPLICABILITY,
                f"parameter occurrence path does not exist: {record['id']}",
            )
            continue
        language = "python" if path.suffix == ".py" else "shell"
        symbol_text = _symbol_text(path, language, record["symbol"])
        if symbol_text is None:
            add(
                RULE_APPLICABILITY,
                f"parameter occurrence symbol does not exist: {record['id']}",
            )
        elif (
            _parameter_occurrences(symbol_text, language, record["parameter"])
            != 1
        ):
            add(
                RULE_APPLICABILITY,
                "parameter occurrence is not exactly once in recorded symbol: "
                f"{record['id']}",
            )
        enrolled = record["disposition"] == "enrolled_member"
        if enrolled and identity not in family_members:
            add(
                RULE_APPLICABILITY,
                f"enrolled parameter occurrence is not a registered member: {record['id']}",
            )
        if not enrolled and identity in family_members:
            add(
                RULE_APPLICABILITY,
                f"non-enrolled disposition masquerades as member: {record['id']}",
            )
    surface_counts = Counter(item["id"] for item in surfaces)
    for identifier, count in sorted(surface_counts.items()):
        if count != 1:
            add(
                RULE_REFERENCE,
                f"surface ID must occur exactly once: {identifier}",
            )
    surface_by_id = {item["id"]: item for item in surfaces}

    for record in surfaces:
        path = root / record["path"]
        if not path.is_file():
            add(
                RULE_REFERENCE,
                f"registered surface path does not exist: {record['path']}",
            )
            continue
        if not _symbol_exists(path, record["language"], record["symbol"]):
            add(
                RULE_REFERENCE,
                "registered surface symbol does not exist: "
                f"{record['path']}::{record['symbol']}",
            )

        expected = AUDIENCE_COMBINATIONS[record["kind"]]
        observed = (
            record["audience"],
            record["availability"],
            record["required_runtime_facts"],
        )
        if observed != expected:
            add(
                RULE_APPLICABILITY,
                f"{record['id']} has an invalid audience/availability/"
                "runtime-facts combination",
            )
        if (
            record["kind"] == "callable_documentation"
            and record["language"] != "python"
        ):
            add(
                RULE_APPLICABILITY,
                f"{record['id']} callable documentation must be Python",
            )
        if (
            record["kind"] == "public_design"
            and record["language"] != "markdown"
        ):
            add(
                RULE_APPLICABILITY,
                f"{record['id']} public design must be Markdown",
            )

    for section in ("examples", "option_realizations", "format_vocabulary"):
        for record in data[section]:
            surface = surface_by_id.get(record["surface_id"])
            if surface is None:
                add(
                    RULE_REFERENCE,
                    f"{section} references unknown surface "
                    f"{record['surface_id']}",
                )
                continue
            if (
                "language" in record
                and record["language"] != surface["language"]
            ):
                add(
                    RULE_APPLICABILITY,
                    f"{section} language differs from surface "
                    f"{record['surface_id']}",
                )

    family_ids = Counter(item["id"] for item in data["parameter_families"])
    parameters = Counter(
        item["parameter"]
        for item in data["parameter_families"]
        if item["status"] == "registered_shared"
    )
    for identifier, count in sorted(family_ids.items()):
        if count != 1:
            add(
                RULE_REFERENCE,
                f"parameter-family ID must occur exactly once: {identifier}",
            )
    for parameter, count in sorted(parameters.items()):
        if count != 1:
            add(
                RULE_APPLICABILITY,
                "registered shared parameter must occur exactly once: "
                f"{parameter}",
            )

    for family in data["parameter_families"]:
        records = (
            family["members"]
            if family["status"] == "registered_shared"
            else family["local_meanings"]
        )
        member_keys: set[tuple[str, str, str, str]] = set()
        for member in records:
            if member["parameter"] != family["parameter"]:
                add(
                    RULE_APPLICABILITY,
                    f"{family['id']} member parameter does not match family",
                )
            key = (
                member["surface_id"],
                member["path"],
                member["symbol"],
                member["parameter"],
            )
            if key in member_keys:
                add(
                    RULE_REFERENCE,
                    f"{family['id']} contains a duplicate member identity",
                )
            member_keys.add(key)
            path = root / member["path"]
            if not path.is_file():
                add(
                    RULE_REFERENCE,
                    f"parameter member path does not exist: {member['path']}",
                )
            elif member["symbol"] != "<file>":
                language = "python" if path.suffix == ".py" else "shell"
                if not _symbol_exists(path, language, member["symbol"]):
                    add(
                        RULE_REFERENCE,
                        "parameter member symbol does not exist: "
                        f"{member['path']}::{member['symbol']}",
                    )
            if member["owner"] != "PARAMETER.DESCRIPTIONS":
                add(
                    RULE_APPLICABILITY,
                    f"{family['id']} member has a different owner",
                )

    example_counts = Counter(item["surface_id"] for item in data["examples"])
    for surface_id, count in sorted(example_counts.items()):
        if count != 1:
            add(
                RULE_REFERENCE,
                f"examples disposition must occur once: {surface_id}",
            )
    for record in data["examples"]:
        expected_minimum = {
            "required_two": 2,
            "required_one": 1,
            "omitted_trivial": 0,
        }[record["disposition"]]
        if record["minimum_count"] != expected_minimum:
            add(
                RULE_APPLICABILITY,
                f"{record['surface_id']} has an inconsistent example minimum",
            )
        if (
            len(record["example_fingerprints"])
            != record["existing_example_count"]
        ):
            add(
                RULE_APPLICABILITY,
                f"{record['surface_id']} example fingerprints/count differ",
            )
        if (
            record["lifecycle"] == "active_enforced"
            and record["existing_example_count"] < expected_minimum
        ):
            add(
                RULE_APPLICABILITY,
                f"{record['surface_id']} active example count is below minimum",
            )
        if (
            record["lifecycle"] == "deferred_migration"
            and record["deferred_record"]
            not in ALLOWED_DEFERRED_EXAMPLE_RECORDS
        ):
            add(
                RULE_APPLICABILITY,
                f"{record['surface_id']} has an unknown deferred record",
            )

    realization_counts = Counter(
        item["surface_id"] for item in data["option_realizations"]
    )
    for surface_id, count in sorted(realization_counts.items()):
        if count != 1:
            add(
                RULE_REFERENCE,
                f"option realization must occur once: {surface_id}",
            )
    for record in data["option_realizations"]:
        order = record["logical_order"]
        if set(record["roles"]) != set(order):
            add(
                RULE_APPLICABILITY,
                f"{record['surface_id']} roles are not total",
            )
        for name, surface_order in record["surface_orders"].items():
            if isinstance(surface_order, list) and not set(
                surface_order,
            ) <= set(order):
                add(
                    RULE_APPLICABILITY,
                    f"{record['surface_id']} {name} contains an unknown option",
                )
        usage_members = [
            member for row in record["usage_rows"] for member in row["members"]
        ]
        if len(usage_members) != len(set(usage_members)):
            add(
                RULE_APPLICABILITY,
                f"{record['surface_id']} usage rows repeat an option",
            )
        if usage_members != order:
            add(
                RULE_APPLICABILITY,
                f"{record['surface_id']} usage rows are not a complete ordered projection",
            )
        for row in record["usage_rows"]:
            roles = {record["roles"].get(member) for member in row["members"]}
            if len(roles) > 1 and not row["rationale"].strip():
                add(
                    RULE_APPLICABILITY,
                    f"{record['surface_id']} multi-role usage row lacks a rationale",
                )
        for relationship in record["relationships"]:
            if not set(relationship["members"]) <= set(order):
                add(
                    RULE_APPLICABILITY,
                    f"{record['surface_id']} relationship has an unknown option",
                )

    consumer_counts = Counter(
        item["checker"] for item in data["checker_consumers"]
    )
    if set(consumer_counts) != EXPECTED_CONSUMERS:
        missing = sorted(EXPECTED_CONSUMERS - set(consumer_counts))
        extra = sorted(set(consumer_counts) - EXPECTED_CONSUMERS)
        add(
            RULE_APPLICABILITY,
            f"checker-consumer matrix differs; missing={missing!r}; "
            f"extra={extra!r}",
        )
    for checker, count in sorted(consumer_counts.items()):
        if count != 1:
            add(
                RULE_APPLICABILITY,
                f"checker consumer must occur exactly once: {checker}",
            )
        if not _module_path(root, checker).is_file():
            add(
                RULE_REFERENCE,
                f"checker consumer does not exist: {checker}",
            )

    try:
        registry = _registry_diagnostics(root)
    except ValueError as error:
        add(RULE_REFERENCE, str(error))
        registry = {}

    assignment_counts = Counter(
        item["diagnostic_id"] for item in data["checker_assignments"]
    )
    consumer_by_id = {
        item["checker"]: item for item in data["checker_consumers"]
    }
    for diagnostic, count in sorted(assignment_counts.items()):
        if count != 1:
            add(
                RULE_APPLICABILITY,
                f"{diagnostic} must have exactly one permitted emitter",
            )

    for record in data["checker_assignments"]:
        diagnostic = record["diagnostic_id"]
        emitter = record["permitted_emitter"]
        registered_owner = registry.get(diagnostic)
        if registered_owner is None:
            add(
                RULE_REFERENCE,
                f"configured diagnostic is absent from registry: {diagnostic}",
            )
        elif registered_owner != record["owner"]:
            add(
                RULE_APPLICABILITY,
                f"{diagnostic} owner differs from standards registry",
            )

        if not _module_path(root, emitter).is_file():
            add(
                RULE_REFERENCE,
                f"permitted emitter does not exist: {emitter}",
            )
        consumer = consumer_by_id.get(emitter)
        if consumer is None:
            add(
                RULE_REFERENCE,
                f"permitted emitter is absent from consumer matrix: {emitter}",
            )
        elif diagnostic not in consumer["diagnostic_scope"]:
            add(
                RULE_APPLICABILITY,
                f"{diagnostic} is outside {emitter} diagnostic scope",
            )
        if emitter in record["prohibited_emitters"]:
            add(
                RULE_APPLICABILITY,
                f"{diagnostic} permits and prohibits the same emitter",
            )
        for prohibited in record["prohibited_emitters"]:
            prohibited_consumer = consumer_by_id.get(prohibited)
            if (
                prohibited_consumer is not None
                and diagnostic
                not in prohibited_consumer["prohibited_diagnostics"]
            ):
                add(
                    RULE_APPLICABILITY,
                    f"{diagnostic} prohibition is not bidirectional for "
                    f"{prohibited}",
                )

    return findings


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse command-line arguments.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="*", type=Path)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("dev/config/help_contracts.json"),
    )
    parser.add_argument(
        "--schema",
        type=Path,
        default=Path("dev/schemas/help_contracts.schema.json"),
    )
    parser.add_argument("--json", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """
    Validate the configured contract without diagnosing delegated defects.
    """

    args = parse_args(argv)
    root = args.root.resolve()
    config = args.config if args.config.is_absolute() else root / args.config
    schema = args.schema if args.schema.is_absolute() else root / args.schema
    data = json.loads(config.read_text(encoding="utf-8"))
    schema_data = json.loads(schema.read_text(encoding="utf-8"))
    findings = validate_contract(root, data, schema_data)

    if args.json:
        print(json.dumps(findings, indent=2, sort_keys=True))
    else:
        for finding in findings:
            print(f"{finding['rule_id']}: {finding['message']}")
    return 1 if findings else 0


if __name__ == "__main__":
    raise SystemExit(main())
