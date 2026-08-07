#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: standards_registry.py
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
Audit traceability between standards, rules, checkers, and fixtures.
"""

from __future__ import annotations

import argparse
import ast
import dataclasses
import json
import re
import tomllib
from collections import Counter
from pathlib import Path

HEADING = re.compile(r"^#{1,6}\s+(.+?)\s*$", re.MULTILINE)
H2 = re.compile(r"^##\s+(.+?)\s*$", re.MULTILINE)
RULE_ID = re.compile(r"`([A-Z][A-Z0-9]*(?:\.[A-Z0-9_]+)+)`")
MARKDOWN_LINK = re.compile(r"\]\(([^)]+\.md)\)")
INDEX_STANDARD_LINK = re.compile(
    r"^\|\s*\[`[^`]+`\]\(([^)]+\.md)\)\s*\|",
    re.MULTILINE,
)
FENCE = re.compile(r"^ {0,3}(?P<marker>`{3,}|~{3,})")
OWNER_FIELD = re.compile(
    r"^\*\*(?P<label>Classification|Scope|Automation|"
    r"Semantic remainder|Exceptions):\*\*\s+(?P<value>.+)$",
)
OWNER_FIELDS = (
    "Classification",
    "Scope",
    "Automation",
    "Semantic remainder",
    "Exceptions",
)
OWNER_CLASSIFICATION = re.compile(r"`(deterministic|advisory|semantic-only)`")
EMITTED_RULE_ID = re.compile(r"^[A-Z][A-Z0-9]*(?:\.[A-Z0-9_]+)+$")
OWNER_CLASSIFICATIONS = {"deterministic", "advisory", "semantic-only"}
EXECUTION_ROLES = {"checker", "evidence_producer", "independent_evidence"}
COVERAGE_RELATIONS = {"exact", "subset", "independent"}
GENERIC_REMAINING_SCOPE = (
    "Owner scope outside the registered execution remains unproved."
)


@dataclasses.dataclass(frozen=True)
class OwnerRule:
    """
    Represent one canonical owner rule and its lifecycle fields.
    """

    rule_id: str
    normative_document: str
    normative_section: str
    owner_classification: str
    owner_scope: str
    automation: str
    semantic_remainder: str
    exceptions: str


def _is_path(value: str) -> bool:
    """
    Return whether a registry value names a repository path.
    """

    return "/" in value and not value.startswith("independent ")


def _owner_documents(
    root: Path,
    rules: list[dict[str, object]],
) -> tuple[list[Path], list[dict[str, str]]]:
    """
    Return indexed owner documents and fail-closed index findings.
    """

    index = root / "docs/standards/README.md"
    findings: list[dict[str, str]] = []

    if index.is_file():
        relatives = INDEX_STANDARD_LINK.findall(
            index.read_text(encoding="utf-8"),
        )

        for relative, count in sorted(Counter(relatives).items()):
            if count > 1:
                findings.append(
                    {
                        "kind": "duplicate_index_entry",
                        "rule_id": "<standards-index>",
                        "detail": f"{relative}: count={count}",
                    },
                )

        documents: list[Path] = []

        for relative in dict.fromkeys(relatives):
            if relative == "README.md" or Path(relative).name.startswith(
                "bak.",
            ):
                continue

            document = index.parent / relative

            if not document.is_file():
                findings.append(
                    {
                        "kind": "missing_index_target",
                        "rule_id": "<standards-index>",
                        "detail": relative,
                    },
                )

                continue

            documents.append(document)

        discovered = {
            path
            for path in index.parent.glob("*.md")
            if path.name != "README.md" and not path.name.startswith("bak.")
        }

        for document in sorted(discovered - set(documents)):
            findings.append(
                {
                    "kind": "unindexed_standard",
                    "rule_id": "<standards-index>",
                    "detail": document.relative_to(root).as_posix(),
                },
            )

        return sorted(set(documents)), findings

    findings.append(
        {
            "kind": "missing_standards_index",
            "rule_id": "<standards-index>",
            "detail": "docs/standards/README.md",
        },
    )
    documents = sorted(
        {
            root / str(rule.get("normative_document", ""))
            for rule in rules
            if (root / str(rule.get("normative_document", ""))).is_file()
        },
    )

    return documents, findings


def _owner_inventory(
    root: Path,
    rules: list[dict[str, object]],
) -> tuple[dict[str, list[OwnerRule]], list[dict[str, str]], list[str]]:
    """
    Return stable IDs mapped to their maintained document and heading.

    Parameters
    ----------
    root : Path
        Repository root containing maintained standards.
    rules : list[dict[str, object]]
        Configured rule registry entries.

    Returns
    -------
    owners, findings, documents : tuple[
        dict[str, list[OwnerRule]], list[dict[str, str]], list[str]
    ]
        Rule owners, structural findings, and maintained standards documents.
    """

    inventory: dict[str, list[OwnerRule]] = {}
    documents, findings = _owner_documents(root, rules)

    for document in documents:
        relative = document.relative_to(root).as_posix()
        lines = document.read_text(encoding="utf-8").splitlines()
        literal: set[int] = set()
        opened: tuple[str, int] | None = None

        for index, line in enumerate(lines):
            fence = FENCE.match(line)

            if fence is not None:
                marker = fence.group("marker")

                if opened is None:
                    opened = (marker[0], len(marker))
                elif marker[0] == opened[0] and len(marker) >= opened[1]:
                    opened = None

                literal.add(index)

                continue

            if opened is not None:
                literal.add(index)

        headings = [
            index
            for index, line in enumerate(lines)
            if index not in literal and line.startswith("## ")
        ]

        for position, start in enumerate(headings):
            heading = lines[start].removeprefix("## ").strip()
            end = (
                headings[position + 1]
                if position + 1 < len(headings)
                else len(lines)
            )
            identifiers = RULE_ID.findall(heading)

            if not identifiers:
                findings.append(
                    {
                        "kind": "owner_section_without_id",
                        "rule_id": "<missing>",
                        "detail": f"{relative}: {heading}",
                    },
                )

                continue

            fields: dict[str, list[tuple[int, str]]] = {
                name: [] for name in OWNER_FIELDS
            }
            ordered: list[str] = []

            for index in range(start + 1, end):
                if index in literal:
                    continue

                field = OWNER_FIELD.match(lines[index])
                if field is None:
                    continue

                label = field.group("label")
                fields[label].append((index, field.group("value").strip()))
                ordered.append(label)

            for label in OWNER_FIELDS:
                if len(fields[label]) != 1:
                    findings.append(
                        {
                            "kind": "owner_field_count",
                            "rule_id": identifiers[0],
                            "detail": (
                                f"{relative}: {heading}: {label} count is "
                                f"{len(fields[label])}"
                            ),
                        },
                    )

            if ordered != list(OWNER_FIELDS):
                findings.append(
                    {
                        "kind": "owner_field_order",
                        "rule_id": identifiers[0],
                        "detail": f"{relative}: {heading}: fields={ordered}",
                    },
                )

            if any(len(fields[label]) != 1 for label in OWNER_FIELDS):
                continue

            classification_text = fields["Classification"][0][1]
            classifications = OWNER_CLASSIFICATION.findall(classification_text)

            if len(classifications) != 1:
                findings.append(
                    {
                        "kind": "owner_classification_invalid",
                        "rule_id": identifiers[0],
                        "detail": (
                            f"{relative}: {heading}: {classification_text}"
                        ),
                    },
                )

                continue

            for rule_id in identifiers:
                inventory.setdefault(rule_id, []).append(
                    OwnerRule(
                        rule_id=rule_id,
                        normative_document=relative,
                        normative_section=heading,
                        owner_classification=classifications[0],
                        owner_scope=fields["Scope"][0][1],
                        automation=fields["Automation"][0][1],
                        semantic_remainder=fields["Semantic remainder"][0][1],
                        exceptions=fields["Exceptions"][0][1],
                    ),
                )

    for rule_id, owners in sorted(inventory.items()):
        if len(owners) > 1:
            findings.append(
                {
                    "kind": "duplicate_owner",
                    "rule_id": rule_id,
                    "detail": f"owners={owners}",
                },
            )

    return (
        inventory,
        findings,
        [document.relative_to(root).as_posix() for document in documents],
    )


def _emitted_rule_ids(root: Path) -> dict[str, list[str]]:
    """
    Return stable IDs emitted by Python audit implementations.

    Parameters
    ----------
    root : Path
        Repository root containing the Python audit package.

    Returns
    -------
    emitted : dict[str, list[str]]
        Source paths indexed by each statically recoverable rule identifier.
    """

    emitted: dict[str, list[str]] = {}
    audit_root = root / "dev/audit"
    if not audit_root.is_dir():
        return emitted

    for path in sorted(audit_root.glob("*.py")):
        try:
            tree = ast.parse(path.read_text(encoding="utf-8"))
        except SyntaxError:
            continue

        relative = path.relative_to(root).as_posix()

        for node in ast.walk(tree):
            values: list[ast.AST] = []

            if isinstance(node, ast.Assign):
                names = [
                    target.id
                    for target in node.targets
                    if isinstance(target, ast.Name)
                ]

                if any(
                    name == "RULE_ID" or name.startswith("RULE_")
                    for name in names
                ):
                    values.append(node.value)
            elif isinstance(node, ast.AnnAssign):
                if (
                    isinstance(node.target, ast.Name)
                    and (
                        node.target.id == "rule_id"
                        or node.target.id == "RULE_ID"
                        or node.target.id.startswith("RULE_")
                    )
                    and node.value is not None
                ):
                    values.append(node.value)
            elif isinstance(node, ast.Dict):
                values.extend(
                    value
                    for key, value in zip(node.keys, node.values, strict=True)
                    if isinstance(key, ast.Constant) and key.value == "rule_id"
                )
            elif isinstance(node, ast.Call):
                function = node.func
                name = (
                    function.id
                    if isinstance(function, ast.Name)
                    else function.attr
                    if isinstance(function, ast.Attribute)
                    else ""
                )

                if name in {"Advisory", "Finding"}:
                    values.extend(node.args)

            for value in values:
                if (
                    isinstance(value, ast.Constant)
                    and isinstance(value.value, str)
                    and EMITTED_RULE_ID.fullmatch(value.value)
                ):
                    location = f"{relative}:{getattr(value, 'lineno', 1)}"
                    emitted.setdefault(value.value, []).append(location)

    return emitted


def _claims_execution(automation: str) -> bool:
    """
    Return whether owner automation text claims executable evidence.
    """

    lowered = automation.lower()
    if any(
        phrase in lowered
        for phrase in (
            "no dedicated registry entry",
            "no dedicated checker",
            "no registry entry",
            "no checker can",
            "no automation",
        )
    ):
        return False

    return bool(
        re.search(r"`[^`]+\.(?:py|sh)`", automation)
        or re.search(
            r"\b(checker|contract|audit|formatter|test suite|unit tests?)\b",
            lowered,
        ),
    )


def _registry_manifest(root: Path) -> dict[str, object]:
    """
    Return the complete static owner and execution reconciliation.

    Parameters
    ----------
    root : Path
        Repository root used to resolve maintained paths.

    Returns
    -------
    manifest : dict[str, object]
        Registry manifest and reconciliation findings.
    """

    config = root / "dev/config/rules.toml"
    data = tomllib.loads(config.read_text(encoding="utf-8"))
    rules = data.get("rule", [])

    inventory, findings, documents = _owner_inventory(root, rules)

    if data.get("schema_version") != 2:
        found_version = data.get("schema_version")
        findings.append(
            {
                "kind": "registry_schema_version",
                "rule_id": "<registry>",
                "detail": f"expected schema_version 2, found {found_version}",
            },
        )

    counts = Counter(rule.get("rule_id", "") for rule in rules)

    for rule_id, count in sorted(counts.items()):
        if not rule_id or count > 1:
            findings.append(
                {
                    "kind": "duplicate_ownership",
                    "rule_id": rule_id or "<missing>",
                    "detail": f"registry count is {count}",
                },
            )

    registry_ids = {str(rule.get("rule_id", "")) for rule in rules}

    for rule_id, owners in sorted(inventory.items()):
        if len(owners) != 1 or rule_id in registry_ids:
            continue

        owner = owners[0]

        if owner.owner_classification == "semantic-only":
            review_text = (
                f"{owner.automation} {owner.semantic_remainder}".lower()
            )

            if not any(
                token in review_text
                for token in (
                    "review",
                    "decide",
                    "determine",
                    "classify",
                    "choose",
                )
            ):
                findings.append(
                    {
                        "kind": "semantic_only_without_review_procedure",
                        "rule_id": rule_id,
                        "detail": (
                            f"{owner.normative_document}: "
                            f"{owner.normative_section}"
                        ),
                    },
                )
        elif owner.owner_classification == "deterministic":
            findings.append(
                {
                    "kind": "implementation_gap",
                    "rule_id": rule_id,
                    "detail": (
                        "deterministic owner has no executable registry entry"
                    ),
                },
            )
        elif _claims_execution(owner.automation):
            findings.append(
                {
                    "kind": "owner_coverage_gap",
                    "rule_id": rule_id,
                    "detail": (
                        "owner names executable evidence but has no registry "
                        "entry"
                    ),
                },
            )

    for rule in rules:
        rule_id = str(rule.get("rule_id", "<missing>"))
        document = root / str(rule.get("normative_document", ""))
        section = str(rule.get("normative_section", ""))

        if not document.is_file():
            findings.append(
                {
                    "kind": "checker_only_rule",
                    "rule_id": rule_id,
                    "detail": f"missing normative document: {document}",
                },
            )
        elif section not in HEADING.findall(
            document.read_text(encoding="utf-8"),
        ):
            findings.append(
                {
                    "kind": "missing_section",
                    "rule_id": rule_id,
                    "detail": f"missing heading: {section}",
                },
            )
        elif rule_id not in RULE_ID.findall(section):
            findings.append(
                {
                    "kind": "section_does_not_own_rule",
                    "rule_id": rule_id,
                    "detail": f"heading does not contain ID: {section}",
                },
            )

        owners = inventory.get(rule_id, [])
        expected = (str(rule.get("normative_document", "")), section)

        if not owners:
            findings.append(
                {
                    "kind": "checker_only_rule",
                    "rule_id": rule_id,
                    "detail": "registry ID has no maintained owner",
                },
            )
        elif not any(
            (owner.normative_document, owner.normative_section) == expected
            for owner in owners
        ):
            findings.append(
                {
                    "kind": "owner_mismatch",
                    "rule_id": rule_id,
                    "detail": f"registry={expected}; owners={owners}",
                },
            )

        required_fields = (
            "owner_classification",
            "execution_role",
            "coverage_relation",
            "covered_scope",
            "remaining_scope",
        )

        for field in required_fields:
            value = rule.get(field)

            if not isinstance(value, str) or not value.strip():
                findings.append(
                    {
                        "kind": "registry_field_missing",
                        "rule_id": rule_id,
                        "detail": field,
                    },
                )

        owner_classification = str(rule.get("owner_classification", ""))
        execution_role = str(rule.get("execution_role", ""))

        coverage_relation = str(rule.get("coverage_relation", ""))
        covered_scope = str(rule.get("covered_scope", "")).strip()
        remaining_scope = str(rule.get("remaining_scope", "")).strip()

        if remaining_scope == GENERIC_REMAINING_SCOPE:
            findings.append(
                {
                    "kind": "generic_remaining_scope",
                    "rule_id": rule_id,
                    "detail": remaining_scope,
                },
            )

        if owner_classification not in OWNER_CLASSIFICATIONS:
            findings.append(
                {
                    "kind": "registry_owner_classification_invalid",
                    "rule_id": rule_id,
                    "detail": owner_classification,
                },
            )
        elif (
            len(owners) == 1
            and owner_classification != owners[0].owner_classification
        ):
            findings.append(
                {
                    "kind": "owner_classification_mismatch",
                    "rule_id": rule_id,
                    "detail": (
                        f"registry={owner_classification}; "
                        f"owner={owners[0].owner_classification}"
                    ),
                },
            )

        if execution_role not in EXECUTION_ROLES:
            findings.append(
                {
                    "kind": "execution_role_invalid",
                    "rule_id": rule_id,
                    "detail": execution_role,
                },
            )

        if coverage_relation not in COVERAGE_RELATIONS:
            findings.append(
                {
                    "kind": "coverage_relation_invalid",
                    "rule_id": rule_id,
                    "detail": coverage_relation,
                },
            )

        if coverage_relation == "exact":
            if owner_classification != "deterministic":
                findings.append(
                    {
                        "kind": "exact_requires_deterministic_owner",
                        "rule_id": rule_id,
                        "detail": owner_classification,
                    },
                )

            if execution_role != "checker":
                findings.append(
                    {
                        "kind": "exact_requires_checker",
                        "rule_id": rule_id,
                        "detail": execution_role,
                    },
                )

            if remaining_scope != "None":
                findings.append(
                    {
                        "kind": "exact_has_remaining_scope",
                        "rule_id": rule_id,
                        "detail": remaining_scope,
                    },
                )
        elif (
            coverage_relation in {"subset", "independent"}
            and remaining_scope == "None"
        ):
            findings.append(
                {
                    "kind": "partial_coverage_without_remainder",
                    "rule_id": rule_id,
                    "detail": coverage_relation,
                },
            )

        if (coverage_relation == "independent") != (
            execution_role == "independent_evidence"
        ):
            findings.append(
                {
                    "kind": "independent_role_relation_mismatch",
                    "rule_id": rule_id,
                    "detail": (
                        f"role={execution_role}; relation={coverage_relation}"
                    ),
                },
            )

        if (
            execution_role == "evidence_producer"
            and coverage_relation == "exact"
        ):
            findings.append(
                {
                    "kind": "evidence_producer_claims_exact",
                    "rule_id": rule_id,
                    "detail": covered_scope,
                },
            )

        if (
            owner_classification == "semantic-only"
            and coverage_relation != "independent"
        ):
            findings.append(
                {
                    "kind": "semantic_only_execution_not_independent",
                    "rule_id": rule_id,
                    "detail": coverage_relation,
                },
            )

        checker = str(rule.get("source_checker", ""))

        if _is_path(checker) and not (root / checker).exists():
            findings.append(
                {
                    "kind": "checker_only_rule",
                    "rule_id": rule_id,
                    "detail": f"missing checker: {checker}",
                },
            )

        fixture = str(rule.get("parity_test", ""))

        if execution_role == "checker" and (
            not fixture
            or (_is_path(fixture) and not (root / fixture).exists())
        ):
            findings.append(
                {
                    "kind": "checker_without_fixtures",
                    "rule_id": rule_id,
                    "detail": fixture or "no parity/fixture path registered",
                },
            )

        gate = str(rule.get("gate", ""))

        if gate and _is_path(gate) and not (root / gate).exists():
            findings.append(
                {
                    "kind": "missing_registered_gate",
                    "rule_id": rule_id,
                    "detail": f"missing gate: {gate}",
                },
            )

        for exception in rule.get("current_exclusions_or_allowlists", []):
            if "owner=" not in str(exception):
                findings.append(
                    {
                        "kind": "unowned_exception",
                        "rule_id": rule_id,
                        "detail": str(exception),
                    },
                )

    diagnostic_owners: dict[str, str] = {}
    migrations: dict[str, str] = {}

    for rule in rules:
        rule_id = str(rule.get("rule_id", "<missing>"))

        for diagnostic in rule.get("diagnostic_ids", []):
            diagnostic_id = str(diagnostic)
            prior = diagnostic_owners.setdefault(diagnostic_id, rule_id)

            if prior != rule_id:
                findings.append(
                    {
                        "kind": "duplicate_diagnostic",
                        "rule_id": diagnostic_id,
                        "detail": f"mapped to {prior} and {rule_id}",
                    },
                )

            if diagnostic_id in inventory:
                findings.append(
                    {
                        "kind": "diagnostic_is_owner",
                        "rule_id": diagnostic_id,
                        "detail": f"also mapped as a facet of {rule_id}",
                    },
                )

        for retired in rule.get("migrates_from", []):
            retired_id = str(retired)
            prior = migrations.setdefault(retired_id, rule_id)

            if prior != rule_id:
                findings.append(
                    {
                        "kind": "duplicate_migration",
                        "rule_id": retired_id,
                        "detail": f"mapped to {prior} and {rule_id}",
                    },
                )

            if retired_id in inventory:
                findings.append(
                    {
                        "kind": "retired_id_still_owned",
                        "rule_id": retired_id,
                        "detail": f"migrates to {rule_id}",
                    },
                )

    known_ids = set(counts) | set(diagnostic_owners) | set(migrations)

    for emitted_id, locations in sorted(_emitted_rule_ids(root).items()):
        if emitted_id not in known_ids:
            findings.append(
                {
                    "kind": "unowned_checker_id",
                    "rule_id": emitted_id,
                    "detail": f"emitted at {sorted(set(locations))}",
                },
            )

    findings = sorted(
        findings,
        key=lambda item: (item["rule_id"], item["kind"], item["detail"]),
    )

    execution_by_owner: dict[str, list[dict[str, object]]] = {}
    manifest_fields = (
        "owner_classification",
        "execution_role",
        "coverage_relation",
        "covered_scope",
        "remaining_scope",
        "source_checker",
        "execution_kind",
        "parity_test",
        "gate",
        "applicable_paths",
        "blocking",
    )

    for rule in rules:
        rule_id = str(rule.get("rule_id", ""))
        execution_by_owner.setdefault(rule_id, []).append(
            {field: rule.get(field) for field in manifest_fields},
        )

    owners_output: list[dict[str, object]] = []

    for rule_id, owners in sorted(inventory.items()):
        if len(owners) != 1:
            continue

        owner = owners[0]
        executions = execution_by_owner.get(rule_id, [])
        owners_output.append(
            {
                **dataclasses.asdict(owner),
                "coverage_status": "registered"
                if executions
                else "review_only",
                "executions": executions,
            },
        )

    owner_classifications = Counter(
        owner["owner_classification"] for owner in owners_output
    )
    execution_roles = Counter(
        str(rule.get("execution_role", "<missing>")) for rule in rules
    )

    coverage_relations = Counter(
        str(rule.get("coverage_relation", "<missing>")) for rule in rules
    )

    partition = {
        "owner_total": len(owners_output),
        "registered_owner_total": sum(
            owner["coverage_status"] == "registered" for owner in owners_output
        ),
        "review_only_owner_total": sum(
            owner["coverage_status"] == "review_only"
            for owner in owners_output
        ),
        "owner_classification": dict(sorted(owner_classifications.items())),
        "execution_role": dict(sorted(execution_roles.items())),
        "coverage_relation": dict(sorted(coverage_relations.items())),
    }

    return {
        "schema_version": 1,
        "maintained_documents": documents,
        "partition": partition,
        "owners": owners_output,
        "finding_count": len(findings),
        "findings": findings,
    }


def audit_registry(root: Path) -> list[dict[str, str]]:
    """
    Return deterministic registry findings below *root*.
    """

    return list(_registry_manifest(root)["findings"])


def main(argv: list[str] | None = None) -> int:
    """
    Run the registry audit and emit JSON.

    Parameters
    ----------
    argv : list[str] | None
        Explicit arguments, or None to read the process arguments.

    Returns
    -------
    status : int
        Zero when registry ownership is valid and one otherwise.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--report-only", action="store_true")
    parser.add_argument("--manifest-output", type=Path)
    args = parser.parse_args(argv)
    manifest = _registry_manifest(args.root.resolve())
    findings = manifest["findings"]

    if args.manifest_output is not None:
        output = args.manifest_output.resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    print(
        json.dumps(
            {
                "finding_count": len(findings),
                "partition": manifest["partition"],
                "findings": findings,
            },
            indent=2,
            sort_keys=True,
        ),
    )

    return 0 if args.report_only or not findings else 1


if __name__ == "__main__":
    raise SystemExit(main())
