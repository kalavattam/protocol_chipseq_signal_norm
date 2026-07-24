#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_standards_registry.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


"""Tests for complete standards-owner and execution-registry reconciliation."""

from __future__ import annotations

from pathlib import Path

from dev.audit.standards_registry import _registry_manifest, audit_registry


def owner_section(
    rule_id: str,
    *,
    classification: str = "deterministic",
    heading: str | None = None,
    automation: str = "`tests/check.py` checks the bounded rule.",
) -> str:
    """Return one canonical owner section."""

    title = heading or f"Owned section (`{rule_id}`)"
    return f"""
## {title}
**Classification:** `{classification}`.

**Scope:** Fixture scope.

Fixture requirement.

**Automation:** {automation}

**Semantic remainder:** Review the fixture when needed.

**Exceptions:** None.
"""


def registry_entry(
    rule_id: str = "ONE.RULE",
    *,
    section: str | None = None,
    owner_classification: str = "deterministic",
    execution_role: str = "checker",
    coverage_relation: str = "subset",
    remaining_scope: str = "Unproved fixture behavior remains.",
    source_checker: str = "tests/check.py",
    parity_test: str = "tests/check.py",
) -> str:
    """Return one schema-v2 execution entry."""

    heading = section or f"Owned section (`{rule_id}`)"
    return f"""
[[rule]]
rule_id = "{rule_id}"
normative_document = "docs/standards/standard.md"
normative_section = "{heading}"
owner_classification = "{owner_classification}"
source_checker = "{source_checker}"
execution_kind = "original"
execution_role = "{execution_role}"
coverage_relation = "{coverage_relation}"
covered_scope = "Fixture-defined execution scope."
remaining_scope = "{remaining_scope}"
parity_test = "{parity_test}"
scope = "repository"
command = ["python3", "tests/check.py"]
applicable_paths = ["docs/**"]
blocking = false
semantic_review = false
current_exclusions_or_allowlists = []
"""


def write_repository(
    root: Path,
    rules: str,
    standard: str,
    *,
    checker: bool = True,
) -> None:
    """Write a minimal indexed standards repository."""

    (root / "dev/config").mkdir(parents=True, exist_ok=True)
    (root / "docs/standards").mkdir(parents=True, exist_ok=True)
    (root / "docs/standards/standard.md").write_text(
        standard,
        encoding="utf-8",
    )
    (root / "docs/standards/README.md").write_text(
        "| Standard | Ownership |\n"
        "| :--- | :--- |\n"
        "| [`standard.md`](standard.md) | Fixture |\n",
        encoding="utf-8",
    )
    (root / "dev/config/rules.toml").write_text(
        "schema_version = 2\n" + rules,
        encoding="utf-8",
    )
    if checker:
        (root / "tests").mkdir(exist_ok=True)
        (root / "tests/check.py").write_text("pass\n", encoding="utf-8")


def test_accepts_complete_subset_trace(tmp_path: Path) -> None:
    """One owner and one truthful bounded execution reconcile cleanly."""

    write_repository(
        tmp_path,
        registry_entry(),
        "# Standard\n" + owner_section("ONE.RULE"),
    )
    assert audit_registry(tmp_path) == []


def test_detects_complete_owner_inventory_gap(tmp_path: Path) -> None:
    """A deterministic unregistered owner cannot disappear from review."""

    write_repository(
        tmp_path,
        "",
        "# Standard\n" + owner_section("ONE.RULE"),
        checker=False,
    )
    findings = audit_registry(tmp_path)
    assert any(item["kind"] == "implementation_gap" for item in findings)
    manifest = _registry_manifest(tmp_path)
    assert manifest["owners"][0]["coverage_status"] == "review_only"


def test_accepts_semantic_only_review_without_execution(tmp_path: Path) -> None:
    """Semantic-only owners remain visible without dummy registry rows."""

    write_repository(
        tmp_path,
        "",
        "# Standard\n"
        + owner_section(
            "ONE.RULE",
            classification="semantic-only",
            automation="No checker can decide this rule; review the evidence.",
        ),
        checker=False,
    )
    assert audit_registry(tmp_path) == []


def test_detects_duplicate_owner_and_registry_rows(tmp_path: Path) -> None:
    """Duplicate owner and execution identities fail closed."""

    write_repository(
        tmp_path,
        registry_entry() + registry_entry(),
        "# Standard\n"
        + owner_section("ONE.RULE")
        + owner_section("ONE.RULE", heading="Other owner (`ONE.RULE`)"),
    )
    kinds = {item["kind"] for item in audit_registry(tmp_path)}
    assert {"duplicate_owner", "duplicate_ownership"} <= kinds


def test_detects_owner_classification_mismatch(tmp_path: Path) -> None:
    """Execution metadata must repeat the canonical owner classification."""

    write_repository(
        tmp_path,
        registry_entry(owner_classification="advisory"),
        "# Standard\n" + owner_section("ONE.RULE"),
    )
    assert any(
        item["kind"] == "owner_classification_mismatch"
        for item in audit_registry(tmp_path)
    )


def test_exact_requires_checker_deterministic_owner_and_no_remainder(
    tmp_path: Path,
) -> None:
    """Exact coverage has a closed, machine-validated meaning."""

    write_repository(
        tmp_path,
        registry_entry(
            owner_classification="advisory",
            execution_role="evidence_producer",
            coverage_relation="exact",
            remaining_scope="Still unproved.",
        ),
        "# Standard\n"
        + owner_section("ONE.RULE", classification="advisory"),
    )
    kinds = {item["kind"] for item in audit_registry(tmp_path)}
    assert {
        "exact_requires_deterministic_owner",
        "exact_requires_checker",
        "exact_has_remaining_scope",
        "evidence_producer_claims_exact",
    } <= kinds


def test_independent_relation_requires_independent_role(tmp_path: Path) -> None:
    """Independent evidence is distinct from subset implementation."""

    write_repository(
        tmp_path,
        registry_entry(
            execution_role="checker",
            coverage_relation="independent",
        ),
        "# Standard\n" + owner_section("ONE.RULE"),
    )
    assert any(
        item["kind"] == "independent_role_relation_mismatch"
        for item in audit_registry(tmp_path)
    )


def test_detects_unowned_checker_id(tmp_path: Path) -> None:
    """A checker-emitted ID must be owned, mapped, or migrated."""

    write_repository(
        tmp_path,
        registry_entry(),
        "# Standard\n" + owner_section("ONE.RULE"),
    )
    (tmp_path / "dev/audit").mkdir(parents=True)
    (tmp_path / "dev/audit/check.py").write_text(
        'RULE_ID = "CHECK.UNOWNED"\n',
        encoding="utf-8",
    )
    findings = audit_registry(tmp_path)
    assert any(
        item["kind"] == "unowned_checker_id"
        and item["rule_id"] == "CHECK.UNOWNED"
        for item in findings
    )


def test_manifest_partition_is_derived_from_finished_owners(tmp_path: Path) -> None:
    """The report calculates its partition instead of enforcing a fixed count."""

    write_repository(
        tmp_path,
        registry_entry(),
        "# Standard\n"
        + owner_section("ONE.RULE")
        + owner_section(
            "REVIEW.RULE",
            classification="semantic-only",
            automation="No checker can decide this rule; review the evidence.",
        ),
    )
    partition = _registry_manifest(tmp_path)["partition"]
    assert partition["owner_total"] == 2
    assert partition["registered_owner_total"] == 1
    assert partition["review_only_owner_total"] == 1


def test_index_missing_target_fails_closed(tmp_path: Path) -> None:
    """An indexed but absent standard remains a reconciliation error."""

    write_repository(
        tmp_path,
        registry_entry(),
        "# Standard\n" + owner_section("ONE.RULE"),
    )
    with (tmp_path / "docs/standards/README.md").open("a", encoding="utf-8") as handle:
        handle.write("| [`missing.md`](missing.md) | Missing |\n")
    assert any(
        item["kind"] == "missing_index_target"
        for item in audit_registry(tmp_path)
    )


def test_duplicate_index_entry_fails_closed(tmp_path: Path) -> None:
    """Duplicate maintained-standard index rows cannot be deduplicated silently."""

    write_repository(
        tmp_path,
        registry_entry(),
        "# Standard\n" + owner_section("ONE.RULE"),
    )
    with (tmp_path / "docs/standards/README.md").open("a", encoding="utf-8") as handle:
        handle.write("| [`standard.md`](standard.md) | Duplicate |\n")
    assert any(
        item["kind"] == "duplicate_index_entry"
        for item in audit_registry(tmp_path)
    )


def test_unindexed_nonbackup_standard_fails_closed(tmp_path: Path) -> None:
    """A non-backup owner document cannot disappear outside the index."""

    write_repository(
        tmp_path,
        registry_entry(),
        "# Standard\n" + owner_section("ONE.RULE"),
    )
    (tmp_path / "docs/standards/omitted.md").write_text(
        "# Omitted\n" + owner_section("OMITTED.RULE"),
        encoding="utf-8",
    )
    assert any(
        item["kind"] == "unindexed_standard"
        for item in audit_registry(tmp_path)
    )


def test_fenced_example_ids_are_not_enumerated_as_owners(tmp_path: Path) -> None:
    """Fence-aware extraction ignores example H2 headings and IDs."""

    write_repository(
        tmp_path,
        registry_entry(),
        "# Standard\n"
        + owner_section("ONE.RULE")
        + "\n```markdown\n## Example (`EXAMPLE.RULE`)\n```\n",
    )
    manifest = _registry_manifest(tmp_path)
    assert manifest["finding_count"] == 0
    assert [owner["rule_id"] for owner in manifest["owners"]] == ["ONE.RULE"]


def test_stock_remaining_scope_placeholder_is_rejected(tmp_path: Path) -> None:
    """A nonempty stock sentence is not precise coverage metadata."""

    write_repository(
        tmp_path,
        registry_entry(
            remaining_scope=(
                "Owner scope outside the registered execution remains unproved."
            )
        ),
        "# Standard\n" + owner_section("ONE.RULE"),
    )
    assert any(
        item["kind"] == "generic_remaining_scope"
        for item in audit_registry(tmp_path)
    )


def test_exception_allowlist_requires_owner_marker(tmp_path: Path) -> None:
    """The registered exception subset accepts owner= and rejects its absence."""

    owned = registry_entry().replace(
        "current_exclusions_or_allowlists = []",
        'current_exclusions_or_allowlists = ["owner=ONE.RULE; fixture"]',
    )
    write_repository(
        tmp_path,
        owned,
        "# Standard\n" + owner_section("ONE.RULE"),
    )
    assert audit_registry(tmp_path) == []
    unowned = owned.replace("owner=ONE.RULE; ", "")
    write_repository(
        tmp_path,
        unowned,
        "# Standard\n" + owner_section("ONE.RULE"),
    )
    assert any(
        item["kind"] == "unowned_exception"
        for item in audit_registry(tmp_path)
    )
