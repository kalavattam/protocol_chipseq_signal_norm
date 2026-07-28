#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_python_source_evidence.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


"""
Test stable nonblocking Python source-policy evidence.
"""

from __future__ import annotations

import json
from pathlib import Path

import jsonschema
from dev.audit.python_source_evidence import (
    _load_config,
    _membership_fingerprint,
    _semantic_record_surfaces,
    produce,
)
from dev.audit.python_source_policy import report

ROOT = Path(__file__).resolve().parents[3]
CONFIG_DIR = ROOT / "dev" / "config" / "pilots"
SCHEMA_DIR = ROOT / "dev" / "schemas"
CONFIG_PATH = CONFIG_DIR / "python_source_policy.json"
SCHEMA_PATH = SCHEMA_DIR / "python_source_candidate.schema.json"
REVIEW_SCHEMA_PATH = SCHEMA_DIR / "python_semantic_review.schema.json"
REVIEW_MANIFEST_PATH = (
    ROOT
    / "artifacts"
    / "tests"
    / "python_source_migration_20260726_corrected"
    / "final"
    / "semantic_review_manifest.json"
)
PILOT_PATHS = [
    "src/protocol_chipseq_signal_norm/cli/parse_metadata_siqchip.py",
    "src/protocol_chipseq_signal_norm/utilities/utils_io.py",
    "dev/tools/markdown_format.py",
    "tests/unit/metadata/test_parse_siqchip.py",
]


def configuration() -> dict[str, object]:
    """
    Return one small evidence configuration.
    """

    return {
        "schema_version": 2,
        "pilot_id": "unit",
        "threshold_definitions": {
            "X": "Multiline call physical-line span.",
            "Y_assignment_span": "Multiline assignment physical-line span.",
            "Y_simple_run_length": "Simple assignment run statement count.",
            "Z": "Attached compact-transfer preparatory physical lines.",
        },
        "pilot_paths": ["src/sample.py"],
        "x_multiline_call_span_thresholds": [5, 6, 7],
        "y_multiline_assignment_span_thresholds": [4, 5, 6],
        "y_simple_assignment_run_statement_thresholds": [4, 5, 6],
        "z_compact_transfer_preparatory_line_limit": 1,
        "line_length": 79,
        "naming_length_thresholds": {
            "production_callable": [40, 48, 56],
            "test_function": [56, 64, 72],
            "local_or_attribute": [24, 32, 40],
            "type": [32, 40, 48],
        },
        "disposition": "unresolved",
        "semantic_candidates_are_nonblocking": True,
    }


def write_source(root: Path) -> str:
    """
    Write one synthetic maintained source and return its relative path.
    """

    path = root / "src" / "sample.py"
    path.parent.mkdir(parents=True)
    path.write_text(
        '''
"""
Provide stable source-policy evidence.
"""


def build(flag: bool, other: bool) -> tuple[object, int]:
    """
    Build one multiline value and one guarded result.
    """

    value = configure(
        "one",
        "two",
        "three",
        "four",
    )
    if flag:
        first = 1
        second = 2
        return value, first + second
    if other:
        result = 1
        return value, result
    return value, 0
''',
        encoding="utf-8",
    )

    return "src/sample.py"


def activate_review(
    root: Path,
    paths: list[str],
    config: dict[str, object],
) -> dict[str, object]:
    """
    Write explicit unit-fixture decisions for every current semantic record.
    """

    unresolved = produce(root, paths, config)
    groups = []
    surfaces = _semantic_record_surfaces(
        unresolved["inventories"],
        unresolved["candidates"],
    )

    for owner, records in sorted(surfaces.items()):
        buckets: dict[str, list[tuple[dict[str, object], str]]] = {}

        for record, _ in records:
            if (owner.startswith("docstring_") and not record["present"]) or (
                owner == "source_headers" and not record["applicable"]
            ):
                disposition = "omitted_by_role"
            elif owner == "source_headers" and record["previously_corrected"]:
                disposition = "refactored"
            elif owner == "raw_width_candidates":
                disposition = "retained_indivisible"
            elif owner == "reusable_error_boundaries":
                disposition = "exception"
            elif owner == "completed_rename_decisions":
                disposition = "refactored"
            else:
                disposition = "retained_coherent"

            key = (
                record["review_keys"][owner]
                if owner in {"docstring_applicability", "docstring_content"}
                else record["review_key"]
            )
            buckets.setdefault(disposition, []).append((record, key))

        for disposition, bucket in sorted(buckets.items()):
            members = sorted(key for _, key in bucket)
            details = ", ".join(
                f"{record.get('path')}:{record.get('line', 1)}"
                for record, _ in bucket
            )
            group_id = f"unit_{owner}_{disposition}"
            groups.append(
                {
                    "group_id": group_id,
                    "owner": owner,
                    "disposition": disposition,
                    "rationale": (
                        f"Every member of {group_id} was reviewed against the "
                        f"bounded fixture source and assertions: {details}."
                    ),
                    "reviewed_facts": {
                        "bounded_fixture_members": details,
                    },
                    "members": members,
                    "membership_fingerprint": _membership_fingerprint(
                        members,
                    ),
                },
            )

    manifest_path = root / "semantic_review.json"
    manifest_path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "review_protocol": {
                    "status": "reviewed",
                    "hash_role": (
                        "Hashes invalidate stale explicit decisions only."
                    ),
                    "membership_rule": (
                        "Every group names exact reviewed fixture members."
                    ),
                    "review_surfaces": ["bounded unit fixture"],
                },
                "reviewed_source_fingerprint": unresolved[
                    "source_fingerprint"
                ],
                "decision_groups": groups,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    config["semantic_review_manifest"] = "semantic_review.json"

    return config


def test_pilot_configuration_has_the_exact_approved_boundary() -> None:
    """
    Keep one source of truth for pilot paths and adjacent values.
    """

    config = _load_config(CONFIG_PATH)

    assert config["pilot_paths"] == PILOT_PATHS
    assert config["x_multiline_call_span_thresholds"] == [5, 6, 7]
    assert config["y_multiline_assignment_span_thresholds"] == [4, 5, 6]
    assert config["y_simple_assignment_run_statement_thresholds"] == [
        4,
        5,
        6,
    ]
    assert set(config["threshold_definitions"]) == {
        "X",
        "Y_assignment_span",
        "Y_simple_run_length",
        "Z",
    }
    assert config["z_compact_transfer_preparatory_line_limit"] == 1
    assert config["line_length"] == 79


def test_explicit_review_manifest_matches_its_schema() -> None:
    """
    Keep decision membership and hash semantics structurally inspectable.
    """

    schema = json.loads(REVIEW_SCHEMA_PATH.read_text(encoding="utf-8"))
    manifest = json.loads(REVIEW_MANIFEST_PATH.read_text(encoding="utf-8"))

    jsonschema.validate(manifest, schema)


def test_configured_pilot_has_no_unresolved_raw_width_candidates() -> None:
    """
    Require project-level physical-width closure independent of Ruff.
    """

    config = _load_config(CONFIG_PATH)

    payload = report(ROOT, config["pilot_paths"])

    assert payload["finding_count"] == 0
    assert {
        fact["path"]: fact["line_length"]
        for fact in payload["facts"]
        if fact["line_length"].get("over_79_raw", 0)
        or fact["line_length"].get("review_candidates", 0)
    } == {}


def test_strict_x_y_and_z_boundaries_are_not_shifted(
    tmp_path: Path,
) -> None:
    """
    Treat X and Y as strict span thresholds and Z as selected at one.
    """

    path = write_source(tmp_path)

    payload = produce(tmp_path, [path], configuration())

    assert payload["inventories"]["x_distributions"] == {
        "5": 1,
        "6": 0,
        "7": 0,
    }
    assert payload["inventories"]["y_distributions"] == {
        "multiline_assignment_span": {
            "4": 1,
            "5": 1,
            "6": 0,
        },
        "simple_assignment_run_length": {
            "4": 0,
            "5": 0,
            "6": 0,
        },
    }
    assert payload["inventories"]["z_distribution"] == {"1": 1}


def test_z_excludes_a_transfer_already_separated_from_preparation(
    tmp_path: Path,
) -> None:
    """
    Do not report a transfer whose immediate source boundary is blank.
    """

    path = tmp_path / "src" / "sample.py"
    path.parent.mkdir(parents=True)
    path.write_text(
        """
def build(flag: bool) -> int:
    if flag:
        first = 1
        second = 2

        return first + second
    return 0
""".removeprefix("\n"),
        encoding="utf-8",
    )

    payload = produce(tmp_path, ["src/sample.py"], configuration())

    assert payload["inventories"]["z_distribution"] == {"1": 0}

    transfer = next(
        record
        for record in payload["inventories"]["transfers"]
        if record["line"] == 6
    )

    assert transfer["blank_line_state"] == "separated"
    assert transfer["predecessor_kind"] == "Assign"
    assert transfer["physical_preparation"] is False


def test_density_is_separate_and_structurally_descriptive(
    tmp_path: Path,
) -> None:
    """
    Record statement, nesting, construct, and physical-region evidence.
    """

    path = write_source(tmp_path)
    density = produce(
        tmp_path,
        [path],
        configuration(),
    )["inventories"]["paragraph_density"]

    assert density["region_count"] > 0
    assert density["physical_span_histogram"]
    assert density["statement_count_histogram"]

    record = density["largest_regions"][0]

    assert {
        "physical_span",
        "logical_line_count",
        "statement_count",
        "nesting",
        "construct_counts",
        "comment_line_count",
    } <= set(record)
    assert all(
        "FunctionDef" not in region["construct_counts"]
        for region in density["largest_regions"]
    )


def test_candidate_schema_signatures_and_output_are_stable(
    tmp_path: Path,
) -> None:
    """
    Validate schema shape, unique signatures, and producer idempotence.
    """

    path = write_source(tmp_path)
    first = produce(tmp_path, [path], configuration())
    second = produce(tmp_path, [path], configuration())

    assert first == second

    schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))

    jsonschema.Draft202012Validator(schema).validate(first)

    signatures = [candidate["signature"] for candidate in first["candidates"]]

    assert len(signatures) == len(set(signatures))

    repeated = first["inventories"]["repeated_statements"]
    transfers = first["inventories"]["transfers"]

    assert len({record["signature"] for record in repeated}) == len(repeated)
    assert len({record["signature"] for record in transfers}) == len(transfers)


def test_naming_evidence_separates_roles_and_conservative_renames(
    tmp_path: Path,
) -> None:
    """
    Keep constants out of locals and leave public surfaces unresolved.
    """

    path = write_source(tmp_path)
    naming = produce(
        tmp_path,
        [path],
        configuration(),
    )["inventories"]["naming"]

    assert naming["length_distributions"]
    assert naming["adjacent_threshold_counts"]
    assert naming["reference_scope"] == "all maintained Python source"
    assert "and" not in naming["abbreviation_candidate_distribution"]
    assert all(
        "direct_rename_status" in record for record in naming["records"]
    )


def test_repeated_runs_expose_every_member_and_split_at_blank_lines(
    tmp_path: Path,
) -> None:
    """
    Keep complete member evidence without bridging distinct source regions.
    """

    path = tmp_path / "src" / "sample.py"
    path.parent.mkdir(parents=True)
    path.write_text(
        """
def build() -> int:
    first = 1
    second = 2
    third = 3
    fourth = 4
    fifth = 5

    sixth = 6
    seventh = 7
    return first + second + third + fourth + fifth + sixth + seventh
""",
        encoding="utf-8",
    )

    payload = produce(tmp_path, ["src/sample.py"], configuration())
    runs = payload["inventories"]["repeated_statements"]

    assert [run["statement_count"] for run in runs] == [5, 2]
    assert [len(run["members"]) for run in runs] == [5, 2]
    assert all(
        member["disposition"] == "unresolved"
        for run in runs
        for member in run["members"]
    )

    simple_y = payload["inventories"]["y_distributions"][
        "simple_assignment_run_length"
    ]

    assert simple_y == {"4": 1, "5": 0, "6": 0}


def test_transfer_inventory_covers_all_suites_and_z_boundaries(
    tmp_path: Path,
) -> None:
    """
    Expose every transfer while applying Z only to eligible primary bodies.
    """

    path = tmp_path / "src" / "sample.py"
    path.parent.mkdir(parents=True)
    path.write_text(
        """
def generate(flag: bool, values: list[int]):
    if flag:
        result = 1
        return result
    else:
        raise ValueError
    for value in values:
        if value < 0:
            continue
        break
    try:
        yield 1
    except ValueError:
        raise
    finally:
        yield from values
    with open("value.txt") as handle:
        return handle.read()
    match values:
        case []:
            return 0
    return 1


def stop() -> None:
    first = 1
    second = 2

    sys.exit(first + second)
""",
        encoding="utf-8",
    )

    payload = produce(tmp_path, ["src/sample.py"], configuration())
    transfers = payload["inventories"]["transfers"]

    assert {
        "break",
        "continue",
        "process_exit",
        "raise",
        "return",
        "yield",
        "yield_from",
    } == {record["transfer_kind"] for record in transfers}
    assert {"body", "orelse", "finalbody"} <= {
        record["suite_field"] for record in transfers
    }
    assert {"ExceptHandler", "With", "match_case"} <= {
        record["suite_owner"] for record in transfers
    }

    compact = next(
        record
        for record in transfers
        if record["transfer_kind"] == "return"
        and record["predecessor_kind"] == "Assign"
    )

    assert compact["z_applies"] is True
    assert compact["preparation_physical_line_count"] == 1
    assert compact["z_exceeded"] is False

    process_exit = next(
        record
        for record in transfers
        if record["transfer_kind"] == "process_exit"
    )

    assert process_exit["blank_line_state"] == "separated"
    assert process_exit["substantive_preceding_phase"] is True
    assert process_exit["preceding_phase_statement_count"] == 2
    assert [
        member["kind"] for member in process_exit["preceding_phase_members"]
    ] == ["Assign", "Assign"]
    assert all(
        member["source_fingerprint"] == process_exit["source_fingerprint"]
        for member in process_exit["preceding_phase_members"]
    )


def test_matching_cohort_hash_cannot_resolve_semantic_records(
    tmp_path: Path,
) -> None:
    """
    Ignore the superseded false-green switch without explicit decisions.
    """

    path = write_source(tmp_path)
    config = configuration()
    baseline = produce(tmp_path, [path], config)
    config.update(
        {
            "semantic_review_complete": True,
            "semantic_review_active": True,
            "reviewed_cohort_fingerprint": baseline["source_fingerprint"],
        },
    )
    false_green_attempt = produce(tmp_path, [path], config)

    assert false_green_attempt["semantic_review"]["complete"] is False
    assert (
        false_green_attempt["semantic_review"]["explicit_decision_count"] == 0
    )
    assert all(
        counts["unresolved"] == counts["total"]
        for counts in false_green_attempt["semantic_owner_counts"].values()
        if counts["total"]
    )


def test_matching_manifest_hash_without_decisions_cannot_complete(
    tmp_path: Path,
) -> None:
    """
    Require explicit decisions even when reviewed source identity is current.
    """

    path = write_source(tmp_path)
    config = configuration()
    baseline = produce(tmp_path, [path], config)
    manifest = {
        "schema_version": 1,
        "review_protocol": {
            "status": "reviewed",
            "hash_role": "Hashes invalidate only.",
            "membership_rule": "Groups must name exact members.",
            "review_surfaces": ["bounded unit fixture"],
        },
        "reviewed_source_fingerprint": baseline["source_fingerprint"],
        "decision_groups": [],
    }
    (tmp_path / "semantic_review.json").write_text(
        json.dumps(manifest),
        encoding="utf-8",
    )
    config["semantic_review_manifest"] = "semantic_review.json"

    false_green_attempt = produce(tmp_path, [path], config)
    review = false_green_attempt["semantic_review"]

    assert review["source_fingerprint_matches"] is True
    assert review["explicit_decision_count"] == 0
    assert review["unresolved_record_count"] > 0
    assert review["complete"] is False


def test_role_aware_docstring_groups_resolve_every_member(
    tmp_path: Path,
) -> None:
    """
    Keep explicit membership for documented and behavior-named test roles.
    """

    path = tmp_path / "tests" / "test_sample.py"
    path.parent.mkdir(parents=True)
    path.write_text(
        """
def test_clear_behavior_name() -> None:
    assert True


def _fixture_value() -> int:
    return 1
""",
        encoding="utf-8",
    )
    config = activate_review(
        tmp_path,
        ["tests/test_sample.py"],
        configuration(),
    )

    inventory = produce(
        tmp_path,
        ["tests/test_sample.py"],
        config,
    )["inventories"]["docstrings"]

    assert inventory["disposition_counts"] == {
        "applicability_disposition": {
            "resolved": 3,
            "unresolved": 0,
        },
        "content_disposition": {
            "resolved": 3,
            "unresolved": 0,
        },
    }
    assert sum(
        group["member_count"] for group in inventory["review_groups"]
    ) == len(inventory["records"])
    assert all(group["members"] for group in inventory["review_groups"])


def test_explicit_group_membership_invalidates_after_source_change(
    tmp_path: Path,
) -> None:
    """
    Invalidate an entire explicit group when a source member changes.
    """

    path = write_source(tmp_path)
    config = activate_review(tmp_path, [path], configuration())
    reviewed = produce(tmp_path, [path], config)

    assert reviewed["semantic_review"]["complete"] is True
    assert reviewed["semantic_review"]["explicit_decision_count"] > 0
    assert all(
        counts["unresolved"] == 0
        for counts in reviewed["semantic_owner_counts"].values()
    )

    source = tmp_path / path
    source.write_text(
        source.read_text(encoding="utf-8") + "\n",
        encoding="utf-8",
    )
    stale = produce(tmp_path, [path], config)

    assert stale["semantic_review"]["complete"] is False
    assert stale["semantic_review"]["stale_group_count"] > 0
    assert stale["semantic_review"]["unresolved_record_count"] > 0
    assert "never create" in stale["semantic_review"]["hash_role"]


def test_raw_width_exceptions_have_candidate_specific_reasons(
    tmp_path: Path,
) -> None:
    """
    Retain only exact behavior-named test identifiers with local reasons.
    """

    path = tmp_path / "tests" / "test_width.py"
    path.parent.mkdir(parents=True)
    identifier = (
        "test_one_behavior_name_that_is_intentionally_long_for_discovery_"
        "and_diagnostics"
    )
    path.write_text(
        f"def {identifier}() -> None:\n    assert True\n",
        encoding="utf-8",
    )
    config = activate_review(
        tmp_path,
        ["tests/test_width.py"],
        configuration(),
    )
    raw_width = produce(
        tmp_path,
        ["tests/test_width.py"],
        config,
    )["inventories"]["raw_width"]
    disposition_reason = raw_width["records"][0]["disposition_reason"]

    assert raw_width["disposition_counts"]["unresolved"] == 0
    assert raw_width["records"][0]["disposition"] == "retained_indivisible"
    assert "tests/test_width.py:1" in disposition_reason


def test_source_header_review_separates_direct_and_import_only_modules(
    tmp_path: Path,
) -> None:
    """
    Require canonical headers only where direct execution is intentional.
    """

    direct = tmp_path / "tests" / "test_direct.py"
    library = tmp_path / "dev" / "library.py"

    direct.parent.mkdir(parents=True)
    library.parent.mkdir(parents=True)

    direct.write_text(
        """
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_direct.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI Codex (GPT-5.6) was used in development.
#
# Distributed under the MIT license.


if __name__ == "__main__":
    raise SystemExit(0)
""".removeprefix("\n"),
        encoding="utf-8",
    )
    library.write_text(
        '"""Provide an import-only fixture module."""\n',
        encoding="utf-8",
    )

    paths = ["dev/library.py", "tests/test_direct.py"]
    config = activate_review(tmp_path, paths, configuration())
    headers = produce(
        tmp_path,
        paths,
        config,
    )["inventories"]["source_headers"]
    by_path = {record["path"]: record for record in headers["records"]}

    assert by_path["tests/test_direct.py"]["structure_valid"] is True
    assert by_path["tests/test_direct.py"]["disposition"] in {
        "refactored",
        "retained_coherent",
    }
    assert by_path["dev/library.py"]["disposition"] == "omitted_by_role"
    assert headers["disposition_counts"]["unresolved"] == 0


def test_source_header_review_can_adopt_audit_script_sources(
    tmp_path: Path,
) -> None:
    """
    Keep an explicit script-source decision separate from guard discovery.
    """

    script = tmp_path / "dev" / "audit" / "generate_prompts.py"
    script.parent.mkdir(parents=True)
    script.write_text(
        """
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: generate_prompts.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI Codex (GPT-5.6) was used in development.
#
# Distributed under the MIT license.


\"\"\"Provide an importable audit-script implementation.\"\"\"
""".removeprefix("\n"),
        encoding="utf-8",
    )
    paths = ["dev/audit/generate_prompts.py"]
    config = activate_review(tmp_path, paths, configuration())

    record = produce(
        tmp_path,
        paths,
        config,
    )["inventories"]["source_headers"]["records"][0]

    assert record["adopted_script_source"] is True
    assert record["applicable"] is True
    assert record["structure_valid"] is True
    assert record["disposition"] == "refactored"
