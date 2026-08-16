#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_dependency_closure.py
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
Focused regressions for configured dependency-closure packages.
"""

from __future__ import annotations

import copy
import json
import subprocess
import tempfile
import unittest
from collections.abc import Callable
from pathlib import Path
from types import SimpleNamespace

from dev.audit.run import load_target_selection


class DependencyClosurePackageTest(unittest.TestCase):
    """
    Exercise family-neutral package metadata and evidence selection.
    """

    def make_tracked_shell(
        self,
        source: str,
    ) -> tuple[tempfile.TemporaryDirectory, Path]:
        temporary = tempfile.TemporaryDirectory()
        root = Path(temporary.name)

        subprocess.run(["git", "init", "-q", str(root)], check=True)
        subprocess.run(
            [
                "git",
                "-C",
                str(root),
                "config",
                "user.email",
                "audit@example.test",
            ],
            check=True,
        )
        subprocess.run(
            ["git", "-C", str(root), "config", "user.name", "Audit Test"],
            check=True,
        )

        (root / "sample.sh").write_text(source, encoding="utf-8")
        subprocess.run(
            ["git", "-C", str(root), "add", "sample.sh"],
            check=True,
        )
        subprocess.run(
            ["git", "-C", str(root), "commit", "-qm", "fixture"],
            check=True,
        )

        return temporary, root

    def synthetic_selection(self) -> dict[str, object]:
        return {
            "targets": [{"path": "sample.sh", "role": "primary"}],
            "context": [],
            "evidence_selection": {
                "target_selector": {
                    "kind": "changed_shell_semantic_units",
                    "selection_reason": "synthetic changed unit",
                },
                "context_selectors": [],
            },
        }

    def minimal_bundle_selection(self) -> dict[str, object]:
        return {
            "package": {
                "bundle_id": "synthetic-package",
                "title": "Synthetic package",
                "semantic_review_path": "semantic_review/synthetic.md",
                "supplied_artifacts": [
                    "semantic_review/synthetic.md",
                    "findings.ndjson",
                    "facts.ndjson",
                    "adapter_limitations.ndjson",
                    "pilot_report.json",
                ],
                "statements": [],
                "semantic_only_topics": [],
                "size_limits": {
                    "semantic_markdown_max_bytes": 1048576,
                    "semantic_markdown_max_lines": 10000,
                },
            },
            "targets": [
                {
                    "path": "sample.sh",
                    "role": "primary",
                    "target_role": "synthetic target",
                    "subcohort": "shared_runtime",
                },
            ],
            "rule_scopes": [],
            "dependency_edges": [],
            "semantic_questions": [],
        }

    def valid_unit_coverage(self) -> dict[str, object]:
        return {
            "configured_semantic_unit_count": 1,
            "function_unit_count": 1,
            "top_level_region_count": 0,
            "changed_block_count": 1,
            "changed_blocks": [
                {"represented": True, "exact_diff_retained": True},
            ],
            "uncovered_changed_blocks": [],
            "overlapping_units": [],
            "changed_top_level_regions": [],
            "segmentation_failures": [],
            "all_changed_blocks_covered": True,
            "evidence_truncated": False,
        }

    def load_child(self, child: str) -> dict[str, object]:
        """
        Load one real linked child with its dry run-group identity.
        """

        root = Path(__file__).resolve().parents[3]
        selection = load_target_selection(
            SimpleNamespace(
                path_values=None,
                paths_from=Path(
                    (
                        "dev/config/pilots/download_fastqs_dependency_closure."
                        "json"
                    ),
                ),
                package_child=child,
                umbrella_run_id="linked-test-run-group",
            ),
            root,
        )

        self.assertIsNotNone(selection)

        return selection

    def load_modified_linked_config(
        self,
        mutate: Callable[[dict[str, object]], None],
    ) -> dict[str, object]:
        """
        Load a temporary external-scratch mutation of the linked config.
        """

        root = Path(__file__).resolve().parents[3]
        source = (
            root / "dev/config/pilots/download_fastqs_dependency_closure.json"
        )
        value = json.loads(source.read_text(encoding="utf-8"))
        mutate(value)

        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            suffix=".json",
            prefix="linked-config-",
        ) as handle:
            json.dump(value, handle)
            handle.flush()

            return load_target_selection(
                SimpleNamespace(
                    path_values=None,
                    paths_from=Path(handle.name),
                    package_child=(
                        "download-fastqs-dependency-closure-runtime-production"
                    ),
                    umbrella_run_id="linked-test-run-group",
                ),
                root,
            )

    def test_nested_standalone_group_closure_does_not_end_function(
        self,
    ) -> None:
        from dev.audit.generate_prompts import _shell_source_units

        source = """
function demo() {
    require_optarg x y || {
        return 1
}
    echo after_group
}
""".removeprefix("\n")

        units = _shell_source_units("sample.sh", source)
        function = next(
            unit for unit in units if unit["unit_kind"] == "shell_function"
        )

        self.assertIn("echo after_group", function["source"])
        self.assertEqual(function["end_line"], 6)

    def test_multiple_nested_groups_remain_in_one_complete_function(
        self,
    ) -> None:
        from dev.audit.generate_prompts import _shell_source_units

        source = """
function demo() {
    first || {
        return 1
}
    second || {
        return 1
}
    echo complete
}
"""
        functions = [
            unit
            for unit in _shell_source_units("sample.sh", source)
            if unit["unit_kind"] == "shell_function"
        ]

        self.assertEqual(len(functions), 1)
        self.assertEqual(functions[0]["source"].count("return 1"), 2)
        self.assertIn("echo complete", functions[0]["source"])

    def test_last_function_excludes_following_top_level_dispatch(self) -> None:
        from dev.audit.generate_prompts import _shell_source_units

        source = """
function main() {
    echo body
}

main "$@"
""".removeprefix("\n")

        units = _shell_source_units("sample.sh", source)
        function = next(
            unit for unit in units if unit["unit_kind"] == "shell_function"
        )
        top_level = next(
            unit for unit in units if unit["unit_kind"] == "top_level_region"
        )

        self.assertNotIn('main "$@"', function["source"])
        self.assertIn('main "$@"', top_level["source"])

    def test_top_level_bootstrap_group_is_one_complete_top_level_region(
        self,
    ) -> None:
        from dev.audit.generate_prompts import _shell_source_units

        source = """
source ./helpers.sh || {
    echo failed
    return 1
}

function main() {
    :
}
"""

        units = _shell_source_units("sample.sh", source)
        bootstrap = units[0]

        self.assertEqual(bootstrap["unit_kind"], "top_level_region")
        self.assertIn("source ./helpers.sh || {", bootstrap["source"])
        self.assertIn("return 1\n}", bootstrap["source"])

    def test_deletion_only_hunk_inside_function_keeps_before_unit(
        self,
    ) -> None:
        from dev.audit.generate_prompts import configured_semantic_units

        before = """
function demo() {
    echo retained
    echo removed_exactly
    echo tail
}
"""

        temporary, root = self.make_tracked_shell(before)
        self.addCleanup(temporary.cleanup)

        (root / "sample.sh").write_text(
            before.replace("    echo removed_exactly\n", ""),
            encoding="utf-8",
        )

        units, coverage = configured_semantic_units(
            root,
            self.synthetic_selection(),
        )
        demo = next(unit for unit in units if unit["unit_name"] == "demo")

        self.assertNotIn("before_source", demo)
        self.assertTrue(
            any(
                "-    echo removed_exactly" in item
                for item in demo["diff_evidence"]
            ),
        )
        self.assertTrue(coverage["all_changed_blocks_covered"])

    def test_deletion_only_top_level_block_keeps_complete_before_region(
        self,
    ) -> None:
        from dev.audit.generate_prompts import configured_semantic_units

        before = """
source ./helpers.sh || {
    echo bootstrap_removed
    return 1
}

function main() {
    :
}
"""

        temporary, root = self.make_tracked_shell(before)
        self.addCleanup(temporary.cleanup)

        (root / "sample.sh").write_text(
            "function main() {\n    :\n}\n",
            encoding="utf-8",
        )

        units, coverage = configured_semantic_units(
            root,
            self.synthetic_selection(),
        )
        before_regions = [
            unit
            for unit in units
            if unit["unit_kind"] == "top_level_region"
            and unit.get("source_state") == "before"
        ]

        self.assertEqual(len(before_regions), 1)
        self.assertIn("bootstrap_removed", before_regions[0]["source"])
        self.assertTrue(coverage["all_changed_blocks_covered"])

    def test_no_residual_function_body_is_emitted_as_top_level(self) -> None:
        from dev.audit.generate_prompts import _shell_source_units

        source = """
function demo() {
    first || {
        return 1
}
    echo must_stay_in_function
}

main "$@"
"""

        units = _shell_source_units("sample.sh", source)
        top_level_source = "".join(
            unit["source"]
            for unit in units
            if unit["unit_kind"] == "top_level_region"
        )

        self.assertNotIn("must_stay_in_function", top_level_source)
        self.assertIn('main "$@"', top_level_source)

    def test_changed_blocks_have_complete_units_and_exact_diffs(
        self,
    ) -> None:
        from dev.audit.generate_prompts import configured_semantic_units

        before = """
echo old_bootstrap

function demo() {
    echo old_body
    echo removed_body
}
"""
        temporary, root = self.make_tracked_shell(before)
        self.addCleanup(temporary.cleanup)
        (root / "sample.sh").write_text(
            """
echo new_bootstrap

function demo() {
    echo new_body
}
""",
            encoding="utf-8",
        )

        _, coverage = configured_semantic_units(
            root,
            self.synthetic_selection(),
        )

        self.assertTrue(coverage["all_changed_blocks_covered"])
        self.assertEqual(coverage["uncovered_changed_blocks"], [])
        self.assertTrue(coverage["changed_blocks"])

        for block in coverage["changed_blocks"]:
            self.assertTrue(block["exact_diff_retained"])
            self.assertTrue(
                block["complete_current_units"]
                or block["complete_before_units"],
            )

    def test_deletion_hunk_coordinates_are_zero_count_and_not_an_invented_line(
        self,
    ) -> None:
        from dev.audit.generate_prompts import _changed_blocks

        before = """
function demo() {
    echo keep
    echo delete_me
}
""".removeprefix("\n")

        temporary, root = self.make_tracked_shell(before)
        self.addCleanup(temporary.cleanup)

        (root / "sample.sh").write_text(
            before.replace("    echo delete_me\n", ""),
            encoding="utf-8",
        )

        blocks = _changed_blocks(root, "sample.sh", 3)
        deletion = next(block for block in blocks if block["new_count"] == 0)

        self.assertEqual(deletion["old_start"], 3)
        self.assertEqual(deletion["new_start"], 2)
        self.assertNotIn("start_line", deletion)
        self.assertIn("-    echo delete_me", deletion["diff_evidence"])

    def test_settled_relocation_reuses_historical_origin(self) -> None:
        from dev.audit.generate_prompts import (
            _changed_blocks,
            _relocation_origin,
        )

        temporary, root = self.make_tracked_shell("echo current\n")
        self.addCleanup(temporary.cleanup)

        legacy = root / "legacy/sample.sh"
        legacy.parent.mkdir()
        legacy.write_text("echo legacy\n", encoding="utf-8")

        subprocess.run(
            ["git", "-C", str(root), "add", "legacy/sample.sh"],
            check=True,
        )
        subprocess.run(
            ["git", "-C", str(root), "commit", "-qm", "add legacy"],
            check=True,
        )

        subprocess.run(
            ["git", "-C", str(root), "rm", "-q", "legacy/sample.sh"],
            check=True,
        )
        subprocess.run(
            ["git", "-C", str(root), "commit", "-qm", "settle relocation"],
            check=True,
        )

        config = root / "dev/config/path_relocations.json"
        config.parent.mkdir(parents=True)
        config.write_text(
            json.dumps(
                {
                    "exact": {
                        "sample.sh": "legacy/sample.sh",
                    },
                },
            ),
            encoding="utf-8",
        )

        origin = _relocation_origin(root, "sample.sh")
        blocks = _changed_blocks(root, "sample.sh", 1, origin)

        self.assertEqual(origin, "legacy/sample.sh")
        self.assertIn("-echo legacy", blocks[0]["diff_evidence"])
        self.assertIn("+echo current", blocks[0]["diff_evidence"])

    def test_representative_real_nested_group_function_is_complete(
        self,
    ) -> None:
        from dev.audit.generate_prompts import _shell_source_units

        root = Path(__file__).resolve().parents[3]
        path = "bin/execute_download_fastqs.sh"
        source = (root / path).read_text(encoding="utf-8")

        units = _shell_source_units(path, source)
        parse_args = next(
            unit for unit in units if unit["unit_name"] == "parse_args"
        )

        self.assertIn(
            'require_optarg "${1}" "${2:-}" "main" || {',
            parse_args["source"],
        )
        self.assertIn('fil_in="${2}"', parse_args["source"])
        self.assertTrue(parse_args["source"].rstrip().endswith("}"))

    def test_selection_loads_package_roles_scopes_questions_edges_and_units(
        self,
    ) -> None:
        root = Path(__file__).resolve().parents[3]

        selection = load_target_selection(
            SimpleNamespace(
                path_values=None,
                paths_from=Path(
                    (
                        "dev/config/pilots/download_fastqs_dependency_closure."
                        "json"
                    ),
                ),
            ),
            root,
        )

        self.assertEqual(selection["package"]["bundle_id"], selection["name"])
        self.assertEqual(len(selection["targets"]), 27)
        self.assertEqual(
            {target["subcohort"] for target in selection["targets"]},
            {
                "shared_runtime",
                "download_production_help",
                "shared_test_infrastructure",
                "download_tests_registration",
            },
        )

        self.assertEqual(
            selection["evidence_selection"]["target_selector"]["kind"],
            "changed_shell_semantic_units",
        )
        self.assertEqual(
            len(selection["evidence_selection"]["context_selectors"]),
            5,
        )

        self.assertEqual(
            {scope["rule_id"] for scope in selection["rule_scopes"]},
            {"SHELL.SYNTAX", "TEST.PROOF.PROPORTIONAL", "SHELL.LINE_LENGTH"},
        )

        self.assertTrue(selection["semantic_questions"])
        self.assertTrue(selection["dependency_edges"])

    def test_configured_units_cover_changes_and_context_functions(
        self,
    ) -> None:
        from dev.audit.generate_prompts import configured_semantic_units

        root = Path(__file__).resolve().parents[3]

        selection = load_target_selection(
            SimpleNamespace(
                path_values=None,
                paths_from=Path(
                    (
                        "dev/config/pilots/download_fastqs_dependency_closure."
                        "json"
                    ),
                ),
            ),
            root,
        )
        units, coverage = configured_semantic_units(root, selection)

        self.assertTrue(coverage["all_changed_blocks_covered"])
        self.assertFalse(coverage["evidence_truncated"])

        wrap_units = [
            unit
            for unit in units
            if unit["path"] == "lib/bash/core/wrap_cmd.sh"
        ]

        self.assertIn(
            "function get_submit_logs()",
            "\n".join(unit["source"] for unit in wrap_units),
        )
        self.assertIn(
            "function print_built_cmd()",
            "\n".join(unit["source"] for unit in wrap_units),
        )

        for unit in wrap_units:
            if unit["unit_name"] in {"get_submit_logs", "print_built_cmd"}:
                self.assertTrue(unit["source"].rstrip().endswith("}"))

        context_units = [
            unit for unit in units if unit["evidence_role"] == "context"
        ]

        self.assertEqual(len(context_units), 5)
        self.assertTrue(
            all(unit["unit_name"] == "run_jobs" for unit in context_units),
        )
        self.assertTrue(
            all("get_submit_logs" in unit["source"] for unit in context_units),
        )

    def test_configured_renderer_uses_complete_dependency_context(
        self,
    ) -> None:
        from dev.audit.generate_prompts import (
            configured_semantic_units,
            write_configured_bundle,
        )

        root = Path(__file__).resolve().parents[3]

        selection = load_target_selection(
            SimpleNamespace(
                path_values=None,
                paths_from=Path(
                    (
                        "dev/config/pilots/download_fastqs_dependency_closure."
                        "json"
                    ),
                ),
            ),
            root,
        )
        units, coverage = configured_semantic_units(root, selection)
        targets = [
            {
                **target,
                "git_state_labels": ["tracked_modified"],
                "checks_run": [],
                "findings_count": 0,
                "content_fingerprint": "sha256:" + "1" * 64,
            }
            for target in selection["targets"]
        ]

        with tempfile.TemporaryDirectory() as temp_dir:
            report = Path(temp_dir)

            record = write_configured_bundle(
                report,
                selection,
                targets,
                [],
                [],
                selection["adapter_limitations"],
                units,
                coverage,
                "configured-test",
            )[0]
            markdown = (
                report / selection["package"]["semantic_review_path"]
            ).read_text(encoding="utf-8")

        self.assertEqual(
            record["bundle_id"],
            "download-fastqs-dependency-closure",
        )
        self.assertIn(selection["package"]["statements"][0], markdown)
        self.assertIn(
            "| Edge kind | Reachability | Runtime status |",
            markdown,
        )
        self.assertIn("## Declarative semantic questions", markdown)

        self.assertEqual(record["semantic_unit_coverage"], coverage)
        self.assertEqual(
            markdown.count("Selection reason:"),
            coverage["configured_semantic_unit_count"],
        )

        self.assertIn("Exact change evidence:", markdown)
        self.assertIn("add_shell_tests", markdown)
        self.assertNotIn("[truncated]", markdown)

    def test_existing_reports_option_accepts_bounded_private_tmp_dry_render(
        self,
    ) -> None:
        from dev.audit.run import resolve_reports_base

        root = Path(__file__).resolve().parents[3]
        dry_base = Path("/private/tmp/dependency-closure-dry-base")

        self.assertEqual(
            resolve_reports_base(root, Path("/private/tmp")),
            Path("/private/tmp"),
        )
        self.assertEqual(resolve_reports_base(root, dry_base), dry_base)
        self.assertEqual(
            resolve_reports_base(root, Path("artifacts/dev/audit")),
            root / "artifacts/dev/audit",
        )

        with self.assertRaisesRegex(
            ValueError,
            "artifacts/dev/audit or /private/tmp",
        ):
            resolve_reports_base(root, root.parent / "unexpected-reports")

    def test_configured_rule_scopes_select_only_declared_rules_and_paths(
        self,
    ) -> None:
        from dev.audit.run import configured_rule_selection

        rules = [
            {"rule_id": "ONE", "applicable_paths": ["**"]},
            {"rule_id": "TWO", "applicable_paths": ["**"]},
            {"rule_id": "UNDECLARED", "applicable_paths": ["**"]},
        ]
        selection = {
            "package": {},
            "rule_scopes": [
                {
                    "rule_id": "ONE",
                    "paths": ["scripts/**/*.sh"],
                    "classification": "strict",
                },
                {
                    "rule_id": "TWO",
                    "paths": ["tests/**/*.sh"],
                    "classification": "advisory",
                },
            ],
        }
        selected, scopes = configured_rule_selection(rules, selection, None)

        self.assertEqual(
            [rule["rule_id"] for rule in selected],
            ["ONE", "TWO"],
        )
        self.assertEqual(
            scopes,
            {"ONE": ["scripts/**/*.sh"], "TWO": ["tests/**/*.sh"]},
        )

    def test_bundle_rejects_every_incomplete_or_truncated_coverage_state(
        self,
    ) -> None:
        from dev.audit.generate_prompts import write_configured_bundle

        target = {
            **self.minimal_bundle_selection()["targets"][0],
            "git_state_labels": ["tracked_modified"],
            "checks_run": [],
            "findings_count": 0,
            "content_fingerprint": "sha256:" + "1" * 64,
        }
        unit = {
            "path": "sample.sh",
            "unit_name": "demo",
            "unit_kind": "shell_function",
            "source_state": "current",
            "evidence_role": "target",
            "start_line": 1,
            "end_line": 3,
            "selection_reason": "synthetic",
            "source": "function demo() {\n    :\n}\n",
            "diff_evidence": ["@@ -2 +2 @@\n-    old\n+    :\n"],
        }
        cases = {
            "coverage flag": {"all_changed_blocks_covered": False},
            "uncovered hunk": {
                "uncovered_changed_blocks": [{"block_id": "one"}],
            },
            "overlap": {"overlapping_units": [{"left": "a", "right": "b"}]},
            "segmentation": {"segmentation_failures": [{"path": "sample.sh"}]},
            "truncation": {"evidence_truncated": True},
        }

        for label, mutation in cases.items():
            coverage = {**self.valid_unit_coverage(), **mutation}

            with (
                self.subTest(label=label),
                tempfile.TemporaryDirectory() as temp_dir,
                self.assertRaisesRegex(
                    ValueError,
                    "incomplete semantic-unit coverage",
                ),
            ):
                write_configured_bundle(
                    Path(temp_dir),
                    self.minimal_bundle_selection(),
                    [target],
                    [],
                    [],
                    [],
                    [unit],
                    coverage,
                    "guard-test",
                )

    def test_dependency_graph_rejects_missing_or_conflicting_edges(
        self,
    ) -> None:
        from dev.audit.run import validate_dependency_graph

        root = Path(__file__).resolve().parents[3]
        raw = __import__("json").loads(
            (
                root
                / "dev/config/pilots/download_fastqs_dependency_closure.json"
            ).read_text(encoding="utf-8"),
        )
        known = {target["path"] for target in raw["targets"]} | {
            item["path"] if isinstance(item, dict) else item
            for item in raw["context"]
        }

        report = validate_dependency_graph(
            raw["dependency_edges"],
            known,
            raw["edge_requirements"],
        )

        self.assertEqual(
            set(report["satisfied_requirement_groups"]),
            {
                "wrapper_roots_to_shared_runtime",
                "download_tests_to_wrapper_roots",
                "download_tests_to_shared_test_infrastructure",
                "runner_to_six_download_tests",
                "fixture_and_runner_ownership",
                "get_submit_logs_current_callers",
            },
        )

        required = raw["edge_requirements"][0]
        missing_identity = (
            required["edges"][0]["from"],
            required["edges"][0]["to"],
        )
        missing = [
            edge
            for edge in raw["dependency_edges"]
            if (edge["from"], edge["to"]) != missing_identity
        ]

        duplicate = raw["dependency_edges"] + [
            copy.deepcopy(raw["dependency_edges"][0]),
        ]

        contradictory = copy.deepcopy(raw["dependency_edges"])
        changed = copy.deepcopy(contradictory[0])
        changed["runtime_status"] = "unknown"
        contradictory.append(changed)

        undeclared = copy.deepcopy(raw["dependency_edges"])
        undeclared[0]["to"] = "undeclared/path.sh"

        cases = {
            "missing": (missing, "missing required dependency edge"),
            "duplicate": (duplicate, "duplicate dependency edge"),
            "contradictory": (
                contradictory,
                "contradictory dependency edge",
            ),
            "undeclared": (undeclared, "undeclared dependency endpoint"),
        }

        for label, (edges, message) in cases.items():
            with (
                self.subTest(label=label),
                self.assertRaisesRegex(ValueError, message),
            ):
                validate_dependency_graph(
                    edges,
                    known,
                    raw["edge_requirements"],
                )

    def test_actual_rule_scopes_report_exact_cohort_matches_and_enforcement(
        self,
    ) -> None:
        from dev.audit.run import configured_rule_scope_report

        root = Path(__file__).resolve().parents[3]

        selection = load_target_selection(
            SimpleNamespace(
                path_values=None,
                paths_from=Path(
                    (
                        "dev/config/pilots/download_fastqs_dependency_closure."
                        "json"
                    ),
                ),
            ),
            root,
        )
        report = configured_rule_scope_report(selection)
        expected = [target["path"] for target in selection["targets"]]

        self.assertEqual(
            set(report),
            {"SHELL.SYNTAX", "TEST.PROOF.PROPORTIONAL", "SHELL.LINE_LENGTH"},
        )

        for rule_id, row in report.items():
            with self.subTest(rule_id=rule_id):
                self.assertEqual(row["matched_targets"], expected)
                self.assertNotIn("not/in/cohort.sh", row["matched_targets"])

        self.assertEqual(report["SHELL.SYNTAX"]["enforcement"], "strict")
        self.assertEqual(
            report["TEST.PROOF.PROPORTIONAL"]["enforcement"],
            "strict",
        )
        self.assertEqual(
            report["SHELL.LINE_LENGTH"]["enforcement"],
            "advisory",
        )

    def test_function_modification_omits_redundant_before_source(
        self,
    ) -> None:
        from dev.audit.generate_prompts import configured_semantic_units

        before = "function demo() {\n    echo old\n}\n"

        temporary, root = self.make_tracked_shell(before)
        self.addCleanup(temporary.cleanup)

        (root / "sample.sh").write_text(
            "function demo() {\n    echo new\n}\n",
            encoding="utf-8",
        )

        units, coverage = configured_semantic_units(
            root,
            self.synthetic_selection(),
        )
        demo = next(unit for unit in units if unit["unit_name"] == "demo")

        self.assertNotIn("before_source", demo)
        self.assertEqual(len(demo["diff_evidence"]), 1)
        self.assertTrue(coverage["all_changed_blocks_covered"])

    def test_whole_deleted_function_retains_complete_before_state_unit(
        self,
    ) -> None:
        from dev.audit.generate_prompts import configured_semantic_units

        before = """
function removed() {
    echo deleted_body
}

function retained() {
    :
}
"""

        temporary, root = self.make_tracked_shell(before)
        self.addCleanup(temporary.cleanup)

        (root / "sample.sh").write_text(
            "function retained() {\n    :\n}\n",
            encoding="utf-8",
        )

        units, coverage = configured_semantic_units(
            root,
            self.synthetic_selection(),
        )
        removed = next(
            unit for unit in units if unit["unit_name"] == "removed"
        )

        self.assertEqual(removed["source_state"], "before")
        self.assertIn("function removed()", removed["source"])
        self.assertIn("deleted_body", removed["source"])
        self.assertTrue(coverage["all_changed_blocks_covered"])

    def test_partial_function_deletion_keeps_exact_lines_only(
        self,
    ) -> None:
        from dev.audit.generate_prompts import configured_semantic_units

        before = """
function demo() {
    echo retained
    echo deleted_line
}
"""

        temporary, root = self.make_tracked_shell(before)
        self.addCleanup(temporary.cleanup)

        (root / "sample.sh").write_text(
            before.replace("    echo deleted_line\n", ""),
            encoding="utf-8",
        )

        units, coverage = configured_semantic_units(
            root,
            self.synthetic_selection(),
        )
        demo = next(unit for unit in units if unit["unit_name"] == "demo")

        self.assertNotIn("before_source", demo)
        self.assertIn("-    echo deleted_line", "".join(demo["diff_evidence"]))
        self.assertTrue(
            all(block["represented"] for block in coverage["changed_blocks"]),
        )

    def test_partial_top_level_deletion_uses_current_region_and_exact_diff(
        self,
    ) -> None:
        from dev.audit.generate_prompts import configured_semantic_units

        before = "echo retained\necho deleted_line\necho tail\n"

        temporary, root = self.make_tracked_shell(before)
        self.addCleanup(temporary.cleanup)

        (root / "sample.sh").write_text(
            before.replace("echo deleted_line\n", ""),
            encoding="utf-8",
        )

        units, coverage = configured_semantic_units(
            root,
            self.synthetic_selection(),
        )
        current = [unit for unit in units if unit["source_state"] == "current"]
        before_units = [
            unit for unit in units if unit["source_state"] == "before"
        ]

        self.assertEqual(len(current), 1)
        self.assertEqual(before_units, [])
        self.assertIn(
            "-echo deleted_line",
            "".join(current[0]["diff_evidence"]),
        )
        self.assertTrue(coverage["all_changed_blocks_covered"])

    def test_coverage_metadata_fingerprints_without_duplicating_exact_diff(
        self,
    ) -> None:
        from dev.audit.generate_prompts import configured_semantic_units

        before = "function demo() {\n    echo old\n}\n"

        temporary, root = self.make_tracked_shell(before)
        self.addCleanup(temporary.cleanup)

        (root / "sample.sh").write_text(
            "function demo() {\n    echo new\n}\n",
            encoding="utf-8",
        )

        units, coverage = configured_semantic_units(
            root,
            self.synthetic_selection(),
        )
        demo = next(unit for unit in units if unit["unit_name"] == "demo")
        block = coverage["changed_blocks"][0]

        self.assertIn("diff_evidence", demo)
        self.assertNotIn("diff_evidence", block)
        self.assertRegex(
            block["diff_evidence_fingerprint"],
            r"^sha256:[0-9a-f]{64}$",
        )
        self.assertTrue(block["exact_diff_retained"])

    def test_markdown_contributions_cover_every_line_and_byte(self) -> None:
        from dev.audit.generate_prompts import (
            configured_semantic_units,
            write_configured_bundle,
        )

        root = Path(__file__).resolve().parents[3]

        selection = load_target_selection(
            SimpleNamespace(
                path_values=None,
                paths_from=Path(
                    (
                        "dev/config/pilots/download_fastqs_dependency_closure."
                        "json"
                    ),
                ),
            ),
            root,
        )
        units, coverage = configured_semantic_units(root, selection)
        targets = [
            {
                **target,
                "git_state_labels": ["tracked_modified"],
                "checks_run": [],
                "findings_count": 0,
                "content_fingerprint": "sha256:" + "2" * 64,
            }
            for target in selection["targets"]
        ]

        with tempfile.TemporaryDirectory() as temp_dir:
            record = write_configured_bundle(
                Path(temp_dir),
                selection,
                targets,
                [],
                [],
                selection["adapter_limitations"],
                units,
                coverage,
                "contribution-test",
            )[0]

        contributions = record["markdown_contributions"]

        self.assertEqual(
            set(contributions),
            {
                "shared_runtime",
                "download_production_help",
                "shared_test_infrastructure",
                "download_tests_registration",
                "context_only_callers",
                "package_overhead",
            },
        )
        self.assertEqual(
            sum(row["lines"] for row in contributions.values()),
            record["markdown_lines"],
        )
        self.assertEqual(
            sum(row["bytes"] for row in contributions.values()),
            record["markdown_bytes"],
        )

    def test_linked_children_filter_exact_ordered_partitions_and_rule_scopes(
        self,
    ) -> None:
        runtime_id = "download-fastqs-dependency-closure-runtime-production"
        tests_id = "download-fastqs-dependency-closure-tests-registration"
        runtime = self.load_child(runtime_id)
        tests = self.load_child(tests_id)

        self.assertEqual(
            (
                len(runtime["targets"]),
                len(runtime["context"]),
                len(runtime["dependency_edges"]),
            ),
            (13, 8, 28),
        )
        self.assertEqual(
            (
                len(tests["targets"]),
                len(tests["context"]),
                len(tests["dependency_edges"]),
            ),
            (14, 4, 24),
        )
        self.assertEqual(
            [item["path"] for item in runtime["targets"] + tests["targets"]],
            [
                item["path"]
                for item in runtime["linked_package"][
                    "umbrella_target_ownership"
                ]
            ],
        )

        combined_edges = {
            (edge["from"], edge["to"]): edge
            for edge in runtime["dependency_edges"] + tests["dependency_edges"]
        }

        self.assertEqual(
            [
                combined_edges[(edge["from"], edge["to"])]
                for edge in runtime["linked_package"][
                    "umbrella_dependency_edges"
                ]
            ],
            runtime["linked_package"]["umbrella_dependency_edges"],
        )
        self.assertEqual(len(runtime["semantic_questions"]), 4)
        self.assertEqual(len(tests["semantic_questions"]), 2)
        self.assertEqual(len(runtime["adapter_limitations"]), 4)
        self.assertEqual(len(tests["adapter_limitations"]), 4)

        for selection, patterns in (
            (runtime, ["bin/**/*.sh", "lib/bash/**/*.sh"]),
            (tests, ["tests/**/*.sh"]),
        ):
            expected = [item["path"] for item in selection["targets"]]

            self.assertEqual(
                selection["rule_scope_report"]["SHELL.SYNTAX"]["paths"],
                patterns,
            )
            self.assertEqual(
                selection["rule_scope_report"]["SHELL.SYNTAX"][
                    "matched_targets"
                ],
                expected,
            )
            self.assertEqual(
                selection["rule_scope_report"]["TEST.PROOF.PROPORTIONAL"][
                    "matched_targets"
                ],
                expected,
            )
            self.assertEqual(
                selection["rule_scope_report"]["SHELL.LINE_LENGTH"][
                    "matched_targets"
                ],
                expected,
            )

    def test_linked_configuration_rejects_child_ceiling_and_ownership_defects(
        self,
    ) -> None:
        def third_child(value: dict[str, object]) -> None:
            value["package"]["children"].append(
                copy.deepcopy(value["package"]["children"][0]),
            )

        with self.assertRaisesRegex(ValueError, "at most two"):
            self.load_modified_linked_config(third_child)

        mutations = {
            "target": lambda value: value["targets"][0].update(
                subcohort="unowned",
            ),
            "context": lambda value: value["context"][0].update(
                child="unknown-child",
            ),
            "edge": lambda value: value["dependency_edges"][0].update(
                child="unknown-child",
            ),
        }

        for label, mutate in mutations.items():
            with (
                self.subTest(label=label),
                self.assertRaisesRegex(ValueError, "owner"),
            ):
                self.load_modified_linked_config(mutate)

        def duplicate_target_owner(value: dict[str, object]) -> None:
            value["package"]["children"][1]["subcohorts"].append(
                "shared_runtime",
            )

        with self.assertRaisesRegex(ValueError, "exactly one"):
            self.load_modified_linked_config(duplicate_target_owner)

    def test_linked_graph_validates_umbrella_endpoints_before_child_filtering(
        self,
    ) -> None:
        tests = self.load_child(
            "download-fastqs-dependency-closure-tests-registration",
        )
        cross_edges = [
            edge
            for edge in tests["dependency_edges"]
            if edge["from"].startswith("tests/")
            and edge["to"].startswith(("bin/", "lib/bash/"))
        ]

        self.assertEqual(len(cross_edges), 7)

        def unknown_endpoint(value: dict[str, object]) -> None:
            value["dependency_edges"][22]["to"] = "bin/not-declared.sh"

        with self.assertRaisesRegex(ValueError, "endpoints must be declared"):
            self.load_modified_linked_config(unknown_endpoint)

    def test_runtime_linked_semantic_units_preserve_final_counts(self) -> None:
        from dev.audit.generate_prompts import configured_semantic_units

        root = Path(__file__).resolve().parents[3]
        runtime = self.load_child(
            "download-fastqs-dependency-closure-runtime-production",
        )

        _, runtime_coverage = configured_semantic_units(root, runtime)

        # Structural coverage holds for the whole cohort in any tree state, as
        # it follows from file content rather than commit status.
        self.assertTrue(runtime_coverage["all_changed_blocks_covered"])
        self.assertEqual(runtime_coverage["uncovered_changed_blocks"], [])
        self.assertEqual(runtime_coverage["overlapping_units"], [])
        self.assertEqual(runtime_coverage["segmentation_failures"], [])
        self.assertFalse(runtime_coverage["evidence_truncated"])

        # Select the targets carrying a frozen historical baseline. Every other
        # target falls back to 'git diff HEAD', so its block count depends on
        # whether the edit is committed and cannot be pinned.
        relocations = json.loads(
            (root / "dev/config/path_relocations.json").read_text(
                encoding="utf-8",
            ),
        )["exact"]
        historical = copy.deepcopy(runtime)
        historical["targets"] = [
            target
            for target in runtime["targets"]
            if str(target["path"]) in relocations
        ]

        # Pin the partition itself, so adding a non-relocated target fails here
        # rather than silently dropping out of the counts below.
        self.assertEqual(len(runtime["targets"]), 13)
        self.assertEqual(len(historical["targets"]), 12)

        historical_units, historical_coverage = configured_semantic_units(
            root,
            historical,
        )

        # Taken over the historical-baseline subset, these counts are identical
        # in a clean checkout and in a dirty working tree.
        self.assertEqual(
            (
                len(historical_units),
                historical_coverage["changed_block_count"],
            ),
            (91, 313),
        )

    def test_tests_linked_semantic_units_report_derived_counts(self) -> None:
        from dev.audit.generate_prompts import configured_semantic_units

        root = Path(__file__).resolve().parents[3]
        tests = self.load_child(
            "download-fastqs-dependency-closure-tests-registration",
        )

        test_units, test_coverage = configured_semantic_units(root, tests)

        self.assertEqual(
            len(test_units),
            test_coverage["configured_semantic_unit_count"],
        )
        self.assertEqual(
            test_coverage["changed_block_count"],
            len(test_coverage["changed_blocks"]),
        )

        self.assertGreater(len(test_units), 0)
        self.assertGreater(test_coverage["changed_block_count"], 0)

        self.assertTrue(test_coverage["all_changed_blocks_covered"])
        self.assertEqual(test_coverage["uncovered_changed_blocks"], [])
        self.assertEqual(test_coverage["overlapping_units"], [])
        self.assertEqual(test_coverage["segmentation_failures"], [])
        self.assertFalse(test_coverage["evidence_truncated"])

    def test_linked_children_share_bounded_final_size_limits(self) -> None:
        for child in (
            "download-fastqs-dependency-closure-runtime-production",
            "download-fastqs-dependency-closure-tests-registration",
        ):
            with self.subTest(child=child):
                limits = self.load_child(child)["package"]["size_limits"]

                self.assertEqual(limits["semantic_markdown_max_lines"], 16000)
                self.assertEqual(
                    limits["semantic_markdown_max_bytes"],
                    1048576,
                )

    def test_linked_renderer_names_sibling_and_shared_graph(self) -> None:
        from dev.audit.generate_prompts import (
            configured_semantic_units,
            write_configured_bundle,
        )

        root = Path(__file__).resolve().parents[3]

        selection = self.load_child(
            "download-fastqs-dependency-closure-runtime-production",
        )
        units, coverage = configured_semantic_units(root, selection)
        targets = [
            {
                **target,
                "git_state_labels": ["tracked_modified"],
                "checks_run": [],
                "findings_count": 0,
                "content_fingerprint": "sha256:" + "3" * 64,
            }
            for target in selection["targets"]
        ]

        with tempfile.TemporaryDirectory() as temporary:
            write_configured_bundle(
                Path(temporary),
                selection,
                targets,
                [],
                [],
                selection["adapter_limitations"],
                units,
                coverage,
                "linked-render-test",
            )

            markdown = (
                Path(temporary) / selection["package"]["semantic_review_path"]
            ).read_text(encoding="utf-8")

        self.assertIn("not the entire closure", markdown)
        self.assertIn(
            selection["linked_package"]["sibling_package_id"],
            markdown,
        )
        self.assertIn(
            selection["linked_package"]["graph_ownership_fingerprint"],
            markdown,
        )

        production_edges = markdown.count("→ `bin/") + markdown.count(
            "→ `lib/bash/",
        )

        self.assertEqual(production_edges, 7)

    def test_linked_header_validation_rejects_identity_mismatches(
        self,
    ) -> None:
        from dev.audit.run import validate_linked_pair_headers

        runtime = self.load_child(
            "download-fastqs-dependency-closure-runtime-production",
        )["linked_package"]
        tests = self.load_child(
            "download-fastqs-dependency-closure-tests-registration",
        )["linked_package"]

        pilots = [
            {
                "package_id": runtime["child_package_id"],
                "sibling_package_id": runtime["sibling_package_id"],
                "bundle_id": runtime["bundle_id"],
                "umbrella_run_id": "linked-test-run-group",
                "linked_package": runtime,
            },
            {
                "package_id": tests["child_package_id"],
                "sibling_package_id": tests["sibling_package_id"],
                "bundle_id": tests["bundle_id"],
                "umbrella_run_id": "linked-test-run-group",
                "linked_package": tests,
            },
        ]

        validate_linked_pair_headers(pilots)

        cases = {
            "bundle": (0, "bundle_id", "wrong-bundle"),
            "run group": (1, "umbrella_run_id", "wrong-run-group"),
            "child": (1, "package_id", pilots[0]["package_id"]),
            "config fingerprint": (1, "config_fingerprint", "sha256:wrong"),
            "graph fingerprint": (
                1,
                "graph_ownership_fingerprint",
                "sha256:wrong",
            ),
        }

        for label, (index, field, replacement) in cases.items():
            malformed = copy.deepcopy(pilots)

            if "fingerprint" in field:
                malformed[index]["linked_package"][field] = replacement
            else:
                malformed[index][field] = replacement

            with self.subTest(label=label), self.assertRaises(ValueError):
                validate_linked_pair_headers(malformed)

    def test_linked_partition_union_rejects_gaps_overlaps_and_contradictions(
        self,
    ) -> None:
        from dev.audit.run import validate_linked_partition_union

        runtime = self.load_child(
            "download-fastqs-dependency-closure-runtime-production",
        )
        tests = self.load_child(
            "download-fastqs-dependency-closure-tests-registration",
        )

        linkage = runtime["linked_package"]
        targets = [
            item["path"] for item in runtime["targets"] + tests["targets"]
        ]
        contexts = list(runtime["context"] + tests["context"])
        edges = list(runtime["dependency_edges"] + tests["dependency_edges"])

        validate_linked_partition_union(
            linkage,
            targets,
            contexts,
            edges,
            ["hunk-a", "hunk-b"],
        )

        cases = {
            "target gap": (targets[:-1], contexts, edges, ["hunk-a"]),
            "target overlap": (
                [*targets[:-1], targets[0]],
                contexts,
                edges,
                ["hunk-a"],
            ),
            "context gap": (targets, contexts[:-1], edges, ["hunk-a"]),
            "edge gap": (targets, contexts, edges[:-1], ["hunk-a"]),
            "edge overlap": (
                targets,
                contexts,
                [*edges, copy.deepcopy(edges[0])],
                ["hunk-a"],
            ),
            "edge contradiction": (
                targets,
                contexts,
                [{**edges[0], "runtime_status": "exercised"}, *edges[1:]],
                ["hunk-a"],
            ),
            "hunk overlap": (targets, contexts, edges, ["hunk-a", "hunk-a"]),
        }

        for label, arguments in cases.items():
            with self.subTest(label=label), self.assertRaises(ValueError):
                validate_linked_partition_union(linkage, *arguments)


if __name__ == "__main__":
    unittest.main()
