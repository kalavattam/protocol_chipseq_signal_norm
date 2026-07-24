#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_help_examples.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


"""Focused regressions for strict shell-help Examples documents."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

AUDIT_DIR = Path(__file__).resolve().parents[3] / "dev" / "audit"
REPO_ROOT = AUDIT_DIR.parents[1]
sys.path.insert(0, str(AUDIT_DIR))

from help_examples import (
    RepositoryResult,
    accepted_function_aliases,
    analyze_help_document,
    classify_wrapper_source,
    compliance_summary,
    invocation_facts,
    repository_crosswalk,
    scan_repository,
    short_help_advertises_details,
    undispositioned_review_findings,
)

USAGE = """Usage
-----
  demo
    [--help] [--mode <mode>] [--dry_run] value

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  --mode : {'one', 'two'}
    Processing mode.

  -dr, --dry_run : flag
    Run in dry-run mode.
"""

def examples(*entries: str, final: str = "") -> str:
    """Build one rendered document with a canonical Examples section."""

    body = f"{USAGE}\nExamples\n--------\n" + "\n\n".join(entries)
    if final:
        body += f"\n\n{final}"
    return body


def entry(number: int, description: str, *code: str) -> str:
    """Build one canonical numbered entry with one Bash code block."""

    rendered = "\n".join(f"    {line}" for line in code)
    return (
        f"  {number}. {description}\n"
        "    '''bash\n"
        f"{rendered}\n"
        "    '''"
    )


def analyze(body: str):
    """Analyze one strict in-memory owner with settled alias visibility."""

    return analyze_help_document(
        body,
        owner="demo",
        accepted_aliases={"-h", "--hlp", "--help", "--mode", "-dr", "--dry_run"},
        public_aliases={"-h", "--help", "--mode", "-dr", "--dry_run"},
        hidden_aliases={"--hlp"},
    )


def rules(body: str) -> set[str]:
    """Return strict rule identifiers for one fixture."""

    return {finding.rule_id for finding in analyze(body).findings}


class HelpExamplesTest(unittest.TestCase):
    """Prove strict structure, parser validity, and distinctness."""

    def test_missing_examples_fails_for_strict_owner(self) -> None:
        self.assertIn("HELP.EXAMPLES.REQUIRED", rules(USAGE))

    def test_one_example_fails(self) -> None:
        body = examples(entry(1, "Use mode one.", "demo --mode one value"))
        self.assertIn("HELP.EXAMPLES.COUNT", rules(body))

    def test_two_numbered_examples_pass(self) -> None:
        body = examples(
            entry(1, "Use mode one.", "demo --mode one value"),
            entry(2, "Use mode two in dry-run mode.", "demo --mode two --dry_run value"),
        )
        self.assertEqual(analyze(body).findings, [])

    def test_unnumbered_composite_examples_fail(self) -> None:
        body = f"{USAGE}\nExamples\n--------\n  '''bash\n  demo --mode one value\n  demo --mode two value\n  '''"
        self.assertIn("HELP.EXAMPLES.ENTRY", rules(body))

    def test_numbering_must_begin_at_one(self) -> None:
        body = examples(
            entry(2, "Use mode one.", "demo --mode one value"),
            entry(3, "Use mode two.", "demo --mode two value"),
        )
        self.assertIn("HELP.EXAMPLES.NUMBERING", rules(body))

    def test_skipped_or_duplicate_numbering_fails(self) -> None:
        skipped = examples(
            entry(1, "Use mode one.", "demo --mode one value"),
            entry(3, "Use mode two.", "demo --mode two value"),
        )
        duplicated = examples(
            entry(1, "Use mode one.", "demo --mode one value"),
            entry(1, "Use mode two.", "demo --mode two value"),
        )
        self.assertIn("HELP.EXAMPLES.NUMBERING", rules(skipped))
        self.assertIn("HELP.EXAMPLES.NUMBERING", rules(duplicated))

    def test_examples_must_be_final(self) -> None:
        body = examples(
            entry(1, "Use mode one.", "demo --mode one value"),
            entry(2, "Use mode two.", "demo --mode two value"),
            final="Notes\n-----\n  Late notes.",
        )
        self.assertIn("HELP.EXAMPLES.FINAL", rules(body))

    def test_more_than_one_examples_section_fails(self) -> None:
        first = examples(
            entry(1, "Use mode one.", "demo --mode one value"),
            entry(2, "Use mode two.", "demo --mode two value"),
        )
        body = f"{first}\n\nExamples\n--------\n{entry(1, 'Repeat.', 'demo --mode one value')}"
        self.assertIn("HELP.EXAMPLES.SECTION_COUNT", rules(body))

    def test_entry_without_code_block_fails(self) -> None:
        body = examples(
            "  1. Use mode one.",
            entry(2, "Use mode two.", "demo --mode two value"),
        )
        self.assertIn("HELP.EXAMPLES.CODE_BLOCK", rules(body))

    def test_entry_with_multiple_code_blocks_fails(self) -> None:
        double = (
            entry(1, "Use mode one.", "demo --mode one value")
            + "\n    '''bash\n    demo --dry_run value\n    '''"
        )
        body = examples(double, entry(2, "Use mode two.", "demo --mode two value"))
        self.assertIn("HELP.EXAMPLES.CODE_BLOCK", rules(body))

    def test_entry_without_owner_invocation_fails(self) -> None:
        body = examples(
            entry(1, "Use mode one.", "other --mode one value"),
            entry(2, "Use mode two.", "demo --mode two value"),
        )
        self.assertIn("HELP.EXAMPLES.OWNER_INVOCATION", rules(body))

    def test_hidden_alias_in_example_fails(self) -> None:
        body = examples(
            entry(1, "Request hidden help.", "demo --hlp"),
            entry(2, "Use mode two.", "demo --mode two value"),
        )
        self.assertIn("HELP.EXAMPLES.ALIAS_VISIBILITY", rules(body))

    def test_retired_or_rejected_alias_in_example_fails(self) -> None:
        body = examples(
            entry(1, "Use a retired spelling.", "demo --old_mode one value"),
            entry(2, "Use mode two.", "demo --mode two value"),
        )
        self.assertIn("HELP.EXAMPLES.ALIAS_ACCEPTANCE", rules(body))

    def test_transposed_alias_in_example_fails(self) -> None:
        body = examples(
            entry(1, "Use a transposed spelling.", "demo --run_dry value"),
            entry(2, "Use mode two.", "demo --mode two value"),
        )
        self.assertIn("HELP.EXAMPLES.ALIAS_ACCEPTANCE", rules(body))

    def test_public_short_alias_in_example_passes(self) -> None:
        body = examples(
            entry(1, "Display help.", "demo -h"),
            entry(2, "Use public dry-run spelling.", "demo -dr value"),
        )
        self.assertEqual(analyze(body).findings, [])

    def test_canonical_long_alias_in_example_passes(self) -> None:
        body = examples(
            entry(1, "Use mode one.", "demo --mode one value"),
            entry(2, "Use canonical dry-run spelling.", "demo --dry_run value"),
        )
        self.assertEqual(analyze(body).findings, [])

    def test_setup_options_remain_ignored_with_runtime_requirement_bullets(self) -> None:
        source = """function demo() {
    local show_help
    show_help=$(cat << EOM
Usage
-----
  demo value

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - rm (when cleaning temporary output)
EOM
    )
}
"""
        accepted = accepted_function_aliases(source, "demo")
        body = examples(
            entry(1, "Create one output directory.", "mkdir -p first", "demo first"),
            entry(2, "Create another output directory.", "mkdir -p second", "demo second third"),
        )
        analysis = analyze_help_document(
            body,
            owner="demo",
            accepted_aliases=accepted,
            public_aliases=set(),
            hidden_aliases=set(),
        )
        self.assertEqual(accepted, set())
        self.assertEqual(analysis.findings, [])

    def test_exact_duplicate_invocations_fail(self) -> None:
        body = examples(
            entry(1, "First.", "demo --mode one value"),
            entry(2, "Second.", "demo --mode one value"),
        )
        self.assertIn("HELP.EXAMPLES.DUPLICATE", rules(body))

    def test_cosmetic_only_duplicate_signatures_fail(self) -> None:
        body = examples(
            entry(1, "First.", "demo   --mode one   value"),
            entry(2, "Second.", "demo --mode one value  # same branch"),
        )
        self.assertIn("HELP.EXAMPLES.SIGNATURE_DUPLICATE", rules(body))

    def test_material_distinctness_candidate_is_emitted(self) -> None:
        body = examples(
            entry(1, "Use mode one.", "demo --mode one value"),
            entry(2, "Use mode two.", "demo --mode two value"),
        )
        analysis = analyze(body)
        self.assertTrue(
            any(review.rule_id == "HELP.EXAMPLES.REVIEW.MATERIAL_DISTINCTNESS" for review in analysis.reviews)
        )

    def test_explicit_success_and_rejection_outcomes_are_distinct(self) -> None:
        body = examples(
            entry(1, "Accept valid streamed input.", "printf '%s\\n' valid | demo"),
            entry(
                2,
                "Handle rejected streamed input.",
                "if ! printf '%s\\n' invalid | demo; then",
                "    printf '%s\\n' 'rejected as expected'",
                "fi",
            ),
        )
        analysis = analyze(body)
        self.assertEqual(analysis.findings, [])
        self.assertEqual(analysis.reviews, [])

    def test_semantic_candidate_generation_is_deterministic(self) -> None:
        body = examples(
            entry(1, "Use mode one.", "demo --mode one value"),
            entry(2, "Use mode two.", "demo --mode two value"),
        )
        first = [review.as_dict() for review in analyze(body).reviews]
        second = [review.as_dict() for review in analyze(body).reviews]
        self.assertEqual(first, second)

    def test_signature_preserves_harness_lifecycle_shapes(self) -> None:
        code = (
            'tmp="$(mktemp -d)"',
            "trap 'rm -rf -- \"${tmp}\"' EXIT",
            "if ! demo value; then",
            "    record_skip 'dependency unavailable'",
            "fi",
        )
        _, signature, _, _, outcome = invocation_facts(code, "demo")
        self.assertEqual(signature[7], ("assignment",))
        self.assertEqual(signature[8], ("if",))
        self.assertEqual(outcome, ("guarded_rejection",))
        self.assertEqual(signature[10], ("temporary_path",))
        self.assertEqual(signature[11], ("skip", "expected_failure"))
        self.assertEqual(signature[12], ("exit_trap", "path_cleanup"))

    def test_signature_ignores_shellcheck_source_comments(self) -> None:
        code = (
            "# shellcheck source=tests/support/test_helpers.sh",
            "source tests/support/test_helpers.sh",
            "demo BAM input.bam",
        )
        _, signature, _, _, outcome = invocation_facts(code, "demo")
        self.assertEqual(signature[1], ("value", "path"))
        self.assertEqual(signature[2], ("BAM",))
        self.assertEqual(signature[7], ("source",))
        self.assertEqual(outcome, ("direct",))

    def test_signature_ignores_help_rendering_continuation_backslashes(self) -> None:
        compact = (
            "demo --mode one value",
        )
        multiline_source = (
            r"demo --mode one \\",
            "    value",
        )
        rendered = (
            "demo --mode one \\",
            "    value",
        )
        self.assertEqual(
            invocation_facts(compact, "demo")[1],
            invocation_facts(multiline_source, "demo")[1],
        )
        self.assertEqual(
            invocation_facts(compact, "demo")[1],
            invocation_facts(rendered, "demo")[1],
        )

    def test_owner_with_candidate_cannot_be_strict_green(self) -> None:
        analysis = analyze(
            examples(
                entry(1, "Use mode one.", "demo --mode one value"),
                entry(2, "Use mode two.", "demo --mode two value"),
            )
        )
        pending = undispositioned_review_findings(
            "lib/bash/core/demo.sh::demo", analysis.reviews
        )
        self.assertTrue(pending)
        self.assertTrue(all(item.rule_id.endswith("UNDISPOSITIONED") for item in pending))

    def test_details_surface_and_short_advertisement_classify(self) -> None:
        source = """if [[ \"${1:-}\" == --details ]]; then
    detail_execute_demo
elif [[ \"${1:-}\" == --help ]]; then
    help_execute_demo
fi
"""
        ownership = classify_wrapper_source("bin/execute_demo.sh", source)
        self.assertEqual(ownership.classification, "details_full_document")
        self.assertEqual(ownership.full_help_owner, "detail_execute_demo")
        documented = """function help_execute_demo() {
    cat << EOM
Usage
-----
  execute_demo.sh
    [--help] [--details]

Parameters
----------
  -h, --help : flag
    Display concise help.

  -d, --details : flag
    Display full help.
EOM
}
"""
        self.assertTrue(
            short_help_advertises_details(documented, "help_execute_demo")
        )
        self.assertFalse(
            short_help_advertises_details(
                documented.replace("  -d, --details : flag", "  --verbose : flag"),
                "help_execute_demo",
            )
        )

    def test_compatibility_delegation_requires_unambiguous_target(self) -> None:
        valid = 'exec "${BASH}" "${dir_scr}/execute_demo.sh" "$@"\n'
        ambiguous = 'exec "${BASH}" "${dir_scr}/${target}" "$@"\n'
        self.assertEqual(
            classify_wrapper_source("bin/execute_old.sh", valid).classification,
            "compatibility_delegate",
        )
        self.assertEqual(
            classify_wrapper_source("bin/execute_old.sh", ambiguous).classification,
            "no_valid_full_help_owner",
        )

    def test_crosswalk_is_disjoint_and_names_shared_usage_identity(self) -> None:
        report = repository_crosswalk(REPO_ROOT)
        categories = {row["category"]: row for row in report["categories"]}
        shared = categories["shared_heredoc_nonowners"]
        self.assertIn(
            "lib/bash/help/help_execute_compute_signal.sh::<file>@16",
            shared["identities"],
        )

    def test_missing_examples_always_fail_direct_enforcement(self) -> None:
        analysis = analyze_help_document(
            USAGE,
            owner="demo",
            accepted_aliases={"-h", "--hlp", "--help", "--mode", "-dr", "--dry_run"},
            public_aliases={"-h", "--help", "--mode", "-dr", "--dry_run"},
            hidden_aliases={"--hlp"},
        )
        self.assertEqual(analysis.status, "strict_violation")
        self.assertIn("HELP.EXAMPLES.REQUIRED", {row.rule_id for row in analysis.findings})

    def test_serialized_entry_structure_resolves_se_pe_distinction(self) -> None:
        body = examples(
            entry(
                1,
                "Submit one single-end entry.",
                "submit_trim_fastqs.sh --csv_fil_in reads/sample.fastq.gz",
            ),
            entry(
                2,
                "Submit two paired-end entries.",
                "submit_trim_fastqs.sh --csv_fil_in reads/a_R1.fastq.gz,reads/a_R2.fastq.gz;reads/b_R1.fastq.gz,reads/b_R2.fastq.gz",
            ),
        )
        analysis = analyze_help_document(
            body,
            owner="submit_trim_fastqs.sh",
            accepted_aliases={"--csv_fil_in"},
            public_aliases={"--csv_fil_in"},
            hidden_aliases=set(),
        )
        self.assertEqual(analysis.findings, [])
        self.assertEqual(analysis.reviews, [])
        self.assertNotEqual(
            analysis.examples[0].signature[5],
            analysis.examples[1].signature[5],
        )

    def test_repository_owner_universe_is_directly_strict(self) -> None:
        result = scan_repository(REPO_ROOT)
        functions = {
            row["identity"]: row
            for row in result.inventory
            if row.get("kind") == "function"
        }
        scripts = {
            row["identity"]: row
            for row in result.inventory
            if row.get("kind") == "script"
        }
        categories = {
            row["category"]: row
            for row in repository_crosswalk(REPO_ROOT)["categories"]
        }
        self.assertEqual(result.findings, [])
        self.assertEqual(result.reviews, [])
        self.assertEqual(result.alias_findings, [])
        self.assertEqual(
            set(functions),
            set(categories["shell_function_heredoc_owners"]["identities"]),
        )
        self.assertEqual(
            set(scripts),
            set(categories["top_level_shell_script_help_owners"]["identities"]),
        )
        self.assertTrue(all("migrated" not in row for row in functions.values()))
        self.assertTrue(all(row["status"] == "strict_green" for row in functions.values()))
        self.assertTrue(all(row["status"] == "strict_green" for row in scripts.values()))

    def test_global_compliance_is_source_derived_positive(self) -> None:
        summary = compliance_summary(scan_repository(REPO_ROOT))
        self.assertTrue(summary["global_compliance"])
        self.assertEqual(summary["remaining"], 0)
        self.assertEqual(summary["script_green"], summary["script_total"])
        self.assertEqual(summary["function_green"], summary["function_total"])

    def test_global_compliance_fails_for_non_green_owner(self) -> None:
        result = scan_repository(REPO_ROOT)
        inventory = [dict(row) for row in result.inventory]
        owner = next(row for row in inventory if row.get("kind") == "function")
        owner["status"] = "strict_violation"
        modified = RepositoryResult(
            result.findings,
            result.reviews,
            inventory,
            result.ownership,
            result.alias_findings,
        )
        self.assertFalse(compliance_summary(modified)["global_compliance"])

    def test_global_compliance_fails_for_alias_finding(self) -> None:
        result = scan_repository(REPO_ROOT)
        modified = RepositoryResult(
            result.findings,
            result.reviews,
            result.inventory,
            result.ownership,
            [object()],
        )
        self.assertFalse(compliance_summary(modified)["global_compliance"])

    def test_global_compliance_fails_for_semantic_candidate(self) -> None:
        analysis = analyze(
            examples(
                entry(1, "Use mode one.", "demo --mode one value"),
                entry(2, "Use mode two.", "demo --mode two value"),
            )
        )
        result = scan_repository(REPO_ROOT)
        modified = RepositoryResult(
            result.findings,
            analysis.reviews,
            result.inventory,
            result.ownership,
            result.alias_findings,
        )
        self.assertFalse(compliance_summary(modified)["global_compliance"])


if __name__ == "__main__":
    unittest.main()
