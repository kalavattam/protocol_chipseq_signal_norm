"""Focused tests for target-family pilot state and review evidence."""

from __future__ import annotations

import copy
import hashlib
import json
import subprocess
import sys
import tempfile
import tomllib
import unittest
from pathlib import Path
from types import SimpleNamespace

AUDIT_DIR = Path(__file__).resolve().parents[3] / "dev" / "audit"
sys.path.insert(0, str(AUDIT_DIR))

from generate_prompts import (
    markdown_table,
    validate_bundle_tables,
    validate_markdown_table,
    write_pilot_bundle,
)
from parse_findings import CommandResult, parse_result
from run import (
    anchor_evidence_windows,
    bounded_markdown_provision,
    controlled_smoke_command,
    controlled_smoke_sources,
    controlled_smoke_target,
    entry_fingerprint,
    focused_diff_excerpt,
    load_target_selection,
    main,
    markdown_fences_balanced,
    normalize_controlled_smoke,
    pilot_report_self_hash,
    run_command,
    selected_git_state,
    validate_controlled_smoke_evidence,
    validate_manifest,
    validate_runtime_evidence,
    verify_runtime_fact_rows,
)
from shell_help_pilot import (
    alias_documentation_findings,
    assignment_facts,
    command_reference_facts,
    documented_alias_facts,
    line_length_facts,
    parser_alias_facts,
    resolve_alias_facts,
    source_style,
    supporting_alignment,
    validate_command_registry,
    wrapper_contract,
)


class ShellHelpPilotTest(unittest.TestCase):
    """Exercise untracked coverage and bounded new-file review evidence."""

    def make_repo(self, root: Path) -> Path:
        subprocess.run(["git", "init", "-q", str(root)], check=True)
        subprocess.run(["git", "-C", str(root), "config", "user.email", "audit@example.test"], check=True)
        subprocess.run(["git", "-C", str(root), "config", "user.name", "Audit Test"], check=True)
        (root / "tracked.txt").write_text("initial\n", encoding="utf-8")
        subprocess.run(["git", "-C", str(root), "add", "tracked.txt"], check=True)
        subprocess.run(["git", "-C", str(root), "commit", "-qm", "initial"], check=True)
        return root

    def no_index(self, root: Path, path: str) -> CommandResult:
        return run_command(["git", "diff", "--no-index", "--check", "--no-ext-diff", "/dev/null", str(root / path)], root)

    def test_untracked_clean_file_is_covered_without_finding(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = self.make_repo(Path(tmp) / "repo")
            (root / "new.txt").write_text("clean\n", encoding="utf-8")
            result = self.no_index(root, "new.txt")
            self.assertEqual(result.returncode, 1)
            self.assertEqual(result.stdout, "")
            self.assertEqual(result.stderr, "")

    def test_untracked_trailing_whitespace_produces_finding(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = self.make_repo(Path(tmp) / "repo")
            (root / "new.txt").write_text("bad \n", encoding="utf-8")
            result = self.no_index(root, "new.txt")
            self.assertNotEqual(result.returncode, 1)
            parsed = parse_result("git_diff_check", result, "new.txt")
            self.assertEqual(parsed[0]["line"], 1)
            self.assertIn("trailing whitespace", str(parsed[0]["message"]))

    def test_target_states_identify_untracked_and_staged_worktree(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = self.make_repo(Path(tmp) / "repo")
            (root / "tracked.txt").write_text("staged\n", encoding="utf-8")
            subprocess.run(["git", "-C", str(root), "add", "tracked.txt"], check=True)
            (root / "tracked.txt").write_text("staged and worktree\n", encoding="utf-8")
            (root / "new.txt").write_text("new\n", encoding="utf-8")
            status = {line[3:]: {"path": line[3:], "git_status": line[:2]} for line in subprocess.run(["git", "-C", str(root), "status", "--porcelain"], text=True, capture_output=True, check=True).stdout.splitlines()}
            self.assertEqual(selected_git_state(root, "new.txt", status)["git_state_labels"], ["untracked"])
            labels = selected_git_state(root, "tracked.txt", status)["git_state_labels"]
            self.assertIn("tracked_staged", labels)
            self.assertIn("tracked_modified", labels)

    def test_untracked_primary_is_explicit_new_file_evidence_and_stale(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            report = Path(tmp) / "report"
            report.mkdir()
            target = {"schema_version": 1, "path": "new.sh", "role": "primary", "git_state_labels": ["untracked"], "content_fingerprint": "sha256:" + "a" * 64}
            records = write_pilot_bundle(report, [target], [], [], [], [], {}, lambda _: "#!/usr/bin/env bash\n", 80)
            body = (report / records[0]["path"]).read_text(encoding="utf-8")
            self.assertIn("git_state: untracked", body)
            self.assertIn("evidence_kind: new_file_source", body)
            self.assertIn(target["content_fingerprint"], body)

    def test_help_module_is_not_misclassified_as_a_top_level_wrapper(self) -> None:
        facts: list[dict[str, object]] = []
        findings: list[dict[str, object]] = []
        wrapper_contract("lib/bash/help/help_submit_download_fastqs.sh", "function help_submit_download_fastqs() { :; }", facts, findings)
        self.assertEqual(facts, [])
        self.assertEqual(findings, [])

    def test_help_body_accepts_owned_comments_before_heredoc(self) -> None:
        from shell_help_pilot import help_body

        text = """function help_tool_alpha() {
    # The owning wrapper initializes the interpolated default.
    # shellcheck disable=SC2154
    cat << EOM
Usage
-----
EOM
}
"""
        self.assertEqual(help_body(text), "Usage\n-----")

    def test_parser_alias_observations_are_policy_neutral(self) -> None:
        text = '''function resolve_dir_scr() {
    case "${1}" in
        -ds|--dir[_-]scr)
            :
            ;;
    esac
}
function parse_args() {
    case "${1}" in
        -h|--help)
            :
            ;;
        -*)
            :
            ;;
    esac
}
- env_protocol (or renamed equivalent)
'''
        facts: list[dict[str, object]] = []
        parser_alias_facts("bin/submit_download_fastqs.sh", text, facts)
        self.assertEqual(facts[0]["topic"], "parser_alias_observations")
        self.assertEqual({item["pattern"] for item in facts[0]["value"]}, {"-ds|--dir[_-]scr", "-h|--help", "-*"})
        self.assertNotIn("visibility", str(facts[0]))
        resolve_alias_facts(facts)
        resolved = {item["alias"]: item for item in facts[-1]["value"]}
        self.assertEqual(resolved["--dir-scr"]["visibility"], "hidden_systematic_compatibility")
        self.assertEqual(resolved["--dir-scr"]["retention"], "indefinite")
        self.assertEqual(resolved["-h"]["visibility"], "public")
        self.assertEqual(resolved["-h"]["retention"], "canonical_or_supported")
        self.assertEqual(resolved["-* ".strip()]["visibility"], "unsupported_fallback")

    def alias_comparison(
        self,
        parser_text: str,
        help_text: str,
        removed: list[str] | None = None,
        owner: str = "bin/tool_alpha.sh",
        source: str = "docs/tool_alpha_help.sh",
    ) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
        """Return generic resolved facts and documentation findings for one interface."""

        facts: list[dict[str, object]] = []
        findings: list[dict[str, object]] = []
        parser_alias_facts(owner, parser_text, facts)
        documented_alias_facts(source, help_text, facts)
        resolve_alias_facts(facts)
        alias_documentation_findings(
            facts,
            [{"path": owner, "documentation_source": source, "removed_aliases": removed or []}],
            findings,
        )
        return facts, findings

    def test_public_short_alias_is_required_in_help(self) -> None:
        facts, findings = self.alias_comparison(
            "function parse_args() {\n    case x in\n        -q|--queue_name)\n            :\n            ;;\n    esac\n}\n",
            "function help_tool_alpha() {\ncat << EOM\nParameters\n----------\n  --queue_name : str\nEOM\n}\n",
        )
        resolved = {row["alias"]: row for item in facts if item["topic"] == "resolved_alias_classifications" for row in item["value"]}
        self.assertEqual(resolved["-q"]["visibility"], "public")
        self.assertEqual([item["evidence"] for item in findings], ["-q"])

    def test_documented_canonical_long_alias_has_no_finding(self) -> None:
        _, findings = self.alias_comparison(
            "function parse_args() {\n    case x in\n        --queue_name)\n            :\n            ;;\n    esac\n}\n",
            "function help_tool_alpha() {\ncat << EOM\nParameters\n----------\n  --queue_name : str\nEOM\n}\n",
        )
        self.assertEqual(findings, [])

    def test_all_public_parameter_aliases_have_no_finding(self) -> None:
        _, findings = self.alias_comparison(
            "function parse_args() {\n    case x in\n        -q|--queue|--queue_name)\n            :\n            ;;\n    esac\n}\n",
            "function help_tool_alpha() {\ncat << EOM\nParameters\n----------\n  -q, --queue, --queue_name : str\nEOM\n}\n",
        )
        self.assertEqual(findings, [])

    def test_long_only_help_row_omits_public_short_alias(self) -> None:
        facts, findings = self.alias_comparison(
            "function parse_args() {\n    case x in\n        -q|--queue_name)\n            :\n            ;;\n    esac\n}\n",
            "function help_tool_alpha() {\ncat << EOM\nParameters\n----------\n  --queue_name : str\nEOM\n}\n",
        )
        documented = next(item for item in facts if item["topic"] == "documented_aliases")
        self.assertEqual(documented["value"], ["--queue_name"])
        self.assertEqual([item["evidence"] for item in findings], ["-q"])

    def test_help_metavariable_does_not_prevent_exact_alias_matching(self) -> None:
        _, findings = self.alias_comparison(
            "function parse_args() {\n    case x in\n        -q|--queue_name)\n            :\n            ;;\n    esac\n}\n",
            "function help_tool_alpha() {\ncat << EOM\nParameters\n----------\n  -q, --queue_name FILE : file\nEOM\n}\n",
        )
        self.assertEqual(findings, [])

    def test_externalized_help_is_bound_to_owning_wrapper(self) -> None:
        facts, findings = self.alias_comparison(
            "function parse_args() {\n    case x in\n        --queue_name)\n            :\n            ;;\n    esac\n}\n",
            "function help_tool_alpha() {\ncat << EOM\nParameters\n----------\n  --queue_name : str\nEOM\n}\n",
        )
        association = next(item for item in facts if item["topic"] == "documentation_source_associations")
        self.assertEqual(association["value"]["documentation_source"], "docs/tool_alpha_help.sh")
        self.assertEqual(findings, [])

    def test_hidden_systematic_hyphen_alias_is_not_required_in_help(self) -> None:
        facts, findings = self.alias_comparison(
            "function parse_args() {\n    case x in\n        --queue[_-]name)\n            :\n            ;;\n    esac\n}\n",
            "function help_tool_alpha() {\ncat << EOM\nParameters\n----------\n  --queue_name : str\nEOM\n}\n",
        )
        resolved = {row["alias"]: row for item in facts if item["topic"] == "resolved_alias_classifications" for row in item["value"]}
        self.assertEqual(resolved["--queue-name"]["visibility"], "hidden_systematic_compatibility")
        self.assertEqual(findings, [])

    def test_hidden_long_help_alias_is_not_required_in_help(self) -> None:
        facts, findings = self.alias_comparison(
            "function parse_args() {\n    case x in\n        -h|--help|--hlp)\n            :\n            ;;\n    esac\n}\n",
            "function help_tool_alpha() {\ncat << EOM\nParameters\n----------\n  -h, --help : flag\nEOM\n}\n",
        )
        resolved = {row["alias"]: row for item in facts if item["topic"] == "resolved_alias_classifications" for row in item["value"]}
        self.assertEqual(resolved["-h"]["visibility"], "public")
        self.assertEqual(resolved["--hlp"]["visibility"], "hidden_legacy_compatibility")
        self.assertEqual(findings, [])

    def test_documented_hidden_alias_fails(self) -> None:
        _, findings = self.alias_comparison(
            "function parse_args() {\n    case x in\n        --queue[_-]name)\n            :\n            ;;\n    esac\n}\n",
            "function help_tool_alpha() {\ncat << EOM\nParameters\n----------\n  --queue_name, --queue-name : str\nEOM\n}\n",
        )
        self.assertEqual([item["evidence"] for item in findings], ["--queue-name"])

    def test_documented_removed_alias_fails(self) -> None:
        _, findings = self.alias_comparison(
            "function parse_args() {\n    case x in\n        --queue_name)\n            :\n            ;;\n    esac\n}\n",
            "function help_tool_alpha() {\ncat << EOM\nParameters\n----------\n  -i, --queue_name : str\nEOM\n}\n",
            removed=["-i"],
        )
        self.assertEqual([item["evidence"] for item in findings], ["-i"])

    def test_duplicate_documented_alias_is_reported(self) -> None:
        _, findings = self.alias_comparison(
            "function parse_args() {\n    case x in\n        -q|--queue_name)\n            :\n            ;;\n    esac\n}\n",
            "function help_tool_alpha() {\ncat << EOM\nParameters\n----------\n  -q, -q, --queue_name : str\nEOM\n}\n",
        )
        self.assertTrue(any("duplicate" in str(item["message"]) for item in findings))

    def test_removed_aliases_are_recorded_as_rejected_policy(self) -> None:
        facts, findings = self.alias_comparison(
            "function parse_args() {\n    case x in\n        --queue)\n            :\n            ;;\n    esac\n}\n",
            "function help_tool_alpha() {\ncat << EOM\nParameters\n----------\n  --queue : str\nEOM\n}\n",
            removed=["-i", "--infile", "-eo"],
        )
        removed = next(item for item in facts if item["topic"] == "removed_alias_classifications")
        self.assertEqual(
            removed["value"],
            [
                {"alias": "--infile", "acceptance": "rejected", "retention": "removed"},
                {"alias": "-eo", "acceptance": "rejected", "retention": "removed"},
                {"alias": "-i", "acceptance": "rejected", "retention": "removed"},
            ],
        )
        resolved = {row["alias"] for item in facts if item["topic"] == "resolved_alias_classifications" for row in item["value"]}
        self.assertTrue({row["alias"] for row in removed["value"]}.isdisjoint(resolved))
        self.assertEqual(findings, [])

    def test_missing_public_alias_still_produces_a_finding(self) -> None:
        _, findings = self.alias_comparison(
            "function parse_args() {\n    case x in\n        --queue|--public_extra)\n            :\n            ;;\n    esac\n}\n",
            "function help_tool_alpha() {\ncat << EOM\nParameters\n----------\n  --queue : str\nEOM\n}\n",
        )
        self.assertEqual([item["evidence"] for item in findings], ["--public_extra"])

    def test_alias_adapter_is_filename_and_family_independent(self) -> None:
        facts, findings = self.alias_comparison(
            "function parse_args() {\n    case x in\n        -z|--zeta_value)\n            :\n            ;;\n    esac\n}\n",
            "function help_completely_unrelated() {\ncat << EOM\nParameters\n----------\n  -z, --zeta_value VALUE : str\nEOM\n}\n",
            owner="utilities/arbitrary_command",
            source="reference/independent_documentation",
        )
        association = next(item for item in facts if item["topic"] == "documentation_source_associations")
        self.assertEqual(association["path"], "utilities/arbitrary_command")
        self.assertEqual(findings, [])

    def test_current_download_pilot_has_no_canonical_or_hidden_documentation_findings(self) -> None:
        root = Path(__file__).resolve().parents[3]
        config = json.loads((root / "dev/config/pilots/download_fastqs_shell_help.json").read_text(encoding="utf-8"))
        interface_paths = {item["path"] for item in config["interfaces"]}
        source_paths = {item["documentation_source"] for item in config["interfaces"]}
        snapshot = {
            "schema_version": 1,
            "targets": [{"path": path, "role": "primary", "content": (root / path).read_text(encoding="utf-8")} for path in sorted(interface_paths)],
            "context": [],
            "documentation_sources": [{"path": path, "content": (root / path).read_text(encoding="utf-8")} for path in sorted(source_paths)],
            "adapter_config": {"interfaces": config["interfaces"]},
        }
        with tempfile.TemporaryDirectory() as tmp:
            snapshot_path = Path(tmp) / "snapshot.json"
            snapshot_path.write_text(json.dumps(snapshot), encoding="utf-8")
            proc = subprocess.run([sys.executable, str(AUDIT_DIR / "shell_help_pilot.py"), "--mode", "interface_facts", "--snapshot", str(snapshot_path)], text=True, capture_output=True, check=True)
        payload = json.loads(proc.stdout)
        alias_findings = [item for item in payload["findings"] if item["topic"] == "aliases"]
        self.assertEqual(alias_findings, [])
        classifications = {
            row["path"]: {item["alias"]: item["acceptance"] for item in row["value"]}
            for row in payload["facts"]
            if row["topic"] == "removed_alias_classifications"
        }
        self.assertEqual(classifications["bin/execute_download_fastqs.sh"], {"-i": "rejected", "--infile": "rejected", "-eo": "rejected"})

    def test_interface_assignments_exclude_runtime_expression(self) -> None:
        text = '''function init_arg_defs() {
    env_nam="env_protocol"
    threads=1
}
function parse_args() {
    env_nam="${2}"
}
function run_jobs() {
    tmp="$(mktemp)"
}
'''
        facts: list[dict[str, object]] = []
        assignment_facts("bin/execute_download_fastqs.sh", text, facts)
        rows = facts[0]["value"]
        self.assertEqual({row["name"] for row in rows}, {"env_nam", "threads"})
        self.assertNotIn("$(mktemp)", str(rows))

    def test_bundle_fences_source_comments_and_has_human_sections(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            report = Path(tmp) / "report"
            report.mkdir()
            target = {"schema_version": 1, "path": "new.sh", "role": "primary", "git_state_labels": ["untracked"], "content_fingerprint": "sha256:" + "a" * 64, "checks_run": [], "findings_count": 0}
            facts = [{"topic": "line_length_candidates", "path": "new.sh", "value": [{"line": 1, "length": 99, "excerpt": "1: # not a heading", "exception_cues": [], "classification_status": "semantic/manual"}]}]
            records = write_pilot_bundle(report, [target], [], facts, [], [], {}, lambda _: "# not a heading\necho hi\n", 80, {"aliases": "## Rule heading\n# comment"}, {})
            body = (report / records[0]["path"]).read_text(encoding="utf-8")
            self.assertIn("## Executive summary", body)
            self.assertIn("## Supporting-test alignment", body)
            self.assertIn("```bash\n# not a heading", body)
            in_fence = False
            outside = []
            for line in body.splitlines():
                if line.startswith("```"):
                    in_fence = not in_fence
                    continue
                if not in_fence:
                    outside.append(line)
            self.assertNotIn("# not a heading", outside)

    def test_wrapper_examples_and_registry_runtime_case_are_findings(self) -> None:
        text = """function help_submit_download_fastqs() {
cat << EOM
Usage
-----
Parameters
----------
Notes
-----
  - Bash >= 4.4
Examples
--------
  1. One.
    '''bash
    submit_download_fastqs.sh --dir_scr VALUE a b c d e f g h
    '''
EOM
}
"""
        facts: list[dict[str, object]] = []
        findings: list[dict[str, object]] = []
        source_style("lib/bash/help/help_submit_download_fastqs.sh", text, facts, findings)
        self.assertTrue(any("two materially" in item["message"] for item in findings))
        callables, concepts = validate_command_registry({"schema_version": 1, "commands": [{"callable": "bash", "conceptual_names": ["Bash"]}]})
        questions: list[dict[str, object]] = []
        command_reference_facts("lib/bash/help/help_submit_download_fastqs.sh", text, callables, concepts, {"runtime_version", "external_programs"}, facts, findings, questions)
        self.assertTrue(any(item.get("topic") == "runtime_requirements" for item in findings))
        self.assertEqual(questions, [])

    def test_scoped_supporting_options_exclude_unrelated_bind(self) -> None:
        text = '''TEST_NAME="x"
python -m http.server --bind 127.0.0.1
bash bin/execute_download_fastqs.sh --env_nam env --fil_in x --dir_out y --dir_sym z
then
'''
        facts: list[dict[str, object]] = []
        supporting_alignment("tests/integration/local/download_fastqs/test_execute_download_fastqs_se_local.sh", text, True, facts)
        self.assertNotIn("--bind", facts[0]["value"]["public_options_invoked"])

    def test_line_length_heredoc_is_not_no_cue_candidate(self) -> None:
        text = '''function help_x() {
cat << EOM
This user-facing help line is deliberately long enough to exceed the configured source preference threshold.
EOM
}
'''
        facts: list[dict[str, object]] = []
        line_length_facts("lib/bash/help/help_x.sh", text, facts)
        candidate = facts[0]["value"][0]
        self.assertEqual(candidate["location_kind"], "user_facing_help_heredoc")
        self.assertEqual(candidate["source_checker_relationship"], "excluded_by_documented_source_checker_policy")

    def test_verify_marks_an_untracked_target_stale(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            public = self.make_repo(Path(tmp) / "public")
            private = self.make_repo(Path(tmp) / "private")
            (public / "new.txt").write_text("before\n", encoding="utf-8")
            (public / "baseline.json").write_text("{}\n", encoding="utf-8")
            (public / "rules.toml").write_text("schema_version = 1\n", encoding="utf-8")
            fingerprint, _ = entry_fingerprint(public, "new.txt")
            base_fp, _ = entry_fingerprint(public, "baseline.json")
            rules_fp, _ = entry_fingerprint(public, "rules.toml")
            report = Path(tmp) / "report"
            report.mkdir()
            for name in ("ledger.ndjson", "artifacts.ndjson", "findings.ndjson", "checks.ndjson", "prompts.ndjson", "facts.ndjson", "policy_questions.ndjson", "adapter_limitations.ndjson"):
                (report / name).write_text("", encoding="utf-8")
            (report / "cohort_progress.json").write_text("{}\n", encoding="utf-8")
            (report / "rule_manifest.json").write_text("{}\n", encoding="utf-8")
            (report / "targets.ndjson").write_text(json.dumps({"schema_version": 1, "run_id": "test", "path": "new.txt", "role": "primary", "content_fingerprint": fingerprint, "git_state_labels": ["untracked"], "whitespace_coverage": []}) + "\n", encoding="utf-8")
            (report / "semantic_review.ndjson").write_text(json.dumps({"schema_version": 1, "run_id": "test", "bundle_id": "test", "path": "semantic_review/test.md", "target_fingerprints": {"new.txt": fingerprint}}) + "\n", encoding="utf-8")
            (report / "pilot_report.json").write_text(json.dumps({"schema_version": 1}) + "\n", encoding="utf-8")
            (report / "run.json").write_text(json.dumps({"status": "completed", "report_format_version": 2, "baseline_path": "baseline.json", "baseline_fingerprint": base_fp, "rule_manifest_path": "rules.toml", "rule_manifest_fingerprint": rules_fp, "consumed_evidence": []}), encoding="utf-8")
            args = ["--public-root", str(public), "--private-root", str(private), "--verify", str(report)]
            self.assertEqual(main(args), 0)
            (public / "new.txt").write_text("after\n", encoding="utf-8")
            self.assertEqual(main(args), 3)

    def test_renderer_uses_artifact_findings_and_omits_resolved_decisions(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            report = Path(tmp) / "report"
            report.mkdir()
            target = {"schema_version": 1, "path": "bin/submit_download_fastqs.sh", "role": "primary", "git_state_labels": ["tracked_modified"], "content_fingerprint": "sha256:" + "b" * 64, "checks_run": [], "findings_count": 1}
            findings = [{"topic": "types_and_placeholders", "path": target["path"], "message": "noncanonical placeholder", "line": 1}]
            questions = [{"topic": "other", "question": "Resolve a remaining question.", "paths": [target["path"]]}]
            record = write_pilot_bundle(report, [target], [], [], findings, questions, {}, lambda _: "echo x\n", 80)[0]
            body = (report / record["path"]).read_text(encoding="utf-8")
            self.assertIn("1 deterministic findings:\n- 1 noncanonical placeholder:", body)
            self.assertIn("Resolve a remaining question.", body)
            self.assertNotIn("Resolve runtime-requirements capitalization policy", body)
            self.assertNotIn("Decide whether submit help benefits", body)

    def test_injected_findings_change_each_rendered_topic_result(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            report = Path(tmp) / "report"
            report.mkdir()
            target = {"schema_version": 1, "path": "bin/execute_download_fastqs.sh", "role": "primary", "git_state_labels": ["tracked_modified"], "content_fingerprint": "sha256:" + "d" * 64, "checks_run": [], "findings_count": 10}
            topics = ("help_ownership", "aliases", "examples", "types_and_placeholders", "requiredness_defaults", "runtime_requirements", "dir_scr", "environment_handling", "stale_naming")
            findings = [{"topic": topic, "path": target["path"], "message": f"injected-{topic}", "line": None} for topic in topics]
            record = write_pilot_bundle(report, [target], [], [], findings, [], {}, lambda _: "echo x\n", 80)[0]
            body = (report / record["path"]).read_text(encoding="utf-8")
            for topic in topics:
                self.assertIn(f"1 injected-{topic}:", body)
            self.assertIn("Not applicable; this topic is semantic-review only.", body)

    def test_alias_comparison_keeps_fallback_and_compatibility_separate(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            report = Path(tmp) / "report"
            report.mkdir()
            targets = [{"schema_version": 1, "path": "bin/submit_download_fastqs.sh", "role": "primary", "git_state_labels": ["tracked_modified"], "content_fingerprint": "sha256:" + "c" * 64, "checks_run": [], "findings_count": 0}]
            facts = [
                {"topic": "resolved_alias_classifications", "path": targets[0]["path"], "value": [{"alias": "-*", "visibility": "unsupported_fallback", "retention": "not_an_option", "locations": [{"function": "parse_args", "line": 2}]}, {"alias": "-ds", "visibility": "hidden_short_compatibility", "retention": "indefinite", "locations": [{"function": "parse_args", "line": 1}]}, {"alias": "--verbose", "visibility": "public", "retention": "canonical_or_supported", "locations": [{"function": "parse_args", "line": 1}]}, {"alias": "--hlp", "visibility": "hidden_legacy_compatibility", "retention": "indefinite", "locations": [{"function": "main", "line": 1}]}]},
                {"topic": "documented_aliases", "path": "lib/bash/help/help_submit_download_fastqs.sh", "value": ["-ds"]},
            ]
            record = write_pilot_bundle(report, targets, [], facts, [], [], {}, lambda _: "echo x\n", 80)[0]
            body = (report / record["path"]).read_text(encoding="utf-8")
            self.assertIn("parser fallback; not a supported option", body)
            self.assertIn("hidden compatibility; intentionally undocumented", body)
            self.assertIn("public accepted but undocumented", body)

    def test_alias_results_retain_public_and_hidden_aliases_without_retirement(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            report = Path(tmp) / "report"
            report.mkdir()
            targets = [
                {"schema_version": 1, "path": "bin/execute_download_fastqs.sh", "role": "primary", "git_state_labels": ["tracked_modified"], "content_fingerprint": "sha256:" + "e" * 64, "checks_run": [], "findings_count": 0},
                {"schema_version": 1, "path": "bin/submit_download_fastqs.sh", "role": "primary", "git_state_labels": ["tracked_modified"], "content_fingerprint": "sha256:" + "f" * 64, "checks_run": [], "findings_count": 0},
            ]
            facts = [
                {"topic": "resolved_alias_classifications", "path": targets[0]["path"], "value": [{"alias": "--dry-run", "visibility": "hidden_systematic_compatibility", "retention": "indefinite", "locations": [{"function": "parse_args", "line": 1}]}, {"alias": "--verbose", "visibility": "public", "retention": "canonical_or_supported", "locations": [{"function": "parse_args", "line": 1}]}, {"alias": "--hlp", "visibility": "hidden_legacy_compatibility", "retention": "indefinite", "locations": [{"function": "main", "line": 1}]}]},
                {"topic": "resolved_alias_classifications", "path": targets[1]["path"], "value": [{"alias": "--hlp", "visibility": "hidden_legacy_compatibility", "retention": "indefinite", "locations": [{"function": "main", "line": 1}]}]},
                {"topic": "documented_aliases", "path": "lib/bash/help/help_execute_download_fastqs.sh", "value": ["-h", "--help"]},
            ]
            record = write_pilot_bundle(report, targets, [], facts, [], [], {}, lambda _: "echo x\n", 80)[0]
            body = (report / record["path"]).read_text(encoding="utf-8")
            self.assertIn("hidden compatibility; intentionally undocumented", body)
            self.assertIn("public accepted but undocumented", body)
            self.assertNotIn("proposed retirement", body)
            self.assertNotIn("retire legacy hidden compatibility alias", body)

    def test_untracked_help_evidence_contains_complete_help_body(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            report = Path(tmp) / "report"
            report.mkdir()
            target = {"schema_version": 1, "path": "lib/bash/help/help_submit_x.sh", "role": "primary", "git_state_labels": ["untracked"], "content_fingerprint": "sha256:" + "1" * 64, "checks_run": [], "findings_count": 0}
            source = "function help_submit_x() {\n    cat << EOM\nUsage\n-----\n  submit_x -x <str>\n\nParameters\n----------\n  -x, --example : str\n    Example parameter.\n\nNotes\n-----\n  - Keep this note.\nEOM\n}\n"
            record = write_pilot_bundle(report, [target], [], [], [], [], {}, lambda _: source, 80)[0]
            body = (report / record["path"]).read_text(encoding="utf-8")
            self.assertIn("evidence_kind: new_file_help_heredoc", body)
            self.assertIn("-x, --example : str", body)
            self.assertIn("Notes", body)

    def test_runtime_topic_renders_exact_callable_rule(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            report = Path(tmp) / "report"
            report.mkdir()
            target = {"schema_version": 1, "path": "bin/execute_download_fastqs.sh", "role": "primary", "git_state_labels": ["tracked_modified"], "content_fingerprint": "sha256:" + "2" * 64, "checks_run": [], "findings_count": 0}
            rule = "[HELP.RUNTIME.EXACT_CALLABLE]\nUse exact callable spellings such as `bash`, `bamCompare`, and `multiBigwigSummary`; conceptual prose may use project names."
            record = write_pilot_bundle(report, [target], [], [], [], [], {}, lambda _: "echo x\n", 80, {"runtime requirements": rule})[0]
            body = (report / record["path"]).read_text(encoding="utf-8")
            self.assertIn("[HELP.RUNTIME.EXACT_CALLABLE]", body)
            self.assertIn("`bamCompare`, and `multiBigwigSummary`", body)

    def test_registry_preserves_mixed_case_and_rejects_duplicate_callables(self) -> None:
        registry = {
            "schema_version": 1,
            "commands": [
                {"callable": "bamCompare", "conceptual_names": ["deepTools"]},
                {"callable": "multiBigwigSummary", "conceptual_names": ["deepTools"]},
            ],
        }
        callables, concepts = validate_command_registry(registry)
        self.assertEqual(callables, {"bamCompare", "multiBigwigSummary"})
        self.assertEqual(concepts["deepTools"], {"bamCompare", "multiBigwigSummary"})
        duplicate = copy.deepcopy(registry)
        duplicate["commands"].append({"callable": "bamCompare", "conceptual_names": []})
        with self.assertRaisesRegex(ValueError, "duplicate command registry callable"):
            validate_command_registry(duplicate)

    def test_public_registry_and_rule_have_current_bounded_consumption(self) -> None:
        root = Path(__file__).resolve().parents[3]
        registry = json.loads((root / "dev/config/command_names.json").read_text(encoding="utf-8"))
        callables, _ = validate_command_registry(registry)
        self.assertIn("bamCompare", callables)
        self.assertIn("multiBigwigSummary", callables)
        self.assertNotIn("bamcompare", callables)
        config = tomllib.loads((root / "dev/config/rules.toml").read_text(encoding="utf-8"))
        rules = validate_manifest(config)
        rule = next(
            item
            for item in rules
            if item["rule_id"] == "HELP.RUNTIME.REQUIREMENTS"
        )
        self.assertEqual(
            rule["source_checker"],
            "dev/audit/help_runtime_requirements.py",
        )
        self.assertEqual(rule["coverage_relation"], "subset")
        self.assertEqual(rule["scope"], "repository")
        self.assertEqual(
            rule["parity_test"],
            "tests/unit/dev_audit/test_help_runtime_requirements.py",
        )
        self.assertEqual(
            set(rule["applicable_paths"]),
            {"bin/**/*.sh", "lib/bash/**/*.sh"},
        )

    def test_generic_command_resolution_does_not_guess_unknown_or_ambiguous_names(self) -> None:
        callables, concepts = validate_command_registry({
            "schema_version": 1,
            "commands": [
                {"callable": "bash", "conceptual_names": ["Bash"]},
                {"callable": "bamCompare", "conceptual_names": ["deepTools"]},
                {"callable": "multiBigwigSummary", "conceptual_names": ["deepTools"]},
            ],
        })
        text = """function help_example() {
cat << EOM
Notes
-----
  Bash is the implementation language.
  Runtime requirements:
    External programs:
      - Bash >= 4.4
      - bamCompare
      - deepTools
      - mysteryTool
EOM
}
"""
        facts: list[dict[str, object]] = []
        findings: list[dict[str, object]] = []
        questions: list[dict[str, object]] = []
        command_reference_facts("bin/example.sh", text, callables, concepts, {"runtime_version", "external_programs"}, facts, findings, questions)
        self.assertTrue(any("'bash' instead of conceptual name 'Bash'" in item["message"] for item in findings))
        self.assertFalse(any("bamCompare" in item["message"] for item in findings))
        self.assertTrue(any("ambiguous runtime-command reference 'deepTools'" in item["question"] for item in questions))
        self.assertTrue(any("unknown runtime-command reference 'mysteryTool'" in item["question"] for item in questions))
        resolutions = next(
            item for item in facts
            if item["topic"] == "command_reference_resolutions"
        )["value"]
        self.assertIn("exact_callable", {item["status"] for item in resolutions})
        self.assertNotIn("Bash is the implementation language.", str(facts))

    def test_command_reference_adapter_consumes_snapshot_registry(self) -> None:
        registry = {"schema_version": 1, "commands": [{"callable": "bamCompare", "conceptual_names": ["deepTools bamCompare"]}]}
        source = """function help_example() {
cat << EOM
Notes
-----
  Runtime requirements:
    External programs:
      - deepTools bamCompare
EOM
}
"""
        snapshot = {
            "schema_version": 1,
            "targets": [{"path": "bin/example.sh", "role": "primary", "content": source}],
            "context": [{"path": "registry.json", "content": json.dumps(registry)}],
            "adapter_config": {"registry_path": "registry.json", "reference_scopes": ["external_programs"]},
        }
        with tempfile.TemporaryDirectory() as tmp:
            snapshot_path = Path(tmp) / "snapshot.json"
            snapshot_path.write_text(json.dumps(snapshot), encoding="utf-8")
            proc = subprocess.run(
                [sys.executable, str(AUDIT_DIR / "shell_help_pilot.py"), "--mode", "command_references", "--snapshot", str(snapshot_path)],
                text=True,
                capture_output=True,
                check=True,
            )
        payload = json.loads(proc.stdout)
        self.assertEqual(payload["findings"][0]["message"], "callable reference must use exact registry spelling 'bamCompare' instead of conceptual name 'deepTools bamCompare'")
        self.assertEqual(payload["policy_questions"], [])

    def test_markdown_table_validation_and_escaped_pipes(self) -> None:
        table = markdown_table(["Topic", "Path", "Finding"], [["aliases", "bin/x.sh", "value a | b"]])
        self.assertIn("a \\| b", table)
        validate_markdown_table(table.splitlines())
        with self.assertRaises(ValueError):
            validate_markdown_table(["| A | B | C |", "| --- | --- |", "| 1 | 2 | 3 |"])
        with self.assertRaises(ValueError):
            validate_markdown_table(["| A | B | C |", "| --- | --- | --- |", "| 1 | 2 |"])
        with self.assertRaises(ValueError):
            validate_markdown_table(["| A | B | C |", "| --- | --- | --- |", "| 1 | 2 | 3 | 4 |"])

    def test_bundle_validates_all_rendered_tables_and_groups_line_candidates(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            report = Path(tmp) / "report"
            report.mkdir()
            target = {"schema_version": 1, "path": "bin/execute_download_fastqs.sh", "role": "primary", "git_state_labels": ["tracked_modified"], "content_fingerprint": "sha256:" + "3" * 64, "checks_run": [], "findings_count": 0}
            facts = [{"topic": "line_length_candidates", "path": "tests/smoke/test_x.sh", "value": [{"line": 10, "length": 81, "excerpt": "9: x\n10: log_out_a=\"long\"\n11: x", "exception_cues": [], "location_kind": "ordinary_shell_code", "source_checker_relationship": "independent_pilot_only_candidate", "classification_status": "semantic/manual"}, {"line": 11, "length": 82, "excerpt": "10: x\n11: log_err_a=\"long\"\n12: x", "exception_cues": [], "location_kind": "ordinary_shell_code", "source_checker_relationship": "independent_pilot_only_candidate", "classification_status": "semantic/manual"}]}]
            record = write_pilot_bundle(report, [target], [], facts, [], [], {}, lambda _: "echo x\n", 80, {"--dir_scr": "```txt\nnested fence\n```"})[0]
            body = (report / record["path"]).read_text(encoding="utf-8")
            self.assertGreaterEqual(validate_bundle_tables(body), 5)
            self.assertEqual(record["line_length_rendering"]["detailed_candidate_count"], 2)
            self.assertEqual(record["line_length_rendering"]["grouped_candidate_count"], 1)
            self.assertIn("Not applicable; this topic is semantic-review only.", body)

    def test_supporting_gap_requires_a_real_submit_command_boundary(self) -> None:
        quoted = "TEST_NAME=\"x\"\nassert_pattern_found cfg \"${TEST_BASH} ${ROOT_REPO}/bin/submit_download_fastqs.sh\"\nbash bin/execute_download_fastqs.sh --fil_in x --dir_out y --dir_sym z\nthen\n"
        facts: list[dict[str, object]] = []
        supporting_alignment("tests/integration/parallel/download_fastqs/test_execute_download_fastqs_parallel.sh", quoted, True, facts)
        self.assertFalse(facts[0]["value"]["submit_wrapper_directly_invoked"])
        self.assertEqual(facts[0]["value"]["apparent_gap_or_uncertainty"], "No direct submit-wrapper assertion.")
        direct = "TEST_NAME=\"x\"\nrun_capture x log \"${TEST_BASH}\" \"${ROOT_REPO}/bin/submit_download_fastqs.sh\" arg\nbash bin/execute_download_fastqs.sh --fil_in x --dir_out y --dir_sym z\nthen\n"
        facts = []
        supporting_alignment("tests/integration/parallel/download_fastqs/test_execute_download_fastqs_parallel.sh", direct, True, facts)
        self.assertTrue(facts[0]["value"]["submit_wrapper_directly_invoked"])
        self.assertIn("line 2", facts[0]["value"]["apparent_gap_or_uncertainty"])

    def test_focused_diff_uses_complete_line_windows_for_oversized_hunks(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = self.make_repo(Path(tmp) / "repo")
            (root / "tracked.txt").write_text("\n".join(f"before-{index}" for index in range(80)) + "\n", encoding="utf-8")
            subprocess.run(["git", "-C", str(root), "add", "tracked.txt"], check=True)
            subprocess.run(["git", "-C", str(root), "commit", "-qm", "large baseline"], check=True)
            changed = [f"after-{index}" for index in range(80)]
            changed[40] = "dir_scr=changed"
            (root / "tracked.txt").write_text("\n".join(changed) + "\n", encoding="utf-8")
            excerpt = focused_diff_excerpt(root, "tracked.txt", 400)
            self.assertIn("@@", excerpt)
            self.assertIn("+dir_scr=changed\n", excerpt)
            self.assertIn("[earlier lines omitted from this hunk]", excerpt)
            self.assertIn("[later lines omitted from this hunk]", excerpt)
            self.assertTrue(excerpt.endswith("\n"))

    def test_focused_diff_reports_when_no_complete_line_window_fits(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = self.make_repo(Path(tmp) / "repo")
            (root / "tracked.txt").write_text("dir_scr=" + "x" * 800 + "\n", encoding="utf-8")
            excerpt = focused_diff_excerpt(root, "tracked.txt", 300)
            self.assertIn("[relevant hunk omitted: no complete-line window fits excerpt budget]", excerpt)
            self.assertTrue(excerpt.endswith("\n"))

    def test_execute_anchor_windows_require_environment_and_final_lifecycle_source(self) -> None:
        root = Path(__file__).resolve().parents[3]
        target_content = {
            path: (root / path).read_text(encoding="utf-8")
            for path in (
                "bin/execute_download_fastqs.sh",
                "bin/submit_download_fastqs.sh",
                "lib/bash/help/help_execute_download_fastqs.sh",
                "lib/bash/help/help_submit_download_fastqs.sh",
            )
        }
        windows = anchor_evidence_windows(root, target_content)["bin/execute_download_fastqs.sh"]
        self.assertIn('handle_env "${env_nam}"', windows["environment_resolution"]["rendered_evidence"])
        self.assertIn("setup_env               || return 1", windows["final_dispatch"]["rendered_evidence"])
        self.assertIn("run_jobs                || return 1", windows["final_dispatch"]["rendered_evidence"])

    def test_submit_and_environment_anchor_windows_render_named_behavior(self) -> None:
        root = Path(__file__).resolve().parents[3]
        target_content = {
            path: (root / path).read_text(encoding="utf-8")
            for path in (
                "bin/execute_download_fastqs.sh",
                "bin/submit_download_fastqs.sh",
                "lib/bash/help/help_execute_download_fastqs.sh",
                "lib/bash/help/help_submit_download_fastqs.sh",
            )
        }
        windows = anchor_evidence_windows(root, target_content)
        submit = windows["bin/submit_download_fastqs.sh"]
        self.assertIn("-h|--hlp|--help)", submit["early_help"]["rendered_evidence"])
        self.assertIn("help_submit_download_fastqs >&2", submit["early_help"]["rendered_evidence"])
        self.assertIn("return 0", submit["early_help"]["rendered_evidence"])
        self.assertIn("-ds|--dir[_-]scr)", submit["bootstrap_dir_scr"]["rendered_evidence"])
        self.assertIn('printf "%s\\n" "${args[i + 1]}"', submit["bootstrap_dir_scr"]["rendered_evidence"])
        self.assertIn("required option '--dir_scr' was not supplied.", submit["bootstrap_dir_scr"]["rendered_evidence"])
        self.assertIn("if [[ ${#args_pos[@]} -ne 8 ]]; then", submit["positional_validation"]["rendered_evidence"])
        self.assertIn('nam_job="${args_pos[7]}"', submit["positional_validation"]["rendered_evidence"])
        for call in ("source_helpers_submit", "parse_args", "validate_args", "run_downloads"):
            self.assertIn(call, submit["worker_after_bootstrap"]["rendered_evidence"])
        environment = windows["lib/bash/help/help_execute_download_fastqs.sh"]["environment_documentation"]
        self.assertIn("Runtime requirements:", environment["rendered_evidence"])
        self.assertIn(
            "A compatible Conda environment providing the listed tools",
            environment["rendered_evidence"],
        )
        self.assertEqual(environment["evidence_kind"], "complete_line_source_window")

    def test_anchor_validation_rejects_deleted_and_descriptive_substitutes(self) -> None:
        root = Path(__file__).resolve().parents[3]
        target_content = {
            path: (root / path).read_text(encoding="utf-8")
            for path in (
                "bin/execute_download_fastqs.sh",
                "bin/submit_download_fastqs.sh",
                "lib/bash/help/help_execute_download_fastqs.sh",
                "lib/bash/help/help_submit_download_fastqs.sh",
            )
        }
        all_evidence = anchor_evidence_windows(root, target_content)

        def target(path: str) -> dict[str, object]:
            return {"schema_version": 1, "path": path, "role": "primary", "git_state_labels": ["tracked_modified"], "content_fingerprint": "sha256:" + "5" * 64, "checks_run": [], "findings_count": 0}

        with tempfile.TemporaryDirectory() as tmp:
            report = Path(tmp) / "report"
            report.mkdir()
            submit_path = "bin/submit_download_fastqs.sh"
            deleted_only = copy.deepcopy(all_evidence[submit_path])
            deleted_only["worker_after_bootstrap"] = {
                "evidence_kind": "complete_line_diff_window",
                "rendered_evidence": "\n".join(f"-{marker}" for marker in (
                    'source_helpers_submit "${0##*/}" "${dir_scr}" || return 1',
                    'parse_args "$@"',
                    "validate_args",
                    "run_downloads",
                )),
            }
            with self.assertRaisesRegex(ValueError, "requires current source evidence"):
                write_pilot_bundle(report, [target(submit_path)], [], [], [], [], {}, target_content.__getitem__, 12000, anchor_evidence={submit_path: deleted_only})

            descriptive_only = copy.deepcopy(all_evidence[submit_path])
            descriptive_only["positional_validation"] = {
                "evidence_kind": "complete_line_source_window",
                "rendered_evidence": "The submit wrapper requires 8 positional arguments.",
            }
            with self.assertRaisesRegex(ValueError, "missing required behavioral anchor group positional_validation"):
                write_pilot_bundle(report, [target(submit_path)], [], [], [], [], {}, target_content.__getitem__, 12000, anchor_evidence={submit_path: descriptive_only})

            help_path = "lib/bash/help/help_execute_download_fastqs.sh"
            insufficient_environment = copy.deepcopy(all_evidence[help_path])
            insufficient_environment["environment_documentation"] = {
                "evidence_kind": "complete_line_source_window",
                "rendered_evidence": "Conda/Mamba environment to activate.",
            }
            with self.assertRaisesRegex(ValueError, "missing required behavioral anchor group environment_documentation"):
                write_pilot_bundle(report, [target(help_path)], [], [], [], [], {}, target_content.__getitem__, 12000, anchor_evidence={help_path: insufficient_environment})

            execute_path = "bin/execute_download_fastqs.sh"
            deleted_parser = copy.deepcopy(all_evidence[execute_path])
            deleted_parser["argument_parsing"] = {
                "evidence_kind": "complete_line_diff_window",
                "rendered_evidence": "-function parse_args() {\n",
            }
            with self.assertRaisesRegex(ValueError, "missing required behavioral anchor group argument_parsing"):
                write_pilot_bundle(report, [target(execute_path)], [], [], [], [], {}, target_content.__getitem__, 12000, anchor_evidence={execute_path: deleted_parser})

    def test_anchor_metadata_uses_its_exact_rendered_evidence(self) -> None:
        root = Path(__file__).resolve().parents[3]
        target_content = {
            path: (root / path).read_text(encoding="utf-8")
            for path in (
                "bin/execute_download_fastqs.sh",
                "bin/submit_download_fastqs.sh",
                "lib/bash/help/help_execute_download_fastqs.sh",
                "lib/bash/help/help_submit_download_fastqs.sh",
            )
        }
        targets = [
            {"schema_version": 1, "path": path, "role": "primary", "git_state_labels": ["untracked"], "content_fingerprint": "sha256:" + str(index) * 64, "checks_run": [], "findings_count": 0}
            for index, path in enumerate(target_content, 1)
        ]
        with tempfile.TemporaryDirectory() as tmp:
            report = Path(tmp) / "report"
            report.mkdir()
            evidence = anchor_evidence_windows(root, target_content)
            record = write_pilot_bundle(report, targets, [], [], [], [], {}, target_content.__getitem__, 12000, anchor_evidence=evidence)[0]
            expected = evidence["bin/execute_download_fastqs.sh"]["environment_resolution"]["rendered_evidence"]
            metadata = next(item for item in record["behavioral_anchor_coverage"] if item["path"] == "bin/execute_download_fastqs.sh" and item["anchor_group"] == "environment_resolution")
            self.assertEqual(metadata["rendered_evidence"], expected)
            self.assertNotEqual(metadata["rendered_evidence"], focused_diff_excerpt(root, "bin/execute_download_fastqs.sh", 12000))

    def test_missing_execute_anchor_windows_fail_bundle_generation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            report = Path(tmp) / "report"
            report.mkdir()
            target = {"schema_version": 1, "path": "bin/execute_download_fastqs.sh", "role": "primary", "git_state_labels": ["tracked_modified"], "content_fingerprint": "sha256:" + "4" * 64, "checks_run": [], "findings_count": 0}
            with self.assertRaisesRegex(ValueError, "missing required behavioral anchor group"):
                write_pilot_bundle(report, [target], [], [], [], [], {}, lambda _: "echo x\n", 80, anchor_evidence={target["path"]: {}})

    def test_environment_provision_fences_are_balanced_after_bounding(self) -> None:
        text = """### Recommended environment:
Use the project environment.

```txt
env_protocol
```

Use the documented installer.

### External programs:
Do not include this provision.
"""
        provision = bounded_markdown_provision(text, "Recommended environment:", 4000)
        self.assertTrue(markdown_fences_balanced(provision))
        self.assertNotIn("External programs", provision)
        bounded = bounded_markdown_provision(text, "Recommended environment:", len(b"### Recommended environment:\nUse the project environment.\n\n"))
        self.assertTrue(markdown_fences_balanced(bounded))
        self.assertNotIn("```txt", bounded)
        self.assertFalse(markdown_fences_balanced("```txt\nopen\n"))

    def test_pilot_report_self_hash_excludes_its_own_hash_field(self) -> None:
        report = {"schema_version": 1, "run_id": "run", "package_provenance": {"artifact_hashes": {"facts.ndjson": "sha256:facts"}, "pilot_report_self_hash": {"algorithm": "sha256", "canonicalization": "utf8_json_sorted_keys_compact_separators", "excluded_json_pointer": "/package_provenance/pilot_report_self_hash/value", "value": ""}}}
        digest = pilot_report_self_hash(report)
        report["package_provenance"]["pilot_report_self_hash"]["value"] = digest
        self.assertEqual(pilot_report_self_hash(report), digest)

    def controlled_smoke_selection(self) -> tuple[Path, dict[str, object], dict[str, object]]:
        """Load the repository's declarative direct-submit smoke target."""

        root = Path(__file__).resolve().parents[3]
        args = SimpleNamespace(path_values=None, paths_from=Path("dev/config/pilots/download_fastqs_shell_help.json"))
        selection = load_target_selection(args, root)
        return root, selection, controlled_smoke_target(selection)

    def controlled_smoke_capture(self, root: Path, selection: dict[str, object], target: dict[str, object], include_observations: bool = True, **overrides: object) -> dict[str, object]:
        """Build one bounded captured result without launching a smoke test."""

        evidence = target["evidence"]
        messages = [message for group in evidence["required_pass_groups"].values() for message in group]
        if include_observations:
            messages.extend(observation["pass_message"] for observation in evidence.get("optional_observations", {}).values())
        raw = {"stdout": "\n".join([*(f"PASS: {message}" for message in messages), "Summary for submit download-fastqs local: pass=40 fail=0 warn=0 skip=0"]), "stderr": ""}
        raw["digest"] = "sha256:" + hashlib.sha256(json.dumps({"stdout": raw["stdout"], "stderr": raw["stderr"]}, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()
        capture = {"schema_version": 1, "run_id": "runtime-test", "manifest_path": selection["selection_path"], "manifest_fingerprint": entry_fingerprint(root, selection["selection_path"])[0], "target_path": target["path"], "command": controlled_smoke_command(target), "outer_context": {"project_environment_active": False}, "outer_context_ok": True, "source_fingerprints_before": controlled_smoke_sources(root, target), "source_fingerprints_after": controlled_smoke_sources(root, target), "raw": raw, "result": {"exit_status": 0, "timed_out": False, "launch_error": None}}
        capture.update(overrides)
        return capture

    def test_declarative_controlled_smoke_source_coverage_and_runtime_evidence(self) -> None:
        root, selection, target = self.controlled_smoke_selection()
        text = (root / target["path"]).read_text(encoding="utf-8")
        facts: list[dict[str, object]] = []
        supporting_alignment(target["path"], text, True, facts, target["evidence"])
        coverage = facts[0]["value"]
        self.assertEqual(facts[0]["topic"], "controlled_smoke_source_coverage")
        self.assertTrue(all(group["present"] for group in coverage["coverage_groups"].values()))
        with tempfile.TemporaryDirectory() as tmp:
            capture_path = Path(tmp) / "capture.json"
            evidence_path = Path(tmp) / "evidence.json"
            capture_path.write_text(json.dumps(self.controlled_smoke_capture(root, selection, target)), encoding="utf-8")
            args = SimpleNamespace(run_id="runtime-test", normalize_controlled_smoke=capture_path, runtime_evidence_out=evidence_path)
            self.assertEqual(normalize_controlled_smoke(args, root, selection), 0)
            normalized = json.loads(evidence_path.read_text(encoding="utf-8"))
            self.assertEqual(normalized["status"], "passed")
            self.assertEqual(normalized["environment_selection"], {"submit_worker": "not_directly_evidenced"})
            self.assertNotIn("selected_environment", normalized["environment_selection"])
            self.assertEqual(normalized["worker_execution"]["status"], "passed")
            self.assertEqual(normalized["worker_execution"]["environment_evidence"], "not_directly_evidenced")
            self.assertEqual(
                normalized["observations"],
                {
                    name: {"value": item["value"], "evidence": {"kind": "pass_message", "message": item["pass_message"]}}
                    for name, item in target["evidence"]["optional_observations"].items()
                },
            )
            self.assertEqual(validate_runtime_evidence(normalized, root, selection, "runtime-test"), normalized)

    def test_generator_environment_cannot_become_submit_worker_environment_evidence(self) -> None:
        root, selection, target = self.controlled_smoke_selection()
        capture = self.controlled_smoke_capture(root, selection, target)
        capture["outer_context"] = {"conda_default_env": "env_protocol", "project_environment_active": True}
        with tempfile.TemporaryDirectory() as tmp:
            capture_path = Path(tmp) / "capture.json"
            evidence_path = Path(tmp) / "evidence.json"
            capture_path.write_text(json.dumps(capture), encoding="utf-8")
            args = SimpleNamespace(run_id="runtime-test", normalize_controlled_smoke=capture_path, runtime_evidence_out=evidence_path)
            self.assertEqual(normalize_controlled_smoke(args, root, selection), 0)
            normalized = json.loads(evidence_path.read_text(encoding="utf-8"))
            self.assertEqual(normalized["environment_selection"], {"submit_worker": "not_directly_evidenced"})
            self.assertNotIn("selected_environment", normalized["environment_selection"])
            self.assertNotIn("outer_invocation_context", normalized)

    def test_optional_observations_are_not_required_for_valid_core_evidence(self) -> None:
        root, selection, target = self.controlled_smoke_selection()
        with tempfile.TemporaryDirectory() as tmp:
            capture_path = Path(tmp) / "capture.json"
            evidence_path = Path(tmp) / "evidence.json"
            capture_path.write_text(json.dumps(self.controlled_smoke_capture(root, selection, target, include_observations=False)), encoding="utf-8")
            args = SimpleNamespace(run_id="runtime-test", normalize_controlled_smoke=capture_path, runtime_evidence_out=evidence_path)
            self.assertEqual(normalize_controlled_smoke(args, root, selection), 0)
            normalized = json.loads(evidence_path.read_text(encoding="utf-8"))
            self.assertEqual(normalized["status"], "passed")
            self.assertNotIn("observations", normalized)
            self.assertEqual(validate_runtime_evidence(normalized, root, selection, "runtime-test"), normalized)

    def test_generic_controlled_smoke_schema_does_not_require_observations(self) -> None:
        _, _, target = self.controlled_smoke_selection()
        evidence = copy.deepcopy(target["evidence"])
        evidence.pop("optional_observations")
        self.assertEqual(validate_controlled_smoke_evidence(evidence), evidence)

    def test_controlled_smoke_core_evidence_fails_closed(self) -> None:
        root, selection, target = self.controlled_smoke_selection()
        cases: list[tuple[str, dict[str, object]]] = []
        failed = self.controlled_smoke_capture(root, selection, target)
        failed["result"] = {"exit_status": 1, "timed_out": False, "launch_error": None}
        cases.append(("failed", failed))
        skipped = self.controlled_smoke_capture(root, selection, target)
        skipped_raw = {"stdout": "Summary for submit download-fastqs local: pass=1 fail=0 warn=0 skip=1", "stderr": ""}
        skipped_raw["digest"] = "sha256:" + hashlib.sha256(json.dumps(skipped_raw, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()
        cases.append(("skipped", {**skipped, "raw": skipped_raw}))
        stale = self.controlled_smoke_capture(root, selection, target)
        stale["manifest_fingerprint"] = "sha256:stale"
        cases.append(("stale", stale))
        malformed = self.controlled_smoke_capture(root, selection, target)
        malformed["result"] = []
        cases.append(("malformed", malformed))
        source_mismatched = self.controlled_smoke_capture(root, selection, target)
        source_mismatched["source_fingerprints_after"] = []
        cases.append(("source-mismatched", source_mismatched))
        incomplete = self.controlled_smoke_capture(root, selection, target)
        incomplete.pop("result")
        cases.append(("incomplete", incomplete))
        with tempfile.TemporaryDirectory() as tmp:
            for name, capture in cases:
                capture_path = Path(tmp) / f"capture-{name}.json"
                evidence_path = Path(tmp) / f"evidence-{name}.json"
                capture_path.write_text(json.dumps(capture), encoding="utf-8")
                args = SimpleNamespace(run_id="runtime-test", normalize_controlled_smoke=capture_path, runtime_evidence_out=evidence_path)
                with self.subTest(name=name):
                    self.assertEqual(normalize_controlled_smoke(args, root, selection), 4)
                    self.assertEqual(json.loads(evidence_path.read_text(encoding="utf-8"))["status"], "invalid")

    def test_runtime_fact_verification_remains_optional_for_historical_reports(self) -> None:
        self.assertEqual(verify_runtime_fact_rows([], Path(__file__).resolve().parents[3]), [])

    def test_five_file_supplied_package_contract_is_unchanged(self) -> None:
        source = (AUDIT_DIR / "run.py").read_text(encoding="utf-8")
        self.assertIn('supplied = ("semantic_review/download-fastqs-shell-help-pilot.md", "findings.ndjson", "facts.ndjson", "adapter_limitations.ndjson", "pilot_report.json")', source)


if __name__ == "__main__":
    unittest.main()
