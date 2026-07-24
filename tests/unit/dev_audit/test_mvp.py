"""Focused isolated tests for the cleanup-audit MVP."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

AUDIT_DIR = Path(__file__).resolve().parents[3] / "dev" / "audit"
sys.path.insert(0, str(AUDIT_DIR))

from parse_findings import CommandResult, parse_result
from run import (
    BASELINE_COUNT,
    BaselineValidationError,
    EvidenceReader,
    build_inventory,
    canonical_report_manifest,
    cohort_progress,
    entry_fingerprint,
    executable_ledger_paths,
    load_path_relocations,
    main,
    path_matches,
    run_command,
    run_shell_syntax_batch,
    validate_baseline,
    validate_canonical_report_manifest,
    validate_manifest,
)
from shell_validation import (
    MINIMUM_BASH,
    POSIX_BOOTSTRAP,
    resolve_env_bash,
    syntax_command,
)


class AuditMvpTest(unittest.TestCase):
    """Exercise audit data behavior without touching the working repository."""

    def make_repo(self, directory: Path) -> Path:
        """Initialize a small repository with one committed file."""

        subprocess.run(["git", "init", "-q", str(directory)], check=True)
        subprocess.run(["git", "-C", str(directory), "config", "user.email", "audit@example.test"], check=True)
        subprocess.run(["git", "-C", str(directory), "config", "user.name", "Audit Test"], check=True)
        (directory / "tracked.txt").write_text("initial\n", encoding="utf-8")
        subprocess.run(["git", "-C", str(directory), "add", "tracked.txt"], check=True)
        subprocess.run(["git", "-C", str(directory), "commit", "-qm", "initial"], check=True)
        return directory

    def baseline(self) -> dict[str, object]:
        """Return a valid synthetic immutable cohort payload."""

        return {
            "schema_version": 1,
            "cohort": "baseline_cleanup",
            "authored_path_count": BASELINE_COUNT,
            "paths": [
                {"path": f"src/file_{index:03d}.txt", "initial_git_status": "??"}
                for index in range(BASELINE_COUNT)
            ],
        }

    def artifact_config(self) -> dict[str, object]:
        """Return the minimum artifact configuration used by baseline tests."""

        return {
            "artifacts": {
                "bytecode_globs": ["**/*.pyc", "*.pyc"],
                "other_generated_globs": ["artifacts/tests/**"],
            }
        }

    def manifest_rule(self) -> dict[str, object]:
        """Return one complete manifest rule for validation tests."""

        return {
            "rule_id": "RULE",
            "description": "x",
            "normative_document": "x",
            "normative_section": "x",
            "owner_classification": "advisory",
            "source_checker": "x",
            "execution_kind": "independent",
            "execution_role": "independent_evidence",
            "coverage_relation": "independent",
            "covered_scope": "Independent evidence for the declared scope.",
            "remaining_scope": "Owner interpretation remains review-owned.",
            "scope": "per_path",
            "applicable_paths": ["**"],
            "default_severity": "error",
            "blocking": False,
            "semantic_review": False,
            "required_environment": "bash",
            "required_commands": ["bash"],
            "known_side_effects": ["none"],
            "output_paths": [],
            "output_parser": "x",
            "strict_availability": "safe_adapter",
            "current_exclusions_or_allowlists": [],
        }

    def test_parser_and_staged_whitespace_evidence(self) -> None:
        """Normalize both standard checker output forms."""

        bash = parse_result(
            "bash_stderr",
            CommandResult(["bash"], "", "bad.sh: line 7: syntax error", 2),
            "bad.sh",
        )
        self.assertEqual(bash[0]["line"], 7)
        staged = parse_result(
            "git_diff_check",
            CommandResult(["git", "diff", "--cached", "--check"], "tracked.txt:2: trailing whitespace.\n", "", 2),
            ".",
        )
        self.assertEqual(staged[0]["path"], "tracked.txt")
        self.assertEqual(staged[0]["line"], 2)

    def test_private_evidence_detects_content_change(self) -> None:
        """Fingerprint only the private file actually consumed."""

        with tempfile.TemporaryDirectory() as tmp:
            public = self.make_repo(Path(tmp) / "public")
            private = self.make_repo(Path(tmp) / "private")
            (private / "evidence.txt").write_text("before\n", encoding="utf-8")
            reader = EvidenceReader(public, private)
            self.assertEqual(reader.read("private", "evidence.txt", "test"), "before\n")
            (private / "evidence.txt").write_text("after\n", encoding="utf-8")
            self.assertEqual(reader.verify()[0]["verification_status"], "stale")

    def test_cohort_reconciliation_preserves_renamed_and_missing_paths(self) -> None:
        """Keep baseline identity visible when current paths diverge."""

        with tempfile.TemporaryDirectory() as tmp:
            root = self.make_repo(Path(tmp) / "public")
            (root / "current.txt").write_text("current\n", encoding="utf-8")
            baseline = {"authored_path_count": 2, "paths": [{"path": "old.txt"}, {"path": "gone.txt"}]}
            ledger = [{"path": "current.txt", "cohort_origin_path": "old.txt", "review_status": "unreviewed", "disposition": "unresolved", "scope_cohort": "baseline_cleanup"}]
            progress = cohort_progress(root, baseline, ledger)
            states = {row["baseline_path"]: row["state"] for row in progress["baseline_paths"]}
            self.assertEqual(states["old.txt"], "renamed")
            self.assertEqual(states["gone.txt"], "missing_or_unmatched")

    def test_path_relocations_preserve_baseline_identity(self) -> None:
        """Map an explicitly moved dirty path to its baseline identity."""

        with tempfile.TemporaryDirectory() as tmp:
            root = self.make_repo(Path(tmp) / "public")
            target = root / "tests/unit/example/test_worker.py"
            target.parent.mkdir(parents=True)
            target.write_text("# moved test\n", encoding="utf-8")
            config = root / "dev/config/path_relocations.json"
            config.parent.mkdir(parents=True)
            config.write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "exact": {
                            "tests/unit/example/test_worker.py":
                                "dev/test_worker.py"
                        },
                        "basename_origins": [],
                    }
                ),
                encoding="utf-8",
            )
            self.assertEqual(
                load_path_relocations(root, {"dev/test_worker.py"}),
                {
                    "tests/unit/example/test_worker.py":
                        "dev/test_worker.py"
                },
            )

    def test_inventory_separates_baseline_and_audit_files(self) -> None:
        """Classify current dirty paths without putting bytecode in the ledger."""

        with tempfile.TemporaryDirectory() as tmp:
            root = self.make_repo(Path(tmp) / "public")
            (root / "tracked.txt").write_text("changed\n", encoding="utf-8")
            (root / "dev/audit").mkdir(parents=True)
            (root / "dev/audit/tool.py").write_text("x = 1\n", encoding="utf-8")
            (root / "cache.pyc").write_bytes(b"bytecode")
            (root / "dev/.DS_Store").write_bytes(b"finder")
            subprocess.run(["git", "-C", str(root), "add", "-f", "dev/.DS_Store"], check=True)
            config = {
                "artifacts": {"bytecode_globs": ["**/*.pyc", "*.pyc"], "other_generated_globs": ["**/.DS_Store"]},
                "audit_infrastructure": {"paths": ["dev/audit/**"], "exclude": ["artifacts/dev/audit/**"]},
            }
            baseline = {"paths": [{"path": "tracked.txt"}]}
            rules: list[dict[str, object]] = []
            ledger, artifacts, _ = build_inventory(root, config, baseline, rules)
            by_path = {row["path"]: row for row in ledger}
            self.assertEqual(by_path["tracked.txt"]["scope_cohort"], "baseline_cleanup")
            self.assertEqual(by_path["dev/audit/tool.py"]["scope_cohort"], "audit_infrastructure")
            self.assertEqual({row["artifact_type"] for row in artifacts}, {"generated_bytecode", "other_generated"})

    def test_inventory_does_not_execute_checkers_on_deleted_later_paths(self) -> None:
        """Retain deleted later paths without treating them as executable source."""

        with tempfile.TemporaryDirectory() as tmp:
            root = self.make_repo(Path(tmp) / "public")
            legacy = root / "legacy.sh"
            legacy.write_text("#!/usr/bin/env bash\n", encoding="utf-8")
            subprocess.run(["git", "-C", str(root), "add", "legacy.sh"], check=True)
            subprocess.run(["git", "-C", str(root), "commit", "-qm", "legacy"], check=True)
            legacy.unlink()
            config = {
                "artifacts": {"bytecode_globs": [], "other_generated_globs": []},
                "audit_infrastructure": {"paths": ["dev/audit/**"], "exclude": []},
            }
            baseline = {"paths": [{"path": "tracked.txt"}]}
            rules = [{
                "rule_id": "SHELL.SYNTAX",
                "applicable_paths": ["**/*.sh"],
                "semantic_review": False,
            }]
            ledger, _, _ = build_inventory(root, config, baseline, rules)
            by_path = {row["path"]: row for row in ledger}
            self.assertEqual(by_path["legacy.sh"]["applicable_rules"], [])
            self.assertEqual(by_path["legacy.sh"]["content_source"], "index_deleted")
            self.assertNotIn("legacy.sh", executable_ledger_paths(ledger))

    def test_manifest_rejects_unproven_exact_adapter(self) -> None:
        """Require declared parity coverage before claiming an exact adapter."""

        rule = self.manifest_rule()
        rule.update({"rule_id": "EXACT", "execution_kind": "adapter", "coverage_relation": "exact"})
        with self.assertRaises(ValueError):
            validate_manifest({"schema_version": 2, "rule": [rule]})

    def test_manifest_rejects_strict_safe_output_paths(self) -> None:
        """Forbid strict-safe rules from silently declaring repository output."""

        rule = self.manifest_rule()
        rule["output_paths"] = ["artifacts/tests/**"]
        with self.assertRaisesRegex(ValueError, "must not declare"):
            validate_manifest({"schema_version": 2, "rule": [rule]})

    def test_recursive_glob_includes_direct_children(self) -> None:
        """Treat manifest **/ globs as zero-or-more directory components."""

        self.assertTrue(path_matches("bin/entry.sh", ["bin/**/*.sh"]))
        self.assertTrue(path_matches("dev/schemas/finding.schema.json", ["dev/schemas/**"]))
        self.assertTrue(
            path_matches(
                "tests/unit/dev_audit/test_mvp.py",
                ["tests/unit/dev_audit/**"],
            )
        )

    def test_canonical_report_integrity_ignores_only_ds_store(self) -> None:
        """Ignore Finder metadata while retaining every report artifact."""

        with tempfile.TemporaryDirectory() as tmp:
            report = Path(tmp) / "report"
            nested = report / "semantic_review"
            nested.mkdir(parents=True)
            artifact = report / "run.json"
            artifact.write_text('{"status":"completed"}\n', encoding="utf-8")
            (nested / "review.md").write_text("canonical review\n", encoding="utf-8")
            expected = canonical_report_manifest(report)

            finder_metadata = report / ".DS_Store"
            finder_metadata.write_bytes(b"first Finder state")
            (nested / ".DS_Store").write_bytes(b"nested Finder state")
            validate_canonical_report_manifest(report, expected)
            finder_metadata.write_bytes(b"changed Finder state")
            validate_canonical_report_manifest(report, expected)
            finder_metadata.unlink()
            (nested / ".DS_Store").unlink()
            validate_canonical_report_manifest(report, expected)

            artifact.write_text('{"status":"changed"}\n', encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "canonical report integrity mismatch"):
                validate_canonical_report_manifest(report, expected)

    def test_verify_reports_stale_fingerprints(self) -> None:
        """Mark an existing ledger and prompt stale after its target changes."""

        with tempfile.TemporaryDirectory() as tmp:
            public = self.make_repo(Path(tmp) / "public")
            private = self.make_repo(Path(tmp) / "private")
            (public / "tracked.txt").write_text("dirty\n", encoding="utf-8")
            fingerprint, _ = entry_fingerprint(public, "tracked.txt")
            (public / "baseline.json").write_text("{}\n", encoding="utf-8")
            (public / "rules.toml").write_text("schema_version = 1\n", encoding="utf-8")
            baseline_fingerprint, _ = entry_fingerprint(public, "baseline.json")
            rules_fingerprint, _ = entry_fingerprint(public, "rules.toml")
            report = Path(tmp) / "report"
            report.mkdir()
            row = {"path": "tracked.txt", "content_fingerprint": fingerprint}
            (report / "ledger.ndjson").write_text(json.dumps(row) + "\n", encoding="utf-8")
            (report / "findings.ndjson").write_text(json.dumps(row) + "\n", encoding="utf-8")
            (report / "prompts.ndjson").write_text(json.dumps({"batch_id": "one", "target_fingerprints": {"tracked.txt": fingerprint}}) + "\n", encoding="utf-8")
            (report / "artifacts.ndjson").write_text("", encoding="utf-8")
            (report / "checks.ndjson").write_text("", encoding="utf-8")
            (report / "cohort_progress.json").write_text("{}\n", encoding="utf-8")
            (report / "rule_manifest.json").write_text("{}\n", encoding="utf-8")
            run = {
                "status": "completed",
                "baseline_path": "baseline.json",
                "baseline_fingerprint": baseline_fingerprint,
                "rule_manifest_path": "rules.toml",
                "rule_manifest_fingerprint": rules_fingerprint,
                "consumed_evidence": [],
            }
            (report / "run.json").write_text(json.dumps(run) + "\n", encoding="utf-8")
            args = ["--public-root", str(public), "--private-root", str(private), "--verify", str(report)]
            self.assertEqual(main(args), 0)
            (public / "tracked.txt").write_text("changed again\n", encoding="utf-8")
            self.assertEqual(main(args), 3)
            with (report / "staleness.ndjson").open(encoding="utf-8") as handle:
                self.assertEqual(sum(1 for _ in handle), 3)

    def test_verify_requires_completed_unless_partial_is_allowed(self) -> None:
        """Do not describe an aborted audit as successfully verified by default."""

        with tempfile.TemporaryDirectory() as tmp:
            public = self.make_repo(Path(tmp) / "public")
            private = self.make_repo(Path(tmp) / "private")
            (public / "baseline.json").write_text("{}\n", encoding="utf-8")
            (public / "rules.toml").write_text("schema_version = 1\n", encoding="utf-8")
            base_fp, _ = entry_fingerprint(public, "baseline.json")
            rules_fp, _ = entry_fingerprint(public, "rules.toml")
            report = Path(tmp) / "report"
            report.mkdir()
            for name in ("ledger.ndjson", "artifacts.ndjson", "findings.ndjson", "checks.ndjson", "prompts.ndjson"):
                (report / name).write_text("", encoding="utf-8")
            (report / "cohort_progress.json").write_text("{}\n", encoding="utf-8")
            (report / "rule_manifest.json").write_text("{}\n", encoding="utf-8")
            (report / "run.json").write_text(json.dumps({"status": "aborted", "baseline_path": "baseline.json", "baseline_fingerprint": base_fp, "rule_manifest_path": "rules.toml", "rule_manifest_fingerprint": rules_fp, "consumed_evidence": []}), encoding="utf-8")
            base_args = ["--public-root", str(public), "--private-root", str(private), "--verify", str(report)]
            self.assertEqual(main(base_args), 4)
            self.assertEqual(main([*base_args, "--allow-partial"]), 0)

    def test_verify_marks_missing_artifacts_incomplete(self) -> None:
        """Distinguish an incomplete bundle from malformed report content."""

        with tempfile.TemporaryDirectory() as tmp:
            public = self.make_repo(Path(tmp) / "public")
            private = self.make_repo(Path(tmp) / "private")
            report = Path(tmp) / "report"
            report.mkdir()
            self.assertEqual(main(["--public-root", str(public), "--private-root", str(private), "--verify", str(report)]), 4)
            summary = json.loads((report / "verification.json").read_text(encoding="utf-8"))
            self.assertEqual(summary["status"], "error_or_incomplete_report")

    def test_baseline_validation_rejects_requested_invariants(self) -> None:
        """Report duplicate, order, path, count, and malformed cohort defects."""

        cases: list[tuple[str, callable]] = [
            ("duplicate", lambda value: value["paths"].__setitem__(1, dict(value["paths"][0]))),
            ("unsorted", lambda value: value["paths"].reverse()),
            ("bytecode", lambda value: value["paths"][0].update({"path": "cache.pyc"})),
            ("absolute", lambda value: value["paths"][0].update({"path": "/absolute.py"})),
            ("parent", lambda value: value["paths"][0].update({"path": "src/../bad.py"})),
            ("wrong_schema", lambda value: value.update({"schema_version": 2})),
            ("wrong_cohort", lambda value: value.update({"cohort": "other"})),
            ("wrong_count", lambda value: value.update({"authored_path_count": 3})),
            ("malformed", lambda value: value["paths"].__setitem__(0, "not-an-entry")),
        ]
        for name, mutate in cases:
            with self.subTest(name=name):
                value = self.baseline()
                mutate(value)
                with self.assertRaises(BaselineValidationError):
                    validate_baseline(value, self.artifact_config())

    def test_timeout_is_a_controlled_command_result(self) -> None:
        """Return a structured timeout rather than an uncaught subprocess error."""

        with tempfile.TemporaryDirectory() as tmp:
            root = self.make_repo(Path(tmp) / "public")
            result = run_command([sys.executable, "-c", "import time; time.sleep(1)"], root, timeout_seconds=0.01)
            self.assertTrue(result.timed_out)
            self.assertIsNone(result.returncode)

    def test_shell_syntax_batch_returns_per_path_subresults(self) -> None:
        """Represent the atomic adapter's individual path outcomes explicitly."""

        with tempfile.TemporaryDirectory() as tmp:
            root = self.make_repo(Path(tmp) / "public")
            (root / "good.sh").write_text("#!/usr/bin/env bash\necho good\n", encoding="utf-8")
            (root / "bad.sh").write_text("if then\n", encoding="utf-8")
            result, subresults = run_shell_syntax_batch(["good.sh", "bad.sh"], root, 10, 4096)
            self.assertEqual(result.returncode, 1)
            self.assertEqual([row["status"] for row in subresults], ["completed", "finding"])

    def test_authoritative_bash_satisfies_minimum_and_is_reported(self) -> None:
        """Resolve Bash from env_protocol rather than the invoking shell."""

        resolved = resolve_env_bash()
        self.assertTrue(resolved.path.is_absolute())
        self.assertGreaterEqual(resolved.version_tuple, MINIMUM_BASH)
        self.assertTrue(resolved.satisfies_minimum)
        self.assertIn(str(resolved.path), resolved.label)
        self.assertIn(resolved.version, resolved.label)

    def test_syntax_routing_uses_env_bash_except_for_posix_bootstrap(self) -> None:
        """Keep the install entrypoint under POSIX sh while Bash owns .sh gates."""

        resolved = resolve_env_bash()
        bash_command = syntax_command("bin/submit_compute_signal.sh", resolved)
        posix_command = syntax_command(POSIX_BOOTSTRAP, resolved)
        self.assertEqual(bash_command[:2], [str(resolved.path), "-n"])
        self.assertEqual(posix_command[:2], ["sh", "-n"])


if __name__ == "__main__":
    unittest.main()
