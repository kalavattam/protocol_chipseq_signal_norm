#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_coordinator.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Focused tests for the portable Slurm source-and-results coordinator.
"""

from __future__ import annotations

import argparse
import contextlib
import importlib.util
import io
import os
import subprocess
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest import mock

COORDINATOR = (
    Path(__file__).resolve().parents[2]
    / "integration"
    / "slurm"
    / "coordinator.py"
)
SPEC = importlib.util.spec_from_file_location("slurm_coordinator", COORDINATOR)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot load coordinator: {COORDINATOR}")

bundle = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = bundle
SPEC.loader.exec_module(bundle)


class TemporaryRepository:
    """
    Create a small dirty Git repository with relevant exclusion cases.
    """

    def __init__(self, root: Path) -> None:
        self.root = root

    def run(self, *args: str) -> None:
        """
        Run a required Git command in the temporary repository.
        """

        subprocess.run(
            ["git", *args],
            cwd=self.root,
            check=True,
            capture_output=True,
        )

    def create(self) -> None:
        """
        Populate committed, modified, untracked, and ignored files.
        """

        self.root.mkdir()
        self.run("init", "-q")
        self.run("config", "user.name", "Bundle Test")
        self.run("config", "user.email", "bundle@example.invalid")

        (self.root / ".gitignore").write_text(
            "ignored/\nartifacts/tests/\nartifacts/dev/audit/\n",
            encoding="utf-8",
        )
        (self.root / "bin").mkdir()
        (self.root / "bin" / "needed.sh").write_text(
            "#!/usr/bin/env bash\n# original\n",
            encoding="utf-8",
        )
        self.run("add", ".gitignore", "bin/needed.sh")
        self.run("commit", "-qm", "fixture")

        (self.root / "bin" / "needed.sh").write_text(
            "#!/usr/bin/env bash\n# dirty\n",
            encoding="utf-8",
        )
        required = (
            self.root / "src/protocol_chipseq_signal_norm/required_new.py"
        )
        required.parent.mkdir(parents=True)
        required.write_text("# untracked\n", encoding="utf-8")

        for relative in (
            "ignored/cache.txt",
            "artifacts/tests/generated.txt",
            "artifacts/dev/audit/report.json",
            "pkg/__pycache__/module.pyc",
            ".DS_Store",
        ):
            path = self.root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("excluded\n", encoding="utf-8")


def fake_executable(directory: Path, name: str, body: str) -> None:
    """
    Write one executable fake command.
    """

    path = directory / name
    path.write_text(
        "#!/usr/bin/env bash\nset -euo pipefail\n" + body,
        encoding="utf-8",
    )
    path.chmod(0o700)


def result_job(key: str, state: str = "COMPLETED") -> dict[str, object]:
    """
    Return one complete verifier job record.
    """

    succeeded = state == "COMPLETED"
    return {
        "job_key": key,
        "job_id": "123",
        "job_name": f"wet-{key}",
        "command": ["sbatch", "job.sh"],
        "requested_resources": {"cpus_per_task": 1},
        "submit_timestamp": "2026-07-16T00:00:00Z",
        "start_timestamp": "2026-07-16T00:00:01Z",
        "finish_timestamp": "2026-07-16T00:00:02Z",
        "final_state": state,
        "exit_code": "0:0" if succeeded else "1:0",
        "stdout_path": f"stdout/{key}.out",
        "stderr_path": f"stderr/{key}.err",
        "assertions": [{"kind": "nonempty", "passed": succeeded}],
        "cleanup_result": "passed",
        "success": succeeded,
    }


def write_result_bundle(
    root: Path,
    run_id: str,
    source_digest: str,
    jobs: list[dict[str, object]],
    required: list[str],
) -> None:
    """
    Write a minimal structurally complete result bundle for verification.
    """

    for name in ("stdout", "stderr", "artifacts"):
        (root / name).mkdir(parents=True, exist_ok=True)

    bundle.write_json(
        root / "preflight.json",
        {"run_id": run_id, "success": True},
    )
    bundle.write_json(
        root / "run_manifest.json",
        {
            "run_id": run_id,
            "source_manifest_sha256": source_digest,
            "required_job_keys": required,
        },
    )
    bundle.write_json(root / "jobs.json", {"run_id": run_id, "jobs": jobs})

    success = len(jobs) == len(required) and all(
        job["success"] for job in jobs
    )
    bundle.write_json(
        root / "summary.json",
        {
            "run_id": run_id,
            "success": success,
            "job_count": len(jobs),
            "required_job_count": len(required),
        },
    )

    (root / "commands.log").write_text("sbatch job.sh\n", encoding="utf-8")
    (root / "summary.txt").write_text(
        f"run_id: {run_id}\nsuccess: {'yes' if success else 'no'}\n",
        encoding="utf-8",
    )
    (root / "exit_status.txt").write_text(
        "0\n" if success else "1\n",
        encoding="utf-8",
    )

    bundle.write_result_checksums(root)


class SourceBundleTests(unittest.TestCase):
    """
    Test source inventory, deterministic archives, and path safety.
    """

    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.base = Path(self.temp.name)
        self.repo = TemporaryRepository(self.base / "repo")
        self.repo.create()

    def tearDown(self) -> None:
        self.temp.cleanup()

    def manifest(self) -> dict[str, object]:
        """
        Build the fixed-time test manifest.
        """

        return bundle.build_source_manifest(
            self.repo.root,
            "wet-test-001",
            "test scope",
            "2026-07-16T00:00:00Z",
        )

    def test_manifest_is_deterministic_and_includes_dirty_untracked_files(
        self,
    ) -> None:
        first = self.manifest()
        second = self.manifest()

        self.assertEqual(first, second)

        paths = {entry["path"] for entry in first["files"]}

        self.assertIn("bin/needed.sh", paths)
        self.assertIn(
            "src/protocol_chipseq_signal_norm/required_new.py",
            paths,
        )

        tracked = next(
            entry
            for entry in first["files"]
            if entry["path"] == "bin/needed.sh"
        )

        self.assertEqual(
            tracked["sha256"],
            bundle.sha256_bytes(b"#!/usr/bin/env bash\n# dirty\n"),
        )
        self.assertEqual(tracked["git_state"], "tracked")

        required_path = "src/protocol_chipseq_signal_norm/required_new.py"
        required = next(
            entry for entry in first["files"] if entry["path"] == required_path
        )

        self.assertEqual(required["git_state"], "required_untracked")
        self.assertTrue(first["porcelain_status"])
        self.assertEqual(len(first["public_git_head"]), 40)

    def test_manifest_excludes_ignored_outputs_caches_reports_and_metadata(
        self,
    ) -> None:
        paths = {entry["path"] for entry in self.manifest()["files"]}

        self.assertNotIn("ignored/cache.txt", paths)
        self.assertNotIn("artifacts/tests/generated.txt", paths)
        self.assertNotIn("artifacts/dev/audit/report.json", paths)
        self.assertNotIn("pkg/__pycache__/module.pyc", paths)
        self.assertNotIn(".DS_Store", paths)

    def test_archive_bytes_and_checksums_are_deterministic(self) -> None:
        manifest = self.manifest()
        first = self.base / "first.tar.gz"
        second = self.base / "second.tar.gz"

        digest_1 = bundle.create_source_archive(
            self.repo.root,
            first,
            manifest,
        )
        digest_2 = bundle.create_source_archive(
            self.repo.root,
            second,
            manifest,
        )

        self.assertEqual(digest_1, digest_2)
        self.assertEqual(first.read_bytes(), second.read_bytes())

        extracted = self.base / "extracted"
        source = bundle.safe_extract_archive(first, extracted)

        self.assertIn("# dirty", (source / "bin" / "needed.sh").read_text())

    def test_archive_path_traversal_is_rejected(self) -> None:
        archive = self.base / "unsafe.tar.gz"
        payload = self.base / "payload"
        payload.write_text("bad\n", encoding="utf-8")

        with tarfile.open(archive, "w:gz") as handle:
            handle.add(payload, arcname="../escape")

        with self.assertRaises(bundle.BundleError):
            bundle.safe_extract_archive(archive, self.base / "unsafe")

    def test_run_id_and_remote_root_validation(self) -> None:
        self.assertEqual(bundle.validate_run_id("wet-7i_01"), "wet-7i_01")

        for value in ("../wet", ".wet", "wet/one", "wet..one", ""):
            with (
                self.subTest(value=value),
                self.assertRaises(bundle.BundleError),
            ):
                bundle.validate_run_id(value)

        self.assertEqual(
            bundle.safe_remote_run_dir("cluster/runs", "wet-1").as_posix(),
            "cluster/runs/wet-1",
        )

        for root in (
            "/",
            "../runs",
            "runs/../other",
            "runs with spaces",
            "x;id",
            "",
        ):
            with (
                self.subTest(root=root),
                self.assertRaises(bundle.BundleError),
            ):
                bundle.safe_remote_run_dir(root, "wet-1")


class InterfaceSafetyTests(unittest.TestCase):
    """
    Test dry runs, command rendering, and cleanup confinement.
    """

    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.base = Path(self.temp.name)

    def tearDown(self) -> None:
        self.temp.cleanup()

    def test_prepare_dry_run_writes_nothing(self) -> None:
        repo = TemporaryRepository(self.base / "repo")
        repo.create()
        bundle_dir = self.base / "bundles"
        args = argparse.Namespace(
            root=str(repo.root),
            run_id="dry-1",
            scope="dry scope",
            created_at="2026-07-16T00:00:00Z",
            bundle_dir=str(bundle_dir),
            dry_run=True,
        )

        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(bundle.prepare(args), 0)

        self.assertFalse(bundle_dir.exists())

    def test_push_pull_and_remote_launch_commands_are_isolated(self) -> None:
        session = self.base / "wet-1"
        config = {
            "ssh_host": "rhino",
            "ssh_user": "user",
            "remote_run_dir": "cluster/runs/wet-1",
            "remote_python": "/opt/env/bin/python3",
        }
        commands = bundle.push_commands(session, config)
        rendered = "\n".join(bundle.shlex.join(item) for item in commands)

        self.assertIn("cluster/runs/wet-1/incoming", rendered)
        self.assertNotIn("--delete", rendered)

        pull = bundle.pull_command(session, config, self.base / "results")

        self.assertIn("cluster/runs/wet-1/results/", " ".join(pull))

        launch = bundle.remote_launch_command({**config, "run_id": "wet-1"})

        self.assertTrue(
            launch.startswith(
                "PYTHONDONTWRITEBYTECODE=1 "
                "PYTHONPATH=cluster/runs/wet-1/incoming/launcher/src "
                "RUN_SLURM=1 WAIT_SLURM=1 CONFIRM_SLURM_WET=1 "
                "/opt/env/bin/python3 ",
            ),
        )
        self.assertNotIn("--protect-args", rendered)
        self.assertNotIn("--protect-args", pull)

    def test_pull_fails_closed_when_result_destination_exists(self) -> None:
        """
        Refuse to overlay a prior local result tree.
        """

        bundle_dir = self.base / "bundles"
        session = bundle_dir / "wet-1"
        (session / "incoming").mkdir(parents=True)
        bundle.write_json(
            session / "incoming" / "remote_config.json",
            {
                "run_id": "wet-1",
                "ssh_host": "rhino",
                "ssh_user": "user",
                "remote_run_dir": "cluster/runs/wet-1",
                "result_destination": str(self.base / "results"),
            },
        )

        destination = self.base / "results" / "wet-1"
        destination.mkdir(parents=True)

        args = argparse.Namespace(
            bundle_dir=str(bundle_dir),
            run_id="wet-1",
            result_destination=None,
            dry_run=False,
        )

        with self.assertRaisesRegex(bundle.BundleError, "already exists"):
            bundle.pull(args)

    def test_clean_requires_marker_and_remains_confined(self) -> None:
        base = self.base / "bundles"
        session = base / "wet-1"
        session.mkdir(parents=True)

        outside = self.base / "outside.txt"
        outside.write_text("keep\n", encoding="utf-8")

        args = argparse.Namespace(
            bundle_dir=str(base),
            run_id="wet-1",
            dry_run=False,
        )

        with self.assertRaises(bundle.BundleError):
            bundle.clean(args)

        bundle.write_json(
            session / "local_state.json",
            {
                "run_id": "wet-1",
                "session_created_by": "slurm_bundle.py/1.0.0",
            },
        )

        self.assertEqual(bundle.clean(args), 0)
        self.assertFalse(session.exists())
        self.assertTrue(outside.exists())

        args.run_id = "../outside.txt"

        with self.assertRaises(bundle.BundleError):
            bundle.clean(args)


class ResultVerifierTests(unittest.TestCase):
    """
    Test result schema, checksums, job completeness, and malformed input.
    """

    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.base = Path(self.temp.name)
        self.run_id = "wet-results-1"
        self.digest = "a" * 64

    def tearDown(self) -> None:
        self.temp.cleanup()

    def valid(self) -> Path:
        """
        Create and return one valid two-job result tree.
        """

        result = self.base / "valid"
        required = ["align_bowtie2_se", "compute_signal_bam_se"]
        jobs = [result_job(key) for key in required]
        write_result_bundle(result, self.run_id, self.digest, jobs, required)

        return result

    def test_complete_results_and_checksums_verify(self) -> None:
        report = bundle.verify_results(self.valid(), self.run_id, self.digest)

        self.assertTrue(report["success"])
        self.assertEqual(report["jobs"], 2)

    def test_checksum_damage_is_rejected(self) -> None:
        result = self.valid()

        (result / "summary.txt").write_text("damaged\n", encoding="utf-8")

        with self.assertRaisesRegex(bundle.BundleError, "checksum mismatch"):
            bundle.verify_results(result, self.run_id, self.digest)

    def test_missing_and_failed_jobs_are_rejected(self) -> None:
        required = ["align_bowtie2_se", "compute_signal_bam_se"]

        missing = self.base / "missing"
        write_result_bundle(
            missing,
            self.run_id,
            self.digest,
            [result_job("align_bowtie2_se")],
            required,
        )

        with self.assertRaisesRegex(bundle.BundleError, "silently omitted"):
            bundle.verify_results(missing, self.run_id, self.digest)

        failed = self.base / "failed"
        jobs = [result_job(required[0]), result_job(required[1], "FAILED")]
        write_result_bundle(failed, self.run_id, self.digest, jobs, required)

        with self.assertRaises(bundle.BundleError):
            bundle.verify_results(failed, self.run_id, self.digest)

    def test_malformed_result_json_and_run_id_mismatch_are_rejected(
        self,
    ) -> None:
        malformed = self.valid()
        (malformed / "jobs.json").write_text("{bad\n", encoding="utf-8")
        bundle.write_result_checksums(malformed)

        with self.assertRaisesRegex(bundle.BundleError, "malformed JSON"):
            bundle.verify_results(malformed, self.run_id, self.digest)

        mismatch = self.valid()

        with self.assertRaisesRegex(bundle.BundleError, "run ID mismatch"):
            bundle.verify_results(mismatch, "another-run", self.digest)


class FakeSlurmTests(unittest.TestCase):
    """
    Exercise gate, preflight, accounting parsing, and failure handling.
    """

    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.base = Path(self.temp.name)
        self.bin = self.base / "bin"
        self.bin.mkdir()

        fake_executable(
            self.bin,
            "sbatch",
            'if [[ -n "${FAKE_SBATCH_TOUCHED:-}" ]]; then\n'
            '  : > "${FAKE_SBATCH_TOUCHED}"\nfi\n'
            'if [[ "${1:-}" == "--test-only" ]]; then exit 0; fi\n'
            'counter="${FAKE_SLURM_COUNTER}"\n'
            'value=100\n[[ -f "${counter}" ]] && value="$(< "${counter}")"\n'
            'value=$(( value + 1 ))\nprintf "%s\\n" "${value}" > '
            '"${counter}"\n'
            'printf "%s\\n" "${value}"\n',
        )

        fake_executable(self.bin, "squeue", "exit 0\n")
        fake_executable(
            self.bin,
            "sacct",
            'job=""\nwhile [[ "$#" -gt 0 ]]; do\n'
            '  if [[ "$1" == "-j" ]]; then job="$2"; shift 2; else shift; fi\n'
            'done\nstate="${FAKE_SLURM_STATE:-COMPLETED}"\n'
            'code="0:0"\n[[ "${state}" == "COMPLETED" ]] || code="1:0"\n'
            'printf "%s|%s|%s|submit|start|end\\n" "${job}" '
            '"${state}" "${code}"\n',
        )
        fake_executable(self.bin, "scontrol", "exit 0\n")

        fake_executable(
            self.bin,
            "conda",
            'if [[ "${1:-}" == "env" && "${2:-}" == "list" ]]; then\n'
            '  printf "env_protocol /fake/env_protocol\\n"\n  exit 0\n'
            "fi\nexit 1\n",
        )

    def tearDown(self) -> None:
        self.temp.cleanup()

    def runtime_config(self, name: str) -> tuple[dict[str, object], Path]:
        """
        Create a minimal verified runtime tree for fake scheduler tests.
        """

        session = self.base / name
        source = session / "source_bundle" / "source"
        incoming = session / "incoming"
        result = session / "results"

        source.mkdir(parents=True)
        incoming.mkdir()
        manifest = {"files": []}
        bundle.write_json(source.parent / "manifest.json", manifest)

        archive = incoming / "source.tar.gz"
        archive.write_bytes(b"archive\n")

        config: dict[str, object] = {
            "run_id": name,
            "session_dir": str(session),
            "source_dir": str(source),
            "result_dir": str(result),
            "source_archive": archive.name,
            "source_archive_sha256": bundle.sha256_file(archive),
            "source_manifest_sha256": bundle.sha256_file(
                source.parent / "manifest.json",
            ),
            "partition": "test-partition",
            "account": "test-account",
            "environment_name": "env_protocol",
            "declared_validation_scope": "fake wet scope",
            "resources": {
                "job_count": 2,
                "cpus_per_task": 1,
                "memory": "256M",
                "time_limit": "00:01:00",
                "poll_seconds": 1,
                "wait_seconds": 60,
            },
        }

        config_path = session / "runtime_config.json"
        bundle.write_json(config_path, config)

        return config, config_path

    def definitions(
        self,
        config: dict[str, object],
    ) -> list[dict[str, object]]:
        """
        Return two fake-worker definitions with stable passing assertions.
        """

        marker = Path(str(config["session_dir"])) / "runtime_config.json"
        return [
            {
                "job_key": key,
                "job_name": f"wet-{key}",
                "worker": Path("/bin/true"),
                "arguments": [],
                "assertions": [{"kind": "nonempty", "path": str(marker)}],
                "output_directories": [
                    Path(str(config["result_dir"])) / "artifacts" / key,
                    Path(str(config["result_dir"]))
                    / "artifacts"
                    / key
                    / "logs",
                ],
            }
            for key in ("job_one", "job_two")
        ]

    def wet_environment(self, state: str = "COMPLETED") -> dict[str, str]:
        """
        Return the exact wet gates and fake command PATH.
        """

        return {
            **os.environ,
            "PATH": f"{self.bin}:{os.environ['PATH']}",
            "RUN_SLURM": "1",
            "WAIT_SLURM": "1",
            "CONFIRM_SLURM_WET": "1",
            "FAKE_SLURM_COUNTER": str(self.base / "counter"),
            "FAKE_SLURM_STATE": state,
        }

    def test_real_job_definitions_declare_each_worker_output_directory(
        self,
    ) -> None:
        config = {
            "run_id": "directory-contract",
            "source_dir": str(self.base / "source"),
            "result_dir": str(self.base / "results"),
            "environment_name": "env_protocol",
        }
        definitions = bundle.job_definitions(config)

        self.assertEqual(len(definitions), 2)

        for definition in definitions:
            arguments = definition["arguments"]
            configured = (
                {
                    Path(arguments[arguments.index("--dir_out") + 1]),
                    Path(arguments[arguments.index("--dir_eo") + 1]),
                }
                if "--dir_out" in arguments
                else {
                    Path(
                        arguments[arguments.index("--csv_fil_out") + 1],
                    ).parent,
                    Path(arguments[arguments.index("--dir_eo") + 1]),
                }
            )

            self.assertEqual(set(definition["output_directories"]), configured)

    def test_worker_output_directories_exist_before_sbatch_test_only(
        self,
    ) -> None:
        config, _ = self.runtime_config("fake-output-directories")
        definitions = self.definitions(config)
        observed: list[bool] = []

        def inspect(command: list[str], **_: object) -> object:
            if command[0] == "sbatch" and "--test-only" in command:
                observed.append(
                    all(
                        Path(directory).is_dir()
                        for definition in definitions
                        for directory in definition["output_directories"]
                    ),
                )

            return bundle.CommandResult(tuple(command), 0, "", "")

        self.assertEqual(
            bundle.prepare_job_output_directories(definitions),
            (True, "all worker output directories prepared"),
        )

        inspect(["sbatch", "--test-only", "job.sh"])

        self.assertEqual(observed, [True])

    def test_fake_slurm_success_records_both_terminal_jobs(self) -> None:
        config, config_path = self.runtime_config("fake-success")
        args = argparse.Namespace(config=str(config_path))

        with (
            mock.patch.dict(os.environ, self.wet_environment(), clear=True),
            mock.patch.object(
                bundle,
                "job_definitions",
                return_value=self.definitions(config),
            ),
            mock.patch.object(
                bundle,
                "prepare_remote_fixtures",
                return_value=(True, "fake fixtures ready"),
            ),
        ):
            self.assertEqual(bundle.wet_run(args), 0)

        jobs = bundle.load_json(Path(str(config["result_dir"])) / "jobs.json")[
            "jobs"
        ]

        self.assertEqual(
            [job["final_state"] for job in jobs],
            ["COMPLETED"] * 2,
        )
        self.assertTrue(all(job["exit_code"] == "0:0" for job in jobs))

        for definition in self.definitions(config):
            for directory in definition["output_directories"]:
                self.assertTrue(directory.is_dir())

        commands = (
            (Path(str(config["result_dir"])) / "commands.log")
            .read_text(encoding="utf-8")
            .splitlines()
        )
        test_only = [
            index for index, row in enumerate(commands) if "--test-only" in row
        ]
        submissions = [
            index
            for index, row in enumerate(commands)
            if row.startswith("sbatch ") and "--test-only" not in row
        ]

        self.assertEqual(len(test_only), 2)
        self.assertEqual(len(submissions), 2)
        self.assertTrue(
            all(
                test < submit
                for test, submit in zip(test_only, submissions, strict=True)
            ),
        )

    def test_fake_failed_job_is_not_submission_success(self) -> None:
        config, config_path = self.runtime_config("fake-failed")
        args = argparse.Namespace(config=str(config_path))

        with (
            mock.patch.dict(
                os.environ,
                self.wet_environment("FAILED"),
                clear=True,
            ),
            mock.patch.object(
                bundle,
                "job_definitions",
                return_value=self.definitions(config),
            ),
            mock.patch.object(
                bundle,
                "prepare_remote_fixtures",
                return_value=(True, "fake fixtures ready"),
            ),
        ):
            self.assertEqual(bundle.wet_run(args), 1)

        summary = bundle.load_json(
            Path(str(config["result_dir"])) / "summary.json",
        )

        self.assertFalse(summary["success"])
        self.assertEqual(summary["passed_job_count"], 0)

    def test_missing_gate_prevents_fake_sbatch_call(self) -> None:
        _, config_path = self.runtime_config("fake-gate")
        args = argparse.Namespace(config=str(config_path))
        environment = self.wet_environment()
        environment.pop("CONFIRM_SLURM_WET")

        with (
            mock.patch.dict(os.environ, environment, clear=True),
            self.assertRaisesRegex(bundle.BundleError, "CONFIRM_SLURM_WET=1"),
        ):
            bundle.wet_run(args)

        self.assertFalse((self.base / "counter").exists())

    def test_retired_slurm_wait_spelling_cannot_enable_submission(
        self,
    ) -> None:
        _, config_path = self.runtime_config("fake-retired-gate")
        args = argparse.Namespace(config=str(config_path))

        environment = self.wet_environment()
        environment.pop("WAIT_SLURM")
        retired_name = "SLURM" + "_WAIT"
        environment[retired_name] = "1"

        with (
            mock.patch.dict(os.environ, environment, clear=True),
            self.assertRaisesRegex(bundle.BundleError, "WAIT_SLURM=1"),
        ):
            bundle.wet_run(args)

        self.assertFalse((self.base / "counter").exists())

    def test_shell_runner_missing_confirmation_never_calls_fake_sbatch(
        self,
    ) -> None:
        config = self.base / "shell_config.json"
        config.write_text("{}\n", encoding="utf-8")
        touched = self.base / "sbatch_touched"

        environment = {
            **os.environ,
            "PATH": f"{self.bin}:{os.environ['PATH']}",
            "RUN_SLURM": "1",
            "WAIT_SLURM": "1",
            "FAKE_SBATCH_TOUCHED": str(touched),
        }
        runner = COORDINATOR.parent / "run_wet_tests.sh"

        completed = subprocess.run(
            ["bash", str(runner), "--config", str(config)],
            env=environment,
            text=True,
            capture_output=True,
            check=False,
        )

        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("CONFIRM_SLURM_WET=1", completed.stderr)
        self.assertFalse(touched.exists())

    def test_sacct_parser_ignores_steps_and_parses_terminal_row(self) -> None:
        text = "123.batch|COMPLETED|0:0|s|a|e\n123|FAILED|1:0|s|a|e\n"

        self.assertEqual(bundle.parse_sacct(text, "123")["state"], "FAILED")


if __name__ == "__main__":
    unittest.main()
