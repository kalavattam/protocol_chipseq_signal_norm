#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: run.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Run the strict, read-only cleanup audit MVP.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import tomllib
import traceback
from collections import Counter
from dataclasses import dataclass
from dataclasses import field as dataclass_field
from pathlib import Path, PurePosixPath
from typing import Any

from dev.audit.generate_prompts import (
    REQUIRED_BEHAVIORAL_ANCHORS,
    SOURCE_ONLY_ANCHORS,
    build_batches,
    configured_semantic_units,
    write_configured_bundle,
    write_pilot_bundle,
    write_prompt_bundle,
)
from dev.audit.parse_findings import (
    CommandResult,
    parse_pilot_payload,
    parse_result,
)
from dev.audit.shell_validation import resolve_env_bash, syntax_command

BASELINE_SCHEMA_VERSION = 1
BASELINE_COHORT = "baseline_cleanup"
BASELINE_COUNT = 184
BASELINE_STATUSES = {" M", "??"}
EXIT_GUARD_ABORTED = 2
EXIT_VERIFY_STALE = 3
EXIT_VERIFY_NONCOMPLETED = 4
EXIT_VERIFY_MALFORMED = 5
EXIT_RUNTIME_ERROR = 6
PILOT_ARTIFACT_SCHEMA_VERSION = 1
CANONICAL_REPORT_METADATA_EXCLUSIONS = frozenset({".DS_Store"})
PILOT_ARTIFACT_REQUIRED_FIELDS = {
    "targets.ndjson": {
        "schema_version",
        "run_id",
        "path",
        "role",
        "content_fingerprint",
        "git_state_labels",
        "whitespace_coverage",
    },
    "facts.ndjson": {
        "schema_version",
        "run_id",
        "rule_id",
        "topic",
        "path",
        "content_fingerprint",
        "certainty",
        "value",
    },
    "policy_questions.ndjson": {
        "schema_version",
        "run_id",
        "rule_id",
        "topic",
        "question",
        "paths",
    },
    "adapter_limitations.ndjson": {
        "schema_version",
        "run_id",
        "rule_id",
        "adapter",
        "limitation",
    },
    "semantic_review.ndjson": {
        "schema_version",
        "run_id",
        "bundle_id",
        "path",
        "target_fingerprints",
    },
}
CONTROLLED_SMOKE_KIND = "controlled_smoke_summary_v1"
CONTROLLED_SMOKE_SUMMARY = re.compile(
    r"^Summary for (?P<test>.+): pass=(?P<pass>\d+) fail=(?P<fail>\d+) "
    r"warn=(?P<warn>\d+) skip=(?P<skip>\d+)$",
    re.MULTILINE,
)


class BaselineValidationError(ValueError):
    """
    Raised when the immutable cleanup cohort violates an invariant.
    """


class IncompleteReportError(ValueError):
    """
    Raised when a report bundle lacks a required artifact.
    """


def sha256_bytes(value: bytes) -> str:
    """
    Return a prefixed SHA-256 digest.
    """

    return "sha256:" + hashlib.sha256(value).hexdigest()


def json_bytes(value: Any) -> bytes:
    """
    Serialize one value deterministically for hashing.
    """

    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    ).encode("utf-8")


def file_sha256(path: Path) -> str:
    """
    Return a SHA-256 fingerprint of one generated artifact file.
    """

    return sha256_bytes(path.read_bytes())


def canonical_report_manifest(report: Path) -> tuple[tuple[str, str], ...]:
    """
    Hash canonical report files in relative-path order, excluding OS metadata.
    """

    if not report.is_dir():
        raise ValueError(f"report directory does not exist: {report}")

    paths = sorted(
        (
            path
            for path in report.rglob("*")
            if path.is_file()
            and path.name not in CANONICAL_REPORT_METADATA_EXCLUSIONS
        ),
        key=lambda path: path.relative_to(report).as_posix(),
    )

    return tuple(
        (path.relative_to(report).as_posix(), file_sha256(path))
        for path in paths
    )


def validate_canonical_report_manifest(
    report: Path,
    expected: tuple[tuple[str, str], ...],
) -> None:
    """
    Require exact canonical report path and byte-hash equality.
    """

    if canonical_report_manifest(report) != expected:
        raise ValueError("canonical report integrity mismatch")


def pilot_report_self_hash(report: dict[str, object]) -> str:
    """
    Hash canonical report JSON while omitting its documented self field.
    """

    provenance = dict(report.get("package_provenance", {}))
    self_hash = dict(provenance.get("pilot_report_self_hash", {}))
    self_hash.pop("value", None)
    provenance["pilot_report_self_hash"] = self_hash
    value = dict(report)
    value["package_provenance"] = provenance

    return sha256_bytes(json_bytes(value))


def bounded_text(value: str | bytes | None, maximum: int) -> tuple[str, bool]:
    """
    Decode and bound captured subprocess output.
    """

    if value is None:
        return "", False

    text = (
        value.decode("utf-8", "replace") if isinstance(value, bytes) else value
    )
    if len(text.encode("utf-8")) <= maximum:
        return text, False

    return text.encode("utf-8")[:maximum].decode(
        "utf-8",
        "ignore",
    ) + "\n[truncated]", True


def run_command(
    argv: list[str],
    cwd: Path,
    env: dict[str, str] | None = None,
    timeout_seconds: float = 120.0,
    max_output_bytes: int = 65536,
) -> CommandResult:
    """
    Run argv without a shell, timeout safely, and capture bounded output.
    """

    started = time.monotonic()

    try:
        proc = subprocess.run(
            argv,
            cwd=cwd,
            env=env,
            text=True,
            capture_output=True,
            check=False,
            timeout=timeout_seconds,
        )
        stdout, _ = bounded_text(proc.stdout, max_output_bytes)
        stderr, _ = bounded_text(proc.stderr, max_output_bytes)

        return CommandResult(
            argv,
            stdout,
            stderr,
            proc.returncode,
            round(time.monotonic() - started, 6),
        )
    except subprocess.TimeoutExpired as exc:
        stdout, _ = bounded_text(exc.stdout, max_output_bytes)
        stderr, _ = bounded_text(exc.stderr, max_output_bytes)

        return CommandResult(
            argv,
            stdout,
            stderr,
            None,
            round(time.monotonic() - started, 6),
            True,
        )
    except OSError as exc:
        return CommandResult(
            argv,
            "",
            "",
            None,
            round(time.monotonic() - started, 6),
            False,
            f"{type(exc).__name__}: {exc}",
        )


def run_shell_syntax_batch(
    paths: list[str],
    root: Path,
    timeout_seconds: float,
    max_output_bytes: int,
) -> tuple[CommandResult, list[dict[str, object]]]:
    """
    Run one atomic shell-syntax adapter with structured path subresults.
    """

    started = time.monotonic()
    interpreter = resolve_env_bash()
    if not interpreter.satisfies_minimum:
        raise RuntimeError(interpreter.label)

    subresults: list[dict[str, object]] = []
    stdout_parts: list[str] = []
    stderr_parts: list[str] = []

    for index, path in enumerate(paths):
        remaining = timeout_seconds - (time.monotonic() - started)

        if remaining <= 0:
            subresults.extend(
                {"path": pending, "status": "not_run_due_timeout"}
                for pending in paths[index:]
            )
            aggregate = CommandResult(
                [
                    "audit-shell-syntax-batch",
                    str(interpreter.path),
                    "-n",
                    *paths,
                ],
                "\n".join(stdout_parts),
                "\n".join(stderr_parts),
                None,
                round(time.monotonic() - started, 6),
                True,
            )

            return aggregate, subresults

        command = syntax_command(path, interpreter)
        command[-1] = str(root / path)

        result = run_command(
            command,
            root,
            dict(os.environ, PYTHONDONTWRITEBYTECODE="1"),
            remaining,
            max_output_bytes,
        )
        subresults.append(
            {
                "path": path,
                "interpreter_path": command[0],
                "interpreter_version": interpreter.version
                if path != "install/scripts/install_envs_entrypoint.sh"
                else "POSIX sh",
                "status": "completed"
                if result.returncode == 0
                else "finding"
                if result.returncode
                else "infrastructure_error",
                "exit_status": result.returncode,
                "stdout": result.stdout,
                "stderr": result.stderr,
                "timed_out": result.timed_out,
                "launch_error": result.launch_error,
            },
        )
        stdout_parts.append(result.stdout)
        stderr_parts.append(result.stderr)

        if result.timed_out or result.launch_error:
            subresults.extend(
                {
                    "path": pending,
                    "status": "not_run_due_timeout"
                    if result.timed_out
                    else "not_run_due_infrastructure_error",
                }
                for pending in paths[index + 1 :]
            )
            aggregate = CommandResult(
                [
                    "audit-shell-syntax-batch",
                    str(interpreter.path),
                    "-n",
                    *paths,
                ],
                interpreter.label + "\n" + "\n".join(stdout_parts),
                "\n".join(stderr_parts),
                None,
                round(time.monotonic() - started, 6),
                result.timed_out,
                result.launch_error,
            )

            return aggregate, subresults

    failed = any(item["status"] == "finding" for item in subresults)

    return CommandResult(
        ["audit-shell-syntax-batch", str(interpreter.path), "-n", *paths],
        interpreter.label + "\n" + "\n".join(stdout_parts),
        "\n".join(stderr_parts),
        1 if failed else 0,
        round(time.monotonic() - started, 6),
    ), subresults


def git(
    root: Path,
    *args: str,
    text: bool = True,
) -> subprocess.CompletedProcess[Any]:
    """
    Run a read-only Git query rooted at root.
    """

    return subprocess.run(
        ["git", *args],
        cwd=root,
        text=text,
        capture_output=True,
        check=False,
    )


def resolve_repo(path: Path, label: str) -> Path:
    """
    Resolve and validate a repository root.
    """

    root = path.resolve()
    proc = git(root, "rev-parse", "--show-toplevel")
    if proc.returncode or Path(proc.stdout.strip()).resolve() != root:
        raise ValueError(f"{label} is not a Git repository root: {root}")

    return root


def resolve_reports_base(root: Path, reports_dir: Path) -> Path:
    """
    Resolve a bounded report output directory.

    Formal reports belong under 'artifacts/dev/audit'; dry renders may use
    '/private/tmp'.

    Parameters
    ----------
    root : Path
        Repository root that owns the formal report directory.
    reports_dir : Path
        Requested absolute or repository-relative report directory.

    Returns
    -------
    path : Path
        Resolved, permitted report directory.

    Raises
    ------
    ValueError
        If the requested directory is outside the permitted report roots.
    """

    resolved = (
        (root / reports_dir).resolve()
        if not reports_dir.is_absolute()
        else reports_dir.resolve()
    )
    private_temporary_root = Path("/private/tmp").resolve()

    if (
        resolved == root / "artifacts/dev/audit"
        or resolved == private_temporary_root
        or private_temporary_root in resolved.parents
    ):
        return resolved

    raise ValueError(
        "reports_dir must be artifacts/dev/audit or /private/tmp subtree",
    )


def parse_status(root: Path) -> tuple[bytes, list[dict[str, str]]]:
    """
    Parse porcelain v1 -z entries, preserving rename origins.
    """

    proc = git(root, "status", "--porcelain=v1", "-z", "-uall", text=False)
    if proc.returncode:
        raise RuntimeError(proc.stderr.decode("utf-8", "replace"))

    raw = proc.stdout
    records: list[dict[str, str]] = []
    fields = raw.split(b"\0")
    index = 0

    while index < len(fields):
        field = fields[index]
        index += 1
        if not field:
            continue

        status = field[:2].decode("ascii")
        path = field[3:].decode("utf-8", "surrogateescape")
        record = {"git_status": status, "path": path}

        if "R" in status or "C" in status:
            if index >= len(fields):
                raise ValueError("truncated porcelain rename record")

            record["rename_from"] = fields[index].decode(
                "utf-8",
                "surrogateescape",
            )
            index += 1

        records.append(record)

    return raw, records


def relative_path(root: Path, path: Path) -> str:
    """
    Return a normalized repository-relative POSIX path.
    """

    return path.resolve().relative_to(root).as_posix()


def entry_fingerprint(root: Path, path: str) -> tuple[str, str]:
    """
    Fingerprint current content, a symlink target, or deleted Git content.
    """

    supplied = Path(path)
    candidate = supplied if supplied.is_absolute() else root / supplied
    if candidate.is_symlink():
        return sha256_bytes(
            path.encode() + b"\0symlink\0" + os.readlink(candidate).encode(),
        ), "symlink"

    if candidate.is_file():
        return sha256_bytes(
            path.encode() + b"\0file\0" + candidate.read_bytes(),
        ), "worktree"

    for spec, source in (
        (f":{path}", "index_deleted"),
        (f"HEAD:{path}", "head_deleted"),
    ):
        proc = git(root, "show", spec, text=False)
        if proc.returncode == 0:
            return sha256_bytes(
                path.encode() + b"\0deleted\0" + proc.stdout,
            ), source

    return sha256_bytes(path.encode() + b"\0missing\0"), "missing"


def path_matches(path: str, patterns: list[str]) -> bool:
    """
    Return whether a repository-relative path matches a configured glob.
    """

    for pattern in patterns:
        pieces: list[str] = []
        index = 0

        while index < len(pattern):
            if pattern.startswith("**/", index):
                pieces.append("(?:.*/)?")
                index += 3
            elif pattern.startswith("**", index):
                pieces.append(".*")
                index += 2
            elif pattern[index] == "*":
                pieces.append("[^/]*")
                index += 1
            elif pattern[index] == "?":
                pieces.append("[^/]")
                index += 1
            else:
                pieces.append(re.escape(pattern[index]))
                index += 1

        if re.fullmatch("".join(pieces), path):
            return True

    return False


def validate_target_path(path: object) -> str:
    """
    Validate one explicitly selected repository-relative path.
    """

    if not isinstance(path, str) or not path:
        raise ValueError("selected target path must be a non-empty string")

    pure = PurePosixPath(path)
    if (
        pure.is_absolute()
        or "\\" in path
        or pure.as_posix() != path
        or any(part in {".", ".."} for part in pure.parts)
    ):
        raise ValueError(f"invalid selected target path: {path!r}")

    return path


def validate_controlled_smoke_evidence(value: object) -> dict[str, object]:
    """
    Validate one narrow manifest-declared controlled smoke specification.

    Parameters
    ----------
    value : object
        Candidate controlled-smoke configuration from the rule manifest.

    Returns
    -------
    evidence : dict[str, object]
        Validated command, source, assertion, observation, and limit fields.

    Raises
    ------
    ValueError
        If any required field, type, marker group, or command form is invalid.
    """

    if (
        not isinstance(value, dict)
        or value.get("kind") != CONTROLLED_SMOKE_KIND
    ):
        raise ValueError(
            "supporting evidence must use controlled_smoke_summary_v1",
        )

    command = value.get("command")
    expected_command_fields = {
        "environment",
        "argv",
    }
    command_fields_valid = (
        isinstance(command, dict) and set(command) == expected_command_fields
    )

    if not command_fields_valid:
        raise ValueError("controlled smoke command needs environment and argv")

    environment, argv = command["environment"], command["argv"]
    if not isinstance(environment, dict) or not all(
        isinstance(key, str) and isinstance(item, str)
        for key, item in environment.items()
    ):
        raise ValueError(
            "controlled smoke command environment must map strings",
        )

    if not isinstance(argv, list) or argv != ["bash", "{path}"]:
        raise ValueError("controlled smoke command must be bash {path}")

    if (
        not isinstance(value.get("test_identifier"), str)
        or not value["test_identifier"]
    ):
        raise ValueError("controlled smoke evidence needs test_identifier")

    context = value.get("outer_context")
    expected_context_fields = {
        "require_wget_absent",
        "require_conda_available",
        "require_inactive_project_environment",
    }
    if (
        not isinstance(context, dict)
        or set(context) != expected_context_fields
        or not all(isinstance(item, bool) for item in context.values())
    ):
        raise ValueError("controlled smoke evidence has invalid outer_context")

    dependencies = value.get("source_dependencies")
    if not isinstance(dependencies, list) or not dependencies:
        raise ValueError("controlled smoke evidence needs source_dependencies")

    normalized_dependencies = [
        validate_target_path(path) for path in dependencies
    ]
    if len(normalized_dependencies) != len(set(normalized_dependencies)):
        raise ValueError("controlled smoke source_dependencies must be unique")

    groups = value.get("source_assertion_groups")
    required_groups = value.get("required_pass_groups")

    for label, collection in (
        ("source_assertion_groups", groups),
        ("required_pass_groups", required_groups),
    ):
        if (
            not isinstance(collection, dict)
            or not collection
            or not all(
                isinstance(name, str)
                and isinstance(markers, list)
                and markers
                and all(
                    isinstance(marker, str) and marker for marker in markers
                )
                for name, markers in collection.items()
            )
        ):
            raise ValueError(f"controlled smoke evidence has invalid {label}")

    observations = value.get("optional_observations", {})
    if not isinstance(observations, dict):
        raise ValueError(
            "controlled smoke optional_observations must be an object",
        )

    for name, observation in observations.items():
        if (
            not isinstance(name, str)
            or not name
            or not isinstance(observation, dict)
            or set(observation) != {"pass_message", "value"}
            or not isinstance(observation["pass_message"], str)
            or not observation["pass_message"]
            or not isinstance(observation["value"], (str, bool, int, float))
        ):
            raise ValueError(
                "controlled smoke evidence has invalid optional_observations",
            )

    if not isinstance(value.get("limitation"), str) or not value["limitation"]:
        raise ValueError("controlled smoke evidence needs limitation")

    return value


def validate_dependency_graph(
    edges: object,
    known_paths: set[str],
    requirements: object,
) -> dict[str, object]:
    """
    Validate declared edge identities and exact required relationships.

    Parameters
    ----------
    edges : object
        Candidate dependency-edge rows.
    known_paths : set[str]
        Declared endpoints that edges may reference.
    requirements : object
        Named exact edge sets that the graph must satisfy.

    Returns
    -------
    graph : dict[str, object]
        Canonical edge inventory and satisfied requirement names.

    Raises
    ------
    ValueError
        If rows, endpoints, identities, or required relationships are invalid.
    """

    if not isinstance(edges, list):
        raise ValueError("dependency_edges must be a list")

    identities: dict[tuple[str, str], dict[str, object]] = {}

    for edge in edges:
        if not isinstance(edge, dict):
            raise ValueError("dependency edge must be an object")

        source, target = edge.get("from"), edge.get("to")
        if source not in known_paths or target not in known_paths:
            raise ValueError(
                f"undeclared dependency endpoint: {source!r} -> {target!r}",
            )

        identity = (str(source), str(target))
        previous = identities.get(identity)

        if previous is not None:
            if previous == edge:
                raise ValueError(
                    (
                        f"duplicate dependency edge: {identity[0]} -> "
                        f"{identity[1]}"
                    ),
                )

            raise ValueError(
                (
                    f"contradictory dependency edge: {identity[0]} -> "
                    f"{identity[1]}"
                ),
            )

        identities[identity] = edge

    if not isinstance(requirements, list) or not requirements:
        raise ValueError("edge_requirements must be a nonempty list")

    satisfied: list[str] = []
    seen_names: set[str] = set()

    for requirement in requirements:
        if (
            not isinstance(requirement, dict)
            or not {"name", "edges"} <= set(requirement)
            or set(requirement) - {"name", "edges", "child"}
        ):
            raise ValueError(
                "dependency edge requirement needs name and edges",
            )

        name = requirement["name"]
        required_edges = requirement["edges"]
        if not isinstance(name, str) or not name or name in seen_names:
            raise ValueError(
                (
                    "dependency edge requirement names must be unique and "
                    "nonempty"
                ),
            )

        if not isinstance(required_edges, list) or not required_edges:
            raise ValueError(
                f"dependency edge requirement must contain edges: {name}",
            )

        for required in required_edges:
            expected_edge_fields = {
                "from",
                "to",
            }
            edge_fields_valid = (
                isinstance(required, dict)
                and set(required) == expected_edge_fields
            )

            if not edge_fields_valid:
                raise ValueError(f"invalid required dependency edge: {name}")

            identity = (required["from"], required["to"])
            if identity not in identities:
                raise ValueError(
                    (
                        f"missing required dependency edge ({name}): "
                        f"{identity[0]} -> "
                        f"{identity[1]}"
                    ),
                )

        seen_names.add(name)
        satisfied.append(name)

    return {
        "edge_count": len(edges),
        "satisfied_requirement_groups": satisfied,
    }


def validate_supplied_package(package: dict[str, object]) -> None:
    """
    Validate one configured five-artifact package identity.
    """

    semantic_path = validate_target_path(package["semantic_review_path"])
    if not semantic_path.startswith(
        "semantic_review/",
    ) or not semantic_path.endswith(".md"):
        raise ValueError(
            "package semantic_review_path must name semantic_review/*.md",
        )

    supplied = package["supplied_artifacts"]
    expected = {
        semantic_path,
        "findings.ndjson",
        "facts.ndjson",
        "adapter_limitations.ndjson",
        "pilot_report.json",
    }
    if (
        not isinstance(supplied, list)
        or len(supplied) != 5
        or set(supplied) != expected
    ):
        raise ValueError(
            (
                "package supplied_artifacts must preserve the five-artifact "
                "package"
            ),
        )


def linked_configuration_metadata(
    value: dict[str, object],
    selection_path: str,
) -> dict[str, object]:
    """
    Return canonical linked-package ownership and graph metadata.
    """

    package = value["package"]
    children = package["children"]
    target_ownership = [
        {
            "path": target["path"],
            "child": target.get("child")
            or next(
                child["package_id"]
                for child in children
                if target.get("subcohort") in child["subcohorts"]
            ),
        }
        for target in value["targets"]
    ]
    context_ownership = [
        {"path": context["path"], "child": context["child"]}
        for context in value.get("context", [])
    ]
    graph = list(value["dependency_edges"])

    return {
        "bundle_id": package["bundle_id"],
        "selection_path": selection_path,
        "config_fingerprint": sha256_bytes(json_bytes(value)),
        "target_ownership_fingerprint": sha256_bytes(
            json_bytes(target_ownership),
        ),
        "graph_ownership_fingerprint": sha256_bytes(json_bytes(graph)),
        "umbrella_target_ownership": target_ownership,
        "umbrella_context_ownership": context_ownership,
        "umbrella_dependency_edges": graph,
        "umbrella_target_count": len(target_ownership),
        "umbrella_context_count": len(context_ownership),
        "umbrella_edge_count": len(graph),
        "child_ids": [child["package_id"] for child in children],
    }


def _validate_selection_package(
    value: dict[str, Any],
    package_child: str | None,
    umbrella_run_id: str | None,
) -> tuple[dict[str, Any] | None, list[dict[str, Any]], set[str], bool]:
    """
    Validate package identity, limits, and optional linked children.

    Parameters
    ----------
    value : dict[str, Any]
        Parsed target-selection payload containing optional package metadata.
    package_child : str | None
        Requested linked child package, when one is selected.
    umbrella_run_id : str | None
        Run identity supplied for linked child generation.

    Returns
    -------
    package, children, child_ids, was_selected : tuple[
        dict[str, Any] | None, list[dict[str, Any]], set[str], bool
    ]
        Package metadata, validated children, child identifiers, and whether
        the package is linked.

    Raises
    ------
    ValueError
        If package fields, limits, children, or linked options are invalid.
    """

    package = value.get("package")
    children: list[dict[str, Any]] = []
    child_ids: set[str] = set()
    linked = isinstance(package, dict) and "children" in package
    if package is None:
        return None, children, child_ids, linked

    if not isinstance(package, dict):
        raise ValueError("package configuration must be an object")

    required_package = {
        "bundle_id",
        "title",
        "semantic_review_path",
        "supplied_artifacts",
        "statements",
        "semantic_only_topics",
        "size_limits",
    }
    allowed_package = required_package | ({"children"} if linked else set())
    if set(package) != allowed_package:
        raise ValueError(
            "package configuration has missing or unknown fields",
        )

    for field in ("bundle_id", "title"):
        if not isinstance(package[field], str) or not package[field]:
            raise ValueError(f"package {field} must be a nonempty string")

    for field in ("statements", "semantic_only_topics"):
        if not isinstance(package[field], list) or not all(
            isinstance(item, str) and item for item in package[field]
        ):
            raise ValueError(
                f"package {field} must contain nonempty strings",
            )

    limits = package["size_limits"]
    expected_limit_fields = {
        "semantic_markdown_max_bytes",
        "semantic_markdown_max_lines",
    }
    if (
        not isinstance(limits, dict)
        or set(limits) != expected_limit_fields
        or not all(type(item) is int and item > 0 for item in limits.values())
    ):
        raise ValueError(
            "package size_limits must contain positive byte and line limits",
        )

    validate_supplied_package(package)

    if not linked:
        if package_child is not None or umbrella_run_id is not None:
            raise ValueError(
                "linked-package options require package.children",
            )

        return package, children, child_ids, linked

    raw_children = package["children"]
    if (
        not isinstance(raw_children, list)
        or not raw_children
        or len(raw_children) > 2
    ):
        raise ValueError(
            "package.children must contain at most two child definitions",
        )

    child_fields = {
        "package_id",
        "title",
        "semantic_review_path",
        "supplied_artifacts",
        "subcohorts",
    }

    for child in raw_children:
        if not isinstance(child, dict) or set(child) != child_fields:
            raise ValueError("linked child has missing or unknown fields")

        if not all(
            isinstance(child[field], str) and child[field]
            for field in ("package_id", "title")
        ):
            raise ValueError(
                "linked child identity fields must be nonempty strings",
            )

        if (
            not isinstance(child["subcohorts"], list)
            or not child["subcohorts"]
            or not all(
                isinstance(item, str) and item for item in child["subcohorts"]
            )
        ):
            raise ValueError(
                "linked child subcohorts must be nonempty strings",
            )

        validate_supplied_package(child)

        children.append(dict(child))

    child_ids = {str(child["package_id"]) for child in children}
    if len(child_ids) != len(children):
        raise ValueError("linked child package ids must be unique")

    if package_child is not None and package_child not in child_ids:
        raise ValueError(f"unknown package child: {package_child}")

    return package, children, child_ids, linked


def _normalize_selection_targets(
    value: dict[str, Any],
    linked: bool,
    children: list[dict[str, Any]],
    child_ids: set[str],
) -> list[dict[str, Any]]:
    """
    Validate target roles, optional evidence, and linked ownership.

    Parameters
    ----------
    value : dict[str, Any]
        Parsed target-selection payload.
    linked : bool
        Whether every target must belong to one linked child.
    children : list[dict[str, Any]]
        Validated linked child definitions.
    child_ids : set[str]
        Accepted linked child identifiers.

    Returns
    -------
    targets : list[dict[str, Any]]
        Normalized unique targets with validated roles and ownership.

    Raises
    ------
    ValueError
        If target structure, roles, evidence, or ownership are invalid.
    """

    targets = value.get("targets")
    if not isinstance(targets, list):
        raise ValueError("target selection targets and context must be lists")

    normalized: list[dict[str, Any]] = []

    for item in targets:
        allowed_roles = {
            "primary",
            "supporting",
        }
        item_has_allowed_role = (
            isinstance(item, dict) and item.get("role") in allowed_roles
        )

        if not item_has_allowed_role:
            raise ValueError(
                "target selection entries need primary or supporting role",
            )

        allowed = {"path", "role", "evidence", "target_role", "subcohort"} | (
            {"child"} if linked else set()
        )
        if set(item) - allowed:
            raise ValueError("target selection entry has unknown fields")

        target = {
            "path": validate_target_path(item.get("path")),
            "role": str(item["role"]),
        }

        for field in ("target_role", "subcohort"):
            if field not in item:
                continue

            if not isinstance(item[field], str) or not item[field]:
                raise ValueError(
                    f"target selection {field} must be a nonempty string",
                )

            target[field] = item[field]

        if "evidence" in item:
            if target["role"] != "supporting":
                raise ValueError(
                    "controlled smoke evidence requires a supporting target",
                )

            target["evidence"] = validate_controlled_smoke_evidence(
                item["evidence"],
            )

        if linked:
            owners = [
                str(child["package_id"])
                for child in children
                if item.get("subcohort") in child["subcohorts"]
            ]
            owner = item.get("child")

            if owner is None and len(owners) == 1:
                owner = owners[0]

            if owner not in child_ids or owners != [owner]:
                raise ValueError(
                    (
                        "every linked target must have exactly one known "
                        "child owner"
                    ),
                )

            target["child"] = owner

        normalized.append(target)

    if not normalized or len({item["path"] for item in normalized}) != len(
        normalized,
    ):
        raise ValueError("target selection must contain unique targets")

    return normalized


def _normalize_selection_interfaces(
    value: dict[str, Any],
    target_paths: set[str],
) -> list[dict[str, Any]]:
    """
    Validate documented interface mappings and removed aliases.

    Parameters
    ----------
    value : dict[str, Any]
        Parsed target-selection payload.
    target_paths : set[str]
        Validated target paths available to interface mappings.

    Returns
    -------
    interfaces : list[dict[str, Any]]
        Normalized unique interface mappings.

    Raises
    ------
    ValueError
        If mappings, paths, or removed aliases violate the schema.
    """

    interfaces = value.get("interfaces", [])
    if not isinstance(interfaces, list):
        raise ValueError("target selection targets and context must be lists")

    normalized: list[dict[str, Any]] = []

    for item in interfaces:
        expected_interface_fields = {
            "path",
            "documentation_source",
            "removed_aliases",
        }
        interface_fields_valid = (
            isinstance(item, dict) and set(item) == expected_interface_fields
        )

        if not interface_fields_valid:
            raise ValueError(
                (
                    "interface entries need path, documentation_source, and "
                    "removed_aliases"
                ),
            )

        path = validate_target_path(item["path"])
        documentation_source = validate_target_path(
            item["documentation_source"],
        )
        removed_aliases = item["removed_aliases"]

        if (
            path not in target_paths
            or documentation_source not in target_paths
        ):
            raise ValueError(
                (
                    "interface paths and documentation sources must be "
                    "selected targets"
                ),
            )

        if not isinstance(removed_aliases, list) or not all(
            isinstance(alias, str)
            and re.fullmatch(r"--?[A-Za-z][A-Za-z0-9_-]*", alias)
            for alias in removed_aliases
        ):
            raise ValueError(
                (
                    "interface removed_aliases must contain exact option "
                    "spellings"
                ),
            )

        if len(removed_aliases) != len(set(removed_aliases)):
            raise ValueError("interface removed_aliases must be unique")

        normalized.append(
            {
                "path": path,
                "documentation_source": documentation_source,
                "removed_aliases": removed_aliases,
            },
        )

    if len({str(item["path"]) for item in normalized}) != len(normalized):
        raise ValueError("interface paths must be unique")

    return normalized


def _normalize_selection_context(
    value: dict[str, Any],
    linked: bool,
    child_ids: set[str],
    target_paths: set[str],
) -> tuple[list[str], list[dict[str, str]]]:
    """
    Validate context paths, uniqueness, and linked ownership.

    Parameters
    ----------
    value : dict[str, Any]
        Parsed target-selection payload.
    linked : bool
        Whether each context path must name one linked child.
    child_ids : set[str]
        Accepted linked child identifiers.
    target_paths : set[str]
        Validated targets from which context must remain distinct.

    Returns
    -------
    paths, records : tuple[list[str], list[dict[str, str]]]
        Normalized context paths and their explicit ownership records.

    Raises
    ------
    ValueError
        If context structure, paths, uniqueness, or ownership are invalid.
    """

    raw_context = value.get("context", [])
    if not isinstance(raw_context, list):
        raise ValueError("target selection targets and context must be lists")

    records: list[dict[str, str]] = []

    for item in raw_context:
        if linked:
            if (
                not isinstance(item, dict)
                or set(item) != {"path", "child"}
                or item.get("child") not in child_ids
            ):
                raise ValueError(
                    (
                        "every linked context path must have exactly one "
                        "known child owner"
                    ),
                )

            records.append(
                {
                    "path": validate_target_path(item["path"]),
                    "child": str(item["child"]),
                },
            )
        else:
            records.append(
                {"path": validate_target_path(item), "child": ""},
            )

    paths = [item["path"] for item in records]
    if len(paths) != len(set(paths)) or set(paths) & target_paths:
        raise ValueError(
            "context paths must be unique and distinct from targets",
        )

    return paths, records


def _validate_evidence_selection(
    value: dict[str, Any],
    linked: bool,
    normalized_context: list[str],
    context_records: list[dict[str, str]],
) -> tuple[dict[str, Any], dict[str, Any], list[dict[str, Any]]]:
    """
    Validate configured target and context semantic-unit selectors.

    Parameters
    ----------
    value : dict[str, Any]
        Parsed target-selection payload.
    linked : bool
        Whether context selectors require child ownership.
    normalized_context : list[str]
        Validated context paths.
    context_records : list[dict[str, str]]
        Context paths paired with their linked owners.

    Returns
    -------
    evidence_selection, target_selector, context_selectors : tuple[
        dict[str, Any], dict[str, Any], list[dict[str, Any]]
    ]
        Complete evidence selection, target selector, and context selectors.

    Raises
    ------
    ValueError
        If selector shape, paths, names, or ownership are invalid.
    """

    evidence_selection = value.get("evidence_selection")

    if not isinstance(evidence_selection, dict) or set(
        evidence_selection,
    ) != {"target_selector", "context_selectors"}:
        raise ValueError(
            "configured package needs target and context evidence selectors",
        )

    target_selector = evidence_selection["target_selector"]
    if (
        not isinstance(target_selector, dict)
        or set(target_selector) != {"kind", "selection_reason"}
        or target_selector.get("kind") != "changed_shell_semantic_units"
        or not isinstance(target_selector.get("selection_reason"), str)
        or not target_selector["selection_reason"]
    ):
        raise ValueError("invalid configured target semantic-unit selector")

    context_selectors = evidence_selection["context_selectors"]
    if not isinstance(context_selectors, list):
        raise ValueError("context_selectors must be a list")

    for selector in context_selectors:
        expected_fields = {
            "path",
            "kind",
            "name",
            "selection_reason",
        } | ({"child"} if linked else set())
        if (
            not isinstance(selector, dict)
            or set(selector) != expected_fields
            or selector.get("kind") != "shell_function"
        ):
            raise ValueError(
                "invalid configured context semantic-unit selector",
            )

        selector_path = validate_target_path(selector["path"])
        if selector_path not in normalized_context:
            raise ValueError(
                (
                    "context semantic-unit selector path is not declared "
                    "context"
                ),
            )

        if not all(
            isinstance(selector.get(field), str) and selector[field]
            for field in ("name", "selection_reason")
        ):
            raise ValueError(
                (
                    "context semantic-unit selector needs name and "
                    "selection_reason"
                ),
            )

        if linked:
            owner = next(
                item["child"]
                for item in context_records
                if item["path"] == selector["path"]
            )
            if selector.get("child") != owner:
                raise ValueError(
                    "context selector ownership must match its context path",
                )

    return evidence_selection, target_selector, context_selectors


def _validate_rule_scopes(value: dict[str, Any]) -> list[dict[str, Any]]:
    """
    Validate configured rule identities, patterns, and enforcement.

    Parameters
    ----------
    value : dict[str, Any]
        Parsed target-selection payload.

    Returns
    -------
    scopes : list[dict[str, Any]]
        Unique validated rule-scope definitions.

    Raises
    ------
    ValueError
        If scope identities, paths, classifications, or enforcement fail.
    """

    rule_scopes = value.get("rule_scopes")
    if not isinstance(rule_scopes, list) or not rule_scopes:
        raise ValueError("configured package needs rule_scopes")

    for scope in rule_scopes:
        expected_scope_fields = {
            "rule_id",
            "paths",
            "classification",
            "enforcement",
        }
        scope_fields_valid = (
            isinstance(scope, dict) and set(scope) == expected_scope_fields
        )

        if not scope_fields_valid:
            raise ValueError("invalid configured rule scope")

        if (
            not isinstance(scope["rule_id"], str)
            or not scope["rule_id"]
            or not isinstance(scope["classification"], str)
            or not scope["classification"]
        ):
            raise ValueError(
                "configured rule scope needs rule_id and classification",
            )

        if (
            not isinstance(scope["paths"], list)
            or not scope["paths"]
            or not all(
                isinstance(pattern, str) and pattern
                for pattern in scope["paths"]
            )
        ):
            raise ValueError(
                "configured rule scope paths must be nonempty patterns",
            )

        if scope["enforcement"] not in {"strict", "advisory"}:
            raise ValueError(
                (
                    "configured rule scope enforcement must be strict or "
                    "advisory"
                ),
            )

    if len({scope["rule_id"] for scope in rule_scopes}) != len(rule_scopes):
        raise ValueError("configured rule scope ids must be unique")

    return rule_scopes


def _validate_semantic_questions(
    value: dict[str, Any],
    linked: bool,
    known_paths: set[str],
    child_ids: set[str],
) -> list[dict[str, Any]]:
    """
    Validate semantic questions and their declared applicability.

    Parameters
    ----------
    value : dict[str, Any]
        Parsed target-selection payload.
    linked : bool
        Whether questions must declare linked-child applicability.
    known_paths : set[str]
        Target and context paths allowed in each question.
    child_ids : set[str]
        Accepted linked child identifiers.

    Returns
    -------
    questions : list[dict[str, Any]]
        Validated semantic questions.

    Raises
    ------
    ValueError
        If question text, paths, or child applicability are invalid.
    """

    questions = value.get("semantic_questions")
    if not isinstance(questions, list) or not questions:
        raise ValueError("configured package needs semantic_questions")

    for question in questions:
        expected_fields = {
            "rule_id",
            "topic",
            "question",
            "paths",
        } | ({"children"} if linked else set())
        if not isinstance(question, dict) or set(question) != expected_fields:
            raise ValueError("invalid configured semantic question")

        if not all(
            isinstance(question[field], str) and question[field]
            for field in ("rule_id", "topic", "question")
        ):
            raise ValueError(
                "semantic question text fields must be nonempty",
            )

        if (
            not isinstance(question["paths"], list)
            or not question["paths"]
            or not set(question["paths"]) <= known_paths
        ):
            raise ValueError(
                (
                    "semantic question paths must be declared targets or "
                    "context"
                ),
            )

        if linked and (
            not isinstance(question["children"], list)
            or not question["children"]
            or not set(question["children"]) <= child_ids
            or len(question["children"]) != len(set(question["children"]))
        ):
            raise ValueError(
                (
                    "semantic question child applicability must name known "
                    "children"
                ),
            )

    return questions


def _validate_dependency_edges(
    value: dict[str, Any],
    linked: bool,
    known_paths: set[str],
    child_ids: set[str],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]] | None]:
    """
    Validate dependency edges, classifications, and child requirements.

    Parameters
    ----------
    value : dict[str, Any]
        Parsed target-selection payload.
    linked : bool
        Whether edges and requirements need linked-child ownership.
    known_paths : set[str]
        Target and context paths allowed as edge endpoints.
    child_ids : set[str]
        Accepted linked child identifiers.

    Returns
    -------
    edges, edge_requirements : tuple[
        list[dict[str, Any]], list[dict[str, Any]] | None
    ]
        Validated dependency edges and optional edge requirements.

    Raises
    ------
    ValueError
        If edge shape, classification, evidence, or ownership are invalid.
    """

    edges = value.get("dependency_edges")
    edge_requirements = value.get("edge_requirements")
    edge_kinds = {
        "source",
        "bootstrap source",
        "direct call",
        "transitive call",
        "test-only",
    }
    reachability = {
        "statically observed",
        "conditional",
        "dynamic or uncertain",
        "loaded only",
    }
    runtime_statuses = {"exercised", "not observed exercised", "unknown"}

    if not isinstance(edges, list) or not edges:
        raise ValueError("configured package needs dependency_edges")

    for edge in edges:
        required = {
            "from",
            "to",
            "kind",
            "reachability",
            "runtime_status",
        } | ({"child"} if linked else set())
        if (
            not isinstance(edge, dict)
            or not required <= set(edge)
            or set(edge) - (required | {"evidence_reference"})
        ):
            raise ValueError("invalid configured dependency edge")

        if edge["from"] not in known_paths or edge["to"] not in known_paths:
            raise ValueError(
                (
                    "dependency edge endpoints must be declared targets or "
                    "context"
                ),
            )

        if (
            edge["kind"] not in edge_kinds
            or edge["reachability"] not in reachability
            or edge["runtime_status"] not in runtime_statuses
        ):
            raise ValueError(
                "dependency edge has invalid separate classification",
            )

        if edge["runtime_status"] == "exercised" and not isinstance(
            edge.get("evidence_reference"),
            str,
        ):
            raise ValueError(
                "exercised dependency edge needs evidence_reference",
            )

        if linked and edge.get("child") not in child_ids:
            raise ValueError(
                (
                    "every linked dependency edge must have exactly one known "
                    "child owner"
                ),
            )

    if linked:
        if not isinstance(edge_requirements, list):
            raise ValueError("edge_requirements must be a list")

        for requirement in edge_requirements:
            if (
                not isinstance(requirement, dict)
                or requirement.get("child") not in child_ids
            ):
                raise ValueError(
                    (
                        "every linked edge requirement must have exactly one "
                        "known child owner"
                    ),
                )

    return edges, edge_requirements


def _validate_adapter_limitations(
    value: dict[str, Any],
    linked: bool,
    child_ids: set[str],
) -> list[dict[str, Any]]:
    """
    Validate adapter limitation text and linked child applicability.

    Parameters
    ----------
    value : dict[str, Any]
        Parsed target-selection payload.
    linked : bool
        Whether each limitation must name linked children.
    child_ids : set[str]
        Accepted linked child identifiers.

    Returns
    -------
    limitations : list[dict[str, Any]]
        Validated adapter limitation records.

    Raises
    ------
    ValueError
        If limitation fields or child applicability are invalid.
    """

    limitations = value.get("adapter_limitations")
    if not isinstance(limitations, list):
        raise ValueError("configured adapter_limitations must be a list")

    for limitation in limitations:
        expected_fields = {
            "rule_id",
            "adapter",
            "limitation",
        } | ({"children"} if linked else set())
        if (
            not isinstance(limitation, dict)
            or set(limitation) != expected_fields
        ):
            raise ValueError("invalid configured adapter limitation")

        if not all(
            isinstance(limitation.get(field), str) and limitation[field]
            for field in ("rule_id", "adapter", "limitation")
        ):
            raise ValueError("invalid configured adapter limitation")

        if linked:
            limitation_children = limitation["children"]
            if (
                not isinstance(limitation_children, list)
                or not limitation_children
            ):
                raise ValueError(
                    "adapter limitation must name child applicability",
                )

            unique_children = set(limitation_children)
            unknown_children = unique_children - child_ids
            has_duplicates = len(limitation_children) != len(unique_children)
            if unknown_children or has_duplicates:
                raise ValueError(
                    (
                        "adapter limitation child applicability must name "
                        "known children"
                    ),
                )

    return limitations


@dataclass
class _LinkedSelectionContext:
    """
    Carry validated umbrella state into linked-child derivation.
    """

    selection: dict[str, Any]
    value: dict[str, Any]
    selection_relative_path: str
    package: dict[str, Any]
    children: list[dict[str, Any]]
    package_child: str | None
    umbrella_run_id: str | None
    targets: list[dict[str, Any]]
    context_records: list[dict[str, str]]
    target_selector: dict[str, Any]
    context_selectors: list[dict[str, Any]]
    questions: list[dict[str, Any]]
    edges: list[dict[str, Any]]
    edge_requirements: list[dict[str, Any]] | None
    limitations: list[dict[str, Any]]
    rule_scopes: list[dict[str, Any]]
    known_paths: set[str]


def _apply_linked_child_selection(
    context: _LinkedSelectionContext,
) -> None:
    """
    Derive one linked child package from the validated umbrella selection.

    Parameters
    ----------
    context : _LinkedSelectionContext
        Validated umbrella selection and requested child identity.
    """

    selection = context.selection
    value = context.value
    package = context.package
    children = context.children
    package_child = context.package_child

    targets = context.targets
    context_records = context.context_records
    target_selector = context.target_selector
    context_selectors = context.context_selectors

    questions = context.questions
    edges = context.edges
    edge_requirements = context.edge_requirements

    limitations = context.limitations
    rule_scopes = context.rule_scopes
    known_paths = context.known_paths

    linkage = linked_configuration_metadata(
        value,
        context.selection_relative_path,
    )
    selection["linked_umbrella"] = linkage
    if package_child is None:
        return

    child = next(
        item for item in children if item["package_id"] == package_child
    )
    sibling = next(
        item for item in children if item["package_id"] != package_child
    )

    child_targets = [
        target for target in targets if target["child"] == package_child
    ]
    child_context_records = [
        item for item in context_records if item["child"] == package_child
    ]
    child_context = [item["path"] for item in child_context_records]

    child_edges = [edge for edge in edges if edge["child"] == package_child]
    child_requirements = [
        item
        for item in (edge_requirements or [])
        if item["child"] == package_child
    ]
    child_questions = [
        item for item in questions if package_child in item["children"]
    ]
    child_limitations = [
        item for item in limitations if package_child in item["children"]
    ]

    owned_paths = {str(target["path"]) for target in child_targets}
    owned_paths.update(child_context)
    sibling_cross_edges = []

    for edge in edges:
        if edge["child"] == package_child:
            continue

        from_owned = edge["from"] in owned_paths
        to_owned = edge["to"] in owned_paths

        if from_owned != to_owned:
            sibling_cross_edges.append(edge)

    selection["targets"] = child_targets
    selection["context"] = child_context
    selection["evidence_selection"] = {
        "target_selector": target_selector,
        "context_selectors": [
            selector
            for selector in context_selectors
            if selector["child"] == package_child
        ],
    }

    selection["semantic_questions"] = child_questions
    selection["dependency_edges"] = child_edges
    selection["edge_requirements"] = child_requirements
    selection["dependency_graph_validation"] = validate_dependency_graph(
        child_edges,
        known_paths,
        child_requirements,
    )
    selection["adapter_limitations"] = child_limitations

    selection["rule_scopes"] = [
        {
            **scope,
            "paths": [
                pattern
                for pattern in scope["paths"]
                if any(
                    path_matches(str(target["path"]), [pattern])
                    for target in child_targets
                )
            ],
        }
        for scope in rule_scopes
    ]

    selection["package"] = {
        "bundle_id": package["bundle_id"],
        "title": child["title"],
        "semantic_review_path": child["semantic_review_path"],
        "supplied_artifacts": child["supplied_artifacts"],
        "statements": package["statements"],
        "semantic_only_topics": package["semantic_only_topics"],
        "size_limits": package["size_limits"],
    }

    selection["linked_package"] = {
        **linkage,
        "child_package_id": package_child,
        "sibling_package_id": sibling["package_id"],
        "umbrella_run_id": context.umbrella_run_id,
        "child_target_count": len(child_targets),
        "child_context_count": len(child_context),
        "child_edge_count": len(child_edges),
        "child_primary_paths": [target["path"] for target in child_targets],
        "child_context_paths": child_context,
        "child_dependency_edges": child_edges,
        "sibling_cross_package_edges": sibling_cross_edges,
    }


def load_target_selection(
    args: argparse.Namespace,
    root: Path,
) -> dict[str, object] | None:
    """
    Load explicit target roles and declared read-only context.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments or an explicit argument sequence.
    root : Path
        Repository root used to resolve maintained paths.

    Returns
    -------
    selection : dict[str, object] | None
        Validated target selection, or None when none was requested.

    Raises
    ------
    ValueError
        If a supplied value violates the validated contract.
    """

    path_values = getattr(args, "path_values", None)
    paths_from = getattr(args, "paths_from", None)
    package_child = getattr(args, "package_child", None)
    umbrella_run_id = getattr(args, "umbrella_run_id", None)

    if path_values and paths_from:
        raise ValueError("--path and --paths-from are mutually exclusive")

    if path_values:
        paths = [validate_target_path(path) for path in path_values]
        if len(paths) != len(set(paths)):
            raise ValueError("selected target paths must be unique")

        return {
            "schema_version": 1,
            "name": "explicit-path-selection",
            "targets": [{"path": path, "role": "primary"} for path in paths],
            "context": [],
        }

    if not paths_from:
        return None

    selection_path = (
        (root / paths_from).resolve()
        if not paths_from.is_absolute()
        else paths_from.resolve()
    )
    value = json.loads(selection_path.read_text(encoding="utf-8"))
    if (
        not isinstance(value, dict)
        or value.get("schema_version") != 1
        or not isinstance(value.get("name"), str)
    ):
        raise ValueError(
            "target selection must be a schema_version 1 object with name",
        )

    selection_relative_path = (
        selection_path.as_posix()
        if selection_path.is_absolute() and root not in selection_path.parents
        else relative_path(root, selection_path)
    )
    package, children, child_ids, linked = _validate_selection_package(
        value,
        package_child,
        umbrella_run_id,
    )

    normalized = _normalize_selection_targets(
        value,
        linked,
        children,
        child_ids,
    )
    target_paths = {str(item["path"]) for item in normalized}
    normalized_interfaces = _normalize_selection_interfaces(
        value,
        target_paths,
    )

    normalized_context, normalized_context_records = (
        _normalize_selection_context(
            value,
            linked,
            child_ids,
            target_paths,
        )
    )

    selection = {
        "schema_version": 1,
        "name": value["name"],
        "targets": normalized,
        "context": normalized_context,
        "interfaces": normalized_interfaces,
        "selection_path": selection_relative_path,
    }

    if package is not None:
        selection["package"] = dict(package)

        (
            evidence_selection,
            target_selector,
            context_selectors,
        ) = _validate_evidence_selection(
            value,
            linked,
            normalized_context,
            normalized_context_records,
        )

        selection["evidence_selection"] = evidence_selection

        rule_scopes = _validate_rule_scopes(value)
        selection["rule_scopes"] = rule_scopes

        known_paths = target_paths | set(normalized_context)
        questions = _validate_semantic_questions(
            value,
            linked,
            known_paths,
            child_ids,
        )
        selection["semantic_questions"] = questions

        edges, edge_requirements = _validate_dependency_edges(
            value,
            linked,
            known_paths,
            child_ids,
        )
        selection["dependency_graph_validation"] = validate_dependency_graph(
            edges,
            known_paths,
            edge_requirements,
        )
        selection["dependency_edges"] = edges
        selection["edge_requirements"] = edge_requirements

        configured_limitations = _validate_adapter_limitations(
            value,
            linked,
            child_ids,
        )
        selection["adapter_limitations"] = configured_limitations

        if linked:
            _apply_linked_child_selection(
                _LinkedSelectionContext(
                    selection=selection,
                    value=value,
                    selection_relative_path=selection_relative_path,
                    package=package,
                    children=children,
                    package_child=package_child,
                    umbrella_run_id=umbrella_run_id,
                    targets=normalized,
                    context_records=normalized_context_records,
                    target_selector=target_selector,
                    context_selectors=context_selectors,
                    questions=questions,
                    edges=edges,
                    edge_requirements=edge_requirements,
                    limitations=configured_limitations,
                    rule_scopes=rule_scopes,
                    known_paths=known_paths,
                ),
            )

        selection["rule_scope_report"] = configured_rule_scope_report(
            selection,
        )

    return selection


def selected_git_state(
    root: Path,
    path: str,
    status_rows: dict[str, dict[str, str]],
) -> dict[str, object]:
    """
    Classify a selected path without hiding staged/worktree combinations.
    """

    record = status_rows.get(path)
    labels: list[str] = []

    if record is None:
        labels.append("clean_context")
        status = "  "
    else:
        status = record["git_status"]

        if status == "??":
            labels.append("untracked")
        else:
            if status[0] != " ":
                labels.append("tracked_staged")

            if status[1] == "M":
                labels.append("tracked_modified")

            if "D" in status:
                labels.append("tracked_deleted")

            if "R" in status or "C" in status:
                labels.append("tracked_renamed")

    return {
        "git_status": status,
        "git_state_labels": labels,
        "rename_from": record.get("rename_from") if record else None,
    }


def build_target_records(
    root: Path,
    selection: dict[str, object],
    status_rows: dict[str, dict[str, str]],
) -> list[dict[str, object]]:
    """
    Return selected targets with current content fingerprints and state.
    """

    records: list[dict[str, object]] = []

    for item in selection["targets"]:
        path = str(item["path"])
        state = selected_git_state(root, path, status_rows)
        candidate = root / path
        if not candidate.is_file() or candidate.is_symlink():
            raise ValueError(
                f"selected target is not a regular readable file: {path}",
            )

        fingerprint, source = entry_fingerprint(root, path)
        record = {
            "schema_version": PILOT_ARTIFACT_SCHEMA_VERSION,
            "path": path,
            "role": item["role"],
            "content_fingerprint": fingerprint,
            "content_source": source,
            **state,
            "checks_run": [],
            "findings_count": 0,
            "whitespace_coverage": [],
        }

        for field in ("target_role", "subcohort"):
            if field in item:
                record[field] = item[field]

        if "evidence" in item:
            record["evidence"] = item["evidence"]

        records.append(record)

    return records


def validate_baseline(
    baseline: object,
    config: dict[str, object],
) -> dict[str, object]:
    """
    Validate every immutable baseline cohort invariant with precise errors.

    Parameters
    ----------
    baseline : object
        Baseline cohort payload to validate.
    config : dict[str, object]
        Validated policy or runtime configuration.

    Returns
    -------
    baseline : dict[str, object]
        The validated baseline payload.

    Raises
    ------
    BaselineValidationError
        If any immutable baseline invariant is violated.
    """

    errors: list[str] = []
    if not isinstance(baseline, dict):
        raise BaselineValidationError(
            "baseline cohort validation failed:\n- root: expected JSON object",
        )

    if baseline.get("schema_version") != BASELINE_SCHEMA_VERSION:
        errors.append(
            (
                f"schema_version: expected {BASELINE_SCHEMA_VERSION!r}, got "
                f"{baseline.get('schema_version')!r}"
            ),
        )

    if baseline.get("cohort") != BASELINE_COHORT:
        errors.append(
            (
                f"cohort: expected {BASELINE_COHORT!r}, got "
                f"{baseline.get('cohort')!r}"
            ),
        )

    declared = baseline.get("authored_path_count")

    if type(declared) is not int or declared != BASELINE_COUNT:
        errors.append(
            (
                f"authored_path_count: expected integer {BASELINE_COUNT}, got "
                f"{declared!r}"
            ),
        )

    entries = baseline.get("paths")

    if not isinstance(entries, list):
        errors.append("paths: expected list")
        entries = []

    if len(entries) != BASELINE_COUNT:
        errors.append(
            f"paths length: expected {BASELINE_COUNT}, got {len(entries)}",
        )

    paths: list[str] = []
    artifact_patterns = list(
        config.get("artifacts", {}).get("bytecode_globs", []),
    ) + list(config.get("artifacts", {}).get("other_generated_globs", []))

    for index, entry in enumerate(entries):
        prefix = f"paths[{index}]"

        if not isinstance(entry, dict):
            errors.append(
                f"{prefix}: expected object with path and initial_git_status",
            )

            continue

        path = entry.get("path")
        status = entry.get("initial_git_status")

        if not isinstance(path, str):
            errors.append(
                f"{prefix}.path: expected string, got {type(path).__name__}",
            )
        else:
            pure = PurePosixPath(path)

            if not path:
                errors.append(f"{prefix}.path: must not be empty")

            if "\\" in path:
                errors.append(f"{prefix}.path: must use POSIX separators")

            if pure.is_absolute():
                errors.append(
                    (
                        f"{prefix}.path: must be repository-relative, got "
                        f"absolute path "
                        f"{path!r}"
                    ),
                )

            if any(part in {".", ".."} for part in pure.parts):
                errors.append(
                    f"{prefix}.path: must not contain '.' or '..' traversal",
                )

            if pure.as_posix() != path:
                errors.append(f"{prefix}.path: must be normalized POSIX path")

            if path_matches(path, artifact_patterns):
                errors.append(
                    (
                        f"{prefix}.path: configured generated-artifact path "
                        f"is not allowed: "
                        f"{path!r}"
                    ),
                )

            paths.append(path)

        if not isinstance(status, str) or status not in BASELINE_STATUSES:
            errors.append(
                (
                    f"{prefix}.initial_git_status: expected one of "
                    f"{sorted(BASELINE_STATUSES)!r}, "
                    f"got "
                    f"{status!r}"
                ),
            )

    if len(paths) != len(set(paths)):
        errors.append("paths: duplicate path values are not allowed")

    if paths != sorted(paths):
        errors.append("paths: entries must be sorted by normalized path")

    if type(declared) is int and len(set(paths)) != declared:
        errors.append(
            f"unique path count: declared {declared}, found {len(set(paths))}",
        )

    if errors:
        raise BaselineValidationError(
            "baseline cohort validation failed:\n"
            + "\n".join(f"- {error}" for error in errors),
        )

    return baseline


class EvidenceReader:
    """
    Read explicitly consumed public/private evidence with fingerprints.
    """

    def __init__(self, public_root: Path, private_root: Path) -> None:
        """
        Initialize evidence roots and an empty consumption ledger.

        Parameters
        ----------
        public_root : Path
            Root of the public repository whose evidence may be consumed.
        private_root : Path
            Root of the paired private repository whose evidence may be
            consumed.
        """

        self.public_root = public_root.resolve()
        self.private_root = private_root.resolve()
        self.records: dict[tuple[str, str], dict[str, object]] = {}

    def read(self, repository: str, relative: str, purpose: str) -> str:
        """
        Read one declared evidence file.

        Symlinks and paths resolving outside the selected repository root are
        rejected.

        Parameters
        ----------
        repository : str
            Declared repository name, either 'public' or 'private'.
        relative : str
            Repository-relative evidence path.
        purpose : str
            Inspectable reason for consuming the evidence.

        Returns
        -------
        text : str
            Evidence text decoded with replacement for invalid bytes.

        Raises
        ------
        ValueError
            If the path is a symlink, missing, or outside its repository.
        """

        root = (
            self.public_root if repository == "public" else self.private_root
        )
        candidate = root / relative
        resolved = candidate.resolve()
        if (
            root not in (resolved, *resolved.parents)
            or candidate.is_symlink()
            or not candidate.is_file()
        ):
            raise ValueError(
                f"evidence path is not a regular in-root file: {relative}",
            )

        content = candidate.read_bytes()
        self.records[(repository, relative)] = {
            "repository": repository,
            "path": relative,
            "purpose": purpose,
            "fingerprint": sha256_bytes(
                relative.encode() + b"\0file\0" + content,
            ),
            "byte_count": len(content),
            "fingerprinted_before": True,
        }

        return content.decode("utf-8", "replace")

    def verify(self) -> list[dict[str, object]]:
        """
        Return consumed evidence whose current bytes no longer match.
        """

        stale: list[dict[str, object]] = []

        for (repository, relative), record in self.records.items():
            root = (
                self.public_root
                if repository == "public"
                else self.private_root
            )
            fingerprint, _ = entry_fingerprint(root, relative)

            record["fingerprinted_after"] = True
            record["verification_status"] = (
                "fresh" if fingerprint == record["fingerprint"] else "stale"
            )

            if fingerprint != record["fingerprint"]:
                stale.append(record)

        return stale


def git_visible_state(root: Path) -> dict[str, object]:
    """
    Capture Git-visible repository state without inspecting ignored files.
    """

    status, _ = parse_status(root)

    return {
        "branch": git(root, "branch", "--show-current").stdout.strip(),
        "head": git(root, "rev-parse", "HEAD").stdout.strip(),
        "status_digest": sha256_bytes(status),
        "staged_diff_digest": sha256_bytes(
            git(root, "diff", "--cached", "--binary", text=False).stdout,
        ),
        "worktree_diff_digest": sha256_bytes(
            git(root, "diff", "--binary", text=False).stdout,
        ),
    }


def validate_manifest(config: dict[str, object]) -> list[dict[str, object]]:
    """
    Validate manifest shape and strict-safe output policy.

    Parameters
    ----------
    config : dict[str, object]
        Parsed rule-registry document.

    Returns
    -------
    rules : list[dict[str, object]]
        Validated rule rows in registry order.

    Raises
    ------
    ValueError
        If registry fields, ownership metadata, coverage, or execution policy
        are inconsistent.
    """

    if config.get("schema_version") != 2:
        raise ValueError("rule manifest schema_version must be 2")

    rules = list(config.get("rule", []))
    seen: set[str] = set()
    required = {
        "rule_id",
        "description",
        "normative_document",
        "normative_section",
        "owner_classification",
        "source_checker",
        "execution_kind",
        "execution_role",
        "coverage_relation",
        "covered_scope",
        "remaining_scope",
        "scope",
        "applicable_paths",
        "default_severity",
        "blocking",
        "semantic_review",
        "required_environment",
        "required_commands",
        "known_side_effects",
        "output_paths",
        "output_parser",
        "strict_availability",
        "current_exclusions_or_allowlists",
    }

    for rule in rules:
        if not isinstance(rule, dict):
            raise ValueError("rule manifest contains a non-object rule")

        missing = required - set(rule)
        if missing:
            raise ValueError(
                f"rule {rule.get('rule_id')} missing {sorted(missing)}",
            )

        rule_id = str(rule["rule_id"])
        if rule_id in seen:
            raise ValueError(f"duplicate rule id: {rule_id}")

        seen.add(rule_id)

        if rule["owner_classification"] not in {
            "deterministic",
            "advisory",
            "semantic-only",
        }:
            raise ValueError(f"invalid owner classification for {rule_id}")

        if rule["execution_role"] not in {
            "checker",
            "evidence_producer",
            "independent_evidence",
        }:
            raise ValueError(f"invalid execution role for {rule_id}")

        if rule["coverage_relation"] not in {"exact", "subset", "independent"}:
            raise ValueError(f"invalid coverage relation for {rule_id}")

        if (
            not isinstance(rule["covered_scope"], str)
            or not rule["covered_scope"].strip()
        ):
            raise ValueError(f"covered_scope must be nonempty for {rule_id}")

        if (
            not isinstance(rule["remaining_scope"], str)
            or not rule["remaining_scope"].strip()
        ):
            raise ValueError(f"remaining_scope must be nonempty for {rule_id}")

        if (
            rule["coverage_relation"] == "exact"
            and rule["remaining_scope"].rstrip(".") != "None"
        ):
            raise ValueError(
                f"exact coverage must declare no remaining scope: {rule_id}",
            )

        if (rule["coverage_relation"] == "independent") != (
            rule["execution_role"] == "independent_evidence"
        ):
            raise ValueError(
                (
                    f"independent relation and execution role must agree: "
                    f"{rule_id}"
                ),
            )

        if (
            rule["execution_kind"] == "adapter"
            and rule["coverage_relation"] == "exact"
            and not rule.get("parity_test")
        ):
            raise ValueError(f"exact adapter requires parity_test: {rule_id}")

        if (
            rule["scope"] == "repository"
            and "commands" not in rule
            and "command" not in rule
        ):
            raise ValueError(f"repository rule needs command: {rule_id}")

        if (
            rule["strict_availability"] == "safe_adapter"
            and rule["output_paths"]
        ):
            raise ValueError(
                (
                    f"strict-safe rule {rule_id} must not declare repository "
                    f"output_paths"
                ),
            )

        registry_path = rule.get("adapter_registry_path")
        reference_scopes = rule.get("adapter_reference_scopes")

        if registry_path is not None:
            validate_target_path(registry_path)

            if (
                rule.get("adapter") != "shell_help_pilot"
                or rule.get("adapter_mode") != "command_references"
            ):
                raise ValueError(
                    (
                        f"registry-backed adapter configuration is invalid for"
                        f" "
                        f"{rule_id}"
                    ),
                )

            if (
                not isinstance(reference_scopes, list)
                or not reference_scopes
                or not all(
                    isinstance(scope, str) and scope
                    for scope in reference_scopes
                )
            ):
                raise ValueError(
                    (
                        f"registry-backed adapter needs nonempty reference "
                        f"scopes: "
                        f"{rule_id}"
                    ),
                )

    return rules


def configured_rule_selection(
    all_rules: list[dict[str, object]],
    selection: dict[str, object] | None,
    requested_ids: list[str] | None,
) -> tuple[list[dict[str, object]], dict[str, list[str]]]:
    """
    Select declared package rules and return their additional path scopes.

    Parameters
    ----------
    all_rules : list[dict[str, object]]
        Validated registry rules available to the audit.
    selection : dict[str, object] | None
        Optional package selection with per-rule path scopes.
    requested_ids : list[str] | None
        Optional caller-requested subset of rule identifiers.

    Returns
    -------
    rules, configured_scopes : tuple[
        list[dict[str, object]], dict[str, list[str]]
    ]
        Selected rules in registry order and package-specific path scopes.

    Raises
    ------
    ValueError
        If package or caller selections name unavailable rule identifiers.
    """

    by_id = {str(rule["rule_id"]): rule for rule in all_rules}
    configured_scopes: dict[str, list[str]] = {}

    if selection and "package" in selection:
        configured_scopes = {
            str(scope["rule_id"]): list(scope["paths"])
            for scope in selection["rule_scopes"]
        }
        unknown = set(configured_scopes) - set(by_id)
        if unknown:
            raise ValueError(
                f"configured package has unknown rule ids: {sorted(unknown)}",
            )

        selected_ids = set(configured_scopes)
    else:
        selected_ids = set(by_id)

    if requested_ids:
        requested = set(requested_ids)
        unknown = requested - set(by_id)
        if unknown:
            raise ValueError(f"unknown rule ids: {sorted(unknown)}")

        outside = requested - selected_ids
        if outside:
            raise ValueError(
                (
                    f"requested rule ids are outside configured package scope:"
                    f" "
                    f"{sorted(outside)}"
                ),
            )

        selected_ids = requested

    return [
        rule for rule in all_rules if str(rule["rule_id"]) in selected_ids
    ], configured_scopes


def configured_rule_scope_report(
    selection: dict[str, object],
) -> dict[str, dict[str, object]]:
    """
    Return exact cohort targets matched by each configured rule scope.
    """

    targets = [str(target["path"]) for target in selection["targets"]]
    report: dict[str, dict[str, object]] = {}

    for scope in selection.get("rule_scopes", []):
        rule_id = str(scope["rule_id"])
        matched = [
            path
            for path in targets
            if path_matches(path, list(scope["paths"]))
        ]
        report[rule_id] = {
            "paths": list(scope["paths"]),
            "matched_targets": matched,
            "matched_target_count": len(matched),
            "classification": scope["classification"],
            "enforcement": scope["enforcement"],
        }

    return report


def classify_artifact(path: str, config: dict[str, object]) -> str | None:
    """
    Return generated artifact category, if configured.
    """

    artifacts = config["artifacts"]
    if path_matches(path, list(artifacts["bytecode_globs"])):
        return "generated_bytecode"

    if path_matches(path, list(artifacts["other_generated_globs"])):
        return "other_generated"

    return None


def subsystem_and_risk(path: str) -> tuple[str, str]:
    """
    Assign deliberately coarse review routing metadata.
    """

    if path == "src/protocol_chipseq_signal_norm/cli/compute_signal.py":
        return "compute_signal", "high"

    if path.startswith("lib/bash/help/"):
        return "shell_help", "medium"

    if path.startswith("bin/") or path.startswith("install/scripts/"):
        return "workflow_shell", "high"

    if path.startswith("lib/bash/"):
        return "shared_shell", "high"

    if path.startswith("src/protocol_chipseq_signal_norm/"):
        return "python_package", "high"

    if path.startswith("tests/"):
        return "tests", "medium"

    if path.startswith("docs/") or path.endswith(".md"):
        return "documentation", "low"

    if path.startswith("dev/"):
        return "audit", "medium"

    return "repository", "medium"


def load_path_relocations(
    root: Path,
    baseline_paths: set[str],
) -> dict[str, str]:
    """
    Load declared current-to-baseline path identities for repository moves.

    Parameters
    ----------
    root : Path
        Repository root containing relocation configuration and Git state.
    baseline_paths : set[str]
        Paths belonging to the committed cleanup cohort.

    Returns
    -------
    relocations : dict[str, str]
        Current paths mapped to unique historical cohort paths.

    Raises
    ------
    ValueError
        If relocation configuration is malformed, ambiguous, or duplicated.
    """

    source = root / "dev/config/path_relocations.json"
    if not source.is_file():
        return {}

    value = json.loads(source.read_text(encoding="utf-8"))
    if value.get("schema_version") != 1:
        raise ValueError("path relocation schema_version must be 1")

    exact = value.get("exact")
    origins = value.get("basename_origins")
    if not isinstance(exact, dict) or not isinstance(origins, list):
        raise ValueError("path relocations require exact and basename_origins")

    relocations = {
        str(target): str(origin) for target, origin in exact.items()
    }

    dirty_paths = {record["path"] for record in parse_status(root)[1]}

    for target in dirty_paths:
        if target in relocations or target in baseline_paths:
            continue

        matches = [
            f"{prefix}/{PurePosixPath(target).name}"
            for prefix in origins
            if f"{prefix}/{PurePosixPath(target).name}" in baseline_paths
        ]

        if len(matches) == 1:
            relocations[target] = matches[0]
        elif len(matches) > 1:
            raise ValueError(f"ambiguous basename relocation for {target}")

    unknown = sorted(set(relocations.values()) - baseline_paths)
    if unknown:
        raise ValueError(
            f"path relocations name unknown baseline paths: {unknown}",
        )

    duplicate_origins = [
        origin
        for origin, count in Counter(relocations.values()).items()
        if count > 1
    ]
    if duplicate_origins:
        raise ValueError(f"duplicate relocation origins: {duplicate_origins}")

    return relocations


def build_inventory(
    root: Path,
    config: dict[str, object],
    baseline: dict[str, object],
    rules: list[dict[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, str]]:
    """
    Return authored ledger drafts, generated artifacts, and dirty fingerprints.

    Parameters
    ----------
    root : Path
        Repository root whose dirty paths are inventoried.
    config : dict[str, object]
        Audit infrastructure and artifact-classification configuration.
    baseline : dict[str, object]
        Committed cleanup cohort and path identities.
    rules : list[dict[str, object]]
        Validated rules used to determine per-path applicability.

    Returns
    -------
    authored, generated, fingerprints : tuple[
        list[dict[str, object]], list[dict[str, object]], dict[str, str]
    ]
        Authored ledger rows, generated-artifact rows, and path fingerprints.
    """

    _, records = parse_status(root)
    baseline_paths = {item["path"] for item in baseline["paths"]}
    relocations = load_path_relocations(root, baseline_paths)
    relocated_origins = set(relocations.values())
    infrastructure = config["audit_infrastructure"]

    ledger: list[dict[str, object]] = []
    artifacts: list[dict[str, object]] = []
    fingerprints: dict[str, str] = {}

    for status in records:
        path = status["path"]
        if path in relocated_origins and not (root / path).exists():
            continue

        fingerprint, source = entry_fingerprint(root, path)

        fingerprints[path] = fingerprint
        generated = classify_artifact(path, config)

        if generated:
            artifacts.append(
                {
                    "path": path,
                    "git_status": status["git_status"],
                    "artifact_type": generated,
                    "classification_rule": generated,
                    "content_fingerprint": fingerprint,
                    "content_source": source,
                },
            )

            continue

        origin = status.get("rename_from")

        if path in baseline_paths:
            cohort, resolution, cohort_origin = (
                "baseline_cleanup",
                "direct_path",
                path,
            )
        elif path in relocations:
            cohort = "baseline_cleanup"
            resolution = "declared_relocation"
            cohort_origin = relocations[path]
        elif origin in baseline_paths:
            cohort, resolution, cohort_origin = (
                "baseline_cleanup",
                "git_rename",
                origin,
            )
        elif path_matches(
            path,
            list(infrastructure["paths"]),
        ) and not path_matches(path, list(infrastructure["exclude"])):
            cohort, resolution, cohort_origin = (
                "audit_infrastructure",
                "infrastructure_path",
                None,
            )
        else:
            cohort, resolution, cohort_origin = "later_added", "new_path", None

        subsystem, risk = subsystem_and_risk(path)
        deleted = "D" in status["git_status"] and not (root / path).exists()
        applicable = (
            []
            if deleted
            else [
                str(rule["rule_id"])
                for rule in rules
                if path_matches(path, list(rule["applicable_paths"]))
            ]
        )
        ledger.append(
            {
                "path": path,
                "git_status": status["git_status"],
                "generated": False,
                "subsystem": subsystem,
                "risk": risk,
                "applicable_rules": applicable,
                "checks_run": [],
                "findings_count": 0,
                "semantic_review_required": risk == "high"
                or any(
                    bool(rule["semantic_review"])
                    for rule in rules
                    if str(rule["rule_id"]) in applicable
                ),
                "human_review_required": True,
                "review_status": "unreviewed",
                "disposition": "unresolved",
                "content_fingerprint": fingerprint,
                "content_source": source,
                "scope_cohort": cohort,
                "cohort_origin_path": cohort_origin,
                "cohort_resolution": resolution,
            },
        )

    return (
        sorted(ledger, key=lambda row: str(row["path"])),
        sorted(artifacts, key=lambda row: str(row["path"])),
        fingerprints,
    )


def executable_ledger_paths(ledger: list[dict[str, object]]) -> list[str]:
    """
    Return current authored paths that strict adapters can read.
    """

    deleted_sources = {"index_deleted", "head_deleted", "missing"}
    return [
        str(row["path"])
        for row in ledger
        if row["content_source"] not in deleted_sources
    ]


def public_guard(
    root: Path,
    fingerprints: dict[str, str],
    evidence: EvidenceReader,
) -> dict[str, object]:
    """
    Capture the per-transaction public strict guard.
    """

    return {
        "git_visible": git_visible_state(root),
        "dirty_content_digest": sha256_bytes(json_bytes(fingerprints)),
        "consumed_public_evidence": sorted(
            (
                record
                for record in evidence.records.values()
                if record["repository"] == "public"
            ),
            key=lambda record: str(record["path"]),
        ),
    }


def write_ndjson(path: Path, rows: list[dict[str, object]]) -> None:
    """
    Write records as UTF-8 NDJSON.
    """

    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(
                json.dumps(row, sort_keys=True, ensure_ascii=False) + "\n",
            )


def write_json(path: Path, value: object) -> None:
    """
    Write stable, readable JSON.
    """

    path.write_text(
        json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def cohort_progress(
    root: Path,
    baseline: dict[str, object],
    ledger: list[dict[str, object]],
) -> dict[str, object]:
    """
    Reconcile all baseline paths without silently losing vanished identities.
    """

    by_origin = {
        str(row.get("cohort_origin_path")): row
        for row in ledger
        if row.get("cohort_origin_path")
    }
    rows: list[dict[str, object]] = []

    for item in baseline["paths"]:
        path = str(item["path"])
        row = by_origin.get(path)

        if row:
            rows.append(
                {
                    "baseline_path": path,
                    "state": "renamed"
                    if row["path"] != path
                    else "current_dirty",
                    "current_path": row["path"],
                    "review_status": row["review_status"],
                    "disposition": row["disposition"],
                },
            )
        elif (root / path).exists() or (root / path).is_symlink():
            rows.append(
                {
                    "baseline_path": path,
                    "state": "clean",
                    "current_path": None,
                    "review_status": "not_applicable",
                    "disposition": "not_currently_dirty",
                },
            )
        else:
            rows.append(
                {
                    "baseline_path": path,
                    "state": "missing_or_unmatched",
                    "current_path": None,
                    "review_status": "unresolved",
                    "disposition": "reconciliation_required",
                },
            )

    scopes = Counter(str(row["scope_cohort"]) for row in ledger)

    return {
        "baseline_authored_path_count": baseline["authored_path_count"],
        "baseline_state_counts": dict(
            Counter(str(row["state"]) for row in rows),
        ),
        "current_authored_counts": dict(scopes),
        "all_current_authored_files": len(ledger),
        "reviewed_current_authored_files": sum(
            row["review_status"] != "unreviewed" for row in ledger
        ),
        "baseline_paths": rows,
    }


def current_fingerprints(
    root: Path,
    config: dict[str, object],
    baseline: dict[str, object],
    all_rules: list[dict[str, object]],
) -> dict[str, str]:
    """
    Rebuild dirty-entry fingerprints for a transaction guard comparison.
    """

    return build_inventory(root, config, baseline, all_rules)[2]


def finding_record(
    run_id: str,
    rule: dict[str, object],
    fragment: dict[str, object],
    fingerprint: str,
    checker: str,
    extra: dict[str, object] | None = None,
) -> dict[str, object]:
    """
    Attach stable audit metadata to one parsed finding fragment.
    """

    evidence: dict[str, object] = {
        "source": checker,
        "text": str(fragment["evidence"]),
    }

    if extra:
        evidence.update(
            {
                key: value
                for key, value in extra.items()
                if key != "covered_paths"
            },
        )

        if "covered_paths" in extra:
            evidence["rule_applicable_target_paths"] = extra["covered_paths"]

        if "selection_context_paths" in extra:
            evidence["selection_context_paths"] = extra[
                "selection_context_paths"
            ]

    key = sha256_bytes(
        json_bytes(
            [
                rule["rule_id"],
                fragment["path"],
                fragment["line"],
                fragment["message"],
                fragment["evidence"],
                fingerprint,
            ],
        ),
    )

    return {
        "run_id": run_id,
        "finding_id": f"{run_id}:{key.removeprefix('sha256:')[:16]}",
        "rule_id": rule["rule_id"],
        "path": fragment["path"],
        "line": fragment["line"],
        "severity": rule["default_severity"],
        "confidence": "high"
        if rule["default_severity"] == "error"
        else "medium",
        "message": fragment["message"],
        "evidence": evidence,
        "checker": checker,
        "content_fingerprint": fingerprint,
        "status": "open",
    }


def executable_info(
    argv: list[str],
    root: Path,
    timeout: float,
    maximum: int,
) -> dict[str, object]:
    """
    Return path, version, and controlled availability for one executable.
    """

    path = shutil.which(argv[0])
    if not path:
        return {
            "available": False,
            "argv": argv,
            "path": None,
            "version": None,
        }

    result = run_command(
        [path, *argv[1:]],
        root,
        timeout_seconds=timeout,
        max_output_bytes=maximum,
    )

    return {
        "available": result.returncode == 0
        and not result.timed_out
        and not result.launch_error,
        "argv": result.argv,
        "path": path,
        "version": (result.stdout or result.stderr).strip(),
        "launch_error": result.launch_error,
        "timed_out": result.timed_out,
    }


def environment_preflight(
    public_root: Path,
    private_root: Path,
    selected_rules: list[dict[str, object]],
    timeout: float,
    maximum: int,
) -> dict[str, object]:
    """
    Record every required executable and project-environment dependency.
    """

    conda_path = shutil.which("conda")
    runner = {"executable": sys.executable, "version": sys.version}
    git_info = executable_info(
        ["git", "--version"],
        public_root,
        timeout,
        maximum,
    )
    bash_info = executable_info(
        ["bash", "--version"],
        public_root,
        timeout,
        maximum,
    )
    conda_info = executable_info(
        ["conda", "--version"],
        public_root,
        timeout,
        maximum,
    )
    project: dict[str, object] = {
        "available": False,
        "python": None,
        "pytest": None,
    }

    if conda_path:
        python_probe = run_command(
            [
                conda_path,
                "run",
                "-n",
                "env_protocol",
                "python",
                "-c",
                "import sys; print(sys.executable); print(sys.version)",
            ],
            public_root,
            dict(os.environ, PYTHONDONTWRITEBYTECODE="1"),
            timeout,
            maximum,
        )
        project["python"] = {
            "available": python_probe.returncode == 0
            and not python_probe.timed_out
            and not python_probe.launch_error,
            "argv": python_probe.argv,
            "stdout": python_probe.stdout,
            "stderr": python_probe.stderr,
            "timed_out": python_probe.timed_out,
            "launch_error": python_probe.launch_error,
        }
        pytest = run_command(
            [
                conda_path,
                "run",
                "-n",
                "env_protocol",
                "python",
                "-c",
                "import pytest; print(pytest.__version__)",
            ],
            public_root,
            dict(os.environ, PYTHONDONTWRITEBYTECODE="1"),
            timeout,
            maximum,
        )
        project["pytest"] = {
            "available": pytest.returncode == 0
            and not pytest.timed_out
            and not pytest.launch_error,
            "argv": pytest.argv,
            "version": pytest.stdout.strip(),
            "stderr": pytest.stderr,
            "timed_out": pytest.timed_out,
            "launch_error": pytest.launch_error,
        }
        project["available"] = bool(project["python"]["available"])

    commands = sorted(
        {
            command
            for rule in selected_rules
            if rule["strict_availability"] == "safe_adapter"
            for command in rule["required_commands"]
        },
    )
    required = {
        command: executable_info(
            [command, "--version"],
            public_root,
            timeout,
            maximum,
        )
        for command in commands
    }

    return {
        "runner_python": runner,
        "git": git_info,
        "bash": bash_info,
        "conda": conda_info,
        "env_protocol": project,
        "required_commands": required,
        "public_git": git_visible_state(public_root),
        "private_git": git_visible_state(private_root),
    }


def rule_available(
    rule: dict[str, object],
    preflight: dict[str, object],
) -> str | None:
    """
    Return controlled unavailability reason for a strict-safe rule.
    """

    for command in rule["required_commands"]:
        if not preflight["required_commands"][command]["available"]:
            return f"required command unavailable: {command}"

    if (
        rule["required_environment"] == "env_protocol"
        and not preflight["env_protocol"]["available"]
    ):
        return "env_protocol Python unavailable"

    return None


def rule_excerpt(
    root: Path,
    evidence: EvidenceReader,
    rule: dict[str, object],
    maximum: int,
) -> str:
    """
    Return the named Markdown section, bounded for prompt use.
    """

    text = evidence.read(
        "public",
        str(rule["normative_document"]),
        f"rule excerpt for {rule['rule_id']}",
    )
    lines = text.splitlines()
    target = str(rule["normative_section"])
    start = next(
        (
            index
            for index, line in enumerate(lines)
            if line.lstrip("#").strip() == target
        ),
        0,
    )
    level = len(lines[start]) - len(lines[start].lstrip("#")) if lines else 1
    end = next(
        (
            index
            for index in range(start + 1, len(lines))
            if lines[index].startswith("#")
            and len(lines[index]) - len(lines[index].lstrip("#")) <= level
        ),
        len(lines),
    )

    return "\n".join(lines[start:end])[:maximum]


def _complete_line_diff_window(
    hunk: str,
    line_index: int,
    budget: int,
) -> str | None:
    """
    Return one complete-line diff hunk window around a changed line.
    """

    lines = hunk.splitlines(keepends=True)
    if not lines or line_index <= 0 or line_index >= len(lines):
        return None

    start, end = max(1, line_index - 3), min(len(lines), line_index + 4)

    while True:
        candidate = [lines[0]]

        if start > 1:
            candidate.append("[earlier lines omitted from this hunk]\n")

        candidate.extend(lines[start:end])

        if end < len(lines):
            candidate.append("[later lines omitted from this hunk]\n")

        rendered = "".join(candidate)
        if len(rendered.encode("utf-8")) <= budget:
            return rendered

        if start < line_index:
            start += 1

            continue

        if end > line_index + 1:
            end -= 1

            continue

        return None


def focused_diff_excerpt(
    root: Path,
    path: str,
    maximum: int,
    anchors: tuple[str, ...] | None = None,
) -> str:
    """
    Return bounded changed hunks relevant to the shell/help pilot review.

    Parameters
    ----------
    root : Path
        Repository root used for the unstaged Git diff.
    path : str
        Repository-relative target path.
    maximum : int
        Maximum UTF-8 byte length of the returned complete-line excerpt.
    anchors : tuple[str, ...] | None
        Optional markers used to select relevant complete hunk windows.

    Returns
    -------
    excerpt : str
        Complete bounded diff windows, or an empty string when none qualify.
    """

    result = run_command(
        ["git", "diff", "--no-ext-diff", "--unified=3", "--", path],
        root,
        dict(os.environ, PYTHONDONTWRITEBYTECODE="1"),
        30,
        65536,
    )

    if result.returncode or not result.stdout:
        return ""

    relevant_anchors = anchors or (
        "dir_scr",
        "parse_args",
        "build_cmd",
        "env_nam",
        "handle_env",
        "source_helpers",
        "help_",
    )

    hunks = re.split(r"(?=^@@ )", result.stdout, flags=re.MULTILINE)
    header, body = hunks[0], hunks[1:]
    selected = [
        hunk
        for hunk in body
        if any(anchor in hunk for anchor in relevant_anchors)
    ] or body[:2]
    kept = header

    for hunk in selected:
        if anchors is not None or len((kept + hunk).encode("utf-8")) > maximum:
            lines = hunk.splitlines(keepends=True)
            indexes = [
                index
                for index, line in enumerate(lines)
                if line.startswith(("+", "-"))
                and any(anchor in line for anchor in relevant_anchors)
            ]

            if not indexes:
                indexes = [
                    index
                    for index, line in enumerate(lines)
                    if line.startswith(("+", "-"))
                ]

            rendered_windows: set[str] = set()
            omitted_windows = False

            for index in indexes:
                remaining = maximum - len(kept.encode("utf-8"))
                excerpt = _complete_line_diff_window(hunk, index, remaining)

                if excerpt is None:
                    omitted_windows = True
                    continue

                if excerpt in rendered_windows:
                    continue

                kept += excerpt
                rendered_windows.add(excerpt)
                if anchors is not None:
                    break

            if not rendered_windows:
                kept += (
                    "[relevant hunk omitted: no complete-line window fits "
                    "excerpt budget]\n"
                )
            elif omitted_windows or len(rendered_windows) < len(set(indexes)):
                kept += (
                    "[additional relevant windows omitted from this hunk]\n"
                )

            continue

        kept += hunk

    return kept


def source_anchor_window(
    source: str,
    markers: tuple[str, ...],
    maximum: int = 4000,
) -> str:
    """
    Return a complete Bash function or local source window containing markers.
    """

    lines = source.splitlines(keepends=True)
    function_starts = [
        index
        for index, line in enumerate(lines)
        if re.match(r"^function\s+[A-Za-z0-9_]+\(\)\s*\{", line)
    ]

    for position, start in enumerate(function_starts):
        end = (
            function_starts[position + 1]
            if position + 1 < len(function_starts)
            else len(lines)
        )
        candidate = "".join(lines[start:end])
        if (
            all(marker in candidate for marker in markers)
            and len(candidate.encode("utf-8")) <= maximum
        ):
            return candidate

    marker_lines = [
        index
        for index, line in enumerate(lines)
        if any(marker in line for marker in markers)
    ]
    if not marker_lines:
        return ""

    first, last = min(marker_lines), max(marker_lines)
    candidate = "".join(lines[max(0, first - 3) : min(len(lines), last + 4)])
    if len(candidate.encode("utf-8")) > maximum:
        return ""

    return candidate


def anchor_evidence_windows(
    root: Path,
    target_content: dict[str, str],
) -> dict[str, dict[str, dict[str, str]]]:
    """
    Build per-anchor complete-line evidence, preferring changed diff windows.
    """

    evidence: dict[str, dict[str, dict[str, str]]] = {}

    for path, groups in REQUIRED_BEHAVIORAL_ANCHORS.items():
        source = target_content.get(path, "")
        if not source:
            raise ValueError(
                f"missing source for behavioral anchor path: {path}",
            )

        evidence[path] = {}

        for group, markers in groups.items():
            diff = (
                ""
                if (path, group) in SOURCE_ONLY_ANCHORS
                else focused_diff_excerpt(root, path, 4000, markers)
            )

            if diff and all(marker in diff for marker in markers):
                rendered, kind = diff, "complete_line_diff_window"
            else:
                rendered, kind = (
                    source_anchor_window(source, markers),
                    "complete_line_source_window",
                )

            if not rendered or not all(
                marker in rendered for marker in markers
            ):
                raise ValueError(
                    (
                        f"missing required behavioral anchor window {group}: "
                        f"{path}"
                    ),
                )

            evidence[path][group] = {
                "evidence_kind": kind,
                "rendered_evidence": rendered,
            }

    return evidence


def markdown_fences_balanced(markdown: str) -> bool:
    """
    Return whether nested Markdown fences have matching delimiters.
    """

    active: tuple[str, int] | None = None

    for line in markdown.splitlines():
        match = re.match(r"^\s*(`{3,}|~{3,})", line)
        if not match:
            continue

        marker = match.group(1)

        if active is None:
            active = (marker[0], len(marker))
        elif marker[0] == active[0] and len(marker) >= active[1]:
            active = None

    return active is None


def bounded_markdown_provision(text: str, heading: str, maximum: int) -> str:
    """
    Extract a bounded Markdown provision without leaving an inner fence open.

    Parameters
    ----------
    text : str
        Complete Markdown document.
    heading : str
        Exact heading text whose section is selected.
    maximum : int
        Maximum character count for the returned provision.

    Returns
    -------
    provision : str
        Complete bounded section text, or an empty string if it cannot fit.

    Raises
    ------
    ValueError
        If the selected provision contains an unmatched fenced block.
    """

    lines = text.splitlines(keepends=True)
    start = next(
        (
            index
            for index, line in enumerate(lines)
            if re.fullmatch(
                r"#{1,6}\s+" + re.escape(heading) + r"\s*",
                line.rstrip(),
            )
        ),
        None,
    )
    if start is None:
        return ""

    level = len(lines[start]) - len(lines[start].lstrip("#"))
    end = next(
        (
            index
            for index in range(start + 1, len(lines))
            if re.match(r"^#{1," + str(level) + r"}\s", lines[index])
            and len(lines[index]) - len(lines[index].lstrip("#")) <= level
        ),
        len(lines),
    )

    blocks: list[str] = []
    index = start

    while index < end:
        block = [lines[index]]
        index += 1

        if re.match(r"^\s*(`{3,}|~{3,})", block[0]):
            marker = re.match(r"^\s*(`{3,}|~{3,})", block[0]).group(1)

            while index < end:
                block.append(lines[index])

                if re.match(
                    r"^\s*"
                    + re.escape(marker[0])
                    + r"{"
                    + str(len(marker))
                    + r",}\s*$",
                    lines[index],
                ):
                    index += 1

                    break

                index += 1
        else:
            while index < end and lines[index].strip():
                block.append(lines[index])
                index += 1

        blocks.append("".join(block))

        while index < end and not lines[index].strip():
            blocks.append(lines[index])
            index += 1

    rendered = ""

    for block in blocks:
        candidate = rendered + block
        if len(candidate.encode("utf-8")) > maximum:
            break

        rendered = candidate

    rendered = rendered.rstrip() + "\n"
    if not markdown_fences_balanced(rendered):
        raise ValueError(
            "bounded Markdown provision contains an unmatched fence",
        )

    return rendered


def read_ndjson(path: Path) -> list[dict[str, object]]:
    """
    Read a required NDJSON report with a precise malformed error.

    Parameters
    ----------
    path : Path
        Required report path whose nonblank lines must be JSON objects.

    Returns
    -------
    rows : list[dict[str, object]]
        Parsed report objects in physical source order.

    Raises
    ------
    IncompleteReportError
        If the required report does not exist.
    ValueError
        If any line is invalid JSON or does not decode to an object.
    """

    if not path.is_file():
        raise IncompleteReportError(f"missing required report: {path.name}")

    rows: list[dict[str, object]] = []

    for number, line in enumerate(
        path.read_text(encoding="utf-8").splitlines(),
        1,
    ):
        try:
            row = json.loads(line)
        except json.JSONDecodeError as exc:
            raise ValueError(
                f"malformed {path.name} line {number}: {exc.msg}",
            ) from exc

        if not isinstance(row, dict):
            raise ValueError(
                f"malformed {path.name} line {number}: expected object",
            )

        rows.append(row)

    return rows


def controlled_smoke_target(
    selection: dict[str, object] | None,
) -> dict[str, object]:
    """
    Return the one manifest-declared controlled-smoke target.
    """

    if selection is None:
        raise ValueError("controlled smoke evidence requires --paths-from")

    matches = [
        item
        for item in selection["targets"]
        if isinstance(item, dict)
        and isinstance(item.get("evidence"), dict)
        and item["evidence"].get("kind") == CONTROLLED_SMOKE_KIND
    ]
    if len(matches) != 1:
        raise ValueError(
            "pilot selection must declare exactly one controlled smoke target",
        )

    return matches[0]


def controlled_smoke_sources(
    root: Path,
    target: dict[str, object],
) -> list[dict[str, str]]:
    """
    Fingerprint the declared test and production dependencies in stable order.
    """

    evidence = target["evidence"]
    paths = [
        str(target["path"]),
        *[str(path) for path in evidence["source_dependencies"]],
    ]
    rows: list[dict[str, str]] = []

    for path in sorted(paths):
        candidate = root / path
        if not candidate.is_file() or candidate.is_symlink():
            raise ValueError(
                f"controlled smoke source is not a regular file: {path}",
            )

        fingerprint, source = entry_fingerprint(root, path)
        rows.append(
            {
                "path": path,
                "content_fingerprint": fingerprint,
                "content_source": source,
            },
        )

    return rows


def controlled_smoke_source_hash(rows: list[dict[str, str]]) -> str:
    """
    Return the canonical fingerprint-set hash for controlled smoke evidence.
    """

    return sha256_bytes(json_bytes(rows))


def controlled_smoke_command(target: dict[str, object]) -> dict[str, object]:
    """
    Return the fixed outer-shell command identity for one controlled target.
    """

    command = target["evidence"]["command"]
    argv = [
        str(part).replace("{path}", str(target["path"]))
        for part in command["argv"]
    ]
    environment = {
        str(key): str(value) for key, value in command["environment"].items()
    }
    display = " ".join(
        [
            *(f"{key}={value}" for key, value in sorted(environment.items())),
            *argv,
        ],
    )

    return {
        "environment": environment,
        "argv": argv,
        "display": display,
        "identity": sha256_bytes(
            json_bytes({"environment": environment, "argv": argv}),
        ),
    }


def outer_smoke_context() -> dict[str, object]:
    """
    Record host-shell facts without changing the command environment.
    """

    conda_environment = os.environ.get("CONDA_DEFAULT_ENV")
    return {
        "path_fingerprint": sha256_bytes(
            os.environ.get("PATH", "").encode("utf-8"),
        ),
        "wget_path": shutil.which("wget"),
        "wget_available": shutil.which("wget") is not None,
        "conda_path": shutil.which("conda"),
        "conda_available": shutil.which("conda") is not None,
        "conda_default_env": conda_environment,
        "project_environment_active": bool(
            conda_environment and conda_environment != "base",
        ),
    }


def capture_controlled_smoke(
    args: argparse.Namespace,
    public_root: Path,
    selection: dict[str, object] | None,
) -> int:
    """
    Capture the exact manifest command from the unreplaced outer environment.
    """

    if not args.run_id:
        raise ValueError("controlled smoke capture requires --run-id")

    target = controlled_smoke_target(selection)
    evidence = target["evidence"]
    command = controlled_smoke_command(target)
    before = controlled_smoke_sources(public_root, target)

    outer = outer_smoke_context()
    requirements = evidence["outer_context"]
    context_ok = (
        (
            not requirements["require_wget_absent"]
            or not outer["wget_available"]
        )
        and (
            not requirements["require_conda_available"]
            or outer["conda_available"]
        )
        and (
            not requirements["require_inactive_project_environment"]
            or not outer["project_environment_active"]
        )
    )

    result = CommandResult(command["argv"], "", "", None, 0.0)

    if context_ok:
        result = run_command(
            command["argv"],
            public_root,
            dict(os.environ, **command["environment"]),
        )

    after = controlled_smoke_sources(public_root, target)
    raw = {"stdout": result.stdout, "stderr": result.stderr}

    capture = {
        "schema_version": 1,
        "run_id": args.run_id,
        "manifest_path": str(selection["selection_path"]),
        "manifest_fingerprint": entry_fingerprint(
            public_root,
            str(selection["selection_path"]),
        )[0],
        "target_path": target["path"],
        "command": command,
        "outer_context": outer,
        "outer_context_ok": context_ok,
        "source_fingerprints_before": before,
        "source_fingerprints_after": after,
        "raw": {**raw, "digest": sha256_bytes(json_bytes(raw))},
        "result": {
            "exit_status": result.returncode,
            "timed_out": result.timed_out,
            "launch_error": result.launch_error,
        },
    }

    write_json(Path(args.capture_controlled_smoke), capture)
    return 0


def normalize_controlled_smoke(
    args: argparse.Namespace,
    public_root: Path,
    selection: dict[str, object] | None,
) -> int:
    """
    Normalize one captured outer-shell test result into fail-closed evidence.

    Parameters
    ----------
    args : argparse.Namespace
        Paths and run identity for captured and normalized evidence.
    public_root : Path
        Repository root used to verify selected source fingerprints.
    selection : dict[str, object] | None
        Validated pilot selection containing the controlled-smoke contract.

    Returns
    -------
    status : int
        Zero for validated passing evidence, or the stable nonpassing status.

    Raises
    ------
    ValueError
        If required run, selection, or controlled-smoke metadata is absent.
    """

    if not args.run_id:
        raise ValueError("controlled smoke normalization requires --run-id")

    target = controlled_smoke_target(selection)
    evidence = target["evidence"]

    try:
        loaded_capture = json.loads(
            Path(args.normalize_controlled_smoke).read_text(encoding="utf-8"),
        )
    except (OSError, json.JSONDecodeError):
        loaded_capture = {}

    capture = loaded_capture if isinstance(loaded_capture, dict) else {}
    command = controlled_smoke_command(target)
    sources = controlled_smoke_sources(public_root, target)

    raw = capture.get("raw", {})
    raw = raw if isinstance(raw, dict) else {}
    stdout, stderr = raw.get("stdout"), raw.get("stderr")
    combined = "\n".join(
        item for item in (stdout, stderr) if isinstance(item, str)
    )

    summaries = list(CONTROLLED_SMOKE_SUMMARY.finditer(combined))
    summary = summaries[-1] if summaries else None
    counts = (
        {
            name: int(summary.group(name))
            for name in ("pass", "fail", "warn", "skip")
        }
        if summary
        else None
    )

    groups = {
        name: {
            "required_messages": messages,
            "passed": all(
                f"PASS: {message}" in combined for message in messages
            ),
        }
        for name, messages in evidence["required_pass_groups"].items()
    }

    source_before = capture.get("source_fingerprints_before")
    source_after = capture.get("source_fingerprints_after")

    result = capture.get("result", {})
    result = result if isinstance(result, dict) else {}

    observations = {
        name: {
            "value": observation["value"],
            "evidence": {
                "kind": "pass_message",
                "message": observation["pass_message"],
            },
        }
        for name, observation in evidence.get(
            "optional_observations",
            {},
        ).items()
        if f"PASS: {observation['pass_message']}" in combined.splitlines()
    }

    expected_manifest = entry_fingerprint(
        public_root,
        str(selection["selection_path"]),
    )[0]
    expected_raw_digest = sha256_bytes(
        json_bytes({"stdout": stdout, "stderr": stderr}),
    )
    checks = {
        "run_id": capture.get("run_id") == args.run_id,
        "manifest": capture.get("manifest_fingerprint") == expected_manifest,
        "target": capture.get("target_path") == target["path"],
        "command": capture.get("command") == command,
        "raw_digest": raw.get("digest") == expected_raw_digest,
        "outer_context": capture.get("outer_context_ok") is True,
        "source_stability": source_before == source_after == sources,
        "summary": summary is not None
        and summary.group("test") == evidence["test_identifier"],
        "exit_status": result.get("exit_status") == 0,
        "counts": counts is not None
        and counts["pass"] > 0
        and counts["fail"] == 0
        and counts["skip"] == 0,
        "pass_groups": all(item["passed"] for item in groups.values()),
    }
    passed = all(checks.values())

    normalized = {
        "schema_version": 1,
        "run_id": args.run_id,
        "status": "passed" if passed else "invalid",
        "manifest_path": str(selection["selection_path"]),
        "manifest_fingerprint": entry_fingerprint(
            public_root,
            str(selection["selection_path"]),
        )[0],
        "target_path": target["path"],
        "command": command,
        "summary": {
            "test_identifier": summary.group("test") if summary else None,
            "counts": counts,
        },
        "exit_status": result.get("exit_status"),
        "pass_groups": groups,
        "source_fingerprints": sources,
        "source_fingerprint_set_hash": controlled_smoke_source_hash(sources),
        "raw_capture_digest": raw.get("digest"),
        "environment_selection": {"submit_worker": "not_directly_evidenced"},
        "worker_execution": {
            "status": "passed" if passed else "not_established",
            "environment_evidence": "not_directly_evidenced",
        },
        "limitation": evidence["limitation"],
        "validation": checks,
    }

    if observations:
        normalized["observations"] = observations

    write_json(Path(args.runtime_evidence_out), normalized)
    return 0 if passed else EXIT_VERIFY_NONCOMPLETED


def validate_runtime_evidence(
    value: object,
    public_root: Path,
    selection: dict[str, object],
    run_id: str,
) -> dict[str, object]:
    """
    Validate a normalized controlled-smoke record against current sources.

    Parameters
    ----------
    value : object
        Candidate normalized runtime-evidence record.
    public_root : Path
        Repository root used for current manifest and source fingerprints.
    selection : dict[str, object]
        Validated pilot selection and controlled-smoke requirements.
    run_id : str
        Required run identity.

    Returns
    -------
    evidence : dict[str, object]
        Validated runtime-evidence record.

    Raises
    ------
    ValueError
        If identity, command, sources, observations, or pass proof is stale or
        malformed.
    """

    target = controlled_smoke_target(selection)
    expected_sources = controlled_smoke_sources(public_root, target)
    expected_command = controlled_smoke_command(target)

    if (
        not isinstance(value, dict)
        or value.get("schema_version") != 1
        or value.get("status") != "passed"
    ):
        raise ValueError(
            "runtime evidence is not a passed schema_version 1 record",
        )

    if (
        value.get("run_id") != run_id
        or value.get("target_path") != target["path"]
    ):
        raise ValueError("runtime evidence run ID or target does not match")

    expected_manifest = entry_fingerprint(
        public_root,
        str(selection["selection_path"]),
    )[0]
    if (
        value.get("manifest_path") != selection["selection_path"]
        or value.get("manifest_fingerprint") != expected_manifest
    ):
        raise ValueError("runtime evidence manifest is stale")

    expected_source_hash = controlled_smoke_source_hash(expected_sources)
    if (
        value.get("command") != expected_command
        or value.get("source_fingerprints") != expected_sources
        or value.get("source_fingerprint_set_hash") != expected_source_hash
    ):
        raise ValueError("runtime evidence command or sources are stale")

    summary = value.get("summary")
    groups = value.get("pass_groups")
    validation = value.get("validation")

    expected_test_identifier = target["evidence"]["test_identifier"]
    if (
        not isinstance(summary, dict)
        or summary.get("test_identifier") != expected_test_identifier
        or not isinstance(summary.get("counts"), dict)
    ):
        raise ValueError("runtime evidence summary is malformed")

    counts = summary["counts"]
    if (
        counts.get("pass", 0) <= 0
        or counts.get("fail") != 0
        or counts.get("skip") != 0
        or value.get("exit_status") != 0
    ):
        raise ValueError(
            "runtime evidence is not a successful controlled result",
        )

    if (
        not isinstance(groups, dict)
        or set(groups) != set(target["evidence"]["required_pass_groups"])
        or not all(
            isinstance(item, dict) and item.get("passed") is True
            for item in groups.values()
        )
    ):
        raise ValueError("runtime evidence is missing required pass groups")

    if (
        not isinstance(validation, dict)
        or not validation
        or not all(item is True for item in validation.values())
    ):
        raise ValueError(
            "runtime evidence did not pass fail-closed validation",
        )

    if (
        not isinstance(value.get("limitation"), str)
        or value["limitation"] != target["evidence"]["limitation"]
    ):
        raise ValueError("runtime evidence limitation does not match manifest")

    if value.get("environment_selection") != {
        "submit_worker": "not_directly_evidenced",
    }:
        raise ValueError(
            "runtime evidence overstates submit-worker environment selection",
        )

    if value.get("worker_execution") != {
        "status": "passed",
        "environment_evidence": "not_directly_evidenced",
    }:
        raise ValueError("runtime evidence worker result is malformed")

    observations = value.get("observations", {})
    configured = target["evidence"].get("optional_observations", {})
    if not isinstance(observations, dict) or not set(observations).issubset(
        configured,
    ):
        raise ValueError("runtime evidence has undeclared observations")

    for name, observation in observations.items():
        specification = configured[name]
        expected = {
            "value": specification["value"],
            "evidence": {
                "kind": "pass_message",
                "message": specification["pass_message"],
            },
        }
        if observation != expected:
            raise ValueError(
                "runtime evidence observation does not match manifest",
            )

    return value


def validate_pilot_artifact_rows(
    name: str,
    rows: list[dict[str, object]],
) -> None:
    """
    Validate version and minimum shape for one versioned pilot record family.
    """

    required = PILOT_ARTIFACT_REQUIRED_FIELDS[name]

    for index, row in enumerate(rows, 1):
        if row.get(
            "schema_version",
        ) != PILOT_ARTIFACT_SCHEMA_VERSION or required - set(row):
            raise ValueError(f"malformed {name} line {index}")


def verify_runtime_fact_rows(
    rows: list[dict[str, object]],
    public_root: Path,
) -> list[dict[str, object]]:
    """
    Return stale source records for normalized controlled-smoke facts.

    Parameters
    ----------
    rows : list[dict[str, object]]
        Fact rows that may contain controlled-smoke execution evidence.
    public_root : Path
        Repository root used to verify recorded source fingerprints.

    Returns
    -------
    rows : list[dict[str, object]]
        Stale source-fingerprint rows with recorded and current values.

    Raises
    ------
    ValueError
        If a controlled-smoke fact lacks successful fail-closed proof.
    """

    stale: list[dict[str, object]] = []

    for row in rows:
        if row.get("topic") != "controlled_smoke_execution_result":
            continue

        value = row.get("value")
        if (
            not isinstance(value, dict)
            or value.get("status") != "passed"
            or value.get("exit_status") != 0
        ):
            raise ValueError("malformed controlled-smoke runtime fact")

        summary = value.get("summary")
        groups = value.get("pass_groups")
        validation = value.get("validation")

        if (
            not isinstance(summary, dict)
            or not isinstance(summary.get("counts"), dict)
            or summary["counts"].get("pass", 0) <= 0
            or summary["counts"].get("fail") != 0
            or summary["counts"].get("skip") != 0
        ):
            raise ValueError("controlled-smoke runtime fact is not successful")

        if (
            not isinstance(groups, dict)
            or not groups
            or not all(
                isinstance(item, dict) and item.get("passed") is True
                for item in groups.values()
            )
        ):
            raise ValueError(
                "controlled-smoke runtime fact is missing pass groups",
            )

        if (
            not isinstance(validation, dict)
            or not validation
            or not all(item is True for item in validation.values())
        ):
            raise ValueError("controlled-smoke runtime fact failed validation")

        if value.get("environment_selection") != {
            "submit_worker": "not_directly_evidenced",
        }:
            raise ValueError(
                (
                    "controlled-smoke runtime fact overstates worker "
                    "environment selection"
                ),
            )

        sources = value.get("source_fingerprints")
        if not isinstance(sources, list) or not all(
            isinstance(item, dict)
            and isinstance(item.get("path"), str)
            and isinstance(item.get("content_fingerprint"), str)
            for item in sources
        ):
            raise ValueError(
                "controlled-smoke runtime fact has malformed sources",
            )

        normalized_sources = [
            {
                "path": str(item["path"]),
                "content_fingerprint": str(item["content_fingerprint"]),
                "content_source": str(item.get("content_source", "")),
            }
            for item in sources
        ]
        if value.get(
            "source_fingerprint_set_hash",
        ) != controlled_smoke_source_hash(normalized_sources):
            raise ValueError(
                "controlled-smoke runtime fact source hash mismatch",
            )

        for source in normalized_sources:
            actual, _ = entry_fingerprint(public_root, source["path"])

            if actual != source["content_fingerprint"]:
                stale.append(
                    {
                        "record_type": "controlled_smoke_runtime",
                        "path": source["path"],
                        "recorded": source["content_fingerprint"],
                        "actual": actual,
                        "status": "stale",
                    },
                )

    return stale


def _load_run_metadata(report: Path) -> dict[str, object]:
    """
    Load and validate the report run metadata.

    Parameters
    ----------
    report : Path
        Report directory under verification.

    Returns
    -------
    metadata : dict[str, object]
        Validated run metadata.

    Raises
    ------
    IncompleteReportError
        If the required metadata file is absent.
    ValueError
        If the run status is malformed.
    """

    run_path = report / "run.json"
    if not run_path.is_file():
        raise IncompleteReportError("missing required report: run.json")

    run = json.loads(run_path.read_text(encoding="utf-8"))
    allowed_statuses = {
        "completed",
        "aborted",
        "error",
    }

    run_status_valid = (
        isinstance(run, dict) and run.get("status") in allowed_statuses
    )

    if not run_status_valid:
        raise ValueError(
            "malformed run.json: status must be completed, aborted, or error",
        )

    return run


def _require_report_artifacts(
    report: Path,
    run: dict[str, object],
) -> None:
    """
    Require every artifact declared by the report format.

    Parameters
    ----------
    report : Path
        Report directory under verification.
    run : dict[str, object]
        Validated run metadata.

    Raises
    ------
    IncompleteReportError
        If a required report artifact is absent.
    """

    required = [
        "ledger.ndjson",
        "artifacts.ndjson",
        "findings.ndjson",
        "checks.ndjson",
        "prompts.ndjson",
        "cohort_progress.json",
        "rule_manifest.json",
    ]

    if run.get("report_format_version") == 2:
        required.extend(
            [
                "targets.ndjson",
                "facts.ndjson",
                "policy_questions.ndjson",
                "adapter_limitations.ndjson",
                "pilot_report.json",
                "semantic_review.ndjson",
            ],
        )

    for name in required:
        if not (report / name).is_file():
            raise IncompleteReportError(f"missing required report: {name}")


def _verify_report_source_fingerprints(
    run: dict[str, object],
    ledger: list[dict[str, object]],
    findings: list[dict[str, object]],
    prompts: list[dict[str, object]],
    public_root: Path,
    private_root: Path,
    stale: list[dict[str, object]],
) -> None:
    """
    Compare primary report source fingerprints with the current repositories.

    Parameters
    ----------
    run : dict[str, object]
        Validated run metadata.
    ledger : list[dict[str, object]]
        Report ledger rows.
    findings : list[dict[str, object]]
        Report finding rows.
    prompts : list[dict[str, object]]
        Report prompt rows.
    public_root : Path
        Public repository root.
    private_root : Path
        Private evidence repository root.
    stale : list[dict[str, object]]
        Mutable staleness records.

    Raises
    ------
    ValueError
        If fingerprint metadata is malformed.
    """

    fingerprint_checks = [
        (
            "baseline",
            run.get("baseline_path"),
            run.get("baseline_fingerprint"),
            public_root,
        ),
        (
            "rule_manifest",
            run.get("rule_manifest_path"),
            run.get("rule_manifest_fingerprint"),
            public_root,
        ),
    ]

    for record_type, path, recorded, root in fingerprint_checks:
        if not isinstance(path, str) or not isinstance(recorded, str):
            raise ValueError(
                (
                    f"malformed run.json: missing {record_type} "
                    f"fingerprint metadata"
                ),
            )

        actual, _ = entry_fingerprint(root, path)

        if actual != recorded:
            stale.append(
                {
                    "record_type": record_type,
                    "path": path,
                    "recorded": recorded,
                    "actual": actual,
                    "status": "stale",
                },
            )

    for evidence in run.get("consumed_evidence", []):
        if (
            not isinstance(evidence, dict)
            or not isinstance(evidence.get("path"), str)
            or not isinstance(evidence.get("fingerprint"), str)
        ):
            raise ValueError(
                "malformed run.json: invalid consumed_evidence record",
            )

        root = (
            public_root
            if evidence.get("repository") == "public"
            else private_root
            if evidence.get("repository") == "private"
            else None
        )
        if root is None:
            raise ValueError(
                (
                    "malformed run.json: consumed evidence has unknown "
                    "repository"
                ),
            )

        actual, _ = entry_fingerprint(root, evidence["path"])

        if actual != evidence["fingerprint"]:
            stale.append(
                {
                    "record_type": "consumed_evidence",
                    "repository": evidence["repository"],
                    "path": evidence["path"],
                    "recorded": evidence["fingerprint"],
                    "actual": actual,
                    "status": "stale",
                },
            )

    for record_type, rows in (("ledger", ledger), ("finding", findings)):
        for row in rows:
            path = row.get("path")
            recorded = row.get("content_fingerprint")
            if path == ".":
                continue

            if not isinstance(path, str) or not isinstance(recorded, str):
                raise ValueError(
                    f"malformed {record_type} record fingerprint",
                )

            actual, _ = entry_fingerprint(public_root, path)

            if actual != recorded:
                stale.append(
                    {
                        "record_type": record_type,
                        "path": path,
                        "recorded": recorded,
                        "actual": actual,
                        "status": "stale",
                    },
                )

    for prompt in prompts:
        targets = prompt.get("target_fingerprints")
        if not isinstance(targets, dict):
            raise ValueError("malformed prompt record target_fingerprints")

        for path, recorded in targets.items():
            actual, _ = entry_fingerprint(public_root, str(path))

            if actual != recorded:
                stale.append(
                    {
                        "record_type": "prompt",
                        "batch_id": prompt.get("batch_id"),
                        "path": path,
                        "recorded": recorded,
                        "actual": actual,
                        "status": "stale",
                    },
                )


def _verify_versioned_report(
    report: Path,
    public_root: Path,
    stale: list[dict[str, object]],
) -> None:
    """
    Validate version-two pilot artifacts and their provenance.

    Parameters
    ----------
    report : Path
        Report directory under verification.
    public_root : Path
        Public repository root.
    stale : list[dict[str, object]]
        Mutable staleness records.

    Raises
    ------
    ValueError
        If versioned artifacts or provenance are malformed.
    """

    targets = read_ndjson(report / "targets.ndjson")
    validate_pilot_artifact_rows("targets.ndjson", targets)

    for target in targets:
        if not isinstance(target.get("path"), str) or not isinstance(
            target.get("content_fingerprint"),
            str,
        ):
            raise ValueError("malformed targets.ndjson record")

        actual, _ = entry_fingerprint(public_root, str(target["path"]))

        if actual != target["content_fingerprint"]:
            stale.append(
                {
                    "record_type": "target",
                    "path": target["path"],
                    "recorded": target["content_fingerprint"],
                    "actual": actual,
                    "status": "stale",
                },
            )

    versioned_rows = {
        name: read_ndjson(report / name)
        for name in (
            "facts.ndjson",
            "policy_questions.ndjson",
            "adapter_limitations.ndjson",
            "semantic_review.ndjson",
        )
    }

    for name, rows in versioned_rows.items():
        validate_pilot_artifact_rows(name, rows)

    pilot_report = json.loads(
        (report / "pilot_report.json").read_text(encoding="utf-8"),
    )
    pilot_schema_version = (
        pilot_report.get("schema_version")
        if isinstance(pilot_report, dict)
        else None
    )
    if (
        not isinstance(pilot_report, dict)
        or pilot_schema_version != PILOT_ARTIFACT_SCHEMA_VERSION
    ):
        raise ValueError("malformed pilot_report.json schema_version")

    provenance = pilot_report.get("package_provenance")

    if provenance is not None:
        hashes = (
            provenance.get("artifact_hashes", {})
            if isinstance(provenance, dict)
            else {}
        )
        supplied = (
            provenance.get("supplied_artifacts", [])
            if isinstance(provenance, dict)
            else []
        )
        if (
            not isinstance(supplied, list)
            or len(supplied) != 5
            or "pilot_report.json" not in supplied
        ):
            raise ValueError(
                "package provenance must declare five supplied artifacts",
            )

        for name in supplied:
            if name == "pilot_report.json":
                continue

            validate_target_path(name)

            if hashes.get(name) != file_sha256(report / name):
                raise ValueError(f"package artifact hash mismatch: {name}")

        self_hash = (
            provenance.get("pilot_report_self_hash", {})
            if isinstance(provenance, dict)
            else {}
        )
        canonicalization = self_hash.get("canonicalization")
        excluded_pointer = self_hash.get("excluded_json_pointer")
        expected_canonicalization = "utf8_json_sorted_keys_compact_separators"
        expected_excluded_pointer = (
            "/package_provenance/pilot_report_self_hash/value"
        )
        expected_self_hash = pilot_report_self_hash(pilot_report)
        if (
            self_hash.get("algorithm") != "sha256"
            or canonicalization != expected_canonicalization
            or excluded_pointer != expected_excluded_pointer
            or self_hash.get("value") != expected_self_hash
        ):
            raise ValueError(
                "package artifact hash mismatch: pilot_report.json",
            )

    runtime_contract = pilot_report.get("runtime_evidence")

    if (
        isinstance(runtime_contract, dict)
        and runtime_contract.get("required") is True
    ):
        runtime_rows = [
            row
            for row in versioned_rows["facts.ndjson"]
            if row.get("topic") == "controlled_smoke_execution_result"
        ]
        if len(runtime_rows) != 1:
            raise ValueError("missing controlled-smoke runtime fact")

        stale.extend(
            verify_runtime_fact_rows(runtime_rows, public_root),
        )


def _verification_outcome(
    status: object,
    stale: list[dict[str, object]],
    allow_partial: bool,
) -> tuple[str, int]:
    """
    Return the public verification status and stable exit code.

    Parameters
    ----------
    status : object
        Validated run status.
    stale : list[dict[str, object]]
        Current staleness records.
    allow_partial : bool
        Whether one fresh aborted report may be accepted.

    Returns
    -------
    status_text, status : tuple[str, int]
        Public status spelling and stable process result.
    """

    if status == "completed":
        return (
            "completed_stale" if stale else "completed_fresh",
            EXIT_VERIFY_STALE if stale else 0,
        )

    if status == "aborted":
        code = 0 if allow_partial and not stale else EXIT_VERIFY_NONCOMPLETED

        return "aborted_stale" if stale else "aborted_fresh", code

    return "error_or_incomplete_report", EXIT_VERIFY_NONCOMPLETED


def verify_report(args: argparse.Namespace) -> int:
    """
    Verify report status and every source fingerprint that determined it.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments or an explicit argument sequence.

    Returns
    -------
    status : int
        Stable report-verification exit status.

    Raises
    ------
    ValueError
        If a supplied value violates the validated contract.
    IncompleteReportError
        If the validated operation cannot be completed.
    """

    read_only = bool(getattr(args, "verify_read_only", None))
    report = Path(
        args.verify_read_only if read_only else args.verify,
    ).resolve()
    summary: dict[str, object] = {
        "report_directory": str(report),
        "allow_partial": args.allow_partial,
        "status": "malformed_report",
        "staleness_count": 0,
    }
    stale: list[dict[str, object]] = []

    try:
        run = _load_run_metadata(report)
        _require_report_artifacts(report, run)

        ledger = read_ndjson(report / "ledger.ndjson")
        findings = read_ndjson(report / "findings.ndjson")
        prompts = read_ndjson(report / "prompts.ndjson")

        public_root = Path(args.public_root).resolve()
        private_root = Path(args.private_root).resolve()

        _verify_report_source_fingerprints(
            run,
            ledger,
            findings,
            prompts,
            public_root,
            private_root,
            stale,
        )

        if run.get("report_format_version") == 2:
            _verify_versioned_report(report, public_root, stale)

        summary["status"], code = _verification_outcome(
            run["status"],
            stale,
            args.allow_partial,
        )
    except IncompleteReportError as exc:
        summary["status"] = "error_or_incomplete_report"
        summary["error"] = f"{type(exc).__name__}: {exc}"
        code = EXIT_VERIFY_NONCOMPLETED
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        summary["status"] = "malformed_report"
        summary["error"] = f"{type(exc).__name__}: {exc}"
        code = EXIT_VERIFY_MALFORMED

    summary["staleness_count"] = len(stale)

    if not read_only:
        write_ndjson(report / "staleness.ndjson", stale)
        write_json(report / "verification.json", summary)

        if code == 0 and (report / "pilot_report.json").is_file():
            pilot_report = json.loads(
                (report / "pilot_report.json").read_text(encoding="utf-8"),
            )

            if isinstance(pilot_report, dict) and isinstance(
                pilot_report.get("package_provenance"),
                dict,
            ):
                pilot_report["verification"]["fresh_verification"] = str(
                    summary["status"],
                )
                pilot_report["package_provenance"]["pilot_report_self_hash"][
                    "value"
                ] = pilot_report_self_hash(pilot_report)
                write_json(report / "pilot_report.json", pilot_report)
    elif code:
        message = (
            f"read-only verification failed: {summary['status']} "
            f"({summary['staleness_count']} stale records): {stale[:3]}"
        )
        print(
            message,
            file=sys.stderr,
        )

    return code


def validate_linked_pair_headers(
    pilots: list[dict[str, object]],
) -> tuple[dict[str, object], dict[str, object], list[str]]:
    """
    Validate shared and reciprocal identities for two linked pilot reports.

    Parameters
    ----------
    pilots : list[dict[str, object]]
        Exactly two child pilot-report payloads.

    Returns
    -------
    left_link, right_link, child_ids : tuple[
        dict[str, object], dict[str, object], list[str]
    ]
        Left linkage, right linkage, and canonical ordered child identities.

    Raises
    ------
    ValueError
        If linkage, ownership fingerprints, identities, or ordering disagree.
    """

    if len(pilots) != 2 or any(
        not isinstance(pilot.get("linked_package"), dict) for pilot in pilots
    ):
        raise ValueError("linked child report is missing linkage metadata")

    left, right = pilots
    left_link = left["linked_package"]
    right_link = right["linked_package"]
    shared_keys = {
        "bundle_id",
        "selection_path",
        "config_fingerprint",
        "target_ownership_fingerprint",
        "graph_ownership_fingerprint",
        "umbrella_target_ownership",
        "umbrella_context_ownership",
        "umbrella_dependency_edges",
        "umbrella_target_count",
        "umbrella_context_count",
        "umbrella_edge_count",
        "child_ids",
    }
    if any(left_link.get(key) != right_link.get(key) for key in shared_keys):
        raise ValueError("linked child umbrella metadata mismatch")

    child_ids = list(left_link["child_ids"])
    actual_children = [left.get("package_id"), right.get("package_id")]
    if actual_children != child_ids or len(set(actual_children)) != 2:
        raise ValueError("linked child identities or ordering mismatch")

    for pilot, linkage in ((left, left_link), (right, right_link)):
        child = pilot["package_id"]
        sibling = next(item for item in child_ids if item != child)
        if (
            pilot.get("bundle_id") != linkage["bundle_id"]
            or linkage.get("child_package_id") != child
            or pilot.get("sibling_package_id") != sibling
            or linkage.get("sibling_package_id") != sibling
            or pilot.get("umbrella_run_id") != linkage.get("umbrella_run_id")
        ):
            raise ValueError("linked child identity metadata mismatch")

    if not left.get("umbrella_run_id") or left.get(
        "umbrella_run_id",
    ) != right.get("umbrella_run_id"):
        raise ValueError("linked child run-group mismatch")

    return left_link, right_link, child_ids


def validate_linked_partition_union(
    linkage: dict[str, object],
    target_union: list[str],
    context_union: list[str],
    edge_union: list[dict[str, object]],
    changed_block_ids: list[str],
) -> None:
    """
    Validate exact disjoint target, context, edge, and hunk partitions.

    Parameters
    ----------
    linkage : dict[str, object]
        Shared umbrella ownership metadata.
    target_union : list[str]
        Ordered child target paths combined across reports.
    context_union : list[str]
        Ordered child context paths combined across reports.
    edge_union : list[dict[str, object]]
        Combined child dependency-edge records.
    changed_block_ids : list[str]
        Combined semantic-unit change-block identities.

    Raises
    ------
    ValueError
        If a combined partition contains a gap, overlap, duplicate, or
        contradiction.
    """

    expected_targets = [
        item["path"] for item in linkage["umbrella_target_ownership"]
    ]
    expected_context = [
        item["path"] for item in linkage["umbrella_context_ownership"]
    ]
    expected_edges = list(linkage["umbrella_dependency_edges"])

    if target_union != expected_targets or len(target_union) != len(
        set(target_union),
    ):
        raise ValueError(
            "linked child target partitions have a gap or overlap",
        )

    if context_union != expected_context or len(context_union) != len(
        set(context_union),
    ):
        raise ValueError(
            "linked child context partitions have a gap or overlap",
        )

    edge_identities = [(edge["from"], edge["to"]) for edge in edge_union]
    if len(edge_identities) != len(set(edge_identities)):
        raise ValueError(
            "linked child edge partitions duplicate an edge identity",
        )

    combined_edges = {(edge["from"], edge["to"]): edge for edge in edge_union}

    try:
        reconstructed_edges = [
            combined_edges[(edge["from"], edge["to"])]
            for edge in expected_edges
        ]
    except KeyError as exc:
        raise ValueError("linked child edge partitions have a gap") from exc

    if reconstructed_edges != expected_edges or len(edge_union) != len(
        expected_edges,
    ):
        raise ValueError(
            (
                "linked child edge partitions have a gap, overlap, or "
                "contradiction"
            ),
        )

    if len(changed_block_ids) != len(set(changed_block_ids)):
        raise ValueError("linked child evidence duplicates a changed block")


@dataclass(frozen=True)
class _LinkedChildExpectations:
    """
    Hold one linked child's validated expected inventories.
    """

    child: str
    targets: list[str]
    context: list[str]
    edges: list[dict[str, object]]
    target_fingerprints: dict[str, object]


def _linked_child_expectations(
    pilot: dict[str, object],
    selection: dict[str, object],
    public_root: Path,
) -> _LinkedChildExpectations:
    """
    Validate child ownership and return its expected inventories.

    Parameters
    ----------
    pilot : dict[str, object]
        Linked child pilot report.
    selection : dict[str, object]
        Current child selection.
    public_root : Path
        Public repository root.

    Returns
    -------
    expectations : _LinkedChildExpectations
        Validated child identity and expected inventories.

    Raises
    ------
    ValueError
        If ownership, scope, or fingerprints disagree.
    """

    child = str(pilot["package_id"])
    expected_targets = [str(item["path"]) for item in selection["targets"]]
    expected_context = list(selection["context"])
    expected_edges = list(selection["dependency_edges"])

    if (
        pilot.get("primary_targets") != expected_targets
        or pilot.get("supporting_targets") != []
    ):
        raise ValueError("linked child target ownership mismatch")

    supplied_context = pilot.get("generator_context_files_not_supplied")
    if supplied_context != expected_context:
        raise ValueError("linked child context ownership mismatch")

    if pilot.get("dependency_edges") != expected_edges:
        raise ValueError("linked child edge ownership mismatch")

    expected_scope = selection["rule_scope_report"]
    if pilot.get("rule_scope_report") != expected_scope:
        raise ValueError("linked child rule-scope report mismatch")

    target_fingerprints = pilot.get("target_fingerprints")
    context_fingerprints = pilot.get("context_fingerprints")

    if not isinstance(target_fingerprints, dict) or set(
        target_fingerprints,
    ) != set(expected_targets):
        raise ValueError("linked child target fingerprints are malformed")

    if not isinstance(context_fingerprints, dict) or set(
        context_fingerprints,
    ) != set(expected_context):
        raise ValueError("linked child context fingerprints are malformed")

    for path, recorded in {
        **target_fingerprints,
        **context_fingerprints,
    }.items():
        actual, _ = entry_fingerprint(public_root, path)
        if actual != recorded:
            raise ValueError(f"linked child fingerprint is stale: {path}")

    return _LinkedChildExpectations(
        child=child,
        targets=expected_targets,
        context=expected_context,
        edges=expected_edges,
        target_fingerprints=target_fingerprints,
    )


def _verify_linked_child_rows(
    report: Path,
    pilot: dict[str, object],
    selection: dict[str, object],
    child: str,
    expected_targets: list[str],
    expected_edges: list[dict[str, object]],
) -> int:
    """
    Validate one child report's facts, findings, and review questions.

    Parameters
    ----------
    report : Path
        Child report directory.
    pilot : dict[str, object]
        Linked child pilot report.
    selection : dict[str, object]
        Current child selection.
    child : str
        Child package identity.
    expected_targets : list[str]
        Ordered child target paths.
    expected_edges : list[dict[str, object]]
        Ordered child dependency edges.

    Returns
    -------
    edge_count : int
        Validated dependency-edge fact count.

    Raises
    ------
    ValueError
        If row identity, ownership, or filtering disagrees.
    """

    facts = read_ndjson(report / "facts.ndjson")
    edge_facts = [
        row for row in facts if row.get("topic") == "dependency_edge"
    ]

    if [row.get("value") for row in edge_facts] != expected_edges:
        raise ValueError(
            "linked child edge facts do not match its partition",
        )

    findings = read_ndjson(report / "findings.ndjson")
    questions = read_ndjson(report / "policy_questions.ndjson")
    limitations = read_ndjson(report / "adapter_limitations.ndjson")

    for rows in (facts, findings, questions, limitations):
        if any(
            row.get("run_id") != pilot["run_id"]
            or row.get("bundle_id") != pilot["bundle_id"]
            or row.get("package_id") != child
            or row.get("umbrella_run_id") != pilot["umbrella_run_id"]
            for row in rows
        ):
            raise ValueError("linked child artifact row identity mismatch")

    if any(row.get("path") not in expected_targets for row in findings):
        raise ValueError(
            "linked child finding references a sibling target",
        )

    expected_questions = [
        dict(item) for item in selection["semantic_questions"]
    ]
    if [
        {
            key: row[key]
            for key in (
                "rule_id",
                "topic",
                "question",
                "paths",
                "children",
            )
        }
        for row in questions
    ] != expected_questions:
        raise ValueError(
            "linked child semantic questions are incorrectly filtered",
        )

    expected_limitations = list(selection["adapter_limitations"])
    if [
        {
            key: row[key]
            for key in ("rule_id", "adapter", "limitation", "children")
        }
        for row in limitations
    ] != expected_limitations:
        raise ValueError(
            "linked child limitations are incorrectly filtered",
        )

    return len(edge_facts)


def _verify_linked_child_coverage(
    pilot: dict[str, object],
    selection: dict[str, object],
    public_root: Path,
) -> dict[str, object]:
    """
    Validate one child report's complete semantic-unit coverage.

    Parameters
    ----------
    pilot : dict[str, object]
        Linked child pilot report.
    selection : dict[str, object]
        Current child selection.
    public_root : Path
        Public repository root.

    Returns
    -------
    coverage : dict[str, object]
        Validated current coverage.

    Raises
    ------
    ValueError
        If coverage is incomplete, stale, or malformed.
    """

    coverage = pilot.get("semantic_unit_coverage")

    if (
        not isinstance(coverage, dict)
        or coverage.get("all_changed_blocks_covered") is not True
        or coverage.get("uncovered_changed_blocks")
        or coverage.get("overlapping_units")
        or coverage.get("segmentation_failures")
        or coverage.get("evidence_truncated") is not False
        or any(
            block.get("represented") is not True
            or block.get("exact_diff_retained") is not True
            for block in coverage.get("changed_blocks", [])
        )
    ):
        raise ValueError("linked child evidence coverage failed closed")

    _, current_coverage = configured_semantic_units(public_root, selection)

    if coverage != current_coverage:
        raise ValueError(
            "linked child semantic-unit coverage is stale or malformed",
        )

    return coverage


def _verify_linked_child_semantic_report(
    report: Path,
    pilot: dict[str, object],
    child: str,
    expected_targets: list[str],
    target_fingerprints: dict[str, object],
    coverage: dict[str, object],
    edge_fact_count: int,
) -> None:
    """
    Validate one child's semantic-review metadata and rendered size.

    Parameters
    ----------
    report : Path
        Child report directory.
    pilot : dict[str, object]
        Linked child pilot report.
    child : str
        Child package identity.
    expected_targets : list[str]
        Ordered child target paths.
    target_fingerprints : dict[str, object]
        Exact child target fingerprints.
    coverage : dict[str, object]
        Validated semantic-unit coverage.
    edge_fact_count : int
        Validated dependency-edge fact count.

    Raises
    ------
    ValueError
        If semantic metadata or rendered-size evidence disagrees.
    """

    semantic_rows = read_ndjson(report / "semantic_review.ndjson")
    if len(semantic_rows) != 1:
        raise ValueError(
            "linked child must contain one semantic-review record",
        )

    semantic = semantic_rows[0]
    if (
        semantic.get("primary_paths") != expected_targets
        or semantic.get("target_fingerprints") != target_fingerprints
        or semantic.get("semantic_unit_coverage") != coverage
        or semantic.get("child_package_id") != child
    ):
        raise ValueError("linked child semantic-review metadata mismatch")

    markdown_path = report / str(semantic["path"])
    markdown = markdown_path.read_text(encoding="utf-8")
    markdown_lines = len(markdown.splitlines())
    markdown_bytes = len(markdown.encode("utf-8"))
    semantic_limits = semantic.get("size_limits")
    expected_limits = pilot.get("size_limits")
    byte_limit_exceeded = semantic["size_limits"].get(
        "semantic_markdown_exceeds_byte_limit",
    )
    line_limit_exceeded = semantic["size_limits"].get(
        "semantic_markdown_exceeds_line_limit",
    )
    if (
        semantic.get("markdown_lines") != markdown_lines
        or semantic.get("markdown_bytes") != markdown_bytes
        or semantic_limits != expected_limits
        or byte_limit_exceeded is not False
        or line_limit_exceeded is not False
        or semantic.get("fact_count") != edge_fact_count
    ):
        raise ValueError(
            "linked child Markdown size or fact metadata mismatch",
        )


def _verify_linked_child(
    report: Path,
    pilot: dict[str, object],
    selection: dict[str, object],
    public_root: Path,
) -> tuple[list[str], list[str], list[dict[str, object]], list[str]]:
    """
    Validate one complete child report and return its partition members.

    Parameters
    ----------
    report : Path
        Child report directory.
    pilot : dict[str, object]
        Linked child pilot report.
    selection : dict[str, object]
        Current child selection.
    public_root : Path
        Public repository root.

    Returns
    -------
    targets, context, edges, changed_blocks : tuple[
        list[str], list[str], list[dict[str, object]], list[str]
    ]
        Targets, context, edges, and changed-block identities.
    """

    expected = _linked_child_expectations(
        pilot,
        selection,
        public_root,
    )
    edge_fact_count = _verify_linked_child_rows(
        report,
        pilot,
        selection,
        expected.child,
        expected.targets,
        expected.edges,
    )
    coverage = _verify_linked_child_coverage(
        pilot,
        selection,
        public_root,
    )

    _verify_linked_child_semantic_report(
        report,
        pilot,
        expected.child,
        expected.targets,
        expected.target_fingerprints,
        coverage,
        edge_fact_count,
    )

    changed_block_ids = [
        str(block["block_id"]) for block in coverage["changed_blocks"]
    ]

    return (
        expected.targets,
        expected.context,
        expected.edges,
        changed_block_ids,
    )


def verify_linked_pair_read_only(args: argparse.Namespace) -> int:
    """
    Verify two linked child reports without changing either report tree.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments or an explicit argument sequence.

    Returns
    -------
    status : int
        Verification summary for the linked pair.

    Raises
    ------
    ValueError
        If a supplied value violates the validated contract.
    """

    reports = [
        Path(path).resolve() for path in args.verify_linked_pair_read_only
    ]

    if len(set(reports)) != 2:
        return EXIT_VERIFY_MALFORMED

    for report in reports:
        child_args = argparse.Namespace(**vars(args))
        child_args.verify = None
        child_args.verify_read_only = report
        child_args.verify_linked_pair_read_only = None
        if verify_report(child_args):
            return EXIT_VERIFY_MALFORMED

    try:
        pilots = [
            json.loads(
                (report / "pilot_report.json").read_text(encoding="utf-8"),
            )
            for report in reports
        ]

        for report, pilot in zip(reports, pilots, strict=True):
            if (report / "verification.json").exists() or (
                report / "staleness.ndjson"
            ).exists():
                raise ValueError(
                    (
                        "linked read-only verification requires no persistent "
                        "verification artifacts"
                    ),
                )

            verification = pilot.get("verification", {})
            fresh_verification = verification.get("fresh_verification")
            if fresh_verification != "pending":
                raise ValueError(
                    (
                        "linked read-only verification requires pending "
                        "fresh_verification"
                    ),
                )

        left, _ = pilots
        left_link, _, child_ids = validate_linked_pair_headers(pilots)
        shared_keys = {
            "bundle_id",
            "selection_path",
            "config_fingerprint",
            "target_ownership_fingerprint",
            "graph_ownership_fingerprint",
            "umbrella_target_ownership",
            "umbrella_context_ownership",
            "umbrella_dependency_edges",
            "umbrella_target_count",
            "umbrella_context_count",
            "umbrella_edge_count",
            "child_ids",
        }

        public_root = Path(args.public_root).resolve()
        config_path = public_root / str(left_link["selection_path"])
        raw = json.loads(config_path.read_text(encoding="utf-8"))
        current_link = linked_configuration_metadata(
            raw,
            str(left_link["selection_path"]),
        )
        if any(
            current_link.get(key) != left_link.get(key) for key in shared_keys
        ):
            raise ValueError(
                (
                    "linked report configuration or ownership fingerprint is "
                    "stale"
                ),
            )

        selections: list[dict[str, object]] = []

        for child in child_ids:
            selection = load_target_selection(
                argparse.Namespace(
                    path_values=None,
                    paths_from=Path(str(left_link["selection_path"])),
                    package_child=child,
                    umbrella_run_id=left["umbrella_run_id"],
                ),
                public_root,
            )

            if selection is None:
                raise ValueError("linked child selection did not load")

            selections.append(selection)

        target_union: list[str] = []
        context_union: list[str] = []
        edge_union: list[dict[str, object]] = []
        changed_block_ids: list[str] = []

        for report, pilot, selection in zip(
            reports,
            pilots,
            selections,
            strict=True,
        ):
            (
                child_targets,
                child_context,
                child_edges,
                child_block_ids,
            ) = _verify_linked_child(
                report,
                pilot,
                selection,
                public_root,
            )
            target_union.extend(child_targets)
            context_union.extend(child_context)
            edge_union.extend(child_edges)
            changed_block_ids.extend(child_block_ids)

        validate_linked_partition_union(
            left_link,
            target_union,
            context_union,
            edge_union,
            changed_block_ids,
        )
    except (
        OSError,
        ValueError,
        KeyError,
        TypeError,
        json.JSONDecodeError,
    ) as exc:
        print(f"linked pair verification failed: {exc}", file=sys.stderr)

        return EXIT_VERIFY_MALFORMED

    return 0


def stage_verify_promote(args: argparse.Namespace) -> int:
    """
    Verify final staged package bytes, then atomically promote supplied files.
    """

    report = Path(args.stage_verify_promote).resolve()
    supplied = (
        "semantic_review/download-fastqs-shell-help-pilot.md",
        "findings.ndjson",
        "facts.ndjson",
        "adapter_limitations.ndjson",
        "pilot_report.json",
    )

    with tempfile.TemporaryDirectory(
        prefix="pilot-package-stage-",
        dir=report.parent,
    ) as temporary:
        stage = Path(temporary) / "report"
        shutil.copytree(report, stage)

        pilot_path = stage / "pilot_report.json"
        pilot = json.loads(pilot_path.read_text(encoding="utf-8"))
        pilot["verification"]["staged_exact_byte_verification"] = "passed"
        pilot["package_provenance"]["artifact_hashes"] = {
            name: file_sha256(stage / name)
            for name in supplied
            if name != "pilot_report.json"
        }
        pilot["package_provenance"]["pilot_report_self_hash"]["value"] = (
            pilot_report_self_hash(pilot)
        )
        write_json(pilot_path, pilot)

        stage_args = argparse.Namespace(**vars(args))
        stage_args.verify = None
        stage_args.verify_read_only = stage
        code = verify_report(stage_args)
        if code:
            return code

        for name in supplied:
            destination = report / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            os.replace(stage / name, destination)

    return 0


def coverage_summary(
    all_rules: list[dict[str, object]],
    selected: list[dict[str, object]],
    ledger: list[dict[str, object]],
    checks: list[dict[str, object]],
    findings: list[dict[str, object]],
) -> dict[str, object]:
    """
    Describe what strict-safe execution did and did not cover.
    """

    selected_ids = {str(rule["rule_id"]) for rule in selected}
    availability = {
        str(rule["rule_id"]): str(rule["strict_availability"])
        for rule in all_rules
    }
    executed = {
        str(check["rule_id"])
        for check in checks
        if check.get("record_type") == "transaction"
        and check.get("status") == "completed"
    }
    unavailable_statuses = {
        "unavailable_environment",
        "unavailable_command",
        "infrastructure_error",
    }
    unavailable = {
        str(check["rule_id"])
        for check in checks
        if check.get("status") in unavailable_statuses
    }
    covered = only_unavailable = applicable_no_execution = no_rules = 0

    for row in ledger:
        applicable = set(row["applicable_rules"])
        ran = set(row["checks_run"])

        if ran:
            covered += 1
        elif not applicable:
            no_rules += 1
        elif applicable and all(
            availability.get(rule_id) == "unavailable"
            for rule_id in applicable
        ):
            only_unavailable += 1
        else:
            applicable_no_execution += 1

    findings_by_rule: dict[str, int] = Counter(
        str(row["rule_id"]) for row in findings
    )
    findings_by_severity: dict[str, int] = Counter(
        str(row["severity"]) for row in findings
    )

    return {
        "manifest_rule_total": len(all_rules),
        "selected_rule_total": len(selected),
        "strict_safe_rules_executed": sorted(executed),
        "strict_safe_rules_unavailable": sorted(unavailable),
        "rules_unavailable_in_strict_mode_by_design": sorted(
            str(rule["rule_id"])
            for rule in all_rules
            if rule["strict_availability"] == "unavailable"
        ),
        "executable_check_transactions": sum(
            1 for check in checks if check.get("record_type") == "transaction"
        ),
        "authored_files_covered_by_executed_strict_safe_rule": covered,
        "authored_files_covered_only_by_applicable_unavailable_rules": (
            only_unavailable
        ),
        "authored_files_with_applicable_rules_but_no_executed_rule": (
            applicable_no_execution
        ),
        "authored_files_with_no_applicable_rule": no_rules,
        "findings_total": len(findings),
        "findings_by_rule": dict(findings_by_rule),
        "findings_by_severity": dict(findings_by_severity),
        "selected_rule_ids": sorted(selected_ids),
        "guard_policy": {
            "public": "before and after every executable transaction",
            "private": (
                "run start, before/after each rule, and final completion"
            ),
        },
    }


_RuleTransaction = tuple[str, str, list[str], dict[str, Any]]


@dataclass
class _RuleExecutionContext:
    """
    Carry shared state for deterministic rule transactions.
    """

    public_root: Path
    private_root: Path
    temporary_directory: Path
    selected: list[dict[str, Any]]
    preflight: dict[str, Any]
    execution_paths: list[str]
    configured_scopes: dict[str, list[str]]
    selection: dict[str, Any] | None
    target_records: list[dict[str, Any]]
    consumed_context: list[dict[str, Any]]
    evidence: EvidenceReader
    checks: list[dict[str, Any]]
    findings: list[dict[str, Any]]
    facts: list[dict[str, Any]]
    policy_questions: list[dict[str, Any]]
    limitations: list[dict[str, Any]]
    config: dict[str, Any]
    baseline: dict[str, Any]
    all_rules: list[dict[str, Any]]
    fingerprints: dict[str, str]
    ledger_by_path: dict[str, dict[str, Any]]
    timeout: float
    maximum: int
    run_id: str


@dataclass
class _RuleExecutionState:
    """
    Track abort and diagnostic identity across rule transactions.
    """

    last_rule: str | None = None
    last_transaction: str | None = None
    aborted: bool = False


def _shell_help_transactions(
    context: _RuleExecutionContext,
    rule: dict[str, Any],
    rule_id: str,
    matched: list[str],
) -> list[_RuleTransaction]:
    """
    Build one snapshot-backed Shell-help pilot transaction.

    Parameters
    ----------
    context : _RuleExecutionContext
        Shared roots, evidence, targets, and temporary workspace.
    rule : dict[str, Any]
        Validated Shell-help adapter rule.
    rule_id : str
        Stable rule identifier.
    matched : list[str]
        Selected target paths covered by the transaction.

    Returns
    -------
    transactions : list[_RuleTransaction]
        One snapshot-backed adapter transaction.

    Raises
    ------
    ValueError
        If required adapter context was not declared.
    """

    snapshot_path = context.temporary_directory / f"{rule_id}.snapshot.json"
    adapter_config: dict[str, object] = {}
    registry_path = rule.get("adapter_registry_path")

    if registry_path is not None:
        context_paths = {str(row["path"]) for row in context.consumed_context}
        if context.selection and registry_path not in context_paths:
            raise ValueError(
                (
                    f"adapter registry is not declared pilot context: "
                    f"{registry_path}"
                ),
            )

        if not context.selection and registry_path not in context_paths:
            content = context.evidence.read(
                "public",
                registry_path,
                "manifest-declared adapter registry",
            )
            record = dict(
                context.evidence.records[("public", registry_path)],
            )
            record["content"] = content
            context.consumed_context.append(record)

        adapter_config["registry_path"] = registry_path
        adapter_config["reference_scopes"] = list(
            rule["adapter_reference_scopes"],
        )

    documentation_sources: list[dict[str, str]] = []

    if rule.get("adapter_mode") == "interface_facts" and context.selection:
        interfaces = list(context.selection.get("interfaces", []))
        adapter_config["interfaces"] = interfaces
        source_paths = {
            str(item["documentation_source"]) for item in interfaces
        }
        documentation_sources = [
            {
                "path": str(row["path"]),
                "content": str(row["content"]),
            }
            for row in context.target_records
            if row["path"] in source_paths
        ]

    snapshot = {
        "schema_version": 1,
        "targets": [
            {
                "path": row["path"],
                "role": row["role"],
                "content": row["content"],
                "evidence": row.get("evidence"),
            }
            for row in context.target_records
            if row["path"] in matched
        ],
        "context": [
            {"path": row["path"], "content": row["content"]}
            for row in context.consumed_context
        ],
        "documentation_sources": documentation_sources,
        "adapter_config": adapter_config,
    }
    write_json(snapshot_path, snapshot)
    command = [str(part) for part in rule["command"]] + [
        "--mode",
        str(rule["adapter_mode"]),
        "--snapshot",
        str(snapshot_path),
    ]

    return [
        (
            f"{rule_id}:pilot",
            ".",
            command,
            {
                "transaction_kind": "pilot_snapshot",
                "covered_paths": matched,
                "selection_context_paths": [
                    str(row["path"]) for row in context.consumed_context
                ],
            },
        ),
    ]


def _proportional_proof_transactions(
    context: _RuleExecutionContext,
    rule_id: str,
    matched: list[str],
) -> list[_RuleTransaction]:
    """
    Build exact Git whitespace transactions for selected pilot paths.

    Parameters
    ----------
    context : _RuleExecutionContext
        Shared roots, targets, and Git-state evidence.
    rule_id : str
        Stable proportional-proof rule identifier.
    matched : list[str]
        Selected paths requiring exact whitespace checks.

    Returns
    -------
    transactions : list[_RuleTransaction]
        Git diff transactions appropriate to each path state.
    """

    transactions: list[_RuleTransaction] = []

    for path in matched:
        target = next(
            row for row in context.target_records if row["path"] == path
        )
        labels = set(target["git_state_labels"])

        if "untracked" in labels or "clean_context" in labels:
            transactions.append(
                (
                    f"{rule_id}:no-index:{path}",
                    path,
                    [
                        "git",
                        "diff",
                        "--no-index",
                        "--check",
                        "--no-ext-diff",
                        "/dev/null",
                        str(context.public_root / path),
                    ],
                    {
                        "transaction_kind": "single_path",
                        "covered_paths": [path],
                        "whitespace_mechanism": "git_diff_no_index_check",
                        "whitespace_coverage_relation": "independent",
                        "normal_exit_statuses": [1],
                    },
                ),
            )

            continue

        for command_id, command_spec in (
            ("unstaged", ["git", "diff", "--check", "--", path]),
            (
                "staged",
                ["git", "diff", "--cached", "--check", "--", path],
            ),
        ):
            transactions.append(
                (
                    f"{rule_id}:{command_id}:{path}",
                    path,
                    command_spec,
                    {
                        "transaction_kind": "single_path",
                        "covered_paths": [path],
                        "whitespace_mechanism": (
                            f"git_diff_{command_id}_check"
                        ),
                        "whitespace_coverage_relation": (
                            "exact_existing_git_diff"
                        ),
                    },
                ),
            )

    return transactions


def _rule_transactions(
    context: _RuleExecutionContext,
    rule: dict[str, Any],
    rule_id: str,
    matched: list[str],
) -> list[_RuleTransaction]:
    """
    Build the declared transactions for one available rule.

    Parameters
    ----------
    context : _RuleExecutionContext
        Shared execution and evidence state.
    rule : dict[str, Any]
        Validated rule definition.
    rule_id : str
        Stable rule identifier.
    matched : list[str]
        Paths matched by the rule and configured scope.

    Returns
    -------
    transactions : list[_RuleTransaction]
        Adapter, batch, per-path, or repository transactions for the rule.
    """

    if rule.get("adapter") == "shell_help_pilot":
        return _shell_help_transactions(context, rule, rule_id, matched)

    if rule["scope"] == "per_path" and rule.get("batch_paths"):
        return [
            (
                f"{rule_id}:batch",
                ".",
                [
                    str(part)
                    for part in rule["command"]
                    if str(part) != "{path}"
                ]
                + [str(context.public_root / path) for path in matched],
                {
                    "transaction_kind": "atomic_batch",
                    "covered_paths": matched,
                },
            ),
        ]

    if rule["scope"] == "per_path":
        return [
            (
                f"{rule_id}:{path}",
                path,
                [
                    str(part).replace(
                        "{path}",
                        str(context.public_root / path),
                    )
                    for part in rule["command"]
                ],
                {
                    "transaction_kind": "single_path",
                    "covered_paths": [path],
                },
            )
            for path in matched
        ]

    if rule_id == "TEST.PROOF.PROPORTIONAL" and context.selection:
        return _proportional_proof_transactions(context, rule_id, matched)

    return [
        (
            f"{rule_id}:{command['id']}",
            ".",
            list(command["argv"]),
            {
                "transaction_kind": "single_command",
                "diff_scope": command["id"],
                "covered_paths": matched,
            },
        )
        for command in rule.get(
            "commands",
            [{"id": "default", "argv": rule.get("command", [])}],
        )
    ]


def _record_transaction_coverage(
    context: _RuleExecutionContext,
    rule_id: str,
    extra: dict[str, Any],
) -> None:
    """
    Record completed checks and any path-specific whitespace mechanism.

    Parameters
    ----------
    context : _RuleExecutionContext
        Mutable ledger and target evidence.
    rule_id : str
        Rule completed for the covered paths.
    extra : dict[str, Any]
        Transaction coverage and optional whitespace metadata.
    """

    for path in extra["covered_paths"]:
        if path in context.ledger_by_path:
            context.ledger_by_path[path]["checks_run"].append(rule_id)

        for target in context.target_records:
            if target["path"] != path:
                continue

            target["checks_run"].append(rule_id)

            if "whitespace_mechanism" in extra:
                target["whitespace_coverage"].append(
                    {
                        "mechanism": extra["whitespace_mechanism"],
                        "coverage_relation": extra[
                            "whitespace_coverage_relation"
                        ],
                        "status": "completed",
                    },
                )


def _record_pilot_payload(
    context: _RuleExecutionContext,
    rule: dict[str, Any],
    rule_id: str,
    result: CommandResult,
) -> dict[str, list[dict[str, object]]] | None:
    """
    Validate and record facts, questions, and limitations from one adapter.

    Parameters
    ----------
    context : _RuleExecutionContext
        Mutable pilot evidence and selected paths.
    rule : dict[str, Any]
        Rule whose adapter produced the result.
    rule_id : str
        Stable rule identifier added to normalized evidence.
    result : CommandResult
        Captured adapter command result.

    Returns
    -------
    payload : dict[str, list[dict[str, object]]] | None
        Parsed pilot payload, or None for a nonpilot rule.

    Raises
    ------
    ValueError
        If a pilot fact references a path outside the selection.
    """

    if rule.get("adapter") != "shell_help_pilot":
        return None

    pilot_payload = parse_pilot_payload(result)

    for fact in pilot_payload["facts"]:
        path = str(fact.get("path", ""))
        if path not in context.execution_paths:
            raise ValueError(
                f"pilot fact referenced non-target path: {path}",
            )

        fingerprint = next(
            row["content_fingerprint"]
            for row in context.target_records
            if row["path"] == path
        )
        context.facts.append(
            {
                "schema_version": PILOT_ARTIFACT_SCHEMA_VERSION,
                "run_id": context.run_id,
                "rule_id": rule_id,
                "content_fingerprint": fingerprint,
                **fact,
            },
        )

    for question in pilot_payload["policy_questions"]:
        context.policy_questions.append(
            {
                "schema_version": PILOT_ARTIFACT_SCHEMA_VERSION,
                "run_id": context.run_id,
                "rule_id": rule_id,
                **question,
            },
        )

    for limitation in pilot_payload["limitations"]:
        context.limitations.append(
            {
                "schema_version": PILOT_ARTIFACT_SCHEMA_VERSION,
                "run_id": context.run_id,
                "rule_id": rule_id,
                **limitation,
            },
        )

    return pilot_payload


def _transaction_fragments(
    rule: dict[str, Any],
    result: CommandResult,
    target: str,
    extra: dict[str, Any],
    batch_subresults: list[dict[str, object]] | None,
    pilot_payload: dict[str, list[dict[str, object]]] | None,
) -> list[dict[str, object]]:
    """
    Normalize ordinary, pilot, and Shell-batch findings.

    Parameters
    ----------
    rule : dict[str, Any]
        Rule that owns output parsing.
    result : CommandResult
        Captured command result.
    target : str
        Transaction target used by ordinary parsers.
    extra : dict[str, Any]
        Transaction metadata including accepted nonzero statuses.
    batch_subresults : list[dict[str, object]] | None
        Optional per-path Shell syntax results.
    pilot_payload : dict[str, list[dict[str, object]]] | None
        Optional normalized pilot payload.

    Returns
    -------
    fragments : list[dict[str, object]]
        Normalized finding fragments from every applicable source.
    """

    normal_exit = result.returncode in extra.get(
        "normal_exit_statuses",
        [],
    )
    fragments = (
        parse_result(str(rule["output_parser"]), result, target)
        if result.returncode
        and batch_subresults is None
        and not normal_exit
        and pilot_payload is None
        else []
    )

    if pilot_payload is not None:
        fragments.extend(pilot_payload["findings"])

    if batch_subresults is not None:
        for subresult in batch_subresults:
            if subresult["status"] != "finding":
                continue

            child = CommandResult(
                ["bash", "-n", str(subresult["path"])],
                str(subresult["stdout"]),
                str(subresult["stderr"]),
                int(subresult["exit_status"]),
            )
            fragments.extend(
                parse_result(
                    str(rule["output_parser"]),
                    child,
                    str(subresult["path"]),
                ),
            )

    return fragments


def _record_transaction_findings(
    context: _RuleExecutionContext,
    rule: dict[str, Any],
    fragments: list[dict[str, object]],
    target: str,
    extra: dict[str, Any],
) -> set[str]:
    """
    Normalize finding paths and append evidence-backed finding records.

    Parameters
    ----------
    context : _RuleExecutionContext
        Mutable findings, targets, and source fingerprints.
    rule : dict[str, Any]
        Rule that produced the fragments.
    fragments : list[dict[str, object]]
        Parsed finding fragments.
    target : str
        Fallback target for findings outside the public repository.
    extra : dict[str, Any]
        Transaction metadata copied to finding evidence.

    Returns
    -------
    paths : set[str]
        Repository-relative paths with at least one finding.
    """

    paths_with_findings: set[str] = set()

    for fragment in fragments:
        path = str(fragment["path"])
        candidate = Path(path)

        if candidate.is_absolute():
            try:
                path = relative_path(context.public_root, candidate)
            except ValueError:
                path = target

            fragment["path"] = path

        paths_with_findings.add(path)

        record = finding_record(
            context.run_id,
            rule,
            fragment,
            context.fingerprints.get(path, sha256_bytes(b"repository")),
            str(rule["source_checker"]),
            extra,
        )

        if "topic" in fragment:
            record["topic"] = fragment["topic"]

        context.findings.append(record)

        for target_record in context.target_records:
            if target_record["path"] == path:
                target_record["findings_count"] += 1

    return paths_with_findings


def _run_rule_transaction(
    context: _RuleExecutionContext,
    state: _RuleExecutionState,
    rule: dict[str, Any],
    rule_id: str,
    transaction: _RuleTransaction,
) -> bool:
    """
    Execute, guard, parse, and record one rule transaction.

    Parameters
    ----------
    context : _RuleExecutionContext
        Shared execution, repository, and evidence state.
    state : _RuleExecutionState
        Mutable abort and diagnostic identity.
    rule : dict[str, Any]
        Rule that owns the transaction.
    rule_id : str
        Stable rule identifier.
    transaction : _RuleTransaction
        Command, target, identity, and coverage metadata.

    Returns
    -------
    timed_out : bool
        True when the transaction timed out or failed to launch.
    """

    transaction_id, target, command, extra = transaction
    state.last_transaction = transaction_id
    before = public_guard(
        context.public_root,
        context.fingerprints,
        context.evidence,
    )
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")

    if rule["required_environment"] == "env_protocol":
        env["PYTHONPYCACHEPREFIX"] = str(
            context.temporary_directory / "pycache",
        )

    batch_subresults: list[dict[str, object]] | None = None

    if rule.get("adapter") == "shell_syntax_batch":
        result, batch_subresults = run_shell_syntax_batch(
            extra["covered_paths"],
            context.public_root,
            context.timeout,
            context.maximum,
        )
    else:
        result = run_command(
            command,
            context.public_root,
            env,
            context.timeout,
            context.maximum,
        )

    current = current_fingerprints(
        context.public_root,
        context.config,
        context.baseline,
        context.all_rules,
    )
    after = public_guard(
        context.public_root,
        current,
        context.evidence,
    )

    evidence_stale = context.evidence.verify()
    check = {
        "record_type": "transaction",
        "transaction_id": transaction_id,
        "rule_id": rule_id,
        "target": target,
        "argv": result.argv,
        "timeout_seconds": context.timeout,
        "elapsed_seconds": result.duration_seconds,
        "stdout": result.stdout,
        "stderr": result.stderr,
        "exit_status": result.returncode,
        "timed_out": result.timed_out,
        "launch_error": result.launch_error,
        "public_state_before": before,
        "public_state_after": after,
        **extra,
    }

    if batch_subresults is not None:
        check["subresults"] = batch_subresults

    context.checks.append(check)

    if before != after or evidence_stale:
        check["status"] = "aborted_state_change"
        state.aborted = True

        return False

    if result.timed_out or result.launch_error:
        check["status"] = "infrastructure_error"
        check["infrastructure_error"] = (
            "timeout" if result.timed_out else result.launch_error
        )

        return True

    check["status"] = "completed"
    _record_transaction_coverage(context, rule_id, extra)

    pilot_payload = _record_pilot_payload(
        context,
        rule,
        rule_id,
        result,
    )
    fragments = _transaction_fragments(
        rule,
        result,
        target,
        extra,
        batch_subresults,
        pilot_payload,
    )

    paths_with_findings = _record_transaction_findings(
        context,
        rule,
        fragments,
        target,
        extra,
    )

    if (
        extra["transaction_kind"] == "atomic_batch"
        and batch_subresults is None
    ):
        check["subresults"] = [
            {
                "path": path,
                "status": (
                    "finding"
                    if path in paths_with_findings
                    else "completed_without_reported_failure"
                ),
            }
            for path in extra["covered_paths"]
        ]

    return False


def _execute_rules(
    context: _RuleExecutionContext,
) -> _RuleExecutionState:
    """
    Execute selected rules with per-rule and per-transaction guards.

    Parameters
    ----------
    context : _RuleExecutionContext
        Selected rules and mutable audit evidence.

    Returns
    -------
    state : _RuleExecutionState
        Last completed identities and whether repository mutation aborted work.

    Raises
    ------
    RuntimeError
        If a transaction times out or cannot be launched.
    """

    state = _RuleExecutionState()

    for rule in context.selected:
        rule_id = str(rule["rule_id"])
        state.last_rule = rule_id

        if rule["strict_availability"] != "safe_adapter":
            context.checks.append(
                {
                    "record_type": "availability",
                    "rule_id": rule_id,
                    "status": "unavailable_in_strict_mode",
                    "reason": rule["known_side_effects"],
                },
            )

            continue

        unavailable = rule_available(rule, context.preflight)

        if unavailable:
            context.checks.append(
                {
                    "record_type": "availability",
                    "rule_id": rule_id,
                    "status": "unavailable_environment",
                    "reason": unavailable,
                },
            )

            continue

        private_before = git_visible_state(context.private_root)
        context.checks.append(
            {
                "record_type": "rule_private_guard",
                "rule_id": rule_id,
                "cadence": "per_rule",
                "state_before": private_before,
            },
        )
        matched = [
            path
            for path in context.execution_paths
            if path_matches(path, list(rule["applicable_paths"]))
            and (
                rule_id not in context.configured_scopes
                or path_matches(
                    path,
                    context.configured_scopes[rule_id],
                )
            )
        ]

        if rule["scope"] == "per_path" and not matched:
            context.checks.append(
                {
                    "record_type": "availability",
                    "rule_id": rule_id,
                    "status": "not_applicable_to_selected_targets",
                },
            )

            continue

        transactions = _rule_transactions(
            context,
            rule,
            rule_id,
            matched,
        )
        rule_timeout = False

        for transaction in transactions:
            rule_timeout = _run_rule_transaction(
                context,
                state,
                rule,
                rule_id,
                transaction,
            )
            if state.aborted or rule_timeout:
                break

        private_after = git_visible_state(context.private_root)
        context.checks.append(
            {
                "record_type": "rule_private_guard",
                "rule_id": rule_id,
                "cadence": "per_rule",
                "state_after": private_after,
                "changed": private_after != private_before,
            },
        )

        if private_after != private_before:
            state.aborted = True

        if state.aborted:
            break

        if rule_timeout:
            raise RuntimeError(
                f"infrastructure timeout or launch error in rule {rule_id}",
            )

    return state


@dataclass
class _PromptGenerationContext:
    """
    Carry the shared state for one prompt-generation phase.
    """

    public_root: Path
    report_dir: Path
    selection: dict[str, Any] | None
    target_records: list[dict[str, Any]]
    consumed_context: list[dict[str, Any]]
    facts: list[dict[str, Any]]
    findings: list[dict[str, Any]]
    policy_questions: list[dict[str, Any]]
    limitations: list[dict[str, Any]]
    config: dict[str, Any]
    baseline: dict[str, Any]
    all_rules: list[dict[str, Any]]
    evidence: EvidenceReader
    run_id: str
    runtime_evidence: dict[str, Any] | None


def _write_semantic_artifacts(
    context: _PromptGenerationContext,
    semantic_records: list[dict[str, Any]],
) -> None:
    """
    Write common pilot artifacts after removing in-memory source content.

    Parameters
    ----------
    context : _PromptGenerationContext
        Mutable target, context, and artifact state.
    semantic_records : list[dict[str, Any]]
        Semantic-review records rendered for the selected package.
    """

    for target in context.target_records:
        target.pop("content", None)
        target["run_id"] = context.run_id

    for consumed in context.consumed_context:
        consumed.pop("content", None)

    write_ndjson(
        context.report_dir / "targets.ndjson",
        context.target_records,
    )
    write_ndjson(context.report_dir / "facts.ndjson", context.facts)
    write_ndjson(
        context.report_dir / "policy_questions.ndjson",
        context.policy_questions,
    )
    write_ndjson(
        context.report_dir / "adapter_limitations.ndjson",
        context.limitations,
    )
    write_ndjson(
        context.report_dir / "semantic_review.ndjson",
        semantic_records,
    )
    write_json(
        context.report_dir / "false_positive_status.json",
        {
            "schema_version": PILOT_ARTIFACT_SCHEMA_VERSION,
            "status": "unassessed_pending_independent_review",
        },
    )


def _write_hashed_pilot_report(
    report_dir: Path,
    package_files: list[str],
    pilot_report: dict[str, Any],
) -> None:
    """
    Write one pilot report and then add its artifact and self hashes.

    Parameters
    ----------
    report_dir : Path
        Directory containing the pilot package.
    package_files : list[str]
        Relative package artifacts included in provenance.
    pilot_report : dict[str, Any]
        Mutable report receiving artifact and self hashes.
    """

    write_json(report_dir / "pilot_report.json", pilot_report)
    provenance = pilot_report["package_provenance"]
    provenance["artifact_hashes"] = {
        name: file_sha256(report_dir / name)
        for name in package_files
        if name != "pilot_report.json"
    }
    provenance["pilot_report_self_hash"]["value"] = pilot_report_self_hash(
        pilot_report,
    )
    write_json(report_dir / "pilot_report.json", pilot_report)


def _public_prompt_state(
    context: _PromptGenerationContext,
) -> dict[str, object]:
    """
    Capture the guarded public state used by prompt generation.

    Parameters
    ----------
    context : _PromptGenerationContext
        Public root, inventory configuration, and consumed evidence.

    Returns
    -------
    state : dict[str, object]
        Current guarded public repository state.
    """

    return public_guard(
        context.public_root,
        current_fingerprints(
            context.public_root,
            context.config,
            context.baseline,
            context.all_rules,
        ),
        context.evidence,
    )


def _write_configured_prompt_package(
    context: _PromptGenerationContext,
) -> bool:
    """
    Render one configured semantic package and its provenance report.

    Parameters
    ----------
    context : _PromptGenerationContext
        Validated package selection and mutable prompt artifacts.

    Returns
    -------
    changed : bool
        Whether public source or consumed evidence changed while rendering.

    Raises
    ------
    ValueError
        If package selection is absent or configured units are invalid.
    """

    selection = context.selection
    if selection is None or "package" not in selection:
        raise ValueError(
            "configured prompt package requires package selection",
        )

    fingerprints_by_path = {
        str(row["path"]): str(row["content_fingerprint"])
        for row in context.target_records
    }
    fingerprints_by_path.update(
        {
            str(row["path"]): str(row["fingerprint"])
            for row in context.consumed_context
        },
    )

    for edge in selection["dependency_edges"]:
        context.facts.append(
            {
                "schema_version": PILOT_ARTIFACT_SCHEMA_VERSION,
                "run_id": context.run_id,
                "rule_id": "LAYOUT.PRODUCTION.ROOTS",
                "topic": "dependency_edge",
                "path": edge["from"],
                "content_fingerprint": fingerprints_by_path[str(edge["from"])],
                "certainty": "declared",
                "value": edge,
            },
        )

    context.policy_questions.extend(
        {
            "schema_version": PILOT_ARTIFACT_SCHEMA_VERSION,
            "run_id": context.run_id,
            **question,
        }
        for question in selection["semantic_questions"]
    )
    context.limitations.extend(
        {
            "schema_version": PILOT_ARTIFACT_SCHEMA_VERSION,
            "run_id": context.run_id,
            **limitation,
        }
        for limitation in selection["adapter_limitations"]
    )

    linkage = selection.get("linked_package")

    if linkage:
        row_identity = {
            "bundle_id": linkage["bundle_id"],
            "package_id": linkage["child_package_id"],
            "umbrella_run_id": linkage["umbrella_run_id"],
        }

        for rows in (
            context.facts,
            context.findings,
            context.policy_questions,
            context.limitations,
        ):
            for row in rows:
                row.update(row_identity)

    semantic_units, unit_coverage = configured_semantic_units(
        context.public_root,
        selection,
    )
    prompt_before = _public_prompt_state(context)

    semantic_records = write_configured_bundle(
        context.report_dir,
        selection,
        context.target_records,
        context.facts,
        context.findings,
        context.limitations,
        semantic_units,
        unit_coverage,
        context.run_id,
    )
    prompt_after = _public_prompt_state(context)

    _write_semantic_artifacts(context, semantic_records)

    package_files = list(selection["package"]["supplied_artifacts"])
    pilot_report = {
        "schema_version": PILOT_ARTIFACT_SCHEMA_VERSION,
        "run_id": context.run_id,
        "package_timestamp": time.strftime(
            "%Y-%m-%dT%H:%M:%SZ",
            time.gmtime(),
        ),
        "pilot_outcome": "dry_render_unreviewed",
        "verification": {
            "fresh_verification": "pending",
            "staged_exact_byte_verification": "not_requested",
        },
        "selection": selection["name"],
        "bundle_id": selection["package"]["bundle_id"],
        "primary_targets": [
            row["path"]
            for row in context.target_records
            if row["role"] == "primary"
        ],
        "supporting_targets": [
            row["path"]
            for row in context.target_records
            if row["role"] == "supporting"
        ],
        "generator_context_files_not_supplied": [
            row["path"] for row in context.consumed_context
        ],
        "rule_scopes": selection["rule_scope_report"],
        "rule_scope_report": selection["rule_scope_report"],
        "semantic_only_topics": selection["package"]["semantic_only_topics"],
        "policy_questions": context.policy_questions,
        "dependency_edges": selection["dependency_edges"],
        "dependency_graph_validation": selection[
            "dependency_graph_validation"
        ],
        "semantic_unit_coverage": unit_coverage,
        "markdown_contributions": semantic_records[0][
            "markdown_contributions"
        ],
        "size_limits": semantic_records[0]["size_limits"],
        "runtime_evidence": {
            "required": False,
            "status": "not_declared",
            "source_fingerprint_set_hash": None,
        },
        "proposed_production_changes_not_applied": [],
        "package_provenance": {
            "supplied_artifacts": package_files,
            "artifact_counts": {
                "findings": len(context.findings),
                "facts": len(context.facts),
                "adapter_limitations": len(context.limitations),
                "semantic_review_markdown": 1,
                "pilot_report": 1,
            },
            "artifact_hashes": {},
            "pilot_report_self_hash": {
                "algorithm": "sha256",
                "canonicalization": (
                    "utf8_json_sorted_keys_compact_separators"
                ),
                "excluded_json_pointer": (
                    "/package_provenance/pilot_report_self_hash/value"
                ),
                "value": "",
            },
        },
    }

    if linkage:
        pilot_report.update(
            {
                "package_id": linkage["child_package_id"],
                "sibling_package_id": linkage["sibling_package_id"],
                "umbrella_run_id": linkage["umbrella_run_id"],
                "linked_package": linkage,
                "target_fingerprints": {
                    str(row["path"]): str(row["content_fingerprint"])
                    for row in context.target_records
                },
                "context_fingerprints": {
                    str(row["path"]): str(row["fingerprint"])
                    for row in context.consumed_context
                },
            },
        )

    _write_hashed_pilot_report(
        context.report_dir,
        package_files,
        pilot_report,
    )

    return prompt_before != prompt_after


def _owner_excerpt(
    context: _PromptGenerationContext,
    content_by_path: dict[str, str],
    path: str,
    rule_id: str,
) -> str:
    """
    Return one bounded standards-owner section for a pilot prompt.

    Parameters
    ----------
    context : _PromptGenerationContext
        Prompt limits and repository evidence.
    content_by_path : dict[str, str]
        Standards content keyed by repository-relative path.
    path : str
        Standard containing the requested owner.
    rule_id : str
        Rule identifier whose owner section is required.

    Returns
    -------
    excerpt : str
        Bounded nonempty owner provision.

    Raises
    ------
    ValueError
        If the owner section is missing or empty.
    """

    text = content_by_path[path]
    heading = next(
        (
            line.removeprefix("## ").strip()
            for line in text.splitlines()
            if line.startswith("## ") and f"`{rule_id}`" in line
        ),
        None,
    )
    if heading is None:
        raise ValueError(f"missing owner section {rule_id} in {path}")

    excerpt = bounded_markdown_provision(
        text,
        heading,
        int(context.config["prompting"]["max_rule_chars"]),
    )
    if not excerpt:
        raise ValueError(f"empty owner section {rule_id} in {path}")

    return excerpt


def _write_shell_help_prompt_package(
    context: _PromptGenerationContext,
) -> bool:
    """
    Render one Shell-help pilot package and its provenance report.

    Parameters
    ----------
    context : _PromptGenerationContext
        Validated Shell-help selection and mutable prompt artifacts.

    Returns
    -------
    changed : bool
        Whether public source or consumed evidence changed while rendering.

    Raises
    ------
    ValueError
        If target selection or required owner provisions are absent.
    """

    selection = context.selection
    if selection is None:
        raise ValueError("Shell-help prompt package requires target selection")

    target_content = {
        str(row["path"]): str(row["content"]) for row in context.target_records
    }
    coverage = {
        "target_selection": selection["name"],
        "topics": {
            "type_placeholder_vocabulary_and_syntax": "fully_deterministic",
            "type_placeholder_semantic_appropriateness": (
                "deterministic_extraction_plus_semantic_review"
            ),
            "runtime_requirement_case": "fully_deterministic",
        },
        "whitespace": {
            row["path"]: row["whitespace_coverage"]
            for row in context.target_records
        },
    }
    context_content = {
        str(row["path"]): str(row["content"])
        for row in context.consumed_context
    }
    environment_provision = _owner_excerpt(
        context,
        context_content,
        "docs/standards/help.md",
        "HELP.RUNTIME.REQUIREMENTS",
    )
    review_rules = {
        "help ownership": _owner_excerpt(
            context,
            context_content,
            "docs/standards/help.md",
            "HELP.AUDIENCE",
        ),
        "aliases": _owner_excerpt(
            context,
            context_content,
            "docs/standards/help.md",
            "HELP.ALIAS.PUBLIC",
        ),
        "headings and required sections": _owner_excerpt(
            context,
            context_content,
            "docs/standards/help.md",
            "HELP.SECTION.SCHEMA",
        )
        + _owner_excerpt(
            context,
            context_content,
            "docs/standards/help.md",
            "HELP.EXAMPLES",
        ),
        "types and placeholders": _owner_excerpt(
            context,
            context_content,
            "docs/standards/help.md",
            "HELP.PARAMETER.VOCABULARY",
        ),
        "requiredness/defaults": _owner_excerpt(
            context,
            context_content,
            "docs/standards/help.md",
            "HELP.PARAMETER.VOCABULARY",
        ),
        "runtime requirements": environment_provision,
        "--dir_scr": _owner_excerpt(
            context,
            context_content,
            "docs/standards/shell.md",
            "SHELL.SUBMIT.BOOTSTRAP",
        ),
        "environment handling": environment_provision,
        "stale naming": _owner_excerpt(
            context,
            context_content,
            "docs/standards/help.md",
            "PARAMETER.DESCRIPTIONS",
        ),
        "behavior-sensitive diff review": _owner_excerpt(
            context,
            context_content,
            "docs/standards/shell.md",
            "SHELL.WRAPPER_TOPOLOGY",
        ),
    }
    diff_excerpts = {
        str(row["path"]): focused_diff_excerpt(
            context.public_root,
            str(row["path"]),
            3000,
        )
        for row in context.target_records
        if row["role"] == "primary"
        and "tracked_modified" in row["git_state_labels"]
    }

    anchor_evidence = anchor_evidence_windows(
        context.public_root,
        target_content,
    )

    prompt_before = _public_prompt_state(context)
    semantic_records = write_pilot_bundle(
        context.report_dir,
        context.target_records,
        context.consumed_context,
        context.facts,
        context.findings,
        context.policy_questions,
        coverage,
        lambda path: target_content[path],
        int(context.config["prompting"]["max_excerpt_chars"]),
        review_rules,
        diff_excerpts,
        context.run_id,
        anchor_evidence,
    )
    prompt_after = _public_prompt_state(context)

    _write_semantic_artifacts(context, semantic_records)

    proposed_actions = [
        action
        for row in context.facts
        if row.get("topic") == "proposed_production_actions"
        for action in row.get("value", [])
    ]
    package_files = [
        "semantic_review/download-fastqs-shell-help-pilot.md",
        "findings.ndjson",
        "facts.ndjson",
        "adapter_limitations.ndjson",
        "pilot_report.json",
    ]
    runtime_evidence = context.runtime_evidence

    pilot_report = {
        "schema_version": PILOT_ARTIFACT_SCHEMA_VERSION,
        "run_id": context.run_id,
        "package_timestamp": time.strftime(
            "%Y-%m-%dT%H:%M:%SZ",
            time.gmtime(),
        ),
        "pilot_outcome": "pass_with_revisions",
        "verification": {
            "fresh_verification": "pending",
            "staged_exact_byte_verification": "pending",
        },
        "selection": selection["name"],
        "primary_targets": [
            row["path"]
            for row in context.target_records
            if row["role"] == "primary"
        ],
        "supporting_targets": [
            row["path"]
            for row in context.target_records
            if row["role"] == "supporting"
        ],
        "generator_context_files_not_supplied": [
            row["path"] for row in context.consumed_context
        ],
        "deterministic_topics": coverage["topics"],
        "semantic_only_topics": [
            "topology",
            "environment_handling",
            "description_consistency",
            "line_length",
            "behavior_sensitive_diff",
        ],
        "policy_questions": context.policy_questions,
        "runtime_evidence": {
            "required": runtime_evidence is not None,
            "status": (
                runtime_evidence["status"]
                if runtime_evidence
                else "not_declared"
            ),
            "source_fingerprint_set_hash": (
                runtime_evidence["source_fingerprint_set_hash"]
                if runtime_evidence
                else None
            ),
        },
        "proposed_production_changes_not_applied": proposed_actions,
        "adapter_dispositions": {
            "HELP.SOURCE_STYLE.PILOT": "migrated_to:HELP.SOURCE_STYLE",
            "HELP.COMMAND_REFERENCES.PILOT": (
                "migrated_to:HELP.RUNTIME.REQUIREMENTS"
            ),
            "SHELL.WRAPPER_TOPOLOGY.PILOT": (
                "migrated_to:SHELL.WRAPPER_TOPOLOGY"
            ),
            "DOWNLOAD.INTERFACE.FACTS.PILOT": (
                "migrated_to:HELP.ALIAS.PUBLIC"
            ),
        },
        "package_provenance": {
            "supplied_artifacts": package_files,
            "artifact_counts": {
                "findings": len(context.findings),
                "facts": len(context.facts),
                "adapter_limitations": len(context.limitations),
                "semantic_review_markdown": 1,
                "pilot_report": 1,
            },
            "artifact_hashes": {},
            "pilot_report_self_hash": {
                "algorithm": "sha256",
                "canonicalization": (
                    "utf8_json_sorted_keys_compact_separators"
                ),
                "excluded_json_pointer": (
                    "/package_provenance/pilot_report_self_hash/value"
                ),
                "value": "",
            },
        },
    }

    _write_hashed_pilot_report(
        context.report_dir,
        package_files,
        pilot_report,
    )

    return prompt_before != prompt_after


def _write_general_prompt_package(
    context: _PromptGenerationContext,
    ledger_by_path: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    """
    Render the ordinary finding-oriented prompt package.

    Parameters
    ----------
    context : _PromptGenerationContext
        Findings, rules, prompt limits, and source evidence.
    ledger_by_path : dict[str, dict[str, Any]]
        Audit ledger indexed by repository-relative path.

    Returns
    -------
    records : list[dict[str, Any]]
        Prompt batch records written to the report directory.
    """

    prompt_rules = {str(rule["rule_id"]): rule for rule in context.all_rules}
    batches = build_batches(
        context.findings,
        ledger_by_path,
        int(context.config["prompting"]["max_paths"]),
        int(context.config["prompting"]["max_findings"]),
    )

    def source_excerpt(path: str, line: int | None) -> str:
        if path == ".":
            return (
                "Repository-level finding; inspect recorded checker evidence."
            )

        content = context.evidence.read(
            "public",
            path,
            f"prompt source excerpt for {path}",
        )
        lines = content.splitlines()
        if line is None:
            return "\n".join(lines[:160])

        start = max(0, line - 21)

        return "\n".join(
            f"{number + 1}: {value}"
            for number, value in enumerate(
                lines[start : start + 80],
                start,
            )
        )

    return write_prompt_bundle(
        context.report_dir,
        batches,
        ledger_by_path,
        prompt_rules,
        lambda rule: rule_excerpt(
            context.public_root,
            context.evidence,
            rule,
            int(context.config["prompting"]["max_rule_chars"]),
        ),
        source_excerpt,
        int(context.config["prompting"]["max_excerpt_chars"]),
    )


def _parse_run_args(
    argv: list[str] | None,
) -> tuple[argparse.ArgumentParser, argparse.Namespace]:
    """
    Build the audit parser and parse one argument vector.

    Parameters
    ----------
    argv : list[str] | None
        Optional argument vector; use process arguments when omitted.

    Returns
    -------
    parser, arguments : tuple[argparse.ArgumentParser, argparse.Namespace]
        Configured parser and parsed arguments.
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--public-root", required=True, type=Path)
    parser.add_argument("--private-root", required=True, type=Path)
    parser.add_argument(
        "--rules",
        type=Path,
        default=Path("dev/config/rules.toml"),
    )
    parser.add_argument(
        "--reports-dir",
        type=Path,
        default=Path("artifacts/dev/audit"),
    )
    parser.add_argument("--run-id")

    parser.add_argument("--rule", action="append", dest="rule_ids")
    parser.add_argument("--path", action="append", dest="path_values")
    parser.add_argument("--paths-from", type=Path)
    parser.add_argument("--package-child")
    parser.add_argument("--umbrella-run-id")
    parser.add_argument("--timeout-seconds", type=float)

    parser.add_argument("--verify", type=Path)
    parser.add_argument("--verify-read-only", type=Path)
    parser.add_argument("--verify-linked-pair-read-only", nargs=2, type=Path)
    parser.add_argument("--stage-verify-promote", type=Path)

    parser.add_argument("--capture-controlled-smoke", type=Path)
    parser.add_argument("--normalize-controlled-smoke", type=Path)
    parser.add_argument("--runtime-evidence-out", type=Path)
    parser.add_argument("--runtime-evidence", type=Path)

    parser.add_argument("--allow-partial", action="store_true")

    return parser, parser.parse_args(argv)


@dataclass
class _AuditRunContext:
    """
    Carry mutable state through one complete audit run.

    Parameters
    ----------
    args : argparse.Namespace
        Validated command-line arguments.
    public_root : Path
        Public repository root.
    private_root : Path
        Private evidence repository root.
    config_path : Path
        Absolute rule-manifest path.
    reports_base : Path
        Parent directory for audit reports.
    selection : dict[str, Any] | None
        Optional bounded target selection.
    run_id : str
        Stable identifier for this run.
    """

    args: argparse.Namespace
    public_root: Path
    private_root: Path
    config_path: Path
    reports_base: Path
    selection: dict[str, Any] | None
    run_id: str
    report_dir: Path | None = None
    temporary_directory: Path | None = None
    phase: str = "initialization"
    last_rule: str | None = None
    last_transaction: str | None = None
    aborted: bool = False
    evidence: EvidenceReader | None = None
    runtime_evidence: dict[str, Any] | None = None
    config: dict[str, Any] = dataclass_field(default_factory=dict)
    baseline: dict[str, Any] = dataclass_field(default_factory=dict)
    all_rules: list[dict[str, Any]] = dataclass_field(default_factory=list)
    selected: list[dict[str, Any]] = dataclass_field(default_factory=list)
    configured_scopes: dict[str, list[str]] = dataclass_field(
        default_factory=dict,
    )
    execution_paths: list[str] = dataclass_field(default_factory=list)
    preflight: dict[str, Any] = dataclass_field(default_factory=dict)
    fingerprints: dict[str, str] = dataclass_field(default_factory=dict)
    ledger_by_path: dict[str, dict[str, Any]] = dataclass_field(
        default_factory=dict,
    )
    timeout: float = 0.0
    maximum: int = 0
    checks: list[dict[str, Any]] = dataclass_field(default_factory=list)
    findings: list[dict[str, Any]] = dataclass_field(default_factory=list)
    ledger: list[dict[str, Any]] = dataclass_field(default_factory=list)
    artifacts: list[dict[str, Any]] = dataclass_field(default_factory=list)
    prompt_records: list[dict[str, Any]] = dataclass_field(
        default_factory=list,
    )
    target_records: list[dict[str, Any]] = dataclass_field(
        default_factory=list,
    )
    facts: list[dict[str, Any]] = dataclass_field(default_factory=list)
    policy_questions: list[dict[str, Any]] = dataclass_field(
        default_factory=list,
    )
    limitations: list[dict[str, Any]] = dataclass_field(default_factory=list)
    consumed_context: list[dict[str, Any]] = dataclass_field(
        default_factory=list,
    )


def _dispatch_special_mode(
    parser: argparse.ArgumentParser,
    args: argparse.Namespace,
    public_root: Path,
    selection: dict[str, Any] | None,
) -> int | None:
    """
    Run a requested verification or controlled-smoke mode.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        Parser used for stable command-line diagnostics.
    args : argparse.Namespace
        Parsed command-line arguments.
    public_root : Path
        Public repository root.
    selection : dict[str, Any] | None
        Optional bounded target selection.

    Returns
    -------
    status : int | None
        Mode status, or 'None' when a normal audit should continue.
    """

    if args.stage_verify_promote:
        return stage_verify_promote(args)

    if args.verify_linked_pair_read_only:
        return verify_linked_pair_read_only(args)

    if args.verify or args.verify_read_only:
        return verify_report(args)

    if args.capture_controlled_smoke:
        return capture_controlled_smoke(args, public_root, selection)

    if not args.normalize_controlled_smoke:
        return None

    if args.runtime_evidence_out is None:
        parser.error(
            "--normalize-controlled-smoke requires --runtime-evidence-out",
        )

    return normalize_controlled_smoke(args, public_root, selection)


def _audit_run_context(
    args: argparse.Namespace,
    public_root: Path,
    private_root: Path,
    selection: dict[str, Any] | None,
) -> _AuditRunContext:
    """
    Build and validate the immutable request portion of an audit run.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.
    public_root : Path
        Public repository root.
    private_root : Path
        Private evidence repository root.
    selection : dict[str, Any] | None
        Optional bounded target selection.

    Returns
    -------
    context : _AuditRunContext
        Initialized mutable audit state.

    Raises
    ------
    ValueError
        If the default report location is not ignored.
    """

    config_path = (
        (public_root / args.rules).resolve()
        if not args.rules.is_absolute()
        else args.rules.resolve()
    )
    reports_base = resolve_reports_base(public_root, args.reports_dir)
    ignore_status = git(
        public_root,
        "check-ignore",
        "-q",
        "artifacts/dev/audit/audit-sentinel",
    ).returncode

    if (
        reports_base == public_root / "artifacts/dev/audit"
        and ignore_status != 0
    ):
        raise ValueError(
            (
                "artifacts/dev/audit must be ignored before an audit can "
                "write reports"
            ),
        )

    return _AuditRunContext(
        args=args,
        public_root=public_root,
        private_root=private_root,
        config_path=config_path,
        reports_base=reports_base,
        selection=selection,
        run_id=args.run_id or time.strftime("%Y%m%dT%H%M%SZ", time.gmtime()),
    )


def _initialize_audit_evidence(context: _AuditRunContext) -> None:
    """
    Create the report and load every immutable audit input.

    Parameters
    ----------
    context : _AuditRunContext
        Mutable audit state.

    Raises
    ------
    ValueError
        If report identity, runtime evidence, or runner limits are invalid.
    """

    context.report_dir = context.reports_base / context.run_id
    if context.report_dir.exists():
        raise ValueError(
            f"report directory already exists: {context.report_dir}",
        )

    selection = context.selection
    args = context.args

    if selection and any(
        isinstance(item, dict) and "evidence" in item
        for item in selection["targets"]
    ):
        if args.runtime_evidence is None:
            raise ValueError("pilot selection requires --runtime-evidence")

        context.runtime_evidence = validate_runtime_evidence(
            json.loads(args.runtime_evidence.read_text(encoding="utf-8")),
            context.public_root,
            selection,
            context.run_id,
        )
    elif args.runtime_evidence is not None:
        raise ValueError(
            "--runtime-evidence requires a declared controlled smoke target",
        )

    context.report_dir.mkdir(parents=True)
    context.evidence = EvidenceReader(
        context.public_root,
        context.private_root,
    )

    if selection and selection.get("selection_path"):
        context.evidence.read(
            "public",
            str(selection["selection_path"]),
            "target selection manifest",
        )

    if context.runtime_evidence is not None:
        for source in context.runtime_evidence["source_fingerprints"]:
            context.evidence.read(
                "public",
                str(source["path"]),
                "controlled smoke runtime evidence",
            )

    context.phase = "manifest_loading"
    context.config = tomllib.loads(
        context.evidence.read(
            "public",
            relative_path(context.public_root, context.config_path),
            "rule manifest",
        ),
    )
    configured_timeout = context.config.get("runner", {}).get(
        "timeout_seconds",
        120,
    )
    context.timeout = (
        args.timeout_seconds
        if args.timeout_seconds is not None
        else float(configured_timeout)
    )
    context.maximum = int(
        context.config.get("runner", {}).get("max_output_bytes", 65536),
    )
    if context.timeout <= 0 or context.maximum <= 0:
        raise ValueError(
            "timeout_seconds and max_output_bytes must be positive",
        )

    baseline_path = "dev/config/baseline_cleanup_cohort.json"
    context.phase = "baseline_validation"
    context.baseline = validate_baseline(
        json.loads(
            context.evidence.read(
                "public",
                baseline_path,
                "baseline cohort",
            ),
        ),
        context.config,
    )

    context.phase = "manifest_validation"
    context.all_rules = validate_manifest(context.config)
    (
        context.selected,
        context.configured_scopes,
    ) = configured_rule_selection(
        context.all_rules,
        selection,
        args.rule_ids,
    )


def _prepare_audit_execution(context: _AuditRunContext) -> None:
    """
    Build inventories, selected targets, and environment preflight.

    Parameters
    ----------
    context : _AuditRunContext
        Initialized audit state.

    Raises
    ------
    ValueError
        If required evidence or report state is not initialized.
    """

    if context.evidence is None or context.report_dir is None:
        raise ValueError("audit evidence is not initialized")

    context.phase = "inventory"
    (
        context.ledger,
        context.artifacts,
        context.fingerprints,
    ) = build_inventory(
        context.public_root,
        context.config,
        context.baseline,
        context.all_rules,
    )
    context.ledger_by_path = {str(row["path"]): row for row in context.ledger}

    _, status_rows_raw = parse_status(context.public_root)
    status_rows = {row["path"]: row for row in status_rows_raw}

    if context.selection:
        context.target_records = build_target_records(
            context.public_root,
            context.selection,
            status_rows,
        )
        context.execution_paths = [
            str(row["path"]) for row in context.target_records
        ]

        for target in context.target_records:
            target["content"] = context.evidence.read(
                "public",
                str(target["path"]),
                "selected pilot target",
            )

        for path in context.selection["context"]:
            content = context.evidence.read(
                "public",
                str(path),
                "declared pilot context",
            )
            record = dict(
                context.evidence.records[("public", str(path))],
            )
            record["content"] = content
            context.consumed_context.append(record)
    else:
        context.execution_paths = executable_ledger_paths(context.ledger)

    context.phase = "preflight"
    context.preflight = environment_preflight(
        context.public_root,
        context.private_root,
        context.selected,
        context.timeout,
        context.maximum,
    )
    write_json(
        context.report_dir / "rule_manifest.json",
        {
            "schema_version": context.config["schema_version"],
            "rules": context.all_rules,
            "selected_rule_ids": [
                rule["rule_id"] for rule in context.selected
            ],
        },
    )
    context.temporary_directory = Path(
        tempfile.mkdtemp(prefix="cleanup-audit-"),
    )


def _execute_audit(context: _AuditRunContext) -> None:
    """
    Execute selected rules and reconcile their per-path finding counts.

    Parameters
    ----------
    context : _AuditRunContext
        Prepared audit state.

    Raises
    ------
    ValueError
        If required evidence or temporary execution state is absent.
    """

    if context.evidence is None or context.temporary_directory is None:
        raise ValueError("audit execution is not prepared")

    context.phase = "execution"
    execution_context = _RuleExecutionContext(
        public_root=context.public_root,
        private_root=context.private_root,
        temporary_directory=context.temporary_directory,
        selected=context.selected,
        preflight=context.preflight,
        execution_paths=context.execution_paths,
        configured_scopes=context.configured_scopes,
        selection=context.selection,
        target_records=context.target_records,
        consumed_context=context.consumed_context,
        evidence=context.evidence,
        checks=context.checks,
        findings=context.findings,
        facts=context.facts,
        policy_questions=context.policy_questions,
        limitations=context.limitations,
        config=context.config,
        baseline=context.baseline,
        all_rules=context.all_rules,
        fingerprints=context.fingerprints,
        ledger_by_path=context.ledger_by_path,
        timeout=context.timeout,
        maximum=context.maximum,
        run_id=context.run_id,
    )
    execution_state = _execute_rules(execution_context)
    context.last_rule = execution_state.last_rule
    context.last_transaction = execution_state.last_transaction
    context.aborted = execution_state.aborted

    for finding in context.findings:
        path = str(finding["path"])

        if path in context.ledger_by_path:
            context.ledger_by_path[path]["findings_count"] += 1

    if context.runtime_evidence is None:
        return

    target = next(
        row
        for row in context.target_records
        if row["path"] == context.runtime_evidence["target_path"]
    )
    context.facts.append(
        {
            "schema_version": PILOT_ARTIFACT_SCHEMA_VERSION,
            "run_id": context.run_id,
            "rule_id": "TEST.PROOF.PROPORTIONAL",
            "content_fingerprint": target["content_fingerprint"],
            "topic": "controlled_smoke_execution_result",
            "path": target["path"],
            "certainty": "deterministic",
            "value": context.runtime_evidence,
        },
    )


def _generate_audit_prompts(context: _AuditRunContext) -> None:
    """
    Generate the applicable prompt package and verify consumed evidence.

    Parameters
    ----------
    context : _AuditRunContext
        Executed audit state.
    """

    if (
        context.report_dir is None
        or context.evidence is None
        or context.aborted
    ):
        return

    context.phase = "prompt_generation"
    prompt_context = _PromptGenerationContext(
        public_root=context.public_root,
        report_dir=context.report_dir,
        selection=context.selection,
        target_records=context.target_records,
        consumed_context=context.consumed_context,
        facts=context.facts,
        findings=context.findings,
        policy_questions=context.policy_questions,
        limitations=context.limitations,
        config=context.config,
        baseline=context.baseline,
        all_rules=context.all_rules,
        evidence=context.evidence,
        run_id=context.run_id,
        runtime_evidence=context.runtime_evidence,
    )

    if context.selection and "package" in context.selection:
        prompt_changed = _write_configured_prompt_package(prompt_context)
        context.aborted = context.aborted or prompt_changed
    elif context.selection:
        prompt_changed = _write_shell_help_prompt_package(prompt_context)
        context.aborted = context.aborted or prompt_changed
    else:
        context.prompt_records = _write_general_prompt_package(
            prompt_context,
            context.ledger_by_path,
        )

    if context.evidence.verify():
        context.aborted = True


def _write_completed_audit(context: _AuditRunContext) -> int:
    """
    Write canonical artifacts and the final completed-or-aborted run record.

    Parameters
    ----------
    context : _AuditRunContext
        Executed audit state.

    Returns
    -------
    status : int
        Stable completed or guard-aborted process status.

    Raises
    ------
    ValueError
        If required evidence or report state is not initialized.
    """

    if context.report_dir is None or context.evidence is None:
        raise ValueError("audit report is not initialized")

    context.phase = "report_writing"
    write_ndjson(context.report_dir / "checks.ndjson", context.checks)
    write_ndjson(context.report_dir / "ledger.ndjson", context.ledger)
    write_ndjson(context.report_dir / "artifacts.ndjson", context.artifacts)
    write_ndjson(context.report_dir / "findings.ndjson", context.findings)
    write_json(
        context.report_dir / "cohort_progress.json",
        cohort_progress(
            context.public_root,
            context.baseline,
            context.ledger,
        ),
    )

    _generate_audit_prompts(context)
    write_ndjson(
        context.report_dir / "prompts.ndjson",
        context.prompt_records,
    )

    baseline_path = "dev/config/baseline_cleanup_cohort.json"
    run = {
        "run_id": context.run_id,
        "report_format_version": 2 if context.selection else 1,
        "target_selection": (
            context.selection["name"] if context.selection else None
        ),
        "status": "aborted" if context.aborted else "completed",
        "abort_reason": (
            "repository_or_consumed_evidence_changed"
            if context.aborted
            else None
        ),
        "preflight": context.preflight,
        "baseline_path": baseline_path,
        "baseline_fingerprint": entry_fingerprint(
            context.public_root,
            baseline_path,
        )[0],
        "rule_manifest_path": relative_path(
            context.public_root,
            context.config_path,
        ),
        "rule_manifest_fingerprint": entry_fingerprint(
            context.public_root,
            relative_path(context.public_root, context.config_path),
        )[0],
        "consumed_evidence": list(context.evidence.records.values()),
        "private_evidence": [
            record
            for record in context.evidence.records.values()
            if record["repository"] == "private"
        ],
        "private_final_git": git_visible_state(context.private_root),
        "coverage_summary": coverage_summary(
            context.all_rules,
            context.selected,
            context.ledger,
            context.checks,
            context.findings,
        ),
        "partial_reports": [
            "rule_manifest.json",
            "checks.ndjson",
            "ledger.ndjson",
            "artifacts.ndjson",
            "findings.ndjson",
            "cohort_progress.json",
            "prompts.ndjson",
        ],
    }

    if context.selection and context.selection.get("linked_package"):
        run["linked_package"] = {
            key: context.selection["linked_package"][key]
            for key in (
                "bundle_id",
                "child_package_id",
                "sibling_package_id",
                "umbrella_run_id",
                "config_fingerprint",
                "target_ownership_fingerprint",
                "graph_ownership_fingerprint",
            )
        }

    write_json(context.report_dir / "run.json", run)

    return EXIT_GUARD_ABORTED if context.aborted else 0


def _write_failed_audit(
    context: _AuditRunContext,
    exc: Exception,
) -> int:
    """
    Preserve available partial artifacts and write a stable error record.

    Parameters
    ----------
    context : _AuditRunContext
        Audit state at the failed phase.
    exc : Exception
        Exception that terminated the run.

    Returns
    -------
    status : int
        Stable runtime-error process status.

    Raises
    ------
    Exception
        Re-raises 'exc' when no report directory is available.
    """

    if context.report_dir is None:
        raise exc

    partial: list[str] = []

    try:
        if context.ledger:
            write_ndjson(
                context.report_dir / "ledger.ndjson",
                context.ledger,
            )
            partial.append("ledger.ndjson")

        if context.artifacts:
            write_ndjson(
                context.report_dir / "artifacts.ndjson",
                context.artifacts,
            )
            partial.append("artifacts.ndjson")

        for name, records in (
            ("checks.ndjson", context.checks),
            ("findings.ndjson", context.findings),
            ("prompts.ndjson", context.prompt_records),
        ):
            write_ndjson(context.report_dir / name, records)
            partial.append(name)
    except Exception:
        pass

    write_json(
        context.report_dir / "run.json",
        {
            "run_id": context.run_id,
            "status": "error",
            "abort_reason": None,
            "error_type": type(exc).__name__,
            "error_message": str(exc),
            "error_traceback": traceback.format_exc(limit=12),
            "failed_phase": context.phase,
            "last_completed_rule": context.last_rule,
            "last_completed_check_transaction": context.last_transaction,
            "partial_reports_available": bool(partial),
            "partial_reports": partial,
            "exit_code_category": "runtime_error",
        },
    )

    return EXIT_RUNTIME_ERROR


def main(argv: list[str] | None = None) -> int:
    """
    Run one audit or verify one existing report bundle.

    Parameters
    ----------
    argv : list[str] | None
        Optional argument vector; use process arguments when omitted.

    Returns
    -------
    status : int
        Stable audit or verification process status.

    Raises
    ------
    ValueError
        If roots, modes, evidence, or report state violate audit contracts.
    """

    parser, args = _parse_run_args(argv)

    public_root = resolve_repo(args.public_root, "public root")
    private_root = resolve_repo(args.private_root, "private root")
    if public_root == private_root:
        raise ValueError("public and private roots must differ")

    args.public_root, args.private_root = public_root, private_root

    modes = (
        args.verify,
        args.verify_read_only,
        args.verify_linked_pair_read_only,
        args.stage_verify_promote,
        args.capture_controlled_smoke,
        args.normalize_controlled_smoke,
    )

    if sum(value is not None for value in modes) > 1:
        parser.error(
            "verification and controlled-smoke modes are mutually exclusive",
        )

    selection = load_target_selection(args, public_root)

    if (
        selection
        and selection.get("linked_umbrella")
        and not selection.get("linked_package")
    ):
        parser.error("linked package generation requires --package-child")

    if (
        selection
        and selection.get("linked_package")
        and not args.umbrella_run_id
    ):
        parser.error("linked package generation requires --umbrella-run-id")

    mode_status = _dispatch_special_mode(
        parser,
        args,
        public_root,
        selection,
    )

    if mode_status is not None:
        return mode_status

    context = _audit_run_context(
        args,
        public_root,
        private_root,
        selection,
    )

    try:
        _initialize_audit_evidence(context)
        _prepare_audit_execution(context)
        _execute_audit(context)

        return _write_completed_audit(context)
    except Exception as exc:
        return _write_failed_audit(context, exc)
    finally:
        if context.temporary_directory is not None:
            shutil.rmtree(
                context.temporary_directory,
                ignore_errors=True,
            )


if __name__ == "__main__":
    raise SystemExit(main())
