#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: python_version_policy.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


"""
Inventory and enforce the repository-wide Python version policy.
"""

from __future__ import annotations

import argparse
import ast
import dataclasses
import json
import re
import sys
import tomllib
from collections import Counter
from pathlib import Path

MINIMUM_PYTHON = (3, 11)
POLICY_LABEL = "Python >= 3.11"
CANONICAL_GUARD = "sys.version_info >= (3, 11)"
CANONICAL_MESSAGE = "Python >= 3.11 required."
STALE_MESSAGE = "Python >= " + "3.10 required."

RULE_ID = "PY.VERSION.FLOOR"

VERSION_GUARD = re.compile(
    r"sys\.version_info\s*>=\s*\(\s*(?P<major>\d+)\s*,\s*(?P<minor>\d+)\s*\)",
)
PYTHON_REQUIREMENT = re.compile(
    r"^\s*-\s+python(?:=|>=)(?P<major>\d+)\.(?P<minor>\d+)"
    r"(?:\.\d+)?(?:\s|$)",
    re.M,
)


@dataclasses.dataclass(frozen=True)
class PythonSource:
    """
    One maintained Python source covered by the shared policy.
    """

    path: str
    status: str
    role: str
    minimum: str
    shebang: str
    local_guard: str
    policy_coverage: str
    direct_execution_tests: tuple[str, ...]

    def as_dict(self) -> dict[str, object]:
        """
        Return a deterministic JSON-ready record.
        """

        return dataclasses.asdict(self)


@dataclasses.dataclass(frozen=True)
class ExcludedPythonSource:
    """
    One Python source deliberately outside maintained policy scope.
    """

    path: str
    status: str
    reason: str

    def as_dict(self) -> dict[str, object]:
        """
        Return a deterministic JSON-ready record.
        """

        return dataclasses.asdict(self)


@dataclasses.dataclass(frozen=True)
class Finding:
    """
    One deterministic Python-policy violation.
    """

    rule_id: str
    path: str
    line: int
    message: str

    def as_dict(self) -> dict[str, object]:
        """
        Return a deterministic JSON-ready record.
        """

        return dataclasses.asdict(self)

    def format(self) -> str:
        """
        Render one line-oriented diagnostic.
        """

        return f"{self.rule_id}: {self.path}:{self.line}: {self.message}"


def is_supported_version(version: tuple[int, ...]) -> bool:
    """
    Return whether one version tuple satisfies the repository floor.
    """

    return version[:2] >= MINIMUM_PYTHON


def environment_python_version(text: str) -> tuple[int, int] | None:
    """
    Return the environment's declared Python major/minor version.
    """

    match = PYTHON_REQUIREMENT.search(text)
    if match is None:
        return None

    return int(match.group("major")), int(match.group("minor"))


def parse_floor_syntax(text: str, filename: str = "<unknown>") -> ast.Module:
    """
    Parse source using the repository's Python 3.11 grammar floor.
    """

    return ast.parse(
        text,
        filename=filename,
        feature_version=MINIMUM_PYTHON,
    )


def maintained_python_paths(root: Path) -> list[Path]:
    """
    Return the source-derived maintained Python file universe.
    """

    paths = {
        *(root / "src").rglob("*.py"),
        *(root / "dev").rglob("*.py"),
        *(root / "tests").rglob("*.py"),
    }

    return sorted(
        path
        for path in paths
        if "__pycache__" not in path.parts
        and "artifacts/tests" not in path.as_posix()
    )


def excluded_python_paths(root: Path) -> list[ExcludedPythonSource]:
    """
    Return explicit Python exclusions from maintained-source discovery.
    """

    return []


def has_entrypoint(tree: ast.Module) -> bool:
    """
    Return whether a module has a top-level '__main__' guard.
    """

    return any(
        isinstance(node, ast.If) and "__main__" in ast.unparse(node.test)
        for node in tree.body
    )


def direct_tests(path: str, role: str) -> tuple[str, ...]:
    """
    Return the established execution gate covering one source role.
    """

    if (
        path.startswith("src/protocol_chipseq_signal_norm/cli/")
        and role == "entry_point"
    ):
        return ("tests/contract/repository/test_python_startup.sh",)

    if path.startswith("src/protocol_chipseq_signal_norm/"):
        return (
            "tests/unit",
            "tests/contract/repository/test_python_startup.sh",
        )

    if path.startswith("tests/unit/dev_audit/"):
        return ("dev/audit unit suite",)

    if path.startswith("dev/audit/"):
        return ("dev/audit unit suite", "complete audit report")

    return ("tests/unit",)


def local_guard(text: str) -> str:
    """
    Classify one source's optional local version guard.
    """

    tree = parse_floor_syntax(text)

    for node in ast.walk(tree):
        if not isinstance(node, ast.Assert):
            continue

        expression = ast.get_source_segment(text, node.test) or ast.unparse(
            node.test,
        )
        match = VERSION_GUARD.search(expression)
        if match is not None:
            return f"Python >= {match.group('major')}.{match.group('minor')}"

    return "none; covered centrally"


def inventory_repository(root: Path) -> list[PythonSource]:
    """
    Return every maintained Python file with source-derived policy coverage.
    """

    root = root.resolve()
    rows: list[PythonSource] = []

    for source in maintained_python_paths(root):
        text = source.read_text(encoding="utf-8")
        tree = parse_floor_syntax(text, filename=str(source))
        role = "entry_point" if has_entrypoint(tree) else "imported_module"
        path = str(source.relative_to(root))
        lines = text.splitlines()
        rows.append(
            PythonSource(
                path=path,
                status="maintained",
                role=role,
                minimum=POLICY_LABEL,
                shebang=lines[0]
                if lines and lines[0].startswith("#!")
                else "none",
                local_guard=local_guard(text),
                policy_coverage="dev/audit/python_version_policy.py",
                direct_execution_tests=direct_tests(path, role),
            ),
        )

    return rows


def guard_findings(path: str, text: str) -> list[Finding]:
    """
    Return stale local-version-guard findings for one Python source.
    """

    findings: list[Finding] = []

    try:
        tree = parse_floor_syntax(text, filename=path)
    except SyntaxError:
        return findings

    for node in ast.walk(tree):
        if not isinstance(node, ast.Assert):
            continue

        expression = ast.get_source_segment(text, node.test) or ast.unparse(
            node.test,
        )
        match = VERSION_GUARD.search(expression)
        if match is None:
            continue

        version = (int(match.group("major")), int(match.group("minor")))

        if version != MINIMUM_PYTHON:
            findings.append(
                Finding(
                    RULE_ID,
                    path,
                    node.lineno,
                    f"facet=guard; local guard must use '{CANONICAL_GUARD}'",
                ),
            )

        message = (
            node.msg.value
            if isinstance(node.msg, ast.Constant)
            and isinstance(node.msg.value, str)
            else None
        )

        if message != CANONICAL_MESSAGE:
            findings.append(
                Finding(
                    RULE_ID,
                    path,
                    node.lineno,
                    (
                        f"facet=guard; local guard message must use "
                        f"'{CANONICAL_MESSAGE}'"
                    ),
                ),
            )

    return findings


def scan_repository(root: Path) -> tuple[list[Finding], list[PythonSource]]:
    """
    Audit interpreter, documentation, environment, and source guards.

    Parameters
    ----------
    root : Path
        Repository root containing policy, environment, and Python sources.

    Returns
    -------
    findings, sources : tuple[list[Finding], list[PythonSource]]
        Policy diagnostics and the complete maintained source inventory.
    """

    root = root.resolve()
    findings: list[Finding] = []

    if not is_supported_version(sys.version_info[:2]):
        findings.append(
            Finding(
                RULE_ID,
                "<interpreter>",
                1,
                (
                    f"facet=interpreter; repository checks require "
                    f"{POLICY_LABEL}; observed "
                    f"{sys.version.split()[0]}"
                ),
            ),
        )

    style = root / "docs/standards/python.md"
    style_text = style.read_text(encoding="utf-8") if style.is_file() else ""

    if POLICY_LABEL not in style_text:
        findings.append(
            Finding(
                RULE_ID,
                "docs/standards/python.md",
                1,
                f"facet=documentation; owner must state '{POLICY_LABEL}'",
            ),
        )

    environment = root / "install/envs/env_protocol.yml"
    environment_text = (
        environment.read_text(encoding="utf-8")
        if environment.is_file()
        else ""
    )
    environment_version = environment_python_version(environment_text)

    if environment_version is None or not is_supported_version(
        environment_version,
    ):
        findings.append(
            Finding(
                RULE_ID,
                "install/envs/env_protocol.yml",
                1,
                (
                    "facet=environment; authoritative environment must select "
                    "Python 3.11 or newer"
                ),
            ),
        )

    metadata = root / "pyproject.toml"
    metadata_data = (
        tomllib.loads(metadata.read_text(encoding="utf-8"))
        if metadata.is_file()
        else {}
    )

    if metadata_data.get("project", {}).get("requires-python") != ">=3.11":
        findings.append(
            Finding(
                RULE_ID,
                "pyproject.toml",
                1,
                "facet=metadata; requires-python must be exactly '>=3.11'",
            ),
        )

    ruff_target = (
        metadata_data.get("tool", {}).get("ruff", {}).get("target-version")
    )

    if ruff_target != "py311":
        findings.append(
            Finding(
                RULE_ID,
                "pyproject.toml",
                1,
                "facet=ruff_target; Ruff target-version must be 'py311'",
            ),
        )

    inventory: list[PythonSource] = []

    for source in maintained_python_paths(root):
        path = str(source.relative_to(root))
        text = source.read_text(encoding="utf-8")

        try:
            parse_floor_syntax(text, filename=str(source))
        except SyntaxError as exc:
            findings.append(
                Finding(
                    RULE_ID,
                    path,
                    exc.lineno or 1,
                    f"facet=source; syntax error: {exc.msg}",
                ),
            )

            continue

        findings.extend(guard_findings(path, text))

    if not any("facet=source" in finding.message for finding in findings):
        inventory = inventory_repository(root)

    return sorted(
        findings,
        key=lambda item: (item.path, item.line, item.rule_id, item.message),
    ), inventory


def write_json(path: Path, value: object) -> None:
    """
    Write deterministic, newline-terminated JSON.
    """

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse Python-policy audit arguments.

    Parameters
    ----------
    argv : list[str] | None
        Explicit arguments, or None to read the process arguments.

    Returns
    -------
    arguments : argparse.Namespace
        Parsed repository and JSON-output options.
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--inventory-output", type=Path)
    parser.add_argument("--strict", action="store_true")

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """
    Print findings and deterministic maintained-source summary counts.

    Parameters
    ----------
    argv : list[str] | None
        Explicit arguments, or None to read the process arguments.

    Returns
    -------
    status : int
        Zero when version policy passes and one when findings remain.
    """

    args = parse_args(argv)
    root = args.root.resolve()

    findings, inventory = scan_repository(root)
    excluded = excluded_python_paths(root)

    if args.inventory_output:
        write_json(
            args.inventory_output,
            {
                "excluded": [row.as_dict() for row in excluded],
                "findings": [finding.as_dict() for finding in findings],
                "maintained": [row.as_dict() for row in inventory],
                "minimum": POLICY_LABEL,
            },
        )

    for finding in findings:
        print(finding.format())

    roles = Counter(row.role for row in inventory)
    print(
        f"Python policy: {len(inventory)} maintained file(s), "
        f"{roles['entry_point']} entry point(s), "
        f"{roles['imported_module']} imported/test module(s), "
        f"{len(excluded)} explicit exclusion(s), "
        f"{len(findings)} finding(s)",
    )

    return 1 if args.strict and findings else 0


if __name__ == "__main__":
    raise SystemExit(main())
