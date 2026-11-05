#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: doc_coverage.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6);
# - Anthropic Claude Code (Opus 5, Opus 5.5).
#
# Distributed under the MIT license.


"""
Find maintained Python and Shell functions that lack needed documentation.

Python modules under 'src/protocol_chipseq_signal_norm' need a module
docstring, and their public and nontrivial functions need docstrings. Shell
functions under 'bin/', 'lib/bash/', 'install/scripts/', and 'tests/support/'
need local help when they are public or complex, and a no-op heredoc is
reported rather than taken for documentation. Run from the repository root.
Prints one '<severity><tab><path><tab><line><tab><message>' row per finding,
with severity 'FAIL' or 'WARN', and exits 0; the calling test reads the rows.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

PYTHON_ROOTS = (Path("src/protocol_chipseq_signal_norm"),)
ROOTS_SH = (
    Path("bin"),
    Path("lib/bash"),
    Path("install/scripts"),
    Path("tests/support"),
)

FUNC_RE = re.compile(
    r"^(?:function\s+([A-Za-z_][A-Za-z0-9_]*)\s*\(\)|"
    r"([A-Za-z_][A-Za-z0-9_]*)\s*\(\)\s*\{)",
)
LIFECYCLE_RE = re.compile(
    r"^(?:source_helpers.*|resolve_dir_scr|init_args_hardcoded|init_arg_defs|"
    r"init_defs|parse_args|canonicalize_args|validate_args|prepare_vecs|"
    r"validate_vecs|config_exec|setup_env|check_tools|print_.*debug|"
    r"report_plan|report_results|run_jobs|main)$",
)
SECTION_RE = re.compile(
    r"(?m)^[ \t]*(?:Usage|Parameters|Returns|Notes|References|See Also|"
    r"Examples)\n^[ \t]*-+$",
)
NOOP_HEREDOC_RE = re.compile(r"(?m)^[ \t]*:[ \t]*<<['\"]?EOM['\"]?")
COMPLEX_SHELL_RE = re.compile(
    r"\b(?:getopts|case |while |for |if |awk |perl |python|samtools|"
    r"curl|wget|sbatch|parallel|mk"
    r"temp|rm |mv |cp |ln |printf )\b|"
    r"declare -g|local -n|arr_cmd_bld|Expected globals|Generated globals",
)


def emit(severity: str, path: Path, line: int, message: str) -> None:
    print(f"{severity}\t{path}\t{line}\t{message}")


def list_python_files() -> list[Path]:
    files = {
        path
        for root in PYTHON_ROOTS
        if root.exists()
        for path in root.rglob("*.py")
        if path.name != "__init__.py"
    }

    return sorted(files)


def shell_files() -> list[Path]:
    files: set[Path] = set()

    for root in ROOTS_SH:
        if not root.exists():
            continue

        files.update(root.rglob("*.sh"))

    return sorted(files)


def is_public_python(node: ast.AST) -> bool:
    name = getattr(node, "name", "")
    return bool(name and not name.startswith("_"))


def is_tiny_python_helper(node: ast.AST) -> bool:
    body = getattr(node, "body", [])
    branch_kinds = (ast.If, ast.For, ast.While, ast.Try, ast.With, ast.Match)
    statements = [stmt for stmt in body if not isinstance(stmt, ast.Expr)]
    if len(statements) <= 1:
        return True

    return len(body) <= 3 and not any(
        isinstance(child, branch_kinds)
        for child in ast.walk(node)
        if child is not node
    )


def check_python() -> None:
    for path in list_python_files():
        tree = ast.parse(path.read_text(), filename=str(path))

        if not ast.get_docstring(tree):
            emit("FAIL", path, 1, "module docstring is required")

        node_kinds = (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)

        for node in tree.body:
            if not isinstance(node, node_kinds):
                continue

            if is_public_python(node) and not ast.get_docstring(node):
                emit(
                    "FAIL",
                    path,
                    node.lineno,
                    f"{node.name} needs a docstring",
                )

            for child in ast.walk(node):
                if child is node:
                    continue

                if not isinstance(child, node_kinds):
                    continue

                if ast.get_docstring(child) or is_tiny_python_helper(child):
                    continue

                emit(
                    "WARN",
                    path,
                    child.lineno,
                    f"{child.name} is a nontrivial nested/private helper "
                    "without docs",
                )


def function_blocks(path: Path) -> list[tuple[int, str, str]]:
    lines = path.read_text(errors="ignore").splitlines()
    starts: list[tuple[int, str]] = []

    for idx, line in enumerate(lines):
        match = FUNC_RE.match(line)

        if match:
            starts.append((idx, match.group(1) or match.group(2)))

    blocks: list[tuple[int, str, str]] = []

    for idx, (start, name) in enumerate(starts):
        end = starts[idx + 1][0] if idx + 1 < len(starts) else len(lines)
        blocks.append((start + 1, name, "\n".join(lines[start:end])))

    return blocks


def is_tiny_shell_helper(body: str) -> bool:
    code_lines = [
        line
        for line in body.splitlines()[1:]
        if line.strip()
        and not line.lstrip().startswith("#")
        and line.strip() not in {"}", ";;"}
    ]

    return len(code_lines) <= 5


def has_local_help(body: str) -> bool:
    return bool(SECTION_RE.search(body))


def check_shell() -> None:
    for path in shell_files():
        if NOOP_HEREDOC_RE.search(path.read_text(errors="ignore")):
            emit(
                "FAIL",
                path,
                1,
                "do not use no-op heredoc blocks as documentation coverage",
            )

        for line, name, body in function_blocks(path):
            if has_local_help(body):
                continue

            if path.match("lib/bash/help/*.sh"):
                continue

            if LIFECYCLE_RE.match(name):
                continue

            if is_tiny_shell_helper(body):
                continue

            private = name.startswith("_")
            public_reusable = path.match("lib/bash/**/*.sh") and not private
            shared_test = path.match("tests/support/*.sh") and not private
            complex_private = private and COMPLEX_SHELL_RE.search(body)

            if public_reusable:
                emit(
                    "WARN",
                    path,
                    line,
                    f"{name} public reusable helper lacks local help",
                )
            elif shared_test:
                emit(
                    "WARN",
                    path,
                    line,
                    f"{name} shared test helper lacks local help",
                )
            elif complex_private:
                emit(
                    "WARN",
                    path,
                    line,
                    f"{name} complex private helper lacks local help",
                )


def main() -> int:
    """
    Print the Python findings, then the Shell findings.
    """

    check_python()
    check_shell()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
