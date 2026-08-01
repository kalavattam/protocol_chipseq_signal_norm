#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: numeric_emission_inventory.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.

"""
Inventory bounded Python and Shell emission syntax without semantic choices.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import re
from pathlib import Path
from typing import Any

PYTHON_CLASSES = [
    "python.direct_print",
    "python.sys_stream_write",
    "python.bound_text_stream_write",
    "python.csv_writer_row",
    "python.json_dump",
    "python.path_write_text",
    "python.same_sink_formatting_expression",
]
SHELL_CLASSES = [
    "shell.direct_printf_or_builtin_printf",
    "shell.direct_echo",
    "shell.static_redirection",
    "shell.literal_heredoc_or_herestring",
]
UNSUPPORTED = [
    "dynamic_or_indirect_call_target",
    "function_or_command_variable_dispatch",
    "eval_or_eval_like_construct",
    "sourced_code_outside_scanned_statement",
    "generated_command_or_source",
    "unresolved_command_substitution",
    "embedded_language",
    "dynamic_file_descriptor_or_stream_alias",
    "unrecognized_writer_factory",
    "otherwise_unsupported_construct",
]
NUMERIC_HINT = re.compile(
    r"(?i)(?:\bdp\b|round|format|decimal|precision|float|int|numeric|"
    r"number|value|coef|ratio|depth|pseudo|floor|eps|\\.\\d*f|%[.0-9*]*[fdeg])",
)


def maintained_paths(root: Path, languages: set[str]) -> list[Path]:
    """
    Return exact maintained source candidates for active languages.
    """

    paths: set[Path] = set()
    if "python" in languages:
        for directory in ("src", "dev", "tests"):
            paths.update((root / directory).rglob("*.py"))
    if "shell" in languages:
        for directory in ("bin", "lib", "install", "tests"):
            paths.update((root / directory).rglob("*.sh"))
    return sorted(path for path in paths if "__pycache__" not in path.parts)


def _record(
    *,
    relative: str,
    language: str,
    line: int,
    column: int,
    symbol: str,
    syntactic_class: str,
    source: str,
    unsupported_reason: str | None = None,
) -> dict[str, Any]:
    """
    Build one deterministic unresolved evidence record.
    """

    payload = f"{relative}:{line}:{column}:{symbol}:{syntactic_class}:{source}"
    digest = hashlib.sha256(payload.encode()).hexdigest()
    base: dict[str, Any] = {
        "site_id": digest[:16],
        "fingerprint": digest,
        "language": language,
        "path": relative,
        "line": line,
        "column": column,
        "symbol": symbol,
        "syntactic_class": syntactic_class,
        "source_expression": source,
        "value_source": "unresolved",
        "dp_source": "unresolved",
        "dependencies": [],
        "portability": "unresolved",
        "nonfinite_reachability": "unknown",
        "applicability": "unresolved",
        "exception": "unresolved",
        "review_status": "evidence_only",
    }
    if unsupported_reason is None:
        formatting = []
        if NUMERIC_HINT.search(source):
            formatting.append("same_sink_numeric_or_formatting_lexeme")
        base.update(
            {
                "numeric_relevance": (
                    "candidate" if formatting else "not_statically_apparent"
                ),
                "formatting_evidence": formatting,
            },
        )
    else:
        base["unsupported_reason"] = unsupported_reason
    return base


def _enclosing_symbols(tree: ast.AST) -> dict[ast.AST, str]:
    parents = {
        child: parent
        for parent in ast.walk(tree)
        for child in ast.iter_child_nodes(parent)
    }
    result: dict[ast.AST, str] = {}
    for node in ast.walk(tree):
        parent = node
        symbol = "<module>"
        while parent in parents:
            parent = parents[parent]
            if isinstance(parent, (ast.FunctionDef, ast.AsyncFunctionDef)):
                symbol = parent.name
                break
        result[node] = symbol
    return result


def _python_bindings(tree: ast.AST) -> tuple[set[str], set[str]]:
    streams: set[str] = set()
    csv_writers: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Assign, ast.AnnAssign)):
            targets = (
                node.targets if isinstance(node, ast.Assign) else [node.target]
            )
            value = node.value
            names = {
                target.id for target in targets if isinstance(target, ast.Name)
            }
            if isinstance(value, ast.Call):
                if (
                    isinstance(value.func, ast.Name)
                    and value.func.id == "open"
                ):
                    streams.update(names)
                if (
                    isinstance(value.func, ast.Attribute)
                    and isinstance(value.func.value, ast.Name)
                    and value.func.value.id == "csv"
                    and value.func.attr == "writer"
                ):
                    csv_writers.update(names)
        if isinstance(node, ast.withitem) and isinstance(
            node.optional_vars,
            ast.Name,
        ):
            context = node.context_expr
            if isinstance(context, ast.Call) and (
                (
                    isinstance(context.func, ast.Name)
                    and context.func.id == "open"
                )
                or (
                    isinstance(context.func, ast.Attribute)
                    and context.func.attr == "open"
                )
            ):
                streams.add(node.optional_vars.id)
    return streams, csv_writers


def python_inventory(
    path: Path,
    relative: str,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """
    Enumerate recognized Python sinks and unsupported review candidates.
    """

    source = path.read_text(encoding="utf-8")
    try:
        tree = ast.parse(source)
    except SyntaxError as error:
        return [], [
            _record(
                relative=relative,
                language="python",
                line=error.lineno or 1,
                column=error.offset or 0,
                symbol="<module>",
                syntactic_class="python.syntax_error",
                source=error.msg,
                unsupported_reason="otherwise_unsupported_construct",
            ),
        ]
    symbols = _enclosing_symbols(tree)
    streams, csv_writers = _python_bindings(tree)
    sites: list[dict[str, Any]] = []
    reviews: list[dict[str, Any]] = []
    lines = source.splitlines()
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        expression = lines[node.lineno - 1].strip()
        kind: str | None = None
        unsupported: str | None = None
        func = node.func
        if isinstance(func, ast.Name) and func.id == "print":
            kind = "python.direct_print"
        elif (
            isinstance(func, ast.Attribute)
            and isinstance(func.value, ast.Attribute)
            and isinstance(func.value.value, ast.Name)
            and func.value.value.id == "sys"
            and func.value.attr in {"stdout", "stderr"}
            and func.attr in {"write", "writelines"}
        ):
            kind = "python.sys_stream_write"
        elif (
            isinstance(func, ast.Attribute)
            and isinstance(func.value, ast.Name)
            and func.value.id in streams
            and func.attr in {"write", "writelines"}
        ):
            kind = "python.bound_text_stream_write"
        elif (
            isinstance(func, ast.Attribute)
            and isinstance(func.value, ast.Name)
            and func.value.id in csv_writers
            and func.attr in {"writerow", "writerows"}
        ):
            kind = "python.csv_writer_row"
        elif (
            isinstance(func, ast.Attribute)
            and isinstance(func.value, ast.Name)
            and func.value.id == "json"
            and func.attr == "dump"
        ):
            kind = "python.json_dump"
        elif isinstance(func, ast.Attribute) and func.attr == "write_text":
            kind = "python.path_write_text"
        elif isinstance(func, ast.Name) and func.id in {"eval", "exec"}:
            unsupported = "eval_or_eval_like_construct"
        elif isinstance(func, (ast.Subscript, ast.Call)):
            unsupported = "dynamic_or_indirect_call_target"
        elif isinstance(func, ast.Attribute) and func.attr in {
            "write",
            "writelines",
            "writerow",
            "writerows",
        }:
            unsupported = "unrecognized_writer_factory"
        elif (
            isinstance(func, ast.Name)
            and func.id
            not in {
                "abs",
                "all",
                "any",
                "bool",
                "dict",
                "enumerate",
                "float",
                "int",
                "len",
                "list",
                "max",
                "min",
                "range",
                "round",
                "set",
                "sorted",
                "str",
                "sum",
                "tuple",
                "zip",
            }
            and NUMERIC_HINT.search(expression)
        ):
            unsupported = "function_or_command_variable_dispatch"

        if kind is not None:
            sites.append(
                _record(
                    relative=relative,
                    language="python",
                    line=node.lineno,
                    column=node.col_offset,
                    symbol=symbols[node],
                    syntactic_class=kind,
                    source=expression,
                ),
            )
        elif unsupported is not None:
            reviews.append(
                _record(
                    relative=relative,
                    language="python",
                    line=node.lineno,
                    column=node.col_offset,
                    symbol=symbols[node],
                    syntactic_class="python.unsupported",
                    source=expression,
                    unsupported_reason=unsupported,
                ),
            )
    return sites, reviews


def _shell_symbol(lines: list[str], index: int) -> str:
    for line in reversed(lines[: index + 1]):
        match = re.match(
            r"\s*(?:function\s+)?([A-Za-z_][A-Za-z0-9_]*)"
            r"\s*(?:\(\))?\s*\{",
            line,
        )
        if match:
            return match.group(1)
    return "<file>"


def shell_inventory(
    path: Path,
    relative: str,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """
    Enumerate bounded direct Shell sinks and unsupported constructs.
    """

    lines = path.read_text(encoding="utf-8").splitlines()
    sites: list[dict[str, Any]] = []
    reviews: list[dict[str, Any]] = []
    heredoc = re.compile(
        r"^\s*(?:cat|tee)\b.*?(?:<<-?|<<<)\s*['\"]?[A-Za-z0-9_]+",
    )
    direct = re.compile(
        r"^\s*(?:builtin\s+)?(?P<command>printf|echo)\b",
    )
    for index, line in enumerate(lines):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        symbol = _shell_symbol(lines, index)
        kind: str | None = None
        unsupported: str | None = None
        match = direct.match(line)
        if heredoc.search(line):
            kind = "shell.literal_heredoc_or_herestring"
        elif match:
            kind = (
                "shell.direct_printf_or_builtin_printf"
                if match.group("command") == "printf"
                else "shell.direct_echo"
            )
            if re.search(r"(?:^|\s)(?:>{1,2}|[0-9]+>{1,2})", line):
                kind = "shell.static_redirection"
            if "$(" in line or "`" in line:
                unsupported = "unresolved_command_substitution"
        elif re.search(r"\b(?:eval|source)\b|^\s*\.\s+", line):
            unsupported = (
                "eval_or_eval_like_construct"
                if re.search(r"\beval\b", line)
                else "sourced_code_outside_scanned_statement"
            )
        elif re.search(r"\b(?:awk|perl|python|python3|Rscript)\b", line):
            unsupported = "embedded_language"
        elif re.search(r"\$\{?[A-Za-z_][A-Za-z0-9_]*\}?\s+", line):
            unsupported = "function_or_command_variable_dispatch"
        elif re.search(r"\bexec\s+[0-9]*>&|>&\$\{?", line):
            unsupported = "dynamic_file_descriptor_or_stream_alias"

        if kind is not None:
            sites.append(
                _record(
                    relative=relative,
                    language="shell",
                    line=index + 1,
                    column=max(0, len(line) - len(line.lstrip())),
                    symbol=symbol,
                    syntactic_class=kind,
                    source=stripped,
                ),
            )
        if unsupported is not None:
            reviews.append(
                _record(
                    relative=relative,
                    language="shell",
                    line=index + 1,
                    column=max(0, len(line) - len(line.lstrip())),
                    symbol=symbol,
                    syntactic_class="shell.unsupported",
                    source=stripped,
                    unsupported_reason=unsupported,
                ),
            )
    return sites, reviews


def inventory(root: Path, languages: set[str]) -> dict[str, Any]:
    """
    Return a complete inventory only within recognized syntactic classes.
    """

    sites: list[dict[str, Any]] = []
    reviews: list[dict[str, Any]] = []
    for path in maintained_paths(root, languages):
        relative = path.relative_to(root).as_posix()
        if path.suffix == ".py":
            found, candidates = python_inventory(path, relative)
        else:
            found, candidates = shell_inventory(path, relative)
        sites.extend(found)
        reviews.extend(candidates)

    def key(item: dict[str, Any]) -> tuple[object, ...]:
        return (
            item["path"],
            item["line"],
            item["column"],
            item["fingerprint"],
        )

    return {
        "schema_version": 1,
        "inventory_scope": {
            "claim": (
                "complete_only_within_recognized_static_syntactic_classes"
            ),
            "python_syntactic_classes": PYTHON_CLASSES,
            "shell_syntactic_classes": SHELL_CLASSES,
            "unsupported_constructs": UNSUPPORTED,
        },
        "sites": sorted(sites, key=key),
        "review_candidates": sorted(reviews, key=key),
    }


def main(argv: list[str] | None = None) -> int:
    """
    Emit the read-only inventory to stdout or one explicit artifact path.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--languages",
        nargs="+",
        choices=("python", "shell"),
        default=("python", "shell"),
    )
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args(argv)
    result = inventory(args.root.resolve(), set(args.languages))
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.json_out is None:
        print(rendered, end="")
    else:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(rendered, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
