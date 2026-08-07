#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: help_option_order.py
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
Validate reviewed option order against actual registered help surfaces.
"""

from __future__ import annotations

import argparse
import ast
import contextlib
import inspect
import io
import json
import os
import re
import runpy
import subprocess
import sys
from pathlib import Path
from typing import Any

PARAMETER_ROW = re.compile(r"^\s+(?P<aliases>-[^:]+?)\s+:\s+", re.MULTILINE)
OPTION_TOKEN = re.compile(r"--[A-Za-z][A-Za-z0-9_-]*")
SECTION = re.compile(
    r"(?ms)^(?P<name>[A-Z][A-Za-z ]+)\n-+\n"
    r"(?P<body>.*?)(?=^[A-Z][A-Za-z ]+\n-+\n|\Z)",
)


def _logical(alias: str) -> str:
    return alias.lstrip("-").replace("-", "_")


def _deduplicate(values: list[str]) -> list[str]:
    return list(dict.fromkeys(values))


def python_parser_order(text: str) -> list[str]:
    """
    Extract actual visible argparse registration order from source.
    """

    tree = ast.parse(text)
    calls = [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr == "add_argument"
    ]
    calls.sort(key=lambda node: (node.lineno, node.col_offset))
    result = ["help"] if "add_help_cap(" in text else []

    for node in calls:
        aliases = [
            value.value
            for value in node.args
            if isinstance(value, ast.Constant)
            and isinstance(value.value, str)
            and value.value.startswith("--")
        ]

        if aliases:
            result.append(_logical(aliases[-1]))

    return _deduplicate(result)


def rendered_section(text: str, names: set[str]) -> str:
    """
    Return one actual rendered help section body.
    """

    matches = [
        match.group("body")
        for match in SECTION.finditer(text)
        if match.group("name") in names
    ]

    return "\n".join(matches)


def rendered_usage_order(text: str) -> list[str]:
    """
    Extract exact logical option order from rendered Usage.
    """

    return _deduplicate(
        [
            _logical(token)
            for token in OPTION_TOKEN.findall(
                rendered_section(text, {"Usage"}),
            )
        ],
    )


def rendered_usage_rows(text: str) -> list[list[str]]:
    """
    Extract logical option groups from physical Usage continuation rows.
    """

    rows: list[list[str]] = []

    for line in rendered_section(text, {"Usage"}).splitlines():
        if not line.startswith("    "):
            continue

        members = _deduplicate(
            [_logical(token) for token in OPTION_TOKEN.findall(line)],
        )

        if members:
            rows.append(members)

    return rows


def rendered_parameter_order(text: str) -> list[str]:
    """
    Extract exact logical option order from rendered Options/Parameters rows.
    """

    body = rendered_section(text, {"Options", "Parameters"})
    result: list[str] = []

    for line in body.splitlines():
        if not re.match(r"^  -(?!\s)", line):
            continue

        aliases = OPTION_TOKEN.findall(line)

        if aliases:
            result.append(_logical(aliases[-1]))

    return _deduplicate(result)


def rendered_alias_map(text: str) -> dict[str, str]:
    """
    Map every displayed long alias to its row's canonical logical option.
    """

    body = rendered_section(text, {"Options", "Parameters"})
    result: dict[str, str] = {}

    for line in body.splitlines():
        if not re.match(r"^  -(?!\s)", line):
            continue

        aliases = OPTION_TOKEN.findall(line)
        if not aliases:
            continue

        canonical = _logical(aliases[-1])

        for alias in aliases:
            result[_logical(alias)] = canonical

    return result


def _shell_pattern_logical(pattern: str) -> str:
    """
    Normalize one bounded Shell case-pattern long option.
    """

    value = pattern.lstrip("-")
    value = value.replace("[_-]", "_")
    value = value.replace("[e]?", "e")

    return value.replace("-", "_")


def shell_function_text(text: str, symbol: str) -> str:
    """
    Return the bounded body of one named Shell function.
    """

    match = re.search(
        rf"(?ms)^function {re.escape(symbol)}\(\) \{{\n"
        r"(?P<body>.*?)^}\n",
        text,
    )
    if match is None:
        raise RuntimeError(f"Shell function not found: {symbol}")

    return match.group("body")


def shell_parser_order(
    text: str,
    alias_map: dict[str, str],
) -> list[str]:
    """
    Extract ordered public logical options from bounded Shell parse_args.
    """

    result: list[str] = []

    for line in shell_function_text(text, "parse_args").splitlines():
        match = re.match(r"^\s*(?P<patterns>-[^)]*)\)\s*$", line)
        if match is None:
            continue

        for value in re.findall(
            r"--[A-Za-z][A-Za-z0-9_\[\]?-]*",
            match.group("patterns"),
        ):
            logical = _shell_pattern_logical(value)
            result.append(alias_map.get(logical, logical))

    if re.search(r"\^\(-h\|--h\[e\]\?lp\)\$", text):
        result.insert(0, "help")

    return _deduplicate(result)


def shell_reporting_order(
    text: str,
    symbol: str,
    logical_options: set[str],
) -> list[str]:
    """
    Extract registered argument rows from one bounded Shell reporter.
    """

    result: list[str] = []

    for line in shell_function_text(text, symbol).splitlines():
        match = re.match(
            r'^\s*echo "(?P<label>[A-Za-z][A-Za-z0-9_]*)=',
            line,
        )

        if match is not None and match.group("label") in logical_options:
            result.append(match.group("label"))

    return result


def render_help(root: Path, command: list[str]) -> str:
    """
    Execute one approved read-only help command.
    """

    resolved = [
        sys.executable if item == "{python}" else item for item in command
    ]
    environment = {
        **os.environ,
        "PYTHONDONTWRITEBYTECODE": "1",
    }

    result = subprocess.run(
        resolved,
        cwd=root,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )

    if result.returncode:
        raise RuntimeError(
            f"help command failed ({result.returncode}): {resolved!r}",
        )

    return result.stdout + result.stderr


def python_reporting_order(
    root: Path,
    parser_path: Path,
    symbol: str,
    probes: list[list[str]],
    call_arguments: list[dict[str, Any]],
) -> list[str]:
    """
    Probe the actual Python verbose argument reporter without running work.
    """

    source = (root / parser_path).read_text(encoding="utf-8")
    tree = ast.parse(source)
    function = next(
        node
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name == symbol
    )
    static: list[str] = []

    for node in sorted(
        (
            node
            for node in ast.walk(function)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "print"
            and node.args
        ),
        key=lambda node: (node.lineno, node.col_offset),
    ):
        argument = node.args[0]
        prefix = ""

        if isinstance(argument, ast.Constant) and isinstance(
            argument.value,
            str,
        ):
            prefix = argument.value
        elif isinstance(argument, ast.JoinedStr) and argument.values:
            first = argument.values[0]

            if isinstance(first, ast.Constant) and isinstance(
                first.value, str
            ):
                prefix = first.value

        tokens = OPTION_TOKEN.findall(prefix)

        if tokens:
            static.append(_logical(tokens[0]))

    static = _deduplicate(static)

    namespace = runpy.run_path(
        str(root / parser_path),
        run_name="phase_c_help_option_order_probe",
    )
    parse_args = namespace["parse_args"]
    reporter = namespace[symbol]
    if not callable(reporter):
        raise RuntimeError(f"verbose reporter is not callable: {symbol}")

    signature = inspect.signature(reporter)
    binding_kinds: list[tuple[str, Any]] = []
    preflight: list[Any] = []

    for binding in call_arguments:
        if not isinstance(binding, dict):
            raise RuntimeError("verbose reporter binding must be an object")

        if set(binding) == {"source"} and binding["source"] == "parsed_args":
            binding_kinds.append(("parsed_args", None))
            preflight.append(object())

            continue

        if set(binding) == {"literal"}:
            try:
                literal = json.loads(
                    json.dumps(binding["literal"], allow_nan=False),
                )
            except (TypeError, ValueError) as error:
                raise RuntimeError(
                    "verbose reporter literal must be JSON-only data",
                ) from error

            binding_kinds.append(("literal", literal))
            preflight.append(literal)

            continue

        raise RuntimeError(
            "verbose reporter binding is not a closed data binding",
        )

    try:
        signature.bind(*preflight)
    except TypeError as error:
        raise RuntimeError(
            f"verbose reporter bindings do not match {symbol}{signature}",
        ) from error

    observed: list[list[str]] = []

    for arguments in probes:
        parsed = parse_args(arguments)
        resolved = [
            parsed if kind == "parsed_args" else value
            for kind, value in binding_kinds
        ]
        stream = io.StringIO()

        with contextlib.redirect_stderr(stream):
            reporter(*resolved)

        sequence = [
            _logical(token)
            for token in OPTION_TOKEN.findall(stream.getvalue())
        ]
        positions = [static.index(item) for item in sequence]
        if positions != sorted(positions):
            raise RuntimeError(
                "verbose probe order differs from reporting source order",
            )

        observed.append(sequence)

    emitted = {item for sequence in observed for item in sequence}

    return [item for item in static if item in emitted]


def actual_surfaces(
    root: Path,
    record: dict[str, Any],
    surface: dict[str, Any],
    overrides: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """
    Collect actual parser, rendered, and reporting surfaces.
    """

    adapters = record["adapters"]
    overrides = overrides or {}
    parser_path = Path(adapters["parser_path"])
    parser_text = overrides.get("parser_text")

    if parser_text is None:
        parser_text = (root / parser_path).read_text(encoding="utf-8")

    rendered = overrides.get("rendered_help")

    if rendered is None:
        rendered = render_help(root, adapters["help_command"])

    usage = rendered_usage_order(rendered)
    usage_rows = rendered_usage_rows(rendered)
    parameters = rendered_parameter_order(rendered)

    if record["language"] == "python":
        parser = python_parser_order(parser_text)
    else:
        alias_map = rendered_alias_map(rendered)

        if adapters["parser"] == "shell_bounded_parse_args_order":
            parser = shell_parser_order(parser_text, alias_map)
        else:
            accepted = set(shell_parser_order(parser_text, alias_map))
            expected_set = set(record["logical_order"])
            parser = (
                parameters if accepted == expected_set else sorted(accepted)
            )

    reporting_expected = record["surface_orders"]["reporting_order"]

    if isinstance(reporting_expected, dict):
        reporting: list[str] | dict[str, str] = reporting_expected
    elif "reporting_texts" in overrides:
        reporting = _deduplicate(
            [
                _logical(token)
                for text in overrides["reporting_texts"]
                for token in OPTION_TOKEN.findall(text)
            ],
        )
    elif record["language"] == "shell":
        reporting = shell_reporting_order(
            overrides.get("reporting_text", parser_text),
            adapters["reporting_symbol"],
            set(reporting_expected),
        )
    else:
        reporting = python_reporting_order(
            root,
            parser_path,
            adapters["reporting_symbol"],
            adapters["reporting_probes"],
            adapters["reporting_call_arguments"],
        )

    return {
        "parser_registration": parser,
        "usage": usage,
        "usage_rows": usage_rows,
        "options_parameters": parameters,
        "reporting_order": reporting,
    }


def validate_order(
    root: Path,
    data: dict[str, Any],
    surface_overrides: dict[str, dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    """
    Return role, relationship, and actual-surface parity findings.
    """

    root = root.resolve()
    surfaces = {item["id"]: item for item in data["surfaces"]}
    category_rank = {
        category: position
        for position, category in enumerate(
            data["option_category_precedence"],
        )
    }
    findings: list[dict[str, Any]] = []

    def add(rule_id: str, surface: str, message: str) -> None:
        findings.append(
            {"rule_id": rule_id, "surface_id": surface, "message": message},
        )

    for record in data["option_realizations"]:
        surface_id = record["surface_id"]
        order = record["logical_order"]
        roles = record["roles"]

        if len(order) != len(set(order)):
            add(
                "HELP.OPTION.ORDER.SURFACE_PARITY",
                surface_id,
                "logical order contains duplicates",
            )

        if not record.get("roles_reviewed"):
            add(
                "HELP.OPTION.ORDER.ROLE_UNREVIEWED",
                surface_id,
                "semantic roles require explicit review",
            )

        if set(roles) != set(order) or any(
            role not in category_rank for role in roles.values()
        ):
            add(
                "HELP.OPTION.ORDER.ROLE_UNREVIEWED",
                surface_id,
                "every option requires one explicit recognized role",
            )
        else:
            observed_categories = [
                roles[option]
                for index, option in enumerate(order)
                if index == 0 or roles[option] != roles[order[index - 1]]
            ]
            ranks = [
                category_rank[category] for category in observed_categories
            ]

            if ranks != sorted(ranks):
                add(
                    "HELP.OPTION.ORDER.CATEGORY",
                    surface_id,
                    "category order differs without a reviewed deviation",
                )

        for relationship in record.get("relationships", []):
            members = relationship["members"]

            try:
                positions = [order.index(member) for member in members]
            except ValueError:
                add(
                    "HELP.OPTION.ORDER.GROUP",
                    surface_id,
                    "relationship references an absent option",
                )

                continue

            kind = relationship["kind"]

            if kind == "adjacent" and positions != list(
                range(positions[0], positions[0] + len(positions)),
            ):
                add(
                    "HELP.OPTION.ORDER.ADJACENCY",
                    surface_id,
                    "adjacent members are not immediately consecutive",
                )
            elif kind == "contiguous" and positions != list(
                range(min(positions), max(positions) + 1),
            ):
                add(
                    "HELP.OPTION.ORDER.CONTIGUITY",
                    surface_id,
                    "semantic-group members are not contiguous",
                )
            elif (
                kind == "same_group"
                and len({roles.get(member) for member in members}) != 1
            ):
                add(
                    "HELP.OPTION.ORDER.GROUP",
                    surface_id,
                    "same-group members do not share one reviewed role",
                )

        expected_usage_rows = [row["members"] for row in record["usage_rows"]]

        overrides = (surface_overrides or {}).get(surface_id)

        try:
            actual = actual_surfaces(
                root,
                record,
                surfaces[surface_id],
                overrides,
            )
        except (KeyError, OSError, RuntimeError, SystemExit) as error:
            add(
                "HELP.OPTION.ORDER.SURFACE_PARITY",
                surface_id,
                f"actual surface probe failed: {error}",
            )

            continue

        for name, expected in record["surface_orders"].items():
            observed = actual[name]

            if isinstance(expected, dict):
                if observed != expected:
                    add(
                        "HELP.OPTION.ORDER.SURFACE_PARITY",
                        surface_id,
                        f"{name} non-applicability differs",
                    )

                continue

            if observed != expected:
                add(
                    "HELP.OPTION.ORDER.SURFACE_PARITY",
                    surface_id,
                    f"{name} differs: expected={expected!r}; "
                    f"actual={observed!r}",
                )

        if actual["usage_rows"] != expected_usage_rows:
            add(
                "HELP.OPTION.ORDER.GROUP",
                surface_id,
                "usage rows differ: expected="
                f"{expected_usage_rows!r}; actual={actual['usage_rows']!r}",
            )

    return findings


def main(argv: list[str] | None = None) -> int:
    """
    Run the read-only option-order checker.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("dev/config/help_contracts.json"),
    )
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args(argv)
    root = args.root.resolve()
    config = args.config if args.config.is_absolute() else root / args.config

    findings = validate_order(
        root,
        json.loads(config.read_text(encoding="utf-8")),
    )

    if args.json:
        print(json.dumps(findings, indent=2, sort_keys=True))
    else:
        for finding in findings:
            print(
                f"{finding['rule_id']}: {finding['surface_id']}: "
                f"{finding['message']}",
            )

    return 1 if findings else 0


if __name__ == "__main__":
    raise SystemExit(main())
