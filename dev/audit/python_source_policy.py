#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: python_source_policy.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


"""
Check deterministic portions of the bounded Python source-policy pilot.
"""

from __future__ import annotations

import argparse
import ast
import dataclasses
import hashlib
import io
import json
import re
import tokenize
from collections import Counter
from collections.abc import Iterable, Iterator
from itertools import pairwise
from pathlib import Path

from dev.audit.python_naming_vocabulary import prohibited_internal_segments
from dev.audit.python_version_policy import maintained_python_paths

VERSION = "1.10"
DEFAULT_CONFIG = (
    Path(__file__).resolve().parents[1]
    / "config"
    / "pilots"
    / "python_source_policy.json"
)

RULE_DOC_LAYOUT = "PY.DOCSTRING.LAYOUT"
RULE_DOC_NUMPY = "PY.DOCSTRING.NUMPY"
RULE_ANNOTATIONS = "PY.TYPE.ANNOTATIONS"
RULE_STRINGS = "PY.STRING.QUOTES"
RULE_CLI_HELP_LAYOUT = "PY.CLI.HELP.LAYOUT"
RULE_HELP_SENTENCES = "HELP.PROSE.SENTENCES"
RULE_COMMENTS = "PY.COMMENT.FORM"
RULE_NAMING = "PY.NAMING.IDENTIFIERS"
RULE_TOPOLOGY = "PY.SOURCE.LAYOUT"
RULE_MULTILINE = "SOURCE.DELIMITED.MULTILINE"
RULE_LINE_LENGTH = "PY.FORMAT.LINE_LENGTH"
RULE_IDS = (
    RULE_DOC_LAYOUT,
    RULE_DOC_NUMPY,
    RULE_ANNOTATIONS,
    RULE_STRINGS,
    RULE_CLI_HELP_LAYOUT,
    RULE_HELP_SENTENCES,
    RULE_COMMENTS,
    RULE_NAMING,
    RULE_TOPOLOGY,
    RULE_MULTILINE,
)

NUMPY_SECTION_ORDER = (
    "Parameters",
    "Returns",
    "Yields",
    "Raises",
    "Warns",
    "Warnings",
    "See Also",
    "Notes",
    "References",
    "Examples",
)
NUMPY_SECTION_INDEX = {
    name: index for index, name in enumerate(NUMPY_SECTION_ORDER)
}
PSEUDO_SECTION = re.compile(
    r"^(?:Exits|Side effects|Yields):$",
    re.IGNORECASE,
)
SECTION_UNDERLINE = re.compile(r"^-+$")
NUMPY_TYPED_ENTRY = re.compile(
    r"^(?P<names>\*{0,2}[A-Za-z_][A-Za-z0-9_]*(?:,\s*"
    r"\*{0,2}[A-Za-z_][A-Za-z0-9_]*)*)\s+:\s+(?P<type>\S.*)$",
)
NUMPY_TYPE_CLOSER = re.compile(
    r"^[\]\)}]+(?:\s*(?:\||,)\s*\S.*)?$",
)

STRING_FORM = re.compile(
    r"(?is)^(?P<prefix>[rubf]*)(?P<quote>'''|\"\"\"|'|\")",
)
SNAKE_CASE = re.compile(r"^_*[a-z][a-z0-9]*(?:_[a-z0-9]+)*_?$")
UPPER_CAMEL = re.compile(r"^[A-Z][A-Za-z0-9]*$")
DIRECTIVE = re.compile(
    r"#\s*(?:fmt:|isort:|mypy:|noqa\b|nosec\b|pragma:|"
    r"pyright:|pylint:|ruff:|type:)",
)
COMMENT_LITERAL_LABEL = re.compile(r"^(?:Example|Examples|Usage):(?:\s|$)")
COMMENT_CODE_FRAGMENT = re.compile(
    r"^(?:"
    r"[A-Za-z_][A-Za-z0-9_.]*(?:\[[^\]]+\])?\s*"
    r"(?:=|==|!=|<=|>=|<|>|\+=|-=|\*=|/=)"
    r"|[-+*/%<>=|&^~]+(?:\s|$))",
)
OPAQUE_IDENTIFIER_SEGMENTS = prohibited_internal_segments()
DEFINITION = (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)
CALLABLE = (ast.FunctionDef, ast.AsyncFunctionDef)
TRANSFER = (ast.Return, ast.Raise, ast.Break, ast.Continue)
DELIMITED = (ast.Call, ast.List, ast.Tuple, ast.Dict, ast.Set)


@dataclasses.dataclass(frozen=True, order=True)
class Finding:
    """
    Represent one deterministic Python source-policy finding.
    """

    path: str
    line: int
    column: int
    rule_id: str
    message: str

    def format(self) -> str:
        """
        Render one stable line-oriented finding.
        """

        return (
            f"{self.rule_id}: {self.path}:{self.line}:{self.column}: "
            f"{self.message}"
        )


@dataclasses.dataclass(frozen=True)
class Analysis:
    """
    Hold deterministic findings and stable per-file source facts.
    """

    findings: tuple[Finding, ...]
    facts: dict[str, object]


def source_fingerprint(text: str) -> str:
    """
    Return one stable SHA-256 identity for source text.
    """

    return f"sha256:{hashlib.sha256(text.encode('utf-8')).hexdigest()}"


def _tokens(text: str) -> list[tokenize.TokenInfo]:
    """
    Return Python tokens, or an empty list for un-tokenizable source.
    """

    try:
        return list(tokenize.generate_tokens(io.StringIO(text).readline))
    except (IndentationError, tokenize.TokenError):
        return []


def _object_name(
    node: ast.AST,
    parents: dict[ast.AST, ast.AST],
) -> str:
    """
    Return one dotted lexical object name.
    """

    names: list[str] = []
    current: ast.AST | None = node

    while current is not None:
        if isinstance(current, ast.Module):
            names.append("<module>")

            break

        name = getattr(current, "name", None)

        if isinstance(name, str):
            names.append(name)

        current = parents.get(current)

    return ".".join(reversed(names))


def _parents(tree: ast.AST) -> dict[ast.AST, ast.AST]:
    """
    Return a parent map for one AST.
    """

    return {
        child: parent
        for parent in ast.walk(tree)
        for child in ast.iter_child_nodes(parent)
    }


def _doc_objects(tree: ast.Module) -> list[ast.AST]:
    """
    Return module, class, and callable objects in source order.
    """

    return sorted(
        (
            node
            for node in ast.walk(tree)
            if isinstance(node, (ast.Module, *DEFINITION))
        ),
        key=lambda node: (
            getattr(node, "lineno", 0),
            getattr(node, "col_offset", 0),
        ),
    )


def _doc_token(
    node: ast.AST,
    tokens: list[tokenize.TokenInfo],
) -> tuple[ast.Expr, tuple[tokenize.TokenInfo, ...]] | None:
    """
    Return an actual first-statement docstring and its source tokens.
    """

    body = getattr(node, "body", None)
    if (
        not body
        or not isinstance(body[0], ast.Expr)
        or not isinstance(body[0].value, ast.Constant)
        or not isinstance(body[0].value.value, str)
    ):
        return None

    expression = body[0]
    string_tokens = tuple(
        token
        for token in tokens
        if token.type == tokenize.STRING
        and _position_within(token, expression)
    )

    return (expression, string_tokens) if string_tokens else None


def _finding(
    findings: list[Finding],
    path: str,
    token: tokenize.TokenInfo,
    rule_id: str,
    message: str,
) -> None:
    """
    Append one token-anchored finding.
    """

    findings.append(
        Finding(
            path,
            token.start[0],
            token.start[1] + 1,
            rule_id,
            message,
        ),
    )


def _finding_at(
    findings: list[Finding],
    path: str,
    line: int,
    column: int,
    rule_id: str,
    message: str,
) -> None:
    """
    Append one source-position finding.
    """

    findings.append(
        Finding(
            path,
            line,
            column,
            rule_id,
            message,
        ),
    )


def _indent_width(line: str) -> int:
    """
    Return the expanded indentation width of one physical line.
    """

    prefix = line[: len(line) - len(line.lstrip(" \t"))]
    return len(prefix.expandtabs(8))


def _numpy_doc_body(
    node: ast.AST,
    tokens: list[tokenize.TokenInfo],
) -> tuple[tokenize.TokenInfo, int, list[str]] | None:
    """
    Return one canonical multiline NumPy-docstring body.

    Parameters
    ----------
    node : ast.AST
        Candidate module, class, or callable.
    tokens : list[tokenize.TokenInfo]
        Token stream retaining physical string boundaries.

    Returns
    -------
    doc_body : tuple[tokenize.TokenInfo, int, list[str]] | None
        String token, expected content indentation, and physical body lines,
        or 'None' when the candidate uses another recognized form.
    """

    located = _doc_token(node, tokens)
    if located is None:
        return None

    _, doc_tokens = located
    if len(doc_tokens) != 1:
        return None

    token = doc_tokens[0]
    form = STRING_FORM.match(token.string)
    if form is None or form.group("quote") != '"""':
        return None

    physical_lines = token.string.splitlines()
    if len(physical_lines) < 3:
        return None

    expected_indent = (
        0 if isinstance(node, ast.Module) else node.col_offset + 4
    )

    return token, expected_indent, physical_lines[1:-1]


def _numpy_section_rows(
    body: list[str],
    token: tokenize.TokenInfo,
    expected_indent: int,
    path: str,
    findings: list[Finding],
    counts: Counter[str],
) -> list[tuple[int, str]]:
    """
    Validate NumPy section headings and return their body indexes.

    Parameters
    ----------
    body : list[str]
        Physical docstring body lines.
    token : tokenize.TokenInfo
        Owning docstring token.
    expected_indent : int
        Required heading and body indentation.
    path : str
        Repository-relative diagnostic path.
    findings : list[Finding]
        Mutable deterministic finding collection.
    counts : Counter[str]
        Mutable NumPy-docstring fact counts.

    Returns
    -------
    rows : list[tuple[int, str]]
        Canonical section indexes and names in source order.
    """

    section_rows: list[tuple[int, str]] = []

    for index, line in enumerate(body):
        stripped = line.strip()
        source_line = token.start[0] + index + 1

        if PSEUDO_SECTION.fullmatch(stripped):
            _finding_at(
                findings,
                path,
                source_line,
                _indent_width(line) + 1,
                RULE_DOC_NUMPY,
                f"noncanonical NumPy pseudo-section '{stripped}'",
            )
            counts["pseudo_section"] += 1

            continue

        if stripped not in NUMPY_SECTION_INDEX:
            if SECTION_UNDERLINE.fullmatch(stripped):
                previous = body[index - 1].strip() if index else ""

                if previous not in NUMPY_SECTION_INDEX:
                    _finding_at(
                        findings,
                        path,
                        source_line,
                        _indent_width(line) + 1,
                        RULE_DOC_NUMPY,
                        "orphan or noncanonical NumPy section underline",
                    )
                    counts["orphan_underline"] += 1

            continue

        section_rows.append((index, stripped))

        previous = body[index - 1] if index else ""

        if previous.strip():
            _finding_at(
                findings,
                path,
                source_line,
                expected_indent + 1,
                RULE_DOC_NUMPY,
                f"NumPy section '{stripped}' requires a blank line before "
                "its heading",
            )

            counts["missing_section_boundary"] += 1

        if _indent_width(line) != expected_indent:
            _finding_at(
                findings,
                path,
                source_line,
                _indent_width(line) + 1,
                RULE_DOC_NUMPY,
                "NumPy section heading must align with docstring content",
            )
            counts["heading_indentation"] += 1

        expected_underline = " " * expected_indent + "-" * len(stripped)
        following = body[index + 1] if index + 1 < len(body) else ""

        if following != expected_underline:
            _finding_at(
                findings,
                path,
                source_line,
                expected_indent + 1,
                RULE_DOC_NUMPY,
                f"NumPy section '{stripped}' requires an exact underline",
            )
            counts["malformed_underline"] += 1

    return section_rows


def _check_numpy_typed_section_entries(
    node: ast.AST,
    name: str,
    body: list[str],
    body_start: int,
    body_end: int,
    token: tokenize.TokenInfo,
    expected_indent: int,
    path: str,
    findings: list[Finding],
    counts: Counter[str],
) -> set[str]:
    """
    Validate typed entries in one parameter, return, or yield section.

    Parameters
    ----------
    node : ast.AST
        Callable that owns the documented contract.
    name : str
        Canonical NumPy section name.
    body : list[str]
        Physical docstring body lines.
    body_start : int
        Zero-based first line after the section underline.
    body_end : int
        Exclusive zero-based section boundary.
    token : tokenize.TokenInfo
        Owning docstring token.
    expected_indent : int
        Required entry indentation.
    path : str
        Repository-relative diagnostic path.
    findings : list[Finding]
        Mutable deterministic finding collection.
    counts : Counter[str]
        Mutable NumPy-docstring fact counts.

    Returns
    -------
    documented_parameters : set[str]
        Normalized parameter identities documented by this section.
    """

    documented_parameters: set[str] = set()
    continuation_rows: set[int] = set()

    for content_index in range(body_start, body_end):
        if content_index in continuation_rows:
            continue

        line = body[content_index]

        if not line.strip() or _indent_width(line) != expected_indent:
            continue

        entry = NUMPY_TYPED_ENTRY.fullmatch(line.strip())

        if entry is None:
            _finding_at(
                findings,
                path,
                token.start[0] + content_index + 1,
                expected_indent + 1,
                RULE_DOC_NUMPY,
                f"NumPy section '{name}' entry requires 'name : type'",
            )
            counts["untyped_or_unnamed_entry"] += 1

            continue

        if name == "Parameters":
            documented_parameters.update(
                parameter.strip().lstrip("*")
                for parameter in entry.group("names").split(",")
            )

        documented_type, closing_index = _numpy_entry_type(
            body,
            content_index,
            body_end,
            expected_indent,
            entry,
        )

        if closing_index is not None:
            continuation_rows.update(
                range(content_index + 1, closing_index + 1),
            )
            closing_line = body[closing_index]

            if (
                NUMPY_TYPE_CLOSER.fullmatch(closing_line.strip())
                and _indent_width(closing_line) != expected_indent
            ):
                _finding_at(
                    findings,
                    path,
                    token.start[0] + closing_index + 1,
                    _indent_width(closing_line) + 1,
                    RULE_DOC_NUMPY,
                    "multiline NumPy type closing delimiter must align with "
                    "the entry line",
                )
                counts["type_closer_indentation"] += 1

        _check_numpy_entry_annotation(
            node,
            name,
            entry,
            documented_type,
            path,
            token.start[0] + content_index + 1,
            expected_indent,
            findings,
            counts,
        )

    return documented_parameters


def _check_missing_parameter_entries(
    node: ast.AST,
    documented_parameters: set[str],
    source_line: int,
    expected_indent: int,
    path: str,
    findings: list[Finding],
    counts: Counter[str],
) -> None:
    """
    Report annotated callable parameters omitted from a Parameters section.

    Parameters
    ----------
    node : ast.AST
        Callable that owns the documented contract.
    documented_parameters : set[str]
        Normalized names present in the section.
    source_line : int
        One-based section heading line.
    expected_indent : int
        Required section indentation.
    path : str
        Repository-relative diagnostic path.
    findings : list[Finding]
        Mutable deterministic finding collection.
    counts : Counter[str]
        Mutable NumPy-docstring fact counts.
    """

    if not isinstance(node, CALLABLE):
        return

    expected_parameters = set(
        _callable_parameter_annotations(node),
    ) - {"self", "cls"}
    missing_parameters = sorted(
        expected_parameters - documented_parameters,
    )

    if missing_parameters:
        _finding_at(
            findings,
            path,
            source_line,
            expected_indent + 1,
            RULE_DOC_NUMPY,
            "NumPy section 'Parameters' omits annotated parameter(s): "
            + ", ".join(missing_parameters),
        )
        counts["missing_parameter_entry"] += 1


def _check_numpy_section_content(
    node: ast.AST,
    body: list[str],
    section_rows: list[tuple[int, str]],
    token: tokenize.TokenInfo,
    expected_indent: int,
    path: str,
    findings: list[Finding],
    counts: Counter[str],
) -> None:
    """
    Validate NumPy section order, uniqueness, content, and indentation.

    Parameters
    ----------
    node : ast.AST
        Callable, class, or module that owns the docstring.
    body : list[str]
        Physical docstring body lines.
    section_rows : list[tuple[int, str]]
        Canonical section indexes and names in source order.
    token : tokenize.TokenInfo
        Owning docstring token.
    expected_indent : int
        Required heading and body indentation.
    path : str
        Repository-relative diagnostic path.
    findings : list[Finding]
        Mutable deterministic finding collection.
    counts : Counter[str]
        Mutable NumPy-docstring fact counts.
    """

    previous_order = -1
    seen: set[str] = set()

    for position, (index, name) in enumerate(section_rows):
        source_line = token.start[0] + index + 1

        if name in seen:
            _finding_at(
                findings,
                path,
                source_line,
                expected_indent + 1,
                RULE_DOC_NUMPY,
                f"NumPy section '{name}' must not be repeated",
            )
            counts["duplicate_section"] += 1

        seen.add(name)
        order = NUMPY_SECTION_INDEX[name]

        if order < previous_order:
            _finding_at(
                findings,
                path,
                source_line,
                expected_indent + 1,
                RULE_DOC_NUMPY,
                f"NumPy section '{name}' is out of canonical order",
            )
            counts["section_order"] += 1

        previous_order = max(previous_order, order)
        body_start = index + 2
        body_end = (
            section_rows[position + 1][0]
            if position + 1 < len(section_rows)
            else len(body)
        )
        content = [line for line in body[body_start:body_end] if line.strip()]

        if not content:
            _finding_at(
                findings,
                path,
                source_line,
                expected_indent + 1,
                RULE_DOC_NUMPY,
                f"NumPy section '{name}' requires substantive content",
            )
            counts["empty_section"] += 1

            continue

        if min(_indent_width(line) for line in content) != expected_indent:
            _finding_at(
                findings,
                path,
                token.start[0] + body_start + 1,
                expected_indent + 1,
                RULE_DOC_NUMPY,
                f"NumPy section '{name}' body is overindented",
            )
            counts["body_indentation"] += 1

            continue

        if name not in {"Parameters", "Returns", "Yields"}:
            continue

        documented_parameters = _check_numpy_typed_section_entries(
            node,
            name,
            body,
            body_start,
            body_end,
            token,
            expected_indent,
            path,
            findings,
            counts,
        )

        if name == "Parameters":
            _check_missing_parameter_entries(
                node,
                documented_parameters,
                source_line,
                expected_indent,
                path,
                findings,
                counts,
            )


def _doc_type_expression(value: str) -> str:
    """
    Remove supported default metadata from a documented type.
    """

    return re.sub(
        r"(?:\s*=\s*.+|,\s*(?:optional|default(?:\s*=\s*|\s+).+))$",
        "",
        value.strip(),
    )


def _normalized_doc_type(value: str) -> str:
    """
    Return a comparison form for one documented type expression.
    """

    without_default = _doc_type_expression(value)

    try:
        expression = ast.parse(without_default, mode="eval").body
    except SyntaxError:
        return " ".join(without_default.split())

    return ast.unparse(expression)


def _numpy_entry_type(
    body: list[str],
    content_index: int,
    body_end: int,
    expected_indent: int,
    entry: re.Match[str],
) -> tuple[str, int | None]:
    """
    Return one documented type expression and its final continuation row.
    """

    pieces = [entry.group("type")]
    final_index: int | None = None

    for next_index in range(content_index + 1, body_end):
        candidate = " ".join(pieces)

        try:
            ast.parse(_doc_type_expression(candidate), mode="eval")
        except SyntaxError:
            pass
        else:
            return candidate, final_index

        line = body[next_index]
        stripped = line.strip()
        is_closing_row = NUMPY_TYPE_CLOSER.fullmatch(stripped) is not None

        if (
            stripped
            and _indent_width(line) <= expected_indent
            and not is_closing_row
        ):
            break

        if stripped:
            pieces.append(stripped)
            final_index = next_index

    return " ".join(pieces), final_index


def _callable_parameter_annotations(node: ast.AST) -> dict[str, str]:
    """
    Return source-spelled annotations for one callable's parameters.
    """

    if not isinstance(node, CALLABLE):
        return {}

    arguments = node.args
    parameters = [
        *arguments.posonlyargs,
        *arguments.args,
        *arguments.kwonlyargs,
    ]

    if arguments.vararg is not None:
        parameters.append(arguments.vararg)

    if arguments.kwarg is not None:
        parameters.append(arguments.kwarg)

    return {
        parameter.arg: ast.unparse(parameter.annotation)
        for parameter in parameters
        if parameter.annotation is not None
    }


def _yield_annotation(node: ast.AST) -> str | None:
    """
    Return the yielded-value annotation encoded by an iterator return type.
    """

    if not isinstance(node, CALLABLE) or node.returns is None:
        return None

    annotation = node.returns
    if not isinstance(annotation, ast.Subscript):
        return None

    owner = ast.unparse(annotation.value)
    if owner not in {
        "AsyncGenerator",
        "AsyncIterable",
        "AsyncIterator",
        "Generator",
        "Iterable",
        "Iterator",
    }:
        return None

    if isinstance(annotation.slice, ast.Tuple):
        first = annotation.slice.elts[0]
    else:
        first = annotation.slice

    return ast.unparse(first)


def _direct_transfer_values(
    node: ast.AST,
    section_name: str,
) -> list[ast.expr | None]:
    """
    Return transfer values owned directly by one callable.
    """

    transfer_values: list[ast.expr | None] = []

    def visit(current: ast.AST) -> None:
        if current is not node and isinstance(
            current,
            (*CALLABLE, ast.ClassDef, ast.Lambda),
        ):
            return

        if section_name == "Returns" and isinstance(current, ast.Return):
            transfer_values.append(current.value)

            return

        if section_name == "Yields" and isinstance(
            current,
            (ast.Yield, ast.YieldFrom),
        ):
            value = current.value if isinstance(current, ast.Yield) else None
            transfer_values.append(value)

            return

        for child in ast.iter_child_nodes(current):
            visit(child)

    visit(node)

    return transfer_values


def _direct_value_names(value: ast.expr | None) -> tuple[str, ...] | None:
    """
    Return direct variable names from one transfer value.
    """

    if isinstance(value, ast.Name):
        return (value.id,)

    if isinstance(value, ast.Tuple) and all(
        isinstance(item, ast.Name) for item in value.elts
    ):
        return tuple(
            item.id for item in value.elts if isinstance(item, ast.Name)
        )

    return None


def _stable_transfer_names(
    node: ast.AST,
    section_name: str,
) -> tuple[str, ...] | None:
    """
    Return direct tuple-component names shared by every owned transfer.
    """

    transfers = _direct_transfer_values(node, section_name)
    if not transfers:
        return None

    names = [_direct_value_names(value) for value in transfers]
    if any(value is None for value in names):
        return None

    distinct_names = {value for value in names if value is not None}
    if len(distinct_names) != 1:
        return None

    stable_names = distinct_names.pop()
    if len(stable_names) < 2:
        return None

    return stable_names


def _check_numpy_entry_names(
    node: ast.AST,
    section_name: str,
    entry: re.Match[str],
    path: str,
    source_line: int,
    expected_indent: int,
    findings: list[Finding],
    counts: Counter[str],
) -> None:
    """
    Compare documented tuple names with stable direct transfer names.
    """

    if section_name not in {"Returns", "Yields"}:
        return

    expected_names = _stable_transfer_names(node, section_name)
    if expected_names is None:
        return

    documented_names = tuple(
        name.strip().lstrip("*") for name in entry.group("names").split(",")
    )
    if documented_names == expected_names:
        return

    documented = ", ".join(documented_names)
    expected = ", ".join(expected_names)
    transfer = "returned" if section_name == "Returns" else "yielded"
    _finding_at(
        findings,
        path,
        source_line,
        expected_indent + 1,
        RULE_DOC_NUMPY,
        f"NumPy section '{section_name}' names '{documented}' do not match "
        f"stable {transfer} names '{expected}'",
    )
    counts["transfer_name_mismatch"] += 1


def _check_numpy_entry_annotation(
    node: ast.AST,
    section_name: str,
    entry: re.Match[str],
    documented_type_text: str,
    path: str,
    source_line: int,
    expected_indent: int,
    findings: list[Finding],
    counts: Counter[str],
) -> None:
    """
    Compare a typed NumPy entry with its owning callable annotation.
    """

    documented_type = _normalized_doc_type(documented_type_text)
    expected_types: set[str] = set()

    if section_name == "Parameters":
        if not isinstance(node, CALLABLE):
            return

        annotations = _callable_parameter_annotations(node)
        names = [
            name.strip().lstrip("*")
            for name in entry.group("names").split(",")
        ]
        unknown_names = sorted(
            name for name in names if name not in annotations
        )

        if unknown_names:
            _finding_at(
                findings,
                path,
                source_line,
                expected_indent + 1,
                RULE_DOC_NUMPY,
                "NumPy section 'Parameters' documents unknown parameter(s): "
                + ", ".join(unknown_names),
            )
            counts["unknown_parameter_entry"] += 1

            return

        expected_types = {
            annotations[name] for name in names if name in annotations
        }

        if len(expected_types) != 1:
            return
    elif section_name == "Returns":
        if not isinstance(node, CALLABLE) or node.returns is None:
            return

        expected_types = {ast.unparse(node.returns)}
    else:
        annotation = _yield_annotation(node)
        if annotation is None:
            return

        expected_types = {annotation}

    expected_type = next(iter(expected_types))

    if documented_type == _normalized_doc_type(expected_type):
        _check_numpy_entry_names(
            node,
            section_name,
            entry,
            path,
            source_line,
            expected_indent,
            findings,
            counts,
        )
    else:
        _finding_at(
            findings,
            path,
            source_line,
            expected_indent + 1,
            RULE_DOC_NUMPY,
            f"NumPy section '{section_name}' type '{documented_type}' "
            f"does not match annotation '{expected_type}'",
        )
        counts["annotation_mismatch"] += 1


def _check_numpy_docstrings(
    tree: ast.Module,
    tokens: list[tokenize.TokenInfo],
    path: str,
    findings: list[Finding],
) -> Counter[str]:
    """
    Check safely provable NumPy docstring structure.

    Parameters
    ----------
    tree : ast.Module
        Parsed Python module.
    tokens : list[tokenize.TokenInfo]
        Token stream retaining physical docstring boundaries.
    path : str
        Repository-relative diagnostic path.
    findings : list[Finding]
        Mutable deterministic finding collection.

    Returns
    -------
    counts : Counter[str]
        Structural fact counts for inspected docstrings.
    """

    counts: Counter[str] = Counter()

    for node in _doc_objects(tree):
        doc_body = _numpy_doc_body(node, tokens)
        if doc_body is None:
            continue

        token, expected_indent, body = doc_body
        counts["checked"] += 1

        if body and len(body[0]) > 79:
            _finding_at(
                findings,
                path,
                token.start[0] + 1,
                expected_indent + 1,
                RULE_DOC_NUMPY,
                "docstring summary exceeds 79 columns",
            )
            counts["overlong_summary"] += 1

        section_rows = _numpy_section_rows(
            body,
            token,
            expected_indent,
            path,
            findings,
            counts,
        )
        _check_numpy_section_content(
            node,
            body,
            section_rows,
            token,
            expected_indent,
            path,
            findings,
            counts,
        )

    return counts


def _check_docstring_summary(
    token: tokenize.TokenInfo,
    literal: str,
    expected_indent: int,
    path: str,
    findings: list[Finding],
) -> None:
    """
    Check the opening boundary, summary placement, and paragraph break.

    Parameters
    ----------
    token : tokenize.TokenInfo
        Owning docstring token.
    literal : str
        Docstring value between its delimiters.
    expected_indent : int
        Required content indentation.
    path : str
        Repository-relative diagnostic path.
    findings : list[Finding]
        Mutable deterministic finding collection.
    """

    if not literal.startswith("\n"):
        _finding(
            findings,
            path,
            token,
            RULE_DOC_LAYOUT,
            "opening docstring delimiter must occupy its own line",
        )

    content_lines = literal.splitlines()

    if len(content_lines) < 2 or not content_lines[1].strip():
        _finding(
            findings,
            path,
            token,
            RULE_DOC_LAYOUT,
            "docstring summary must begin on the line after the opener",
        )
    else:
        summary = content_lines[1]

        if summary[:expected_indent] != " " * expected_indent or (
            len(summary) > expected_indent
            and summary[expected_indent].isspace()
        ):
            _finding(
                findings,
                path,
                token,
                RULE_DOC_LAYOUT,
                "docstring summary must align with the opener",
            )

        if not summary.strip().endswith("."):
            _finding(
                findings,
                path,
                token,
                RULE_DOC_LAYOUT,
                "docstring summary must end with a full stop",
            )

    if len(content_lines) <= 2:
        return

    remaining = content_lines[2:]
    substantive = [line for line in remaining if line.strip()]

    if substantive and remaining[0].strip():
        _finding(
            findings,
            path,
            token,
            RULE_DOC_LAYOUT,
            "docstring summary must be followed by a blank line before more "
            "content",
        )


def _check_docstring_physical_lines(
    token: tokenize.TokenInfo,
    expected_indent: int,
    path: str,
    findings: list[Finding],
) -> None:
    """
    Check closing-boundary and physical indentation invariants.

    Parameters
    ----------
    token : tokenize.TokenInfo
        Owning docstring token.
    expected_indent : int
        Required content indentation.
    path : str
        Repository-relative diagnostic path.
    findings : list[Finding]
        Mutable deterministic finding collection.
    """

    physical_lines = token.string.splitlines()
    closing = physical_lines[-1] if physical_lines else ""

    if closing.strip() != '"""':
        _finding(
            findings,
            path,
            token,
            RULE_DOC_LAYOUT,
            "closing docstring delimiter must occupy its own line",
        )
    elif closing != " " * expected_indent + '"""':
        _finding(
            findings,
            path,
            token,
            RULE_DOC_LAYOUT,
            "closing docstring delimiter must align with the opener",
        )

    for content_line in physical_lines[1:-1]:
        if not content_line:
            continue

        if not content_line.strip():
            _finding(
                findings,
                path,
                token,
                RULE_DOC_LAYOUT,
                "blank docstring lines must contain no whitespace",
            )

            continue

        indent = content_line[
            : len(content_line) - len(content_line.lstrip(" \t"))
        ]

        if "\t" in indent:
            _finding(
                findings,
                path,
                token,
                RULE_DOC_LAYOUT,
                "docstring indentation must use spaces",
            )

        if len(indent.expandtabs(8)) < expected_indent:
            _finding(
                findings,
                path,
                token,
                RULE_DOC_LAYOUT,
                "docstring content must not dedent before its opener",
            )


def _check_docstring_owner_adjacency(
    node: ast.AST,
    token: tokenize.TokenInfo,
    lines: list[str],
    path: str,
    findings: list[Finding],
) -> None:
    """
    Check definition-to-docstring and docstring-to-body adjacency.

    Parameters
    ----------
    node : ast.AST
        Owning class or callable.
    token : tokenize.TokenInfo
        Owning docstring token.
    lines : list[str]
        Physical source lines.
    path : str
        Repository-relative diagnostic path.
    findings : list[Finding]
        Mutable deterministic finding collection.
    """

    if not isinstance(
        node,
        (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef),
    ):
        return

    if token.start[0] > 1:
        previous = lines[token.start[0] - 2].strip()

        if not previous.endswith(":"):
            _finding(
                findings,
                path,
                token,
                RULE_DOC_LAYOUT,
                "definition signature must be immediately followed by its "
                "docstring",
            )

    body = getattr(node, "body", [])
    if len(body) <= 1:
        return

    after = token.end[0]
    first_after = lines[after] if after < len(lines) else None
    second_after = lines[after + 1] if after + 1 < len(lines) else None

    if first_after is None or first_after.strip():
        _finding(
            findings,
            path,
            token,
            RULE_DOC_LAYOUT,
            "one blank line is required after the docstring",
        )
    elif second_after is not None and not second_after.strip():
        _finding(
            findings,
            path,
            token,
            RULE_DOC_LAYOUT,
            "exactly one blank line is allowed after the docstring",
        )


def _check_one_docstring(
    node: ast.AST,
    expression: ast.Expr,
    doc_tokens: tuple[tokenize.TokenInfo, ...],
    lines: list[str],
    path: str,
    findings: list[Finding],
    counts: Counter[str],
) -> None:
    """
    Check one present docstring's complete source-layout contract.

    Parameters
    ----------
    node : ast.AST
        Owning module, class, or callable.
    expression : ast.Expr
        First-statement string expression.
    doc_tokens : tuple[tokenize.TokenInfo, ...]
        Physical tokens contributing to the expression.
    lines : list[str]
        Physical source lines.
    path : str
        Repository-relative diagnostic path.
    findings : list[Finding]
        Mutable deterministic finding collection.
    counts : Counter[str]
        Mutable docstring fact counts.
    """

    token = doc_tokens[0]

    if len(doc_tokens) != 1:
        _finding(
            findings,
            path,
            token,
            RULE_DOC_LAYOUT,
            "docstring must be one triple-double-quoted token",
        )

    form = STRING_FORM.match(token.string)
    if form is None:
        return

    prefix = form.group("prefix")
    quote = form.group("quote")
    is_docstring = isinstance(
        expression.value,
        ast.Constant,
    ) and isinstance(expression.value.value, str)
    counts["present" if is_docstring else "invalid_prefix_value"] += 1

    if quote != '"""':
        _finding(
            findings,
            path,
            token,
            RULE_DOC_LAYOUT,
            "docstring must use triple double quotes",
        )

    if prefix not in {"", "r"}:
        _finding(
            findings,
            path,
            token,
            RULE_DOC_LAYOUT,
            "docstring prefix must be absent or lowercase raw 'r'",
        )

    if quote != '"""':
        return

    expected_indent = (
        0 if isinstance(node, ast.Module) else node.col_offset + 4
    )

    if token.start[1] != expected_indent:
        _finding(
            findings,
            path,
            token,
            RULE_DOC_LAYOUT,
            "docstring opener must align exactly one indentation level inside "
            "its owner",
        )

    literal = token.string[len(prefix) + 3 : -3]

    if prefix == "r" and "\\" not in literal:
        _finding(
            findings,
            path,
            token,
            RULE_DOC_LAYOUT,
            "raw docstring prefix requires literal backslash content",
        )

    _check_docstring_summary(
        token,
        literal,
        expected_indent,
        path,
        findings,
    )
    _check_docstring_physical_lines(
        token,
        expected_indent,
        path,
        findings,
    )
    _check_docstring_owner_adjacency(
        node,
        token,
        lines,
        path,
        findings,
    )


def _check_docstrings(
    tree: ast.Module,
    tokens: list[tokenize.TokenInfo],
    lines: list[str],
    path: str,
    findings: list[Finding],
) -> tuple[set[tuple[int, int]], Counter[str]]:
    """
    Check exact layout for every recognized docstring-like first value.

    Parameters
    ----------
    tree : ast.Module
        Parsed syntax tree for the maintained source.
    tokens : list[tokenize.TokenInfo]
        Token stream paired with the parsed source.
    lines : list[str]
        Physical source lines used for location-aware checks.
    path : str
        Repository-relative path associated with the source.
    findings : list[Finding]
        Mutable finding collection populated by the check.

    Returns
    -------
    positions, counts : tuple[set[tuple[int, int]], Counter[str]]
        Docstring token positions and measured source facts.
    """

    positions: set[tuple[int, int]] = set()
    counts: Counter[str] = Counter()

    for node in _doc_objects(tree):
        located = _doc_token(node, tokens)

        if located is None:
            counts["missing"] += 1

            continue

        expression, doc_tokens = located
        positions.update(token.start for token in doc_tokens)
        _check_one_docstring(
            node,
            expression,
            doc_tokens,
            lines,
            path,
            findings,
            counts,
        )

    return positions, counts


def _check_annotations(
    tree: ast.Module,
    parents: dict[ast.AST, ast.AST],
    path: str,
    findings: list[Finding],
) -> Counter[str]:
    """
    Check annotation presence on maintained callable boundaries.

    Parameters
    ----------
    tree : ast.Module
        Parsed module whose callables are inspected.
    parents : dict[ast.AST, ast.AST]
        Parent map used to distinguish method receivers.
    path : str
        Repository-relative path used in diagnostics.
    findings : list[Finding]
        Mutable deterministic finding collection.

    Returns
    -------
    counts : Counter[str]
        Annotation-presence facts for the inspected callables.
    """

    counts: Counter[str] = Counter()

    for node in ast.walk(tree):
        if not isinstance(node, CALLABLE):
            continue

        counts["functions"] += 1
        object_name = _object_name(node, parents)
        arguments = [
            *node.args.posonlyargs,
            *node.args.args,
            *node.args.kwonlyargs,
        ]

        if node.args.vararg is not None:
            arguments.append(node.args.vararg)

        if node.args.kwarg is not None:
            arguments.append(node.args.kwarg)

        parent = parents.get(node)
        positional = (*node.args.posonlyargs, *node.args.args)
        receiver: ast.arg | None = None

        if isinstance(parent, ast.ClassDef) and positional:
            decorators = {
                decorator.id
                if isinstance(decorator, ast.Name)
                else decorator.attr
                if isinstance(decorator, ast.Attribute)
                else ""
                for decorator in node.decorator_list
            }
            first = positional[0]

            if ("classmethod" in decorators and first.arg == "cls") or (
                "staticmethod" not in decorators
                and "classmethod" not in decorators
                and first.arg == "self"
            ):
                receiver = first

        for argument in arguments:
            if argument is receiver:
                counts["object_role_exclusions"] += 1

                continue

            if argument.annotation is None:
                counts["missing_parameters"] += 1
                findings.append(
                    Finding(
                        path,
                        argument.lineno,
                        argument.col_offset + 1,
                        RULE_ANNOTATIONS,
                        f"{object_name} parameter '{argument.arg}' "
                        "requires an annotation",
                    ),
                )
            else:
                counts["annotated_parameters"] += 1

        if node.returns is None:
            counts["missing_returns"] += 1
            findings.append(
                Finding(
                    path,
                    node.lineno,
                    node.col_offset + 1,
                    RULE_ANNOTATIONS,
                    f"{object_name} requires a return annotation",
                ),
            )
        else:
            counts["annotated_returns"] += 1

    return counts


def _string_content(token: tokenize.TokenInfo) -> tuple[str, str, str] | None:
    """
    Return prefix, quote, and lexical content for one string token.
    """

    form = STRING_FORM.match(token.string)
    if form is None:
        return None

    prefix = form.group("prefix")
    quote = form.group("quote")

    return prefix, quote, token.string[len(prefix) + len(quote) : -len(quote)]


def _check_strings(
    tokens: list[tokenize.TokenInfo],
    doc_positions: set[tuple[int, int]],
    path: str,
    findings: list[Finding],
) -> Counter[str]:
    """
    Check bounded ordinary-string quote selection.

    Parameters
    ----------
    tokens : list[tokenize.TokenInfo]
        Complete token stream for the maintained source.
    doc_positions : set[tuple[int, int]]
        Token positions owned by docstrings and excluded from this check.
    path : str
        Repository-relative path used in diagnostics.
    findings : list[Finding]
        Mutable deterministic finding collection.

    Returns
    -------
    counts : Counter[str]
        Quote-form and explicit-exception facts.
    """

    counts: Counter[str] = Counter()
    fstring_depth = 0
    fstring_start = getattr(tokenize, "FSTRING_START", -1)
    fstring_end = getattr(tokenize, "FSTRING_END", -1)

    for token in tokens:
        if token.type == fstring_start:
            fstring_depth += 1

            continue

        if token.type == fstring_end:
            fstring_depth = max(0, fstring_depth - 1)
            continue

        if token.type != tokenize.STRING or token.start in doc_positions:
            continue

        parsed = _string_content(token)
        if parsed is None:
            continue

        _, quote, content = parsed
        multiline = token.start[0] != token.end[0]
        counts["total"] += 1
        counts["multiline" if multiline else "one_line"] += 1

        if quote == "'''":
            if multiline and '"""' in content:
                counts["triple_single_literal_exceptions"] += 1
            else:
                _finding(
                    findings,
                    path,
                    token,
                    RULE_STRINGS,
                    "ordinary multiline string should use triple double quotes"
                    "",
                )
        elif quote == '"""':
            if not multiline:
                _finding(
                    findings,
                    path,
                    token,
                    RULE_STRINGS,
                    "ordinary one-line string should use double quotes",
                )
        elif quote == "'":
            if fstring_depth:
                counts["fstring_expression_quote_exceptions"] += 1
            elif '"' in content:
                counts["single_quote_literal_exceptions"] += 1
            else:
                _finding(
                    findings,
                    path,
                    token,
                    RULE_STRINGS,
                    "ordinary one-line string should use double quotes",
                )

    return counts


def _position_within(
    token: tokenize.TokenInfo,
    node: ast.AST,
) -> bool:
    """
    Return whether one token begins within an AST node's source span.
    """

    end_line = getattr(node, "end_lineno", None)
    end_column = getattr(node, "end_col_offset", None)
    if end_line is None or end_column is None:
        return False

    start = (node.lineno, node.col_offset)
    end = (end_line, end_column)

    return start <= token.start < end


def _inside_parse_args(
    node: ast.AST,
    parents: dict[ast.AST, ast.AST],
) -> bool:
    """
    Return whether one node belongs to a function named parse_args.
    """

    current: ast.AST | None = node

    while current is not None:
        if isinstance(current, CALLABLE):
            return current.name == "parse_args"

        current = parents.get(current)

    return False


def _decoded_string(token: tokenize.TokenInfo) -> str | None:
    """
    Return one literal string token's value without executing code.
    """

    try:
        value = ast.literal_eval(token.string)
    except (SyntaxError, ValueError):
        return None

    return value if isinstance(value, str) else None


def _static_string_piece(token: tokenize.TokenInfo) -> str | None:
    """
    Return one string token's value when it has no dynamic interpolation.
    """

    try:
        expression = ast.parse(token.string, mode="eval").body
    except SyntaxError:
        return None

    if isinstance(expression, ast.Constant) and isinstance(
        expression.value,
        str,
    ):
        return expression.value

    if isinstance(expression, ast.JoinedStr) and all(
        isinstance(value, ast.Constant) and isinstance(value.value, str)
        for value in expression.values
    ):
        return "".join(value.value for value in expression.values)

    return None


def _next_help_piece(value: str) -> str:
    """
    Return the first movable prose word and its following separator.
    """

    match = re.match(r"^[ \t]*(\S+)([ \t]+(?=\S)|)", value)
    if match is None:
        return ""

    return match.group(0)


def _is_help_prose(paragraph: str) -> bool:
    """
    Return whether one help paragraph is prose rather than structured text.
    """

    lines = [line.strip() for line in paragraph.splitlines() if line.strip()]
    if not lines:
        return False

    structured = re.compile(
        r"^(?:[-+*]\s|"
        r"[A-Za-z_][A-Za-z0-9_]*\s*(?::=|=)|[\[(].*(?::=|=))"
        r"",
    )

    return not any(structured.match(line) for line in lines)


def _check_adjacent_prose_literals(
    tree: ast.Module,
    parents: dict[ast.AST, ast.AST],
    tokens: list[tokenize.TokenInfo],
    doc_positions: set[tuple[int, int]],
    lines: list[str],
    path: str,
    findings: list[Finding],
) -> Counter[str]:
    """
    Reject safely provable premature breaks in adjacent prose literals.

    Parameters
    ----------
    tree : ast.Module
        Parsed module used to locate constructed prose expressions.
    parents : dict[ast.AST, ast.AST]
        Parent map used to identify CLI help and expression ownership.
    tokens : list[tokenize.TokenInfo]
        Complete token stream retaining adjacent literal boundaries.
    doc_positions : set[tuple[int, int]]
        Token positions owned by docstrings and excluded from this check.
    lines : list[str]
        Physical source lines used for exact diagnostic columns.
    path : str
        Repository-relative path used in diagnostics.
    findings : list[Finding]
        Mutable deterministic finding collection.

    Returns
    -------
    counts : Counter[str]
        Inspected prose groups, exclusions, and movable-break findings.
    """

    counts: Counter[str] = Counter()

    for node in ast.walk(tree):
        if not isinstance(node, (ast.Constant, ast.JoinedStr)):
            continue

        if isinstance(node, ast.Constant) and not isinstance(node.value, str):
            continue

        if isinstance(parents.get(node), ast.JoinedStr):
            continue

        if node.lineno == (node.end_lineno or node.lineno):
            continue

        if (
            isinstance(node, ast.Constant)
            and any(separator in node.value for separator in ("\n", "\r"))
        ) or (
            isinstance(node, ast.JoinedStr)
            and any(
                isinstance(value, ast.Constant)
                and isinstance(value.value, str)
                and any(separator in value.value for separator in ("\n", "\r"))
                for value in node.values
            )
        ):
            counts["explicit_newline_exclusions"] += 1

            continue

        parent = parents.get(node)

        if (
            isinstance(parent, ast.keyword)
            and parent.arg == "help"
            and _inside_parse_args(parent, parents)
        ):
            counts["cli_help_exclusions"] += 1

            continue

        literal_tokens = [
            token
            for token in tokens
            if token.type == tokenize.STRING
            and token.start not in doc_positions
            and _position_within(token, node)
        ]
        if len(literal_tokens) < 2:
            continue

        counts["adjacent_groups"] += 1

        for current, following in pairwise(literal_tokens):
            if current.start[0] == following.start[0]:
                continue

            current_value = _static_string_piece(current)
            following_value = _static_string_piece(following)
            if (
                current_value is None
                or following_value is None
                or current_value.endswith(("\n ", "\r "))
            ):
                continue

            piece = _next_help_piece(following_value)
            stripped_piece = piece.strip()
            word_piece = (
                bool(stripped_piece)
                and stripped_piece[:1].isalpha()
                and current_value.endswith((" ", "\t"))
            )
            punctuation_piece = bool(stripped_piece) and all(
                not character.isalnum() and not character.isspace()
                for character in stripped_piece
            )
            if not piece or not (word_piece or punctuation_piece):
                continue

            source_line = lines[current.start[0] - 1]

            if len(source_line) + len(piece) <= 79:
                _finding(
                    findings,
                    path,
                    current,
                    RULE_MULTILINE,
                    "adjacent prose literal breaks before a word or "
                    "punctuation fragment that would still fit within 79 "
                    "columns",
                )
                counts["premature_prose_breaks"] += 1

    return counts


def _check_cli_help_layout(
    tree: ast.Module,
    parents: dict[ast.AST, ast.AST],
    tokens: list[tokenize.TokenInfo],
    lines: list[str],
    path: str,
    findings: list[Finding],
) -> Counter[str]:
    """
    Check greedy 79-column wrapping for parse_args help literals.

    Parameters
    ----------
    tree : ast.Module
        Parsed Python module.
    parents : dict[ast.AST, ast.AST]
        Parent relationships for parsed nodes.
    tokens : list[tokenize.TokenInfo]
        Token stream retaining literal source boundaries.
    lines : list[str]
        Physical source lines.
    path : str
        Repository-relative diagnostic path.
    findings : list[Finding]
        Mutable deterministic finding collection.

    Returns
    -------
    counts : Counter[str]
        Help-literal inspection and exclusion counts.
    """

    counts: Counter[str] = Counter()

    for node in ast.walk(tree):
        if (
            not isinstance(node, ast.keyword)
            or node.arg != "help"
            or not _inside_parse_args(node, parents)
        ):
            continue

        if not (
            isinstance(node.value, ast.Constant)
            and isinstance(node.value.value, str)
        ):
            counts["dynamic_exclusions"] += 1

            continue

        for paragraph in re.split(r"\n[ \t]*\n", node.value.value):
            prose = paragraph.strip()

            if not prose or not _is_help_prose(prose):
                counts["structured_help_exclusions"] += 1

                continue

            after_literal = re.sub(r"^'[^']+'[ \t]+", "", prose)
            first_letter = re.search(r"[A-Za-z]", after_literal)
            first_word = after_literal.split(maxsplit=1)[0]

            if (
                first_letter is not None
                and not first_letter.group(0).isupper()
                and not any(char.isupper() for char in first_word[1:])
            ):
                findings.append(
                    Finding(
                        path,
                        node.value.lineno,
                        node.value.col_offset + 1,
                        RULE_HELP_SENTENCES,
                        "help prose must begin with sentence capitalization",
                    ),
                )

            if prose[-1] not in ".?!":
                findings.append(
                    Finding(
                        path,
                        node.value.lineno,
                        node.value.col_offset + 1,
                        RULE_HELP_SENTENCES,
                        "help prose must end with terminal punctuation",
                    ),
                )

        literal_tokens = [
            token
            for token in tokens
            if token.type == tokenize.STRING
            and _position_within(token, node.value)
        ]

        if len(literal_tokens) == 1:
            counts["single_static_literals"] += 1
            token = literal_tokens[0]

            if token.start[0] != token.end[0]:
                _finding(
                    findings,
                    path,
                    token,
                    RULE_CLI_HELP_LAYOUT,
                    "parse_args help prose must use adjacent one-line literals"
                    "",
                )
            elif len(lines[token.start[0] - 1]) > 79:
                _finding(
                    findings,
                    path,
                    token,
                    RULE_CLI_HELP_LAYOUT,
                    "splittable static parse_args help source line exceeds "
                    "79 columns",
                )

            continue

        if not literal_tokens:
            counts["unrecognized_static_literals"] += 1

            continue

        counts["multiline_literal_groups"] += 1

        for token in literal_tokens:
            if token.start[0] != token.end[0]:
                _finding(
                    findings,
                    path,
                    token,
                    RULE_CLI_HELP_LAYOUT,
                    "parse_args help prose must use adjacent one-line literals"
                    "",
                )

            source_line = lines[token.start[0] - 1]

            if len(source_line) > 79:
                _finding(
                    findings,
                    path,
                    token,
                    RULE_CLI_HELP_LAYOUT,
                    "parse_args help source line exceeds 79 columns",
                )

        for current, following in pairwise(literal_tokens):
            if current.start[0] == current.end[0]:
                current_value = _decoded_string(current)
                following_value = _decoded_string(following)
            else:
                current_value = None
                following_value = None

            if (
                current_value is None
                or following_value is None
                or current_value.endswith(("\n", "\r"))
            ):
                continue

            piece = _next_help_piece(following_value)
            if not piece:
                continue

            source_line = lines[current.start[0] - 1]

            if len(source_line) + len(piece) <= 79:
                _finding(
                    findings,
                    path,
                    current,
                    RULE_CLI_HELP_LAYOUT,
                    "parse_args help line breaks before a prose word that "
                    "would still fit within 79 columns",
                )

    return counts


def _is_header_comment(
    line: int,
    value: str,
    lines: list[str],
) -> bool:
    """
    Return whether one comment belongs to the governed source header.
    """

    preamble = all(
        not source_line.strip() or source_line.lstrip().startswith("#")
        for source_line in lines[:line]
    )

    return preamble


def _comment_prose(value: str) -> str | None:
    """
    Return normalized ordinary prose, or None for a non-prose fragment.
    """

    content = value.removeprefix("# ").strip()

    if (
        not content
        or COMMENT_LITERAL_LABEL.match(content)
        or COMMENT_CODE_FRAGMENT.match(content)
        or not re.search(r"[A-Za-z]", content)
        or not re.search(r"\s", content)
    ):
        return None

    return content


def _comment_paragraphs(
    tokens: list[tokenize.TokenInfo],
    lines: list[str],
) -> Iterator[tuple[tokenize.TokenInfo, tokenize.TokenInfo, str]]:
    """
    Yield complete ordinary full-line and inline comment paragraphs.

    Parameters
    ----------
    tokens : list[tokenize.TokenInfo]
        Complete token stream containing comment markers.
    lines : list[str]
        Physical source lines used to recover inline comment context.

    Yields
    ------
    first, last, paragraph : tuple[
        tokenize.TokenInfo, tokenize.TokenInfo, str
    ]
        Boundary tokens and normalized prose for each complete paragraph.
    """

    paragraph: list[tokenize.TokenInfo] = []

    def flush() -> tuple[tokenize.TokenInfo, tokenize.TokenInfo, str] | None:
        if not paragraph:
            return None

        content = " ".join(
            token.string.removeprefix("# ").strip() for token in paragraph
        )
        prose = _comment_prose(f"# {content}")

        if prose is None:
            return None

        return paragraph[0], paragraph[-1], prose

    for token in tokens:
        if token.type != tokenize.COMMENT:
            continue

        line = token.start[0]
        value = token.string

        if _is_header_comment(line, value, lines) or DIRECTIVE.match(value):
            continue

        before = lines[line - 1][: token.start[1]]
        inline = bool(before.strip())

        if inline:
            completed = flush()

            if completed is not None:
                yield completed
                paragraph = []

            prose = _comment_prose(value)
            if prose is not None:
                yield token, token, prose

            continue

        continues = bool(
            paragraph
            and line == paragraph[-1].start[0] + 1
            and token.start[1] == paragraph[-1].start[1]
            and value != "#",
        )

        if value == "#":
            completed = flush()
            if completed is not None:
                yield completed

            paragraph = []
            continue

        if not continues:
            completed = flush()
            if completed is not None:
                yield completed

            paragraph = []

        paragraph.append(token)

    completed = flush()
    if completed is not None:
        yield completed


def _check_comment_prose(
    tokens: list[tokenize.TokenInfo],
    lines: list[str],
    path: str,
    findings: list[Finding],
) -> Counter[str]:
    """
    Check sentence capitalization and punctuation by whole comment paragraph.
    """

    counts: Counter[str] = Counter()

    for first, last, prose in _comment_paragraphs(tokens, lines):
        counts["prose_paragraphs"] += 1
        sentence = re.sub(r"^(?:TODO|NOTE|MAYBE):\s*", "", prose)
        first_letter = re.search(r"[A-Za-z]", sentence)

        if first_letter is not None and first_letter.group().islower():
            _finding(
                findings,
                path,
                first,
                RULE_COMMENTS,
                "ordinary comment prose must begin with a capital letter",
            )

        if re.search(r"[.!?][\"')\]]? {2,}(?=[A-Z])", sentence):
            _finding(
                findings,
                path,
                first,
                RULE_COMMENTS,
                "ordinary comment prose must use one space between sentences",
            )
            counts["repeated_sentence_spacing"] += 1

        if not prose.endswith((".", "?", "!")):
            _finding(
                findings,
                path,
                last,
                RULE_COMMENTS,
                "ordinary comment prose must end with terminal punctuation",
            )

    return counts


def _check_comments(
    tokens: list[tokenize.TokenInfo],
    lines: list[str],
    path: str,
    findings: list[Finding],
) -> Counter[str]:
    """
    Check ordinary comment markers, inline spacing, and source width.

    Parameters
    ----------
    tokens : list[tokenize.TokenInfo]
        Complete token stream containing comments.
    lines : list[str]
        Physical source lines used for spacing and width checks.
    path : str
        Repository-relative path used in diagnostics.
    findings : list[Finding]
        Mutable deterministic finding collection.

    Returns
    -------
    counts : Counter[str]
        Comment-form, prose, separator, and width facts.
    """

    counts: Counter[str] = Counter()
    ordinary_by_line = {
        token.start[0]: token
        for token in tokens
        if token.type == tokenize.COMMENT
        and not _is_header_comment(
            token.start[0],
            token.string,
            lines,
        )
        and not DIRECTIVE.match(token.string)
    }

    for token in tokens:
        if token.type != tokenize.COMMENT:
            continue

        value = token.string
        line = token.start[0]

        if _is_header_comment(line, value, lines):
            counts["headers"] += 1

            continue

        if DIRECTIVE.match(value):
            counts["directives"] += 1

            continue

        counts["ordinary"] += 1
        before = lines[line - 1][: token.start[1]]
        inline = bool(before.strip())

        if value == "#":
            prior = ordinary_by_line.get(line - 1)
            following = ordinary_by_line.get(line + 1)
            valid_separator = (
                not inline
                and prior is not None
                and following is not None
                and prior.string.startswith("# ")
                and following.string.startswith("# ")
                and not lines[prior.start[0] - 1][: prior.start[1]].strip()
                and not lines[following.start[0] - 1][
                    : following.start[1]
                ].strip()
            )

            if not valid_separator:
                _finding(
                    findings,
                    path,
                    token,
                    RULE_COMMENTS,
                    "bare '#' is allowed only between full-line comment "
                    "paragraphs",
                )
        elif not value.startswith("# ") or value.startswith("#  "):
            _finding(
                findings,
                path,
                token,
                RULE_COMMENTS,
                "ordinary comment marker must be '# ' or the separator '#'",
            )

        if inline and len(before) - len(before.rstrip(" ")) != 2:
            _finding(
                findings,
                path,
                token,
                RULE_COMMENTS,
                "inline comment requires exactly two preceding spaces",
            )

        if len(lines[line - 1]) > 79:
            _finding(
                findings,
                path,
                token,
                RULE_COMMENTS,
                "ordinary comment source line exceeds 79 columns",
            )

        if lines[line - 1].endswith((" ", "\t")):
            _finding(
                findings,
                path,
                token,
                RULE_COMMENTS,
                "ordinary comment source line has trailing whitespace",
            )

    counts.update(
        _check_comment_prose(
            tokens,
            lines,
            path,
            findings,
        ),
    )

    return counts


def _canonical_function_name(name: str) -> bool:
    """
    Return whether one function-like name has canonical casing.
    """

    if name in {"setUp", "tearDown"}:
        return True

    if name.startswith("__") and name.endswith("__"):
        return True

    return bool(SNAKE_CASE.fullmatch(name))


def _opaque_identifier_segments(name: str) -> list[str]:
    """
    Return prohibited private-shorthand segments from one identifier.
    """

    return sorted(
        {
            segment.casefold()
            for segment in name.strip("_").split("_")
            if segment.casefold() in OPAQUE_IDENTIFIER_SEGMENTS
        },
    )


def _check_opaque_identifier(
    name: str,
    path: str,
    line: int,
    column: int,
    findings: list[Finding],
    counts: Counter[str],
) -> None:
    """
    Reject implementation shorthand with an established descriptive form.
    """

    opaque_segments = _opaque_identifier_segments(name)
    if not opaque_segments:
        return

    findings.append(
        Finding(
            path,
            line,
            column,
            RULE_NAMING,
            f"identifier '{name}' contains opaque shorthand segment(s): "
            + ", ".join(opaque_segments),
        ),
    )
    counts["opaque_shorthand"] += 1


def _check_naming(
    tree: ast.Module,
    path: str,
    findings: list[Finding],
) -> Counter[str]:
    """
    Check bounded role-aware casing and visibility forms.

    Parameters
    ----------
    tree : ast.Module
        Parsed module whose bound names are inspected.
    path : str
        Repository-relative path used in diagnostics.
    findings : list[Finding]
        Mutable deterministic finding collection.

    Returns
    -------
    counts : Counter[str]
        Name-role counts and casing findings.
    """

    counts: Counter[str] = Counter()
    module_name = Path(path).stem

    if (
        path.endswith(".py")
        and module_name != "__init__"
        and not SNAKE_CASE.fullmatch(module_name)
    ):
        findings.append(
            Finding(
                path,
                1,
                1,
                RULE_NAMING,
                f"module name '{module_name}' is not lowercase snake_case",
            ),
        )

    for node in ast.walk(tree):
        if isinstance(node, CALLABLE):
            counts["functions"] += 1

            if not _canonical_function_name(node.name):
                findings.append(
                    Finding(
                        path,
                        node.lineno,
                        node.col_offset + 1,
                        RULE_NAMING,
                        f"function or method '{node.name}' is not snake_case",
                    ),
                )

            _check_opaque_identifier(
                node.name,
                path,
                node.lineno,
                node.col_offset + 1,
                findings,
                counts,
            )

            for argument in (
                *node.args.posonlyargs,
                *node.args.args,
                *node.args.kwonlyargs,
                node.args.vararg,
                node.args.kwarg,
            ):
                if argument is None:
                    continue

                counts["parameters"] += 1

                if argument.arg == "_":
                    counts["unused_parameter_exclusions"] += 1

                    continue

                if not _canonical_function_name(argument.arg):
                    findings.append(
                        Finding(
                            path,
                            argument.lineno,
                            argument.col_offset + 1,
                            RULE_NAMING,
                            f"parameter '{argument.arg}' is not snake_case",
                        ),
                    )

                _check_opaque_identifier(
                    argument.arg,
                    path,
                    argument.lineno,
                    argument.col_offset + 1,
                    findings,
                    counts,
                )
        elif isinstance(node, ast.ClassDef):
            counts["classes"] += 1
            type_name = node.name.lstrip("_")

            if not type_name or not UPPER_CAMEL.fullmatch(type_name):
                findings.append(
                    Finding(
                        path,
                        node.lineno,
                        node.col_offset + 1,
                        RULE_NAMING,
                        f"class '{node.name}' is not UpperCamelCase",
                    ),
                )

            _check_opaque_identifier(
                node.name,
                path,
                node.lineno,
                node.col_offset + 1,
                findings,
                counts,
            )
        elif isinstance(node, ast.Name) and isinstance(node.ctx, ast.Store):
            _check_opaque_identifier(
                node.id,
                path,
                node.lineno,
                node.col_offset + 1,
                findings,
                counts,
            )
        elif isinstance(node, ast.Attribute) and isinstance(
            node.ctx,
            ast.Store,
        ):
            _check_opaque_identifier(
                node.attr,
                path,
                node.lineno,
                max(1, (node.end_col_offset or 1) - len(node.attr) + 1),
                findings,
                counts,
            )

    return counts


def _blank_count_before(lines: list[str], line: int) -> int:
    """
    Return consecutive blank physical lines before one 1-based line.
    """

    count = 0
    index = line - 2

    while index >= 0 and not lines[index].strip():
        count += 1
        index -= 1

    return count


def _start_line(node: ast.AST) -> int:
    """
    Return the first decorator or definition line.
    """

    decorators = getattr(node, "decorator_list", [])
    return min(
        [getattr(node, "lineno", 1)]
        + [decorator.lineno for decorator in decorators],
    )


def _attached_comment_start(lines: list[str], line: int) -> int:
    """
    Return the start of a contiguous comment block attached to a definition.
    """

    index = line - 2

    while index >= 0 and lines[index].lstrip().startswith("#"):
        line = index + 1
        index -= 1

    return line


def _statement_suites(
    tree: ast.Module,
) -> Iterator[tuple[ast.AST, str, list[ast.stmt]]]:
    """
    Yield every nonempty statement suite with its owner and field.
    """

    for owner in ast.walk(tree):
        for field, value in ast.iter_fields(owner):
            if (
                isinstance(value, list)
                and value
                and all(isinstance(item, ast.stmt) for item in value)
            ):
                yield owner, field, value


def _is_compact_guard(node: ast.stmt) -> bool:
    """
    Return whether one if statement owns one direct guarded transfer.
    """

    return (
        isinstance(node, ast.If)
        and not node.orelse
        and len(node.body) == 1
        and not _is_control_flow(node.body[0])
        and _transfer_kind(node.body[0]) is not None
    )


def _needs_control_boundary(node: ast.stmt) -> bool:
    """
    Return whether one compound statement needs a visible sibling boundary.
    """

    return isinstance(
        node,
        (
            ast.AsyncFor,
            ast.AsyncWith,
            ast.For,
            ast.If,
            ast.Match,
            ast.Try,
            ast.While,
            ast.With,
        ),
    ) and not _is_compact_guard(node)


def _is_control_flow(node: ast.stmt) -> bool:
    """
    Return whether one statement owns a branch or iterative suite.
    """

    return isinstance(
        node,
        (
            ast.AsyncFor,
            ast.AsyncWith,
            ast.For,
            ast.If,
            ast.Match,
            ast.Try,
            ast.While,
            ast.With,
        ),
    )


def _qualified_call_name(node: ast.AST) -> str:
    """
    Return one statically visible dotted call name.
    """

    if isinstance(node, ast.Name):
        return node.id

    if isinstance(node, ast.Attribute):
        owner = _qualified_call_name(node.value)
        return f"{owner}.{node.attr}" if owner else node.attr

    return ""


def _call_leaf(node: ast.AST) -> str:
    """
    Return the final component of one statically qualified call name.
    """

    return _qualified_call_name(node).rsplit(".", maxsplit=1)[-1]


def _transfer_kind(node: ast.stmt) -> str | None:
    """
    Return the recognized transfer kind for one statement, when applicable.
    """

    if isinstance(node, ast.Return):
        return "return"

    if isinstance(node, ast.Raise):
        return "raise"

    if isinstance(node, ast.Break):
        return "break"

    if isinstance(node, ast.Continue):
        return "continue"

    if isinstance(node, ast.Expr) and isinstance(
        node.value,
        (ast.Yield, ast.YieldFrom),
    ):
        if isinstance(node.value, ast.YieldFrom):
            return "yield_from"

        return "yield"

    if not isinstance(node, ast.Expr) or not isinstance(node.value, ast.Call):
        return None

    if _qualified_call_name(node.value.func) in {
        "exit",
        "os._exit",
        "pytest.fail",
        "quit",
        "sys.exit",
    }:
        return "process_exit"

    return None


def _is_direct_side_effect(node: ast.stmt) -> bool:
    """
    Return whether one direct statement visibly mutates or emits output.
    """

    if isinstance(node, (ast.AugAssign, ast.Delete)):
        return True

    if not isinstance(node, ast.Expr) or not isinstance(node.value, ast.Call):
        return False

    name = _qualified_call_name(node.value.func)
    leaf = name.rsplit(".", maxsplit=1)[-1]

    return name == "print" or leaf in {
        "add",
        "append",
        "clear",
        "close",
        "discard",
        "extend",
        "flush",
        "insert",
        "mkdir",
        "remove",
        "rename",
        "replace",
        "reverse",
        "rmdir",
        "sort",
        "touch",
        "unlink",
        "update",
        "write",
        "writelines",
    }


def _attached_preparation(
    suite: list[ast.stmt],
    index: int,
    lines: list[str],
) -> tuple[ast.stmt, ...]:
    """
    Return statements physically attached before one suite member.
    """

    start = index

    while start > 0:
        current_start = _attached_comment_start(
            lines,
            suite[start].lineno,
        )
        if _blank_count_before(lines, current_start):
            break

        start -= 1

    return tuple(suite[start:index])


def _requires_transfer_boundary(
    suite: list[ast.stmt],
    index: int,
    lines: list[str],
) -> bool:
    """
    Return whether an attached transfer follows a substantive source phase.
    """

    preparation = _attached_preparation(suite, index, lines)
    if not preparation:
        return False

    if len(preparation) > 1:
        return True

    predecessor = preparation[0]
    physical_span = (
        (predecessor.end_lineno or predecessor.lineno) - predecessor.lineno + 1
    )

    return (
        _needs_control_boundary(predecessor)
        or _is_direct_side_effect(predecessor)
        or physical_span > 1
    )


def _is_test_assertion(node: ast.stmt) -> bool:
    """
    Return whether one direct test statement performs an assertion.
    """

    if isinstance(node, ast.Assert):
        return True

    if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call):
        function = node.value.func
        name = (
            function.id
            if isinstance(function, ast.Name)
            else function.attr
            if isinstance(function, ast.Attribute)
            else ""
        )

        return name.startswith("assert") or name == "fail"

    if not isinstance(node, (ast.With, ast.AsyncWith)):
        return False

    return any(
        isinstance(item.context_expr, ast.Call)
        and (
            (
                isinstance(item.context_expr.func, ast.Name)
                and item.context_expr.func.id in {"raises", "warns"}
            )
            or (
                isinstance(item.context_expr.func, ast.Attribute)
                and item.context_expr.func.attr in {"raises", "warns"}
            )
        )
        for item in node.items
    )


def _semantic_phase_kind(node: ast.stmt) -> str:
    """
    Return one safely recognized statement phase.
    """

    if _is_test_assertion(node):
        return "assertion"

    if isinstance(node, (ast.Import, ast.ImportFrom)):
        return "local_import"

    if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call):
        name = _qualified_call_name(node.value.func)
        leaf = name.rsplit(".", maxsplit=1)[-1]
        if leaf.startswith(("validate", "verify")):
            return "verification_action"

    if _needs_control_boundary(node):
        return "control_flow"

    return "setup"


def _assigned_names(node: ast.stmt) -> set[str]:
    """
    Return direct names bound by one assignment statement.
    """

    def bound_names(target: ast.AST) -> set[str]:
        if isinstance(target, ast.Name):
            return {target.id}

        if isinstance(target, ast.Starred):
            return bound_names(target.value)

        if isinstance(target, (ast.List, ast.Tuple)):
            return {name for item in target.elts for name in bound_names(item)}

        return set()

    if isinstance(node, ast.Assign):
        targets = node.targets
    elif isinstance(node, ast.AnnAssign):
        targets = [node.target]
    else:
        return set()

    return {name for target in targets for name in bound_names(target)}


ACTION_RESULT_WORDS = {
    "analysis",
    "analyses",
    "binding",
    "bindings",
    "capture",
    "captures",
    "coverage",
    "decision",
    "decisions",
    "diagnostic",
    "diagnostics",
    "digest",
    "evidence",
    "fact",
    "facts",
    "finding",
    "findings",
    "fingerprint",
    "fingerprints",
    "fp",
    "index",
    "indices",
    "inventory",
    "ledger",
    "markdown",
    "message",
    "messages",
    "observed",
    "output",
    "outputs",
    "payload",
    "progress",
    "record",
    "records",
    "report",
    "reports",
    "result",
    "results",
    "selection",
    "selections",
    "status",
    "summary",
    "unit",
    "units",
    "warning",
    "warnings",
}

CONSTRUCTION_CALLS = {
    "Counter",
    "Namespace",
    "Path",
    "PurePath",
    "PurePosixPath",
    "bool",
    "bytearray",
    "bytes",
    "defaultdict",
    "dict",
    "float",
    "int",
    "list",
    "set",
    "str",
    "tuple",
}


def _assignment_value(node: ast.stmt) -> ast.AST | None:
    """
    Return the value owned by one direct assignment statement.
    """

    if isinstance(node, ast.Assign):
        return node.value

    if isinstance(node, ast.AnnAssign):
        return node.value

    return None


def _contains_call(node: ast.AST | None) -> bool:
    """
    Return whether one expression evaluates a call.
    """

    return node is not None and any(
        isinstance(child, ast.Call) for child in ast.walk(node)
    )


def _operation_calls(node: ast.AST | None) -> set[str]:
    """
    Return non-construction call names evaluated by one expression.
    """

    if node is None:
        return set()

    return {
        _call_leaf(child.func)
        for child in ast.walk(node)
        if isinstance(child, ast.Call)
        and _call_leaf(child.func) not in CONSTRUCTION_CALLS
    }


def _loaded_names(node: ast.AST) -> set[str]:
    """
    Return names read by one statement.
    """

    return {
        child.id
        for child in ast.walk(node)
        if isinstance(child, ast.Name) and isinstance(child.ctx, ast.Load)
    }


def _is_action_result_assignment(node: ast.stmt) -> bool:
    """
    Return whether an assignment visibly captures an operation result.
    """

    value = _assignment_value(node)
    if not _contains_call(value) or not _operation_calls(value):
        return False

    return any(
        word in ACTION_RESULT_WORDS
        for name in _assigned_names(node)
        for word in name.casefold().split("_")
    )


def _owning_callable_name(
    owner: ast.AST,
    parents: dict[ast.AST, ast.AST],
) -> str:
    """
    Return the nearest callable name enclosing one statement suite.
    """

    current: ast.AST | None = owner

    while current is not None:
        if isinstance(current, CALLABLE):
            return current.name

        current = parents.get(current)

    return ""


def _data_loop_transition(
    previous: ast.stmt,
    current: ast.stmt,
) -> bool:
    """
    Return whether a multiline case inventory feeds an attached loop.
    """

    if not isinstance(current, (ast.For, ast.AsyncFor)):
        return False

    physical_span = (
        (previous.end_lineno or previous.lineno) - previous.lineno + 1
    )
    if physical_span <= 6:
        return False

    assigned = _assigned_names(previous)
    consumed = {
        child.id
        for child in ast.walk(current.iter)
        if isinstance(child, ast.Name)
    }

    return bool(assigned & consumed)


def _semantic_phase_transition(
    previous: ast.stmt,
    current: ast.stmt,
) -> bool:
    """
    Return whether two recognized statements need a visible boundary.
    """

    previous_kind = _semantic_phase_kind(previous)
    current_kind = _semantic_phase_kind(current)
    if previous_kind == current_kind:
        return False

    explicit_phases = {"assertion", "local_import"}
    if {previous_kind, current_kind} & explicit_phases:
        return True

    if (
        previous_kind == "verification_action"
        and current_kind != "verification_action"
    ):
        return True

    return _data_loop_transition(previous, current)


def _action_result_transition(
    suite: list[ast.stmt],
    index: int,
    lines: list[str],
    owner: ast.AST,
    parents: dict[ast.AST, ast.AST],
) -> bool:
    """
    Return whether an operation-result assignment changes the visible phase.
    """

    previous = suite[index - 1]
    current = suite[index]
    previous_action = _is_action_result_assignment(previous)
    current_action = _is_action_result_assignment(current)
    if previous_action == current_action:
        return False

    previous_names = _assigned_names(previous)
    if (
        previous_action
        and _is_direct_side_effect(current)
        and previous_names & _loaded_names(current)
    ):
        return False

    if previous_action and _assigned_names(current):
        return False

    callable_name = _owning_callable_name(owner, parents)
    if callable_name.startswith("test"):
        return True

    if previous_action:
        return True

    preparation = _attached_preparation(suite, index, lines)
    previous_end = previous.end_lineno or previous.lineno
    previous_span = previous_end - previous.lineno + 1

    return (
        len(preparation) > 1
        or _is_direct_side_effect(previous)
        or previous_span > 6
    )


def _inside_callable(
    owner: ast.AST,
    parents: dict[ast.AST, ast.AST],
) -> bool:
    """
    Return whether one suite belongs to a callable.
    """

    current: ast.AST | None = owner

    while current is not None:
        if isinstance(current, CALLABLE):
            return True

        current = parents.get(current)

    return False


def _check_topology(
    tree: ast.Module,
    parents: dict[ast.AST, ast.AST],
    lines: list[str],
    path: str,
    findings: list[Finding],
) -> Counter[str]:
    """
    Check bounded module and definition blank-line topology.

    Parameters
    ----------
    tree : ast.Module
        Parsed Python module.
    parents : dict[ast.AST, ast.AST]
        Parent lookup for callable and suite ownership.
    lines : list[str]
        Physical source lines.
    path : str
        Repository-relative diagnostic path.
    findings : list[Finding]
        Mutable deterministic finding collection.

    Returns
    -------
    counts : Counter[str]
        Topology and semantic-layout inspection counts.
    """

    counts: Counter[str] = Counter()
    seen_executable_statement = False

    for index, node in enumerate(tree.body):
        module_docstring = (
            index == 0
            and isinstance(node, ast.Expr)
            and isinstance(node.value, ast.Constant)
            and isinstance(node.value.value, str)
        )
        if module_docstring:
            continue

        if isinstance(node, (ast.Import, ast.ImportFrom)):
            counts["module_imports"] += 1

            if seen_executable_statement:
                findings.append(
                    Finding(
                        path,
                        node.lineno,
                        node.col_offset + 1,
                        RULE_TOPOLOGY,
                        "module import must precede executable module "
                        "statements",
                    ),
                )
                counts["late_module_imports"] += 1

            continue

        seen_executable_statement = True

    for node in tree.body:
        if not isinstance(node, DEFINITION):
            continue

        counts["top_level_definitions"] += 1
        line = _start_line(node)
        section_line = _attached_comment_start(lines, line)

        if _blank_count_before(lines, section_line) != 2:
            findings.append(
                Finding(
                    path,
                    line,
                    1,
                    RULE_TOPOLOGY,
                    "top-level definition requires exactly two blank lines "
                    "before it",
                ),
            )

        if not isinstance(node, ast.ClassDef):
            continue

        prior_method: ast.AST | None = None

        for member in node.body:
            if not isinstance(member, CALLABLE):
                continue

            counts["methods"] += 1

            if prior_method is not None:
                member_line = _start_line(member)

                if _blank_count_before(lines, member_line) != 1:
                    findings.append(
                        Finding(
                            path,
                            member_line,
                            member.col_offset + 1,
                            RULE_TOPOLOGY,
                            "methods require exactly one blank line between "
                            "definitions",
                        ),
                    )

            prior_method = member

    for owner, _field, suite in _statement_suites(tree):
        for index, current in enumerate(suite):
            if index == 0:
                continue

            previous = suite[index - 1]
            missing_control_boundary = _needs_control_boundary(
                previous,
            ) or _needs_control_boundary(current)
            completed_compact_guard = _is_compact_guard(previous)
            transfer_kind = _transfer_kind(current)
            missing_transfer_boundary = (
                transfer_kind is not None
                and _requires_transfer_boundary(suite, index, lines)
            )

            inside_callable = _inside_callable(owner, parents)
            previous_phase_kind = _semantic_phase_kind(previous)
            previous_end = previous.end_lineno or previous.lineno
            previous_span = previous_end - previous.lineno + 1
            validation_after_setup = (
                inside_callable
                and previous_phase_kind == "setup"
                and _semantic_phase_kind(current) == "verification_action"
                and (
                    len(_attached_preparation(suite, index, lines)) > 1
                    or previous_span > 6
                )
            )
            semantic_phase_transition = (
                inside_callable
                and _semantic_phase_transition(previous, current)
            ) or validation_after_setup

            action_result_transition = (
                inside_callable
                and _action_result_transition(
                    suite,
                    index,
                    lines,
                    owner,
                    parents,
                )
            )

            if not (
                missing_control_boundary
                or completed_compact_guard
                or missing_transfer_boundary
                or semantic_phase_transition
                or action_result_transition
            ):
                continue

            current_start = _attached_comment_start(lines, current.lineno)
            if _blank_count_before(lines, current_start):
                continue

            if completed_compact_guard:
                message = (
                    "completed compact transfer guard requires a blank-line "
                    "boundary from the following statement"
                )
                count_key = "missing_compact_guard_boundaries"
            elif semantic_phase_transition or action_result_transition:
                message = (
                    "semantic setup, action, validation, and assertion phases "
                    "require a blank-line boundary"
                )
                count_key = "missing_semantic_phase_boundaries"
            elif missing_transfer_boundary:
                message = (
                    "terminal transfer after a substantive phase requires a "
                    "blank-line boundary"
                )
                count_key = "missing_transfer_boundaries"
            else:
                message = (
                    "noncompact control-flow phase requires a blank-line "
                    "boundary from its sibling statement"
                )
                count_key = "missing_control_flow_boundaries"

            findings.append(
                Finding(
                    path,
                    current.lineno,
                    current.col_offset + 1,
                    RULE_TOPOLOGY,
                    message,
                ),
            )
            counts[count_key] += 1

    return counts


def _significant(
    tokens: list[tokenize.TokenInfo],
) -> list[tokenize.TokenInfo]:
    """
    Return tokens relevant to delimited source-form checks.
    """

    ignored = {
        tokenize.ENCODING,
        tokenize.ENDMARKER,
        tokenize.INDENT,
        tokenize.DEDENT,
        tokenize.NEWLINE,
        tokenize.NL,
        tokenize.COMMENT,
    }

    return [token for token in tokens if token.type not in ignored]


def _matched_delimiters(
    significant: list[tokenize.TokenInfo],
) -> tuple[dict[int, int], dict[int, int]]:
    """
    Return opening-to-closing and closing-to-opening token indexes.
    """

    stack: list[tuple[int, tokenize.TokenInfo]] = []
    opening_to_closing: dict[int, int] = {}
    closing_to_opening: dict[int, int] = {}
    pairs = {")": "(", "]": "[", "}": "{"}

    for index, token in enumerate(significant):
        if token.string in {"(", "[", "{"}:
            stack.append((index, token))

            continue

        if token.string not in pairs or not stack:
            continue

        open_index, opening = stack.pop()

        if opening.string != pairs[token.string]:
            continue

        opening_to_closing[open_index] = index
        closing_to_opening[index] = open_index

    return opening_to_closing, closing_to_opening


def _check_detached_comparisons(
    significant: list[tokenize.TokenInfo],
    path: str,
    findings: list[Finding],
) -> Counter[str]:
    """
    Reject comparison operators detached from their left operands.
    """

    counts: Counter[str] = Counter()
    symbolic = {"==", "!=", "<", "<=", ">", ">="}

    for index, token in enumerate(significant):
        if index == 0:
            continue

        is_comparison = token.string in symbolic | {"in", "is"}
        is_not_in = (
            token.string == "not"
            and index + 1 < len(significant)
            and significant[index + 1].string == "in"
        )
        if not (is_comparison or is_not_in):
            continue

        counts["comparisons"] += 1
        left = significant[index - 1]
        if left.end[0] == token.start[0]:
            continue

        _finding(
            findings,
            path,
            token,
            RULE_MULTILINE,
            "comparison operator must remain with its left operand",
        )

    return counts


def _check_block_literal_boundaries(
    tokens: list[tokenize.TokenInfo],
    doc_positions: set[tuple[int, int]],
    lines: list[str],
    path: str,
    findings: list[Finding],
) -> Counter[str]:
    """
    Check content boundaries for every multiline triple-quoted literal.

    Parameters
    ----------
    tokens : list[tokenize.TokenInfo]
        Complete token stream retaining triple-quote boundaries.
    doc_positions : set[tuple[int, int]]
        Token positions owned by docstrings and excluded from this check.
    lines : list[str]
        Physical source lines used to identify delimiter placement.
    path : str
        Repository-relative path used in diagnostics.
    findings : list[Finding]
        Mutable deterministic finding collection.

    Returns
    -------
    counts : Counter[str]
        Checked block literals and opening or closing boundary findings.
    """

    counts: Counter[str] = Counter()
    fstring_start = getattr(tokenize, "FSTRING_START", -1)
    fstring_end = getattr(tokenize, "FSTRING_END", -1)

    for token in tokens:
        if token.start in doc_positions:
            continue

        if token.type == fstring_start:
            if not token.string.endswith(('"""', "'''")):
                continue

            counts["recognized_block_literals"] += 1
            remainder = lines[token.end[0] - 1][token.end[1] :]
            if not remainder:
                continue

            _finding(
                findings,
                path,
                token,
                RULE_MULTILINE,
                "multiline string opener must be on its own line",
            )

            continue

        if token.type == fstring_end:
            if token.string not in {'"""', "'''"}:
                continue

            prefix = lines[token.start[0] - 1][: token.start[1]]
            if not prefix.strip():
                continue

            _finding(
                findings,
                path,
                token,
                RULE_MULTILINE,
                "multiline string closer must be on its own line",
            )

            continue

        if token.type != tokenize.STRING or token.start[0] == token.end[0]:
            continue

        parsed = _string_content(token)
        if parsed is None:
            continue

        _, quote, content = parsed
        if quote not in {'"""', "'''"}:
            continue

        counts["recognized_block_literals"] += 1

        if not content.startswith(("\n", "\r")):
            _finding(
                findings,
                path,
                token,
                RULE_MULTILINE,
                "multiline string opener must be on its own line",
            )

        if not content.endswith(("\n", "\r")):
            _finding(
                findings,
                path,
                token,
                RULE_MULTILINE,
                "multiline string closer must be on its own line",
            )

    return counts


def _delimited_items(node: ast.AST) -> list[tuple[ast.AST, int]]:
    """
    Return directly governed items and their source columns.
    """

    if isinstance(node, ast.Call):
        return [
            (item, item.col_offset) for item in [*node.args, *node.keywords]
        ]

    if isinstance(node, (ast.List, ast.Tuple, ast.Set)):
        return [(item, item.col_offset) for item in node.elts]

    if isinstance(node, ast.Dict):
        return [
            (
                (key, key.col_offset)
                if key is not None
                else (value, value.col_offset - 2)
            )
            for key, value in zip(
                node.keys,
                node.values,
                strict=True,
            )
        ]

    return []


def _compact_call_argument_row(
    node: ast.AST,
    items: list[tuple[ast.AST, int]],
    opening: tokenize.TokenInfo,
    closing: tokenize.TokenInfo,
    lines: list[str],
) -> int | None:
    """
    Return the source row for one bounded compact call payload.
    """

    if not isinstance(node, ast.Call) or len(items) < 2:
        return None

    row = items[0][0].lineno
    if not opening.start[0] < row < closing.start[0]:
        return None

    if any(item.lineno != row or item.end_lineno != row for item, _ in items):
        return None

    if len(lines[row - 1]) > 79:
        return None

    return row


def _check_parameter_lists(
    tree: ast.Module,
    significant: list[tokenize.TokenInfo],
    lines: list[str],
    path: str,
    findings: list[Finding],
    opening_to_closing: dict[int, int],
) -> Counter[str]:
    """
    Check multiline function and method parameter-list source form.

    Parameters
    ----------
    tree : ast.Module
        Parsed Python module.
    significant : list[tokenize.TokenInfo]
        Nontrivia tokens used to locate delimited structures.
    lines : list[str]
        Physical source lines.
    path : str
        Repository-relative diagnostic path.
    findings : list[Finding]
        Mutable deterministic finding collection.
    opening_to_closing : dict[int, int]
        Matching delimiter-token indices.

    Returns
    -------
    counts : Counter[str]
        Parameter-list inspection and exclusion counts.
    """

    counts: Counter[str] = Counter()

    for node in ast.walk(tree):
        if not isinstance(node, CALLABLE):
            continue

        name_index = next(
            (
                index
                for index, token in enumerate(significant)
                if token.type == tokenize.NAME
                and token.string == node.name
                and token.start[0] == node.lineno
                and token.start[1] >= node.col_offset
            ),
            None,
        )

        if name_index is None:
            continue

        open_index = next(
            (
                index
                for index in range(name_index + 1, len(significant))
                if significant[index].string == "("
                and index in opening_to_closing
            ),
            None,
        )

        if open_index is None:
            continue

        close_index = opening_to_closing[open_index]
        opening = significant[open_index]
        closing = significant[close_index]
        if opening.start[0] == closing.start[0]:
            continue

        counts["parameter_lists"] += 1
        inner = significant[open_index + 1 : close_index]
        if not inner:
            continue

        if inner[0].start[0] == opening.start[0]:
            _finding(
                findings,
                path,
                opening,
                RULE_MULTILINE,
                "multiline parameter list must break after the opening "
                "delimiter",
            )

        if inner[-1].end[0] == closing.start[0]:
            _finding(
                findings,
                path,
                closing,
                RULE_MULTILINE,
                "multiline parameter list must break before the closing "
                "delimiter",
            )

        segments: list[list[tokenize.TokenInfo]] = []
        segment: list[tokenize.TokenInfo] = []
        depth = 0

        for token in inner:
            if token.string in {"(", "[", "{"}:
                depth += 1
            elif token.string in {")", "]", "}"}:
                depth -= 1

            if token.string == "," and depth == 0:
                if segment:
                    segments.append(segment)
                    segment = []

                continue

            segment.append(token)

        if segment:
            segments.append(segment)
            _finding(
                findings,
                path,
                closing,
                RULE_MULTILINE,
                "multiline parameter list requires a trailing comma",
            )

        expected_indent = node.col_offset + 4
        prior_line: int | None = None

        for item in segments:
            first = item[0]

            if first.start[0] == prior_line:
                _finding(
                    findings,
                    path,
                    first,
                    RULE_MULTILINE,
                    "multiline parameter list requires one item per line",
                )

            prior_line = first.start[0]

            if first.start[1] != expected_indent:
                _finding(
                    findings,
                    path,
                    first,
                    RULE_MULTILINE,
                    "multiline parameter item requires one continuation indent"
                    "",
                )

        if closing.start[1] != node.col_offset:
            _finding(
                findings,
                path,
                closing,
                RULE_MULTILINE,
                "parameter-list closing delimiter must align with the "
                "definition line",
            )

    return counts


def _check_simple_return_annotations(
    tree: ast.Module,
    significant: list[tokenize.TokenInfo],
    path: str,
    findings: list[Finding],
    opening_to_closing: dict[int, int],
) -> Counter[str]:
    """
    Reject multiline parentheses around one simple return annotation.
    """

    counts: Counter[str] = Counter()

    for node in ast.walk(tree):
        if not isinstance(node, CALLABLE) or node.returns is None:
            continue

        arrow_index = next(
            (
                index
                for index, token in enumerate(significant)
                if token.string == "->"
                and node.lineno <= token.start[0] <= node.returns.lineno
                and token.start[1] >= node.col_offset
            ),
            None,
        )

        if arrow_index is None or arrow_index + 1 >= len(significant):
            continue

        open_index = arrow_index + 1
        opening = significant[open_index]
        if opening.string != "(" or open_index not in opening_to_closing:
            continue

        close_index = opening_to_closing[open_index]
        closing = significant[close_index]
        inner = significant[open_index + 1 : close_index]
        if (
            len(inner) != 1
            or inner[0].type != tokenize.NAME
            or opening.start[0] == closing.start[0]
        ):
            continue

        counts["multiline_simple_return_annotations"] += 1
        _finding(
            findings,
            path,
            opening,
            RULE_MULTILINE,
            "simple return annotation must remain on the definition line",
        )

    return counts


def _check_compound_condition_readability(
    tree: ast.Module,
    path: str,
    findings: list[Finding],
) -> Counter[str]:
    """
    Require named predicates when a compound condition expands a collection.
    """

    counts: Counter[str] = Counter()

    for node in ast.walk(tree):
        if not isinstance(node, (ast.If, ast.While)):
            continue

        if not isinstance(node.test, ast.BoolOp):
            continue

        collection_is_multiline = any(
            isinstance(part, (ast.List, ast.Tuple, ast.Set, ast.Dict))
            and part.end_lineno is not None
            and part.lineno != part.end_lineno
            for part in ast.walk(node.test)
        )

        if not collection_is_multiline:
            continue

        counts["compound_conditions_with_multiline_collections"] += 1
        findings.append(
            Finding(
                path,
                node.test.lineno,
                node.test.col_offset + 1,
                RULE_MULTILINE,
                "compound condition with a multiline collection requires "
                "named predicates",
            ),
        )

    return counts


def _multiline_token_pair(
    node: ast.AST,
    significant: list[tokenize.TokenInfo],
    matched: dict[int, int],
    checked_pairs: set[tuple[int, int]],
) -> tuple[int, int] | None:
    """
    Return one unreviewed multiline delimiter pair for an AST node.
    """

    closing_delimiters = {")", "]", "}"}
    close_index = next(
        (
            index
            for index, token in enumerate(significant)
            if token.string in closing_delimiters
            and token.end == (node.end_lineno, node.end_col_offset)
        ),
        None,
    )

    if close_index is None or close_index not in matched:
        return None

    open_index = matched[close_index]
    pair = (open_index, close_index)

    if pair in checked_pairs:
        return None

    checked_pairs.add(pair)
    opening = significant[open_index]
    closing = significant[close_index]

    return pair if opening.start[0] != closing.start[0] else None


def _normalized_multiline_item_position(
    item: ast.AST,
    item_offset: int,
    expected_offset: int,
    lines: list[str],
) -> tuple[int, int]:
    """
    Normalize one parenthesized item's physical line and indentation.
    """

    item_line = item.lineno
    same_line_prefix = lines[item_line - 1][
        expected_offset:item_offset
    ].strip()

    if same_line_prefix and set(same_line_prefix) == {"("}:
        return item_line, expected_offset

    if item_offset <= expected_offset:
        return item_line, item_offset

    prior_index = item_line - 2

    while prior_index >= 0 and not lines[prior_index].strip():
        prior_index -= 1

    prior_indent = (
        len(lines[prior_index]) - len(lines[prior_index].lstrip(" "))
        if prior_index >= 0
        else None
    )
    parenthesized_item = (
        prior_index >= 0
        and lines[prior_index].strip() == "("
        and prior_indent == expected_offset
    )

    if parenthesized_item:
        return prior_index + 1, expected_offset

    return item_line, item_offset


def _check_multiline_item_layout(
    items: list[tuple[ast.AST, int]],
    expected_offset: int,
    lines: list[str],
    path: str,
    findings: list[Finding],
) -> None:
    """
    Check expanded multiline item separation and continuation indentation.
    """

    prior_line: int | None = None

    for item, item_offset in items:
        item_line, item_offset = _normalized_multiline_item_position(
            item,
            item_offset,
            expected_offset,
            lines,
        )

        if prior_line == item_line:
            findings.append(
                Finding(
                    path,
                    item_line,
                    item_offset + 1,
                    RULE_MULTILINE,
                    "multiline structure requires one item per line",
                ),
            )

        prior_line = item_line

        if item_offset != expected_offset:
            findings.append(
                Finding(
                    path,
                    item_line,
                    item_offset + 1,
                    RULE_MULTILINE,
                    "multiline item requires one continuation indent",
                ),
            )


def _check_one_multiline_structure(
    node: ast.AST,
    open_index: int,
    close_index: int,
    significant: list[tokenize.TokenInfo],
    lines: list[str],
    path: str,
    findings: list[Finding],
) -> Counter[str]:
    """
    Check delimiters and items for one recognized multiline structure.

    Parameters
    ----------
    node : ast.AST
        Collection, call, signature, or subscript owning the token pair.
    open_index : int
        Index of the opening delimiter in the significant token stream.
    close_index : int
        Index of the matching closing delimiter.
    significant : list[tokenize.TokenInfo]
        Nontrivia tokens for the complete source.
    lines : list[str]
        Physical source lines used for indentation checks.
    path : str
        Repository-relative path used in diagnostics.
    findings : list[Finding]
        Mutable deterministic finding collection.

    Returns
    -------
    counts : Counter[str]
        Delimiter, indentation, item, and trailing-comma facts.
    """

    counts: Counter[str] = Counter(recognized=1)
    opening = significant[open_index]
    closing = significant[close_index]
    inner = significant[open_index + 1 : close_index]

    if not inner:
        return counts

    if inner[0].start[0] == opening.start[0]:
        _finding(
            findings,
            path,
            opening,
            RULE_MULTILINE,
            "multiline structure must break after the opening delimiter",
        )

    if inner[-1].end[0] == closing.start[0]:
        _finding(
            findings,
            path,
            closing,
            RULE_MULTILINE,
            "multiline structure must break before the closing delimiter",
        )

    items = _delimited_items(node)
    compact_row = _compact_call_argument_row(
        node,
        items,
        opening,
        closing,
        lines,
    )

    if compact_row is None and inner[-1].string != ",":
        _finding(
            findings,
            path,
            closing,
            RULE_MULTILINE,
            "recognized multiline structure requires a trailing comma",
        )
    elif compact_row is not None and inner[-1].string == ",":
        _finding(
            findings,
            path,
            closing,
            RULE_MULTILINE,
            "compact multiline call row must omit its trailing comma",
        )

    construction_line = lines[opening.start[0] - 1]
    construction_indent = len(construction_line) - len(
        construction_line.lstrip(" "),
    )

    if closing.start[1] != construction_indent:
        _finding(
            findings,
            path,
            closing,
            RULE_MULTILINE,
            "closing delimiter must align with the construction line",
        )

    expected_offset = construction_indent + 4

    if compact_row is not None:
        first_offset = items[0][1]

        if first_offset != expected_offset:
            findings.append(
                Finding(
                    path,
                    compact_row,
                    first_offset + 1,
                    RULE_MULTILINE,
                    "compact call row requires one continuation indent",
                ),
            )

        counts["compact_call_rows"] += 1

        return counts

    _check_multiline_item_layout(
        items,
        expected_offset,
        lines,
        path,
        findings,
    )

    return counts


def _check_multiline(
    tree: ast.Module,
    tokens: list[tokenize.TokenInfo],
    lines: list[str],
    path: str,
    findings: list[Finding],
) -> Counter[str]:
    """
    Check recognized multiline calls and collection displays.

    Parameters
    ----------
    tree : ast.Module
        Parsed syntax tree for the maintained source.
    tokens : list[tokenize.TokenInfo]
        Token stream paired with the parsed source.
    lines : list[str]
        Physical source lines used for location-aware checks.
    path : str
        Repository-relative path associated with the source.
    findings : list[Finding]
        Mutable finding collection populated by the check.

    Returns
    -------
    counts : Counter[str]
        Counts for recognized multiline constructions.
    """

    counts: Counter[str] = Counter()
    significant = _significant(tokens)
    opening_to_closing, matched = _matched_delimiters(significant)
    counts.update(
        _check_detached_comparisons(
            significant,
            path,
            findings,
        ),
    )
    counts.update(
        _check_parameter_lists(
            tree,
            significant,
            lines,
            path,
            findings,
            opening_to_closing,
        ),
    )
    counts.update(
        _check_simple_return_annotations(
            tree,
            significant,
            path,
            findings,
            opening_to_closing,
        ),
    )
    counts.update(
        _check_compound_condition_readability(
            tree,
            path,
            findings,
        ),
    )
    checked_pairs: set[tuple[int, int]] = set()

    for node in ast.walk(tree):
        if (
            not isinstance(node, DELIMITED)
            or node.lineno == node.end_lineno
            or node.end_lineno is None
            or node.end_col_offset is None
        ):
            continue

        if (
            isinstance(node, ast.Call)
            and len(node.args) == 1
            and not node.keywords
            and isinstance(node.args[0], ast.GeneratorExp)
        ):
            counts["generator_call_exclusions"] += 1

            continue

        pair = _multiline_token_pair(
            node,
            significant,
            matched,
            checked_pairs,
        )

        if pair is None:
            continue

        open_index, close_index = pair
        counts.update(
            _check_one_multiline_structure(
                node,
                open_index,
                close_index,
                significant,
                lines,
                path,
                findings,
            ),
        )

    return counts


def _check_line_length(
    lines: list[str],
) -> Counter[str]:
    """
    Inventory physical source lines above 79 columns.
    """

    counts: Counter[str] = Counter()

    for line in lines:
        if len(line) <= 79:
            continue

        counts["over_79_raw"] += 1

        if "http://" in line or "https://" in line:
            counts["url_exceptions"] += 1

            continue

        counts["review_candidates"] += 1

    return counts


def analyze_text(text: str, path: str = "<memory>") -> Analysis:
    """
    Analyze one Python source string without changing it.
    """

    tree = ast.parse(text, filename=path, type_comments=True)
    tokens = _tokens(text)
    lines = text.splitlines()
    parents = _parents(tree)
    findings: list[Finding] = []

    doc_positions, docstrings = _check_docstrings(
        tree,
        tokens,
        lines,
        path,
        findings,
    )
    numpy_docstrings = _check_numpy_docstrings(
        tree,
        tokens,
        path,
        findings,
    )

    facts = {
        "path": path,
        "source_fingerprint": source_fingerprint(text),
        "lines": len(lines),
        "docstrings": dict(sorted(docstrings.items())),
        "numpy_docstrings": dict(sorted(numpy_docstrings.items())),
        "annotations": dict(
            sorted(
                _check_annotations(
                    tree,
                    parents,
                    path,
                    findings,
                ).items(),
            ),
        ),
        "strings": dict(
            sorted(
                _check_strings(
                    tokens,
                    doc_positions,
                    path,
                    findings,
                ).items(),
            ),
        ),
        "cli_help_layout": dict(
            sorted(
                _check_cli_help_layout(
                    tree,
                    parents,
                    tokens,
                    lines,
                    path,
                    findings,
                ).items(),
            ),
        ),
        "comments": dict(
            sorted(
                _check_comments(
                    tokens,
                    lines,
                    path,
                    findings,
                ).items(),
            ),
        ),
        "naming": dict(sorted(_check_naming(tree, path, findings).items())),
        "topology": dict(
            sorted(
                _check_topology(
                    tree,
                    parents,
                    lines,
                    path,
                    findings,
                ).items(),
            ),
        ),
        "multiline": dict(
            sorted(
                (
                    _check_multiline(
                        tree,
                        tokens,
                        lines,
                        path,
                        findings,
                    )
                    + _check_block_literal_boundaries(
                        tokens,
                        doc_positions,
                        lines,
                        path,
                        findings,
                    )
                    + _check_adjacent_prose_literals(
                        tree,
                        parents,
                        tokens,
                        doc_positions,
                        lines,
                        path,
                        findings,
                    )
                ).items(),
            ),
        ),
        "line_length": dict(sorted(_check_line_length(lines).items())),
    }

    return Analysis(tuple(sorted(set(findings))), facts)


def scan(
    root: Path,
    paths: Iterable[str] | None = None,
) -> tuple[list[Finding], list[dict[str, object]]]:
    """
    Analyze explicit paths or the complete maintained Python universe.
    """

    root = root.resolve()
    selected = (
        sorted(set(paths))
        if paths is not None
        else [
            path.relative_to(root).as_posix()
            for path in maintained_python_paths(root)
        ]
    )
    findings: list[Finding] = []
    facts: list[dict[str, object]] = []

    for path in selected:
        target = root / path
        if not target.is_file():
            raise FileNotFoundError(path)

        analysis = analyze_text(target.read_text(encoding="utf-8"), path)
        findings.extend(analysis.findings)
        facts.append(analysis.facts)

    return sorted(set(findings)), facts


def report(
    root: Path,
    paths: Iterable[str] | None = None,
    rule_id: str | None = None,
) -> dict[str, object]:
    """
    Return one stable JSON-ready checker report.
    """

    findings, facts = scan(root, paths)

    if rule_id is not None:
        findings = [
            finding for finding in findings if finding.rule_id == rule_id
        ]

    return {
        "schema_version": 1,
        "producer": "python_source_policy",
        "producer_version": VERSION,
        "paths": [fact["path"] for fact in facts],
        "facts": facts,
        "finding_count": len(findings),
        "findings": [dataclasses.asdict(item) for item in findings],
    }


def _pilot_paths(path: Path) -> tuple[str, ...]:
    """
    Return the configured bounded pilot paths.
    """

    value = json.loads(path.read_text(encoding="utf-8"))
    paths = value.get("pilot_paths")
    if (
        not isinstance(paths, list)
        or not paths
        or not all(isinstance(item, str) and item for item in paths)
        or len(paths) != len(set(paths))
    ):
        raise ValueError("pilot_paths must be a nonempty unique string list")

    return tuple(paths)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse command-line arguments.

    Parameters
    ----------
    argv : list[str] | None
        Explicit arguments, or None to read the process arguments.

    Returns
    -------
    arguments : argparse.Namespace
        Parsed source-policy cohort and rule selectors.
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--json", action="store_true")
    cohort = parser.add_mutually_exclusive_group()
    cohort.add_argument("--all-maintained", action="store_true")
    cohort.add_argument("--pilot", action="store_true")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--rule", choices=RULE_IDS)
    parser.add_argument("--output", type=Path)
    parser.add_argument("paths", nargs="*")

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """
    Run the deterministic Python source-policy checker.

    Parameters
    ----------
    argv : list[str] | None
        Explicit arguments, or None to read the process arguments.

    Returns
    -------
    status : int
        Zero when selected policy checks pass and one otherwise.
    """

    args = parse_args(argv)

    if args.paths and (args.all_maintained or args.pilot):
        raise ValueError(
            "explicit paths cannot be combined with a cohort selector",
        )

    if args.paths:
        paths: Iterable[str] | None = args.paths
    elif args.pilot:
        paths = _pilot_paths(args.config)
    else:
        paths = None

    payload = report(args.root, paths, args.rule)

    if args.json or args.output is not None:
        rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"

        if args.output is not None:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(rendered, encoding="utf-8")
        else:
            print(rendered, end="")
    else:
        for finding in payload["findings"]:
            item = Finding(**finding)
            print(item.format())

        print(
            f"PY.SOURCE.POLICY: {len(payload['paths'])} path(s), "
            f"{payload['finding_count']} deterministic finding(s)",
        )

    return 1 if payload["finding_count"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
