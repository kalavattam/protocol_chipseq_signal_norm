#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: boolean_contracts.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Inventory and enforce source-derived Boolean-like contracts.
"""

from __future__ import annotations

import argparse
import ast
import dataclasses
import json
import re
from collections import Counter
from pathlib import Path

from dev.audit.help_aliases import function_body_spans, parameter_rows
from dev.audit.help_heredoc_reflow import extract_help_heredocs
from dev.audit.help_style import PARAMETER_ENTRY, sections

TRUE_TOKENS = ("true", "t", "yes", "y", "1")
FALSE_TOKENS = ("false", "f", "no", "n", "0")
POLICY_TEXT = (
    "Accepted true-like values: true, t, yes, y, 1. Accepted false-like "
    "values: false, f, no, n, 0. Matching is case-insensitive; surrounding "
    "whitespace and empty values are invalid. Successful normalization emits "
    "'true' or 'false'."
)
ENV_POLICY_TEXT = (
    "Boolean test gates accept true, t, yes, y, 1 and false, f, no, n, 0 "
    "case-insensitively. Unset or empty gates are disabled; surrounding "
    "whitespace and other nonempty values are invalid."
)

RULE_MANUAL = "BOOLEAN.NORMALIZATION.MANUAL"
RULE_REQUIRED = "BOOLEAN.NORMALIZATION.REQUIRED"
RULE_DOCUMENTATION = "BOOLEAN.DOCUMENTATION.EXACT_TOKENS"
RULE_ENVIRONMENT = "BOOLEAN.ENVIRONMENT.CANONICAL"
RULE_ENVIRONMENT_DOC = "BOOLEAN.ENVIRONMENT.DOCUMENTATION"
RULE_FLAG = "BOOLEAN.PRESENCE_FLAG.VALUE_FREE"
RULE_HELPER = "BOOLEAN.NORMALIZATION.HELPER_CONTRACT"
RULE_WHITESPACE = "BOOLEAN.NORMALIZATION.NO_TRIM"

MANUAL_TRUE = "true|t|yes|y|1"
MANUAL_FALSE = "false|f|no|n|0"
ENV_NAME = re.compile(r"\b(?:RUN_[A-Z0-9_]+|WAIT_SLURM)\b")
RAW_ENV_BOOL = re.compile(
    r"\$\{(?P<name>(?:RUN_[A-Z0-9_]+|WAIT_SLURM)):-0\}"
    r'[^\n]*(?:==|!=)\s*"?1"?',
)
NORMALIZED_ENV = re.compile(
    r"\bnormalize_test_gate\s+"
    r"(?P<name>(?:RUN_[A-Z0-9_]+|WAIT_SLURM))\b",
)
NORMALIZE_CALL = re.compile(
    r"\bnormalize_bool\s+"
    r'"?\$\{(?P<name>[A-Za-z_][A-Za-z0-9_]*)\}"?',
)
POSITIONAL_HEAD = re.compile(
    r"(?P<ordinal>\d+)\+?\s+(?P<name>[A-Za-z_][A-Za-z0-9_]*)$",
)


@dataclasses.dataclass(frozen=True)
class Finding:
    """
    One stable Boolean-contract violation.
    """

    rule_id: str
    path: str
    line: int
    owner: str
    message: str

    def format(self) -> str:
        """
        Render one deterministic line-oriented diagnostic.
        """

        return (
            f"{self.rule_id}: {self.path}:{self.line}: "
            f"owner={self.owner}; {self.message}"
        )


@dataclasses.dataclass(frozen=True)
class Contract:
    """
    One source-derived Boolean-like contract or reviewed adjacent mode.
    """

    owner_identity: str
    path: str
    line: int
    public_status: str
    input_surface: str
    name: str
    classification: str
    value_form: str
    accepted_true_like: tuple[str, ...]
    accepted_false_like: tuple[str, ...]
    case_sensitivity: str
    whitespace_behavior: str
    empty_unset_behavior: str
    invalid_value_behavior: str
    canonical_output: str
    error_return_contract: str
    normalization_helper: str
    direct_callers: tuple[str, ...]
    current_help_text: str
    relevant_tests: tuple[str, ...]
    proposed_disposition: str

    def as_dict(self) -> dict[str, object]:
        """
        Return a JSON-ready inventory record.
        """

        return dataclasses.asdict(self)


def shell_sources(root: Path) -> list[Path]:
    """
    Return maintained shell production and test/helper sources.
    """

    candidates = [
        *(root / "bin").glob("*.sh"),
        *(root / "lib/bash").rglob("*.sh"),
        *(root / "install/scripts").glob("*.sh"),
        *(root / "tests").rglob("*.sh"),
    ]

    return sorted(
        path
        for path in set(candidates)
        if "artifacts/tests" not in path.as_posix()
    )


def python_sources(root: Path) -> list[Path]:
    """
    Return maintained Python production and test/helper sources.
    """

    candidates = [
        *(root / "src").rglob("*.py"),
        *(root / "dev").rglob("*.py"),
        *(root / "tests").rglob("*.py"),
    ]

    return sorted(
        path
        for path in set(candidates)
        if "__pycache__" not in path.parts
        and "artifacts/tests" not in path.as_posix()
    )


def descriptions_by_row(
    text: str,
) -> dict[tuple[str, int], tuple[str, str, str]]:
    """
    Reuse recognized help sections to map typed rows to descriptions.
    """

    result: dict[tuple[str, int], tuple[str, str, str]] = {}

    for heredoc in extract_help_heredocs(text):
        for section in sections(heredoc):
            if section.name not in {"Parameters", "Expected globals"}:
                continue

            body = list(section.lines)

            for index, (number, line) in enumerate(body):
                match = PARAMETER_ENTRY.match(line)

                if match is None and section.name == "Expected globals":
                    match = re.match(
                        r"^(?P<indent> *)(?P<head>[A-Za-z_][A-Za-z0-9_]*(?:, "
                        r"[A-Za-z_][A-Za-z0-9_]*)*)\s+:\s+"
                        r"(?P<type>\S.*)$",
                        line,
                    )

                if match is None:
                    continue

                description: list[str] = []

                for _, candidate in body[index + 1 :]:
                    if not candidate.strip():
                        break

                    if PARAMETER_ENTRY.match(candidate):
                        break

                    description.append(candidate.strip())

                result[(heredoc.owner, number)] = (
                    section.name,
                    match.group("head").strip(),
                    " ".join(description),
                )

    return result


def lexical_owner(text: str, line: int) -> str:
    """
    Return the bounded function owner for one source line.
    """

    candidates = [
        (name, start, end)
        for name, (start, end, _) in function_body_spans(text).items()
        if start <= line <= end
    ]

    return candidates[-1][0] if candidates else "<file>"


def function_normalizes(
    owner: str,
    bodies: dict[str, str],
    seen: frozenset[str] = frozenset(),
) -> bool:
    """
    Return whether an owner reaches the canonical helper in one file.
    """

    if owner == "normalize_bool":
        return True

    if owner in seen or owner not in bodies:
        return False

    body = bodies[owner]
    if re.search(r"\bnormalize_bool\b", body):
        return True

    nested_seen = seen | {owner}

    return any(
        re.search(rf"(?m)^\s*{re.escape(candidate)}\b", body)
        and function_normalizes(candidate, bodies, nested_seen)
        for candidate in bodies
        if candidate != owner
    )


def public_status(path: str, owner: str) -> str:
    """
    Classify repository interface visibility from its maintained location.
    """

    if owner.startswith("_"):
        return "internal"

    if path.startswith("lib/bash/") or path.startswith("tests/support/"):
        return "public"

    return "internal"


def call_sites(
    root: Path,
    sources: list[Path],
    owner: str,
    own_identity: str,
) -> tuple[str, ...]:
    """
    Return direct lexical callers of one shell function name.
    """

    callers: set[str] = set()
    pattern = re.compile(rf"(?m)^\s*{re.escape(owner)}(?:\s|$)")

    for source in sources:
        text = source.read_text(encoding="utf-8")
        rel = str(source.relative_to(root))

        for match in pattern.finditer(text):
            line = text.count("\n", 0, match.start()) + 1
            identity = f"{rel}::{lexical_owner(text, line)}"

            if identity != own_identity:
                callers.add(identity)

    return tuple(sorted(callers))


def build_call_index(
    root: Path,
    sources: list[Path],
) -> dict[str, tuple[str, ...]]:
    """
    Index direct shell callers once for deterministic inventory reuse.
    """

    cache = {
        str(source.relative_to(root)): source.read_text(encoding="utf-8")
        for source in sources
    }
    known = sorted(
        {
            owner
            for text in cache.values()
            for owner in function_body_spans(text)
        },
    )
    callers: dict[str, set[str]] = {owner: set() for owner in known}
    known_set = set(known)

    for path, text in cache.items():
        spans = function_body_spans(text)
        regions = [
            ("<file>", text),
            *((owner, body) for owner, (_, _, body) in spans.items()),
        ]

        for lexical, body in regions:
            identity = f"{path}::{lexical}"
            invoked = {
                match.group("name")
                for match in re.finditer(
                    r"(?m)^\s*(?P<name>[A-Za-z_][A-Za-z0-9_]*)\b",
                    body,
                )
                if match.group("name") in known_set
            }

            for owner in invoked - {lexical}:
                callers[owner].add(identity)

    return {
        owner: tuple(sorted(identities))
        for owner, identities in callers.items()
    }


def canonical_fields(classification: str) -> dict[str, object]:
    """
    Return common contract fields for one classification.
    """

    if classification == "presence-only flag":
        return {
            "accepted_true_like": (),
            "accepted_false_like": (),
            "case_sensitivity": "not applicable",
            "whitespace_behavior": "no value consumed",
            "empty_unset_behavior": "omitted means false/default",
            "invalid_value_behavior": "following token is not a Boolean value",
            "canonical_output": "present means true behavior",
            "error_return_contract": "normal option-parser contract",
            "normalization_helper": "none",
            "proposed_disposition": "retain as a value-free presence flag",
        }

    if classification == "internal already-canonical Boolean":
        return {
            "accepted_true_like": ("true",),
            "accepted_false_like": ("false",),
            "case_sensitivity": "canonical lowercase only",
            "whitespace_behavior": "not accepted by internal contract",
            "empty_unset_behavior": "owner-specific default or invalid",
            "invalid_value_behavior": "not a public normalization surface",
            "canonical_output": "true or false",
            "error_return_contract": "caller establishes canonical value",
            "normalization_helper": "upstream presence/parser contract",
            "proposed_disposition": "retain as internal canonical transport",
        }

    if classification == "Boolean environment variable":
        return {
            "accepted_true_like": TRUE_TOKENS,
            "accepted_false_like": FALSE_TOKENS,
            "case_sensitivity": "case-insensitive",
            "whitespace_behavior": "surrounding whitespace invalid",
            "empty_unset_behavior": "disabled false state",
            "invalid_value_behavior": "diagnostic and nonzero termination",
            "canonical_output": "true or false",
            "error_return_contract": "invalid nonempty value returns nonzero",
            "normalization_helper": "normalize_test_gate -> normalize_bool",
            "proposed_disposition": "normalize nonempty values canonically",
        }

    return {
        "accepted_true_like": TRUE_TOKENS,
        "accepted_false_like": FALSE_TOKENS,
        "case_sensitivity": "case-insensitive",
        "whitespace_behavior": "surrounding whitespace invalid",
        "empty_unset_behavior": "empty required value invalid",
        "invalid_value_behavior": "precise diagnostic and nonzero return",
        "canonical_output": "true or false",
        "error_return_contract": (
            "invalid value returns nonzero without fallback"
        ),
        "normalization_helper": "normalize_bool",
        "proposed_disposition": "route through canonical normalization",
    }


def shell_inventory(root: Path, sources: list[Path]) -> list[Contract]:
    """
    Build shell contracts from existing help and function parsers.

    Parameters
    ----------
    root : Path
        Repository root used to resolve source identities.
    sources : list[Path]
        Maintained Shell sources to inspect.

    Returns
    -------
    contracts : list[Contract]
        Source-derived Boolean contract records.
    """

    inventory: list[Contract] = []
    explicit_keys: set[tuple[str, str]] = set()

    caller_index = build_call_index(root, sources)

    for source in sources:
        path = str(source.relative_to(root))
        text = source.read_text(encoding="utf-8")
        bodies = {
            owner: body
            for owner, (_, _, body) in function_body_spans(text).items()
        }
        descriptions = descriptions_by_row(text)

        for row in parameter_rows(text):
            if row.type_name.strip() != "flag":
                continue

            owner_identity = f"{path}::{row.owner}"
            help_text = descriptions.get((row.owner, row.line), ("", "", ""))[
                2
            ]
            fields = canonical_fields("presence-only flag")
            inventory.append(
                Contract(
                    owner_identity=owner_identity,
                    path=path,
                    line=row.line,
                    public_status=public_status(path, row.owner),
                    input_surface="shell option",
                    name=",".join(row.aliases),
                    classification="presence-only flag",
                    value_form="no-value flag",
                    direct_callers=caller_index.get(row.owner, ()),
                    current_help_text=help_text,
                    relevant_tests=(
                        "tests/contract/repository/test_boolean_contracts.sh",
                    ),
                    **fields,
                ),
            )

        for heredoc in extract_help_heredocs(text):
            for section in sections(heredoc):
                if section.name not in {"Parameters", "Expected globals"}:
                    continue

                body = list(section.lines)

                for _index, (line, source_line) in enumerate(body):
                    match = PARAMETER_ENTRY.match(source_line)

                    if match is None and section.name == "Expected globals":
                        match = re.match(
                            r"^(?P<indent> *)(?P<head>[A-Za-z_][A-Za-z0-9_]*"
                            r"(?:, [A-Za-z_][A-Za-z0-9_]*)*)\s+:\s+"
                            r"(?P<type>\S.*)$",
                            source_line,
                        )

                    if match is None:
                        continue

                    boolean_type = match.group("type").strip() in {
                        "bool",
                        "flag",
                    }

                    if not boolean_type:
                        continue

                    head = match.group("head").strip()
                    if head.startswith("-"):
                        continue

                    description = descriptions.get(
                        (heredoc.owner, line),
                        ("", "", ""),
                    )[2]
                    positional = POSITIONAL_HEAD.match(head)
                    names = (
                        (positional.group("name"),)
                        if positional
                        else tuple(part.strip() for part in head.split(","))
                    )

                    normalizes = function_normalizes(heredoc.owner, bodies)
                    visibility = public_status(path, heredoc.owner)
                    explicit = (
                        match.group("type").strip() == "bool"
                        and section.name == "Parameters"
                        and (visibility == "public" or normalizes)
                    )
                    classification = (
                        "Boolean positional argument"
                        if explicit
                        else "internal already-canonical Boolean"
                    )

                    for name in names:
                        owner_identity = f"{path}::{heredoc.owner}"
                        fields = canonical_fields(classification)
                        inventory.append(
                            Contract(
                                owner_identity=owner_identity,
                                path=path,
                                line=line,
                                public_status=visibility,
                                input_surface=(
                                    "function positional argument"
                                    if section.name == "Parameters"
                                    else "expected global"
                                ),
                                name=name,
                                classification=classification,
                                value_form="explicit Boolean value",
                                direct_callers=caller_index.get(
                                    heredoc.owner,
                                    (),
                                ),
                                current_help_text=description,
                                relevant_tests=(
                                    (
                                        "tests/contract/repository/test_boolea"
                                        "n_contracts.sh"
                                    ),
                                ),
                                **fields,
                            ),
                        )

                        if explicit:
                            explicit_keys.add((owner_identity, name))

        for owner, body in bodies.items():
            if owner in {"normalize_bool", "normalize_test_gate"}:
                continue

            for match in NORMALIZE_CALL.finditer(body):
                name = match.group("name")
                owner_identity = f"{path}::{owner}"
                if (owner_identity, name) in explicit_keys:
                    continue

                line = (
                    function_body_spans(text)[owner][0]
                    + body.count("\n", 0, match.start())
                    + 1
                )
                fields = canonical_fields("Boolean positional argument")
                inventory.append(
                    Contract(
                        owner_identity=owner_identity,
                        path=path,
                        line=line,
                        public_status=public_status(path, owner),
                        input_surface=(
                            "conditional function positional argument"
                        ),
                        name=name,
                        classification="Boolean positional argument",
                        value_form="explicit Boolean value",
                        direct_callers=caller_index.get(owner, ()),
                        current_help_text=(
                            "Boolean-like branch documented in owner help"
                        ),
                        relevant_tests=(
                            (
                                "tests/contract/repository/test_boolean_contra"
                                "cts.sh"
                            ),
                        ),
                        **fields,
                    ),
                )
                explicit_keys.add((owner_identity, name))

    return inventory


def python_flag_inventory(root: Path) -> list[Contract]:
    """
    Inventory argparse presence flags without treating them as values.
    """

    inventory: list[Contract] = []

    for source in python_sources(root):
        path = str(source.relative_to(root))

        try:
            tree = ast.parse(source.read_text(encoding="utf-8"))
        except SyntaxError:
            continue

        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue

            if (
                not isinstance(node.func, ast.Attribute)
                or node.func.attr != "add_argument"
            ):
                continue

            action = next(
                (
                    keyword.value.value
                    for keyword in node.keywords
                    if keyword.arg == "action"
                    and isinstance(keyword.value, ast.Constant)
                ),
                None,
            )
            if action not in {"store_true", "store_false"}:
                continue

            aliases = tuple(
                value.value
                for value in node.args
                if isinstance(value, ast.Constant)
                and isinstance(value.value, str)
            )
            fields = canonical_fields("presence-only flag")
            inventory.append(
                Contract(
                    owner_identity=f"{path}::parse_args",
                    path=path,
                    line=node.lineno,
                    public_status="public",
                    input_surface="Python argparse option",
                    name=",".join(aliases),
                    classification="presence-only flag",
                    value_form="no-value flag",
                    direct_callers=(f"{path}::main",),
                    current_help_text="argparse action=" + str(action),
                    relevant_tests=(
                        "tests/contract/repository/test_boolean_contracts.sh",
                    ),
                    **fields,
                ),
            )

    return inventory


def environment_inventory(root: Path, sources: list[Path]) -> list[Contract]:
    """
    Inventory live Boolean test gates and their disabled empty state.
    """

    observed: dict[str, tuple[str, int, str]] = {}

    for source in sources:
        path = str(source.relative_to(root))
        text = source.read_text(encoding="utf-8")

        for pattern in (RAW_ENV_BOOL, NORMALIZED_ENV):
            for match in pattern.finditer(text):
                name = match.group("name")
                if (
                    name == "WAIT_SLURM"
                    and path == "tests/integration/slurm/run_wet_tests.sh"
                ):
                    continue

                line = text.count("\n", 0, match.start()) + 1
                if text.splitlines()[line - 1].lstrip().startswith("#"):
                    continue

                observed.setdefault(
                    name,
                    (path, line, lexical_owner(text, line)),
                )

    inventory: list[Contract] = []
    fields = canonical_fields("Boolean environment variable")

    for name, (path, line, owner) in sorted(observed.items()):
        inventory.append(
            Contract(
                owner_identity=f"{path}::{owner}",
                path=path,
                line=line,
                public_status="public",
                input_surface="environment variable",
                name=name,
                classification="Boolean environment variable",
                value_form="optional explicit Boolean value",
                direct_callers=(),
                current_help_text="documented test gate",
                relevant_tests=(
                    "tests/contract/repository/test_boolean_contracts.sh",
                ),
                **fields,
            ),
        )

    return inventory


def adjacent_mode_inventory(root: Path) -> list[Contract]:
    """
    Disposition Boolean-adjacent generic selectors as non-Boolean modes.
    """

    path = root / "lib/bash/workflows/process_tables.sh"
    if not path.is_file():
        return []

    text = path.read_text(encoding="utf-8")
    rows: list[Contract] = []

    for owner, name, marker in (
        ("check_table_scaling_factor", "type", "str|string"),
        ("_validate_arg_csv", "valid", "IFS=',' read"),
    ):
        span = function_body_spans(text).get(owner)
        if span is None or marker not in span[2]:
            continue

        rows.append(
            Contract(
                owner_identity=(
                    f"lib/bash/workflows/process_tables.sh::{owner}"
                ),
                path="lib/bash/workflows/process_tables.sh",
                line=span[0],
                public_status=public_status(
                    "lib/bash/workflows/process_tables.sh",
                    owner,
                ),
                input_surface="function positional argument",
                name=name,
                classification="non-Boolean enum or mode",
                value_form="choice selector",
                accepted_true_like=(),
                accepted_false_like=(),
                case_sensitivity="owner-specific case-insensitive choices",
                whitespace_behavior="not a Boolean contract",
                empty_unset_behavior="required by owner",
                invalid_value_behavior="owner-specific choice-set diagnostic",
                canonical_output="owner-specific mode",
                error_return_contract="invalid choice returns nonzero",
                normalization_helper="none",
                direct_callers=(),
                current_help_text="generic selector, not a Boolean value",
                relevant_tests=(),
                proposed_disposition="exclude from Boolean normalization",
            ),
        )

    return rows


def inventory_repository(root: Path) -> list[Contract]:
    """
    Return the complete deterministic source-derived inventory.
    """

    root = root.resolve()
    sources = shell_sources(root)
    rows = [
        *shell_inventory(root, sources),
        *python_flag_inventory(root),
        *environment_inventory(root, sources),
        *adjacent_mode_inventory(root),
    ]
    unique = {
        (
            row.owner_identity,
            row.line,
            row.name,
            row.classification,
        ): row
        for row in rows
    }

    return sorted(
        unique.values(),
        key=lambda row: (
            row.path,
            row.line,
            row.owner_identity,
            row.name,
            row.classification,
        ),
    )


def scan_repository(root: Path) -> tuple[list[Finding], list[Contract]]:
    """
    Audit normalization, documentation, and classification contracts.

    Parameters
    ----------
    root : Path
        Repository root used to resolve maintained paths.

    Returns
    -------
    ordered, inventory : tuple[list[Finding], list[Contract]]
        Repository findings and supporting inventory.
    """

    root = root.resolve()
    sources = shell_sources(root)

    inventory = inventory_repository(root)
    findings: list[Finding] = []

    for source in sources:
        path = str(source.relative_to(root))
        text = source.read_text(encoding="utf-8")
        spans = function_body_spans(text)
        bodies = {owner: body for owner, (_, _, body) in spans.items()}

        for owner, (start, _, body) in spans.items():
            if owner == "normalize_bool":
                continue

            if MANUAL_TRUE in body or MANUAL_FALSE in body:
                offset = min(
                    value
                    for value in (
                        body.find(MANUAL_TRUE),
                        body.find(MANUAL_FALSE),
                    )
                    if value >= 0
                )
                line = start + body.count("\n", 0, offset) + 1
                findings.append(
                    Finding(
                        RULE_MANUAL,
                        path,
                        line,
                        owner,
                        "manual Boolean token cases must use normalize_bool",
                    ),
                )

        for line, source_line in enumerate(text.splitlines(), 1):
            if "normalize_bool" not in source_line:
                continue

            if re.search(r"\b(?:xargs|sed|awk|strip|trim)\b", source_line):
                findings.append(
                    Finding(
                        RULE_WHITESPACE,
                        path,
                        line,
                        lexical_owner(text, line),
                        (
                            "do not trim or transform whitespace before "
                            "normalize_bool"
                        ),
                    ),
                )

        for match in RAW_ENV_BOOL.finditer(text):
            line = text.count("\n", 0, match.start()) + 1
            if text.splitlines()[line - 1].lstrip().startswith("#"):
                continue

            findings.append(
                Finding(
                    RULE_ENVIRONMENT,
                    path,
                    line,
                    lexical_owner(text, line),
                    (
                        f"{match.group('name')} accepts only numeric 1; use "
                        f"normalize_test_gate"
                    ),
                ),
            )

        for heredoc in extract_help_heredocs(text):
            body = "\n".join(line for _, line in heredoc.lines)
            bool_rows = [
                contract
                for contract in inventory
                if contract.owner_identity == f"{path}::{heredoc.owner}"
                and contract.classification == "Boolean positional argument"
            ]

            if bool_rows and POLICY_TEXT not in body:
                findings.append(
                    Finding(
                        RULE_DOCUMENTATION,
                        path,
                        bool_rows[0].line,
                        heredoc.owner,
                        (
                            "explicit Boolean help must state the exact "
                            "canonical policy"
                        ),
                    ),
                )

        for contract in inventory:
            if contract.path != path:
                continue

            if contract.classification != "Boolean positional argument":
                continue

            owner = contract.owner_identity.rsplit("::", 1)[1]

            if owner in bodies and not function_normalizes(owner, bodies):
                findings.append(
                    Finding(
                        RULE_REQUIRED,
                        path,
                        contract.line,
                        owner,
                        (
                            f"explicit Boolean '{contract.name}' does not "
                            f"reach normalize_bool"
                        ),
                    ),
                )

    helper = root / "lib/bash/core/check_args.sh"
    helper_text = (
        helper.read_text(encoding="utf-8") if helper.is_file() else ""
    )
    helper_span = function_body_spans(helper_text).get("normalize_bool")
    helper_body = helper_span[2] if helper_span else ""
    required_fragments = (
        MANUAL_TRUE,
        MANUAL_FALSE,
        'echo "true"',
        'echo "false"',
        '[[ -z "${val}" ]]',
        'val_lc="${val,,}"',
    )

    for fragment in required_fragments:
        if fragment not in helper_body:
            findings.append(
                Finding(
                    RULE_HELPER,
                    "lib/bash/core/check_args.sh",
                    helper_span[0] if helper_span else 1,
                    "normalize_bool",
                    (
                        f"canonical helper contract missing source fragment: "
                        f"{fragment}"
                    ),
                ),
            )

    if POLICY_TEXT not in helper_text:
        findings.append(
            Finding(
                RULE_DOCUMENTATION,
                "lib/bash/core/check_args.sh",
                helper_span[0] if helper_span else 1,
                "normalize_bool",
                (
                    "canonical helper help must state strict whitespace and "
                    "output policy"
                ),
            ),
        )

    test_readme = root / "tests/README.md"
    readme_text = (
        test_readme.read_text(encoding="utf-8")
        if test_readme.is_file()
        else ""
    )

    if ENV_POLICY_TEXT not in readme_text:
        findings.append(
            Finding(
                RULE_ENVIRONMENT_DOC,
                "tests/README.md",
                1,
                "test gates",
                "Boolean environment gates must document the canonical policy",
            ),
        )

    for contract in inventory:
        if contract.classification != "presence-only flag":
            continue

        if contract.value_form != "no-value flag":
            findings.append(
                Finding(
                    RULE_FLAG,
                    contract.path,
                    contract.line,
                    contract.owner_identity.rsplit("::", 1)[1],
                    "presence-only flag must not consume a Boolean value",
                ),
            )

    unique = {
        (
            finding.rule_id,
            finding.path,
            finding.line,
            finding.owner,
            finding.message,
        ): finding
        for finding in findings
    }
    ordered = sorted(
        unique.values(),
        key=lambda finding: (
            finding.path,
            finding.line,
            finding.rule_id,
            finding.owner,
            finding.message,
        ),
    )

    return ordered, inventory


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse Boolean-contract audit arguments.

    Parameters
    ----------
    argv : list[str] | None
        Explicit arguments, or None to read the process arguments.

    Returns
    -------
    arguments : argparse.Namespace
        Parsed repository and output options.
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--inventory-json", action="store_true")
    parser.add_argument("--inventory-output", type=Path)

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """
    Print stable findings, inventory counts, and optional JSON.

    Parameters
    ----------
    argv : list[str] | None
        Explicit arguments, or None to read the process arguments.

    Returns
    -------
    status : int
        Zero when the audit passes and one when findings remain.
    """

    args = parse_args(argv)
    findings, inventory = scan_repository(args.root)
    payload = [contract.as_dict() for contract in inventory]

    if args.inventory_output:
        args.inventory_output.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    if args.inventory_json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        for finding in findings:
            print(finding.format())

        counts = Counter(contract.classification for contract in inventory)
        print(f"BOOLEAN.CONTRACTS: total={len(inventory)}")

        for classification in sorted(counts):
            print(
                (
                    f"BOOLEAN.CONTRACTS: "
                    f"{classification}={counts[classification]}"
                ),
            )

        print(f"BOOLEAN.CONTRACTS: findings={len(findings)}")

    return 1 if findings else 0


if __name__ == "__main__":
    raise SystemExit(main())
