#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: python_source_evidence.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Produce nonblocking evidence for the bounded Python source-policy pilot.
"""

from __future__ import annotations

import argparse
import ast
import dataclasses
import hashlib
import io
import json
import re
import statistics
import tokenize
from collections import Counter, defaultdict
from collections.abc import Iterable, Sequence
from pathlib import Path

from dev.audit.python_naming_vocabulary import (
    evidence_candidate_segments,
    ordinary_short_words,
)
from dev.audit.python_source_policy import (
    VERSION as CHECKER_VERSION,
)
from dev.audit.python_source_policy import (
    _is_direct_side_effect,
    _needs_control_boundary,
    _transfer_kind,
    analyze_text,
    source_fingerprint,
)
from dev.audit.python_version_policy import maintained_python_paths

VERSION = "3.0"
DEFAULT_CONFIG = (
    Path(__file__).resolve().parents[1]
    / "config"
    / "pilots"
    / "python_source_policy.json"
)
ASSIGNMENT = (ast.Assign, ast.AnnAssign, ast.AugAssign)
CALLABLE = (ast.FunctionDef, ast.AsyncFunctionDef)
DEFINITION = (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)
COMMON_SHORT_WORDS = ordinary_short_words()
KNOWN_ABBREVIATIONS = evidence_candidate_segments()
SECTION = re.compile(r"(?m)^([A-Z][A-Za-z ]+)\n-{3,}$")
RESOLVED_DISPOSITIONS = {
    "exception",
    "omitted_by_role",
    "refactored",
    "retained_coherent",
    "retained_compact_guard",
    "retained_indivisible",
}
SEMANTIC_OWNERS = {
    "annotation_correctness",
    "completed_rename_decisions",
    "density_regions",
    "docstring_applicability",
    "docstring_content",
    "naming_candidates",
    "raw_width_candidates",
    "repeated_run_members",
    "repeated_runs",
    "reusable_error_boundaries",
    "source_headers",
    "suppressions_and_e402",
    "transfers",
    "x_y_z_candidates",
}


@dataclasses.dataclass(frozen=True)
class Candidate:
    """
    Represent one stable source-layout review candidate.
    """

    signature: str
    language: str
    path: str
    enclosing_object: str
    line: int
    end_line: int
    column: int
    construct: str
    reason: str
    parameter: dict[str, int | str]
    metric: str
    unit: str
    measured_value: int
    physical_span: int
    logical_span: int
    item_count: int
    nesting: int
    source_fingerprint: str
    producer_version: str
    configuration: dict[str, object]
    disposition: str
    disposition_reason: str | None
    false_positive_risk: str
    review_key: str
    review_group: str | None


@dataclasses.dataclass(frozen=True)
class ParsedSource:
    """
    Hold one parsed maintained Python source.
    """

    path: str
    text: str
    lines: tuple[str, ...]
    tree: ast.Module
    parents: dict[ast.AST, ast.AST]
    fingerprint: str


@dataclasses.dataclass(frozen=True)
class TransferContext:
    """
    Hold the bounded predecessor and Z evidence for one transfer.
    """

    predecessor: ast.stmt | None
    blank_count: int
    preceding_phase: tuple[ast.stmt, ...]
    preceding_phase_lines: int
    substantive_preceding_phase: bool
    preparation: tuple[ast.stmt, ...]
    preparation_lines: int
    z_applies: bool
    z_exceeded: bool


@dataclasses.dataclass(frozen=True)
class ExplicitDecisionGroup:
    """
    Hold one validated explicit semantic decision group.
    """

    group_id: str
    owner: str
    disposition: str
    rationale: str
    reviewed_facts: dict[str, object]
    members: tuple[str, ...]
    membership_fingerprint: str


def _parents(tree: ast.AST) -> dict[ast.AST, ast.AST]:
    """
    Return a parent map for one AST.
    """

    return {
        child: parent
        for parent in ast.walk(tree)
        for child in ast.iter_child_nodes(parent)
    }


def _object_name(
    node: ast.AST,
    parents: dict[ast.AST, ast.AST],
) -> str:
    """
    Return the nearest dotted lexical object name.
    """

    names: list[str] = []
    current: ast.AST | None = node

    while current is not None:
        name = getattr(current, "name", None)

        if isinstance(name, str):
            names.append(name)

        current = parents.get(current)

    return ".".join(reversed(names)) or "<module>"


def _enclosing_object(
    node: ast.AST,
    parents: dict[ast.AST, ast.AST],
) -> str:
    """
    Return the nearest module, class, or callable object.
    """

    current = parents.get(node)

    while current is not None:
        if isinstance(current, (ast.ClassDef, *CALLABLE)):
            return _object_name(current, parents)

        current = parents.get(current)

    return "<module>"


def _nesting(node: ast.AST, parents: dict[ast.AST, ast.AST]) -> int:
    """
    Return compound-statement nesting depth.
    """

    compound = (
        ast.If,
        ast.For,
        ast.AsyncFor,
        ast.While,
        ast.Try,
        ast.With,
        ast.AsyncWith,
        ast.Match,
        ast.comprehension,
    )
    depth = 0
    current = parents.get(node)

    while current is not None:
        if isinstance(current, compound):
            depth += 1

        current = parents.get(current)

    return depth


def _candidate_signature(record: dict[str, object]) -> str:
    """
    Return a stable signature for one unsigned candidate record.
    """

    encoded = json.dumps(
        record,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")

    return f"sha256:{hashlib.sha256(encoded).hexdigest()}"


def _review_key(owner: str, record: dict[str, object]) -> str:
    """
    Return a stable semantic-review key for one unresolved record.

    The key includes the record's source fingerprint and semantic facts, but
    excludes generated dispositions and evidence-rendering metadata.

    Parameters
    ----------
    owner : str
        Semantic owner whose decision surface contains the record.
    record : dict[str, object]
        Unresolved record with source-derived semantic facts.

    Returns
    -------
    key : str
        SHA-256 review key that changes with owned semantic facts.
    """

    excluded = {
        "configuration",
        "correctness_disposition",
        "disposition",
        "disposition_reason",
        "members",
        "producer_version",
        "review_group",
        "review_key",
        "review_keys",
        "semantic_disposition",
        "semantic_topic_disposition",
        "signature",
    }

    def scrub(value: object) -> object:
        if isinstance(value, dict):
            return {
                str(key): scrub(item)
                for key, item in sorted(value.items())
                if key not in excluded
            }

        if isinstance(value, list):
            return [scrub(item) for item in value]

        return value

    return _candidate_signature(
        {
            "owner": owner,
            "record": scrub(record),
        },
    )


def _membership_fingerprint(members: Iterable[str]) -> str:
    """
    Return the exact identity of one explicit decision-group membership.
    """

    return _candidate_signature({"members": sorted(members)})


def _candidate(
    source: ParsedSource,
    node: ast.AST,
    *,
    end_node: ast.AST | None = None,
    construct: str,
    reason: str,
    parameter_name: str,
    parameter_value: int,
    metric: str,
    unit: str,
    measured_value: int,
    logical_span: int,
    item_count: int,
    false_positive_risk: str,
    configuration: dict[str, object],
) -> Candidate:
    """
    Build one signed candidate record.
    """

    line = getattr(node, "lineno", 1)
    final_node = end_node or node
    end_line = getattr(final_node, "end_lineno", line) or line
    unsigned: dict[str, object] = {
        "language": "python",
        "path": source.path,
        "enclosing_object": _enclosing_object(node, source.parents),
        "line": line,
        "end_line": end_line,
        "column": getattr(node, "col_offset", 0) + 1,
        "construct": construct,
        "reason": reason,
        "parameter": {
            "name": parameter_name,
            "value": parameter_value,
        },
        "metric": metric,
        "unit": unit,
        "measured_value": measured_value,
        "physical_span": end_line - line + 1,
        "logical_span": logical_span,
        "item_count": item_count,
        "nesting": _nesting(node, source.parents),
        "source_fingerprint": source.fingerprint,
        "producer_version": VERSION,
        "configuration": configuration,
        "disposition": "unresolved",
        "disposition_reason": (
            "No fresh explicit semantic decision covers this candidate."
        ),
        "false_positive_risk": false_positive_risk,
    }
    review_key = _review_key("x_y_z_candidates", unsigned)

    return Candidate(
        signature=_candidate_signature(unsigned),
        **unsigned,
        review_key=review_key,
        review_group=None,
    )


def _parse_sources(root: Path, paths: Iterable[str]) -> list[ParsedSource]:
    """
    Read and parse one stable path set.
    """

    sources: list[ParsedSource] = []

    for path in sorted(set(paths)):
        text = (root / path).read_text(encoding="utf-8")
        tree = ast.parse(text, filename=path, type_comments=True)
        sources.append(
            ParsedSource(
                path=path,
                text=text,
                lines=tuple(text.splitlines()),
                tree=tree,
                parents=_parents(tree),
                fingerprint=source_fingerprint(text),
            ),
        )

    return sources


def _combined_fingerprint(sources: Sequence[ParsedSource]) -> str:
    """
    Return an ordered content identity for one source cohort.
    """

    digest = hashlib.sha256()

    for source in sources:
        digest.update(source.path.encode("utf-8"))
        digest.update(b"\0")
        digest.update(source.fingerprint.encode("ascii"))
        digest.update(b"\0")

    return f"sha256:{digest.hexdigest()}"


def _call_item_count(node: ast.Call) -> int:
    """
    Return positional plus keyword item count for one call.
    """

    return len(node.args) + len(node.keywords)


def _x_evidence(
    sources: Sequence[ParsedSource],
    values: Sequence[int],
    configuration: dict[str, object],
) -> tuple[dict[str, int], list[Candidate]]:
    """
    Measure multiline call spans for adjacent X values.
    """

    counts = {str(value): 0 for value in values}
    candidates: list[Candidate] = []

    for source in sources:
        for node in ast.walk(source.tree):
            if not isinstance(node, ast.Call):
                continue

            span = (node.end_lineno or node.lineno) - node.lineno + 1

            for value in values:
                if span <= value:
                    continue

                counts[str(value)] += 1
                call_name = ast.unparse(node.func)
                risk = (
                    "parser-builder repetition may be one tight topic"
                    if call_name.endswith("add_argument")
                    else "literal or configuration arguments may inflate span"
                )
                candidates.append(
                    _candidate(
                        source,
                        node,
                        construct="multiline_call",
                        reason=(
                            f"physical span {span} exceeds X={value}; "
                            "semantic independence is unresolved"
                        ),
                        parameter_name="X",
                        parameter_value=value,
                        metric="multiline_call_span",
                        unit="physical_lines",
                        measured_value=span,
                        logical_span=1,
                        item_count=_call_item_count(node),
                        false_positive_risk=risk,
                        configuration=configuration,
                    ),
                )

    return counts, candidates


def _assignment_item_count(node: ast.AST) -> int:
    """
    Return a bounded item count for one assignment.
    """

    value = getattr(node, "value", None)
    if isinstance(value, (ast.List, ast.Tuple, ast.Set)):
        return len(value.elts)

    if isinstance(value, ast.Dict):
        return len(value.keys)

    if isinstance(value, ast.Call):
        return _call_item_count(value)

    return 1


def _y_evidence(
    sources: Sequence[ParsedSource],
    values: Sequence[int],
    configuration: dict[str, object],
) -> tuple[dict[str, int], list[Candidate]]:
    """
    Measure multiline assignment spans for adjacent Y values.
    """

    counts = {str(value): 0 for value in values}
    candidates: list[Candidate] = []

    for source in sources:
        for node in ast.walk(source.tree):
            if not isinstance(node, ASSIGNMENT):
                continue

            span = (node.end_lineno or node.lineno) - node.lineno + 1

            for value in values:
                if span <= value:
                    continue

                counts[str(value)] += 1
                value_node = getattr(node, "value", None)
                risk = (
                    "literal or declarative configuration may be clearer as "
                    "one assignment"
                    if isinstance(
                        value_node,
                        (ast.Dict, ast.List, ast.Set, ast.Tuple),
                    )
                    else "assignment span may reflect canonical call wrapping"
                )
                candidates.append(
                    _candidate(
                        source,
                        node,
                        construct="multiline_assignment",
                        reason=(
                            f"physical span {span} exceeds Y={value}; "
                            "visual independence is unresolved"
                        ),
                        parameter_name="Y",
                        parameter_value=value,
                        metric="multiline_assignment_span",
                        unit="physical_lines",
                        measured_value=span,
                        logical_span=1,
                        item_count=_assignment_item_count(node),
                        false_positive_risk=risk,
                        configuration=configuration,
                    ),
                )

    return counts, candidates


def _blank_between(
    source: ParsedSource,
    previous: ast.stmt,
    following: ast.stmt,
) -> bool:
    """
    Return whether a blank physical line separates two statements.
    """

    start = (previous.end_lineno or previous.lineno) + 1
    return any(
        not source.lines[number - 1].strip()
        for number in range(start, following.lineno)
    )


def _simple_assignment_runs(
    source: ParsedSource,
) -> Iterable[tuple[ast.AST, str, list[ast.stmt]]]:
    """
    Yield blank-line-delimited runs of simple assignment statements.
    """

    for owner, field, suite in _statement_suites(source.tree):
        run: list[ast.stmt] = []

        for node in [*suite, None]:
            is_simple = (
                isinstance(node, ASSIGNMENT) and node.lineno == node.end_lineno
            )
            separated = bool(
                is_simple and run and _blank_between(source, run[-1], node),
            )

            if is_simple and not separated:
                run.append(node)

                continue

            if len(run) >= 2:
                yield owner, field, run

            run = [node] if is_simple and node is not None else []


def _y_simple_run_evidence(
    sources: Sequence[ParsedSource],
    values: Sequence[int],
    configuration: dict[str, object],
) -> tuple[dict[str, int], list[Candidate]]:
    """
    Measure simple-assignment run lengths for adjacent Y values.
    """

    counts = {str(value): 0 for value in values}
    candidates: list[Candidate] = []

    for source in sources:
        for _, _, run in _simple_assignment_runs(source):
            measured = len(run)
            first = run[0]
            last = run[-1]

            for value in values:
                if measured <= value:
                    continue

                counts[str(value)] += 1
                candidates.append(
                    _candidate(
                        source,
                        first,
                        end_node=last,
                        construct="simple_assignment_run",
                        reason=(
                            f"{measured} structurally consecutive simple "
                            f"assignments exceed Y={value}; topical unity "
                            "requires semantic review"
                        ),
                        parameter_name="Y",
                        parameter_value=value,
                        metric="simple_assignment_run_length",
                        unit="statements",
                        measured_value=measured,
                        logical_span=measured,
                        item_count=measured,
                        false_positive_risk=(
                            "declarative fields or tightly related setup may "
                            "form one readable topic"
                        ),
                        configuration=configuration,
                    ),
                )

    return counts, candidates


def _statement_suites(
    node: ast.AST,
) -> Iterable[tuple[ast.AST, str, list[ast.stmt]]]:
    """
    Yield owned statement-list fields recursively.
    """

    for field, value in ast.iter_fields(node):
        if (
            isinstance(value, list)
            and value
            and all(isinstance(item, ast.stmt) for item in value)
        ):
            yield node, field, value

            for item in value:
                yield from _statement_suites(item)
        elif isinstance(value, ast.AST):
            yield from _statement_suites(value)
        elif isinstance(value, list):
            for item in value:
                if isinstance(item, ast.AST):
                    yield from _statement_suites(item)


def _blank_count_before(
    source: ParsedSource,
    node: ast.stmt,
) -> int:
    """
    Return blank lines before a statement and its attached comments.
    """

    line = node.lineno
    index = line - 2

    while index >= 0 and source.lines[index].lstrip().startswith("#"):
        line = index + 1
        index -= 1

    count = 0
    index = line - 2

    while index >= 0 and not source.lines[index].strip():
        count += 1
        index -= 1

    return count


def _phase_members(
    source: ParsedSource,
    statements: tuple[ast.stmt, ...],
) -> list[dict[str, object]]:
    """
    Return inspectable membership for one transfer's preceding phase.
    """

    members: list[dict[str, object]] = []

    for index, statement in enumerate(statements):
        end_line = statement.end_lineno or statement.lineno
        unsigned: dict[str, object] = {
            "member_index": index,
            "kind": type(statement).__name__,
            "line": statement.lineno,
            "end_line": end_line,
            "physical_span": end_line - statement.lineno + 1,
            "source_fingerprint": source.fingerprint,
        }
        members.append(
            {
                "signature": _candidate_signature(unsigned),
                **unsigned,
            },
        )

    return members


def _transfer_context(
    source: ParsedSource,
    owner: ast.AST,
    field: str,
    suite: list[ast.stmt],
    index: int,
    value: int,
) -> TransferContext:
    """
    Derive predecessor, preparation, and Z facts for one transfer.
    """

    node = suite[index]
    predecessor = suite[index - 1] if index else None
    blank_count = _blank_count_before(source, node)

    phase_start = index

    if predecessor is not None:
        phase_start = index - 1

        while phase_start > 0 and not _blank_between(
            source,
            suite[phase_start - 1],
            suite[phase_start],
        ):
            phase_start -= 1

    preceding_phase = tuple(suite[phase_start:index])
    preceding_phase_lines = (
        (
            (predecessor.end_lineno or predecessor.lineno)
            - preceding_phase[0].lineno
            + 1
        )
        if predecessor is not None and preceding_phase
        else 0
    )
    predecessor_span = (
        (predecessor.end_lineno or predecessor.lineno) - predecessor.lineno + 1
        if predecessor is not None
        else 0
    )
    substantive_preceding_phase = bool(
        predecessor is not None
        and (
            len(preceding_phase) > 1
            or _needs_control_boundary(predecessor)
            or _is_direct_side_effect(predecessor)
            or predecessor_span > 1
        ),
    )
    preparation = preceding_phase if blank_count == 0 else ()
    preparation_lines = preceding_phase_lines if preparation else 0
    z_applies = field == "body" and isinstance(
        owner,
        (
            ast.If,
            ast.For,
            ast.AsyncFor,
            ast.While,
            ast.match_case,
        ),
    )
    z_exceeded = bool(
        z_applies and blank_count == 0 and preparation_lines > value,
    )

    return TransferContext(
        predecessor=predecessor,
        blank_count=blank_count,
        preceding_phase=preceding_phase,
        preceding_phase_lines=preceding_phase_lines,
        substantive_preceding_phase=substantive_preceding_phase,
        preparation=preparation,
        preparation_lines=preparation_lines,
        z_applies=z_applies,
        z_exceeded=z_exceeded,
    )


def _transfer_record(
    source: ParsedSource,
    owner: ast.AST,
    field: str,
    node: ast.stmt,
    context: TransferContext,
    value: int,
    configuration: dict[str, object],
) -> dict[str, object]:
    """
    Build one complete transfer record with its semantic disposition.
    """

    predecessor = context.predecessor
    unsigned: dict[str, object] = {
        "path": source.path,
        "enclosing_object": _enclosing_object(node, source.parents),
        "suite_owner": type(owner).__name__,
        "suite_field": field,
        "line": node.lineno,
        "end_line": node.end_lineno or node.lineno,
        "column": node.col_offset + 1,
        "transfer_kind": _transfer_kind(node),
        "predecessor_kind": (
            type(predecessor).__name__ if predecessor is not None else None
        ),
        "predecessor_line": (
            predecessor.lineno if predecessor is not None else None
        ),
        "predecessor_end_line": (
            predecessor.end_lineno or predecessor.lineno
            if predecessor is not None
            else None
        ),
        "physical_preparation": bool(context.preparation),
        "preparation_statement_count": len(context.preparation),
        "preparation_physical_line_count": context.preparation_lines,
        "preceding_phase_statement_count": len(context.preceding_phase),
        "preceding_phase_physical_line_count": (context.preceding_phase_lines),
        "preceding_phase_members": _phase_members(
            source,
            context.preceding_phase,
        ),
        "substantive_preceding_phase": (context.substantive_preceding_phase),
        "blank_line_count_before": context.blank_count,
        "blank_line_state": (
            "no_predecessor"
            if predecessor is None
            else "separated"
            if context.blank_count
            else "attached"
        ),
        "z_applies": context.z_applies,
        "z_limit": value if context.z_applies else None,
        "z_exceeded": context.z_exceeded,
        "source_fingerprint": source.fingerprint,
        "producer_version": VERSION,
        "disposition": "unresolved",
        "disposition_reason": (
            "No fresh explicit semantic decision covers this transfer."
        ),
    }
    review_key = _review_key("transfers", unsigned)

    return {
        "signature": _candidate_signature(unsigned),
        "review_key": review_key,
        "review_group": None,
        **unsigned,
    }


def _transfer_candidate(
    source: ParsedSource,
    node: ast.stmt,
    context: TransferContext,
    value: int,
    configuration: dict[str, object],
) -> Candidate:
    """
    Build one exceeded-Z candidate from an already derived context.
    """

    return _candidate(
        source,
        node,
        construct="terminal_transfer",
        reason=(
            f"{context.preparation_lines} attached preparatory physical "
            f"lines exceed Z={value}"
        ),
        parameter_name="Z",
        parameter_value=value,
        metric="compact_transfer_preparation",
        unit="physical_lines",
        measured_value=context.preparation_lines,
        logical_span=len(context.preparation),
        item_count=len(context.preparation),
        false_positive_risk=(
            "preparation may be one indivisible guard or generator operation"
        ),
        configuration=configuration,
    )


def _transfer_inventory(
    sources: Sequence[ParsedSource],
    value: int,
    configuration: dict[str, object],
) -> tuple[list[dict[str, object]], int, list[Candidate]]:
    """
    Inventory every suite transfer and derive bounded Z candidates.
    """

    records: list[dict[str, object]] = []
    candidates: list[Candidate] = []

    for source in sources:
        for owner, field, suite in _statement_suites(source.tree):
            for index, node in enumerate(suite):
                if _transfer_kind(node) is None:
                    continue

                context = _transfer_context(
                    source,
                    owner,
                    field,
                    suite,
                    index,
                    value,
                )
                records.append(
                    _transfer_record(
                        source,
                        owner,
                        field,
                        node,
                        context,
                        value,
                        configuration,
                    ),
                )

                if context.z_exceeded:
                    candidates.append(
                        _transfer_candidate(
                            source,
                            node,
                            context,
                            value,
                            configuration,
                        ),
                    )

    return (
        sorted(
            records,
            key=lambda item: (
                str(item["path"]),
                int(item["line"]),
                str(item["transfer_kind"]),
            ),
        ),
        len(candidates),
        candidates,
    )


def _is_docstring_statement(node: ast.stmt) -> bool:
    """
    Return whether one statement is a literal docstring expression.
    """

    return (
        isinstance(node, ast.Expr)
        and isinstance(node.value, ast.Constant)
        and isinstance(node.value.value, str)
    )


def _statement_subject(node: ast.stmt) -> str:
    """
    Return one compact inspectable subject for a density-region member.

    Parameters
    ----------
    node : ast.stmt
        Statement whose semantic subject is summarized.

    Returns
    -------
    subject : str
        Bounded source-like label for the statement.
    """

    if isinstance(node, ast.Assign):
        subject = ", ".join(ast.unparse(target) for target in node.targets)
    elif isinstance(node, (ast.AnnAssign, ast.AugAssign)):
        subject = ast.unparse(node.target)
    elif isinstance(node, ast.Expr) and isinstance(node.value, ast.Call):
        subject = ast.unparse(node.value.func)
    elif isinstance(node, ast.If):
        subject = f"if {ast.unparse(node.test)}"
    elif isinstance(node, (ast.For, ast.AsyncFor)):
        subject = f"for {ast.unparse(node.target)}"
    elif isinstance(node, ast.While):
        subject = f"while {ast.unparse(node.test)}"
    elif isinstance(node, ast.With):
        subject = "with " + ", ".join(
            ast.unparse(item.context_expr) for item in node.items
        )
    elif isinstance(node, ast.Try):
        subject = "try"
    elif isinstance(node, ast.Return):
        subject = "return"
    elif isinstance(node, ast.Raise):
        subject = "raise"
    else:
        subject = type(node).__name__

    return subject if len(subject) <= 120 else subject[:117] + "..."


def _suite_regions(
    source: ParsedSource,
    owner: ast.AST,
    field: str,
    suite: list[ast.stmt],
) -> list[dict[str, object]]:
    """
    Return blank-line-delimited statement regions for one suite.
    """

    runs: list[list[ast.stmt]] = []
    current: list[ast.stmt] = []

    for position, statement in enumerate(suite):
        if isinstance(statement, DEFINITION) or (
            position == 0 and _is_docstring_statement(statement)
        ):
            if current:
                runs.append(current)
                current = []

            continue

        if current and _blank_between(source, current[-1], statement):
            runs.append(current)
            current = []

        current.append(statement)

    if current:
        runs.append(current)

    records: list[dict[str, object]] = []

    for run in runs:
        line = run[0].lineno
        end_line = run[-1].end_lineno or run[-1].lineno
        construct_counts = Counter(type(item).__name__ for item in run)
        records.append(
            {
                "path": source.path,
                "enclosing_object": _enclosing_object(
                    run[0],
                    source.parents,
                ),
                "suite_owner": type(owner).__name__,
                "suite_field": field,
                "line": line,
                "end_line": end_line,
                "physical_span": end_line - line + 1,
                "logical_line_count": len(run),
                "statement_count": len(run),
                "nesting": max(_nesting(item, source.parents) for item in run),
                "construct_counts": dict(sorted(construct_counts.items())),
                "statement_subjects": [
                    _statement_subject(statement) for statement in run
                ],
                "comment_line_count": sum(
                    source.lines[number - 1].lstrip().startswith("#")
                    for number in range(line, end_line + 1)
                ),
                "source_fingerprint": source.fingerprint,
            },
        )

    return records


def _density_review_group(
    region: dict[str, object],
) -> tuple[str, str]:
    """
    Return one explicit semantic-review group and its shared rationale.
    """

    statement_count = int(region["statement_count"])
    owner = str(region["enclosing_object"])
    suite_owner = str(region["suite_owner"])
    path = str(region["path"])
    constructs = set(region["construct_counts"])

    if statement_count == 1:
        return (
            "single_statement_container",
            "The region contains one syntactic statement. Multiline call and "
            "assignment shape is reviewed independently by X and Y, while "
            "nested compound suites are separate density records.",
        )

    if owner == "parse_args" or owner.endswith(".parse_args"):
        return (
            "ordered_cli_schema",
            "The members are consecutive parser or subparser declarations "
            "that form one ordered CLI schema.",
        )

    if constructs <= {"Import", "ImportFrom"}:
        return (
            "module_import_block",
            "The members are one canonical static import block; paragraph "
            "breaks inside it would misstate import grouping.",
        )

    static_declarations = constructs <= {
        "Assign",
        "AnnAssign",
        "Expr",
    }

    if suite_owner == "Module" and static_declarations:
        return (
            "module_declaration_block",
            "The members declare one adjacent module-level constant, schema, "
            "or registration family.",
        )

    if suite_owner == "ClassDef" and static_declarations:
        return (
            "class_declaration_block",
            "The members declare one adjacent class field or test-fixture "
            "family.",
        )

    if path.startswith("tests/"):
        return (
            "test_phase",
            "The members form one blank-line-delimited test setup, action, or "
            "assertion phase; exact statement subjects expose membership.",
        )

    if suite_owner in {
        "For",
        "If",
        "Try",
        "While",
        "With",
    }:
        return (
            "control_flow_phase",
            "The members form one guard, iteration, exception, or managed-"
            "resource phase inside a separately reviewed callable.",
        )

    return (
        "callable_phase",
        "The members form one blank-line-delimited callable phase; exact "
        "statement subjects expose the retained construction.",
    )


def _paragraph_density(
    sources: Sequence[ParsedSource],
) -> dict[str, object]:
    """
    Measure syntactic blank-line regions independently from X and Y.

    Parameters
    ----------
    sources : Sequence[ParsedSource]
        Parsed maintained sources whose statement suites are reviewed.

    Returns
    -------
    evidence : dict[str, object]
        Density distributions, largest regions, and inspectable review groups.
    """

    regions: list[dict[str, object]] = []

    for source in sources:
        for owner, field, suite in _statement_suites(source.tree):
            regions.extend(_suite_regions(source, owner, field, suite))

    physical_spans = [int(region["physical_span"]) for region in regions]
    statement_counts = [int(region["statement_count"]) for region in regions]

    physical_histogram = Counter(physical_spans)
    statement_histogram = Counter(statement_counts)

    physical_p90 = (
        sorted(physical_spans)[
            max(
                0,
                int(len(physical_spans) * 0.9) - 1,
            )
        ]
        if physical_spans
        else 0
    )

    statement_p90 = (
        sorted(statement_counts)[max(0, int(len(statement_counts) * 0.9) - 1)]
        if statement_counts
        else 0
    )

    reviewed_regions = []
    review_groups: defaultdict[str, dict[str, object]] = defaultdict(dict)

    for region in regions:
        if (
            int(region["physical_span"]) < physical_p90
            and int(region["statement_count"]) < statement_p90
        ):
            continue

        record = dict(region)
        group_id, group_rationale = _density_review_group(region)
        record["suggested_review_group"] = group_id
        record["suggested_group_rationale"] = group_rationale
        record["disposition"] = "unresolved"
        record["disposition_reason"] = (
            "No fresh explicit semantic decision covers this density region."
        )
        record["review_key"] = _review_key("density_regions", record)
        record["review_group"] = None
        reviewed_regions.append(record)

    return {
        "region_count": len(regions),
        "physical_span_histogram": {
            str(key): value
            for key, value in sorted(physical_histogram.items())
        },
        "statement_count_histogram": {
            str(key): value
            for key, value in sorted(statement_histogram.items())
        },
        "physical_span_median": (
            statistics.median(physical_spans) if physical_spans else 0
        ),
        "physical_span_p90": physical_p90,
        "physical_span_maximum": max(physical_spans, default=0),
        "statement_count_median": (
            statistics.median(statement_counts) if statement_counts else 0
        ),
        "statement_count_p90": statement_p90,
        "statement_count_maximum": max(statement_counts, default=0),
        "largest_regions": sorted(
            regions,
            key=lambda item: (
                -int(item["physical_span"]),
                -int(item["statement_count"]),
                str(item["path"]),
                int(item["line"]),
            ),
        )[:30],
        "reviewed_region_count": len(reviewed_regions),
        "reviewed_disposition_counts": {
            "resolved": sum(
                region["disposition"] in RESOLVED_DISPOSITIONS
                for region in reviewed_regions
            ),
            "unresolved": sum(
                region["disposition"] == "unresolved"
                for region in reviewed_regions
            ),
        },
        "review_groups": [
            {
                **group,
                "member_count": len(group["members"]),
            }
            for _, group in sorted(review_groups.items())
        ],
        "reviewed_regions": sorted(
            reviewed_regions,
            key=lambda item: (
                str(item["path"]),
                int(item["line"]),
                int(item["end_line"]),
            ),
        ),
    }


def _definition_role(
    node: ast.AST,
    parents: dict[ast.AST, ast.AST],
) -> str:
    """
    Return one bounded identifier role.
    """

    if isinstance(node, ast.ClassDef):
        return "class"

    if isinstance(node, CALLABLE):
        parent = parents.get(node)
        if node.name.startswith("test_"):
            return "test_function"

        if isinstance(parent, ast.ClassDef):
            return "method"

        if node.name.startswith("_"):
            return "private_function"

        return "function"

    if isinstance(node, ast.arg):
        return "parameter"

    return "local"


def _scope_name(
    node: ast.AST,
    parents: dict[ast.AST, ast.AST],
) -> str:
    """
    Return the nearest callable, class, or module scope.
    """

    current: ast.AST | None = node

    while current is not None:
        if isinstance(current, (ast.ClassDef, *CALLABLE)):
            return _object_name(current, parents)

        current = parents.get(current)

    return "<module>"


def _stored_name_role(
    node: ast.Name,
    parents: dict[ast.AST, ast.AST],
) -> str:
    """
    Return a bounded role for one stored name.
    """

    current = parents.get(node)

    while current is not None:
        if isinstance(current, CALLABLE):
            return "local"

        if isinstance(current, ast.ClassDef):
            if node.id.isupper():
                return "class_constant"

            return "class_attribute"

        if isinstance(current, ast.Module):
            return (
                "module_constant" if node.id.isupper() else "module_variable"
            )

        current = parents.get(current)

    return "local"


def _abbreviation_candidate(segment: str) -> bool:
    """
    Return whether one name segment is plausible abbreviation evidence.
    """

    lower = segment.lower()
    if lower in COMMON_SHORT_WORDS:
        return False

    return lower in KNOWN_ABBREVIATIONS


def _reference_counts(
    sources: Sequence[ParsedSource],
) -> tuple[
    Counter[str],
    Counter[str],
    defaultdict[str, set[str]],
    Counter[tuple[str, str, str]],
]:
    """
    Return bounded identifier and literal reference surfaces.
    """

    uses: Counter[str] = Counter()
    string_mentions: Counter[str] = Counter()
    paths_by_name: defaultdict[str, set[str]] = defaultdict(set)
    scoped_uses: Counter[tuple[str, str, str]] = Counter()

    for source in sources:
        for node in ast.walk(source.tree):
            if isinstance(node, ast.Name):
                uses[node.id] += 1
                paths_by_name[node.id].add(source.path)
                scoped_uses[
                    (
                        source.path,
                        _scope_name(node, source.parents),
                        node.id,
                    )
                ] += 1
            elif isinstance(node, ast.Attribute):
                uses[node.attr] += 1
                paths_by_name[node.attr].add(source.path)

            if isinstance(node, ast.Constant) and isinstance(node.value, str):
                for word in re.findall(r"[A-Za-z_][A-Za-z0-9_]*", node.value):
                    string_mentions[word] += 1

    return uses, string_mentions, paths_by_name, scoped_uses


def _name_records(
    sources: Sequence[ParsedSource],
    reference_sources: Sequence[ParsedSource],
    thresholds: dict[str, list[int]],
) -> dict[str, object]:
    """
    Produce role-aware name, abbreviation, and reference evidence.

    Parameters
    ----------
    sources : Sequence[ParsedSource]
        Parsed maintained sources included in the evidence cohort.
    reference_sources : Sequence[ParsedSource]
        Sources searched for naming references.
    thresholds : dict[str, list[int]]
        Role-specific name-length review thresholds.

    Returns
    -------
    evidence : dict[str, object]
        Role-aware naming records and reference distributions.
    """

    (
        uses,
        string_mentions,
        paths_by_name,
        scoped_uses,
    ) = _reference_counts(reference_sources)
    records: list[dict[str, object]] = []

    for source in sources:
        for node in ast.walk(source.tree):
            named: list[tuple[str, str, int, int, str]] = []

            if isinstance(node, (ast.ClassDef, *CALLABLE)):
                named.append(
                    (
                        node.name,
                        _definition_role(node, source.parents),
                        node.lineno,
                        node.col_offset + 1,
                        _enclosing_object(node, source.parents),
                    ),
                )

            if isinstance(node, CALLABLE):
                for argument in (
                    *node.args.posonlyargs,
                    *node.args.args,
                    *node.args.kwonlyargs,
                    node.args.vararg,
                    node.args.kwarg,
                ):
                    if argument is None:
                        continue

                    named.append(
                        (
                            argument.arg,
                            "parameter",
                            argument.lineno,
                            argument.col_offset + 1,
                            _object_name(node, source.parents),
                        ),
                    )

            if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Store):
                named.append(
                    (
                        node.id,
                        _stored_name_role(node, source.parents),
                        node.lineno,
                        node.col_offset + 1,
                        _enclosing_object(node, source.parents),
                    ),
                )

            if isinstance(node, ast.Attribute) and isinstance(
                node.ctx,
                ast.Store,
            ):
                named.append(
                    (
                        node.attr,
                        "attribute",
                        node.lineno,
                        max(
                            1,
                            (node.end_col_offset or 1) - len(node.attr) + 1,
                        ),
                        _enclosing_object(node, source.parents),
                    ),
                )

            for name, role, line, column, enclosing in named:
                segments = [
                    segment
                    for segment in name.strip("_").split("_")
                    if segment
                ]
                abbreviations = [
                    segment
                    for segment in segments
                    if _abbreviation_candidate(segment)
                ]
                records.append(
                    {
                        "path": source.path,
                        "enclosing_object": enclosing,
                        "role": role,
                        "visibility": (
                            "private"
                            if name.startswith("_")
                            and not name.startswith("__")
                            else "public"
                        ),
                        "line": line,
                        "column": column,
                        "spelling": name,
                        "character_count": len(name),
                        "segment_count": len(segments),
                        "abbreviation_tokens": abbreviations,
                    },
                )
                paths_by_name[name].add(source.path)

    deduplicated = {
        (
            record["path"],
            record["line"],
            record["column"],
            record["role"],
            record["spelling"],
        ): record
        for record in records
    }

    records = sorted(
        deduplicated.values(),
        key=lambda item: (
            str(item["path"]),
            int(item["line"]),
            int(item["column"]),
            str(item["role"]),
        ),
    )
    role_lengths: defaultdict[str, list[int]] = defaultdict(list)
    abbreviations: Counter[str] = Counter()

    for record in records:
        name = str(record["spelling"])
        role_lengths[str(record["role"])].append(len(name))
        abbreviations.update(record["abbreviation_tokens"])
        record["use_count"] = uses[name]
        record["scope_use_count"] = scoped_uses[
            (
                str(record["path"]),
                str(record["enclosing_object"]),
                name,
            )
        ]
        record["reference_path_count"] = len(paths_by_name[name])
        record["string_reference_count"] = string_mentions[name]

        role = str(record["role"])

        if role in {"local", "parameter"} and string_mentions[name] == 0:
            rename_status = "candidate_local_scope"
        elif (
            str(record["visibility"]) == "private"
            and len(paths_by_name[name]) == 1
            and string_mentions[name] == 0
        ):
            rename_status = "candidate_single_path"
        else:
            rename_status = "reference_review_required"

        record["direct_rename_status"] = rename_status

        if role in {"function", "method", "private_function"}:
            family = "production_callable"
        elif role == "test_function":
            family = "test_function"
        elif role == "class":
            family = "type"
        elif role in {
            "attribute",
            "class_attribute",
            "local",
            "module_variable",
            "parameter",
        }:
            family = "local_or_attribute"
        else:
            family = ""

        exceeded = (
            [value for value in thresholds[family] if len(name) > value]
            if family
            else []
        )
        record["length_thresholds_exceeded"] = exceeded
        semantic_candidate = bool(
            record["abbreviation_tokens"] or exceeded,
        )
        record["semantic_candidate"] = semantic_candidate

        if not semantic_candidate:
            record["semantic_disposition"] = "not_candidate"
            record["disposition_reason"] = (
                "No configured naming length or abbreviation review threshold "
                "applies."
            )
        else:
            record["semantic_disposition"] = "unresolved"
            record["disposition_reason"] = (
                "No fresh explicit semantic decision covers this naming "
                "candidate."
            )
            record["review_key"] = _review_key(
                "naming_candidates",
                record,
            )
            record["review_group"] = None

    distributions = {
        role: {
            "count": len(lengths),
            "minimum": min(lengths),
            "median": statistics.median(lengths),
            "maximum": max(lengths),
            "histogram": {
                str(length): count
                for length, count in sorted(Counter(lengths).items())
            },
        }
        for role, lengths in sorted(role_lengths.items())
    }
    threshold_counts: dict[str, dict[str, int]] = {}

    for family, values in sorted(thresholds.items()):
        if family == "production_callable":
            roles = {"function", "method", "private_function"}
        elif family == "test_function":
            roles = {"test_function"}
        elif family == "local_or_attribute":
            roles = {
                "attribute",
                "class_attribute",
                "local",
                "module_variable",
                "parameter",
            }
        elif family == "type":
            roles = {"class"}
        else:
            roles = set()

        family_lengths = [
            int(record["character_count"])
            for record in records
            if record["role"] in roles
        ]
        threshold_counts[family] = {
            str(value): sum(length > value for length in family_lengths)
            for value in values
        }

    semantic_records = [
        record for record in records if record["semantic_candidate"]
    ]

    return {
        "record_count": len(records),
        "length_distributions": distributions,
        "adjacent_threshold_counts": threshold_counts,
        "abbreviation_candidate_distribution": dict(
            sorted(
                abbreviations.items(),
                key=lambda item: (-item[1], item[0]),
            ),
        ),
        "reference_scope": "all maintained Python source",
        "semantic_candidate_count": len(semantic_records),
        "semantic_disposition_counts": {
            "resolved": sum(
                record["semantic_disposition"] in RESOLVED_DISPOSITIONS
                for record in semantic_records
            ),
            "unresolved": sum(
                record["semantic_disposition"] == "unresolved"
                for record in semantic_records
            ),
        },
        "records": records,
    }


def _docstring_role(
    node: ast.AST,
    parents: dict[ast.AST, ast.AST],
    path: str,
) -> str:
    """
    Return one object role for documentation applicability review.

    Parameters
    ----------
    node : ast.AST
        Module, class, or callable whose documentation role is classified.
    parents : dict[ast.AST, ast.AST]
        Parent map used to recover lexical and class ownership.
    path : str
        Repository-relative source path used for test and CLI roles.

    Returns
    -------
    role : str
        Role-aware documentation category for the object.
    """

    if isinstance(node, ast.Module):
        return "module"

    if isinstance(node, ast.ClassDef):
        return (
            "test_class"
            if node.name.startswith("Test") or path.startswith("tests/")
            else "class"
        )

    if isinstance(node, CALLABLE):
        current = parents.get(node)
        callable_ancestor = False
        test_class_ancestor = False

        while current is not None:
            if isinstance(current, CALLABLE):
                callable_ancestor = True

            if isinstance(current, ast.ClassDef) and (
                current.name.startswith("Test") or path.startswith("tests/")
            ):
                test_class_ancestor = True

            current = parents.get(current)

        if callable_ancestor:
            return "nested_helper"

        parent = parents.get(node)
        decorators = {
            decorator.id
            if isinstance(decorator, ast.Name)
            else decorator.attr
            if isinstance(decorator, ast.Attribute)
            else ""
            for decorator in node.decorator_list
        }

        if "fixture" in decorators or node.name.endswith("_fixture"):
            return "fixture"

        if node.name.startswith("test_"):
            return (
                "test_method"
                if isinstance(parent, ast.ClassDef)
                else "test_function"
            )

        if test_class_ancestor or path.startswith("tests/"):
            return "test_helper"

        if isinstance(parent, ast.ClassDef):
            if path.startswith("dev/"):
                return "audit_method"

            if "/cli/" in path:
                return "cli_method"

            if parent.name.startswith("_"):
                return "private_method"

            return "method"

        if node.name.startswith("_"):
            return "private_helper"

        if path.startswith("dev/"):
            return "audit_helper"

        if "/cli/" in path:
            return (
                "cli_boundary"
                if node.name in {"main", "parse_args"}
                else "cli_helper"
            )

        return "reusable_function"

    return "unknown"


def _owned_nodes(node: ast.AST) -> list[ast.AST]:
    """
    Return descendants owned by one object, excluding nested definitions.
    """

    descendants: list[ast.AST] = []
    stack = list(ast.iter_child_nodes(node))

    while stack:
        current = stack.pop()
        if isinstance(current, DEFINITION):
            continue

        descendants.append(current)
        stack.extend(ast.iter_child_nodes(current))

    return descendants


def _applicable_docstring_sections(node: ast.AST) -> list[str]:
    """
    Return syntax-supported NumPy sections applicable to one object.
    """

    if not isinstance(node, CALLABLE):
        return []

    positional = (*node.args.posonlyargs, *node.args.args)
    parameters = [
        argument
        for argument in (
            *positional,
            *node.args.kwonlyargs,
            node.args.vararg,
            node.args.kwarg,
        )
        if argument is not None and argument.arg not in {"self", "cls"}
    ]
    owned = _owned_nodes(node)
    sections: list[str] = []

    if parameters:
        sections.append("Parameters")

    if any(isinstance(item, (ast.Yield, ast.YieldFrom)) for item in owned):
        sections.append("Yields")
    elif any(
        isinstance(item, ast.Return) and item.value is not None
        for item in owned
    ):
        sections.append("Returns")

    if any(isinstance(item, ast.Raise) for item in owned):
        sections.append("Raises")

    if any(
        isinstance(item, ast.Call)
        and isinstance(item.func, ast.Attribute)
        and isinstance(item.func.value, ast.Name)
        and item.func.value.id == "warnings"
        and item.func.attr == "warn"
        for item in owned
    ):
        sections.append("Warns")

    return sections


def _docstring_complexity(node: ast.AST) -> dict[str, int | bool]:
    """
    Return bounded syntax facts used for semantic documentation review.
    """

    owned = _owned_nodes(node)
    statement_count = sum(isinstance(item, ast.stmt) for item in owned)
    branch_count = sum(
        isinstance(
            item,
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
        for item in owned
    )
    physical_span = (
        (getattr(node, "end_lineno", getattr(node, "lineno", 1)) or 1)
        - getattr(node, "lineno", 1)
        + 1
    )

    return {
        "owned_statement_count": statement_count,
        "branch_count": branch_count,
        "physical_span": physical_span,
        "behaviorally_complex": (statement_count >= 80 or branch_count >= 20),
    }


def _docstring_inventory(
    sources: Sequence[ParsedSource],
) -> dict[str, object]:
    """
    Inventory applicability and summary-only/expanded source shapes.

    Parameters
    ----------
    sources : Sequence[ParsedSource]
        Parsed maintained sources and their reference surfaces.

    Returns
    -------
    evidence : dict[str, object]
        Role-aware forms, applicable sections, behavioral facts, and review
        keys for every documentable object.
    """

    uses, _, paths_by_name, _ = _reference_counts(sources)
    records: list[dict[str, object]] = []

    for source in sources:
        for node in ast.walk(source.tree):
            if not isinstance(node, (ast.Module, ast.ClassDef, *CALLABLE)):
                continue

            value = ast.get_docstring(node, clean=True)
            nonblank = (
                [line for line in value.splitlines() if line.strip()]
                if value is not None
                else []
            )
            sections = SECTION.findall(value or "")

            role = _docstring_role(node, source.parents, source.path)
            applicable_sections = _applicable_docstring_sections(node)
            missing_sections = [
                section
                for section in applicable_sections
                if section not in sections
            ]

            complexity = _docstring_complexity(node)
            owned = _owned_nodes(node)
            object_name = _object_name(node, source.parents)
            definition_name = getattr(node, "name", "")
            signature = (
                ast.unparse(node.args)
                if isinstance(node, CALLABLE)
                else object_name
            )

            failure_types = sorted(
                {
                    ast.unparse(item.exc)
                    for item in owned
                    if isinstance(item, ast.Raise) and item.exc is not None
                },
            )
            statement_subjects = [
                type(item).__name__
                for item in getattr(node, "body", [])
                if not (
                    isinstance(item, ast.Expr)
                    and isinstance(item.value, ast.Constant)
                    and isinstance(item.value.value, str)
                )
            ]

            form = (
                "missing"
                if value is None
                else "summary_only"
                if len(nonblank) == 1
                else "expanded"
            )

            record: dict[str, object] = {
                "path": source.path,
                "enclosing_object": object_name,
                "role": role,
                "line": getattr(node, "lineno", 1),
                "present": value is not None,
                "form": form,
                "summary": nonblank[0].strip() if nonblank else None,
                "signature": signature,
                "sections": sections,
                "applicable_sections": applicable_sections,
                "missing_applicable_sections": missing_sections,
                **complexity,
                "direct_statement_subjects": statement_subjects,
                "return_value_count": sum(
                    isinstance(item, ast.Return) and item.value is not None
                    for item in owned
                ),
                "yield_count": sum(
                    isinstance(item, (ast.Yield, ast.YieldFrom))
                    for item in owned
                ),
                "failure_types": failure_types,
                "direct_side_effect_count": sum(
                    isinstance(item, ast.stmt) and _is_direct_side_effect(item)
                    for item in owned
                ),
                "maintained_name_use_count": uses[definition_name],
                "maintained_reference_paths": sorted(
                    paths_by_name[definition_name],
                ),
                "has_examples": "Examples" in sections,
                "source_fingerprint": source.fingerprint,
                "suggested_review_group": f"{role}_{form}",
                "review_group": None,
                "applicability_disposition": "unresolved",
                "content_disposition": "unresolved",
                "disposition_reason": (
                    "No fresh explicit semantic decisions cover this "
                    "docstring's applicability and content."
                ),
            }

            record["review_keys"] = {
                owner: _review_key(owner, record)
                for owner in (
                    "docstring_applicability",
                    "docstring_content",
                )
            }
            records.append(record)

    role_counts: defaultdict[str, Counter[str]] = defaultdict(Counter)

    for record in records:
        role_counts[str(record["role"])][str(record["form"])] += 1

    sorted_records = sorted(
        records,
        key=lambda item: (
            str(item["path"]),
            int(item["line"]),
            str(item["enclosing_object"]),
        ),
    )

    group_members: defaultdict[
        str,
        list[dict[str, object]],
    ] = defaultdict(list)

    for record in sorted_records:
        group_members[str(record["review_group"])].append(
            {
                "path": record["path"],
                "line": record["line"],
                "enclosing_object": record["enclosing_object"],
            },
        )

    disposition_counts = {
        owner: {
            "resolved": sum(
                record[owner] in RESOLVED_DISPOSITIONS
                for record in sorted_records
            ),
            "unresolved": sum(
                record[owner] == "unresolved" for record in sorted_records
            ),
        }
        for owner in (
            "applicability_disposition",
            "content_disposition",
        )
    }

    return {
        "records": sorted_records,
        "role_form_counts": {
            role: dict(sorted(counts.items()))
            for role, counts in sorted(role_counts.items())
        },
        "review_groups": [
            {
                "group_id": group,
                "member_count": len(members),
                "members": members,
            }
            for group, members in sorted(group_members.items())
        ],
        "disposition_counts": disposition_counts,
    }


def _annotation_parameters(
    node: ast.FunctionDef | ast.AsyncFunctionDef,
    parent: ast.AST | None,
) -> list[ast.arg]:
    """
    Return callable parameters excluding an established receiver.
    """

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
        class_receiver = "classmethod" in decorators and first.arg == "cls"
        instance_receiver = (
            "staticmethod" not in decorators
            and "classmethod" not in decorators
            and first.arg == "self"
        )

        if class_receiver or instance_receiver:
            receiver = first

    return [
        argument
        for argument in (
            *node.args.posonlyargs,
            *node.args.args,
            *node.args.kwonlyargs,
            node.args.vararg,
            node.args.kwarg,
        )
        if argument is not None and argument is not receiver
    ]


def _annotation_record(
    source: ParsedSource,
    node: ast.FunctionDef | ast.AsyncFunctionDef,
    uses: Counter[str],
    paths_by_name: dict[str, set[str]],
) -> dict[str, object]:
    """
    Build one annotation-review record from runtime-relevant syntax facts.
    """

    parameters = _annotation_parameters(node, source.parents.get(node))
    missing = [
        argument.arg for argument in parameters if argument.annotation is None
    ]

    owned = _owned_nodes(node)
    return_values = [
        item
        for item in owned
        if isinstance(item, ast.Return) and item.value is not None
    ]
    bare_returns = [
        item
        for item in owned
        if isinstance(item, ast.Return) and item.value is None
    ]
    yields = [
        item for item in owned if isinstance(item, (ast.Yield, ast.YieldFrom))
    ]

    record: dict[str, object] = {
        "path": source.path,
        "enclosing_object": _object_name(node, source.parents),
        "line": node.lineno,
        "parameter_count": len(parameters),
        "missing_parameters": missing,
        "parameter_annotations": [
            {
                "name": argument.arg,
                "annotation": (
                    ast.unparse(argument.annotation)
                    if argument.annotation is not None
                    else None
                ),
            }
            for argument in parameters
        ],
        "return_present": node.returns is not None,
        "return_annotation": (
            ast.unparse(node.returns) if node.returns is not None else None
        ),
        "return_value_count": len(return_values),
        "bare_return_count": len(bare_returns),
        "yield_count": len(yields),
        "raise_count": sum(isinstance(item, ast.Raise) for item in owned),
        "callable_kind": (
            "async" if isinstance(node, ast.AsyncFunctionDef) else "sync"
        ),
        "maintained_name_use_count": uses[node.name],
        "maintained_reference_path_count": len(paths_by_name[node.name]),
        "source_fingerprint": source.fingerprint,
        "correctness_disposition": "unresolved",
        "disposition_reason": (
            "No fresh explicit semantic decision covers this callable's "
            "annotations."
        ),
    }

    record["review_key"] = _review_key("annotation_correctness", record)
    record["review_group"] = None

    return record


def _annotation_inventory(
    sources: Sequence[ParsedSource],
) -> dict[str, object]:
    """
    Inventory annotation presence without judging correctness.
    """

    uses, _, paths_by_name, _ = _reference_counts(sources)
    records: list[dict[str, object]] = []

    for source in sources:
        for node in ast.walk(source.tree):
            if not isinstance(node, CALLABLE):
                continue

            records.append(
                _annotation_record(
                    source,
                    node,
                    uses,
                    paths_by_name,
                ),
            )

    sorted_records = sorted(
        records,
        key=lambda item: (
            str(item["path"]),
            int(item["line"]),
        ),
    )

    return {
        "function_count": len(records),
        "missing_parameter_function_count": sum(
            bool(record["missing_parameters"]) for record in records
        ),
        "missing_return_count": sum(
            not bool(record["return_present"]) for record in records
        ),
        "correctness_counts": {
            "resolved": sum(
                record["correctness_disposition"] in RESOLVED_DISPOSITIONS
                for record in sorted_records
            ),
            "unresolved": sum(
                record["correctness_disposition"] == "unresolved"
                for record in sorted_records
            ),
        },
        "records": sorted_records,
    }


@dataclasses.dataclass(frozen=True)
class RepeatedThresholds:
    """
    Store the configured X and Y threshold families.
    """

    x_call_spans: tuple[int, ...]
    y_assignment_spans: tuple[int, ...]
    y_run_lengths: tuple[int, ...]


def _repeated_thresholds(
    configuration: dict[str, object],
) -> RepeatedThresholds:
    """
    Return integer X and Y thresholds from validated configuration.
    """

    return RepeatedThresholds(
        x_call_spans=tuple(
            int(value)
            for value in configuration["x_multiline_call_span_thresholds"]
        ),
        y_assignment_spans=tuple(
            int(value)
            for value in configuration[
                "y_multiline_assignment_span_thresholds"
            ]
        ),
        y_run_lengths=tuple(
            int(value)
            for value in configuration[
                "y_simple_assignment_run_statement_thresholds"
            ]
        ),
    )


def _repeated_statement_kind(node: ast.stmt | None) -> str:
    """
    Classify one repeated assignment or call statement.
    """

    if isinstance(node, ASSIGNMENT):
        return "assignment_run"

    if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call):
        return "call_run"

    return ""


def _repeated_member_record(
    source: ParsedSource,
    owner: ast.AST,
    field: str,
    run: list[ast.stmt],
    run_kind: str,
    member: ast.stmt,
    member_index: int,
    thresholds: RepeatedThresholds,
) -> dict[str, object]:
    """
    Build one explicit repeated-run member record.
    """

    span = (member.end_lineno or member.lineno) - member.lineno + 1

    if isinstance(member, ASSIGNMENT):
        targets = (
            [ast.unparse(target) for target in member.targets]
            if isinstance(member, ast.Assign)
            else [ast.unparse(member.target)]
        )
        callee: str | None = None
    else:
        targets = []
        callee = ast.unparse(member.value.func)

    simple_assignment_run = run_kind == "assignment_run" and all(
        item.lineno == item.end_lineno for item in run
    )
    unsigned: dict[str, object] = {
        "path": source.path,
        "enclosing_object": _enclosing_object(
            member,
            source.parents,
        ),
        "suite_owner": type(owner).__name__,
        "suite_field": field,
        "run_kind": run_kind,
        "member_index": member_index,
        "line": member.lineno,
        "end_line": member.end_lineno or member.lineno,
        "physical_span": span,
        "source_form": "simple" if span == 1 else "multiline",
        "callee": callee,
        "targets": targets,
        "blank_line_count_before": _blank_count_before(source, member),
        "x_thresholds_exceeded": (
            [value for value in thresholds.x_call_spans if span > value]
            if callee is not None
            else []
        ),
        "y_assignment_span_thresholds_exceeded": (
            [value for value in thresholds.y_assignment_spans if span > value]
            if targets
            else []
        ),
        "y_simple_run_thresholds_exceeded": (
            [
                value
                for value in thresholds.y_run_lengths
                if simple_assignment_run and len(run) > value
            ]
            if targets
            else []
        ),
        "disposition": "unresolved",
        "disposition_reason": (
            "No fresh explicit semantic decision covers this repeated-run "
            "member."
        ),
    }

    return {
        "signature": _candidate_signature(unsigned),
        "review_key": _review_key("repeated_run_members", unsigned),
        "review_group": None,
        **unsigned,
    }


def _repeated_run_record(
    source: ParsedSource,
    owner: ast.AST,
    field: str,
    run: list[ast.stmt],
    run_kind: str,
    thresholds: RepeatedThresholds,
) -> dict[str, object] | None:
    """
    Build one repeated-run record and its exact member inventory.
    """

    if len(run) < 2:
        return None

    members = [
        _repeated_member_record(
            source,
            owner,
            field,
            run,
            run_kind,
            member,
            member_index,
            thresholds,
        )
        for member_index, member in enumerate(run, 1)
    ]
    unsigned: dict[str, object] = {
        "path": source.path,
        "enclosing_object": _enclosing_object(
            run[0],
            source.parents,
        ),
        "suite_owner": type(owner).__name__,
        "suite_field": field,
        "line": run[0].lineno,
        "end_line": run[-1].end_lineno or run[-1].lineno,
        "kind": run_kind,
        "statement_count": len(run),
        "source_fingerprint": source.fingerprint,
        "producer_version": VERSION,
        "semantic_topic_disposition": "unresolved",
        "disposition_reason": (
            "No fresh explicit semantic decision covers this repeated "
            "statement run."
        ),
        "members": members,
    }

    return {
        "signature": _candidate_signature(unsigned),
        "review_key": _review_key("repeated_runs", unsigned),
        "review_group": None,
        **unsigned,
    }


def _suite_repeated_runs(
    source: ParsedSource,
    owner: ast.AST,
    field: str,
    suite: list[ast.stmt],
    thresholds: RepeatedThresholds,
) -> list[dict[str, object]]:
    """
    Return every repeated assignment or call run in one statement suite.
    """

    records: list[dict[str, object]] = []
    run: list[ast.stmt] = []
    run_kind = ""

    for node in [*suite, None]:
        kind = _repeated_statement_kind(node)
        separated = bool(
            kind and run and _blank_between(source, run[-1], node),
        )

        if kind and not separated and (not run or kind == run_kind):
            run.append(node)
            run_kind = kind

            continue

        record = _repeated_run_record(
            source,
            owner,
            field,
            run,
            run_kind,
            thresholds,
        )

        if record is not None:
            records.append(record)

        run = [node] if kind and node is not None else []
        run_kind = kind

    return records


def _repeated_statement_inventory(
    sources: Sequence[ParsedSource],
    configuration: dict[str, object],
) -> list[dict[str, object]]:
    """
    Inventory complete blank-line-delimited repeated-statement runs.

    Parameters
    ----------
    sources : Sequence[ParsedSource]
        Parsed maintained sources and statement suites.
    configuration : dict[str, object]
        Current candidate thresholds and review configuration.

    Returns
    -------
    records : list[dict[str, object]]
        Complete run and member records with stable signatures and boundaries.
    """

    thresholds = _repeated_thresholds(configuration)
    records: list[dict[str, object]] = []

    for source in sources:
        for owner, field, suite in _statement_suites(source.tree):
            records.extend(
                _suite_repeated_runs(
                    source,
                    owner,
                    field,
                    suite,
                    thresholds,
                ),
            )

    return sorted(
        records,
        key=lambda item: (
            str(item["path"]),
            int(item["line"]),
        ),
    )


def _token_inventories(
    sources: Sequence[ParsedSource],
) -> dict[str, object]:
    """
    Inventory comment and ordinary-string source forms.

    Parameters
    ----------
    sources : Sequence[ParsedSource]
        Parsed maintained sources whose tokens are classified.

    Returns
    -------
    inventory : dict[str, object]
        Stable comment-marker and string-quote distributions.
    """

    quote_counts: Counter[str] = Counter()
    comment_counts: Counter[str] = Counter()

    for source in sources:
        doc_positions = {
            (
                node.body[0].lineno,
                node.body[0].col_offset,
            )
            for node in ast.walk(source.tree)
            if isinstance(node, (ast.Module, ast.ClassDef, *CALLABLE))
            and node.body
            and _is_docstring_statement(node.body[0])
        }

        analysis = analyze_text(source.text, source.path)
        comment_counts.update(analysis.facts["comments"])
        reader = io.StringIO(source.text).readline

        try:
            tokens = tokenize.generate_tokens(reader)
            fstring_depth = 0

            for token in tokens:
                if token.type == getattr(tokenize, "FSTRING_START", -1):
                    fstring_depth += 1
                elif token.type == getattr(tokenize, "FSTRING_END", -1):
                    fstring_depth = max(0, fstring_depth - 1)
                elif (
                    token.type == tokenize.STRING
                    and token.start not in doc_positions
                ):
                    prefix = re.match(r"(?i)^[rubf]*", token.string)
                    start = len(prefix.group(0)) if prefix else 0

                    if token.string[start:].startswith('"""'):
                        quote_counts["triple_double"] += 1
                    elif token.string[start:].startswith("'''"):
                        quote_counts["triple_single"] += 1
                    elif token.string[start:].startswith('"'):
                        quote_counts["double"] += 1
                    elif token.string[start:].startswith("'"):
                        key = (
                            "single_in_fstring_expression"
                            if fstring_depth
                            else "single"
                        )
                        quote_counts[key] += 1
        except (IndentationError, tokenize.TokenError):
            continue

    return {
        "ordinary_strings": dict(sorted(quote_counts.items())),
        "comments": dict(sorted(comment_counts.items())),
    }


def _raw_width_inventory(
    sources: Sequence[ParsedSource],
) -> dict[str, object]:
    """
    Review every physical source line wider than the advisory limit.
    """

    fixture_reasons = {
        "tests/unit/dev_audit/test_help_heredoc_reflow.py": (
            "The literal fixture intentionally places list prose near the "
            "ordinary source-width boundary so reflow detection sees the "
            "physical source form."
        ),
        "tests/unit/dev_audit/test_shell_help_pilot.py": (
            "The literal fixture intentionally preserves an exact overlength "
            "command or prose candidate for the line-length pilot."
        ),
    }
    records: list[dict[str, object]] = []

    for source in sources:
        for line_number, line in enumerate(source.lines, 1):
            if len(line) <= 79:
                continue

            stripped = line.strip()
            test_match = re.match(r"def (test_[A-Za-z0-9_]+)\(", stripped)

            record: dict[str, object] = {
                "path": source.path,
                "line": line_number,
                "length": len(line),
                "source": line,
                "source_fingerprint": source.fingerprint,
                "suggested_indivisible_test_identifier": (
                    test_match.group(1) if test_match else None
                ),
                "suggested_fixture_reason": fixture_reasons.get(
                    source.path,
                ),
                "disposition": "unresolved",
                "disposition_reason": (
                    "No fresh explicit semantic decision covers this raw "
                    "width candidate."
                ),
            }

            record["review_key"] = _review_key(
                "raw_width_candidates",
                record,
            )
            record["review_group"] = None
            records.append(record)

    return {
        "candidate_count": len(records),
        "disposition_counts": {
            "resolved": sum(
                record["disposition"] in RESOLVED_DISPOSITIONS
                for record in records
            ),
            "unresolved": sum(
                record["disposition"] == "unresolved" for record in records
            ),
        },
        "records": records,
    }


def _exit_boundary_inventory(
    sources: Sequence[ParsedSource],
) -> dict[str, object]:
    """
    Review explicit process exits against reusable-helper ownership.
    """

    records: list[dict[str, object]] = []

    for source in sources:
        for node in ast.walk(source.tree):
            kind: str | None = None

            if isinstance(node, ast.Call):
                if (
                    isinstance(node.func, ast.Attribute)
                    and isinstance(node.func.value, ast.Name)
                    and node.func.value.id == "sys"
                    and node.func.attr == "exit"
                ):
                    kind = "sys.exit"
            elif isinstance(node, ast.Raise):
                exception = node.exc

                if isinstance(exception, ast.Call):
                    exception = exception.func

                if (
                    isinstance(exception, ast.Name)
                    and exception.id == "SystemExit"
                ):
                    kind = "raise SystemExit"

            if kind is None:
                continue

            owner = _enclosing_object(node, source.parents)
            cli_boundary = owner in {"<module>", "main", "parse_args"}
            record: dict[str, object] = {
                "path": source.path,
                "line": node.lineno,
                "enclosing_object": owner,
                "exit_kind": kind,
                "structural_cli_boundary": cli_boundary,
                "source_fingerprint": source.fingerprint,
                "disposition": "unresolved",
                "disposition_reason": (
                    "No fresh explicit semantic decision covers this "
                    "process-exit boundary."
                ),
            }
            record["review_key"] = _review_key(
                "reusable_error_boundaries",
                record,
            )
            record["review_group"] = None
            records.append(record)

    return {
        "record_count": len(records),
        "disposition_counts": {
            "resolved": sum(
                record["disposition"] in RESOLVED_DISPOSITIONS
                for record in records
            ),
            "unresolved": sum(
                record["disposition"] == "unresolved" for record in records
            ),
        },
        "records": sorted(
            records,
            key=lambda item: (
                str(item["path"]),
                int(item["line"]),
            ),
        ),
    }


def _suppression_inventory(
    root: Path,
    sources: Sequence[ParsedSource],
) -> dict[str, object]:
    """
    Inventory inline Ruff and type-checker suppressions plus E402 policy.
    """

    records = []
    marker = re.compile(r"#\s*(?:noqa\b|type:\s*ignore\b)")

    for source in sources:
        for line_number, line in enumerate(source.lines, 1):
            if marker.search(line):
                record: dict[str, object] = {
                    "path": source.path,
                    "line": line_number,
                    "source": line,
                    "source_fingerprint": source.fingerprint,
                    "disposition": "unresolved",
                    "disposition_reason": (
                        "No fresh explicit semantic decision covers this "
                        "inline suppression."
                    ),
                }
                record["review_key"] = _review_key(
                    "suppressions_and_e402",
                    record,
                )
                record["review_group"] = None
                records.append(record)

    pyproject_path = root / "pyproject.toml"
    pyproject = (
        pyproject_path.read_text(encoding="utf-8")
        if pyproject_path.is_file()
        else ""
    )

    blanket = bool(
        re.search(
            r"(?ms)per-file-ignores.*?E402",
            pyproject,
        ),
    )
    e402_record: dict[str, object] = {
        "path": "pyproject.toml",
        "blanket_e402_ignore_present": blanket,
        "source_fingerprint": source_fingerprint(pyproject),
        "disposition": "unresolved",
        "disposition_reason": (
            "No fresh explicit semantic decision covers the E402 "
            "configuration and import architecture."
        ),
    }

    e402_record["review_key"] = _review_key(
        "suppressions_and_e402",
        e402_record,
    )
    e402_record["review_group"] = None

    return {
        "inline_record_count": len(records),
        "inline_records": records,
        "blanket_e402_ignore_present": blanket,
        "e402_record": e402_record,
        "e402_disposition": "unresolved",
        "e402_disposition_reason": e402_record["disposition_reason"],
    }


def _source_header_inventory(
    sources: Sequence[ParsedSource],
) -> dict[str, object]:
    """
    Review shebang and opening-header applicability for every source.
    """

    corrected_paths = {
        "dev/audit/generate_prompts.py",
        "dev/audit/parse_findings.py",
        "tests/unit/dev_audit/test_dependency_closure.py",
        "tests/unit/dev_audit/test_mvp.py",
        "tests/unit/dev_audit/test_shell_help_pilot.py",
    }
    adopted_script_sources = {
        "dev/audit/generate_prompts.py",
        "dev/audit/parse_findings.py",
    }
    records: list[dict[str, object]] = []

    for source in sources:
        lines = source.lines
        has_shebang = bool(
            lines and lines[0] == "#!/usr/bin/env python3",
        )
        direct_guard = any(
            isinstance(node, ast.If)
            and isinstance(node.test, ast.Compare)
            and isinstance(node.test.left, ast.Name)
            and node.test.left.id == "__name__"
            for node in source.tree.body
        )

        adopted_script_source = source.path in adopted_script_sources
        applicable = has_shebang or direct_guard or adopted_script_source

        script_row = f"# Script: {Path(source.path).name}"
        header_window = lines[:24]
        structure_valid = bool(
            has_shebang
            and len(lines) > 12
            and lines[1] == "# -*- coding: utf-8 -*-"
            and script_row in header_window
            and any(
                line.startswith("# Copyright ")
                and line.endswith(" by Kris Alavattam")
                for line in header_window
            )
            and "# Email: kalavattam@gmail.com" in header_window
            and "# Distributed under the MIT license." in header_window,
        )

        record: dict[str, object] = {
            "path": source.path,
            "source_fingerprint": source.fingerprint,
            "has_approved_shebang": has_shebang,
            "has_direct_execution_guard": direct_guard,
            "adopted_script_source": adopted_script_source,
            "applicable": applicable,
            "structure_valid": structure_valid,
            "previously_corrected": source.path in corrected_paths,
            "disposition": "unresolved",
            "disposition_reason": (
                "No fresh explicit semantic decision covers this source "
                "header."
            ),
        }

        record["review_key"] = _review_key("source_headers", record)
        record["review_group"] = None
        records.append(record)

    return {
        "record_count": len(records),
        "applicable_count": sum(
            bool(record["applicable"]) for record in records
        ),
        "import_only_exception_count": sum(
            record["disposition"] == "exception" for record in records
        ),
        "disposition_counts": _disposition_counts(records),
        "records": records,
    }


def _compute_input_floor_decision(
    sources: Sequence[ParsedSource],
) -> dict[str, object]:
    """
    Record the completed naming and reusable-error boundary decision.
    """

    target = next(
        (
            source
            for source in sources
            if source.path.endswith("/cli/compute_input_floor.py")
        ),
        None,
    )
    if target is None:
        return {"applicable": False}

    old_name = "count_aln_bam"
    new_name = "_count_alignment_records"

    old_references = [
        (source.path, node.lineno)
        for source in sources
        for node in ast.walk(source.tree)
        if (
            (isinstance(node, ast.Name) and node.id == old_name)
            or (isinstance(node, ast.Attribute) and node.attr == old_name)
        )
    ]
    new_references = [
        (source.path, node.lineno)
        for source in sources
        for node in ast.walk(source.tree)
        if (
            (isinstance(node, ast.Name) and node.id == new_name)
            or (isinstance(node, ast.Attribute) and node.attr == new_name)
        )
    ]

    record: dict[str, object] = {
        "applicable": True,
        "source_fingerprint": target.fingerprint,
        "old_spelling": old_name,
        "new_spelling": new_name,
        "visibility_decision": "private reusable helper",
        "compatibility_alias": False,
        "alias_removal_condition": (
            "not applicable; no external contract found"
        ),
        "active_old_reference_count": len(old_references),
        "active_new_reference_count": len(new_references),
        "reference_paths": sorted(
            {path for path, _ in new_references},
        ),
        "exception_types": [
            "AlignmentReadError",
            "FlagParseError",
            "InputFloorValidationError",
        ],
        "disposition": "unresolved",
        "disposition_reason": (
            "No fresh explicit semantic decision covers this completed rename "
            "and error-boundary change."
        ),
    }

    record["review_key"] = _review_key(
        "completed_rename_decisions",
        record,
    )
    record["review_group"] = None

    return record


def _configured_type_checker(root: Path) -> dict[str, object]:
    """
    Record configured static-type-checker evidence without selecting one.
    """

    candidates = (
        "mypy.ini",
        ".mypy.ini",
        "pyrightconfig.json",
        "pyproject.toml",
        "setup.cfg",
    )
    configured: list[str] = []
    markers = ("[tool.mypy]", "[mypy", "[tool.pyright]", "pyright")

    for name in candidates:
        path = root / name
        if not path.is_file():
            continue

        text = path.read_text(encoding="utf-8")

        if any(marker in text for marker in markers):
            configured.append(name)

    return {
        "configured": bool(configured),
        "configuration_paths": configured,
        "selection_authorized": False,
        "installation_attempted": False,
        "recommendation": (
            "Use the configured checker as additional non-authoritative "
            "evidence; retain annotation correctness as semantic review."
            if configured
            else "No configured checker is available for this pilot; enforce "
            "presence, retain correctness as semantic review, and revisit "
            "tool selection separately."
        ),
    }


def _checker_inventory(
    sources: Sequence[ParsedSource],
) -> dict[str, object]:
    """
    Aggregate deterministic checker interactions without enforcing them.
    """

    by_rule: Counter[str] = Counter()
    facts: list[dict[str, object]] = []

    for source in sources:
        analysis = analyze_text(source.text, source.path)
        by_rule.update(finding.rule_id for finding in analysis.findings)
        facts.append(analysis.facts)

    return {
        "checker_version": CHECKER_VERSION,
        "finding_count": sum(by_rule.values()),
        "findings_by_rule": dict(sorted(by_rule.items())),
        "facts": facts,
    }


def _semantic_record_surfaces(
    inventories: dict[str, object],
    candidates: list[dict[str, object]],
) -> dict[str, list[tuple[dict[str, object], str]]]:
    """
    Return every semantic owner record and its disposition field.
    """

    docstrings = inventories["docstrings"]["records"]
    annotations = inventories["annotations"]["records"]
    naming = inventories["naming"]["records"]

    repeated = inventories["repeated_statements"]
    transfers = inventories["transfers"]
    density = inventories["paragraph_density"]["reviewed_regions"]

    raw_width = inventories["raw_width"]["records"]
    exits = inventories["reusable_error_boundaries"]["records"]
    suppressions = inventories["suppressions"]

    headers = inventories["source_headers"]["records"]
    floor = inventories["compute_input_floor_decision"]

    return {
        "docstring_applicability": [
            (record, "applicability_disposition") for record in docstrings
        ],
        "docstring_content": [
            (record, "content_disposition") for record in docstrings
        ],
        "annotation_correctness": [
            (record, "correctness_disposition") for record in annotations
        ],
        "naming_candidates": [
            (record, "semantic_disposition")
            for record in naming
            if record["semantic_candidate"]
        ],
        "completed_rename_decisions": (
            [(floor, "disposition")] if floor.get("applicable") else []
        ),
        "x_y_z_candidates": [
            (candidate, "disposition") for candidate in candidates
        ],
        "repeated_runs": [
            (record, "semantic_topic_disposition") for record in repeated
        ],
        "repeated_run_members": [
            (member, "disposition")
            for record in repeated
            for member in record["members"]
        ],
        "transfers": [(record, "disposition") for record in transfers],
        "density_regions": [(record, "disposition") for record in density],
        "raw_width_candidates": [
            (record, "disposition") for record in raw_width
        ],
        "reusable_error_boundaries": [
            (record, "disposition") for record in exits
        ],
        "suppressions_and_e402": [
            *[
                (record, "disposition")
                for record in suppressions["inline_records"]
            ],
            (suppressions["e402_record"], "disposition"),
        ],
        "source_headers": [(record, "disposition") for record in headers],
    }


def _record_review_key(
    owner: str,
    record: dict[str, object],
) -> str:
    """
    Return the owner-specific review key stored on one semantic record.
    """

    if owner in {"docstring_applicability", "docstring_content"}:
        keys = record.get("review_keys")
        if isinstance(keys, dict) and isinstance(keys.get(owner), str):
            return str(keys[owner])

    key = record.get("review_key")
    if not isinstance(key, str):
        raise ValueError(f"semantic record for {owner} lacks a review key")

    return key


def _load_review_manifest(
    root: Path,
    configuration: dict[str, object],
) -> tuple[dict[str, object], str | None, str | None]:
    """
    Load an optional explicit semantic-decision manifest.

    Parameters
    ----------
    root : Path
        Repository root that confines the configured manifest.
    configuration : dict[str, object]
        Evidence configuration containing the optional manifest path.

    Returns
    -------
    manifest, path, fingerprint : tuple[
        dict[str, object], str | None, str | None
    ]
        Parsed manifest, repository-relative path, and content fingerprint.

    Raises
    ------
    ValueError
        If the path escapes the repository or the manifest schema is invalid.
    """

    relative = configuration.get("semantic_review_manifest")
    if relative is None:
        return {"schema_version": 1, "decision_groups": []}, None, None

    if not isinstance(relative, str) or not relative:
        raise ValueError("semantic_review_manifest must be a nonempty path")

    path = (root / relative).resolve()

    try:
        path.relative_to(root)
    except ValueError as exc:
        raise ValueError(
            "semantic_review_manifest must remain under the repository root",
        ) from exc

    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError("semantic review manifest must be an object")

    if value.get("schema_version") != 1:
        raise ValueError("semantic review manifest schema_version must be 1")

    groups = value.get("decision_groups")
    if not isinstance(groups, list):
        raise ValueError("semantic review decision_groups must be a list")

    protocol = value.get("review_protocol")

    if not isinstance(protocol, dict):
        raise ValueError("semantic review manifest needs review_protocol")

    if protocol.get("status") not in {"provisional", "reviewed"}:
        raise ValueError(
            "semantic review status must be provisional or reviewed",
        )

    reviewed_source = value.get("reviewed_source_fingerprint")

    if reviewed_source is not None and (
        not isinstance(reviewed_source, str)
        or not re.fullmatch(r"sha256:[0-9a-f]{64}", reviewed_source)
    ):
        raise ValueError(
            "reviewed_source_fingerprint must be null or a SHA-256 identity",
        )

    return (
        value,
        relative,
        source_fingerprint(
            path.read_text(encoding="utf-8"),
        ),
    )


def _validated_decision_group(
    value: object,
    seen_group_ids: set[str],
) -> ExplicitDecisionGroup:
    """
    Validate and return one explicit semantic decision group.

    Parameters
    ----------
    value : object
        Candidate manifest value to validate.
    seen_group_ids : set[str]
        Group identifiers already accepted from this manifest.

    Returns
    -------
    group : ExplicitDecisionGroup
        Validated decision fields and exact membership.

    Raises
    ------
    ValueError
        If identity, ownership, disposition, rationale, or membership fails.
    """

    if not isinstance(value, dict):
        raise ValueError("each semantic decision group must be an object")

    group_id = value.get("group_id")
    if not isinstance(group_id, str) or not group_id:
        raise ValueError("semantic decision group_id must be nonempty")

    if group_id in seen_group_ids:
        raise ValueError(f"duplicate semantic decision group: {group_id}")

    seen_group_ids.add(group_id)

    owner = value.get("owner")
    disposition = value.get("disposition")
    rationale = value.get("rationale")
    if owner not in SEMANTIC_OWNERS:
        raise ValueError(f"unknown semantic decision owner: {owner}")

    if disposition not in RESOLVED_DISPOSITIONS:
        raise ValueError(
            f"semantic decision {group_id} has unresolved disposition",
        )

    if not isinstance(rationale, str) or not rationale.strip():
        raise ValueError(
            f"semantic decision {group_id} lacks a rationale",
        )

    reviewed_facts = value.get("reviewed_facts")

    if not isinstance(reviewed_facts, dict) or not reviewed_facts:
        raise ValueError(
            f"semantic decision {group_id} lacks inspectable reviewed_facts",
        )

    members = value.get("members")
    if (
        not isinstance(members, list)
        or not members
        or not all(isinstance(member, str) for member in members)
        or len(members) != len(set(members))
    ):
        raise ValueError(
            f"semantic decision {group_id} needs unique members",
        )

    fingerprint = value.get("membership_fingerprint")

    if fingerprint != _membership_fingerprint(members):
        raise ValueError(
            f"semantic decision {group_id} membership hash is stale",
        )

    return ExplicitDecisionGroup(
        group_id=group_id,
        owner=str(owner),
        disposition=str(disposition),
        rationale=rationale,
        reviewed_facts=reviewed_facts,
        members=tuple(members),
        membership_fingerprint=str(fingerprint),
    )


def _apply_semantic_decisions(
    inventories: dict[str, object],
    candidates: list[dict[str, object]],
    manifest: dict[str, object],
    *,
    manifest_path: str | None,
    manifest_fingerprint: str | None,
    source_fingerprint: str,
) -> dict[str, object]:
    """
    Apply only fresh, explicit, exact-membership semantic decisions.

    Parameters
    ----------
    inventories : dict[str, object]
        Current semantic record inventories.
    candidates : list[dict[str, object]]
        Current X, Y, and Z candidate records.
    manifest : dict[str, object]
        Explicit exact-membership decision groups.
    manifest_path : str | None
        Repository-relative manifest path, when configured.
    manifest_fingerprint : str | None
        Current manifest content fingerprint.
    source_fingerprint : str
        Current selected-cohort source identity.

    Returns
    -------
    evidence : dict[str, object]
        Per-owner decision coverage, staleness, and completion facts.

    Raises
    ------
    ValueError
        If current keys duplicate or a decision repeats owner membership.
    """

    surfaces = _semantic_record_surfaces(inventories, candidates)
    by_owner: dict[
        str,
        dict[str, tuple[dict[str, object], str]],
    ] = {}

    for owner, records in surfaces.items():
        keyed: dict[str, tuple[dict[str, object], str]] = {}

        for record, field in records:
            key = _record_review_key(owner, record)
            if key in keyed:
                raise ValueError(
                    f"duplicate current semantic review key for {owner}: "
                    f"{key}",
                )

            keyed[key] = (record, field)

        by_owner[owner] = keyed

    seen_group_ids: set[str] = set()
    assigned: defaultdict[str, set[str]] = defaultdict(set)
    applied_groups: list[dict[str, object]] = []
    stale_groups: list[dict[str, object]] = []
    groups = manifest.get("decision_groups", [])

    assert isinstance(groups, list)

    for value in groups:
        group = _validated_decision_group(value, seen_group_ids)

        duplicates = assigned[group.owner].intersection(group.members)
        if duplicates:
            raise ValueError(
                f"semantic decision {group.group_id} duplicates owner members",
            )

        assigned[group.owner].update(group.members)

        current = by_owner[group.owner]
        missing = sorted(set(group.members) - set(current))

        if missing:
            stale_groups.append(
                {
                    "group_id": group.group_id,
                    "owner": group.owner,
                    "missing_member_count": len(missing),
                    "missing_members": missing,
                },
            )

            continue

        for member in group.members:
            record, field = current[member]
            record[field] = group.disposition

            decisions = record.setdefault("review_decisions", {})

            assert isinstance(decisions, dict)

            decisions[group.owner] = {
                "group_id": group.group_id,
                "disposition": group.disposition,
                "rationale": group.rationale,
                "reviewed_facts": group.reviewed_facts,
                "membership_fingerprint": group.membership_fingerprint,
            }

            if group.owner in {
                "docstring_applicability",
                "docstring_content",
            }:
                review_groups = record.get("review_groups")

                if not isinstance(review_groups, dict):
                    review_groups = {}

                record["review_groups"] = {
                    **review_groups,
                    group.owner: group.group_id,
                }
            else:
                record["review_group"] = group.group_id

            record["disposition_reason"] = group.rationale

        applied_groups.append(
            {
                "group_id": group.group_id,
                "owner": group.owner,
                "disposition": group.disposition,
                "rationale": group.rationale,
                "reviewed_facts": group.reviewed_facts,
                "member_count": len(group.members),
                "members": list(group.members),
                "membership_fingerprint": group.membership_fingerprint,
            },
        )

    current_record_count = sum(len(records) for records in by_owner.values())
    explicit_decision_count = sum(
        len(group["members"]) for group in applied_groups
    )
    unresolved_record_count = sum(
        record[field] == "unresolved"
        for records in surfaces.values()
        for record, field in records
    )

    protocol = manifest.get("review_protocol")
    review_status = (
        protocol.get("status") if isinstance(protocol, dict) else None
    )
    reviewed_source_fingerprint = manifest.get(
        "reviewed_source_fingerprint",
    )
    source_fingerprint_matches = (
        reviewed_source_fingerprint == source_fingerprint
    )
    complete = (
        manifest_path is not None
        and review_status == "reviewed"
        and source_fingerprint_matches
        and not stale_groups
        and unresolved_record_count == 0
    )

    return {
        "manifest_path": manifest_path,
        "manifest_fingerprint": manifest_fingerprint,
        "review_status": review_status,
        "reviewed_source_fingerprint": reviewed_source_fingerprint,
        "source_fingerprint_matches": source_fingerprint_matches,
        "declared_group_count": len(groups),
        "applied_group_count": len(applied_groups),
        "stale_group_count": len(stale_groups),
        "current_record_count": current_record_count,
        "explicit_decision_count": explicit_decision_count,
        "unresolved_record_count": unresolved_record_count,
        "complete": complete,
        "applied_groups": applied_groups,
        "stale_groups": stale_groups,
        "hash_role": (
            "Record and manifest hashes only invalidate stale explicit "
            "decisions; they never create a disposition."
        ),
    }


def _refresh_semantic_inventory_counts(
    inventories: dict[str, object],
) -> None:
    """
    Refresh inventory-local summaries after explicit decisions are applied.
    """

    docstrings = inventories["docstrings"]
    docstrings["disposition_counts"] = {
        field: {
            "resolved": sum(
                record[field] in RESOLVED_DISPOSITIONS
                for record in docstrings["records"]
            ),
            "unresolved": sum(
                record[field] == "unresolved"
                for record in docstrings["records"]
            ),
        }
        for field in (
            "applicability_disposition",
            "content_disposition",
        )
    }
    docstring_groups: defaultdict[
        tuple[str | None, str | None],
        list[dict[str, object]],
    ] = defaultdict(list)

    for record in docstrings["records"]:
        decision_groups = record.get("review_groups", {})

        if not isinstance(decision_groups, dict):
            decision_groups = {}

        group_key = (
            decision_groups.get("docstring_applicability"),
            decision_groups.get("docstring_content"),
        )
        docstring_groups[group_key].append(
            {
                "path": record["path"],
                "line": record["line"],
                "enclosing_object": record["enclosing_object"],
                "role": record["role"],
                "form": record["form"],
            },
        )

    docstrings["review_groups"] = [
        {
            "group_id": (f"applicability={applicability};content={content}"),
            "applicability_group": applicability,
            "content_group": content,
            "member_count": len(members),
            "members": members,
        }
        for (applicability, content), members in sorted(
            docstring_groups.items(),
            key=lambda item: str(item[0]),
        )
    ]

    annotations = inventories["annotations"]
    annotations["correctness_counts"] = {
        "resolved": sum(
            record["correctness_disposition"] in RESOLVED_DISPOSITIONS
            for record in annotations["records"]
        ),
        "unresolved": sum(
            record["correctness_disposition"] == "unresolved"
            for record in annotations["records"]
        ),
    }

    naming = inventories["naming"]
    naming_records = [
        record for record in naming["records"] if record["semantic_candidate"]
    ]
    naming["semantic_disposition_counts"] = {
        "resolved": sum(
            record["semantic_disposition"] in RESOLVED_DISPOSITIONS
            for record in naming_records
        ),
        "unresolved": sum(
            record["semantic_disposition"] == "unresolved"
            for record in naming_records
        ),
    }

    density = inventories["paragraph_density"]
    density["reviewed_disposition_counts"] = {
        "resolved": sum(
            record["disposition"] in RESOLVED_DISPOSITIONS
            for record in density["reviewed_regions"]
        ),
        "unresolved": sum(
            record["disposition"] == "unresolved"
            for record in density["reviewed_regions"]
        ),
    }

    for name in (
        "raw_width",
        "reusable_error_boundaries",
        "source_headers",
    ):
        inventory = inventories[name]
        inventory["disposition_counts"] = _disposition_counts(
            inventory["records"],
        )

    headers = inventories["source_headers"]
    headers["import_only_exception_count"] = sum(
        record["disposition"] == "exception" for record in headers["records"]
    )
    headers["import_only_omission_count"] = sum(
        record["disposition"] == "omitted_by_role"
        for record in headers["records"]
    )


def _disposition_counts(
    records: Iterable[dict[str, object]],
    key: str = "disposition",
) -> dict[str, int]:
    """
    Count resolved and unresolved semantic decisions.
    """

    values = [str(record.get(key, "unresolved")) for record in records]
    retained = {
        "retained",
        "retained_coherent",
        "retained_compact_guard",
        "retained_indivisible",
    }

    return {
        "total": len(values),
        "resolved": sum(value in RESOLVED_DISPOSITIONS for value in values),
        "fixed": sum(value in {"changed", "refactored"} for value in values),
        "retained": sum(value in retained for value in values),
        "omitted_by_role": sum(value == "omitted_by_role" for value in values),
        "exception": sum(value == "exception" for value in values),
        "unresolved": sum(value == "unresolved" for value in values),
    }


def _semantic_owner_counts(
    inventories: dict[str, object],
    candidates: Sequence[dict[str, object]],
) -> dict[str, dict[str, int]]:
    """
    Return separate resolution counts for every semantic owner.
    """

    docstrings = inventories["docstrings"]
    annotations = inventories["annotations"]
    naming = inventories["naming"]

    repeated = inventories["repeated_statements"]
    transfers = inventories["transfers"]
    density = inventories["paragraph_density"]
    raw_width = inventories["raw_width"]

    exits = inventories["reusable_error_boundaries"]
    suppressions = inventories["suppressions"]
    headers = inventories["source_headers"]
    floor_decision = inventories["compute_input_floor_decision"]

    repeated_members = [
        member for record in repeated for member in record["members"]
    ]
    suppression_records = [
        *suppressions["inline_records"],
        suppressions["e402_record"],
    ]

    floor_records = (
        [floor_decision] if floor_decision.get("applicable") else []
    )

    return {
        "docstring_applicability": _disposition_counts(
            docstrings["records"],
            "applicability_disposition",
        ),
        "docstring_content": _disposition_counts(
            docstrings["records"],
            "content_disposition",
        ),
        "annotation_correctness": _disposition_counts(
            annotations["records"],
            "correctness_disposition",
        ),
        "naming_candidates": _disposition_counts(
            [
                record
                for record in naming["records"]
                if record["semantic_candidate"]
            ],
            "semantic_disposition",
        ),
        "completed_rename_decisions": _disposition_counts(floor_records),
        "x_y_z_candidates": _disposition_counts(candidates),
        "repeated_runs": _disposition_counts(
            repeated,
            "semantic_topic_disposition",
        ),
        "repeated_run_members": _disposition_counts(repeated_members),
        "transfers": _disposition_counts(transfers),
        "density_regions": _disposition_counts(
            density["reviewed_regions"],
        ),
        "raw_width_candidates": _disposition_counts(
            raw_width["records"],
        ),
        "reusable_error_boundaries": _disposition_counts(
            exits["records"],
        ),
        "suppressions_and_e402": _disposition_counts(
            suppression_records,
        ),
        "source_headers": _disposition_counts(headers["records"]),
    }


def produce(
    root: Path,
    paths: Iterable[str],
    configuration: dict[str, object],
) -> dict[str, object]:
    """
    Produce one stable no-write evidence payload.

    Parameters
    ----------
    root : Path
        Repository root used to resolve maintained source paths.
    paths : Iterable[str]
        Repository-relative Python paths in the selected cohort.
    configuration : dict[str, object]
        Validated thresholds and exact-hash semantic-review state.

    Returns
    -------
    evidence : dict[str, object]
        Deterministic inventories, candidate dispositions, and owner counts.
    """

    root = root.resolve()
    sources = _parse_sources(root, paths)

    cohort_fingerprint = _combined_fingerprint(sources)
    manifest, manifest_path, manifest_fingerprint = _load_review_manifest(
        root,
        configuration,
    )

    maintained_paths = [
        path.relative_to(root).as_posix()
        for path in maintained_python_paths(root)
    ]
    selected_paths = [source.path for source in sources]
    reference_sources = (
        sources
        if selected_paths == sorted(maintained_paths)
        else _parse_sources(root, maintained_paths)
    )

    x_span_thresholds = [
        int(value)
        for value in configuration["x_multiline_call_span_thresholds"]
    ]
    y_span_thresholds = [
        int(value)
        for value in configuration["y_multiline_assignment_span_thresholds"]
    ]
    y_run_thresholds = [
        int(value)
        for value in configuration[
            "y_simple_assignment_run_statement_thresholds"
        ]
    ]
    z_line_limit = int(
        configuration["z_compact_transfer_preparatory_line_limit"],
    )

    candidate_config = {
        "x_multiline_call_span_thresholds": x_span_thresholds,
        "y_multiline_assignment_span_thresholds": y_span_thresholds,
        "y_simple_assignment_run_statement_thresholds": y_run_thresholds,
        "z_compact_transfer_preparatory_line_limit": z_line_limit,
    }

    x_counts, x_candidates = _x_evidence(
        sources,
        x_span_thresholds,
        candidate_config,
    )
    y_counts, y_candidates = _y_evidence(
        sources,
        y_span_thresholds,
        candidate_config,
    )
    y_run_counts, y_run_candidates = _y_simple_run_evidence(
        sources,
        y_run_thresholds,
        candidate_config,
    )
    transfers, z_count, z_candidates = _transfer_inventory(
        sources,
        z_line_limit,
        candidate_config,
    )

    selected_x = int(configuration.get("selected_x_threshold", 6))
    selected_y = int(configuration.get("selected_y_threshold", 6))
    candidates = sorted(
        [
            *[
                candidate
                for candidate in x_candidates
                if candidate.parameter["value"] == selected_x
            ],
            *[
                candidate
                for candidate in y_candidates
                if candidate.parameter["value"] == selected_y
            ],
            *[
                candidate
                for candidate in y_run_candidates
                if candidate.parameter["value"] == selected_y
            ],
            *z_candidates,
        ],
        key=lambda item: (
            item.path,
            item.line,
            item.construct,
            str(item.parameter["name"]),
            int(item.parameter["value"]),
        ),
    )

    roots = Counter(
        path.split("/", 1)[0] for path in (s.path for s in sources)
    )

    inventories = {
        "maintained_files": {
            "total": len(sources),
            "by_root": dict(sorted(roots.items())),
            "hashes": {source.path: source.fingerprint for source in sources},
        },
        "x_distributions": x_counts,
        "y_distributions": {
            "multiline_assignment_span": y_counts,
            "simple_assignment_run_length": y_run_counts,
        },
        "z_distribution": {str(z_line_limit): z_count},
        "paragraph_density": _paragraph_density(sources),
        "repeated_statements": _repeated_statement_inventory(
            sources,
            candidate_config,
        ),
        "transfers": transfers,
        "naming": _name_records(
            sources,
            reference_sources,
            {
                str(role): [int(value) for value in values]
                for role, values in configuration[
                    "naming_length_thresholds"
                ].items()
            },
        ),
        "docstrings": _docstring_inventory(sources),
        "annotations": _annotation_inventory(sources),
        "raw_width": _raw_width_inventory(sources),
        "reusable_error_boundaries": _exit_boundary_inventory(sources),
        "suppressions": _suppression_inventory(root, sources),
        "source_headers": _source_header_inventory(sources),
        "compute_input_floor_decision": _compute_input_floor_decision(sources),
        "tokens": _token_inventories(sources),
        "checker_interactions": _checker_inventory(sources),
        "type_checker": _configured_type_checker(root),
    }

    candidate_records = [
        dataclasses.asdict(candidate) for candidate in candidates
    ]
    semantic_review = _apply_semantic_decisions(
        inventories,
        candidate_records,
        manifest,
        manifest_path=manifest_path,
        manifest_fingerprint=manifest_fingerprint,
        source_fingerprint=cohort_fingerprint,
    )
    _refresh_semantic_inventory_counts(inventories)
    suppressions = inventories["suppressions"]
    suppressions["e402_disposition"] = suppressions["e402_record"][
        "disposition"
    ]
    suppressions["e402_disposition_reason"] = suppressions["e402_record"][
        "disposition_reason"
    ]
    owner_counts = _semantic_owner_counts(
        inventories,
        candidate_records,
    )

    return {
        "schema_version": 2,
        "producer": "python_source_evidence",
        "producer_version": VERSION,
        "configuration": configuration,
        "paths": [source.path for source in sources],
        "source_fingerprint": cohort_fingerprint,
        "semantic_review": semantic_review,
        "inventories": inventories,
        "semantic_owner_counts": owner_counts,
        "candidates": candidate_records,
        "limitations": [
            "X, both distinct Y metrics, Z, blank-line-delimited source "
            "density, naming length, abbreviations, docstring applicability, "
            "and type correctness remain review evidence rather than "
            "violations.",
            "AST and token facts cannot decide readability, documentation "
            "usefulness, naming quality, or semantic independence.",
            "Dynamic, serialized, reflected, and generated references require "
            "human review before any rename.",
        ],
    }


def _load_config(path: Path) -> dict[str, object]:
    """
    Load and validate the bounded pilot configuration.

    Parameters
    ----------
    path : Path
        Repository-relative path associated with the source.

    Returns
    -------
    configuration : dict[str, object]
        Validated source-policy configuration.

    Raises
    ------
    ValueError
        If a supplied value violates the validated contract.
    """

    value = json.loads(path.read_text(encoding="utf-8"))
    required = {
        "schema_version",
        "pilot_id",
        "threshold_definitions",
        "pilot_paths",
        "x_multiline_call_span_thresholds",
        "y_multiline_assignment_span_thresholds",
        "y_simple_assignment_run_statement_thresholds",
        "z_compact_transfer_preparatory_line_limit",
        "line_length",
        "naming_length_thresholds",
        "semantic_review_manifest",
        "disposition",
        "semantic_candidates_are_nonblocking",
    }
    missing = required - set(value)
    if missing:
        raise ValueError(
            "pilot configuration missing: " + ", ".join(sorted(missing)),
        )

    if value["schema_version"] != 2:
        raise ValueError("pilot schema_version must be 2")

    if not isinstance(value["pilot_id"], str) or not value["pilot_id"]:
        raise ValueError("pilot_id must be a nonempty string")

    paths = value["pilot_paths"]
    if (
        not isinstance(paths, list)
        or not paths
        or not all(isinstance(item, str) and item for item in paths)
        or len(paths) != len(set(paths))
    ):
        raise ValueError("pilot_paths must be a nonempty unique string list")

    for name in (
        "x_multiline_call_span_thresholds",
        "y_multiline_assignment_span_thresholds",
        "y_simple_assignment_run_statement_thresholds",
    ):
        values = value[name]
        if (
            not isinstance(values, list)
            or not values
            or not all(isinstance(item, int) and item > 0 for item in values)
            or values != sorted(set(values))
        ):
            raise ValueError(f"{name} must be a sorted unique positive list")

    definitions = value["threshold_definitions"]
    expected_definition_names = {
        "X",
        "Y_assignment_span",
        "Y_simple_run_length",
        "Z",
    }
    definition_names_valid = (
        isinstance(definitions, dict)
        and set(definitions) == expected_definition_names
    )

    if not definition_names_valid:
        raise ValueError(
            "threshold_definitions must explain X, both Y metrics, and Z",
        )

    if not all(
        isinstance(definition, str) and definition
        for definition in definitions.values()
    ):
        raise ValueError("threshold definitions must be nonempty strings")

    z_limit = value["z_compact_transfer_preparatory_line_limit"]
    if not isinstance(z_limit, int) or z_limit < 0:
        raise ValueError(
            "z_compact_transfer_preparatory_line_limit must be a nonnegative "
            "integer",
        )

    if value["line_length"] != 79:
        raise ValueError("line_length must be the approved value 79")

    manifest = value["semantic_review_manifest"]
    if not isinstance(manifest, str) or not manifest:
        raise ValueError("semantic_review_manifest must be a nonempty path")

    if value.get("selected_x_threshold") != 6:
        raise ValueError("selected_x_threshold must be 6")

    if value.get("selected_y_threshold") != 6:
        raise ValueError("selected_y_threshold must be 6")

    thresholds = value["naming_length_thresholds"]
    expected_roles = {
        "production_callable",
        "test_function",
        "local_or_attribute",
        "type",
    }
    if not isinstance(thresholds, dict) or set(thresholds) != expected_roles:
        raise ValueError("naming_length_thresholds has invalid roles")

    for role, values in thresholds.items():
        if (
            not isinstance(values, list)
            or len(values) != 3
            or not all(isinstance(item, int) and item > 0 for item in values)
            or values != sorted(set(values))
        ):
            raise ValueError(
                f"naming threshold family '{role}' must contain three "
                "sorted unique positive values",
            )

    if value["disposition"] != "unresolved":
        raise ValueError("pilot disposition must remain unresolved")

    if value["semantic_candidates_are_nonblocking"] is not True:
        raise ValueError("semantic candidates must remain nonblocking")

    return value


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
        Parsed evidence-cohort, configuration, and output options.
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)

    selection = parser.add_mutually_exclusive_group()

    selection.add_argument("--all-maintained", action="store_true")
    selection.add_argument("--pilot", action="store_true")
    parser.add_argument("--output", type=Path)
    parser.add_argument("paths", nargs="*")

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """
    Write or print stable nonblocking Python source-policy evidence.

    Parameters
    ----------
    argv : list[str] | None
        Explicit arguments, or None to read the process arguments.

    Returns
    -------
    status : int
        Zero after evidence is written or rendered successfully.

    Raises
    ------
    ValueError
        If cohort selectors conflict or semantic evidence is malformed.
    """

    args = parse_args(argv)
    root = args.root.resolve()
    configuration = _load_config(args.config)
    if args.paths and (args.all_maintained or args.pilot):
        raise ValueError(
            "explicit paths cannot be combined with a cohort flag",
        )

    if args.paths:
        paths = args.paths
    elif args.all_maintained:
        paths = [
            path.relative_to(root).as_posix()
            for path in maintained_python_paths(root)
        ]
    elif args.pilot:
        paths = configuration["pilot_paths"]
    else:
        paths = configuration["pilot_paths"]

    payload = produce(root, paths, configuration)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"

    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
