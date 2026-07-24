#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: help_style.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


"""Check bounded structural conventions in changed shell help documentation."""

from __future__ import annotations

import argparse
import ast
import dataclasses
import re
import sys
from collections.abc import Iterable
from itertools import pairwise
from pathlib import Path

from help_heredoc_reflow import (
    SECTION_NAMES,
    UNDERLINE,
    Heredoc,
    audit_rendered_example_commands,
    changed_lines,
    extract_help_heredocs,
    indentation,
    run_git,
    shell_paths,
)
from shell_help_pilot import validate_command_registry

RULE_CANONICAL_HELP = "HELP.DOC.CANONICAL_HELP"
RULE_USAGE = "HELP.USAGE.STRUCTURE"
RULE_ENTRY_INDENT = "HELP.ENTRY.INDENT"
RULE_ENTRY_BLANK = "HELP.ENTRY.BLANK"
RULE_GLOBAL_TYPE = "HELP.GLOBAL.TYPE"
RULE_RUNTIME_INDENT = "HELP.RUNTIME.INDENT"
RULE_RUNTIME_STRUCTURE = "HELP.RUNTIME.STRUCTURE"
RULE_RUNTIME_PROGRAMS = "HELP.RUNTIME.REDUNDANT_PROGRAMS"
RULE_RUNTIME_CALLABLE = "HELP.RUNTIME.EXACT_CALLABLE"
RULE_RUNTIME_CARDINALITY = "HELP.RUNTIME.CARDINALITY"
RULE_SECTION_INDENT = "HELP.SECTION.INDENT"
RULE_NOTES_PSEUDO_INDENT = "HELP.NOTES.PSEUDO_HEADER_INDENT"
RULE_UNSUPPORTED_DEPENDENCY = "HELP.SECTION.UNSUPPORTED_DEPENDENCY"
RULE_SECTION_DUPLICATE = "HELP.SECTION.DUPLICATE"
RULE_SECTION_ORDER = "HELP.SECTION.ORDER"
RULE_SECTION_REQUIRED = "HELP.SECTION.REQUIRED"
RULE_LONG_OPTIONS = "HELP.DOC.LONG_OPTIONS"
RULE_PARAMETER_ALIAS_DUPLICATE = "HELP.PARAMETER.ALIAS_DUPLICATE"
RULE_GLOBAL_SUBGROUP = "HELP.GLOBAL.SUBGROUP"
RULE_GLOBAL_GROUP_SYNTAX = "HELP.GLOBAL.GROUP_SYNTAX"
RULE_GLOBAL_REPEATED_DESCRIPTION = "HELP.GLOBAL.REPEATED_DESCRIPTION"
RULE_GLOBAL_MAPPING_REVIEW = "HELP.GLOBAL.MAPPING_REVIEW"
RULE_PROSE_QUOTES = "HELP.PROSE.STRAIGHT_QUOTES"
RULE_PARSER_UNKNOWN = "SHELL.PARSER.UNKNOWN_ERROR"
RULE_PYTHON_PROSE_QUOTES = "PYTHON.DOCSTRING.STRAIGHT_QUOTES"

FUNCTION_SECTION_ORDER = (
    "Usage",
    "Parameters",
    "Expected globals",
    "Generated globals",
    "Returns",
    "Output",
    "Notes",
    "Examples",
)
FUNCTION_SECTION_RANK = {
    name: index for index, name in enumerate(FUNCTION_SECTION_ORDER)
}

HIDDEN_HELP = re.compile(r"(?<![-A-Za-z0-9_])-h\b|--hlp\b")
SHORT_OPTION = re.compile(r"(?<![-A-Za-z0-9_])-(?!-|\d)([A-Za-z0-9][A-Za-z0-9]*)")
COMMAND_NAME = re.compile(r"[A-Za-z0-9_./-]+")
PARAMETER_ENTRY = re.compile(
    r"^(?P<indent> *)(?P<head>(?:\d+\+?\s+\S+|-[^:]+?))\s+:\s+(?P<type>\S.*)$"
)
GLOBAL_ENTRY = re.compile(
    r"^(?P<indent> *)(?P<head>[A-Za-z_][A-Za-z0-9_]*(?:\[\])?)\s+:\s+(?P<type>\S.*)$"
)
UNTYPED_GLOBAL = re.compile(
    r"^(?P<indent> *)(?P<head>[A-Za-z_][A-Za-z0-9_]*(?:\s+[A-Za-z_][A-Za-z0-9_]*)*)\s*$"
)
GROUPED_GLOBAL = re.compile(
    r"^(?P<indent> *)(?P<head>[A-Za-z_][A-Za-z0-9_]*(?:, [A-Za-z_][A-Za-z0-9_]*)+)\s+:\s+(?P<type>\S.*)$"
)
RUNTIME_LABEL = re.compile(r"^(?P<indent> *)Runtime requirements:\s*$")
LIST_ITEM = re.compile(r"^(?P<indent> *)(?P<marker>[-+])\s+(?P<text>\S.*)$")
NOTES_PSEUDO_HEADER = re.compile(r"^[A-Z][^.!?]{0,60}:$")
UNSUPPORTED_DEPENDENCY = re.compile(r"^Dependenc(?:y|ies):?$")
GLOBAL_SUBGROUP = re.compile(r"^(?:Read|Write|Required|Optional|Inputs?|Outputs?):$")
PROSE_BACKTICKS = re.compile(r"`[^`\n]+`")
CANONICAL_UNKNOWN = "echo_err \"unknown option/parameter passed: '${1}'.\""
UNKNOWN_MESSAGE = re.compile(
    r"(?:Unknown argument passed|unknown (?:argument|option|option/parameter) passed)",
    re.I,
)
APPROVED_GLOBAL_TYPE = re.compile(
    r"^(?:"
    r"str|int|num|float|bool|flag|file|dir|path|time|choice|array|"
    r"structured string|"
    r"(?:array|list) of (?:str|int|num|float|bool|file|dir|path|structured string)|"
    r"\{[^{}]+\}"
    r")$"
)
CATEGORY_NAMES = {
    "Programs",
    "Environment",
    "Files",
    "Shell scripts",
    "Python scripts",
    "Configuration files",
}
OLD_CATEGORY_NAMES = {
    "Recommended environment",
    "External programs",
    "Shell scripts",
    "Python scripts",
    "Configuration files",
}


@dataclasses.dataclass(frozen=True)
class Finding:
    """One structural shell-help finding."""

    rule_id: str
    path: str
    line: int
    message: str
    related_lines: tuple[int, ...] = ()

    def format(self) -> str:
        """Render one stable smoke-test diagnostic."""

        return f"{self.rule_id}: {self.path}:{self.line}: {self.message}"


@dataclasses.dataclass(frozen=True)
class Section:
    """One underlined section inside a recognized help heredoc."""

    name: str
    heading_line: int
    heading_indent: int
    lines: tuple[tuple[int, str], ...]


@dataclasses.dataclass(frozen=True)
class Advisory:
    """One non-failing bounded style warning."""

    rule_id: str
    path: str
    line: int
    message: str

    def format(self) -> str:
        """Render one stable advisory diagnostic."""

        return f"{self.rule_id}: {self.path}:{self.line}: {self.message}"


def case_arm_spacing_warnings(
    text: str,
    path: str = "<memory>",
) -> list[Advisory]:
    """Warn for a narrow adjacent standard-case-arm shape.

    This intentionally recognizes only a standalone ';;' followed immediately
    by another simple arm header. It does not interpret nested cases, shared
    patterns, fall-through terminators, or generated compact code.
    """

    lines = text.splitlines()
    warnings: list[Advisory] = []
    for index, line in enumerate(lines[:-1]):
        if line.strip() != ";;":
            continue
        following = lines[index + 1]
        if not re.fullmatch(r"\s*[A-Za-z0-9_$\"'{}:/.*?+|_-]+\)\s*", following):
            continue
        warnings.append(
            Advisory(
                rule_id="SHELL.CASE_ARM.SPACING",
                path=path,
                line=index + 2,
                message="prefer exactly one empty line between adjacent case arms",
            )
        )
    return warnings


def sections(heredoc: Heredoc) -> list[Section]:
    """Return bounded underlined section bodies."""

    body = list(heredoc.lines)
    starts: list[tuple[int, str, int]] = []
    for index, (_number, line) in enumerate(body[:-1]):
        name = line.strip()
        if name not in SECTION_NAMES:
            continue
        if not UNDERLINE.fullmatch(body[index + 1][1].strip()):
            continue
        starts.append((index, name, indentation(line)))

    result: list[Section] = []
    for position, (index, name, heading_indent) in enumerate(starts):
        end = starts[position + 1][0] if position + 1 < len(starts) else len(body)
        result.append(
            Section(
                name=name,
                heading_line=body[index][0],
                heading_indent=heading_indent,
                lines=tuple(body[index + 2 : end]),
            )
        )
    return result


def finding_enabled(lines: tuple[int, ...], active_lines: set[int] | None) -> bool:
    """Return whether one finding intersects the requested current lines."""

    return active_lines is None or bool(set(lines) & active_lines)


def add_finding(
    findings: list[Finding],
    rule_id: str,
    path: str,
    line: int,
    message: str,
    active_lines: set[int] | None,
    *related_lines: int,
) -> None:
    """Append one finding when its evidence intersects active lines."""

    evidence = (line, *related_lines)
    if finding_enabled(evidence, active_lines):
        findings.append(
            Finding(
                rule_id=rule_id,
                path=path,
                line=line,
                message=message,
                related_lines=tuple(related_lines),
            )
        )


def nonblank_groups(lines: tuple[tuple[int, str], ...]) -> list[list[tuple[int, str]]]:
    """Split a section body into blank-line-delimited groups."""

    groups: list[list[tuple[int, str]]] = []
    current: list[tuple[int, str]] = []
    for item in lines:
        if item[1].strip():
            current.append(item)
        elif current:
            groups.append(current)
            current = []
    if current:
        groups.append(current)
    return groups


def check_usage(
    section: Section,
    path: str,
    findings: list[Finding],
    active_lines: set[int] | None,
) -> None:
    """Check canonical help spelling and invocation-line structure."""

    groups = nonblank_groups(section.lines)
    for group in groups[:1]:
        first_number, first_line = group[0]
        first_indent = indentation(first_line)
        if (
            first_indent != section.heading_indent + 2
            or not COMMAND_NAME.fullmatch(first_line.strip())
        ):
            add_finding(
                findings,
                RULE_USAGE,
                path,
                first_number,
                "command or script name must appear alone at section base + 2 spaces",
                active_lines,
            )
        for number, line in group[1:]:
            if indentation(line) != first_indent + 2:
                add_finding(
                    findings,
                    RULE_USAGE,
                    path,
                    number,
                    "argument continuation must be exactly two spaces deeper than the command name",
                    active_lines,
                    first_number,
                )
        for number, line in group:
            if HIDDEN_HELP.search(line):
                add_finding(
                    findings,
                    RULE_CANONICAL_HELP,
                    path,
                    number,
                    "public Usage must document only [--help]",
                    active_lines,
                )
            if SHORT_OPTION.search(line):
                add_finding(
                    findings,
                    RULE_LONG_OPTIONS,
                    path,
                    number,
                    "public Usage must advertise canonical long option names only",
                    active_lines,
                )


def entry_match(section_name: str, line: str) -> re.Match[str] | None:
    """Match one typed entry row for a supported structured section."""

    if section_name == "Parameters":
        return PARAMETER_ENTRY.match(line)
    return GLOBAL_ENTRY.match(line)


def check_entries(
    section: Section,
    path: str,
    findings: list[Finding],
    active_lines: set[int] | None,
) -> None:
    """Check entry indentation, descriptions, types, and blank separation."""

    rows: list[tuple[int, int, str, bool]] = []
    globals_by_description: dict[tuple[str, str], list[int]] = {}
    body = list(section.lines)
    expected_row_indent = section.heading_indent + 2

    for index, (number, line) in enumerate(body):
        match = entry_match(section.name, line)
        missing_type = False
        if match is None and section.name in {"Expected globals", "Generated globals"}:
            match = GROUPED_GLOBAL.match(line)
        if match is None and section.name in {"Expected globals", "Generated globals"}:
            match = UNTYPED_GLOBAL.match(line)
            missing_type = match is not None
        if (
            match is None
            and section.name in {"Expected globals", "Generated globals"}
            and indentation(line) == expected_row_indent
            and ("," in line or line.count(" : ") > 1)
        ):
            add_finding(
                findings,
                RULE_GLOBAL_GROUP_SYNTAX,
                path,
                number,
                "grouped globals require valid identifiers separated by comma and one space before one shared type",
                active_lines,
            )
            continue
        if match is None:
            if (
                section.name in {"Expected globals", "Generated globals"}
                and indentation(line) == expected_row_indent
                and GLOBAL_SUBGROUP.fullmatch(line.strip())
            ):
                add_finding(
                    findings,
                    RULE_GLOBAL_SUBGROUP,
                    path,
                    number,
                    f"{section.name} entries must be flat; remove '{line.strip()}'",
                    active_lines,
                )
            continue

        row_indent = len(match.group("indent"))
        rows.append((index, number, line, missing_type))
        if row_indent != expected_row_indent:
            add_finding(
                findings,
                RULE_ENTRY_INDENT,
                path,
                number,
                "entry row must begin at section base + 2 spaces",
                active_lines,
            )
        if missing_type:
            add_finding(
                findings,
                RULE_GLOBAL_TYPE,
                path,
                number,
                f"{section.name} entry requires 'name : type'",
                active_lines,
            )
        elif (
            section.name in {"Expected globals", "Generated globals"}
            and (match.group("type").count(" : ") or re.search(r",\s*[A-Za-z_][A-Za-z0-9_]*\s+:\s+", match.group("type")))
        ):
            add_finding(
                findings,
                RULE_GLOBAL_GROUP_SYNTAX,
                path,
                number,
                "grouped globals require one shared type after the complete identifier list",
                active_lines,
            )
        elif (
            section.name in {"Expected globals", "Generated globals"}
            and not APPROVED_GLOBAL_TYPE.fullmatch(match.group("type").strip())
        ):
            add_finding(
                findings,
                RULE_GLOBAL_TYPE,
                path,
                number,
                f"{section.name} entry uses an unapproved type: '{match.group('type').strip()}'",
                active_lines,
            )
        elif re.search(r"\s{2,}\S", match.group("type")):
            add_finding(
                findings,
                RULE_ENTRY_INDENT,
                path,
                number,
                "entry description must begin on the next line at row indentation + 2 spaces",
                active_lines,
            )
        if section.name == "Parameters" and match.group("head").lstrip().startswith("-"):
            aliases = re.findall(r"--?[A-Za-z0-9][A-Za-z0-9_-]*", match.group("head"))
            if len(aliases) != len(set(aliases)):
                add_finding(
                    findings,
                    RULE_PARAMETER_ALIAS_DUPLICATE,
                    path,
                    number,
                    "a parameter row must list each alias exactly once",
                    active_lines,
                )

    row_indexes = {row[0] for row in rows}
    for index, number, _, _ in rows:
        next_index = index + 1
        if next_index >= len(body) or next_index in row_indexes:
            continue
        next_number, next_line = body[next_index]
        if not next_line.strip():
            continue
        if indentation(next_line) != expected_row_indent + 2:
            add_finding(
                findings,
                RULE_ENTRY_INDENT,
                path,
                next_number,
                "entry description must begin two spaces deeper than its row",
                active_lines,
                number,
            )

    for previous, current in pairwise(rows):
        previous_index, previous_number, _, _ = previous
        current_index, current_number, _, _ = current
        between = body[previous_index + 1 : current_index]
        if not any(not line.strip() for _, line in between):
            add_finding(
                findings,
                RULE_ENTRY_BLANK,
                path,
                current_number,
                "adjacent entries require one blank line between them",
                active_lines,
                previous_number,
            )

    if section.name not in {"Expected globals", "Generated globals"}:
        return

    for index, number, line, missing_type in rows:
        if missing_type:
            continue
        match = GROUPED_GLOBAL.match(line) or GLOBAL_ENTRY.match(line)
        if match is None:
            continue
        description_lines: list[str] = []
        for _, description in body[index + 1 :]:
            if not description.strip():
                break
            if entry_match(section.name, description) or GROUPED_GLOBAL.match(description):
                break
            description_lines.append(description.strip())
        description = " ".join(description_lines)
        type_name = match.group("type").strip()
        if GROUPED_GLOBAL.match(line):
            common_markers = (
                " respectively",
                "both ",
                "each ",
                "shared ",
            )
            lowered = f" {description.lower()}"
            if not any(marker in lowered for marker in common_markers):
                add_finding(
                    findings,
                    RULE_GLOBAL_MAPPING_REVIEW,
                    path,
                    number,
                    "grouped-global description needs 'respectively' or equally explicit shared wording",
                    active_lines,
                )
            continue
        key = (type_name, description)
        globals_by_description.setdefault(key, []).append(number)

    for (type_name, description), line_numbers in globals_by_description.items():
        if description and len(line_numbers) > 1:
            add_finding(
                findings,
                RULE_GLOBAL_REPEATED_DESCRIPTION,
                path,
                line_numbers[1],
                f"separate {section.name.lower()} entries repeat one '{type_name}' description and should be grouped",
                active_lines,
                line_numbers[0],
            )


def next_nonblank(
    body: list[tuple[int, str]],
    start: int,
) -> tuple[int, str] | None:
    """Return the next nonblank line after one body index."""

    for item in body[start + 1 :]:
        if item[1].strip():
            return item
    return None


def exact_callable_message(
    text: str,
    callables: set[str],
    concepts: dict[str, set[str]],
) -> str | None:
    """Return an exact-callable correction for one program requirement."""

    # Slurm allocation describes an execution resource, not the ``sbatch``
    # callable.  Worker-only submit interfaces may require the former without
    # submitting jobs themselves.
    if text.startswith("Slurm allocation"):
        return None

    for concept in sorted(concepts, key=len, reverse=True):
        if re.match(rf"^{re.escape(concept)}(?:\b|\s|[,;(])", text, re.I):
            choices = concepts[concept]
            if len(choices) == 1:
                callable_name = next(iter(choices))
                if not text.startswith(callable_name):
                    return f"use exact callable spelling '{callable_name}', not '{concept}'"

    token = re.match(r"[A-Za-z][A-Za-z0-9-]*", text)
    if token:
        observed = token.group(0)
        for callable_name in callables:
            if observed.lower() == callable_name.lower() and observed != callable_name:
                return f"use exact callable spelling '{callable_name}', not '{observed}'"
    return None


def check_runtime_blocks(
    heredoc: Heredoc,
    path: str,
    findings: list[Finding],
    active_lines: set[int] | None,
    callables: set[str],
    concepts: dict[str, set[str]],
) -> None:
    """Check relative runtime nesting, category shape, and callable spellings."""

    body = list(heredoc.lines)
    for label_index, (label_number, label_line) in enumerate(body):
        label = RUNTIME_LABEL.match(label_line)
        if label is None:
            continue
        label_indent = len(label.group("indent"))
        direct: list[tuple[int, str, bool, bool]] = []
        program_lines: list[tuple[int, str]] = []
        current_category = ""
        saw_requirement = False
        after_blank = False

        for index in range(label_index + 1, len(body)):
            number, line = body[index]
            stripped = line.strip()
            if not stripped:
                after_blank = True
                continue
            if (
                stripped in SECTION_NAMES
                and index + 1 < len(body)
                and UNDERLINE.fullmatch(body[index + 1][1].strip())
                and indentation(line) <= label_indent
            ):
                break

            item = LIST_ITEM.match(line)
            line_indent = indentation(line)
            if item is None:
                old_name = stripped[:-1] if stripped.endswith(":") else ""
                if old_name in OLD_CATEGORY_NAMES and line_indent > label_indent:
                    add_finding(
                        findings,
                        RULE_RUNTIME_STRUCTURE,
                        path,
                        number,
                        "runtime categories must use '- Category' with '+' nested items",
                        active_lines,
                        label_number,
                    )
                    current_category = old_name
                    after_blank = False
                    continue
                if line_indent <= label_indent:
                    break
                if line_indent == label_indent + 2:
                    direct.append((number, stripped, False, False))
                    program_lines.append((number, stripped))
                    saw_requirement = True
                after_blank = False
                continue

            marker = item.group("marker")
            text = item.group("text")
            item_indent = len(item.group("indent"))
            if after_blank and saw_requirement and item_indent == label_indent:
                break
            saw_requirement = True

            if marker == "-":
                is_category = text in CATEGORY_NAMES
                direct.append((number, text, is_category, True))
                current_category = text if is_category else ""
                if item_indent != label_indent + 2:
                    add_finding(
                        findings,
                        RULE_RUNTIME_INDENT,
                        path,
                        number,
                        "direct runtime bullet must be exactly two spaces deeper than its label",
                        active_lines,
                        label_number,
                    )
                if not is_category:
                    program_lines.append((number, text))
            else:
                if not current_category or item_indent != label_indent + 4:
                    add_finding(
                        findings,
                        RULE_RUNTIME_INDENT,
                        path,
                        number,
                        "nested runtime item must be exactly two spaces deeper than its category",
                        active_lines,
                        label_number,
                    )
                if current_category == "Programs":
                    program_lines.append((number, text))
            after_blank = False

        if len(direct) == 1 and direct[0][1] == "Programs" and direct[0][2]:
            add_finding(
                findings,
                RULE_RUNTIME_PROGRAMS,
                path,
                direct[0][0],
                "use a flat runtime list when Programs is the only requirement category",
                active_lines,
                label_number,
            )

        if len(direct) == 1 and direct[0][3] and not direct[0][2]:
            add_finding(
                findings,
                RULE_RUNTIME_CARDINALITY,
                path,
                direct[0][0],
                "a single runtime requirement must be unbulleted",
                active_lines,
                label_number,
            )
        elif len(direct) >= 2:
            for number, _, is_category, bulleted in direct:
                if not bulleted and not is_category:
                    add_finding(
                        findings,
                        RULE_RUNTIME_CARDINALITY,
                        path,
                        number,
                        "two or more peer runtime requirements must use bullets",
                        active_lines,
                        label_number,
                    )

        for number, text in program_lines:
            message = exact_callable_message(text, callables, concepts)
            if message:
                add_finding(
                    findings,
                    RULE_RUNTIME_CALLABLE,
                    path,
                    number,
                    message,
                    active_lines,
                )


def is_function_help(heredoc: Heredoc) -> bool:
    """Return whether a bounded heredoc documents a shell function."""

    return heredoc.owner not in {"<file>", "<rendered>", "main"} and not heredoc.owner.startswith(
        ("help_", "show_help", "detail_")
    )


def usage_documents_arguments(section: Section) -> bool:
    """Return whether Usage contains arguments beyond the help flag."""

    groups = nonblank_groups(section.lines)
    lines = [line.strip() for _, line in groups[0]] if groups else []
    if len(lines) < 2:
        return False
    synopsis = " ".join(lines[1:])
    synopsis = synopsis.replace("[--help]", "").replace("--help", "")
    return bool(synopsis.strip(" []()|"))


def check_document_structure(
    heredoc: Heredoc,
    path: str,
    findings: list[Finding],
    active_lines: set[int] | None,
) -> None:
    """Check headings, relative indentation, Notes labels, and function order."""

    parsed = sections(heredoc)
    body = list(heredoc.lines)

    for section in parsed:
        heading_index = next(
            index
            for index, (number, _) in enumerate(body)
            if number == section.heading_line
        )
        underline_number, underline_line = body[heading_index + 1]
        if section.heading_indent != 0:
            add_finding(
                findings,
                RULE_SECTION_INDENT,
                path,
                section.heading_line,
                "top-level section headings must begin at column zero",
                active_lines,
            )
        if indentation(underline_line) != 0:
            add_finding(
                findings,
                RULE_SECTION_INDENT,
                path,
                underline_number,
                "top-level section underlines must begin at column zero",
                active_lines,
                section.heading_line,
            )
        for number, line in section.lines:
            if line.strip() and indentation(line) < section.heading_indent + 2:
                add_finding(
                    findings,
                    RULE_SECTION_INDENT,
                    path,
                    number,
                    "top-level section content must begin at section base + 2 spaces",
                    active_lines,
                    section.heading_line,
                )

    names: dict[str, list[int]] = {}
    for section in parsed:
        names.setdefault(section.name, []).append(section.heading_line)
    for name, lines in names.items():
        for duplicate in lines[1:]:
            add_finding(
                findings,
                RULE_SECTION_DUPLICATE,
                path,
                duplicate,
                f"top-level heading '{name}' appears more than once in one help unit",
                active_lines,
                lines[0],
            )

    for index, (number, line) in enumerate(body):
        stripped = line.strip()
        if not UNSUPPORTED_DEPENDENCY.fullmatch(stripped):
            continue
        underlined = (
            index + 1 < len(body)
            and UNDERLINE.fullmatch(body[index + 1][1].strip()) is not None
        )
        if stripped.endswith(":") or underlined:
            add_finding(
                findings,
                RULE_UNSUPPORTED_DEPENDENCY,
                path,
                number,
                "replace Dependency/Dependencies with Notes and 'Runtime requirements:'",
                active_lines,
            )

    for section in parsed:
        if section.name != "Notes":
            continue
        lines = list(section.lines)
        for index, (number, line) in enumerate(lines):
            if (
                indentation(line) != section.heading_indent + 2
                or not NOTES_PSEUDO_HEADER.fullmatch(line.strip())
            ):
                continue
            for content_number, content_line in lines[index + 1 :]:
                if not content_line.strip():
                    break
                if indentation(content_line) <= indentation(line):
                    add_finding(
                        findings,
                        RULE_NOTES_PSEUDO_INDENT,
                        path,
                        content_number,
                        "content beneath a top-level Notes pseudo-header must be exactly two spaces deeper",
                        active_lines,
                        number,
                    )
                    break

    if not is_function_help(heredoc):
        return

    first_by_name = {section.name: section for section in parsed}
    for required in ("Usage", "Returns"):
        if required not in first_by_name:
            add_finding(
                findings,
                RULE_SECTION_REQUIRED,
                path,
                heredoc.start_line,
                f"function help requires a '{required}' section",
                active_lines,
                *(number for number, _ in heredoc.lines),
            )
    usage = first_by_name.get("Usage")
    if (
        usage is not None
        and usage_documents_arguments(usage)
        and "Parameters" not in first_by_name
    ):
        add_finding(
            findings,
            RULE_SECTION_REQUIRED,
            path,
            usage.heading_line,
            "function help documenting arguments requires a 'Parameters' section",
            active_lines,
            *(number for number, _ in usage.lines),
        )

    ordered = [section for section in parsed if section.name in FUNCTION_SECTION_RANK]
    for later_index, later in enumerate(ordered):
        later_rank = FUNCTION_SECTION_RANK[later.name]
        conflicts = [
            earlier
            for earlier in ordered[:later_index]
            if FUNCTION_SECTION_RANK[earlier.name] > later_rank
        ]
        if conflicts:
            earlier = conflicts[-1]
            add_finding(
                findings,
                RULE_SECTION_ORDER,
                path,
                later.heading_line,
                f"'{later.name}' must precede '{earlier.name}' in function help",
                active_lines,
                earlier.heading_line,
            )

    examples = first_by_name.get("Examples")
    if examples is not None and parsed[-1].name != "Examples":
        add_finding(
            findings,
            RULE_SECTION_ORDER,
            path,
            parsed[-1].heading_line,
            f"'{parsed[-1].name}' may not appear after 'Examples'; Examples must be final",
            active_lines,
            examples.heading_line,
        )


def check_long_options_and_prose(
    heredoc: Heredoc,
    path: str,
    findings: list[Finding],
    active_lines: set[int] | None,
) -> None:
    """Check example option spellings and clear prose backtick delimiters."""

    for section in sections(heredoc):
        if section.name == "Examples":
            continue
        if section.name == "Usage":
            continue
        for number, line in section.lines:
            if PROSE_BACKTICKS.search(line):
                add_finding(
                    findings,
                    RULE_PROSE_QUOTES,
                    path,
                    number,
                    "use straight single quotes as prose delimiters in shell help",
                    active_lines,
                )


def check_parser_messages(
    text: str,
    path: str,
    findings: list[Finding],
    active_lines: set[int] | None,
) -> None:
    """Reject bounded obsolete unknown-option diagnostics and stale TODOs."""

    for number, line in enumerate(text.splitlines(), 1):
        if (
            re.search(r"\becho(?:_err)?\b", line)
            and UNKNOWN_MESSAGE.search(line)
            and CANONICAL_UNKNOWN not in line
        ):
            add_finding(
                findings,
                RULE_PARSER_UNKNOWN,
                path,
                number,
                f"use exactly: {CANONICAL_UNKNOWN}",
                active_lines,
            )
        elif re.search(r"TODO.*(?:Unknown argument|unknown option|echo_err)", line, re.I):
            add_finding(
                findings,
                RULE_PARSER_UNKNOWN,
                path,
                number,
                "remove the stale TODO for the settled unknown-option diagnostic",
                active_lines,
            )


def check_python_docstrings(
    text: str,
    path: str,
    active_lines: set[int] | None,
) -> list[Finding]:
    """Check clear prose backtick delimiters in changed Python docstrings."""

    try:
        tree = ast.parse(text)
    except SyntaxError:
        return []
    source_lines = text.splitlines()
    findings: list[Finding] = []
    for node in ast.walk(tree):
        if not isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            continue
        body = getattr(node, "body", [])
        if not body:
            continue
        expression = body[0]
        if not (
            isinstance(expression, ast.Expr)
            and isinstance(expression.value, ast.Constant)
            and isinstance(expression.value.value, str)
        ):
            continue
        start = expression.lineno
        end = expression.end_lineno or start
        in_examples = False
        pending_examples = False
        in_fence = False
        for number in range(start, end + 1):
            line = source_lines[number - 1]
            stripped = line.strip().strip('"\'')
            if stripped == "Examples":
                pending_examples = True
                continue
            if pending_examples and UNDERLINE.fullmatch(stripped):
                in_examples = True
                pending_examples = False
                continue
            if in_examples and stripped and number < end:
                next_stripped = source_lines[number].strip().strip('"\'')
                if (
                    stripped in SECTION_NAMES
                    and UNDERLINE.fullmatch(next_stripped)
                ):
                    in_examples = False
            if stripped.startswith(("```", "'''")):
                in_fence = not in_fence
                continue
            if in_examples or in_fence or stripped.startswith(">>>"):
                continue
            if re.search(r"``[^`\n]+``|(?<!`)`[^`\n]+`(?!`)", line):
                add_finding(
                    findings,
                    RULE_PYTHON_PROSE_QUOTES,
                    path,
                    number,
                    "use straight single quotes as prose delimiters in changed Python docstrings",
                    active_lines,
                )
    return findings


def load_registry(registry_path: Path | None) -> tuple[set[str], dict[str, set[str]]]:
    """Load the existing exact-callable registry when supplied."""

    if registry_path is None:
        return set(), {}
    import json

    return validate_command_registry(
        json.loads(registry_path.read_text(encoding="utf-8"))
    )


def check_help_source(
    text: str,
    path: str = "<memory>",
    *,
    registry_path: Path | None = None,
    active_lines: set[int] | None = None,
) -> list[Finding]:
    """Return bounded structural findings for recognized help heredocs."""

    callables, concepts = load_registry(registry_path)
    findings: list[Finding] = []
    for heredoc in extract_help_heredocs(text):
        heredoc_lines = set(range(heredoc.start_line, heredoc.end_line + 1))
        if active_lines is None:
            heredoc_active = None
        elif heredoc_lines & active_lines:
            heredoc_active = heredoc_lines
        else:
            heredoc_active = set()
        check_document_structure(heredoc, path, findings, heredoc_active)
        for section in sections(heredoc):
            if section.name == "Usage":
                check_usage(section, path, findings, heredoc_active)
            elif section.name in {
                "Parameters",
                "Expected globals",
                "Generated globals",
            }:
                check_entries(section, path, findings, heredoc_active)
        check_runtime_blocks(
            heredoc,
            path,
            findings,
            heredoc_active,
            callables,
            concepts,
        )
        check_long_options_and_prose(heredoc, path, findings, heredoc_active)
    check_parser_messages(text, path, findings, active_lines)
    return findings


def check_rendered_help(
    text: str,
    path: str = "<rendered>",
    *,
    registry_path: Path | None = None,
) -> list[Finding]:
    """Check one already-rendered sectioned help document."""

    lines = tuple(enumerate(text.splitlines(), 1))
    heredoc = Heredoc(
        owner="<rendered>",
        delimiter="<rendered>",
        start_line=1,
        end_line=max(1, len(lines)),
        lines=lines,
    )
    callables, concepts = load_registry(registry_path)
    findings: list[Finding] = []
    check_document_structure(heredoc, path, findings, None)
    for section in sections(heredoc):
        if section.name == "Usage":
            check_usage(section, path, findings, None)
        elif section.name in {
            "Parameters",
            "Expected globals",
            "Generated globals",
        }:
            check_entries(section, path, findings, None)
    check_runtime_blocks(heredoc, path, findings, None, callables, concepts)
    check_long_options_and_prose(heredoc, path, findings, None)
    findings.extend(audit_rendered_example_commands(text, path).findings)
    return findings


def scan_repository(
    root: Path,
    paths: Iterable[str] | None = None,
    *,
    registry_path: Path | None = None,
) -> list[Finding]:
    """Scan changed lines in changed and untracked shell help."""

    root = root.resolve()
    if paths is not None:
        selected = list(paths)
    else:
        tracked_py = run_git(
            root,
            ["diff", "--name-only", "--diff-filter=ACMR", "HEAD", "--", "*.py"],
        ).stdout.splitlines()
        untracked_py = run_git(
            root,
            ["ls-files", "--others", "--exclude-standard", "--", "*.py"],
        ).stdout.splitlines()
        selected = shell_paths(root) + tracked_py + untracked_py
    registry = registry_path or root / "dev/config/command_names.json"
    findings: list[Finding] = []
    for path in sorted(set(selected)):
        target = root / path
        if not target.is_file():
            continue
        if path.endswith(".py"):
            findings.extend(
                check_python_docstrings(
                    target.read_text(encoding="utf-8"),
                    path,
                    changed_lines(root, path),
                )
            )
            continue
        findings.extend(
            check_help_source(
                target.read_text(encoding="utf-8"),
                path,
                registry_path=registry,
                active_lines=changed_lines(root, path),
            )
        )
    return findings


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse the bounded structural-checker command line."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--rendered", action="store_true")
    parser.add_argument("--case-spacing", action="store_true")
    parser.add_argument("paths", nargs="*")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """Print findings and return nonzero when structural rules fail."""

    args = parse_args(argv)
    root = args.root.resolve()
    registry = root / "dev/config/command_names.json"
    if args.case_spacing:
        warnings = [
            warning
            for path in (args.paths or shell_paths(root))
            if (root / path).is_file()
            for warning in case_arm_spacing_warnings(
                (root / path).read_text(encoding="utf-8"),
                path,
            )
        ]
        for warning in warnings:
            print(warning.format())
        print(f"SHELL.CASE_ARM.SPACING: {len(warnings)} advisory warning(s)")
        return 0
    if args.rendered:
        findings = [
            finding
            for path in args.paths
            for finding in check_rendered_help(
                Path(path).read_text(encoding="utf-8"),
                path,
                registry_path=registry,
            )
        ]
    else:
        findings = scan_repository(root, args.paths or None, registry_path=registry)
    for finding in findings:
        print(finding.format())
    if findings:
        print(f"HELP.STRUCTURE: {len(findings)} changed/new violation(s)")
        return 1
    print("HELP.STRUCTURE: pass (zero changed/new violations)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
