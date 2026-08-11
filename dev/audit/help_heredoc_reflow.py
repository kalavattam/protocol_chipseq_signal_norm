#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: help_heredoc_reflow.py
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
Audit physical source and rendered structure in shell help heredocs.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import re
import subprocess
import sys
from collections.abc import Iterable
from itertools import pairwise
from pathlib import Path

from dev.audit.fixture_paths import is_fixture_recipe

RULE_ID = "HELP.HEREDOC.SOURCE_REFLOW"
RULE_EXAMPLE_SOURCE_BACKSLASH = "HELP.EXAMPLES.COMMAND.SOURCE_BACKSLASH"
RULE_EXAMPLE_INDENT = "HELP.EXAMPLES.COMMAND.INDENT"
RULE_EXAMPLE_TAB_INDENT = "HELP.EXAMPLES.COMMAND.TAB_INDENT"
RULE_EXAMPLE_TRAILING_WHITESPACE = "HELP.EXAMPLES.COMMAND.TRAILING_WHITESPACE"
RULE_EXAMPLE_FINAL_CONTINUATION = "HELP.EXAMPLES.COMMAND.FINAL_CONTINUATION"
RULE_EXAMPLE_FENCE = "HELP.EXAMPLES.FENCE"
RULE_EXAMPLE_RENDERED_STRUCTURE = "HELP.EXAMPLES.COMMAND.RENDERED_STRUCTURE"
RULE_EXAMPLE_COLLAPSED = "HELP.EXAMPLES.COMMAND.COLLAPSED"
SECTION_NAMES = {
    "Usage",
    "Parameters",
    "Returns",
    "Output",
    "Notes",
    "References",
    "See Also",
    "Examples",
    "Expected globals",
    "Generated globals",
}
HEREDOC_START = re.compile(
    r"<<(?P<strip>-)?[ \t]*(?P<quote>['\"]?)"
    r"(?P<delimiter>[A-Za-z_][A-Za-z0-9_]*)(?P=quote)(?![A-Za-z0-9_])",
)
FUNCTION_START = re.compile(
    r"^\s*(?:function\s+)?(?P<name>[A-Za-z_][A-Za-z0-9_]*)\s*\(\s*\)\s*\{",
)
PARAMETER_ROW = re.compile(r"^\s+(?:\d+\+?\s+\S+|-{1,2}\S(?:.*?))\s+:\s+\S")
LIST_ROW = re.compile(r"^(?P<indent>\s*)(?:[-+*]|\d+[.)])\s+\S")
UNDERLINE = re.compile(r"^-{3,}$")
EXAMPLE_ENTRY = re.compile(r"^(?P<indent> *)\d+\.\s+\S")
FUNCTION_DEFINITION = re.compile(
    r"^(?:function\s+)?[A-Za-z_][A-Za-z0-9_]*\s*\(\s*\)\s*\{",
)
ASSIGNMENT = re.compile(
    r"^(?:export\s+|local\s+|readonly\s+|declare\s+|typeset\s+)?"
    r"[A-Za-z_][A-Za-z0-9_]*(?:\[[^]]+\])?\+?=",
)
COMPLEX_HEAD = re.compile(
    r"^(?:if|elif|else|then|fi|for|while|until|do|done|case|esac|select)\b",
)
COLLAPSED_CONTINUATION = re.compile(r"\\+[ \t]+--?[A-Za-z0-9]")


@dataclasses.dataclass(frozen=True)
class Heredoc:
    """
    One narrowly recognized help heredoc.
    """

    owner: str
    delimiter: str
    start_line: int
    end_line: int
    lines: tuple[tuple[int, str], ...]
    quoted: bool = False
    strips_tabs: bool = False
    opener: str = ""


@dataclasses.dataclass(frozen=True)
class ExampleCommandFinding:
    """
    One stable source or rendered example-command diagnostic.
    """

    rule_id: str
    path: str
    line: int
    owner: str
    message: str

    def format(self) -> str:
        """
        Render one line-precise diagnostic.
        """

        return (
            f"{self.rule_id}: {self.path}:{self.line}: "
            f"owner={self.owner}; {self.message}"
        )


@dataclasses.dataclass(frozen=True)
class ExampleCommandInventory:
    """
    One explicit command-block classification.
    """

    path: str
    owner: str
    heredoc_start_line: int
    example_line: int
    start_line: int
    end_line: int
    classification: str
    reason: str
    representation: str
    quoted_heredoc: bool

    def as_dict(self) -> dict[str, object]:
        """
        Return a deterministic JSON-ready inventory row.
        """

        return dataclasses.asdict(self)


@dataclasses.dataclass(frozen=True)
class ExampleCommandAudit:
    """
    Findings and explicit classifications for one audit scope.
    """

    findings: list[ExampleCommandFinding]
    inventory: list[ExampleCommandInventory]


@dataclasses.dataclass(frozen=True)
class Finding:
    """
    One consecutive ordinary-prose source-wrap finding.
    """

    path: str
    owner: str
    delimiter: str
    boundary_lines: tuple[int, int]
    physical_lines: tuple[int, int]
    offending_lines: tuple[str, str]
    rendered_prose: str

    def format(self) -> str:
        """
        Render a stable, line-oriented strict-gate diagnostic.
        """

        start, end = self.boundary_lines
        first, second = self.physical_lines

        return (
            f"{RULE_ID}: {self.path}:{first}-{second}: owner={self.owner}; "
            f"heredoc={self.delimiter}@{start}-{end}; "
            f"physical={self.offending_lines[0]!r} | "
            f"{self.offending_lines[1]!r}; "
            f"rendered={self.rendered_prose!r}"
        )


def run_git(
    root: Path,
    args: list[str],
    *,
    check: bool = True,
) -> subprocess.CompletedProcess[str]:
    """
    Run one read-only Git query in 'root'.
    """

    return subprocess.run(
        ["git", "-C", str(root), *args],
        text=True,
        capture_output=True,
        check=check,
    )


def shell_paths(root: Path) -> list[str]:
    """
    Return changed tracked and nonignored untracked shell paths.
    """

    tracked = run_git(
        root,
        ["diff", "--name-only", "--diff-filter=ACMR", "HEAD", "--", "*.sh"],
    ).stdout.splitlines()
    untracked = run_git(
        root,
        ["ls-files", "--others", "--exclude-standard", "--", "*.sh"],
    ).stdout.splitlines()

    # A fixture recipe has no help surface of its own. It is tracked, so this
    # discovery reaches it, and its heredoc bodies are fixture payload: a
    # recipe writing a help-surface fixture would otherwise be read as though
    # the recipe itself declared that help, at the indentation the heredoc
    # happens to sit at.
    return sorted(
        path
        for path in set(tracked) | set(untracked)
        if not is_fixture_recipe(path)
    )


def tracked_at_head(root: Path, path: str) -> bool:
    """
    Return whether 'path' exists as a blob at 'HEAD'.
    """

    result = run_git(root, ["cat-file", "-e", f"HEAD:{path}"], check=False)

    return result.returncode == 0


def changed_lines(root: Path, path: str) -> set[int]:
    """
    Return current-side line numbers changed from 'HEAD'.
    """

    if not tracked_at_head(root, path):
        text = (root / path).read_text(encoding="utf-8")
        return set(range(1, len(text.splitlines()) + 1))

    diff = run_git(
        root,
        ["diff", "--unified=0", "--no-ext-diff", "HEAD", "--", path],
    ).stdout
    changed: set[int] = set()

    for match in re.finditer(
        r"^@@ -\d+(?:,\d+)? \+(\d+)(?:,(\d+))? @@",
        diff,
        re.MULTILINE,
    ):
        start = int(match.group(1))
        count = int(match.group(2) or "1")
        changed.update(range(start, start + count))

    return changed


def is_help_body(
    lines: list[tuple[int, str]],
    owner: str,
    start_source: str,
) -> bool:
    """
    Recognize sectioned help without attempting a general shell parse.
    """

    body = [line for _, line in lines]
    sectioned = any(
        line in SECTION_NAMES
        and index + 1 < len(body)
        and bool(UNDERLINE.fullmatch(body[index + 1]))
        for index, line in enumerate(body)
    )
    context_named = (
        owner.startswith(("help_", "show_help", "detail_"))
        or "help" in start_source.lower()
    )

    return sectioned and (context_named or "cat" in start_source)


def extract_help_heredocs(text: str) -> list[Heredoc]:
    """
    Extract bounded, literal-delimiter, sectioned help heredocs.

    Parameters
    ----------
    text : str
        Complete Shell source.

    Returns
    -------
    heredocs : list[Heredoc]
        Recognized literal-delimiter help heredocs with owners and source rows.
    """

    source = text.splitlines()
    owner = "<file>"
    heredocs: list[Heredoc] = []
    index = 0

    while index < len(source):
        function = FUNCTION_START.match(source[index])

        if function:
            owner = function.group("name")

        start = HEREDOC_START.search(source[index])

        if not start:
            index += 1

            continue

        delimiter = start.group("delimiter")
        quoted = bool(start.group("quote"))
        strips_tabs = bool(start.group("strip"))
        end = index + 1
        body: list[tuple[int, str]] = []

        while end < len(source):
            closing = source[end]
            if closing == delimiter or (
                strips_tabs and closing.lstrip("\t") == delimiter
            ):
                break

            body.append((end + 1, source[end]))
            end += 1

        if end >= len(source):
            index += 1

            continue

        if is_help_body(body, owner, source[index]):
            heredocs.append(
                Heredoc(
                    owner=owner,
                    delimiter=delimiter,
                    start_line=index + 1,
                    end_line=end + 1,
                    lines=tuple(body),
                    quoted=quoted,
                    strips_tabs=strips_tabs,
                    opener=source[index],
                ),
            )

        index = end + 1

    return heredocs


def _example_section_rows(
    heredoc: Heredoc,
) -> list[tuple[tuple[int, str], ...]]:
    """
    Return bounded Examples section bodies without owning owner discovery.
    """

    body = list(heredoc.lines)
    starts: list[tuple[int, str]] = []

    for index, (_, line) in enumerate(body[:-1]):
        name = line.strip()

        if name in SECTION_NAMES and UNDERLINE.fullmatch(
            body[index + 1][1].strip(),
        ):
            starts.append((index, name))

    result: list[tuple[tuple[int, str], ...]] = []

    for position, (index, name) in enumerate(starts):
        if name != "Examples":
            continue

        end = (
            starts[position + 1][0]
            if position + 1 < len(starts)
            else len(body)
        )
        result.append(tuple(body[index + 2 : end]))

    return result


def extract_rendered_examples_text(text: str) -> str:
    """
    Return one rendered Examples section, including its heading and underline.
    """

    lines = text.splitlines()
    starts: list[tuple[int, str]] = []

    for index, line in enumerate(lines[:-1]):
        name = line.strip()

        if name in SECTION_NAMES and UNDERLINE.fullmatch(
            lines[index + 1].strip(),
        ):
            starts.append((index, name))

    examples: list[str] = []

    for position, (index, name) in enumerate(starts):
        if name != "Examples":
            continue

        end = (
            starts[position + 1][0]
            if position + 1 < len(starts)
            else len(lines)
        )
        examples.append("\n".join(lines[index:end]).rstrip("\n") + "\n")

    if len(examples) != 1:
        raise ValueError(
            (
                f"expected exactly one rendered Examples section; found "
                f"{len(examples)}"
            ),
        )

    return examples[0]


def _leading_whitespace(line: str) -> str:
    """
    Return the complete leading spaces-and-tabs prefix.
    """

    match = re.match(r"[ \t]*", line)
    return match.group(0) if match else ""


def _terminal_backslashes(line: str) -> tuple[int, str]:
    """
    Return terminal backslash count and whitespace that follows it.
    """

    match = re.search(r"(?P<slashes>\\+)(?P<trailing>[ \t]*)$", line)
    if match is None:
        return 0, ""

    return len(match.group("slashes")), match.group("trailing")


def _complex_reason(code: list[tuple[int, str]]) -> str | None:
    """
    Return a narrow documented exclusion reason for complex shell syntax.

    Parameters
    ----------
    code : list[tuple[int, str]]
        Numbered physical lines from one candidate example block.

    Returns
    -------
    reason : str | None
        The first applicable exclusion reason, or None for simple commands.
    """

    for _, line in code:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue

        if COMPLEX_HEAD.match(stripped):
            return (
                "shell control-flow syntax has separate indentation semantics"
            )

        if FUNCTION_DEFINITION.match(stripped):
            return "function definition has separate indentation semantics"

        if ASSIGNMENT.match(stripped):
            if "=(" in stripped:
                return "array assignment has separate indentation semantics"

            return (
                "assignment/setup snippet is not one ordinary CLI invocation"
            )

        if "<<" in stripped:
            return "nested heredoc has separate physical-line semantics"

        if stripped.startswith(("(", ")")) or "$(" in stripped:
            return "subshell syntax has separate indentation semantics"

        if re.search(r"(?:^|\s)(?:\||\|\||&&)(?:\s|$)", stripped):
            return (
                "structured pipeline or command list is deliberately excluded"
            )

        if re.search(r"(?:^|\s)(?:>>?|<)(?:\s|$)", stripped):
            return "structured redirection snippet is deliberately excluded"

        if re.search(r"(?:^|\s);(?:\s|$)", stripped):
            return "multi-command shell list is deliberately excluded"

    return None


def _split_simple_commands(
    code: list[tuple[int, str]],
) -> list[list[tuple[int, str]]]:
    """
    Split bounded non-complex code into independent simple commands.
    """

    commands: list[list[tuple[int, str]]] = []
    current: list[tuple[int, str]] = []
    base_indent = 0

    for row in code:
        _, line = row
        stripped = line.strip()

        if not stripped or stripped.startswith("#"):
            if current:
                commands.append(current)
                current = []

            continue

        whitespace = _leading_whitespace(line)
        indent = len(whitespace.replace("\t", ""))

        if not current:
            current = [row]
            base_indent = indent

            continue

        previous_slashes, _ = _terminal_backslashes(current[-1][1])
        continuation = (
            previous_slashes > 0
            or stripped.startswith("-")
            or indent > base_indent
            or "\t" in whitespace
        )

        if continuation:
            current.append(row)

            continue

        commands.append(current)
        current = [row]
        base_indent = indent

    if current:
        commands.append(current)

    return commands


def _finding(
    findings: list[ExampleCommandFinding],
    rule_id: str,
    path: str,
    line: int,
    owner: str,
    message: str,
) -> None:
    """
    Append one command diagnostic.
    """

    findings.append(ExampleCommandFinding(rule_id, path, line, owner, message))


def _audit_simple_command(
    command: list[tuple[int, str]],
    *,
    path: str,
    owner: str,
    quoted: bool,
    rendered: bool,
    example_line: int,
    heredoc_start_line: int,
    findings: list[ExampleCommandFinding],
    inventory: list[ExampleCommandInventory],
) -> None:
    """
    Apply exact relative indentation and representation invariants.

    Parameters
    ----------
    command : list[tuple[int, str]]
        Numbered physical lines forming one simple shell command.
    path : str
        Repository-relative path that owns the command.
    owner : str
        Help function or document owner used in diagnostics.
    quoted : bool
        Whether the source uses a quoted heredoc delimiter.
    rendered : bool
        Whether the command represents rendered rather than source text.
    example_line : int
        Source line on which the enclosing example starts.
    heredoc_start_line : int
        Source line on which the enclosing heredoc starts.
    findings : list[ExampleCommandFinding]
        Mutable collection receiving structural diagnostics.
    inventory : list[ExampleCommandInventory]
        Mutable collection receiving the command review record.
    """

    representation = "rendered" if rendered else "source"
    first_number, first_line = command[0]
    last_number = command[-1][0]
    inventory.append(
        ExampleCommandInventory(
            path=path,
            owner=owner,
            heredoc_start_line=heredoc_start_line,
            example_line=example_line,
            start_line=first_number,
            end_line=last_number,
            classification="checked_simple_command",
            reason=(
                "bounded ordinary command with physical argument continuations"
            ),
            representation=representation,
            quoted_heredoc=quoted,
        ),
    )
    base_whitespace = _leading_whitespace(first_line)
    base_indent = len(base_whitespace.replace("\t", ""))
    expected_slashes = 1 if rendered or quoted else 2

    for index, (number, line) in enumerate(command):
        is_final = index == len(command) - 1
        slash_count, trailing = _terminal_backslashes(line)

        if COLLAPSED_CONTINUATION.search(line):
            _finding(
                findings,
                RULE_EXAMPLE_COLLAPSED,
                path,
                number,
                owner,
                (
                    "continuation and following option are collapsed on one "
                    "physical line"
                ),
            )

            if rendered:
                _finding(
                    findings,
                    RULE_EXAMPLE_RENDERED_STRUCTURE,
                    path,
                    number,
                    owner,
                    (
                        "rendered command did not preserve the required "
                        "physical line boundary"
                    ),
                )

        if trailing:
            _finding(
                findings,
                RULE_EXAMPLE_TRAILING_WHITESPACE,
                path,
                number,
                owner,
                "spaces or tabs follow the terminal continuation backslash",
            )

        if is_final:
            if slash_count:
                _finding(
                    findings,
                    RULE_EXAMPLE_FINAL_CONTINUATION,
                    path,
                    number,
                    owner,
                    (
                        "final command line must not end with a continuation "
                        "backslash"
                    ),
                )

            continue

        if slash_count != expected_slashes:
            rule = (
                RULE_EXAMPLE_RENDERED_STRUCTURE
                if rendered
                else RULE_EXAMPLE_SOURCE_BACKSLASH
            )
            _finding(
                findings,
                rule,
                path,
                number,
                owner,
                (
                    f"non-final {representation} line requires exactly "
                    f"{expected_slashes} "
                    f"terminal backslash(es); observed "
                    f"{slash_count}"
                ),
            )

    for number, line in command[1:]:
        whitespace = _leading_whitespace(line)

        if "\t" in whitespace:
            _finding(
                findings,
                RULE_EXAMPLE_TAB_INDENT,
                path,
                number,
                owner,
                "tabs are forbidden in continuation indentation",
            )

            continue

        observed = len(whitespace)

        if observed != base_indent + 4:
            _finding(
                findings,
                RULE_EXAMPLE_INDENT,
                path,
                number,
                owner,
                (
                    f"continuation indentation must be exactly four spaces "
                    f"beyond base "
                    f"{base_indent}; observed "
                    f"{observed}"
                ),
            )


def _audit_examples_in_heredoc(
    heredoc: Heredoc,
    *,
    path: str,
    rendered: bool,
) -> ExampleCommandAudit:
    """
    Audit recognized Examples fences in one source or rendered unit.

    Parameters
    ----------
    heredoc : Heredoc
        Parsed help heredoc containing recognized section rows.
    path : str
        Repository-relative source path for diagnostics.
    rendered : bool
        Whether rows represent rendered help instead of physical Shell source.

    Returns
    -------
    audit : ExampleCommandAudit
        Deterministic command findings and complete command inventory.
    """

    findings: list[ExampleCommandFinding] = []
    inventory: list[ExampleCommandInventory] = []
    representation = "rendered" if rendered else "source"

    for section in _example_section_rows(heredoc):
        starts = [
            index
            for index, (_, line) in enumerate(section)
            if EXAMPLE_ENTRY.match(line)
        ]

        if not starts and any(line.strip() for _, line in section):
            number = next(number for number, line in section if line.strip())
            _finding(
                findings,
                RULE_EXAMPLE_FENCE,
                path,
                number,
                heredoc.owner,
                (
                    "recognized Examples content has no numbered entry/fence "
                    "structure"
                ),
            )
            inventory.append(
                ExampleCommandInventory(
                    path,
                    heredoc.owner,
                    heredoc.start_line,
                    number,
                    number,
                    number,
                    "malformed_or_ambiguous_block",
                    "Examples content has no numbered entry boundary",
                    representation,
                    heredoc.quoted,
                ),
            )

            continue

        for position, start in enumerate(starts):
            end = (
                starts[position + 1]
                if position + 1 < len(starts)
                else len(section)
            )
            group = list(section[start:end])
            example_line = group[0][0]

            if position and section[start - 1][1].strip():
                rule = (
                    RULE_EXAMPLE_RENDERED_STRUCTURE
                    if rendered
                    else RULE_EXAMPLE_FENCE
                )
                _finding(
                    findings,
                    rule,
                    path,
                    example_line,
                    heredoc.owner,
                    (
                        "numbered Examples entries require one separating "
                        "blank line"
                    ),
                )

            malformed_fences = [
                (number, line)
                for number, line in group[1:]
                if line.strip().startswith(("```", "\\```"))
            ]

            openings = [
                (index, number, line)
                for index, (number, line) in enumerate(group[1:], 1)
                if line.strip() == "'''bash"
            ]

            closings = [
                (index, number, line)
                for index, (number, line) in enumerate(group[1:], 1)
                if line.strip() == "'''"
            ]

            opening_indent = (
                _leading_whitespace(openings[0][2]) if openings else None
            )
            closing_indent = (
                _leading_whitespace(closings[0][2]) if closings else None
            )
            valid_pair = (
                len(openings) == 1
                and len(closings) == 1
                and closings[0][0] > openings[0][0]
                and opening_indent == closing_indent
            )

            if malformed_fences or not valid_pair:
                line = (
                    malformed_fences[0][0]
                    if malformed_fences
                    else example_line
                )
                _finding(
                    findings,
                    RULE_EXAMPLE_FENCE,
                    path,
                    line,
                    heredoc.owner,
                    (
                        "each recognized Examples entry requires one literal "
                        "'''bash ... ''' fence pair"
                    ),
                )
                inventory.append(
                    ExampleCommandInventory(
                        path,
                        heredoc.owner,
                        heredoc.start_line,
                        example_line,
                        line,
                        line,
                        "malformed_or_ambiguous_block",
                        (
                            "missing, changed, duplicated, or unterminated "
                            "literal fence"
                        ),
                        representation,
                        heredoc.quoted,
                    ),
                )

                continue

            opening, closing = openings[0], closings[0]
            fence_indent = len(_leading_whitespace(opening[2]))
            code = [
                (
                    number,
                    line[fence_indent:]
                    if line.startswith(" " * fence_indent)
                    else line,
                )
                for number, line in group[opening[0] + 1 : closing[0]]
            ]

            if not any(line.strip() for _, line in code):
                _finding(
                    findings,
                    RULE_EXAMPLE_FENCE,
                    path,
                    example_line,
                    heredoc.owner,
                    "literal Bash fence must contain nonempty code",
                )
                inventory.append(
                    ExampleCommandInventory(
                        path,
                        heredoc.owner,
                        heredoc.start_line,
                        example_line,
                        opening[1],
                        closing[1],
                        "malformed_or_ambiguous_block",
                        "empty literal Bash fence",
                        representation,
                        heredoc.quoted,
                    ),
                )

                continue

            reason = _complex_reason(code)

            if reason is not None:
                inventory.append(
                    ExampleCommandInventory(
                        path,
                        heredoc.owner,
                        heredoc.start_line,
                        example_line,
                        code[0][0],
                        code[-1][0],
                        "deliberately_excluded_complex_snippet",
                        reason,
                        representation,
                        heredoc.quoted,
                    ),
                )

                continue

            commands = _split_simple_commands(code)

            if not commands:
                inventory.append(
                    ExampleCommandInventory(
                        path,
                        heredoc.owner,
                        heredoc.start_line,
                        example_line,
                        opening[1],
                        closing[1],
                        "malformed_or_ambiguous_block",
                        (
                            "no classifiable command remains after comments "
                            "and blanks"
                        ),
                        representation,
                        heredoc.quoted,
                    ),
                )

                continue

            for command in commands:
                _audit_simple_command(
                    command,
                    path=path,
                    owner=heredoc.owner,
                    quoted=heredoc.quoted,
                    rendered=rendered,
                    example_line=example_line,
                    heredoc_start_line=heredoc.start_line,
                    findings=findings,
                    inventory=inventory,
                )

    findings.sort(
        key=lambda item: (item.path, item.line, item.rule_id, item.message),
    )
    inventory.sort(
        key=lambda item: (
            item.path,
            item.heredoc_start_line,
            item.example_line,
            item.start_line,
            item.classification,
        ),
    )

    return ExampleCommandAudit(findings, inventory)


def audit_source_example_commands(
    text: str,
    path: str = "<memory>",
) -> ExampleCommandAudit:
    """
    Audit every recognized source Examples block in one shell file.
    """

    findings: list[ExampleCommandFinding] = []
    inventory: list[ExampleCommandInventory] = []

    for heredoc in extract_help_heredocs(text):
        result = _audit_examples_in_heredoc(
            heredoc,
            path=path,
            rendered=False,
        )
        findings.extend(result.findings)
        inventory.extend(result.inventory)

    findings.sort(
        key=lambda item: (item.path, item.line, item.rule_id, item.message),
    )
    inventory.sort(
        key=lambda item: (
            item.path,
            item.heredoc_start_line,
            item.example_line,
            item.start_line,
            item.classification,
        ),
    )

    return ExampleCommandAudit(findings, inventory)


def audit_rendered_example_commands(
    text: str,
    path: str = "<rendered>",
) -> ExampleCommandAudit:
    """
    Audit physical command lines and literal fences in rendered help.
    """

    lines = tuple(enumerate(text.splitlines(), 1))
    heredoc = Heredoc(
        owner="<rendered>",
        delimiter="<rendered>",
        start_line=1,
        end_line=max(1, len(lines)),
        lines=lines,
        quoted=True,
        opener="<rendered>",
    )

    return _audit_examples_in_heredoc(
        heredoc,
        path=path,
        rendered=True,
    )


def indentation(line: str) -> int:
    """
    Return leading-space indentation for help prose comparison.
    """

    return len(line) - len(line.lstrip(" "))


def structured_lines(heredoc: Heredoc) -> set[int]:
    """
    Classify the small explicit structured-content exclusions.

    Parameters
    ----------
    heredoc : Heredoc
        Parsed help rows whose tables, lists, fences, and examples are scanned.

    Returns
    -------
    line_numbers : set[int]
        Source line numbers whose structure owns their physical wrapping.
    """

    lines = list(heredoc.lines)
    structured: set[int] = set()
    section = ""
    in_fence = False
    list_indent: int | None = None

    for index, (number, line) in enumerate(lines):
        stripped = line.strip()
        next_line = lines[index + 1][1] if index + 1 < len(lines) else ""

        if stripped in {"'''", "'''bash", "```", "```bash"}:
            structured.add(number)
            in_fence = not in_fence if stripped in {"'''", "```"} else True

            continue

        if in_fence:
            structured.add(number)

            continue

        if not stripped:
            structured.add(number)
            list_indent = None

            continue

        if stripped in SECTION_NAMES and UNDERLINE.fullmatch(
            next_line.strip(),
        ):
            section = stripped
            structured.add(number)

            continue

        if UNDERLINE.fullmatch(stripped):
            structured.add(number)

            continue

        if section == "Usage":
            structured.add(number)

            continue

        if PARAMETER_ROW.match(line):
            structured.add(number)

            continue

        list_match = LIST_ROW.match(line)

        if list_match:
            list_indent = len(list_match.group("indent"))
            structured.add(number)

            continue

        if list_indent is not None and indentation(line) > list_indent:
            structured.add(number)

            continue

        list_indent = None

        if "|" in line or "\t" in line or re.search(r"\S\s{2,}\S", stripped):
            structured.add(number)

            continue

        if re.match(
            r"^(?:[$>]\s*)?(?:bash|sh|python\d*|conda|git|printf|echo|cat)\b",
            stripped,
        ):
            structured.add(number)

            continue

        if re.match(r"^[A-Za-z_][A-Za-z0-9_]*=\S", stripped):
            structured.add(number)

            continue

        if len(stripped.split()) < 3 and not re.search(r"[.!?]$", stripped):
            structured.add(number)

    return structured


def normalized(line: str) -> str:
    """
    Collapse whitespace for exact one-line-to-two-line reflow evidence.
    """

    return " ".join(line.split())


def is_wrap_pair(first: str, second: str) -> bool:
    """
    Return whether two nonblank ordinary lines share one indentation.
    """

    return (
        bool(first.strip())
        and bool(second.strip())
        and indentation(first) == indentation(second)
    )


def findings_for_path(root: Path, path: str) -> list[Finding]:
    """
    Return strict changed-content findings for one current shell path.
    """

    target = root / path
    if not target.is_file():
        return []

    text = target.read_text(encoding="utf-8")
    changed = changed_lines(root, path)
    findings: list[Finding] = []

    for heredoc in extract_help_heredocs(text):
        structured = structured_lines(heredoc)
        body = list(heredoc.lines)

        for (first_number, first), (second_number, second) in pairwise(body):
            if not ({first_number, second_number} & changed):
                continue

            if first_number in structured or second_number in structured:
                continue

            if not is_wrap_pair(first, second):
                continue

            findings.append(
                Finding(
                    path=path,
                    owner=heredoc.owner,
                    delimiter=heredoc.delimiter,
                    boundary_lines=(heredoc.start_line, heredoc.end_line),
                    physical_lines=(first_number, second_number),
                    offending_lines=(first, second),
                    rendered_prose=normalized(
                        f"{first.strip()} {second.strip()}",
                    ),
                ),
            )

    return findings


def scan_repository(
    root: Path,
    paths: Iterable[str] | None = None,
) -> list[Finding]:
    """
    Scan changed/new shell paths while excluding untouched historical prose.
    """

    root = root.resolve()
    selected = list(paths) if paths is not None else shell_paths(root)

    return [
        finding
        for path in sorted(set(selected))
        for finding in findings_for_path(root, path)
    ]


def repository_shell_paths(
    root: Path,
    *,
    production_only: bool = False,
) -> list[str]:
    """
    Return maintained Git-visible shell sources in deterministic order.
    """

    tracked = run_git(root, ["ls-files", "--", "*.sh"]).stdout.splitlines()
    untracked = run_git(
        root,
        ["ls-files", "--others", "--exclude-standard", "--", "*.sh"],
    ).stdout.splitlines()
    prefixes = ("bin/", "lib/bash/", "install/scripts/")

    if not production_only:
        prefixes += ("tests/",)

    return sorted(
        {
            path
            for path in (*tracked, *untracked)
            if path.startswith(prefixes)
            and not path.startswith("artifacts/tests/")
            and (root / path).is_file()
        },
    )


def audit_example_command_repository(
    root: Path,
    paths: Iterable[str] | None = None,
    *,
    production_only: bool = False,
) -> ExampleCommandAudit:
    """
    Audit every maintained source help example in the selected scope.
    """

    root = root.resolve()
    selected = (
        sorted(set(paths))
        if paths is not None
        else repository_shell_paths(root, production_only=production_only)
    )
    findings: list[ExampleCommandFinding] = []
    inventory: list[ExampleCommandInventory] = []

    for path in selected:
        target = root / path
        if not target.is_file():
            continue

        audit = audit_source_example_commands(
            target.read_text(encoding="utf-8"),
            path,
        )
        findings.extend(audit.findings)
        inventory.extend(audit.inventory)

    findings.sort(
        key=lambda item: (item.path, item.line, item.rule_id, item.message),
    )
    inventory.sort(
        key=lambda item: (
            item.path,
            item.heredoc_start_line,
            item.example_line,
            item.start_line,
            item.classification,
        ),
    )

    return ExampleCommandAudit(findings, inventory)


def example_inventory_summary(
    audit: ExampleCommandAudit,
) -> dict[str, object]:
    """
    Return deterministic source/rendered command classification totals.
    """

    checked = [
        row
        for row in audit.inventory
        if row.classification == "checked_simple_command"
    ]
    excluded = [
        row
        for row in audit.inventory
        if row.classification == "deliberately_excluded_complex_snippet"
    ]
    ambiguous = [
        row
        for row in audit.inventory
        if row.classification == "malformed_or_ambiguous_block"
    ]
    affected = sorted({finding.path for finding in audit.findings})
    owners = sorted({row.owner for row in audit.inventory})

    return {
        "applicable_simple_commands": len(checked),
        "compliant_simple_commands": len(
            [
                row
                for row in checked
                if not any(
                    finding.path == row.path
                    and row.start_line <= finding.line <= row.end_line
                    for finding in audit.findings
                )
            ],
        ),
        "findings": len(audit.findings),
        "deliberately_excluded_complex_snippets": len(excluded),
        "ambiguous_snippets": len(ambiguous),
        "files_affected": affected,
        "owners_inventory_count": len(owners),
    }


def fix_source_example_commands(
    text: str,
    path: str = "<memory>",
) -> tuple[str, int, tuple[str, ...]]:
    """
    Conservatively normalize only classified simple source commands.

    Parameters
    ----------
    text : str
        Shell source containing help heredocs to inspect.
    path : str
        Diagnostic path associated with the source text.

    Returns
    -------
    corrected, changed_count, refusals : tuple[str, int, tuple[str, ...]]
        Corrected source, number of changed commands, and any refusal reasons.
    """

    audit = audit_source_example_commands(text, path)
    refused_rules = {RULE_EXAMPLE_COLLAPSED, RULE_EXAMPLE_FENCE}
    refused = tuple(
        sorted(
            {
                f"{finding.rule_id}@{finding.line}"
                for finding in audit.findings
                if finding.rule_id in refused_rules
            }
            | {
                f"ambiguous@{row.example_line}"
                for row in audit.inventory
                if row.classification == "malformed_or_ambiguous_block"
            },
        ),
    )
    if refused:
        return text, 0, refused

    lines = text.splitlines()
    changed_commands = 0

    for row in audit.inventory:
        if row.classification != "checked_simple_command":
            continue

        indexes = list(range(row.start_line - 1, row.end_line))
        if not indexes:
            continue

        original = [lines[index] for index in indexes]
        base_whitespace = _leading_whitespace(original[0])
        expected_slashes = 1 if row.quoted_heredoc else 2
        normalized: list[str] = []

        for position, line in enumerate(original):
            content = (
                line.lstrip(" \t")
                if position
                else line[len(base_whitespace) :]
            )
            content = re.sub(r"\\+[ \t]*$", "", content).rstrip(" \t")
            indentation_prefix = (
                base_whitespace
                if position == 0
                else " " * (len(base_whitespace) + 4)
            )

            if position != len(original) - 1:
                continuation = "\\" * expected_slashes
                content = f"{content} {continuation}"

            normalized.append(indentation_prefix + content)

        if normalized == original:
            continue

        for index, line in zip(indexes, normalized, strict=True):
            lines[index] = line

        changed_commands += 1

    suffix = "\n" if text.endswith("\n") else ""

    return "\n".join(lines) + suffix, changed_commands, ()


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse the strict-gate command line.

    Parameters
    ----------
    argv : list[str] | None
        Explicit arguments, or None to read the process arguments.

    Returns
    -------
    arguments : argparse.Namespace
        Parsed heredoc-audit and extraction options.
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--production-only", action="store_true")
    parser.add_argument("--prose-only", action="store_true")
    parser.add_argument("--rendered", action="store_true")
    parser.add_argument("--extract-rendered-examples", type=Path)
    parser.add_argument("--fix-example-commands", action="store_true")
    parser.add_argument("--example-inventory-output", type=Path)
    parser.add_argument("paths", nargs="*")

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """
    Print all strict findings and return nonzero when any exist.

    Parameters
    ----------
    argv : list[str] | None
        Explicit arguments, or None to read the process arguments.

    Returns
    -------
    status : int
        Zero when the selected heredoc operation succeeds and one otherwise.
    """

    args = parse_args(argv)

    if args.extract_rendered_examples is not None:
        try:
            sys.stdout.write(
                extract_rendered_examples_text(
                    args.extract_rendered_examples.read_text(encoding="utf-8"),
                ),
            )
        except ValueError as exc:
            print(
                (
                    f"{RULE_EXAMPLE_FENCE}: {args.extract_rendered_examples}: "
                    f"{exc}"
                ),
            )

            return 1

        return 0

    if args.fix_example_commands:
        root = args.root.resolve()
        selected = (
            sorted(set(args.paths))
            if args.paths
            else repository_shell_paths(
                root,
                production_only=args.production_only,
            )
        )
        refused: list[str] = []

        for path in selected:
            target = root / path
            if not target.is_file():
                continue

            before = target.read_text(encoding="utf-8")
            after, count, reasons = fix_source_example_commands(before, path)

            if reasons:
                refused.extend(f"{path}:{reason}" for reason in reasons)

                continue

            if after != before:
                target.write_text(after, encoding="utf-8")
                print(
                    (
                        f"HELP.EXAMPLES.COMMAND.FIX: {path}: normalized "
                        f"{count} "
                        f"command(s)"
                    ),
                )

        if refused:
            for item in refused:
                print(f"HELP.EXAMPLES.COMMAND.FIX_REFUSED: {item}")

            return 1

    prose_findings = (
        [] if args.rendered else scan_repository(args.root, args.paths or None)
    )
    example_audit = ExampleCommandAudit([], [])

    if args.rendered:
        findings: list[ExampleCommandFinding] = []
        inventory: list[ExampleCommandInventory] = []

        for path in args.paths:
            result = audit_rendered_example_commands(
                Path(path).read_text(encoding="utf-8"),
                path,
            )
            findings.extend(result.findings)
            inventory.extend(result.inventory)

        example_audit = ExampleCommandAudit(
            sorted(
                findings,
                key=lambda item: (
                    item.path,
                    item.line,
                    item.rule_id,
                    item.message,
                ),
            ),
            sorted(
                inventory,
                key=lambda item: (
                    item.path,
                    item.example_line,
                    item.start_line,
                    item.classification,
                ),
            ),
        )
    elif not args.prose_only:
        example_audit = audit_example_command_repository(
            args.root,
            args.paths or None,
            production_only=args.production_only,
        )

    if args.example_inventory_output is not None:
        args.example_inventory_output.write_text(
            json.dumps(
                {
                    "summary": example_inventory_summary(example_audit),
                    "inventory": [
                        row.as_dict() for row in example_audit.inventory
                    ],
                    "findings": [
                        dataclasses.asdict(row)
                        for row in example_audit.findings
                    ],
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )

    for finding in prose_findings:
        print(finding.format())

    for finding in example_audit.findings:
        print(finding.format())

    if prose_findings:
        print(
            (
                f"{RULE_ID}: {len(prose_findings)} newly introduced or "
                f"modified violation(s)"
            ),
        )
    else:
        print(
            f"{RULE_ID}: pass (zero newly introduced or modified violations)",
        )

    summary = example_inventory_summary(example_audit)

    if not args.prose_only:
        applicable = summary["applicable_simple_commands"]
        compliant = summary["compliant_simple_commands"]
        excluded = summary["deliberately_excluded_complex_snippets"]
        ambiguous = summary["ambiguous_snippets"]
        finding_count = summary["findings"]
        print(
            f"HELP.EXAMPLES.COMMAND: {applicable} applicable; "
            f"{compliant} compliant; {excluded} complex excluded; "
            f"{ambiguous} ambiguous; {finding_count} finding(s)",
        )

    if prose_findings or example_audit.findings:
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
