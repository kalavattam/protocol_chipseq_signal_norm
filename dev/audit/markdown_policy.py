#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: markdown_policy.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


"""
Check maintained Markdown against the proposed golden-source policy.
"""

from __future__ import annotations

import argparse
import dataclasses
import html
import json
import re
from pathlib import Path
from urllib.parse import unquote, unquote_to_bytes

from dev.tools.markdown_format import (
    ANCHOR_CANDIDATE,
    ATX,
    BLOCKQUOTE,
    CANONICAL_ANCHOR,
    FENCE,
    LIST,
    HeadingUnit,
    canonical_blockquote,
    fence_errors,
    fenced_indexes,
    format_document,
    format_table,
    heading_source_at,
    heading_units,
    informal_heading,
    maintained_markdown,
    parse_table,
    split_pipe_row,
)

STANDARD_RULE_ID = re.compile(r"`[A-Z][A-Z0-9]*(?:\.[A-Z0-9_]+)+`")
STANDARD_FIELDS = (
    "Classification",
    "Scope",
    "Automation",
    "Semantic remainder",
    "Exceptions",
)
STANDARD_FIELD = re.compile(
    r"^\*\*(?P<label>Classification|Scope|Automation|"
    r"Semantic remainder|Exceptions):\*\*[ \t]+\S.*$",
)
MARKDOWN_LINK = re.compile(
    r"(?<!!)\[[^\]]*\]\((?P<target>[^)\s]+)(?:\s+[^)]*)?\)",
)


@dataclasses.dataclass(frozen=True)
class Finding:
    """
    Represent one stable Markdown policy finding.
    """

    rule_id: str
    line: int
    message: str
    classification: str = "deterministic"


def _blank_count(lines: list[str], start: int) -> tuple[int, int]:
    index = start

    while index < len(lines) and not lines[index].strip():
        index += 1

    return index - start, index


def _fenced_lines(lines: list[str]) -> set[int]:
    """
    Return zero-based line indexes belonging to fenced code blocks.
    """

    return fenced_indexes(lines)


def check_standard_sections(text: str) -> list[Finding]:
    """
    Check canonical fields and IDs in maintained standards sections.

    Parameters
    ----------
    text : str
        Complete maintained standards document.

    Returns
    -------
    findings : list[Finding]
        Missing, duplicate, malformed, or misplaced owner-field findings.
    """

    lines = text.splitlines()
    literal = _fenced_lines(lines)
    findings: list[Finding] = []
    rule_sections: list[int] = []
    boundaries: list[int] = []

    for index, line in enumerate(lines):
        if index in literal:
            continue

        heading = ATX.match(line)
        if heading is None:
            continue

        level = len(heading.group("marks"))

        if level <= 2:
            boundaries.append(index)

        if level == 2:
            if STANDARD_RULE_ID.search(heading.group("text")):
                rule_sections.append(index)
            else:
                findings.append(
                    Finding(
                        "MD.STANDARD.SECTION",
                        index + 1,
                        "maintained standards H2 is missing a rule ID",
                    ),
                )

    for start in rule_sections:
        end = next(
            (value for value in boundaries if value > start),
            len(lines),
        )
        fields: list[tuple[str, int]] = []

        for index in range(start + 1, end):
            if index in literal:
                continue

            field = STANDARD_FIELD.match(lines[index])

            if field is not None:
                fields.append((field.group("label"), index))

        labels = [label for label, _ in fields]

        for required in STANDARD_FIELDS:
            positions = [index for label, index in fields if label == required]

            if not positions:
                findings.append(
                    Finding(
                        "MD.STANDARD.SECTION",
                        start + 1,
                        f"rule section is missing **{required}:**",
                    ),
                )
            elif len(positions) > 1:
                findings.append(
                    Finding(
                        "MD.STANDARD.SECTION",
                        positions[1] + 1,
                        f"rule section repeats **{required}:**",
                    ),
                )

        if (
            all(labels.count(required) == 1 for required in STANDARD_FIELDS)
            and tuple(labels) != STANDARD_FIELDS
        ):
            findings.append(
                Finding(
                    "MD.STANDARD.SECTION",
                    start + 1,
                    "rule fields must appear in canonical order",
                ),
            )

        adjacent = (
            start + 1 < end
            and start + 1 not in literal
            and (field := STANDARD_FIELD.match(lines[start + 1])) is not None
            and field.group("label") == "Classification"
        )

        if not adjacent:
            findings.append(
                Finding(
                    "MD.STANDARD.SECTION",
                    start + 1,
                    "**Classification:** must immediately follow the heading",
                ),
            )

    return findings


def _anchor_findings(
    lines: list[str],
    literal: set[int],
) -> list[Finding]:
    """
    Return canonical-source, attachment, and uniqueness findings.
    """

    findings: list[Finding] = []
    identifiers: dict[str, list[int]] = {}
    candidates = [
        index
        for index, line in enumerate(lines)
        if index not in literal and ANCHOR_CANDIDATE.fullmatch(line)
    ]

    for index in candidates:
        line = lines[index]
        match = CANONICAL_ANCHOR.fullmatch(line)

        if match is None:
            findings.append(
                Finding(
                    "MD.ANCHOR.CANONICAL",
                    index + 1,
                    'heading anchor must use exact <a id="ID"></a> source',
                ),
            )
        else:
            identifier = html.unescape(match.group("id"))

            if not identifier or any(
                character in " \t\r\n\f\v" for character in identifier
            ):
                findings.append(
                    Finding(
                        "MD.ANCHOR.CANONICAL",
                        index + 1,
                        (
                            "decoded anchor ID must be nonempty and contain "
                            "no ASCII whitespace"
                        ),
                    ),
                )

            identifiers.setdefault(identifier, []).append(index)

        following = index + 1
        attached = (
            following < len(lines)
            and following not in literal
            and heading_source_at(lines, following, literal) is not None
        )

        if not attached:
            findings.append(
                Finding(
                    "MD.ANCHOR.CANONICAL",
                    index + 1,
                    (
                        "heading anchor must be immediately followed by one "
                        "heading"
                    ),
                ),
            )

        if following < len(lines) and ANCHOR_CANDIDATE.fullmatch(
            lines[following],
        ):
            findings.append(
                Finding(
                    "MD.ANCHOR.CANONICAL",
                    following + 1,
                    "consecutive heading anchors are prohibited",
                ),
            )

    for identifier, indexes in identifiers.items():
        for index in indexes[1:]:
            findings.append(
                Finding(
                    "MD.ANCHOR.CANONICAL",
                    index + 1,
                    f"decoded anchor ID is duplicated: {identifier}",
                ),
            )

    return findings


def _boundary_tokens(
    lines: list[str],
    start: int,
    end: int,
    literal: set[int],
) -> list[str]:
    """
    Return boundary source with anchor candidates erased.
    """

    tokens: list[str] = []

    for index in range(start, end):
        if (
            index in literal
            or ANCHOR_CANDIDATE.fullmatch(lines[index]) is not None
        ):
            continue

        token = lines[index]
        if token == "" and tokens and tokens[-1] == "":
            continue

        tokens.append(token)

    return tokens


def _last_real_index(
    lines: list[str],
    start: int,
    end: int,
    literal: set[int],
) -> int | None:
    for index in range(end - 1, start - 1, -1):
        if index in literal:
            return index

        line = lines[index]
        if (
            line.strip()
            and line != "<br />"
            and ANCHOR_CANDIDATE.fullmatch(line) is None
        ):
            return index

    return None


def _section_findings(
    lines: list[str],
    literal: set[int],
    units: list[HeadingUnit],
) -> list[Finding]:
    """
    Return one deterministic finding for each malformed boundary.
    """

    findings: list[Finding] = []

    for position, unit in enumerate(units[1:], 1):
        previous = units[position - 1]
        segment_start = previous.heading + 1
        segment_end = unit.start
        last_real = _last_real_index(
            lines,
            segment_start,
            segment_end,
            literal,
        )
        direct_child = unit.level == previous.level + 1 and last_real is None
        first_boundary = (
            _boundary_tokens(
                lines,
                last_real + 1,
                segment_end,
                literal,
            )[:1]
            if last_real is not None
            else []
        )

        if direct_child:
            actual = _boundary_tokens(
                lines,
                segment_start,
                segment_end,
                literal,
            )
            valid = actual == []
            message = (
                "contentless direct parent and child must have no separator "
                "or blank"
            )
        elif (
            last_real is not None
            and lines[last_real].lower() == "</details>"
            and first_boundary == ["<br />"]
        ):
            actual = _boundary_tokens(
                lines,
                last_real + 1,
                segment_end,
                literal,
            )
            valid = actual == ["<br />", ""]
            message = (
                "details-close boundary must use its direct break with no "
                "second separator"
            )
        else:
            boundary_start = (
                last_real + 1 if last_real is not None else segment_start
            )
            actual = _boundary_tokens(
                lines,
                boundary_start,
                segment_end,
                literal,
            )
            valid = actual == ["", "<br />", ""]
            message = (
                "ordinary headed section needs blank, <br />, blank boundary"
            )

        if not valid:
            findings.append(
                Finding(
                    "MD.SECTION.BREAK",
                    unit.title + 1,
                    message,
                ),
            )

    return findings


def _heading_spacing_findings(
    lines: list[str],
    units: list[HeadingUnit],
) -> list[Finding]:
    """
    Return spacing findings for formal and recognized informal headings.
    """

    findings: list[Finding] = []

    for unit in units:
        if unit.heading + 1 >= len(lines):
            continue

        blanks, following = _blank_count(lines, unit.heading + 1)
        section_break_follows = (
            following < len(lines) and lines[following] == "<br />"
        )
        expected = 1 if section_break_follows else 0

        if following < len(lines) and blanks != expected:
            findings.append(
                Finding(
                    "MD.HEADING.SPACING",
                    unit.title + 1,
                    (
                        "heading needs one blank before an ordinary section "
                        "break and none before its first content block"
                    ),
                ),
            )

    return findings


def _ambiguous_informal_findings(
    lines: list[str],
    literal: set[int],
    units: list[HeadingUnit],
) -> list[Finding]:
    """
    Inventory emphasized-only lines not safely recognized as headings.
    """

    recognized = {unit.heading for unit in units if unit.informal}
    return [
        Finding(
            "MD.HEADING.INFORMAL",
            index + 1,
            "review whether standalone emphasis was intended as H7/H8",
            "advisory",
        )
        for index, line in enumerate(lines)
        if index not in literal
        and index not in recognized
        and informal_heading(line) is not None
    ]


def explicit_anchor_ids(text: str) -> dict[str, list[int]]:
    """
    Return decoded canonical anchor IDs and their source lines.
    """

    lines = text.splitlines()
    literal = fenced_indexes(lines)
    identifiers: dict[str, list[int]] = {}

    for index, line in enumerate(lines):
        if index in literal:
            continue

        match = CANONICAL_ANCHOR.fullmatch(line)
        if match is None:
            continue

        identifiers.setdefault(
            html.unescape(match.group("id")),
            [],
        ).append(index + 1)

    return identifiers


def _decode_fragment(fragment: str) -> str | None:
    try:
        return unquote_to_bytes(fragment).decode("utf-8", errors="strict")
    except (UnicodeDecodeError, ValueError):
        return None


def check_explicit_links(
    texts: dict[Path, str],
) -> dict[Path, list[Finding]]:
    """
    Check only links whose fragment names a present explicit ID.
    """

    resolved_texts = {
        file_path.resolve(): text for file_path, text in texts.items()
    }
    identifiers = {
        file_path: explicit_anchor_ids(text)
        for file_path, text in resolved_texts.items()
    }
    findings = {file_path: [] for file_path in resolved_texts}

    for source, text in resolved_texts.items():
        lines = text.splitlines()
        literal = fenced_indexes(lines)

        for index, line in enumerate(lines):
            if index in literal:
                continue

            for match in MARKDOWN_LINK.finditer(line):
                target = match.group("target")
                if (
                    target.startswith(
                        ("http://", "https://", "mailto:", "ftp://"),
                    )
                    or "#" not in target
                ):
                    continue

                relative, fragment = target.split("#", 1)
                decoded = _decode_fragment(fragment)
                if decoded is None:
                    continue

                target_path = (
                    source
                    if not relative
                    else (source.parent / unquote(relative)).resolve()
                )
                target_ids = identifiers.get(target_path)
                if not target_ids or decoded not in target_ids:
                    continue

                if len(target_ids[decoded]) != 1:
                    findings[source].append(
                        Finding(
                            "MD.ANCHOR.CANONICAL",
                            index + 1,
                            f"explicit-ID link is ambiguous: {target}",
                        ),
                    )

    return findings


def check_text(text: str) -> list[Finding]:
    """
    Return deterministic and advisory findings for one Markdown source.

    Parameters
    ----------
    text : str
        Markdown source to inspect.

    Returns
    -------
    findings : list[Finding]
        Deterministic and advisory policy findings.
    """

    lines = text.splitlines()
    literal = fenced_indexes(lines)

    units = heading_units(lines)
    findings: list[Finding] = []

    if not text.startswith("\n") or text.startswith("\n\n"):
        findings.append(
            Finding(
                "MD.FILE.BOUNDARY",
                1,
                "file must begin with exactly one line feed",
            ),
        )
    elif units:
        first_nonempty = next(
            (index for index, line in enumerate(lines) if line.strip()),
            None,
        )

        if first_nonempty != units[0].start:
            findings.append(
                Finding(
                    "MD.FILE.BOUNDARY",
                    1,
                    (
                        "first heading unit must immediately follow the "
                        "leading line feed"
                    ),
                ),
            )

    if not text.endswith("\n") or text.endswith("\n\n"):
        findings.append(
            Finding(
                "MD.FILE.BOUNDARY",
                max(1, len(lines)),
                "file must end with one terminal line feed and no empty line",
            ),
        )

    findings.extend(
        Finding("MD.FENCE.CLOSED", line, message)
        for line, message in fence_errors(text)
    )
    findings.extend(_anchor_findings(lines, literal))
    findings.extend(_heading_spacing_findings(lines, units))
    findings.extend(_section_findings(lines, literal, units))
    findings.extend(_ambiguous_informal_findings(lines, literal, units))

    index = 0

    while index < len(lines):
        line = lines[index]

        if index in literal:
            index += 1

            continue

        if "\t" in line[: len(line) - len(line.lstrip(" \t"))]:
            findings.append(
                Finding(
                    "MD.INDENT.SPACES",
                    index + 1,
                    "structural indentation contains a tab",
                ),
            )

        canonical_quote = canonical_blockquote(line)
        if canonical_quote is not None and canonical_quote != line:
            findings.append(
                Finding(
                    "MD.BLOCKQUOTE.MARKER",
                    index + 1,
                    (
                        "blockquote text and nested markers require one "
                        "space; an empty quoted separator is exactly '>'"
                    ),
                ),
            )

        if line.rstrip().endswith(":"):
            blanks, following = _blank_count(lines, index + 1)

            if following < len(lines):
                list_or_fence = bool(
                    FENCE.match(lines[following])
                    or LIST.match(lines[following])
                    or BLOCKQUOTE.match(lines[following])
                )
                table_follows = False

                if following + 1 < len(lines):
                    table = parse_table(lines[following : following + 2])
                    table_follows = table is not None

                expected = 1 if table_follows else 0

                if (list_or_fence or table_follows) and blanks != expected:
                    findings.append(
                        Finding(
                            "MD.COLON.STRUCTURE",
                            index + 1,
                            (
                                "colon needs one blank before a table and "
                                "none before a list, fence, or blockquote"
                            ),
                        ),
                    )

        list_match = LIST.match(line)

        if list_match is not None:
            indent = list_match.group("indent")

            if len(indent) % 2:
                findings.append(
                    Finding(
                        "MD.LIST.INDENT",
                        index + 1,
                        "list indentation must use two-space levels",
                    ),
                )

            blanks, following = _blank_count(lines, index + 1)

            if (
                following < len(lines)
                and LIST.match(lines[following])
                and blanks != 0
            ):
                findings.append(
                    Finding(
                        "MD.LIST.SPACING",
                        following + 1,
                        "ordinary list peers and nested lists remain tight",
                    ),
                )

        if index + 1 < len(lines) and split_pipe_row(line) is not None:
            end = index + 2

            while end < len(lines) and split_pipe_row(lines[end]) is not None:
                end += 1

            table = parse_table(lines[index:end])

            if table is not None:
                actual = "\n".join(lines[index:end])

                if actual != format_table(table):
                    findings.append(
                        Finding(
                            "MD.TABLE.CANONICAL",
                            index + 1,
                            "GFM table is not in canonical source form",
                        ),
                    )

                index = end

                continue

        index += 1

    if format_document(text) != text:
        findings.append(
            Finding(
                "MD.PROSE.UNWRAP",
                1,
                "explicit formatter would change prose or block spacing",
                "advisory",
            ),
        )

    return findings


def main(argv: list[str] | None = None) -> int:
    """
    Check explicit paths or the maintained repository Markdown scope.

    Parameters
    ----------
    argv : list[str] | None
        Explicit arguments, or None to read the process arguments.

    Returns
    -------
    status : int
        Zero when deterministic Markdown checks pass and one otherwise.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="*", type=Path)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--json", action="store_true")
    parser.add_argument(
        "--fail-on",
        choices=("deterministic", "all", "none"),
        default="deterministic",
        help="select which finding class determines the exit status",
    )
    args = parser.parse_args(argv)
    root = args.root.resolve()
    paths = args.paths or maintained_markdown(root)

    texts: dict[Path, str] = {}

    for file_path in paths:
        absolute = file_path if file_path.is_absolute() else root / file_path
        if (
            absolute.parent == root / "docs/standards"
            and absolute.name.startswith("bak.")
        ):
            continue

        texts[absolute.resolve()] = absolute.read_text(encoding="utf-8")

    link_findings = check_explicit_links(texts)
    output: list[dict[str, object]] = []

    for absolute, text in texts.items():
        relative = absolute.relative_to(root).as_posix()
        findings = check_text(text)
        findings.extend(link_findings[absolute])

        if (
            absolute.parent == root / "docs/standards"
            and not absolute.name.startswith("bak.")
        ):
            findings.extend(check_standard_sections(text))

        for finding in findings:
            output.append({"path": relative, **dataclasses.asdict(finding)})

    if args.json:
        print(json.dumps(output, indent=2, sort_keys=True))
    else:
        for finding in output:
            print(
                f"{finding['rule_id']}: {finding['path']}:"
                f"{finding['line']}: {finding['message']}",
            )

    if args.fail_on == "none":
        return 0

    if args.fail_on == "all":
        return 1 if output else 0

    return (
        1
        if any(item["classification"] == "deterministic" for item in output)
        else 0
    )


if __name__ == "__main__":
    raise SystemExit(main())
