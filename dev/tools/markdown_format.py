#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: markdown_format.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


"""Pure Markdown parsing and explicit canonical-formatting primitives."""

from __future__ import annotations

import argparse
import csv
import html
import io
import re
import subprocess
import unicodedata
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path

ATX = re.compile(r"^(?P<indent> {0,3})(?P<marks>#{1,6})[ \t]+(?P<text>.*)$")
SETEXT = re.compile(r"^ {0,3}(?P<marker>=+|-+)[ \t]*$")
FENCE = re.compile(r"^ {0,3}(?P<marker>`{3,}|~{3,})(?P<info>.*)$")
LIST = re.compile(r"^(?P<indent> *)(?:[-+*]|[0-9]{1,9}[.)])[ \t]+\S")
SEPARATOR = re.compile(r"^:?-{3,}:?$")
DETAILS_OPEN = re.compile(r"^\s*<details(?:\s[^>]*)?>\s*$", re.IGNORECASE)
DETAILS_CLOSE = re.compile(r"^\s*</details>\s*$", re.IGNORECASE)
SUMMARY = re.compile(
    r"^\s*<summary(?:\s[^>]*)?>.*?</summary>\s*$",
    re.IGNORECASE,
)
CANONICAL_ANCHOR = re.compile(
    r'^<a id="(?P<id>[^"\s]+)"></a>$'
)
ANCHOR_CANDIDATE = re.compile(
    r"^\s*<a\b[^>]*>\s*</a>\s*$",
    re.IGNORECASE,
)
INFORMAL_H7 = re.compile(
    r"^(?:\*\*(?P<md>.+)\*\*|"
    r"<(?:strong|b)>(?P<html>.+)</(?:strong|b)>)$",
    re.IGNORECASE,
)
INFORMAL_H8 = re.compile(
    r"^(?:\*(?P<md>[^*].*?)\*|"
    r"<(?:em|i)>(?P<html>.+)</(?:em|i)>)$",
    re.IGNORECASE,
)


@dataclass(frozen=True)
class ParsedTable:
    """Represent one validated GFM pipe table."""

    header: tuple[str, ...]
    alignments: tuple[str, ...]
    body: tuple[tuple[str, ...], ...]
    indent: str = ""


@dataclass(frozen=True)
class HeadingUnit:
    """Represent one formal or safely recognized informal heading unit."""

    start: int
    title: int
    heading: int
    level: int
    text: str
    anchor_candidates: tuple[int, ...] = ()
    informal: bool = False


def _escaped(text: str, index: int) -> bool:
    count = 0
    index -= 1
    while index >= 0 and text[index] == "\\":
        count += 1
        index -= 1
    return count % 2 == 1


def split_pipe_row(line: str) -> tuple[str, ...] | None:
    """Split one GFM table row outside escapes and code spans."""

    text = line.strip()
    cells: list[str] = []
    current: list[str] = []
    code_run: int | None = None
    index = 0
    structural = False
    while index < len(text):
        char = text[index]
        if char == "`" and not _escaped(text, index):
            end = index + 1
            while end < len(text) and text[end] == "`":
                end += 1
            length = end - index
            current.append(text[index:end])
            if code_run is None:
                code_run = length
            elif code_run == length:
                code_run = None
            index = end
            continue
        if char == "|" and code_run is None and not _escaped(text, index):
            cells.append("".join(current))
            current = []
            structural = True
        else:
            current.append(char)
        index += 1
    cells.append("".join(current))
    if not structural:
        return None
    if cells and not cells[0].strip():
        cells = cells[1:]
    if cells and not cells[-1].strip():
        cells = cells[:-1]
    if len(cells) < 2:
        return None
    return tuple(cell.strip() for cell in cells)


def parse_table(lines: Sequence[str]) -> ParsedTable | None:
    """Parse a complete GFM pipe table with equal-width logical rows."""

    if len(lines) < 2 or any(not line.strip() for line in lines):
        return None
    rows = [split_pipe_row(line) for line in lines]
    if any(row is None for row in rows):
        return None
    parsed = [row for row in rows if row is not None]
    width = len(parsed[0])
    if any(len(row) != width for row in parsed):
        return None
    alignments: list[str] = []
    for cell in parsed[1]:
        compact = re.sub(r"\s+", "", cell)
        if not SEPARATOR.fullmatch(compact):
            return None
        left = compact.startswith(":")
        right = compact.endswith(":")
        alignments.append("center" if left and right else "right" if right else "left")
    indent = lines[0][: len(lines[0]) - len(lines[0].lstrip(" \t"))]
    return ParsedTable(
        parsed[0], tuple(alignments), tuple(parsed[2:]), indent
    )


def visual_width(text: str) -> int:
    """Return a deterministic approximation of editor display width."""

    width = 0
    for char in text:
        value = ord(char)
        if value in {0x200B, 0x200C, 0x200D, 0xFE0E, 0xFE0F}:
            continue
        if unicodedata.combining(char):
            continue
        width += 2 if unicodedata.east_asian_width(char) in {"W", "F"} else 1
    return width


def _separator(alignment: str) -> str:
    return {"left": ":---", "right": "---:", "center": ":---:"}[alignment]


def format_table(table: ParsedTable) -> str:
    """Render one table with leading/trailing pipes and aligned source."""

    ordinary = (table.header, *table.body)
    widths = [
        max(
            4,
            visual_width(_separator(table.alignments[column])),
            *(visual_width(row[column]) for row in ordinary),
        )
        for column in range(len(table.header))
    ]

    def pad(text: str, column: int) -> str:
        return text + " " * (widths[column] - visual_width(text))

    def row(values: Sequence[str]) -> str:
        cells = [pad(value, index) for index, value in enumerate(values)]
        return table.indent + "| " + " | ".join(cells) + " |"

    output = [row(table.header)]
    output.append(row([_separator(value) for value in table.alignments]))
    output.extend(row(values) for values in table.body)
    return "\n".join(output)


def convert_delimited(text: str, delimiter: str) -> str:
    """Convert explicit CSV or TSV input to a canonical left-aligned table."""

    if delimiter not in {",", "\t"}:
        raise ValueError("delimiter must be ',' or a tab")
    rows = list(csv.reader(io.StringIO(text), delimiter=delimiter))
    if not rows or not rows[0]:
        raise ValueError("delimited input is empty")
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError("delimited rows must have equal column counts")
    table = ParsedTable(
        tuple(cell.strip() for cell in rows[0]),
        tuple("left" for _ in range(width)),
        tuple(tuple(cell.strip() for cell in row) for row in rows[1:]),
    )
    return format_table(table)


def fence_errors(text: str) -> list[tuple[int, str]]:
    """Return unclosed or mismatched fenced-code diagnostics."""

    opened: tuple[str, int, int] | None = None
    for number, line in enumerate(text.splitlines(), 1):
        match = FENCE.match(line)
        if match is None:
            continue
        marker = match.group("marker")
        if opened is None:
            opened = (marker[0], len(marker), number)
        elif marker[0] == opened[0] and len(marker) >= opened[1]:
            opened = None
    if opened is None:
        return []
    return [(opened[2], "fenced code block is not closed")]


def _plain_prose(line: str) -> bool:
    stripped = line.strip()
    if not stripped or ATX.match(line) or LIST.match(line) or FENCE.match(line):
        return False
    if stripped.startswith((">", "<", "|", "```", "~~~", "[")):
        return False
    return not line.startswith(("    ", "\t"))


def _format_tables(lines: list[str]) -> list[str]:
    output: list[str] = []
    index = 0
    opened: tuple[str, int] | None = None
    while index < len(lines):
        fence = FENCE.match(lines[index])
        if fence is not None:
            marker = fence.group("marker")
            if opened is None:
                opened = (marker[0], len(marker))
            elif marker[0] == opened[0] and len(marker) >= opened[1]:
                opened = None
            output.append(lines[index])
            index += 1
            continue
        if opened is not None:
            output.append(lines[index])
            index += 1
            continue
        if index + 1 >= len(lines) or split_pipe_row(lines[index]) is None:
            output.append(lines[index])
            index += 1
            continue
        end = index + 2
        while end < len(lines) and split_pipe_row(lines[end]) is not None:
            end += 1
        parsed = parse_table(lines[index:end])
        if parsed is None:
            output.append(lines[index])
            index += 1
            continue
        output.extend(format_table(parsed).splitlines())
        index = end
    return output


def fenced_indexes(lines: Sequence[str]) -> set[int]:
    """Return indexes whose source must remain literal."""

    literal: set[int] = set()
    opened: tuple[str, int] | None = None
    for index, line in enumerate(lines):
        fence = FENCE.match(line)
        if opened is None:
            if fence is not None:
                marker = fence.group("marker")
                opened = (marker[0], len(marker))
                literal.add(index)
            continue
        literal.add(index)
        if fence is None:
            continue
        marker = fence.group("marker")
        if marker[0] == opened[0] and len(marker) >= opened[1]:
            opened = None
    return literal


def informal_heading(line: str) -> tuple[int, str] | None:
    """Return the informal heading level and visible source text."""

    if match := INFORMAL_H7.fullmatch(line):
        return 7, (match.group("md") or match.group("html")).strip()
    if match := INFORMAL_H8.fullmatch(line):
        return 8, (match.group("md") or match.group("html")).strip()
    return None


def heading_source_at(
    lines: Sequence[str],
    index: int,
    literal: set[int],
) -> tuple[int, str, int, bool] | None:
    """Return level, text, final line, and informal status at a title line."""

    if index in literal:
        return None
    if formal := ATX.fullmatch(lines[index]):
        return (
            len(formal.group("marks")),
            formal.group("text").strip(),
            index,
            False,
        )
    if (
        index + 1 < len(lines)
        and index + 1 not in literal
        and lines[index].strip()
        and ANCHOR_CANDIDATE.fullmatch(lines[index]) is None
        and (underline := SETEXT.fullmatch(lines[index + 1])) is not None
    ):
        return (
            1 if underline.group("marker").startswith("=") else 2,
            lines[index].strip(),
            index + 1,
            False,
        )
    if informal := informal_heading(lines[index]):
        return informal[0], informal[1], index, True
    return None


def _attached_anchor_candidates(
    lines: Sequence[str],
    heading: int,
    literal: set[int],
) -> tuple[int, ...]:
    indexes: list[int] = []
    index = heading - 1
    while (
        index >= 0
        and index not in literal
        and ANCHOR_CANDIDATE.fullmatch(lines[index])
    ):
        indexes.append(index)
        index -= 1
    return tuple(reversed(indexes))


def _structural_tokens_only(
    lines: Sequence[str],
    start: int,
    end: int,
    literal: set[int],
) -> bool:
    return all(
        index not in literal
        and (
            not lines[index].strip()
            or lines[index] == "<br />"
            or ANCHOR_CANDIDATE.fullmatch(lines[index])
        )
        for index in range(start, end)
    )


def _trailing_separator_signal(
    lines: Sequence[str],
    start: int,
    end: int,
    literal: set[int],
) -> bool:
    """Return whether the trailing structural cluster contains a break."""

    trailing = start
    for index in range(start, end):
        if index in literal:
            trailing = index + 1
            continue
        line = lines[index]
        if (
            line.strip()
            and line != "<br />"
            and ANCHOR_CANDIDATE.fullmatch(line) is None
        ):
            trailing = index + 1
    return any(
        index not in literal and lines[index] == "<br />"
        for index in range(trailing, end)
    )


def _inside_details(
    lines: Sequence[str],
    end: int,
    literal: set[int],
) -> bool:
    """Return whether an index is inside an unclosed details element."""

    depth = 0
    for index in range(end):
        if index in literal:
            continue
        if DETAILS_OPEN.fullmatch(lines[index]):
            depth += 1
        elif DETAILS_CLOSE.fullmatch(lines[index]) and depth:
            depth -= 1
    return depth > 0


def heading_units(lines: Sequence[str]) -> list[HeadingUnit]:
    """Return formal and safely recognized informal heading units."""

    literal = fenced_indexes(lines)
    units: list[HeadingUnit] = []
    for index in range(len(lines)):
        source = heading_source_at(lines, index, literal)
        if source is None:
            continue
        level, text, final_line, is_informal = source
        anchors = _attached_anchor_candidates(lines, index, literal)
        if is_informal:
            first_entry = not units and all(
                not value.strip()
                or ANCHOR_CANDIDATE.fullmatch(value)
                for value in lines[:index]
            )
            separated = (
                not _inside_details(lines, index, literal)
                and _trailing_separator_signal(
                    lines,
                    units[-1].heading + 1 if units else 0,
                    anchors[0] if anchors else index,
                    literal,
                )
            )
            direct_parent = bool(
                units
                and units[-1].level + 1 == level
                and _structural_tokens_only(
                    lines,
                    units[-1].heading + 1,
                    anchors[0] if anchors else index,
                    literal,
                )
            )
            if not (anchors or first_entry or separated or direct_parent):
                continue
        units.append(
            HeadingUnit(
                start=anchors[0] if anchors else index,
                title=index,
                heading=final_line,
                level=level,
                text=text,
                anchor_candidates=anchors,
                informal=is_informal,
            )
        )
    return units


def maintained_markdown(root: Path) -> list[Path]:
    """Return tracked and nonignored untracked maintained Markdown files."""

    commands = (
        ["git", "ls-files", "-z", "--", "*.md"],
        [
            "git",
            "ls-files",
            "-z",
            "--others",
            "--exclude-standard",
            "--",
            "*.md",
        ],
    )
    paths: set[Path] = set()
    for command in commands:
        result = subprocess.run(
            command,
            cwd=root,
            check=True,
            capture_output=True,
        )
        for value in result.stdout.split(b"\0"):
            if value:
                paths.add(root / value.decode("utf-8"))
    fixture_root = root / "tests/fixtures/markdown"
    excluded = {
        fixture_root / "format/input.md",
        *(fixture_root / "rejected").glob("*.md"),
    }
    return sorted(
        file_path
        for file_path in paths
        if file_path.is_file()
        and file_path not in excluded
        and not (
            file_path.parent == root / "docs/standards"
            and file_path.name.startswith("bak.")
        )
    )


def _table_starts(lines: Sequence[str], index: int) -> bool:
    """Return whether a supported GFM table starts at 'index'."""

    return (
        index + 1 < len(lines)
        and parse_table(lines[index : index + 2]) is not None
    )


def _canonicalize_deterministic_gaps(lines: list[str]) -> list[str]:
    """Canonicalize only checker-owned simple block gaps."""

    literal = fenced_indexes(lines)
    heading_ends = {
        unit.heading for unit in heading_units(lines)
    }
    content = [index for index, line in enumerate(lines) if line.strip()]
    output: list[str] = []
    for position, previous_index in enumerate(content):
        output.append(lines[previous_index])
        if position + 1 == len(content):
            continue
        following_index = content[position + 1]
        previous = lines[previous_index]
        following = lines[following_index]
        desired: int | None = None
        if (
            previous_index not in literal
            and (
                following_index not in literal
                or FENCE.match(following) is not None
            )
            and not any(
                item in literal
                for item in range(previous_index + 1, following_index)
            )
        ):
            table_follows = _table_starts(lines, following_index)
            if previous_index in heading_ends:
                desired = (
                    1
                    if table_follows or following == "<br />"
                    else 0
                )
            elif previous.rstrip().endswith(":") and (
                LIST.match(following)
                or FENCE.match(following)
                or table_follows
            ):
                desired = 1 if table_follows else 0
            elif LIST.match(previous) and LIST.match(following):
                desired = 0

        count = (
            following_index - previous_index - 1
            if desired is None
            else desired
        )
        output.extend("" for _ in range(count))
    return output


def _boundary_cluster_start(
    lines: Sequence[str],
    heading: int,
    literal: set[int],
) -> int:
    """Return the first blank, break, or anchor before one heading."""

    index = heading
    while index > 0:
        previous = index - 1
        if previous in literal:
            break
        line = lines[previous]
        if (
            line.strip()
            and line != "<br />"
            and ANCHOR_CANDIDATE.fullmatch(line) is None
        ):
            break
        index = previous
    return index


def _canonicalize_section_breaks(lines: list[str]) -> list[str]:
    """Canonicalize recognized boundary clusters without guessing IDs."""

    literal = fenced_indexes(lines)
    units = heading_units(lines)
    output = list(lines)
    for position in range(len(units) - 1, -1, -1):
        unit = units[position]
        start = _boundary_cluster_start(lines, unit.title, literal)
        candidate_indexes = [
            index
            for index in range(start, unit.title)
            if ANCHOR_CANDIDATE.fullmatch(lines[index])
        ]
        anchor_line: str | None = None
        if candidate_indexes:
            matches = [
                CANONICAL_ANCHOR.fullmatch(lines[index])
                for index in candidate_indexes
            ]
            if any(match is None for match in matches):
                continue
            identifiers = {
                html.unescape(match.group("id"))
                for match in matches
                if match is not None
            }
            if len(identifiers) != 1:
                continue
            anchor_line = lines[candidate_indexes[-1]]

        boundary: list[str]
        if position == 0:
            boundary = []
        else:
            previous = units[position - 1]
            direct_child = (
                unit.level == previous.level + 1
                and _structural_tokens_only(
                    lines,
                    previous.heading + 1,
                    start,
                    literal,
                )
            )
            details_close = (
                start > 0
                and DETAILS_CLOSE.fullmatch(lines[start - 1]) is not None
                and start < unit.heading
                and lines[start] == "<br />"
            )
            if direct_child:
                boundary = []
            elif details_close:
                boundary = ["<br />", ""]
            else:
                boundary = ["", "<br />", ""]
        if anchor_line is not None:
            boundary.append(anchor_line)
        output[start : unit.title] = boundary
    return output


def format_deterministic(text: str) -> str:
    """Apply only approved, checker-recognized deterministic rewrites."""

    lines = _format_tables(text.splitlines())
    while lines and not lines[0].strip():
        lines.pop(0)
    while lines and not lines[-1].strip():
        lines.pop()
    lines = _canonicalize_deterministic_gaps(lines)
    lines = _canonicalize_section_breaks(lines)
    return "\n" + "\n".join(lines) + "\n"


def _rebase_details(lines: list[str]) -> list[str]:
    output = list(lines)
    parent_level: int | None = None
    index = 0
    while index < len(output):
        heading = ATX.match(output[index])
        if heading:
            parent_level = len(heading.group("marks"))
            index += 1
            continue
        if not DETAILS_OPEN.match(output[index]) or parent_level is None:
            index += 1
            continue
        close = next(
            (
                item
                for item in range(index + 1, len(output))
                if DETAILS_CLOSE.match(output[item])
            ),
            None,
        )
        if close is None or not any(
            SUMMARY.match(output[item]) for item in range(index + 1, close)
        ):
            index += 1
            continue
        headings = [
            item
            for item in range(index + 1, close)
            if ATX.match(output[item])
        ]
        if headings:
            shallow = min(len(ATX.match(output[item]).group("marks")) for item in headings)  # type: ignore[union-attr]
            shift = parent_level + 1 - shallow
            if all(
                1 <= len(ATX.match(output[item]).group("marks")) + shift <= 6  # type: ignore[union-attr]
                for item in headings
            ):
                for item in headings:
                    match = ATX.match(output[item])
                    assert match is not None
                    level = len(match.group("marks")) + shift
                    output[item] = "#" * level + " " + match.group("text").strip()
        index = close + 1
    return output


def format_document(text: str) -> str:
    """Return the proposed canonical form without writing a file."""

    lines = _rebase_details(_format_tables(text.splitlines()))
    output: list[str] = []
    in_fence = False
    index = 0
    while index < len(lines):
        line = lines[index]
        fence = FENCE.match(line)
        if fence:
            in_fence = not in_fence
            output.append(line)
            index += 1
            continue
        if in_fence:
            output.append(line)
            index += 1
            continue
        if _plain_prose(line):
            paragraph = [line.strip()]
            while index + 1 < len(lines) and _plain_prose(lines[index + 1]):
                if line.endswith(("  ", "\\")):
                    break
                index += 1
                line = lines[index]
                paragraph.append(line.strip())
            line = " ".join(paragraph)
        output.append(line)
        index += 1

    compact: list[str] = []
    in_fence = False
    for line in output:
        if FENCE.match(line):
            in_fence = not in_fence
            compact.append(line)
            continue
        if not in_fence and not line and compact and not compact[-1]:
            continue
        compact.append(line)

    spaced: list[str] = []
    in_fence = False
    for index, line in enumerate(compact):
        if FENCE.match(line):
            in_fence = not in_fence
            spaced.append(line)
            continue
        if in_fence or line:
            spaced.append(line)
            continue
        previous = compact[index - 1] if index else ""
        following = compact[index + 1] if index + 1 < len(compact) else ""
        following_table = (
            index + 2 < len(compact)
            and parse_table(compact[index + 1 : index + 3]) is not None
        )
        if ATX.match(previous) and not following_table:
            continue
        if previous.rstrip().endswith(":") and (
            FENCE.match(following) or LIST.match(following)
        ):
            continue
        if LIST.match(previous) and LIST.match(following):
            continue
        spaced.append(line)

    while spaced and not spaced[0]:
        spaced.pop(0)
    while spaced and not spaced[-1]:
        spaced.pop()
    return "\n" + "\n".join(spaced) + "\n"


def main(argv: list[str] | None = None) -> int:
    """Format selected files only when explicitly invoked with '--write'."""

    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="*", type=Path)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--write", action="store_true")
    parser.add_argument(
        "--mode",
        choices=("deterministic", "proposed"),
        default="deterministic",
        help=(
            "use the approved deterministic subset (default) or preview "
            "deferred proposal behaviors"
        ),
    )
    args = parser.parse_args(argv)
    if args.mode == "proposed" and args.write:
        parser.error("--mode proposed is preview-only and cannot be written")
    formatter: Callable[[str], str] = (
        format_deterministic
        if args.mode == "deterministic"
        else format_document
    )
    root = args.root.resolve()
    paths = args.paths or maintained_markdown(root)
    changed = False
    for file_path in paths:
        absolute = (
            file_path
            if file_path.is_absolute()
            else root / file_path
        )
        original = absolute.read_text(encoding="utf-8")
        formatted = formatter(original)
        if formatted == original:
            continue
        changed = True
        if args.write:
            absolute.write_text(formatted, encoding="utf-8")
        else:
            try:
                print(absolute.relative_to(root))
            except ValueError:
                print(absolute)
    return 0 if args.write or not changed else 1


if __name__ == "__main__":
    raise SystemExit(main())
