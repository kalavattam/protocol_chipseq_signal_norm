#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: comment_wrap_format.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Refill ordinary Shell and Python comment prose greedily through column 79.
"""

from __future__ import annotations

import argparse
import io
import itertools
import re
import sys
import tokenize
from pathlib import Path

# These are the only copy of the comment-wrap rules. The Shell checker in
# 'dev/audit/shell_source_form.py' imports them, as
# 'dev/audit/markdown_policy.py' imports 'dev/tools/markdown_format.py', so the
# checker and the formatter cannot disagree. 'SOURCE.PROSE.WRAP' owns the rule.
WRAP_COLUMN = 79
HEADER_ROWS = 16
COMMENT_LINE = re.compile(r"^(?P<indent>[ \t]*)# (?P<pad> *)(?P<body>\S.*)$")
COMMENT_FENCE = re.compile(r"^[ \t]*#[ \t]*(?:'''|```|~~~)\S*[ \t]*$")
LIST_ITEM = re.compile(r"^(?:[-*]|\d+[.)]|[A-Z]{2,}:)\s")
INVENTORY_ROW = re.compile(r"^\S*[^\s.,;:!?](?:\s+#\S+)?$")
FORMULA_START = re.compile(r"^[(\d]\S*\s.*\s=\s")
UNIT_CLOSERS = {"(": ")", "[": "]", "{": "}", "'": "'", '"': '"', "`": "`"}
HEREDOC = re.compile(
    r"(?<!<)(?P<operator><<-?)(?!<)(?P<separator>[ \t]*)"
    r"(?:'(?P<single>[A-Za-z_][A-Za-z0-9_]*)'|"
    r'"(?P<double>[A-Za-z_][A-Za-z0-9_]*)"|'
    r"(?P<plain>[A-Za-z_][A-Za-z0-9_]*))",
)

# 'PY.COMMENT.FORM' leaves tool directives, coverage pragmas, and
# type-checker directives to their own syntax.
PYTHON_DIRECTIVE = re.compile(
    r"^#\s*(?:fmt|noqa|type|pragma|pylint|ruff|mypy|isort|pyright)\b",
)


def shell_candidates(lines: list[str]) -> set[int]:
    """
    Return the one-based Shell lines outside heredoc bodies.

    Parameters
    ----------
    lines : list[str]
        Physical source lines.

    Returns
    -------
    candidates : set[int]
        Line numbers whose comments may be ordinary prose.
    """

    candidates: set[int] = set()
    delimiter: str | None = None

    for number, line in enumerate(lines, 1):
        if delimiter is not None:
            if line.strip() == delimiter:
                delimiter = None

            continue

        candidates.add(number)
        heredoc = HEREDOC.search(line)

        if heredoc is not None:
            delimiter = next(
                value
                for name in ("single", "double", "plain")
                if (value := heredoc.group(name)) is not None
            )

    return candidates


def python_candidates(text: str) -> set[int]:
    """
    Return the one-based Python lines that hold a full-line comment.

    Python's tokenizer separates comments from string literals, so a line
    beginning with '#' inside a string is never a candidate.

    Parameters
    ----------
    text : str
        Python source.

    Returns
    -------
    candidates : set[int]
        Line numbers of full-line comments other than tool directives.
    """

    candidates: set[int] = set()
    lines = text.splitlines()

    for token in tokenize.generate_tokens(io.StringIO(text).readline):
        if token.type != tokenize.COMMENT:
            continue

        number, column = token.start
        line = lines[number - 1]

        if line[:column].strip() or PYTHON_DIRECTIVE.match(token.string):
            continue

        candidates.add(number)

    return candidates


def units(text: str) -> list[str]:
    """
    Split comment text into words, keeping bracketed and quoted spans whole.

    Parameters
    ----------
    text : str
        Comment prose with markers removed.

    Returns
    -------
    units : list[str]
        Words, with each multiword bracketed or quoted span joined into one.
    """

    words = text.split()
    result: list[str] = []
    index = 0

    while index < len(words):
        word = words[index]
        closer = UNIT_CLOSERS.get(word[0])
        end = index + 1

        if closer is not None and closer not in word[1:]:
            end = next(
                (
                    stop
                    for stop in range(index + 2, len(words) + 1)
                    if closer in words[stop - 1]
                ),
                index + 1,
            )

        result.append(" ".join(words[index:end]))
        index = end

    return result


def _introduced_by_list(lines: list[str], number: int, indent: str) -> bool:
    """
    Report whether an indented comment line continues a list item.

    Parameters
    ----------
    lines : list[str]
        Physical source lines.
    number : int
        One-based number of an indented comment line.
    indent : str
        Leading whitespace before the comment marker.

    Returns
    -------
    continues : bool
        Whether the nearest unindented comment line above, in the same comment
        block, is a list item.
    """

    for line in reversed(lines[: number - 1]):
        comment = COMMENT_LINE.match(line)

        if comment is None or comment["indent"] != indent:
            return False

        if not comment["pad"]:
            return bool(LIST_ITEM.match(comment["body"]))

    return False


def _follows_list_item(lines: list[str], number: int, indent: str) -> bool:
    """
    Report whether a list item precedes a line in the same comment block.

    Parameters
    ----------
    lines : list[str]
        Physical source lines.
    number : int
        One-based number of an unindented comment line.
    indent : str
        Leading whitespace before the comment marker.

    Returns
    -------
    follows : bool
        Whether a list item appears above the line before the comment block
        ends.
    """

    for line in reversed(lines[: number - 1]):
        comment = COMMENT_LINE.match(line)

        if comment is None or comment["indent"] != indent:
            return False

        if not comment["pad"] and LIST_ITEM.match(comment["body"]):
            return True

    return False


def _inventory_row(lines: list[str], number: int) -> bool:
    """
    Report whether a comment line is a row of a one-identifier inventory.

    A row holds one identifier, optionally followed by a '#' tag, and a
    neighboring line in the same comment block is another such row, as in a
    function inventory. A lone one-word line inside prose is prose.

    Parameters
    ----------
    lines : list[str]
        Physical source lines.
    number : int
        One-based line number.

    Returns
    -------
    row : bool
        Whether the line is an inventory row.
    """

    def single(index: int) -> re.Match[str] | None:
        if not 1 <= index <= len(lines):
            return None

        comment = COMMENT_LINE.match(lines[index - 1])

        if comment is None or comment["pad"]:
            return None

        return INVENTORY_ROW.match(comment["body"])

    return bool(
        single(number) and (single(number - 1) or single(number + 1)),
    )


def premature_breaks(
    lines: list[str],
    candidates: set[int],
) -> list[tuple[int, str]]:
    """
    Return comment lines that break before a unit that still fits.

    Parameters
    ----------
    lines : list[str]
        Physical source lines.
    candidates : set[int]
        One-based line numbers whose comments may be ordinary prose.

    Returns
    -------
    breaks : list[tuple[int, str]]
        One-based line number of each premature break and the unit that fits.
    """

    breaks: list[tuple[int, str]] = []
    fenced = False

    for number, (line, following_line) in enumerate(
        itertools.pairwise(lines),
        1,
    ):
        if COMMENT_FENCE.match(line):
            fenced = not fenced
            continue

        if (
            fenced
            or number <= HEADER_ROWS
            or number not in candidates
            or number + 1 not in candidates
            or "shellcheck" in line
            or "shellcheck" in following_line
            or COMMENT_FENCE.match(following_line)
        ):
            continue

        current = COMMENT_LINE.match(line)
        following = COMMENT_LINE.match(following_line)

        if (
            current is None
            or following is None
            or current["indent"] != following["indent"]
            or _inventory_row(lines, number)
            or LIST_ITEM.match(following["body"])
            or FORMULA_START.match(following["body"])
        ):
            continue

        if current["pad"]:
            in_list = _introduced_by_list(lines, number, current["indent"])

            if not in_list:
                continue
        else:
            in_list = bool(LIST_ITEM.match(current["body"]))

            if not in_list and _follows_list_item(
                lines,
                number,
                current["indent"],
            ):
                continue

        if bool(following["pad"]) != in_list or (
            current["pad"] and current["pad"] != following["pad"]
        ):
            continue

        unit = units(following["body"])[0]
        text = line.rstrip()
        separator = 0 if re.search(r"\S-$", text) else 1

        if len(text) + separator + len(unit) <= WRAP_COLUMN:
            breaks.append((number, unit))

    return breaks


def _paragraph(
    lines: list[str],
    number: int,
    candidates: set[int],
) -> tuple[int, int]:
    """
    Return the first and last line of the paragraph holding a break.

    The paragraph is the list item and its indented continuation lines when the
    break lies in an item, and otherwise the adjacent plain-prose lines at one
    indentation. It never extends past a candidate boundary, a directive, a
    fence, a list item, or a one-identifier inventory row.

    Parameters
    ----------
    lines : list[str]
        Physical source lines.
    number : int
        One-based line number of a premature break.
    candidates : set[int]
        One-based line numbers whose comments may be ordinary prose.

    Returns
    -------
    first, last : tuple[int, int]
        One-based first and last line numbers.
    """

    current = COMMENT_LINE.match(lines[number - 1])

    assert current is not None

    indent = current["indent"]

    def comment(index: int) -> re.Match[str] | None:
        if not 1 <= index <= len(lines) or index not in candidates:
            return None

        line = lines[index - 1]

        if "shellcheck" in line or COMMENT_FENCE.match(line):
            return None

        match = COMMENT_LINE.match(line)

        return match if match and match["indent"] == indent else None

    def prose(index: int) -> bool:
        match = comment(index)

        return bool(
            match
            and not match["pad"]
            and not LIST_ITEM.match(match["body"])
            and not _inventory_row(lines, index),
        )

    def continuation(index: int) -> bool:
        match = comment(index)

        return bool(match and match["pad"])

    if current["pad"] or LIST_ITEM.match(current["body"]):
        first = number

        # Walk up the item's indented continuation lines to the item line.
        while (above := COMMENT_LINE.match(lines[first - 1])) and above["pad"]:
            first -= 1

        last = first

        while continuation(last + 1):
            last += 1

        return first, last

    first = number

    while prose(first - 1):
        first -= 1

    last = number

    while prose(last + 1):
        last += 1

    return first, last


def _refill(lines: list[str], first: int, last: int) -> list[str]:
    """
    Refill one paragraph greedily, keeping its marker and indentation.

    Parameters
    ----------
    lines : list[str]
        Physical source lines.
    first : int
        One-based first line of the paragraph.
    last : int
        One-based last line of the paragraph.

    Returns
    -------
    refilled : list[str]
        The paragraph's new lines.
    """

    comments = [COMMENT_LINE.match(line) for line in lines[first - 1 : last]]
    head = comments[0]

    assert head is not None

    indent = head["indent"]
    text = " ".join(match["body"] for match in comments if match)
    item = LIST_ITEM.match(head["body"])
    lead = f"{indent}# "
    rest = lead + (" " * len(item.group(0)) if item else "")
    out: list[str] = []
    line = ""

    for unit in units(text):
        if not line:
            line = lead + unit
        elif len(line) + 1 + len(unit) <= WRAP_COLUMN:
            line += " " + unit
        else:
            out.append(line)
            line = rest + unit

    out.append(line)

    return out


def format_text(text: str, language: str) -> tuple[str, list[int]]:
    """
    Refill every paragraph that breaks early, and list what needs review.

    A paragraph with a line ending in an attached hyphen is never rejoined: the
    hyphen may break one word or be a suspended hyphen shared by two, and width
    cannot tell them apart.

    Parameters
    ----------
    text : str
        Shell or Python source.
    language : str
        'shell' or 'python'.

    Returns
    -------
    formatted, review : tuple[str, list[int]]
        Source with every eligible paragraph refilled, and the one-based line
        numbers of breaks left for review.
    """

    lines = text.split("\n")
    review: list[int] = []

    for _ in range(len(lines)):
        if language == "python":
            candidates = python_candidates("\n".join(lines))
        else:
            candidates = shell_candidates(lines)

        pending = [
            number
            for number, _ in premature_breaks(lines, candidates)
            if number not in review
        ]

        if not pending:
            break

        first, last = _paragraph(lines, pending[0], candidates)
        block = lines[first - 1 : last]

        if any(re.search(r"\S-$", line.rstrip()) for line in block[:-1]):
            review.append(pending[0])

            continue

        refilled = _refill(lines, first, last)

        # A refill that changes nothing cannot remove the break, so leave the
        # paragraph for review rather than retrying it.
        if refilled == block:
            review.append(pending[0])

            continue

        offset = len(refilled) - len(block)
        lines[first - 1 : last] = refilled
        review = [n + offset if n > last else n for n in review]

    formatted = "\n".join(lines)
    review.sort()

    return formatted, review


def main(argv: list[str] | None = None) -> int:
    """
    Preview by default and write only with an explicit flag.

    Parameters
    ----------
    argv : list[str] | None
        Explicit arguments, or None to read the process arguments.

    Returns
    -------
    status : int
        Zero when writing or when nothing would change, and one otherwise.
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+", type=Path)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args(argv)
    root = args.root.resolve()
    changed = False

    for path in args.paths:
        absolute = path if path.is_absolute() else root / path
        language = {".sh": "shell", ".py": "python"}.get(absolute.suffix)

        if language is None:
            print(f"skipped (not .sh or .py): {path}", file=sys.stderr)

            continue

        original = absolute.read_text(encoding="utf-8")
        formatted, review = format_text(original, language)

        for number in review:
            print(
                f"review: {path}:{number}: a line ends in a hyphen, so the "
                "paragraph is left as written",
                file=sys.stderr,
            )

        if formatted == original:
            continue

        changed = True

        if args.write:
            absolute.write_text(formatted, encoding="utf-8")
        else:
            print(path)

    return 0 if args.write or not changed else 1


if __name__ == "__main__":
    raise SystemExit(main())
