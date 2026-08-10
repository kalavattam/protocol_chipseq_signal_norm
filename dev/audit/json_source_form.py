#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: json_source_form.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Check maintained JSON source against one canonical rendering.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import subprocess
from collections.abc import Callable
from pathlib import Path

RULE_ID = "JSON.SOURCE.FORM"

# The budget decides structure, not line length. A structure fits when its
# complete physical line fits; an indivisible scalar longer than the budget is
# left alone, because JSON has no string continuation.
BUDGET = 79
STEP = 2
# Serialized inventories under 'dev/audit/' are owned by the producer's
# serializer rather than by an author, so an author-facing rule cannot reach
# them. They are excluded by not appearing here, which is an applicability
# boundary rather than an allowlist; a separate exclusion tuple would be
# unreachable code, because no path can start with both.
MAINTAINED_ROOTS = ("dev/config/", "dev/schemas/", "tests/fixtures/")


@dataclasses.dataclass
class Finding:
    """
    Represent one canonical-form finding.
    """

    rule_id: str
    path: str
    line: int
    message: str

    def format(self) -> str:
        """
        Return the one-line report form for this finding.
        """

        return f"{self.rule_id}: {self.path}:{self.line}: {self.message}"


def _reject_duplicate_keys(pairs: list[tuple[str, object]]) -> dict:
    """
    Build one object, refusing to collapse a duplicated key.

    A duplicate key is not a form defect, but `json.loads` keeps only the last
    occurrence. Rewriting such a file would silently delete data, so parsing
    fails instead and the caller reports the file as unreadable.

    Parameters
    ----------
    pairs : list[tuple[str, object]]
        Member pairs in document order.

    Returns
    -------
    value : dict
        The object built from unique keys.

    Raises
    ------
    ValueError
        If one key appears more than once in the same object.
    """

    value: dict = {}

    for key, member in pairs:
        if key in value:
            raise ValueError(f"duplicate key {key!r}")

        value[key] = member

    return value


def render(
    value: object,
    column: int = 0,
    line_indent: int = 0,
    tail: int = 0,
) -> list[str]:
    """
    Render one value in canonical form.

    A container is kept inline when its complete physical line fits the budget,
    and is otherwise expanded one member per line. Member order is preserved
    exactly as authored, because key ordering is not governed.

    Parameters
    ----------
    value : object
        One parsed JSON value.
    column : int
        Column at which the value's first character sits, used for the fit
        test. This includes any key prefix already placed on the line.
    line_indent : int
        Indentation of the line the value begins on. Members indent from this
        rather than from `column`, so a long key never cascades its width into
        the members below it.
    tail : int
        Width of any separator that will follow the value's last line.

    Returns
    -------
    lines : list[str]
        Rendered lines. The first carries no indentation, because the caller
        places it; every later line carries its absolute indentation.
    """

    if not isinstance(value, (dict, list)) or not value:
        return [json.dumps(value, ensure_ascii=False)]

    flat = json.dumps(value, ensure_ascii=False, separators=(", ", ": "))

    if column + len(flat) + tail <= BUDGET:
        return [flat]

    if isinstance(value, dict):
        opening, closing = "{", "}"
        members = [
            (json.dumps(key, ensure_ascii=False) + ": ", member)
            for key, member in value.items()
        ]
    else:
        opening, closing = "[", "]"
        members = [("", member) for member in value]

    inner = line_indent + STEP
    lines = [opening]

    for position, (prefix, member) in enumerate(members):
        final = position == len(members) - 1
        rendered = render(
            member, inner + len(prefix), inner, 0 if final else 1
        )
        rendered[0] = " " * inner + prefix + rendered[0]

        if not final:
            rendered[-1] += ","

        lines.extend(rendered)

    lines.append(" " * line_indent + closing)

    return lines


def canonical(value: object) -> str:
    """
    Return the complete canonical text for one parsed document.
    """

    return "\n".join(render(value)) + "\n"


def _containers(text: str) -> list[tuple[int, int]]:
    """
    Return every container's opening and closing offset, outermost first.

    Parameters
    ----------
    text : str
        Complete JSON source text.

    Returns
    -------
    spans : list[tuple[int, int]]
        Offset pairs sorted by opening offset.
    """

    stack: list[int] = []
    spans: list[tuple[int, int]] = []
    in_string = False
    escaped = False

    for offset, character in enumerate(text):
        if in_string:
            if escaped:
                escaped = False
            elif character == "\\":
                escaped = True
            elif character == '"':
                in_string = False

            continue

        if character == '"':
            in_string = True
        elif character in "[{":
            stack.append(offset)
        elif character in "]}":
            spans.append((stack.pop(), offset))

    return sorted(spans)


def _line_starts(lines: list[str]) -> list[int]:
    """
    Return the offset at which each line begins.
    """

    starts: list[int] = []
    position = 0

    for line in lines:
        starts.append(position)
        position += len(line) + 1

    return starts


def _locator(starts: list[int]) -> Callable[[int], int]:
    """
    Return a function mapping one offset to its zero-based line number.
    """

    def locate(offset: int) -> int:
        low, high = 0, len(starts) - 1

        while low < high:
            middle = (low + high + 1) // 2

            if starts[middle] <= offset:
                low = middle
            else:
                high = middle - 1

        return low

    return locate


def _depth_at_line_start(text: str, count: int) -> list[int]:
    """
    Return the container depth in force at the start of each line.
    """

    depths = [0] * (count + 1)
    depth = 0
    line = 0
    in_string = False
    escaped = False

    for character in text:
        if character == "\n":
            line += 1

            if line <= count:
                depths[line] = depth

            continue

        if in_string:
            if escaped:
                escaped = False
            elif character == "\\":
                escaped = True
            elif character == '"':
                in_string = False

            continue

        if character == '"':
            in_string = True
        elif character in "[{":
            depth += 1
        elif character in "]}":
            depth -= 1

    return depths


def _structural_findings(
    text: str,
    path: str,
    lines: list[str],
) -> list[Finding]:
    """
    Report every recognized departure from the canonical structure.

    Parameters
    ----------
    text : str
        Complete JSON source text.
    path : str
        Repository-relative path used in the report.
    lines : list[str]
        Source lines with the final empty element removed.

    Returns
    -------
    findings : list[Finding]
        Structural findings in source order.
    """

    findings: list[Finding] = []
    starts = _line_starts(lines)
    locate = _locator(starts)
    reported_overflow: set[int] = set()

    for opening, closing in _containers(text):
        first, last = locate(opening), locate(closing)

        if first == last:
            if len(lines[first]) > BUDGET and first not in reported_overflow:
                reported_overflow.add(first)
                findings.append(
                    Finding(
                        RULE_ID,
                        path,
                        first + 1,
                        f"inline structure runs to {len(lines[first])} "
                        f"columns, past the {BUDGET}-column budget; expand it",
                    ),
                )

            continue

        column = opening - starts[first]
        following = text[closing + 1 : closing + 2]
        member = json.loads(
            text[opening : closing + 1],
            object_pairs_hook=_reject_duplicate_keys,
        )
        flat = json.dumps(member, ensure_ascii=False, separators=(", ", ": "))

        if column + len(flat) + (following == ",") <= BUDGET:
            findings.append(
                Finding(
                    RULE_ID,
                    path,
                    first + 1,
                    f"expanded structure fits the {BUDGET}-column budget "
                    f"inline at {column + len(flat)} columns; inline it",
                ),
            )

        trailing = text[opening + 1 : starts[first] + len(lines[first])]

        if trailing.strip():
            findings.append(
                Finding(
                    RULE_ID,
                    path,
                    first + 1,
                    "opening delimiter of an expanded structure must be the "
                    "last content on its line",
                ),
            )

        leading = text[starts[last] : closing]

        if leading.strip():
            findings.append(
                Finding(
                    RULE_ID,
                    path,
                    last + 1,
                    "closing delimiter of an expanded structure must be the "
                    "first content on its line",
                ),
            )

    depths = _depth_at_line_start(text, len(lines))

    for number, line in enumerate(lines, 1):
        if not line.strip():
            findings.append(
                Finding(RULE_ID, path, number, "blank line"),
            )

            continue

        if "\t" in line:
            findings.append(
                Finding(RULE_ID, path, number, "line contains a tab"),
            )

        if line != line.rstrip():
            findings.append(
                Finding(RULE_ID, path, number, "line has trailing whitespace"),
            )

        depth = depths[number - 1]

        if line.lstrip()[0] in "]}":
            depth -= 1

        found = len(line) - len(line.lstrip(" "))

        if found != depth * STEP:
            findings.append(
                Finding(
                    RULE_ID,
                    path,
                    number,
                    f"indentation is {found}, expected {depth * STEP}",
                ),
            )

    return findings


def check_text(text: str, path: str) -> list[Finding]:
    """
    Report every canonical-form finding for one JSON document.

    Parameters
    ----------
    text : str
        Complete JSON source text.
    path : str
        Repository-relative path used in the report.

    Returns
    -------
    findings : list[Finding]
        Findings in source order, empty when the text is already canonical.
    """

    try:
        value = json.loads(text, object_pairs_hook=_reject_duplicate_keys)
    except ValueError as error:
        return [Finding(RULE_ID, path, 1, f"unreadable JSON: {error}")]

    expected = canonical(value)

    if text == expected:
        return []

    lines = text.split("\n")

    if lines and lines[-1] == "":
        lines.pop()
    else:
        return [
            Finding(
                RULE_ID,
                path,
                max(len(lines), 1),
                "file must end with exactly one newline",
            ),
        ]

    if text.endswith("\n\n"):
        return [
            Finding(
                RULE_ID,
                path,
                len(lines),
                "file must end with exactly one newline",
            ),
        ]

    findings = _structural_findings(text, path, lines)

    if findings:
        return findings

    # Every recognized class is clean, yet the text still differs. Reporting
    # the first differing line keeps the checker from passing a file it cannot
    # explain, which is the failure mode a silent canonical comparison invites.
    produced = expected.split("\n")

    paired = zip(lines, produced, strict=False)

    for number, (found, wanted) in enumerate(paired, 1):
        if found != wanted:
            return [
                Finding(
                    RULE_ID,
                    path,
                    number,
                    "line differs from the canonical rendering; expected "
                    f"{wanted.strip()!r}",
                ),
            ]

    return [
        Finding(
            RULE_ID,
            path,
            len(lines),
            f"file has {len(lines)} lines; the canonical rendering has "
            f"{len(produced) - 1}",
        ),
    ]


def maintained_paths(root: Path) -> list[str]:
    """
    Return current maintained JSON paths without serialized inventories.

    Discovery lists tracked and non-ignored files, so generated fixtures stay
    invisible to the checker even while they sit under a maintained root.

    Parameters
    ----------
    root : Path
        Repository root to discover within.

    Returns
    -------
    paths : list[str]
        Sorted repository-relative paths.
    """

    result = subprocess.run(
        [
            "git",
            "ls-files",
            "-z",
            "--cached",
            "--others",
            "--exclude-standard",
            "--",
            "*.json",
        ],
        cwd=root,
        check=True,
        capture_output=True,
    )

    return sorted(
        path
        for item in result.stdout.split(b"\0")
        if item
        and (path := item.decode("utf-8")).startswith(MAINTAINED_ROOTS)
        and (root / path).is_file()
    )


def scan(root: Path, paths: list[str] | None = None) -> list[Finding]:
    """
    Check explicit paths or current maintained JSON source.

    Parameters
    ----------
    root : Path
        Repository root to resolve paths against.
    paths : list[str] | None
        Explicit repository-relative paths, or None to discover them.

    Returns
    -------
    findings : list[Finding]
        Findings across every inspected path.
    """

    root = root.resolve()
    selected = paths if paths is not None else maintained_paths(root)

    return [
        finding
        for path in sorted(set(selected))
        if (root / path).is_file()
        for finding in check_text(
            (root / path).read_text(encoding="utf-8"),
            path,
        )
    ]


def rewrite(root: Path, paths: list[str]) -> list[str]:
    """
    Rewrite each readable path in canonical form and return those changed.

    Parameters
    ----------
    root : Path
        Repository root to resolve paths against.
    paths : list[str]
        Repository-relative paths to rewrite.

    Returns
    -------
    changed : list[str]
        Paths whose bytes changed.
    """

    changed: list[str] = []

    for path in sorted(set(paths)):
        target = root / path
        text = target.read_text(encoding="utf-8")
        value = json.loads(text, object_pairs_hook=_reject_duplicate_keys)
        expected = canonical(value)

        if expected != text:
            target.write_text(expected, encoding="utf-8")
            changed.append(path)

    return changed


def main(argv: list[str] | None = None) -> int:
    """
    Report or repair maintained JSON canonical-form findings.

    Parameters
    ----------
    argv : list[str] | None
        Explicit arguments, or None to read the process arguments.

    Returns
    -------
    status : int
        Zero when the selected JSON source is canonical and one otherwise.
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--fix", action="store_true")
    parser.add_argument("paths", nargs="*")
    args = parser.parse_args(argv)

    root = args.root.resolve()
    selected = args.paths or maintained_paths(root)

    if args.fix:
        for path in rewrite(root, selected):
            print(f"rewrote {path}")

    findings = scan(root, selected)

    if args.json:
        print(
            json.dumps(
                {
                    "paths": sorted(set(selected)),
                    "finding_count": len(findings),
                    "findings": [
                        dataclasses.asdict(item) for item in findings
                    ],
                },
                indent=2,
            ),
        )
    else:
        for finding in findings:
            print(finding.format())

        print(
            f"{RULE_ID}: {len(set(selected))} path(s), "
            f"{len(findings)} finding(s)",
        )

    return 1 if findings else 0


if __name__ == "__main__":
    raise SystemExit(main())
