#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: shell_source_form.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


"""
Check a deliberately narrow subset of maintained Bash source form.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import re
import subprocess
from pathlib import Path

RULE_ID = "SHELL.SOURCE.FORM"
DIAGNOSTIC_RULE_ID = "SHELL.DIAGNOSTIC.FORM"
FUNCTION = re.compile(
    r"^\s*function\s+(?P<name>[A-Za-z_][A-Za-z0-9_]*)\(\)\s*\{",
)
ARRAY = re.compile(
    r"^\s*(?:local|declare)\s+-(?:[A-Za-z]*[aA][A-Za-z]*)\s+"
    r"(?P<name>[A-Za-z_][A-Za-z0-9_]*)",
)
HEREDOC = re.compile(
    r"<<-?\s*(?:'(?P<single>[A-Za-z_][A-Za-z0-9_]*)'|"
    r'"(?P<double>[A-Za-z_][A-Za-z0-9_]*)"|'
    r"(?P<plain>[A-Za-z_][A-Za-z0-9_]*))",
)
MESSAGE_TYPE = r"[a-z][a-z0-9_]*(?:_[a-z0-9]+)*"
DIAGNOSTIC = re.compile(
    r'^(?P<indent>[ \t]*)echo[ \t]+"'
    rf"(?P<prefix>{MESSAGE_TYPE}\(.+\):)"
    r'(?P<message>.*)"(?P<tail>[ \t].*)?$',
)
CONTINUED_DIAGNOSTIC = re.compile(
    rf'^[ \t]+"(?P<prefix>{MESSAGE_TYPE}\(.+\):)'
    r'(?P<message>.*)"(?P<tail>[ \t].*)?$',
)
EMPTY_DIAGNOSTIC_CONTEXT = re.compile(
    rf'^[ \t]*echo[ \t]+"{MESSAGE_TYPE}\(\):',
)
SINGLE_QUOTED_DIAGNOSTIC = re.compile(
    rf"^[ \t]*echo[ \t]+'{MESSAGE_TYPE}\(.+\):",
)
DIAGNOSTIC_CONTINUATION = re.compile(
    r'^(?P<indent>[ \t]+)"(?P<message>.*)"(?P<tail>[ \t].*)?$',
)
ECHO_CONTINUED = re.compile(r"^[ \t]*echo[ \t]+\\$")
SNAKE_CASE = re.compile(r"^[a-z][a-z0-9]*(?:_[a-z0-9]+)*$")
MAINTAINED_ROOTS = ("bin/", "lib/bash/", "install/scripts/", "tests/")
POSIX_BOOTSTRAP = "install/scripts/install_envs_entrypoint.sh"
DIAGNOSTIC_INDENT_MESSAGE = (
    "diagnostic continuation must be indented one additional level"
)
DIAGNOSTIC_PREFIX_MESSAGE = (
    "diagnostic prefix must be the first quoted argument by itself"
)


@dataclasses.dataclass(frozen=True)
class Finding:
    """
    Represent one bounded shell-source finding.
    """

    rule_id: str
    path: str
    line: int
    message: str


def _append_diagnostic_finding(
    findings: list[Finding],
    path: str,
    line: int,
    message: str,
) -> None:
    """
    Append one diagnostic-form finding.
    """

    findings.append(Finding(DIAGNOSTIC_RULE_ID, path, line, message))


def _check_split_diagnostic_echo(
    lines: list[str],
    index: int,
    path: str,
    findings: list[Finding],
) -> bool:
    """
    Check a diagnostic whose echo command precedes its first argument.
    """

    if not ECHO_CONTINUED.fullmatch(lines[index]) or index + 1 >= len(lines):
        return False

    argument = CONTINUED_DIAGNOSTIC.fullmatch(lines[index + 1])
    if argument is None:
        return True

    _append_diagnostic_finding(
        findings,
        path,
        index + 1,
        "diagnostic echo and prefix must share the first source line",
    )

    if argument.group("message"):
        _append_diagnostic_finding(
            findings,
            path,
            index + 2,
            "diagnostic prefix must be the first quoted argument by itself",
        )

    if len(lines[index + 1]) > 79:
        _append_diagnostic_finding(
            findings,
            path,
            index + 2,
            "recognized diagnostic source line exceeds 79 columns",
        )

    return True


def _check_invalid_diagnostic_prefix(
    line: str,
    line_number: int,
    path: str,
    findings: list[Finding],
) -> bool:
    """
    Check empty-context and single-quoted diagnostic prefixes.
    """

    if EMPTY_DIAGNOSTIC_CONTEXT.match(line):
        _append_diagnostic_finding(
            findings,
            path,
            line_number,
            "diagnostic script-or-function context must be nonempty",
        )

        return True

    if not SINGLE_QUOTED_DIAGNOSTIC.match(line):
        return False

    _append_diagnostic_finding(
        findings,
        path,
        line_number,
        "diagnostic prefix must use a double-quoted first argument",
    )

    return True


def _check_diagnostic_continuations(
    lines: list[str],
    index: int,
    diagnostic: re.Match[str],
    path: str,
    findings: list[Finding],
) -> int:
    """
    Check adjacent prose arguments and return the next source index.
    """

    line_number = index + 1
    expected_indent = len(diagnostic.group("indent")) + 4
    continuation = index + 1

    while continuation < len(lines):
        current = lines[continuation]
        current_number = continuation + 1

        if len(current) > 79:
            _append_diagnostic_finding(
                findings,
                path,
                current_number,
                "recognized diagnostic source line exceeds 79 columns",
            )

        argument = DIAGNOSTIC_CONTINUATION.fullmatch(current)

        if argument is None:
            _append_diagnostic_finding(
                findings,
                path,
                current_number,
                "diagnostic continuation must be an adjacent double-quoted "
                "argument",
            )

            return continuation

        if len(argument.group("indent")) != expected_indent:
            _append_diagnostic_finding(
                findings,
                path,
                current_number,
                DIAGNOSTIC_INDENT_MESSAGE,
            )

        message = argument.group("message")

        if not message or message != message.strip():
            _append_diagnostic_finding(
                findings,
                path,
                current_number,
                "diagnostic prose argument must be nonempty without edge "
                "whitespace",
            )

        tail = (argument.group("tail") or "").strip()

        if tail == "\\":
            continuation += 1

            continue

        if tail not in {"", ">&2"}:
            _append_diagnostic_finding(
                findings,
                path,
                current_number,
                "diagnostic final line has an unrecognized source suffix",
            )

        return continuation + 1

    _append_diagnostic_finding(
        findings,
        path,
        line_number,
        "diagnostic continuation is incomplete",
    )

    return len(lines)


def check_diagnostic_forms(lines: list[str], path: str) -> list[Finding]:
    """
    Check fixture-defined direct-echo diagnostic source forms.

    Parameters
    ----------
    lines : list[str]
        Physical source lines used for location-aware checks.
    path : str
        Repository-relative path associated with the source.

    Returns
    -------
    findings : list[Finding]
        Shell diagnostic-form findings.
    """

    findings: list[Finding] = []
    heredoc_delimiter: str | None = None
    index = 0

    while index < len(lines):
        line = lines[index]

        if heredoc_delimiter is not None:
            if line.strip() == heredoc_delimiter:
                heredoc_delimiter = None

            index += 1

            continue

        heredoc = HEREDOC.search(line)

        if heredoc is not None:
            heredoc_delimiter = next(
                value for value in heredoc.groups() if value is not None
            )
            index += 1

            continue

        if _check_split_diagnostic_echo(lines, index, path, findings):
            index += 1

            continue

        if _check_invalid_diagnostic_prefix(
            line,
            index + 1,
            path,
            findings,
        ):
            index += 1

            continue

        diagnostic = DIAGNOSTIC.fullmatch(line)

        if diagnostic is None:
            index += 1

            continue

        line_number = index + 1

        if len(line) > 79:
            _append_diagnostic_finding(
                findings,
                path,
                line_number,
                "recognized diagnostic source line exceeds 79 columns",
            )

        if diagnostic.group("message"):
            _append_diagnostic_finding(
                findings,
                path,
                line_number,
                DIAGNOSTIC_PREFIX_MESSAGE,
            )
            index += 1

            continue

        if (diagnostic.group("tail") or "").strip() != "\\":
            _append_diagnostic_finding(
                findings,
                path,
                line_number,
                "diagnostic prefix line must continue to a quoted prose "
                "argument",
            )
            index += 1

            continue

        index = _check_diagnostic_continuations(
            lines,
            index,
            diagnostic,
            path,
            findings,
        )

    return findings


def check_text(text: str, path: str = "<memory>") -> list[Finding]:
    """
    Check recognized Bash declarations, indentation, and heredoc forms.

    Parameters
    ----------
    text : str
        Shell source to inspect without executing it.
    path : str
        Diagnostic path associated with the source text.

    Returns
    -------
    findings : list[Finding]
        Stable diagnostics for the deliberately bounded source grammar.
    """

    lines = text.splitlines()
    findings = check_diagnostic_forms(lines, path)

    if path == POSIX_BOOTSTRAP:
        return findings

    heredoc_delimiter: str | None = None

    for index, line in enumerate(lines, 1):
        if heredoc_delimiter is not None:
            if line.strip() == heredoc_delimiter:
                heredoc_delimiter = None

            continue

        leading = line[: len(line) - len(line.lstrip(" \t"))]

        if "\t" in leading:
            findings.append(
                Finding(
                    RULE_ID,
                    path,
                    index,
                    "leading indentation contains a tab",
                ),
            )
        elif (
            leading
            and len(leading) % 4
            and not line.lstrip().startswith(("#", ")"))
        ):
            findings.append(
                Finding(
                    RULE_ID,
                    path,
                    index,
                    "recognized code indentation is not a four-space multiple",
                ),
            )

        function = FUNCTION.match(line)

        if function is not None and not SNAKE_CASE.fullmatch(
            function.group("name"),
        ):
            findings.append(
                Finding(
                    RULE_ID,
                    path,
                    index,
                    "recognized function name is not snake_case",
                ),
            )

        array = ARRAY.match(line)

        if array is not None and not array.group("name").startswith("arr_"):
            findings.append(
                Finding(
                    RULE_ID,
                    path,
                    index,
                    "recognized array declaration should use the arr_ prefix",
                ),
            )

        heredoc = HEREDOC.search(line)

        if heredoc is not None:
            delimiter = next(
                value for value in heredoc.groups() if value is not None
            )
            heredoc_delimiter = delimiter

            if delimiter == "EOF":
                findings.append(
                    Finding(
                        RULE_ID,
                        path,
                        index,
                        (
                            "ordinary shell-text heredoc should use EOM "
                            "instead of EOF"
                        ),
                    ),
                )

    return findings


def maintained_paths(root: Path) -> list[str]:
    """
    Return current maintained shell paths without generated evidence.
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
            "*.sh",
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
    Check explicit paths or current maintained Bash source.
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


def main(argv: list[str] | None = None) -> int:
    """
    Report bounded shell source-form findings.

    Parameters
    ----------
    argv : list[str] | None
        Explicit arguments, or None to read the process arguments.

    Returns
    -------
    status : int
        Zero when selected shell source forms pass and one otherwise.
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--json", action="store_true")
    parser.add_argument("paths", nargs="*")
    args = parser.parse_args(argv)

    findings = scan(args.root, args.paths or None)

    if args.json:
        print(
            json.dumps(
                [dataclasses.asdict(item) for item in findings],
                indent=2,
            ),
        )
    else:
        for finding in findings:
            print(
                f"{finding.rule_id}: {finding.path}:{finding.line}: "
                f"{finding.message}",
            )

    return 1 if findings else 0


if __name__ == "__main__":
    raise SystemExit(main())
