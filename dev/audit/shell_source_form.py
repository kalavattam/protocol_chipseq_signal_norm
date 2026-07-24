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


"""Check a deliberately narrow subset of maintained Bash source form."""

from __future__ import annotations

import argparse
import dataclasses
import json
import re
import subprocess
from pathlib import Path

RULE_ID = "SHELL.SOURCE.FORM"
FUNCTION = re.compile(r"^\s*function\s+(?P<name>[A-Za-z_][A-Za-z0-9_]*)\(\)\s*\{")
ARRAY = re.compile(
    r"^\s*(?:local|declare)\s+-(?:[A-Za-z]*[aA][A-Za-z]*)\s+"
    r"(?P<name>[A-Za-z_][A-Za-z0-9_]*)"
)
HEREDOC = re.compile(
    r"<<-?\s*(?:'(?P<single>[A-Za-z_][A-Za-z0-9_]*)'|"
    r'"(?P<double>[A-Za-z_][A-Za-z0-9_]*)"|'
    r"(?P<plain>[A-Za-z_][A-Za-z0-9_]*))"
)
SNAKE_CASE = re.compile(r"^[a-z][a-z0-9]*(?:_[a-z0-9]+)*$")
MAINTAINED_ROOTS = ("bin/", "lib/bash/", "install/scripts/", "tests/")
POSIX_BOOTSTRAP = "install/scripts/install_envs_entrypoint.sh"


@dataclasses.dataclass(frozen=True)
class Finding:
    """Represent one bounded shell-source finding."""

    rule_id: str
    path: str
    line: int
    message: str


def check_text(text: str, path: str = "<memory>") -> list[Finding]:
    """Check recognized Bash declarations, indentation, and heredoc forms."""

    if path == POSIX_BOOTSTRAP:
        return []
    lines = text.splitlines()
    findings: list[Finding] = []
    heredoc_delimiter: str | None = None
    for index, line in enumerate(lines, 1):
        if heredoc_delimiter is not None:
            if line.strip() == heredoc_delimiter:
                heredoc_delimiter = None
            continue
        leading = line[: len(line) - len(line.lstrip(" \t"))]
        if "\t" in leading:
            findings.append(
                Finding(RULE_ID, path, index, "leading indentation contains a tab")
            )
        elif leading and len(leading) % 4 and not line.lstrip().startswith(("#", ")")):
            findings.append(
                Finding(
                    RULE_ID,
                    path,
                    index,
                    "recognized code indentation is not a four-space multiple",
                )
            )
        function = FUNCTION.match(line)
        if function is not None and not SNAKE_CASE.fullmatch(function.group("name")):
            findings.append(
                Finding(
                    RULE_ID,
                    path,
                    index,
                    "recognized function name is not snake_case",
                )
            )
        array = ARRAY.match(line)
        if array is not None and not array.group("name").startswith("arr_"):
            findings.append(
                Finding(
                    RULE_ID,
                    path,
                    index,
                    "recognized array declaration should use the arr_ prefix",
                )
            )
        heredoc = HEREDOC.search(line)
        if heredoc is not None:
            delimiter = next(value for value in heredoc.groups() if value is not None)
            heredoc_delimiter = delimiter
            if delimiter == "EOF":
                findings.append(
                    Finding(
                        RULE_ID,
                        path,
                        index,
                        "ordinary shell-text heredoc should use EOM instead of EOF",
                    )
                )
    return findings


def maintained_paths(root: Path) -> list[str]:
    """Return current maintained shell paths without generated evidence."""

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
    """Check explicit paths or current maintained Bash source."""

    root = root.resolve()
    selected = paths if paths is not None else maintained_paths(root)
    return [
        finding
        for path in sorted(set(selected))
        if (root / path).is_file()
        for finding in check_text((root / path).read_text(encoding="utf-8"), path)
    ]


def main(argv: list[str] | None = None) -> int:
    """Report bounded shell source-form findings."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--json", action="store_true")
    parser.add_argument("paths", nargs="*")
    args = parser.parse_args(argv)
    findings = scan(args.root, args.paths or None)
    if args.json:
        print(json.dumps([dataclasses.asdict(item) for item in findings], indent=2))
    else:
        for finding in findings:
            print(
                f"{finding.rule_id}: {finding.path}:{finding.line}: "
                f"{finding.message}"
            )
    return 1 if findings else 0


if __name__ == "__main__":
    raise SystemExit(main())
