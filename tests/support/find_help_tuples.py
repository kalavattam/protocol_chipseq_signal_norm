#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: find_help_tuples.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Find help text written as a tuple in calls under a source tree.

A tuple passed as 'help', 'description', or 'epilog' prints as a tuple, not as
its joined text. 'metavar' takes a tuple only when 'nargs' is set, so a tuple
'metavar' without 'nargs' is reported too. Prints 'inspected <n>' for the
number of Python files read, then one line per finding, and exits 1 if there is
any.
"""

from __future__ import annotations

import argparse
import ast
import pathlib

# 'metavar' takes a tuple when 'nargs' is set; the rest are always strings.
KEY_STR = ("help", "description", "epilog")


def parse_args(
    argv: list[str] | None = None,
) -> argparse.Namespace:
    """
    Parse the source root.
    """

    parser = argparse.ArgumentParser(
        description="Find help text written as a tuple in calls.",
    )
    parser.add_argument("root", help="Source tree to search.")

    return parser.parse_args(argv)


def find_tuples(path: pathlib.Path) -> list[str]:
    """
    Return one finding per tuple-valued help keyword in a file.
    """

    found = []
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))

    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue

        keys = {kwd.arg for kwd in node.keywords}

        for kwd in node.keywords:
            if not isinstance(kwd.value, ast.Tuple):
                continue

            if kwd.arg in KEY_STR or (
                kwd.arg == "metavar" and "nargs" not in keys
            ):
                found.append(
                    f"{path}:{kwd.value.lineno}: {kwd.arg}= is a tuple"
                )

    return found


def main(argv: list[str] | None = None) -> int:
    """
    Print the file count and each finding; exit 1 if there is any.
    """

    args = parse_args(argv)
    paths = sorted(pathlib.Path(args.root).rglob("*.py"))
    found = [line for path in paths for line in find_tuples(path)]

    print(f"inspected {len(paths)}")

    for line in found:
        print(line)

    return 1 if found else 0


if __name__ == "__main__":
    raise SystemExit(main())
