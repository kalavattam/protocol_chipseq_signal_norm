#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: make_legacy_shellcheck_runner.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Write a copy of 'run_shellcheck.sh' without its empty-array guards.

Bash before 4.4 treats '"${arr[@]}"' of an empty array as unset under 'set -u',
so the runner guards each expansion. The copy drops the four guards, so the
test can show the fault they prevent: under an older '/bin/bash', the copy
stops with 'unbound variable'. Exits with a message when it does not find
exactly four guards, so a change to the runner surfaces here rather than as a
silently identical copy.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

# An empty-array guard: 'if (( ${#arr_<name>[@]} > 0 )); then', the guarded
# lines one level in, an optional 'else' branch for the empty case, and 'fi'.
GUARD = re.compile(
    r"^(?P<indent> *)if \(\( \$\{#arr_\w+\[@\]\} > 0 \)\); then\n"
    r"(?P<body>(?:(?P=indent) {4}.*\n)+?)"
    r"(?:(?P=indent)else\n(?:(?P=indent) {4}.*\n)+?)?"
    r"(?P=indent)fi\n",
    re.MULTILINE,
)

# The runner's guards: two around the path listings, two around the runs.
GUARD_COUNT = 4


def unguard(match: re.Match[str]) -> str:
    """
    Keep the guarded lines, one level out, and drop the rest.
    """

    return re.sub(r"^ {4}", "", match.group("body"), flags=re.MULTILINE)


def parse_args(
    argv: list[str] | None = None,
) -> argparse.Namespace:
    """
    Parse the runner and the destination.
    """

    parser = argparse.ArgumentParser(
        description=(
            "Write a copy of 'run_shellcheck.sh' without its empty-array "
            "guards."
        ),
    )
    parser.add_argument("runner", help="Path of 'run_shellcheck.sh'.")
    parser.add_argument("destination", help="Path of the copy to write.")

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """
    Write the unguarded copy.
    """

    args = parse_args(argv)
    source = Path(args.runner).read_text(encoding="utf-8")

    source, count = GUARD.subn(unguard, source)
    if count != GUARD_COUNT:
        raise SystemExit(
            f"expected {GUARD_COUNT} empty-array guards, found {count}"
        )

    Path(args.destination).write_text(source, encoding="utf-8")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
