#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: check_scaling_factor_order.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Check the order of options in the calculate-scaling-factor sources.

- 'region <file> <begin> <end> <token>...': the tokens appear in the given
  order between the first '<begin>' and the next '<end>' after it.
- 'fixtures <file>...': every command line with '--csv_mip' sets '--aln_typ'
  first, on that line or in the 12 lines before it, or builds on an 'arr_cmd'
  array that does. Prints the number of commands checked.

Each exits 0 when its order holds. Otherwise it exits nonzero with a message on
stderr. Arguments are read by position, not parsed as options.
"""

from __future__ import annotations

import sys
from pathlib import Path

USAGE = (
    "usage: check_scaling_factor_order.py region <file> <begin> <end> "
    "[<token> ...]\n"
    "       check_scaling_factor_order.py fixtures <file> [<file> ...]"
)


def check_region(file: str, begin: str, end: str, tokens: list[str]) -> None:
    """
    Raise SystemExit unless the tokens appear in order within the region.
    """

    text = Path(file).read_text(encoding="utf-8")
    start = text.find(begin)
    if start < 0:
        raise SystemExit(f"missing region start {begin!r}")

    stop = text.find(end, start + len(begin))
    if stop < 0:
        raise SystemExit(f"missing region end {end!r}")

    region = text[start:stop]
    cursor = -1
    for token in tokens:
        position = region.find(token, cursor + 1)
        if position < 0:
            raise SystemExit(f"missing ordered token {token!r}")
        cursor = position


def check_fixtures(files: list[str]) -> None:
    """
    Print the commands checked; raise SystemExit if any sets '--csv_mip' first.
    """

    failures = []
    checked = 0

    for name in files:
        path = Path(name)
        lines = path.read_text(encoding="utf-8").splitlines()
        for index, line in enumerate(lines):
            if "--csv_mip" not in line:
                continue

            prior = "\n".join(lines[max(0, index - 12) : index])
            uses_base_array = "arr_cmd" in prior
            if (
                "--aln_typ" not in prior
                and "--aln_typ" not in line
                and not uses_base_array
            ):
                failures.append(
                    f"{path}:{index + 1}: --csv_mip lacks prior --aln_typ"
                )
            checked += 1

    if failures:
        raise SystemExit("\n".join(failures))
    if checked == 0:
        raise SystemExit(
            "no calculate-scaling-factor fixture commands checked"
        )

    print(checked)


def main(argv: list[str] | None = None) -> int:
    """
    Run the subcommand's check. Arguments are read by position, since tokens
    such as '-at|' and '--aln_typ' look like options.
    """

    args = sys.argv[1:] if argv is None else argv

    if len(args) >= 4 and args[0] == "region":
        check_region(args[1], args[2], args[3], args[4:])
    elif len(args) >= 2 and args[0] == "fixtures":
        check_fixtures(args[1:])
    else:
        raise SystemExit(USAGE)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
