#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: print_conda_pkg.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Print a Conda package's version and build as '<version>-<build>'.

The package is found by its record in '<prefix>/conda-meta/', and the record's
'name' field must match, so a package whose name starts with another's (such as
'bwa-mem2' for 'bwa') is never taken for it. The values come from the record's
fields rather than from its file name. Prints nothing when no record matches.
"""

from __future__ import annotations

import argparse
import glob
import json
import sys


def parse_args(
    argv: list[str] | None = None,
) -> argparse.Namespace:
    """
    Parse the package name and the environment prefix.
    """

    parser = argparse.ArgumentParser(
        description=(
            "Print a Conda package's version and build as "
            "'<version>-<build>', from its package record."
        ),
    )
    parser.add_argument("pkg", help="Package name, such as 'bowtie2'.")
    parser.add_argument("prefix", help="Environment prefix.")

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """
    Print '<version>-<build>' for each record whose name matches.
    """

    args = parse_args(argv)

    try:
        for path in sorted(
            glob.glob(f"{args.prefix}/conda-meta/{args.pkg}-[0-9]*.json")
        ):
            with open(path, encoding="utf-8") as handle:
                rec = json.load(handle)

            if rec["name"] == args.pkg:
                print(f"{rec['version']}-{rec['build']}")
    except (OSError, KeyError, ValueError) as err:
        print(f"error(print_conda_pkg.py): {err}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
