#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: json_value.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Print one value from a JSON file.

The expression is Python evaluated with the parsed file bound to 'value', as in
'value["summary"]["count"]'. It comes from the test that calls this script,
never from input data.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def parse_args(
    argv: list[str] | None = None,
) -> argparse.Namespace:
    """
    Parse the JSON path and the expression.
    """

    parser = argparse.ArgumentParser(
        description="Print one value from a JSON file.",
    )
    parser.add_argument("source", help="JSON file to read.")
    parser.add_argument(
        "expression", help="Python expression over 'value', the parsed file."
    )

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """
    Print the expression's value.
    """

    args = parse_args(argv)
    value = json.loads(Path(args.source).read_text(encoding="utf-8"))

    print(eval(args.expression, {"value": value}))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
