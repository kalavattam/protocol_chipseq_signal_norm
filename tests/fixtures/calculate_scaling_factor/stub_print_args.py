#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: stub_print_args.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Stand in for the Python BAM-to-BED converter.

Write the arguments given to the file named after '--fil_out', so a test can
see what a caller passed the converter.
"""

from __future__ import annotations

import sys


def main(argv: list[str] | None = None) -> int:
    """
    Write the arguments, space-separated, to the file after '--fil_out'.
    """

    args = sys.argv[1:] if argv is None else argv
    fil_out = args[args.index("--fil_out") + 1]

    with open(fil_out, "w", encoding="utf-8") as fil:
        fil.write(" ".join(args) + "\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
