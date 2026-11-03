#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: check_logged_pseudo.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Check the '--pseudo' pair that a ratio step logged against an expected pair.

The log holds a line '--pseudo <A>:<B>', printed as floats. Exits 0 when both
values equal the expected 'A:B' to a relative tolerance of 1e-12, so a value
applied to the wrong pair cannot pass, and 1 otherwise.
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path


def parse_args(
    argv: list[str] | None = None,
) -> argparse.Namespace:
    """
    Parse the log path and the expected pair.
    """

    parser = argparse.ArgumentParser(
        description=(
            "Check the '--pseudo' pair that a ratio step logged against an "
            "expected pair."
        ),
    )
    parser.add_argument("log", help="Log of the ratio step.")
    parser.add_argument("spec", help="Expected pair, as 'A:B'.")

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """
    Exit 0 if the logged pair matches the expected one.
    """

    args = parse_args(argv)
    text = Path(args.log).read_text(encoding="utf-8")
    found = re.search(r"^--pseudo\s+(\S+):(\S+)$", text, re.MULTILINE)
    want = [float(v) for v in args.spec.split(":")]
    got = [float(v) for v in found.groups()] if found else []
    same = (
        len(got) == 2
        and len(want) == 2
        and all(
            math.isclose(g, w, rel_tol=1e-12)
            for g, w in zip(got, want, strict=True)
        )
    )

    return 0 if same else 1


if __name__ == "__main__":
    raise SystemExit(main())
