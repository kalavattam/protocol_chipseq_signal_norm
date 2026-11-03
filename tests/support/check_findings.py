#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: check_findings.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Check that a checker's JSON report holds no deterministic finding.

The report is a JSON list of findings, each with a 'classification'. Exits 1 if
any finding is classified 'deterministic', and 0 otherwise.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def parse_args(
    argv: list[str] | None = None,
) -> argparse.Namespace:
    """
    Parse the report path.
    """

    parser = argparse.ArgumentParser(
        description=(
            "Exit 1 if a checker's JSON report holds a deterministic finding."
        ),
    )
    parser.add_argument("report", help="JSON report of findings.")

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """
    Exit 1 if the report holds a deterministic finding.
    """

    args = parse_args(argv)
    findings = json.loads(Path(args.report).read_text(encoding="utf-8"))

    return int(
        any(item["classification"] == "deterministic" for item in findings),
    )


if __name__ == "__main__":
    raise SystemExit(main())
