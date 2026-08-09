#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: run_numeric_output_applicability.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.6);
# - Anthropic Claude Code (Opus 5).
#
# Distributed under the MIT license.


"""
Run numeric applicability with a temporary, deterministic inventory.
"""

from __future__ import annotations

import argparse
import tempfile
from pathlib import Path

from dev.audit import numeric_emission_inventory, numeric_output_applicability


def main(argv: list[str] | None = None) -> int:
    """
    Build only an ephemeral inventory, then run the registered checker.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path)
    parser.add_argument("--schema", type=Path)
    parser.add_argument("--inventory", type=Path)
    args = parser.parse_args(argv)
    forwarded: list[str] = ["--root", "."]
    if args.config is not None:
        forwarded.extend(["--config", str(args.config)])
    if args.schema is not None:
        forwarded.extend(["--schema", str(args.schema)])
    if args.inventory is not None:
        forwarded.extend(["--inventory", str(args.inventory)])
        return numeric_output_applicability.main(forwarded)
    with tempfile.TemporaryDirectory(
        prefix="numeric_output_applicability."
    ) as tmp:
        inventory = Path(tmp) / "inventory.json"
        result = numeric_emission_inventory.main(
            ["--root", ".", "--json-out", str(inventory)]
        )
        if result:
            return result
        forwarded.extend(["--inventory", str(inventory)])
        return numeric_output_applicability.main(forwarded)


if __name__ == "__main__":
    raise SystemExit(main())
