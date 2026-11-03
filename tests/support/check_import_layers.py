#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: check_import_layers.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Check the import directions between the package's layers.

Under '<root>/src/protocol_chipseq_signal_norm', nothing imports the retired
'scripts' package, a utility never imports the CLI layer, and a CLI never
imports another CLI. Prints one line per violation and exits 1 if there is any.
"""

from __future__ import annotations

import argparse
import ast
from pathlib import Path

PKG = "protocol_chipseq_signal_norm"


def parse_args(
    argv: list[str] | None = None,
) -> argparse.Namespace:
    """
    Parse the repository root.
    """

    parser = argparse.ArgumentParser(
        description=(
            "Check the import directions between the package's layers."
        ),
    )
    parser.add_argument("root", help="Repository root.")

    return parser.parse_args(argv)


def find_violations(root: Path) -> list[str]:
    """
    Return one line per import that crosses a layer the wrong way.
    """

    findings = []

    for path in sorted(root.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        role = path.relative_to(root).parts[0]

        for node in ast.walk(tree):
            names = []

            if isinstance(node, ast.Import):
                names = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom):
                names = [node.module or ""]

            for name in names:
                if name == "scripts" or name.startswith("scripts."):
                    findings.append(
                        f"{path}:{node.lineno}: retired scripts import",
                    )

                if role == "utilities" and ".cli" in name:
                    findings.append(
                        f"{path}:{node.lineno}: utility imports CLI layer",
                    )

                if role == "cli" and name.startswith(f"{PKG}.cli."):
                    findings.append(
                        f"{path}:{node.lineno}: CLI imports peer CLI",
                    )

    return findings


def main(argv: list[str] | None = None) -> int:
    """
    Print each violation and exit 1 if there is any.
    """

    args = parse_args(argv)
    findings = find_violations(Path(args.root) / "src" / PKG)

    for line in findings:
        print(line)

    return 1 if findings else 0


if __name__ == "__main__":
    raise SystemExit(main())
