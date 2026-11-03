#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: check_python_package.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Check how the installed package is laid out, found, and registered.

Each subcommand exits 0 when its check holds and 1 when it does not:
- 'packages <root>': setuptools finds exactly the intended packages under
  '<root>/src'.
- 'import-path': the package imports from its 'src' tree, wherever the caller
  runs it.
- 'entry-points': every console script is registered, and each one loads the
  same 'main' as its module.
"""

from __future__ import annotations

import argparse
from importlib import import_module
from importlib.metadata import entry_points
from pathlib import Path

PKG = "protocol_chipseq_signal_norm"

PACKAGES = {PKG, f"{PKG}.cli", f"{PKG}.utilities"}

CONSOLE_SCRIPTS = {
    "add_coeffs_namespaced",
    "calculate_scaling_factor_siqchip",
    "calculate_scaling_factor_spike",
    "compute_input_floor",
    "compute_pseudo",
    "compute_pseudo_deeptools",
    "compute_signal",
    "compute_signal_ratio",
    "merge_bins_bdg",
    "parse_metadata_siqchip",
    "relativize_scaling_factors",
    "sum_bdg",
}


def parse_args(
    argv: list[str] | None = None,
) -> argparse.Namespace:
    """
    Parse the subcommand and, for 'packages', the repository root.
    """

    parser = argparse.ArgumentParser(
        description=(
            "Check how the installed package is laid out, found, and "
            "registered."
        ),
    )
    sub = parser.add_subparsers(dest="cmd", required=True)
    pkgs = sub.add_parser(
        "packages", help="setuptools finds only the intended packages."
    )
    pkgs.add_argument("root", help="Repository root.")
    sub.add_parser(
        "import-path", help="The package imports from its 'src' tree."
    )
    sub.add_parser(
        "entry-points",
        help="Each console script loads the same 'main' as its module.",
    )

    return parser.parse_args(argv)


def check_packages(root: str) -> bool:
    """
    Return whether setuptools finds exactly the intended packages.
    """

    from setuptools import find_packages

    return set(find_packages(where=str(Path(root) / "src"))) == PACKAGES


def check_import_path() -> bool:
    """
    Return whether the package imports from its 'src' tree.
    """

    path = Path(import_module(PKG).__file__).resolve()

    return path.name == "__init__.py" and "src" in path.parts


def check_entry_points() -> bool:
    """
    Return whether each console script loads its module's 'main'.
    """

    registered = {
        item.name: item
        for item in entry_points(group="console_scripts")
        if item.name in CONSOLE_SCRIPTS
    }
    if set(registered) != CONSOLE_SCRIPTS:
        return False

    return all(
        entry.load() is import_module(f"{PKG}.cli.{name}").main
        for name, entry in registered.items()
    )


def main(argv: list[str] | None = None) -> int:
    """
    Run the subcommand's check and exit 0 if it holds.
    """

    args = parse_args(argv)

    if args.cmd == "packages":
        ok = check_packages(args.root)
    elif args.cmd == "import-path":
        ok = check_import_path()
    else:
        ok = check_entry_points()

    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
