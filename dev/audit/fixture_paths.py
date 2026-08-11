#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: fixture_paths.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Exclude generated test fixtures from source-derived discovery.
"""

from __future__ import annotations

from pathlib import Path

# Generated fixtures live below this prefix and are ignored by Git, so a
# discovery that asks Git never sees them and needs nothing from this module.
# A discovery that globs the filesystem cannot see that they are ignored, and a
# fixture admitted to a maintained universe is treated as source: a
# deliberately malformed fixture becomes a finding, and a fixture help surface
# becomes a surface the repository owes documentation. Both were observed once
# fixtures took ordinary '.py' and '.sh' suffixes.
FIXTURE_ROOT = "tests/fixtures"


def is_fixture_path(path: Path, root: Path) -> bool:
    """
    Return whether a path is a generated fixture below the repository root.

    Parameters
    ----------
    path : Path
        Candidate path, absolute or relative to 'root'.

    root : Path
        Repository root that 'path' is measured against.

    Returns
    -------
    result : bool
        True when the path lies under the generated-fixture root.
    """

    try:
        relative = path.relative_to(root)
    except ValueError:
        relative = path

    return relative.as_posix().startswith(f"{FIXTURE_ROOT}/")


def is_fixture_recipe(path: str) -> bool:
    """
    Return whether a repository-relative path is a fixture recipe.

    Parameters
    ----------
    path : str
        Repository-relative path in POSIX form.

    Returns
    -------
    result : bool
        True when the path is a 'make.sh' below the fixture root.
    """

    return path.startswith(f"{FIXTURE_ROOT}/") and path.endswith("/make.sh")
