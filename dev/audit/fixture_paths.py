#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: fixture_paths.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5, Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Exclude generated test fixtures from source-derived discovery.
"""

from __future__ import annotations

import re
from pathlib import Path

# Generated fixtures live below this prefix. Git ignores them, but a discovery
# that globs the filesystem would treat them as source, making a malformed
# fixture a finding and a fixture's help a surface owed documentation.
FIXTURE_ROOT = "tests/fixtures"

# The tracked sources beside a root's recipe, which 'testing.md' names: the
# recipe, its generators, and its stubs. They are maintained source, not
# generated fixtures.
FIXTURE_SOURCE = re.compile(
    rf"^{FIXTURE_ROOT}/[^/]+/"
    r"(?:make\.sh|make_\w+\.py|stub_\w+\.(?:sh|py))$",
)


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
        True when the path lies under the generated-fixture root and is not one
        of the tracked sources beside a recipe.
    """

    try:
        relative = path.relative_to(root)
    except ValueError:
        relative = path

    posix = relative.as_posix()
    under_root = posix.startswith(f"{FIXTURE_ROOT}/")

    return under_root and not FIXTURE_SOURCE.match(posix)


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
        True when the path is the 'make.sh' of a fixture root. A 'make.sh'
        deeper down is generated fixture data, not a recipe.
    """

    return re.fullmatch(rf"{FIXTURE_ROOT}/[^/]+/make\.sh", path) is not None
