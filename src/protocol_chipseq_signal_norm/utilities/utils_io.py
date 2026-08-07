#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: utils_io.py
#
# Copyright 2025-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5-series models; most recent: GPT-5.6);
# - Anthropic Claude Code (Opus 5).
#
# Distributed under the MIT license.


"""
Input/output helpers for text and compressed workflow files.

Notes
-----
Helpers centralize stdin/stdout handling, gzip-aware open functions, header
prefix parsing, and single-stdin validation.
"""

from __future__ import annotations

import gzip
import sys
from contextlib import AbstractContextManager, nullcontext
from typing import TextIO

assert sys.version_info >= (3, 11), "Python >= 3.11 required."

DEF_SKP_PFX: tuple[str, ...] = ("#", "track", "browser")


def open_in(path: str) -> AbstractContextManager[TextIO]:
    """
    Open a text input stream.

    Parameters
    ----------
    path : str
        File path or '-' for standard input. A '.gz' suffix is detected
        automatically.

    Returns
    -------
    manager : AbstractContextManager[TextIO]
        Context manager yielding a readable text stream. For '-', the manager
        wraps 'sys.stdin' without closing it on exit.
    """

    if path == "-":
        return nullcontext(sys.stdin)

    if path.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")

    return open(path, encoding="utf-8")


def open_out(path: str) -> AbstractContextManager[TextIO]:
    """
    Open a text output stream.

    Parameters
    ----------
    path : str
        File path or '-' for standard output. A '.gz' suffix is detected
        automatically.

    Returns
    -------
    manager : AbstractContextManager[TextIO]
        Context manager yielding a writable text stream. For '-', the manager
        wraps 'sys.stdout' without closing it on exit.

    Raises
    ------
    OSError
        If the file cannot be opened for writing.
    """

    if path == "-":
        return nullcontext(sys.stdout)

    if path.lower().endswith(".gz"):
        return gzip.open(path, "wt", encoding="utf-8", newline="\n")

    return open(path, "w", encoding="utf-8", newline="\n")


def parse_skp_pfx(
    csv: str | None,
    default: tuple[str, ...] = DEF_SKP_PFX,
) -> tuple[str, ...]:
    """
    Parse a comma-separated list of header prefixes.

    Parameters
    ----------
    csv : str | None
        Comma-delimited prefixes such as '#,track,browser'. 'None' and
        '__default__' select 'default'; an empty string disables skipping.
    default : tuple[str, ...], default=DEF_SKP_PFX
        Fallback prefixes.

    Returns
    -------
    prefixes : tuple[str, ...]
        Effective prefixes for header detection.
    """

    if csv is None:
        return default

    s = csv.strip()
    if not s:
        return ()

    if s == "__default__":
        return default

    return tuple(token.strip() for token in s.split(",") if token.strip())


def is_header(line: str, skp_pfx: tuple[str, ...] = DEF_SKP_PFX) -> bool:
    """
    Predicate for header/metadata lines.

    Parameters
    ----------
    line : str
        Raw input line, with or without a trailing newline.
    skp_pfx : tuple[str, ...], default=DEF_SKP_PFX
        Prefixes that mark a line as a header after left-stripping.

    Returns
    -------
    is_header : bool
        'True' for blank lines or lines whose left-stripped content starts with
        any configured prefix; otherwise, 'False'.
    """

    stripped = line.lstrip()

    return not stripped or stripped.startswith(skp_pfx)


def read_data_line(
    handle: TextIO,
    skp_pfx: tuple[str, ...] = DEF_SKP_PFX,
) -> str:
    """
    Return the next nonempty, non-header bedGraph line.

    Lines beginning with a configured prefix are skipped. An empty string is
    returned at end of file.

    Parameters
    ----------
    handle : TextIO
        Open text-mode handle positioned at the current read location.
    skp_pfx : tuple[str, ...], default=DEF_SKP_PFX
        Prefixes to skip as bedGraph header or metadata lines.

    Returns
    -------
    line : str
        Next data line stripped of its trailing newline, or '' at end of file.

    Notes
    -----
    This helper was extracted from 'compute_signal_ratio.py'.
    """

    while True:
        line = handle.readline()
        if not line:
            return ""

        if is_header(line, skp_pfx):
            continue

        stripped = line.strip()
        if not stripped:
            continue

        return stripped


def ensure_single_stdin(paths: list[str]) -> None:
    """
    Enforce at most one '-' (stdin) among provided paths.

    Parameters
    ----------
    paths : list[str]
        Positional path arguments provided by the user.

    Raises
    ------
    ValueError
        If more than one '-' is present.
    """

    if paths.count("-") > 1:
        raise ValueError("At most one '-' (stdin) path is allowed.")
