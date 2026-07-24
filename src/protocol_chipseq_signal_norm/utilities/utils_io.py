#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: utils_io.py
#
# Copyright 2025-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5-series models; most recent: GPT-5.6) were
# used in development and documentation.
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
            File path or '-' for stdin. '.gz' is auto-detected.

    Returns
    -------
        ContextManager[TextIO]
            A context manager yielding a readable text stream. For '-',
            returns a no-op context over 'sys.stdin' (not closed on exit).

    Raises
    ------
        None.
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
            File path or '-' for stdout. '.gz' is auto-detected.

    Returns
    -------
        ContextManager[TextIO]
            A context manager yielding a writable text stream. For '-',
            returns a no-op context over 'sys.stdout' (not closed on exit).

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
    csv: str | None, default: tuple[str, ...] = DEF_SKP_PFX
) -> tuple[str, ...]:
    """
    Parse a comma-separated list of header prefixes.

    Parameters
    ----------
        csv: str | None
            Comma-delimited prefixes (e.g., "#,track,browser"). Special cases:
                - None: use 'default'
                - "" (empty string): disable skipping
                - "__default__": use 'default'
        default : tuple[str, ...] = DEF_SKP_PFX
            Fallback tuple of prefixes when 'csv' is None or "__default__".

    Returns
    -------
        tuple() | default | tuple(tok.strip() for tok in s.split(",")
        if tok.strip())
            The effective prefix tuple to use for header detection.

    Raises
    ------
        None.
    """
    if csv is None:
        return default

    s = csv.strip()
    if not s:
        return ()

    if s == "__default__":
        return default

    return tuple(tok.strip() for tok in s.split(",") if tok.strip())


def is_header(line: str, skp_pfx: tuple[str, ...] = DEF_SKP_PFX) -> bool:
    """
    Predicate for header/metadata lines.

    Parameters
    ----------
        line : str
            Raw input line (with or without trailing newline).
        skp_pfx: tuple[str, ...] = DEF_SKP_PFX
            Prefixes that mark a line as a header after left-stripping.

    Returns
    -------
        bool
            True for blank/whitespace lines or lines whose left-stripped
            content starts with any prefix in 'skp_pfx'; otherwise False.

    Raises
    ------
        None.
    """
    s = line.lstrip()
    return not s or s.startswith(skp_pfx)


def read_data_line(
    fh: TextIO,
    skp_pfx: tuple[str, ...] = DEF_SKP_PFX
) -> str:
    """
    Return the next non-empty, non-header bedGraph line from an open text
    stream. Lines beginning with any prefix in 'skp_pfx' (e.g., '#', 'track',
    'browser') are skipped. Returns '' on EOF.

    Parameters
    ----------
        fh : TextIO
            Open text-mode handle positioned at the current read location.
        skp_pfx : tuple[str, ...]
            Prefixes to skip as bedGraph header/meta lines.

    Returns
    -------
        s : str
            The next data line stripped of trailing newline or, at EOF, ''.

    Raises
    ------
        None.

    Notes
    -----
        Refactored out of 'compute_signal_ratio.py' for modularization.
    """
    while True:
        lin = fh.readline()
        if not lin:
            return ""

        if is_header(lin, skp_pfx):
            continue

        s = lin.strip()
        if not s:
            continue

        return s


def ensure_single_stdin(paths: list[str]) -> None:
    """
    Enforce at most one '-' (stdin) among provided paths.

    Parameters
    ----------
        paths : list[str]
            The positional path arguments provided by the user.

    Returns
    -------
        None.

    Raises
    ------
        ValueError
            If more than one '-' is present.
    """
    if paths.count("-") > 1:
        raise ValueError("At most one '-' (stdin) path is allowed.")
