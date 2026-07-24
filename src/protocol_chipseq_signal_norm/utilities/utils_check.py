#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: utils_check.py
#
# Copyright 2025-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5-series models; most recent: GPT-5.6) were
# used in development and documentation.
#
# Distributed under the MIT license.


"""
Validation helpers shared by Python command-line scripts.

Notes
-----
Helpers raise or report argument and file validation errors consistently across
user-facing Python CLIs.
"""

from __future__ import annotations

import operator
import os
import sys
from collections.abc import Iterable
from typing import Any, Literal

assert sys.version_info >= (3, 11), "Python >= 3.11 required."

OPS: dict[str, tuple] = {
    "gt": (operator.gt, ">"),
    "ge": (operator.ge, ">="),
    "lt": (operator.lt, "<"),
    "le": (operator.le, "<="),
    "eq": (operator.eq, "=="),
    "ne": (operator.ne, "!="),
}

ALLOWED: tuple[str, ...] = ("bedGraph", "bedgraph", "bdg", "bg", "bed")


def fmt_label(label: str) -> str:
    """
    Format CLI-style label for error messages.
    """
    return label if label.startswith("--") else f"--{label}"


def cap_first(s: str) -> str:
    """
    Ensure the first letter in a string is capitalized
    """
    return s[:1].upper() + s[1:]


def as_iter(x: Any) -> tuple[Any, ...]:
    """
    Return a tuple for any non-string, non-bytes, non-PathLike iterable;
    otherwise wrap as a 1-tuple. Guarantees stable container semantics for
    broadcasting and avoids per-character iteration on strings.
    """
    #  Treat obvious scalars as scalars
    if isinstance(x, (str, bytes, os.PathLike)):
        return (x,)

    #  Is it iterable?
    try:
        iter(x)
    except TypeError:
        return (x,)

    #  If it is iterable, coerce to tuple to stabilize broadcasting semantics
    try:
        return tuple(x)
    except TypeError:
        #  Fallback for rare non-reiterable iterables
        return (x,)


def pair_val_thresh(vals: Any, threshs: Any) -> list[tuple[Any, Any]]:
    """
    Broadcast a scalar threshold across values, or pair elementwise when both
    are iterables of the same length. Raises on length mismatch.
    """
    vals = as_iter(vals)
    threshs = as_iter(threshs)

    if len(threshs) == 1:
        threshs = threshs * len(vals)
    elif len(threshs) != len(vals):
        raise ValueError(f"Length mismatch: {len(vals)=} vs {len(threshs)=}.")

    return list(zip(vals, threshs, strict=True))


def check_cmp(
    value: int | float | Iterable[int | float | None],
    comp: str,
    thresh: int | float | Iterable[int | float],
    label: str,
    allow_none: bool = True
) -> None:
    """
    Check that value(s) satisfy a comparison against a threshold. Works with
    scalars or iterables; thresholds can be scalar (broadcast) or iterable
    (paired).

    Parameters
    ----------
        value : int | float | iterable[int|float|None]
            The number (or numbers) to validate. Elements may be None when
            'allow_none' is True; such entries are skipped.
        comp : str
            One of 'gt', 'ge', 'lt', 'le', 'eq', 'ne'.
        thresh : int | float | iterable[int|float]
            The comparison target (scalar or per-element thresholds).
        label : str
            Short name used in error text (e.g., 'eps' becomes "Error: '--eps'
            ...").
        allow_none : bool
            If True, silently skip None values; if False, None triggers an
            error.

    Returns
    -------
        None.

    Raises
    ------
        ValueError
            If any value violates the comparison or if the comparison itself is
            invalid.
    """
    try:
        op, sym = OPS[comp]
    except KeyError as e:
        raise ValueError(
            f"Unknown comparison '{comp}'. Expected one of "
            f"{', '.join(OPS.keys())}."
        ) from e

    lbl = fmt_label(label)

    for val, thr in pair_val_thresh(value, thresh):
        if val is None:
            if allow_none:
                continue
            raise ValueError(f"'{lbl}' must be {sym} {thr}.")
        if not op(val, thr):
            raise ValueError(f"'{lbl}' must be {sym} {thr}.")


def check_parse_fil_out(
    value: os.PathLike[str] | str,
    allowed: Iterable[str] = ALLOWED
) -> tuple[str, str, bool]:
    """
    Check and parse an output path, inferring format and gzip compression from
    the extension.

    Parameters
    ----------
        value : os.PathLike[str] | str
            Output path provided by the user. If it ends with ".gz", output is
            gzip-compressed.
        allowed : Iterable[str] = ALLOWED
            Iterable of allowed base extensions. Exact lowercase spellings are
            accepted, and exact camelCase spellings explicitly present in
            'allowed' (e.g., 'bedGraph') are also accepted. Leading dots are
            ignored.

    Returns
    -------
        fil_out : str
            Validated output filename, preserving the accepted extension
            spelling and optional '.gz'.
        ext : str
            Validated output format token. Returns 'bedGraph' when that exact
            spelling is supplied and accepted; otherwise returns the accepted
            lowercase form ('bedgraph', 'bdg', 'bg', or 'bed').
        is_gz : bool
            True if the filename ends with '.gz', else False.

    Raises
    ------
        ValueError
            If the base extension is not one of the allowed values.

    Notes
    -----
        Refactored out of 'compute_signal.py' for modularization.
    """
    #  Coerce PathLike to str
    value = os.fspath(value)

    #  Preserve exact allowed spellings (e.g., 'bedGraph') while also
    #  recognizing lowercase canonical forms
    allowed_exact = {ext.lstrip(".") for ext in allowed}
    allowed_lower = {ext.lower() for ext in allowed_exact}

    #  Check if the extension is '.gz'
    is_gz = value.endswith(".gz")

    #  If extension is '.gz', remove '.gz' and extract the base extension
    base = value[:-3] if is_gz else value
    stem, ext_raw = os.path.splitext(base)

    ext_pres = ext_raw.lstrip(".")
    ext_low = ext_pres.lower()

    #  Accept either:
    #    1. an exact allowed spelling (e.g., 'bedGraph'), or
    #    2. an all-lowercase spelling whose lowercase form is allowed
    if not (
        ext_pres in allowed_exact or
        (ext_pres == ext_low and ext_low in allowed_lower)
    ):
        allowed_ord = tuple(ext.lstrip(".") for ext in allowed)
        allowed_lst = ", ".join(f".{e}" for e in allowed_ord)
        raise ValueError(
            f"Invalid extension '{ext_pres}'; allowed: {allowed_lst}."
        )

    #  Reconstruct the final output filename while preserving the accepted
    #  extension spelling
    fil_out = (
        f"{stem}.{ext_pres}.gz"
        if is_gz
        else f"{stem}.{ext_pres}"
    )

    ext = "bedGraph" if ext_pres == "bedGraph" else ext_low

    return fil_out, ext, is_gz


def check_exists(
    path: os.PathLike[str] | str,
    kind: Literal["file", "dir", "any"] = "any",
    label: str | None = None
) -> None:
    """
    Ensure 'path' exists, optionally as a specific kind.

    Parameters
    ----------
        path : str
            File or directory path.
        kind : {'file', 'dir', 'any'}, default 'any'
            - 'file': require an existing regular file
            - 'dir' : require an existing directory
            - 'any' : existence check only (file, dir, or other)
        label : str | None = None
            Optional label for clearer error text (e.g., "First file (A)").

    Returns
    -------
        None

    Raises
    ------
        FileNotFoundError
            If the required path does not exist or is not of the requested
            kind.
    """
    p = os.fspath(path)

    if kind == "file":
        ok, want = os.path.isfile(p), "file"
    elif kind == "dir":
        ok, want = os.path.isdir(p), "directory"
    elif kind == "any":
        ok, want = os.path.exists(p), "path"
    else:
        raise ValueError(
            f"Unknown kind: {kind!r} (expected 'file', 'dir', or 'any')."
        )

    if ok:
        return

    #  Decide how to render the noun in the error message
    what = label or want

    #  Special-case bedGraph labels: keep 'bedGraph' lowercase even at start
    #  of the message, e.g. "bedGraph A not found: ..."
    if what.lower().startswith("bedgraph"):
        msg_what = what
    else:
        msg_what = cap_first(what)

    raise FileNotFoundError(f"{msg_what} not found: {p}")


def check_writable(
    path: os.PathLike[str] | str,
    kind: Literal["file", "dir"] = "file",
    must_exist: bool = False,
    label: str | None = None
) -> None:
    """
    Ensure the writability of a file or directory.

    For files, the parent directory exists and is writable/enterable. If
    'must_exist=True' and the file exists, the file itself must be writable.

    For directories, the directory itself must exist and be writable/enterable.

    Parameters
    ----------
        path : str
            Target file path ('kind="file"') or directory path ('kind="dir"').
        kind : {'file', 'dir'}, default 'file'
            What to validate.
        must_exist : bool
            When 'kind="file"' and the file already exists, require the file
            itself to be writable (default: False).
        label : str | None
            Optional label for nicer error text.

    Returns
    -------
        None.

    Raises
    ------
        FileNotFoundError
            If the relevant directory (or the directory itself when
            'kind="dir"') does not exist.
        PermissionError
            If the directory (or file, when 'must_exist=True') is not writable.
        IsADirectoryError
            If the path points to a directory when a file is expected.
    """
    p = os.fspath(path)

    if kind == "dir":
        dir_path = p
        if not os.path.isdir(dir_path):
            what = label or "directory"
            raise FileNotFoundError(
                f"{cap_first(what)} does not exist: {dir_path}"
            )
        #  To create files within directories, check for write and execute
        #  permissions
        if not os.access(dir_path, os.W_OK | os.X_OK):
            what = label or "directory"
            raise PermissionError(
                f"No write permission for {what.lower()}: {dir_path}"
            )
        return

    #  Otherwise, handle when 'kind == "file"'
    if os.path.isdir(p):
        raise IsADirectoryError(
            f"Path points to a directory, not a file: {p}"
        )

    dir_path = os.path.dirname(p) or "."
    if not os.path.isdir(dir_path):
        raise FileNotFoundError(f"Output directory does not exist: {dir_path}")
    if not os.access(dir_path, os.W_OK | os.X_OK):
        raise PermissionError(
            f"No write permission for output directory: {dir_path}"
        )

    if must_exist and os.path.exists(p) and not os.access(p, os.W_OK):
        what = label or "file"
        raise PermissionError(f"{cap_first(what)} is not writable: {p}")
