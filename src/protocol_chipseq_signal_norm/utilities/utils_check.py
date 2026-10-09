#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: utils_check.py
#
# Copyright 2025-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5-series models; most recent: GPT-5.6);
# - Anthropic Claude Code (Opus 5, Opus 5.5).
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
import warnings
from collections.abc import Iterable
from typing import Any, Literal

assert sys.version_info >= (3, 11), "Python >= 3.11 required."

COMPARISON_OPERATORS: dict[str, tuple] = {
    "gt": (operator.gt, ">"),
    "ge": (operator.ge, ">="),
    "lt": (operator.lt, "<"),
    "le": (operator.le, "<="),
    "eq": (operator.eq, "=="),
    "ne": (operator.ne, "!="),
}

# Formats the tools read and write, by canonical name. A format name or file
# suffix is matched in any letter case and resolves to its canonical spelling.
FORMATS: tuple[str, ...] = ("bam", "cram", "bed", "bedGraph")

# Formats whose files may end in '.gz'.
FORMATS_GZIP: frozenset[str] = frozenset({"bed", "bedGraph"})

# Short bedGraph spellings, refused in favor of 'bedGraph'.
BEDGRAPH_SHORT: frozenset[str] = frozenset({"bdg", "bg"})

ALLOWED_OUTPUT_FORMATS: tuple[str, ...] = ("bedGraph", "bed")

_FORMAT_CANON = {name.casefold(): name for name in FORMATS}


def format_label(label: str) -> str:
    """
    Format CLI-style label for error messages.

    Parameters
    ----------
    label : str
        Human-readable diagnostic label.

    Returns
    -------
    formatted : str
        A display label suitable for diagnostics.
    """

    return label if label.startswith("--") else f"--{label}"


def capitalize_first(value: str) -> str:
    """
    Ensure the first letter in a string is capitalized.

    Parameters
    ----------
    value : str
        String value to capitalize.

    Returns
    -------
    capitalized : str
        The value with its first character capitalized.
    """

    return value[:1].upper() + value[1:]


def as_tuple(value: Any) -> tuple[Any, ...]:
    """
    Normalize a value to stable tuple semantics.

    Strings, bytes, and path-like values remain scalar rather than being
    iterated character by character.

    Parameters
    ----------
    value : Any
        Scalar or iterable value to normalize.

    Returns
    -------
    values : tuple[Any, ...]
        A tuple view of the supplied scalar or iterable value.
    """

    if isinstance(value, (str, bytes, os.PathLike)):
        return (value,)

    try:
        iter(value)
    except TypeError:
        return (value,)

    try:
        return tuple(value)
    except TypeError:
        return (value,)


def pair_values_and_thresholds(
    values: Any,
    thresholds: Any,
) -> list[tuple[Any, Any]]:
    """
    Pair values with scalar or elementwise thresholds.

    A single threshold is broadcast; other threshold collections must match the
    number of values.

    Parameters
    ----------
    values : Any
        Values to pair with one or more thresholds.
    thresholds : Any
        Scalar or elementwise thresholds paired with the values.

    Returns
    -------
    pairs : list[tuple[Any, Any]]
        Values paired with their applicable thresholds.

    Raises
    ------
    ValueError
        If a supplied value violates the validated contract.
    """

    normalized_values = as_tuple(values)
    normalized_thresholds = as_tuple(thresholds)

    if len(normalized_thresholds) == 1:
        normalized_thresholds *= len(normalized_values)
    elif len(normalized_thresholds) != len(normalized_values):
        raise ValueError(
            "Length mismatch: "
            f"values={len(normalized_values)} vs "
            f"thresholds={len(normalized_thresholds)}.",
        )

    return list(
        zip(
            normalized_values,
            normalized_thresholds,
            strict=True,
        ),
    )


def validate_comparison(
    value: int | float | Iterable[int | float | None],
    comparison: str,
    threshold: int | float | Iterable[int | float],
    label: str,
    allow_none: bool = True,
) -> None:
    """
    Check that value(s) satisfy a comparison against a threshold.

    Scalar thresholds are broadcast. Iterable thresholds are paired.

    Parameters
    ----------
    value : int | float | Iterable[int | float | None]
        The number (or numbers) to validate. Elements may be None when
        'allow_none' is True; such entries are skipped.
    comparison : str
        One of 'gt', 'ge', 'lt', 'le', 'eq', 'ne'.
    threshold : int | float | Iterable[int | float]
        The comparison target (scalar or per-element thresholds).
    label : str
        Short name used in error text (e.g., 'eps' becomes "Error: '--eps'
        ...").
    allow_none : bool
        If True, silently skip None values; if False, None triggers an error.

    Raises
    ------
    ValueError
        If any value violates the comparison or if the comparison itself is
        invalid.
    """

    try:
        operation, symbol = COMPARISON_OPERATORS[comparison]
    except KeyError as error:
        raise ValueError(
            f"Unknown comparison '{comparison}'. Expected one of "
            f"{', '.join(COMPARISON_OPERATORS)}.",
        ) from error

    formatted_label = format_label(label)

    for candidate_value, threshold_value in pair_values_and_thresholds(
        value,
        threshold,
    ):
        if candidate_value is None:
            if allow_none:
                continue

            raise ValueError(
                f"'{formatted_label}' must be {symbol} {threshold_value}.",
            )

        if not operation(candidate_value, threshold_value):
            raise ValueError(
                f"'{formatted_label}' must be {symbol} {threshold_value}.",
            )


def canonicalize_format(name: str, allowed: Iterable[str] = FORMATS) -> str:
    """
    Resolve a format name, in any letter case, to its canonical spelling.

    Parameters
    ----------
    name : str
        Format name as supplied, such as 'bedgraph' or 'CRAM'.
    allowed : Iterable[str]
        Canonical formats accepted here (default: 'FORMATS').

    Returns
    -------
    fmt : str
        Canonical format name, such as 'bedGraph' or 'cram'.

    Raises
    ------
    ValueError
        If 'name' is a short bedGraph spelling ('bdg', 'bg') or not one of
        'allowed'.
    """

    allowed = tuple(allowed)
    key = name.strip().casefold()

    if key in BEDGRAPH_SHORT:
        raise ValueError(f"Format '{name}' is not accepted; use 'bedGraph'.")

    fmt = _FORMAT_CANON.get(key)

    if fmt is None or fmt not in allowed:
        quoted = [f"'{entry}'" for entry in allowed]
        choices = (
            quoted[0]
            if len(quoted) == 1
            else f"{', '.join(quoted[:-1])}, or {quoted[-1]}"
        )

        raise ValueError(f"Invalid format '{name}'; choose {choices}.")

    return fmt


def _split_suffix(path: os.PathLike[str] | str) -> tuple[str, str]:
    """
    Split a path's format suffix from an optional trailing '.gz'.

    Parameters
    ----------
    path : os.PathLike[str] | str
        File path.

    Returns
    -------
    ext, gz : tuple[str, str]
        Suffix before any '.gz', without its dot and as typed; and the '.gz' as
        typed, in any letter case as 'open_in()' and 'open_out()' accept it,
        otherwise ''.
    """

    value = os.fspath(path)
    gz = value[-3:] if value.lower().endswith(".gz") else ""
    base = value[: -len(gz)] if gz else value

    return os.path.splitext(base)[1].lstrip("."), gz


def _refuse_short(path: os.PathLike[str] | str, ext: str, gz: str) -> None:
    """
    Refuse a short bedGraph suffix, naming the accepted one.

    Parameters
    ----------
    path : os.PathLike[str] | str
        File path.
    ext, gz : str
        The path's suffix and optional '.gz', from '_split_suffix()'.

    Raises
    ------
    ValueError
        If 'ext' is a short bedGraph spelling ('bdg', 'bg').
    """

    if ext.casefold() in BEDGRAPH_SHORT:
        want = ".bedGraph.gz" if gz else ".bedGraph"

        raise ValueError(
            f"'{os.fspath(path)}': '.{ext}{gz}' is not accepted; the file "
            f"name must end in '{want}'.",
        )


def format_from_path(path: os.PathLike[str] | str) -> str | None:
    """
    Recognize a file's format from its suffix, in any letter case.

    Parameters
    ----------
    path : os.PathLike[str] | str
        File path, such as 'x.cram' or 'x.BEDGRAPH.gz'.

    Returns
    -------
    fmt : str | None
        Canonical format, or None when the suffix names no format; '.gz' is
        recognized only after a BED or bedGraph suffix.

    Raises
    ------
    ValueError
        If the suffix is a short bedGraph spelling ('.bdg', '.bg').

    Notes
    -----
    The file is never renamed; only its format is recognized.
    """

    ext, gz = _split_suffix(path)
    _refuse_short(path, ext, gz)
    fmt = _FORMAT_CANON.get(ext.casefold())

    if fmt is None or (gz and fmt not in FORMATS_GZIP):
        return None

    return fmt


def check_bedgraph_path(path: os.PathLike[str] | str) -> None:
    """
    Refuse a bedGraph input whose name uses a short suffix.

    Parameters
    ----------
    path : os.PathLike[str] | str
        bedGraph input path, or '-' for standard input.

    Raises
    ------
    ValueError
        If the path ends in '.bdg' or '.bg', with or without '.gz'.

    Notes
    -----
    Any other name is read as bedGraph, so piped inputs and temporary files
    keep working; only the short spellings are refused.
    """

    ext, gz = _split_suffix(path)
    _refuse_short(path, ext, gz)


def validate_output_path(
    value: os.PathLike[str] | str,
    allowed: Iterable[str] = ALLOWED_OUTPUT_FORMATS,
) -> tuple[str, str, bool]:
    """
    Validate an output path and infer its format and compression.

    Parameters
    ----------
    value : os.PathLike[str] | str
        Output path provided by the user. If it ends with '.gz', output is
        gzip-compressed.
    allowed : Iterable[str]
        Canonical formats accepted here (default: 'ALLOWED_OUTPUT_FORMATS').

    Returns
    -------
    output_path, fmt, is_compressed : tuple[str, str, bool]
        The output path as given; its canonical format, such as 'bedGraph'; and
        whether it ends with '.gz'.

    Raises
    ------
    ValueError
        If the suffix is a short bedGraph spelling, or names no format in
        'allowed'.

    Warns
    -----
    UserWarning
        If the suffix differs from the canonical one only in letter case, such
        as '.BEDGRAPH.gz'; the path is still written as given.
    """

    value = os.fspath(value)
    allowed = tuple(allowed)
    ext, gz = _split_suffix(value)
    _refuse_short(value, ext, gz)
    fmt = _FORMAT_CANON.get(ext.casefold())

    if fmt is None or fmt not in allowed or (gz and fmt not in FORMATS_GZIP):
        allowed_text = ", ".join(f"'.{entry}'" for entry in allowed)

        raise ValueError(
            f"Invalid extension '.{ext}'; allowed: {allowed_text}, each "
            "optionally followed by '.gz'.",
        )

    canonical = f".{fmt}.gz" if gz else f".{fmt}"

    if f".{ext}{gz}" != canonical:
        warnings.warn(
            f"writing '{value}' as given; the canonical suffix is "
            f"'{canonical}'.",
            stacklevel=2,
        )

    return value, fmt, bool(gz)


def check_exists(
    path: os.PathLike[str] | str,
    kind: Literal["file", "dir", "any"] = "any",
    label: str | None = None,
) -> None:
    """
    Ensure 'path' exists, optionally as a specific kind.

    Parameters
    ----------
    path : os.PathLike[str] | str
        File or directory path.
    kind : Literal["file", "dir", "any"], default "any"
        - 'file': require an existing regular file
        - 'dir' : require an existing directory
        - 'any' : existence check only (file, dir, or other)
    label : str | None = None
        Optional label for clearer error text (e.g., "First file (A)").

    Raises
    ------
    FileNotFoundError
        If the required path does not exist or is not of the requested kind.
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
            f"Unknown kind: {kind!r} (expected 'file', 'dir', or 'any').",
        )

    if ok:
        return

    what = label or want

    # Keep `bedGraph` lowercase when it begins a diagnostic.
    if what.lower().startswith("bedgraph"):
        target_description = what
    else:
        target_description = capitalize_first(what)

    raise FileNotFoundError(f"{target_description} not found: {p}")


def check_writable(
    path: os.PathLike[str] | str,
    kind: Literal["file", "dir"] = "file",
    must_exist: bool = False,
    label: str | None = None,
) -> None:
    """
    Ensure the writability of a file or directory.

    For files, the parent directory exists and is writable/enterable. If
    'must_exist=True' and the file exists, the file itself must be writable.

    For directories, the directory itself must exist and be writable/enterable.

    Parameters
    ----------
    path : os.PathLike[str] | str
        Target file path ('kind="file"') or directory path ('kind="dir"').
    kind : Literal["file", "dir"], default "file"
        What to validate.
    must_exist : bool
        When 'kind="file"' and the file already exists, require the file itself
        to be writable (default: False). Ignored with a warning when
        'kind="dir"'.
    label : str | None
        Optional label for nicer error text.

    Raises
    ------
    FileNotFoundError
        If the relevant directory (or the directory itself when 'kind="dir"')
        does not exist.
    PermissionError
        If the directory (or file, when 'must_exist=True') is not writable.
    IsADirectoryError
        If the path points to a directory when a file is expected.
    """

    p = os.fspath(path)

    if kind == "dir":
        if must_exist:
            warnings.warn(
                "'must_exist' has no effect with kind='dir' and is ignored.",
                stacklevel=2,
            )

        dir_path = p

        if not os.path.isdir(dir_path):
            what = label or "directory"
            raise FileNotFoundError(
                f"{capitalize_first(what)} does not exist: {dir_path}",
            )

        if not os.access(dir_path, os.W_OK | os.X_OK):
            what = label or "directory"
            raise PermissionError(
                f"No write permission for {what.lower()}: {dir_path}",
            )

        return

    if os.path.isdir(p):
        raise IsADirectoryError(f"Path points to a directory, not a file: {p}")

    dir_path = os.path.dirname(p) or "."
    if not os.path.isdir(dir_path):
        raise FileNotFoundError(f"Output directory does not exist: {dir_path}")

    if not os.access(dir_path, os.W_OK | os.X_OK):
        raise PermissionError(
            f"No write permission for output directory: {dir_path}",
        )

    if must_exist and os.path.exists(p) and not os.access(p, os.W_OK):
        what = label or "file"
        raise PermissionError(f"{capitalize_first(what)} is not writable: {p}")
