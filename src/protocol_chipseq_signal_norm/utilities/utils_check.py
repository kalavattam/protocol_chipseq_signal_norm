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

COMPARISON_OPERATORS: dict[str, tuple] = {
    "gt": (operator.gt, ">"),
    "ge": (operator.ge, ">="),
    "lt": (operator.lt, "<"),
    "le": (operator.le, "<="),
    "eq": (operator.eq, "=="),
    "ne": (operator.ne, "!="),
}

ALLOWED_OUTPUT_FORMATS: tuple[str, ...] = (
    "bedGraph",
    "bedgraph",
    "bdg",
    "bg",
    "bed",
)


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

    A single threshold is broadcast; other threshold collections must match
    the number of values.

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
        If True, silently skip None values; if False, None triggers an
        error.

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


def validate_output_path(
    value: os.PathLike[str] | str,
    allowed: Iterable[str] = ALLOWED_OUTPUT_FORMATS,
) -> tuple[str, str, bool]:
    """
    Validate an output path and infer its format and compression.

    Parameters
    ----------
    value : os.PathLike[str] | str
        Output path provided by the user. If it ends with ".gz", output is
        gzip-compressed.
    allowed : Iterable[str]
        Iterable of allowed base extensions. Exact lowercase spellings are
        accepted, and exact camelCase spellings explicitly present in
        'allowed' (e.g., 'bedGraph') are also accepted. Leading dots are
        ignored. Defaults to 'ALLOWED_OUTPUT_FORMATS'.

    Returns
    -------
    output_path, extension, is_compressed : tuple[str, str, bool]
        Validated output filename, preserving the accepted extension
        spelling and optional '.gz'; validated output format token, which is
        'bedGraph' when that exact
        spelling is supplied and accepted; otherwise returns the accepted
        lowercase form ('bedgraph', 'bdg', 'bg', or 'bed'); and whether the
        filename ends with '.gz'.

    Raises
    ------
    ValueError
        If the base extension is not one of the allowed values.

    Notes
    -----
    Refactored out of 'compute_signal.py' for modularization.
    """

    value = os.fspath(value)

    # Preserve exact spellings such as `bedGraph` while recognizing aliases.
    allowed_exact = {extension.lstrip(".") for extension in allowed}
    allowed_lower = {extension.lower() for extension in allowed_exact}

    is_compressed = value.endswith(".gz")

    base = value[:-3] if is_compressed else value
    stem, raw_extension = os.path.splitext(base)

    preserved_extension = raw_extension.lstrip(".")
    lowercase_extension = preserved_extension.lower()

    if not (
        preserved_extension in allowed_exact
        or (
            preserved_extension == lowercase_extension
            and lowercase_extension in allowed_lower
        )
    ):
        allowed_ordered = tuple(extension.lstrip(".") for extension in allowed)
        allowed_text = ", ".join(
            f".{extension}" for extension in allowed_ordered
        )

        raise ValueError(
            f"Invalid extension '{preserved_extension}'; allowed: "
            f"{allowed_text}.",
        )

    output_path = (
        f"{stem}.{preserved_extension}.gz"
        if is_compressed
        else f"{stem}.{preserved_extension}"
    )

    extension = (
        "bedGraph"
        if preserved_extension == "bedGraph"
        else lowercase_extension
    )

    return output_path, extension, is_compressed


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
        When 'kind="file"' and the file already exists, require the file
        itself to be writable (default: False).
    label : str | None
        Optional label for nicer error text.

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
