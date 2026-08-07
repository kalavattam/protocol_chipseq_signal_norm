#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: add_coeffs_namespaced.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.4, GPT-5.5, GPT-5.6) were used in design,
# development, and documentation, with all output reviewed, edited, and
# approved by the author.
#
# Distributed under the MIT license.


"""
Augment coefficient tables with namespaced raw and normalized columns.

The CLI accepts input and output paths, in-place mode, and decimal precision.
It writes a headered TSV or TSV.GZ containing original coefficient values and
values divided by reference statistics.

Notes
-----
The script loads the coefficient table into memory and computes reference
statistics across samples.

Examples
--------
python add_coeffs_namespaced.py <file> [--in_place | --fil_out <file>]
"""

from __future__ import annotations

import argparse
import gzip
import math
import os
import sys
from typing import TextIO

from protocol_chipseq_signal_norm.utilities.utils_cli import (
    CapArgumentParser,
    add_help_cap,
)

REFERENCE_KEYS = ("min", "median", "mean", "gmean", "hmean", "max")
COEFFICIENT_NAMES = (
    "pair_s",
    "pair_alpha_rxinput",
    "ip_alpha_ip",
    "in_alpha_in",
)


# TODO: Integrate this script with the shared table/formatting modules where
# possible, then delete or move duplicated helper functions below.
def open_maybe_gzip_input(path: str) -> TextIO:
    """
    Open a text input file, transparently handling .gz paths.
    """

    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")

    return open(path, encoding="utf-8", newline="")


def open_maybe_gzip_output(path: str) -> TextIO:
    """
    Open a text output file, transparently handling .gz paths.
    """

    if path.endswith(".gz"):
        return gzip.open(path, "wt", encoding="utf-8", newline="")

    return open(path, "w", encoding="utf-8", newline="")


def median(values: list[float]) -> float:
    """
    Return the median of a non-empty list of floats.
    """

    ordered_values = sorted(values)
    count = len(ordered_values)
    if count == 0:
        raise ValueError("median(): empty list")

    midpoint = count // 2
    if count % 2 == 1:
        return ordered_values[midpoint]

    return 0.5 * (ordered_values[midpoint - 1] + ordered_values[midpoint])


def gmean(values: list[float]) -> float:
    """
    Return the geometric mean of positive floats.
    """

    if any(value <= 0 for value in values):
        raise ValueError("gmean() requires all values > 0")

    return math.exp(
        sum(math.log(value) for value in values) / len(values),
    )


def hmean(values: list[float]) -> float:
    """
    Return the harmonic mean of positive floats.
    """

    if any(value <= 0 for value in values):
        raise ValueError("hmean() requires all values > 0")

    return len(values) / sum(1.0 / value for value in values)


def safe_refs(values: list[float], label: str) -> dict[str, float]:
    """
    Return conventional reference statistics for one float series.
    """

    refs: dict[str, float] = {}
    refs["min"] = min(values)
    refs["median"] = median(values)
    refs["mean"] = sum(values) / len(values)
    refs["max"] = max(values)

    try:
        refs["gmean"] = gmean(values)
    except ValueError as e:
        print(f"warning: {label}: gmean unavailable: {e}", file=sys.stderr)

    try:
        refs["hmean"] = hmean(values)
    except ValueError as e:
        print(f"warning: {label}: hmean unavailable: {e}", file=sys.stderr)

    return refs


def _parse_tsv_rows(
    handle: TextIO,
) -> tuple[list[str], list[dict[str, str]]]:
    """
    Parse a headered TSV and reject malformed rows.

    Parameters
    ----------
    handle : TextIO
        Text stream positioned at the TSV header.

    Returns
    -------
    fields, rows : tuple[list[str], list[dict[str, str]]]
        Header fields and nonblank data rows keyed by those fields.

    Raises
    ------
    ValueError
        If the header is absent or a data row has the wrong field count.
    """

    header_line = handle.readline()

    if not header_line:
        raise ValueError("Input has no header.")

    header_line = header_line.rstrip("\n").rstrip("\r")
    fields = header_line.split("\t")

    if not fields or any(field == "" for field in fields):
        raise ValueError("Invalid header: empty field name found.")

    rows: list[dict[str, str]] = []

    for line in handle:
        line = line.rstrip("\n").rstrip("\r")

        if line == "":
            continue

        parts = line.split("\t")

        if len(parts) != len(fields):
            raise ValueError(
                f"Malformed TSV row: expected {len(fields)} fileds, but "
                f"got {len(parts)}: {line!r}.",
            )

        rows.append(dict(zip(fields, parts, strict=True)))

    return fields, rows


def _coefficient_columns(fields: list[str]) -> dict[str, str]:
    """
    Resolve canonical coefficient names to input columns.
    """

    field_names = set(fields)

    if "sample" not in field_names:
        raise ValueError("Missing required column: 'sample'.")

    columns = {
        "pair_s": "pair_s" if "pair_s" in field_names else "s",
        "pair_alpha_rxinput": (
            "pair_alpha_rxinput"
            if "pair_alpha_rxinput" in field_names
            else "alpha_rxinput"
        ),
        "ip_alpha_ip": (
            "ip_alpha_ip" if "ip_alpha_ip" in field_names else "alpha_ip"
        ),
        "in_alpha_in": (
            "in_alpha_in" if "in_alpha_in" in field_names else "alpha_in"
        ),
    }

    for column in columns.values():
        if column not in field_names:
            raise ValueError(f"Missing required column: '{column}'.")

    return columns


def _column_as_floats(
    rows: list[dict[str, str]],
    column: str,
) -> list[float]:
    """
    Return one whitespace-trimmed input column as floats.
    """

    return [float(row[column].strip()) for row in rows]


def _reference_statistics(
    rows: list[dict[str, str]],
    columns: dict[str, str],
) -> dict[str, dict[str, float]]:
    """
    Calculate available reference statistics for every coefficient family.
    """

    return {
        name: safe_refs(_column_as_floats(rows, columns[name]), name)
        for name in COEFFICIENT_NAMES
    }


def _output_fields() -> list[str]:
    """
    Return canonical output fields in rendering order.
    """

    fields = ["sample"]

    for name in COEFFICIENT_NAMES:
        fields.append(name)
        fields.extend(f"{name}_{key}" for key in REFERENCE_KEYS)

    return fields


def _adjusted_value(
    value: float,
    references: dict[str, float],
    key: str,
    format_spec: str,
) -> str:
    """
    Format one reference-adjusted value or an empty unavailable field.
    """

    reference = references.get(key)
    if reference is None or reference == 0:
        return ""

    return format_spec.format(value / reference)


def _augmented_rows(
    rows: list[dict[str, str]],
    columns: dict[str, str],
    references: dict[str, dict[str, float]],
    decimal_places: int,
) -> list[dict[str, str]]:
    """
    Return namespaced raw and reference-adjusted output rows.
    """

    format_spec = f"{{:.{decimal_places}f}}"
    output_rows: list[dict[str, str]] = []

    for row in rows:
        values = {
            name: float(row[columns[name]]) for name in COEFFICIENT_NAMES
        }
        output_row: dict[str, str] = {"sample": row["sample"]}

        for name in COEFFICIENT_NAMES:
            output_row[name] = format_spec.format(values[name])

            for key in REFERENCE_KEYS:
                output_row[f"{name}_{key}"] = _adjusted_value(
                    values[name],
                    references[name],
                    key,
                    format_spec,
                )

        output_rows.append(output_row)

    return output_rows


def _write_tsv(
    handle: TextIO,
    fields: list[str],
    rows: list[dict[str, str]],
) -> None:
    """
    Write headered rows as a field-ordered TSV.
    """

    handle.write("\t".join(fields) + "\n")

    for row in rows:
        handle.write("\t".join(row.get(field, "") for field in fields) + "\n")


def _write_output(
    args: argparse.Namespace,
    input_path: str,
    fields: list[str],
    rows: list[dict[str, str]],
) -> None:
    """
    Write an augmented table to stdout, a selected path, or the input path.
    """

    if not args.in_place and args.fil_out is None:
        _write_tsv(sys.stdout, fields, rows)

        return

    output_path = input_path if args.in_place else args.fil_out
    if output_path is None:
        raise AssertionError("an output path is required outside stdout mode")

    if args.in_place:
        temporary_path = (
            f"{input_path}.tmp.gz"
            if input_path.endswith(".gz")
            else f"{input_path}.tmp"
        )
    else:
        temporary_path = output_path

    with open_maybe_gzip_output(temporary_path) as handle:
        _write_tsv(handle, fields, rows)

    if args.in_place:
        os.replace(temporary_path, output_path)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse CLI args for reading a coefficients TSV and writing an augmented TSV.

    Parameters
    ----------
    argv : list[str] | None
        Explicit arguments, or None to read the process arguments.

    Returns
    -------
    arguments : argparse.Namespace
        Parsed coefficient-input and output options.

    Raises
    ------
    SystemExit
        If argument parsing fails or help is requested.
    """

    p = CapArgumentParser()
    add_help_cap(p)
    p.add_argument(
        "fil_in",
        help="Input file path: coefficients.all.tsv[.gz].",
    )
    p.add_argument(
        "-fo",
        "--fil_out",
        dest="fil_out",
        default=None,
        help=(
            "Output file path. Output path (supports .gz). If omitted, write "
            "to stdout."
        ),
    )
    p.add_argument(
        "-ip",
        "--in_place",
        dest="in_place",
        action="store_true",
        default=False,
        help="Overwrite fil_in (writes a temp file then replaces).",
    )
    p.add_argument(
        "-dp",
        "--dp",
        dest="dp",
        type=int,
        default=24,
        help=(
            "Maximum number of decimal places retained for finite emitted "
            "values. Decimal precision for emitted numeric fields (default: "
            "%(default)s)."
        ),
    )

    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """
    Read coefficient rows, add namespaced columns, and return an exit status.

    Parameters
    ----------
    argv : list[str] | None
        Explicit arguments, or None to read the process arguments.

    Returns
    -------
    status : int
        Zero after the augmented table is written successfully.

    Raises
    ------
    ValueError
        If the coefficient table or requested reference data are invalid.
    """

    argv_parse = sys.argv[1:] if argv is None else argv

    if not argv_parse:
        parse_args(["-h"])

        return 0

    args = parse_args(argv_parse)

    if args.dp < 0:
        raise ValueError("'--dp' must be >= 0.")

    in_path = args.fil_in

    if args.in_place and args.fil_out is not None:
        raise ValueError("Use either '--in_place' or '--fil_out', not both.")

    with open_maybe_gzip_input(in_path) as handle:
        input_fields, rows = _parse_tsv_rows(handle)

    columns = _coefficient_columns(input_fields)
    references = _reference_statistics(rows, columns)

    output_fields = _output_fields()
    output_rows = _augmented_rows(rows, columns, references, args.dp)

    _write_output(args, in_path, output_fields, output_rows)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
