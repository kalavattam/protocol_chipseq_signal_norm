#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: add_coeffs_namespaced.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.4, GPT-5.5, GPT-5.6) were used in development
# and documentation.
#
# Distributed under the MIT license.


"""
Augment coefficient tables with namespaced raw and normalized columns.

Usage
-----
python add_coeffs_namespaced.py <file> [--in_place | --fil_out <file>]

Parameters
----------
Input table path, output path, in-place mode, and decimal precision are parsed
by 'parse_args()'.

Returns
-------
Writes a headered TSV or TSV.GZ containing original coefficient values and
values divided by reference statistics.

Notes
-----
The script loads the coefficient table into memory and computes reference
statistics across samples.
"""

from __future__ import annotations

import argparse
import gzip
import math
import os
import sys

from protocol_chipseq_signal_norm.utilities.utils_cli import (
    CapArgumentParser,
    add_help_cap,
)


# TODO: Integrate this script with the shared table/formatting modules where
#       possible, then delete or move duplicated helper functions below.
def open_maybe_gz_in(path: str):
    """
    Open a text input file, transparently handling .gz paths.
    """
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")

    return open(path, encoding="utf-8", newline="")


def open_maybe_gz_out(path: str):
    """
    Open a text output file, transparently handling .gz paths.
    """
    if path.endswith(".gz"):
        return gzip.open(path, "wt", encoding="utf-8", newline="")

    return open(path, "w", encoding="utf-8", newline="")


def median(xs: list[float]) -> float:
    """
    Return the median of a non-empty list of floats.
    """
    ys = sorted(xs)
    n = len(ys)
    if n == 0:
        raise ValueError("median(): empty list")
    mid = n // 2
    if n % 2 == 1:
        return ys[mid]

    return 0.5 * (ys[mid - 1] + ys[mid])


def gmean(xs: list[float]) -> float:
    """
    Return the geometric mean of positive floats.
    """
    if any(x <= 0 for x in xs):
        raise ValueError("gmean() requires all values > 0")

    return math.exp(sum(math.log(x) for x in xs) / len(xs))


def hmean(xs: list[float]) -> float:
    """
    Return the harmonic mean of positive floats.
    """
    if any(x <= 0 for x in xs):
        raise ValueError("hmean() requires all values > 0")

    return len(xs) / sum(1.0 / x for x in xs)


def safe_refs(xs: list[float], label: str) -> dict[str, float]:
    """
    Compute reference stats (min, median, mean, and max, plus optional gmean
    and hmean).
    """
    refs: dict[str, float] = {}
    refs["min"] = min(xs)
    refs["median"] = median(xs)
    refs["mean"] = sum(xs) / len(xs)
    refs["max"] = max(xs)

    try:
        refs["gmean"] = gmean(xs)
    except ValueError as e:
        print(f"warning: {label}: gmean unavailable: {e}", file=sys.stderr)

    try:
        refs["hmean"] = hmean(xs)
    except ValueError as e:
        print(f"warning: {label}: hmean unavailable: {e}", file=sys.stderr)

    return refs


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse CLI args for reading a coefficients TSV and writing an augmented TSV.
    """
    p = CapArgumentParser()
    add_help_cap(p)
    p.add_argument("fil_in", help="Input file path. coefficients.all.tsv[.gz]")
    p.add_argument(
        "-fo", "--fil_out",
        dest="fil_out",
        default=None,
        help=(
            "Output file path. Output path (supports .gz). If omitted, write "
            "to stdout."
        )
    )
    p.add_argument(
        "-ip", "--in_place",
        dest="in_place",
        action="store_true",
        default=False,
        help="Overwrite fil_in (writes a temp file then replaces)."
    )
    p.add_argument(
        "-dp", "--dp",
        dest="dp",
        type=int,
        default=24,
        help=(
            "Maximum number of decimal places retained for finite emitted "
            "values. Decimal precision for emitted numeric fields (default: "
            "%(default)s)."
        )
    )
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """
    Read coefficient rows, add namespaced columns, and return an exit status.
    """

    #  Resolve argv: default to CLI args (sys.argv[1:]) unless caller supplied
    argv_parse = sys.argv[1:] if argv is None else argv
    if not argv_parse:
        #  No args: show help and exit cleanly
        parse_args(["-h"])
        return 0

    #  Parse CLI options into an argparse.Namespace
    args = parse_args(argv_parse)

    #  Validate assignment to args.dp
    if args.dp < 0:
        raise ValueError("'--dp' must be >= 0.")

    #  Cache input path for readability, reuse
    in_path = args.fil_in

    #  Disallow ambiguous output mode: can't request both in-place and fil_out
    if args.in_place and args.fil_out is not None:
        raise ValueError("Use either '--in_place' or '--fil_out', not both.")

    #  If no fil_out and not in-place, write to stdout
    write_stdout = (not args.in_place) and (args.fil_out is None)

    def parse_tsv_rows(handle) -> tuple[list[str], list[dict[str, str]]]:
        """
        Parse a headered TSV into (fields, rows) dicts; rejects malformed rows.
        """
        #  Read the header line (first row) to get column names
        header_line = handle.readline()
        if not header_line:
            #  Error if empty file (or missing header)
            raise ValueError("Input has no header.")

        #  Strip trailing newline(s) so split() doesn’t keep them
        header_line = header_line.rstrip("\n").rstrip("\r")

        #  Split header into field names on tabs
        fields = header_line.split("\t")

        #  Do header validation: must have names, none can be empty
        if not fields or any(f == "" for f in fields):
            raise ValueError("Invalid header: empty field name found.")

        #  Accumulate data rows as dicts keyed by header field names
        rows: list[dict[str, str]] = []
        for line in handle:
            #  Strip trailing newline(s)
            line = line.rstrip("\n").rstrip("\r")

            #  Skip blank lines (defensive)
            if line == "":
                continue

            #  Split row into tab-separated fields
            parts = line.split("\t")

            #  Enforce rectangular TSV: every row must have exactly len(fields)
            #  cols
            if len(parts) != len(fields):
                raise ValueError(
                    f"Malformed TSV row: expected {len(fields)} fileds, but "
                    f"got {len(parts)}: {line!r}."
                )

            #  Convert the row into a dict: {field_name: value_string, ...}
            rows.append(dict(zip(fields, parts, strict=True)))

        #  Return header fields and all parsed row dicts
        return fields, rows

    #  Read and parse the entire input TSV into memory (header and row dicts)
    with open_maybe_gz_in(in_path) as fh:
        in_fields, rows = parse_tsv_rows(fh)

    #  Convert header fields to a set for fast membership checks
    fns = set(in_fields)

    #  Ensure TSV has column 'sample', which is needed for indexing (row keys)
    if "sample" not in fns:
        raise ValueError("Missing required column: 'sample'.")

    #  Prefer namespaced columns if present, else fall back to legacy names
    col_s = "pair_s" if "pair_s" in fns else "s"

    #  Choose rxinput column name (namespaced preferred; legacy fallback)
    if "pair_alpha_rxinput" in fns:
        col_arx = "pair_alpha_rxinput"
    else:
        col_arx = "alpha_rxinput"

    #  Choose IP/input alpha columns (namespaced preferred; legacy fallback)
    col_aip = "ip_alpha_ip" if "ip_alpha_ip" in fns else "alpha_ip"
    col_ain = "in_alpha_in" if "in_alpha_in" in fns else "alpha_in"

    #  Fail fast if any required coefficient column is missing
    for col in (col_s, col_arx, col_aip, col_ain):
        if col not in fns:
            raise ValueError(f"Missing required column: '{col}'.")

    def col_as_floats(col: str) -> list[float]:
        """
        Extract one column from rows as floats (whitespace-trimmed).
        """
        out: list[float] = []
        for r in rows:
            v = r[col].strip()
            out.append(float(v))
        return out

    #  Pull each coefficient column as a float vector (one value per sample)
    s_vals = col_as_floats(col_s)
    arx_vals = col_as_floats(col_arx)
    aip_vals = col_as_floats(col_aip)
    ain_vals = col_as_floats(col_ain)

    #  Compute reference stats per coefficient (min, median, means, and max,
    #  etc.)
    refs_pair_s = safe_refs(s_vals, "pair_s")
    refs_pair_arx = safe_refs(arx_vals, "pair_alpha_rxinput")
    refs_ip_aip = safe_refs(aip_vals, "ip_alpha_ip")
    refs_in_ain = safe_refs(ain_vals, "in_alpha_in")

    #  Define which reference stats to emit (and in what order) per coefficient
    ref_keys = ["min", "median", "mean", "gmean", "hmean", "max"]

    #  Build the output header fields in “like-with-like” grouped order
    new_fields: list[str] = []
    for name in (
        "pair_s",
        "pair_alpha_rxinput",
        "ip_alpha_ip",
        "in_alpha_in",
    ):
        #  Raw coefficient first, then its per-stat scaled variants
        new_fields.append(name)
        for k in ref_keys:
            new_fields.append(f"{name}_{k}")

    #  Final output header: sample id + all generated fields
    out_fields = ["sample", *new_fields]

    #  Cache decimal precision and build a reusable fixed-point format string
    dp = args.dp
    fmt = f"{{:.{dp}f}}"

    def adj(x: float, ref: dict[str, float], key: str) -> str:
        """
        Return x/ref[key] as a formatted string, or "" if key missing/invalid.
        """
        if key not in ref:
            return ""
        r = ref[key]
        if r == 0:
            return ""
        return fmt.format(x / r)

    def write_tsv(
        handle, fields: list[str], out_rows: list[dict[str, str]]
    ) -> None:
        """
        Write header and rows as TSV in 'fields' order (missing keys: "").
        """
        handle.write("\t".join(fields) + "\n")
        for rr in out_rows:
            handle.write("\t".join(rr.get(f, "") for f in fields) + "\n")

    #  Build output rows: one output dict per input sample
    out_rows: list[dict[str, str]] = []
    for r in rows:
        #  Parse the four coefficient values for given sample row
        s = float(r[col_s])
        arx = float(r[col_arx])
        aip = float(r[col_aip])
        ain = float(r[col_ain])

        #  Map coefficient names to values and to reference-stat dicts
        vals = {
            "pair_s": s,
            "pair_alpha_rxinput": arx,
            "ip_alpha_ip": aip,
            "in_alpha_in": ain,
        }
        refs = {
            "pair_s": refs_pair_s,
            "pair_alpha_rxinput": refs_pair_arx,
            "ip_alpha_ip": refs_ip_aip,
            "in_alpha_in": refs_in_ain,
        }

        #  Fill one output row: sample id, raw coeffs, scaled-by-stat values
        out_row: dict[str, str] = {"sample": r["sample"]}
        for name in (
            "pair_s",
            "pair_alpha_rxinput",
            "ip_alpha_ip",
            "in_alpha_in",
        ):
            out_row[name] = fmt.format(vals[name])
            for k in ref_keys:
                out_row[f"{name}_{k}"] = adj(vals[name], refs[name], k)

        #  Accumulate row dicts for writing
        out_rows.append(out_row)

    #  Write to stdout or to a file, optionally via a temp file for in-place
    if write_stdout:
        write_tsv(sys.stdout, out_fields, out_rows)
    else:
        #  Pick output path: in-place uses a temp file then atomic replace
        if args.in_place:
            if in_path.endswith(".gz"):
                tmp_path = f"{in_path}.tmp.gz"
            else:
                tmp_path = f"{in_path}.tmp"
            out_path = in_path
        else:
            tmp_path = args.fil_out  # args.fil_out guaranteed non-'None' here
            out_path = args.fil_out

        #  Write the TSV to the chosen destination (gzip if tmp_path endswith
        #  .gz)
        with open_maybe_gz_out(tmp_path) as out_fh:
            write_tsv(out_fh, out_fields, out_rows)

        #  If in-place mode, atomically replace the original with the temp file
        if args.in_place:
            os.replace(tmp_path, out_path)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
