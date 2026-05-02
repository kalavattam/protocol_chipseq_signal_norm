#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Distributed under the MIT license.

"""
Script
------
add_coeffs_namespaced.py


Description
-----------
Read a tab-delimited coefficients table (optionally gzip-compressed), compute
reference statistics across all samples for a small set of coefficient columns,
and write an augmented table containing the following:

    (1) the original per-sample coefficient values (namespaced), and
    (2) per-sample values divided by each reference statistic (e.g., x/min,
        x/median, x/mean, ...).

This is intended to produce a “wide” table of per-sample scaling factors that
can be consumed downstream (e.g., to scale raw bedGraph tracks).


Usage
-----
python add_coeffs_namespaced.py \\
    <infile.tsv[.gz]> [--inplace | -o/--outfile <out.tsv[.gz]>] [-d/--dp <dp>]


Arguments
---------
infile
    Input TSV path ('.tsv' or '.tsv.gz'). Must be headered and tab-delimited.

-o, --outfile
    Output TSV path (supports '.gz'). If omitted and '--inplace' is not used,
    output is written to stdout.

--inplace
    Overwrite the input file. The script writes to a temporary file first and
    then atomically replaces the original input path.

-dp, --dp, --rnd, --round, --decimals, --digits
    Decimal precision for emitted numeric fields (default: 24).


Output
------
A headered TSV (or TSV.GZ) with columns:

    sample
    pair_s, pair_s_min, pair_s_median, pair_s_mean, pair_s_gmean,
        pair_s_hmean, pair_s_max
    pair_alpha_rxinput, pair_alpha_rxinput_min, ...
    ip_alpha_ip, ip_alpha_ip_min, ...
    in_alpha_in, in_alpha_in_min, ...

For each coefficient name, the script emits the raw value (formatted with
'--dp' decimal places), followed by values normalized by each reference
statistic (min/median/mean/gmean/hmean/max). If a reference statistic is
unavailable (e.g., geometric mean with non-positive values), that column is
emitted as an empty string for all rows.


Examples
--------
(1) Write augmented TSV to stdout:
'''bash
python add_coeffs_namespaced.py \\
    coefficients.Brn1.tsv > coefficients.Brn1_plus.tsv
'''

(2) Write gzip-compressed output:
'''bash
python add_coeffs_namespaced.py \\
    coefficients.Hmo1.tsv \\
    --outfile coefficients.Hmo1_plus.tsv.gz
'''

(3) Overwrite the input file in-place:
'''bash
python add_coeffs_namespaced.py coefficients.all.tsv --inplace
'''


General notes
-------------
- Input must be tab-delimited with a header row; all data rows must have the
  same number of fields as the header.
- The script requires a 'sample' column (used as the per-row identifier).
- Coefficient columns are chosen via a namespaced-first policy:
      + pair_s             (fallback: s)
      + pair_alpha_rxinput (fallback: alpha_rxinput)
      + ip_alpha_ip        (fallback: alpha_ip)
      + in_alpha_in        (fallback: alpha_in)
- Output column order is grouped “like-with-like”: for each coefficient name,
  the raw value is emitted first, then its /min, /median, /mean, /gmean,
  /hmean, and /max normalized variants.
- Numeric formatting uses fixed-point with '--dp' decimal places (default: 24).


Performance notes
-----------------
- The script loads the entire TSV into memory (rows as dicts) and computes
  reference statistics across all samples; memory use scales with number of
  rows × number of columns.
- For typical coefficient tables (tens to thousands of samples), runtime is
  dominated by parsing and string formatting rather than the statistics.
- In-place mode uses a temporary file plus an atomic rename (os.replace) to
  avoid leaving a partially written file on failure.


#TODO
-----
- Integrate with already-written modules (can delete/move many functions
  below).
"""

from __future__ import annotations

import argparse
import gzip
import math
import os
import sys
from typing import Dict, List


def open_maybe_gz_in(path: str):
    """
    Open a text input file, transparently handling .gz paths.
    """
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")

    return open(path, "rt", encoding="utf-8", newline="")


def open_maybe_gz_out(path: str):
    """
    Open a text output file, transparently handling .gz paths.
    """
    if path.endswith(".gz"):
        return gzip.open(path, "wt", encoding="utf-8", newline="")

    return open(path, "wt", encoding="utf-8", newline="")


def median(xs: List[float]) -> float:
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


def gmean(xs: List[float]) -> float:
    """
    Return the geometric mean of positive floats.
    """
    if any(x <= 0 for x in xs):
        raise ValueError("gmean() requires all values > 0")

    return math.exp(sum(math.log(x) for x in xs) / len(xs))


def hmean(xs: List[float]) -> float:
    """
    Return the harmonic mean of positive floats.
    """
    if any(x <= 0 for x in xs):
        raise ValueError("hmean() requires all values > 0")

    return len(xs) / sum(1.0 / x for x in xs)


def safe_refs(xs: List[float], label: str) -> Dict[str, float]:
    """
    Compute reference stats (min, median, mean, and max, plus optional gmean
    and hmean).
    """
    refs: Dict[str, float] = {}
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


def parse_args(argv: List[str] | None = None) -> argparse.Namespace:
    """
    Parse CLI args for reading a coefficients TSV and writing an augmented TSV.
    """
    p = argparse.ArgumentParser()
    p.add_argument("infile", help="coefficients.all.tsv[.gz]")
    p.add_argument(
        "-o", "--outfile",
        default=None,
        help="Output path (supports .gz). If omitted, write to stdout."
    )
    p.add_argument(
        "--inplace",
        action="store_true",
        help="Overwrite infile (writes a temp file then replaces)."
    )
    p.add_argument(
        "-dp", "--dp", "--rnd", "--round", "--decimals", "--digits",
        type=int,
        default=24,
        help=(
            "Decimal precision for emitted numeric fields (default: "
            "%(default)s)."
        )
    )
    return p.parse_args(argv)


def main(argv: List[str] | None = None) -> int:
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
    in_path = args.infile

    #  Disallow ambiguous output mode: can't request both inplace and outfile
    if args.inplace and args.outfile is not None:
        raise ValueError("Use either '--inplace' or '--outfile', not both.")

    #  If no outfile and not inplace, write to stdout
    write_stdout = (not args.inplace) and (args.outfile is None)

    def parse_tsv_rows(handle) -> tuple[List[str], List[Dict[str, str]]]:
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
        rows: List[Dict[str, str]] = []
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
            rows.append(dict(zip(fields, parts)))

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

    def col_as_floats(col: str) -> List[float]:
        """
        Extract one column from rows as floats (whitespace-trimmed).
        """
        out: List[float] = []
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
    new_fields: List[str] = []
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
    out_fields = ["sample"] + new_fields

    #  Cache decimal precision and build a reusable fixed-point format string
    dp = args.dp
    fmt = f"{{:.{dp}f}}"

    def adj(x: float, ref: Dict[str, float], key: str) -> str:
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
        handle, fields: List[str], out_rows: List[Dict[str, str]]
    ) -> None:
        """
        Write header and rows as TSV in 'fields' order (missing keys: "").
        """
        handle.write("\t".join(fields) + "\n")
        for rr in out_rows:
            handle.write("\t".join(rr.get(f, "") for f in fields) + "\n")

    #  Build output rows: one output dict per input sample
    out_rows: List[Dict[str, str]] = []
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
        out_row: Dict[str, str] = {"sample": r["sample"]}
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

    #  Write to stdout or to a file, optionally via a temp file for inplace
    if write_stdout:
        write_tsv(sys.stdout, out_fields, out_rows)
    else:
        #  Pick output path: inplace uses a temp file then atomic replace
        if args.inplace:
            if in_path.endswith(".gz"):
                tmp_path = f"{in_path}.tmp.gz"
            else:
                tmp_path = f"{in_path}.tmp"
            out_path = in_path
        else:
            tmp_path = args.outfile  # args.outfile guaranteed non-'None' here
            out_path = args.outfile

        #  Write the TSV to the chosen destination (gzip if tmp_path endswith
        #  .gz)
        with open_maybe_gz_out(tmp_path) as out_fh:
            write_tsv(out_fh, out_fields, out_rows)

        #  If inplace mode, atomically replace the original with the temp file
        if args.inplace:
            os.replace(tmp_path, out_path)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
