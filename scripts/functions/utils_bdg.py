#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: utils_bdg.py
#
# Copyright 2025-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5-series models) was used in development.
#
# Distributed under the MIT license.

"""
Script
------
utils_bdg.py


Description
-----------
Helper functions for bedGraph-style text parsing and numeric handling.


Functions
---------
canon_nonfinite()
    Canonicalize NaN/±inf spellings to {'nan','inf','-inf'}.

try_float()
    Return float() with None on failure.

iter_bdg_rows()
    Stream bedGraph rows, yielding '(chrom, start, end, val_tok, val_num)',
    where 'val_num' is None for non-finite tokens (via 'canon_nonfinite') or
    unparsable values; otherwise, 'val_num' is float (via 'try_float').

check_bin_size()
    Compare the first few data bins in two BED or bedGraph files and error if
    their bin sizes differ (uses 'open_in' and 'read_data_line').

key_bin()
    Sort key for binned signals: first by yeast-aware chromosome order (via
    'sort_chrom'), then by bin start.

write_bdg()
    Stream bedGraph output from a dict of {(chrom, bin_start): value}
    using 'open_out' and 'key_bin'; rounds to at most 'rnd' decimals and
    strips non-informative trailing zeros.

generate_track_name()
    Insert '.track' before the main extension of an output filename,
    preserving any trailing '.gz'. Current use: '.track' sidecars optionally
    output by 'compute_signal_ratio.py'.


Example
-------
'''python
from scripts.functions.utils_bdg import (
    iter_bdg_rows, canon_nonfinite, try_float
)
from scripts.functions.utils_io import is_header, open_in

with open_in("track.bdg.gz") as fh:
    skp_prd = lambda line: is_header(line, ("#", "track", "browser"))
    for chrom, s, e, tok, num in iter_bdg_rows(fh, skp_prd=skp_prd):
        ...
'''


Notes
-----
- BED and bedGraph are treated as 0-based, half-open intervals: [start, end).
- Header detection/skip is caller-controlled via 'skp_prd'; see
  'utils_io.is_header'.
- Sorting uses 'utils_chrom.sort_chrom' to respect S. cerevisiae/S. pombe
  Roman-numeral ordering.
- I/O is line-streamed; no clamping to chromosome bounds is performed.
- bedGraph numeric formatting uses at most the requested decimal precision:
  trailing zeros after the decimal point, and a trailing decimal point itself,
  are stripped from finite values.
- Non-finite values ('nan', 'inf', '-inf') are surfaced via 'canon_nonfinite';
  'iter_bdg_rows' reports their numeric as None (caller decides how to handle).
"""

from __future__ import annotations

from .utils_chrom import sort_chrom
from .utils_io import DEF_SKP_PFX, open_in, open_out, read_data_line
from typing import Iterator, Callable

import os
import sys

assert sys.version_info >= (3, 10), "Python >= 3.10 required."


def canon_nonfinite(tok: str) -> str | None:
    """
    Return a canonical non-finite token or None.

    Args:
        tok : str
            Raw string token.

    Returns:
        'nan' | 'inf' | '-inf'
            If the token case-insensitively matches any NaN/±inf spelling
            (including '+inf'); otherwise None.

    Raises:
        None.
    """
    t = tok.strip().lower()
    if t in ("nan", "inf", "-inf", "+inf"):
        if t == "nan":
            return "nan"
        else:
            if t.startswith("-"):
                return "-inf"
            else:
                return "inf"

    return None


def try_float(tok: str) -> float | None:
    """
    Parse a float token.

    Args:
        tok : str
            Raw string token.

    Returns:
        float | None
            On success, the parsed token as a float; on failure, None.

    Raises:
        None.

    Notes:
        Non-finite tokens like 'nan'/'inf' are not canonicalized here; combine
        with 'canon_nonfinite' to treat them specially.
    """
    try:
        return float(tok)
    except ValueError:
        return None


def iter_bdg_rows(
    stream, skp_prd: Callable[[str], bool]
) -> Iterator[tuple[str, int, int, str, float | None]]:
    """
    Stream bedGraph-like rows.

    Yields:
        '(chrom, start, end, val_tok, val_num)', where 'val_num' is None for
        non-finite tokens or unparsable values.

    Raises:
        None.

    Notes:
        - Skips lines where 'skp_prd(line)' is True (headers/metadata/blank).
        - Skips lines with < 4 fields, non-integer coords, or non-positive
          width.
        - Coordinates are not clamped; caller may clamp to chrom bounds.
    """
    for raw in stream:
        if skp_prd(raw):
            continue

        parts = raw.rstrip("\n").split()
        if len(parts) < 4:
            continue

        chrom, s_str, e_str, v_str = parts[:4]

        try:
            s = int(s_str)
            e = int(e_str)
        except ValueError:
            continue

        if e <= s:
            continue

        nf = canon_nonfinite(v_str)
        v_num = None if nf is not None else try_float(v_str)

        yield chrom, s, e, v_str, v_num


def check_bin_size(
    fil_A: str,
    fil_B: str,
    skp_pfx: tuple[str, ...] = DEF_SKP_PFX
) -> None:
    """
    Verify that the bin sizes in two BED or bedGraph files are identical.

    Args:
        fil_A : str
            Path to first BED or bedGraph input file (file A; e.g., IP).
        fil_B : str
            Path to second BED or bedGraph input file (file B; e.g., input).
        skp_pfx : tuple[str, ...]
            Prefixes to skip as BED or bedGraph header/meta lines.

    Returns:
        None.

    Raises:
        ValueError
            If the bin sizes differ, if lines are malformed, or if positions
            are not parseable as integers.

    Notes:
        Refactored out of 'compute_signal_ratio.py' for modularization.
    """
    #  Compare the first few bins from each file to ensure they use the same
    #  bin size (if not, fast fail)
    with open_in(fil_A) as opn_A, open_in(fil_B) as opn_B:
        seen = 0
        while seen < 5:
            lin_A = read_data_line(fh=opn_A, skp_pfx=skp_pfx)
            lin_B = read_data_line(fh=opn_B, skp_pfx=skp_pfx)
            if not lin_A or not lin_B:
                break

            fld_A = lin_A.split()
            fld_B = lin_B.split()

            #  Fail-fast: Require at least chrom, start, end
            if len(fld_A) < 3:
                raise ValueError(
                    "Malformed BED or bedGraph line in first file (file A; "
                    f"e.g., IP) during bin-size check: {lin_A!r}"
                )
            if len(fld_B) < 3:
                raise ValueError(
                    "Malformed BED or bedGraph line in second file (file B; "
                    f"e.g., input) during bin-size check: {lin_B!r}"
                )

            try:
                bin_A = int(fld_A[2]) - int(fld_A[1])
                bin_B = int(fld_B[2]) - int(fld_B[1])
            except ValueError as e:
                raise ValueError(
                    "Non-integer start/end encountered during bin-size "
                    "check:\n"
                    f"  First file (file A; e.g., IP):     {lin_A!r}\n"
                    f"  Second file (file B; e.g., input): {lin_B!r}"
                ) from e

            if bin_A != bin_B:
                raise ValueError(
                    f"Error: Mismatched bin sizes detected.\n"
                    f"  - {fil_A}: {bin_A}-bp bins\n"
                    f"  - {fil_B}: {bin_B}-bp bins\n\n"
                    f"Ensure both files are binned identically."
                )
            seen += 1


def key_bin(chrom: str, start: int) -> tuple[tuple[int, int, str], int]:
    #  Use Roman-order first via sort_chrom(), then start coordinate
    #  (refactored out of 'compute_signal_ratio.py' for modularization)
    return (sort_chrom(chrom), start)


def write_bdg(cvg, fil_out, siz_bin, rnd):
    """
    Write binned signal data to a bedGraph file.

    Args:
        cvg : dict
            Binned signal data, where keys are (chrom, bin_start) and float
            values.
        fil_out: str
            Path to the output file: '.bedGraph', '.bedgraph', '.bdg', or
            '.bg', optionally with '.gz'.
        siz_bin: int
            Bin size in base pairs.
        rnd : int
            Maximum number of decimal places retained for signal values.

    Returns:
        None. Writes interval records to the output file.

    Raises:
        ValueError
            If 'siz_bin' <= 0 or 'rnd' < 0.
        OSError
            On filesystem/permission issues while writing.

    Notes:
        - Signal values are rounded to at most 'rnd' decimal places.
        - After rounding, non-informative trailing zeros are stripped from
          finite decimal representations, and any trailing decimal point is
          removed.
        - Negative zero is emitted as '0'.
        - Output is automatically gzip-compressed if 'fil_out' ends with '.gz'.
        - bedGraph format: 'chrom<tab>start<tab>end<tab>signal<newline>'.
        - Refactored out of 'compute_signal.py' for modularization.
    """
    #  Guard against negative precision or a zero bin
    if siz_bin <= 0:
        raise ValueError("'siz_bin' must be > 0.")
    if rnd < 0:
        raise ValueError("'rnd' must be >= 0.")

    #  (Updated approach: More memory-efficient stream-write)
    with open_out(fil_out) as fh:
        #  Iterate over the signal data, sorted by chromosome (using
        #  Roman-to-Arabic numeral conversion) and bin start position
        for (chrom, bin_start), value in sorted(
            cvg.items(),
            key=lambda kv: key_bin(kv[0][0], kv[0][1])
        ):
            v = round(value, rnd)
            if v == 0.0:
                v = 0.0  # Collapse “-0”

            #  Format with at most 'rnd' decimal places, then strip only
            #  non-informative trailing zeros and any trailing decimal point
            s_val = f"{v:.{rnd}f}"
            if "." in s_val:
                s_val = s_val.rstrip("0").rstrip(".")
            if s_val == "-0":
                s_val = "0"

            #  Write chrom, start, end, and signal to the bedGraph file
            fh.write(
                f"{chrom}\t{bin_start}\t{bin_start + siz_bin}\t{s_val}\n"
            )


def generate_track_name(fil_out: str) -> str:
    """
    Generate a track file name by inserting '.track' before the main extension.

    Args:
        fil_out : str
            Original output filename.

    Returns:
        nam_trk : str
            Track filename with '.track' inserted before the main extension.

    Raises:
        None.

    Notes:
        - Preserves the original main extension spelling (for example,
          '.bedGraph' remains '.track.bedGraph').
        - Preserves any trailing '.gz'.
        - Refactored out of 'compute_signal_ratio.py' for modularization.
    """
    #  Check if the file is gzipped
    is_gz = fil_out.endswith(".gz")

    #  Remove .gz temporarily
    bas = fil_out[:-3] if is_gz else fil_out

    #  Find the last extension (e.g., '.bdg', '.bedgraph')
    nam_bas, ext = os.path.splitext(bas)

    #  Construct the new track filename
    nam_trk = f"{nam_bas}.track{ext}"

    #  Re-add '.gz' if necessary
    if is_gz:
        nam_trk += ".gz"

    return nam_trk
