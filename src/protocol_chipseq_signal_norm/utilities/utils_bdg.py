#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: utils_bdg.py
#
# Copyright 2025-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5-series models; most recent: GPT-5.6) were
# used in development and documentation.
#
# Distributed under the MIT license.


"""
Utilities for bedGraph parsing, validation, sorting, and output.

Notes
-----
Helpers are importable and do not parse command-line arguments. BED and
bedGraph intervals are treated as 0-based, half-open coordinates.
"""

from __future__ import annotations

import os
import sys
from collections.abc import Callable, Iterator

from .utils_chrom import sort_chrom
from .utils_io import DEF_SKP_PFX, is_header, open_in, open_out, read_data_line

assert sys.version_info >= (3, 11), "Python >= 3.11 required."


def canon_nonfinite(tok: str) -> str | None:
    """
    Return a canonical non-finite token or None.

    Parameters
    ----------
        tok : str
            Raw string token.

    Returns
    -------
        'nan' | 'inf' | '-inf'
            If the token case-insensitively matches any NaN/±inf spelling
            (including '+inf'); otherwise None.

    Raises
    ------
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

    Parameters
    ----------
        tok : str
            Raw string token.

    Returns
    -------
        float | None
            On success, the parsed token as a float; on failure, None.

    Raises
    ------
        None.

    Notes
    -----
        Non-finite tokens like 'nan'/'inf' are not canonicalized here; combine
        with 'canon_nonfinite' to treat them specially.
    """
    try:
        return float(tok)
    except ValueError:
        return None


def iter_rows_bdg(
    stream, skp_prd: Callable[[str], bool]
) -> Iterator[tuple[str, int, int, str, float | None]]:
    """
    Stream bedGraph-like rows.

    Yields:
        '(chrom, start, end, val_tok, val_num)', where 'val_num' is None for
        non-finite tokens or unparsable values.

    Raises
    ------
        None.

    Notes
    -----
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


def check_size_bin(
    fil_A: str,
    fil_B: str,
    skp_pfx: tuple[str, ...] = DEF_SKP_PFX,
    max_records: int | None = 5
) -> None:
    """
    Verify that the bin sizes in two BED or bedGraph files are identical.

    Parameters
    ----------
        fil_A : str
            Path to first BED or bedGraph input file (file A; e.g., IP).
        fil_B : str
            Path to second BED or bedGraph input file (file B; e.g., input).
        skp_pfx : tuple[str, ...]
            Prefixes to skip as BED or bedGraph header/meta lines.
        max_records : int | None
            Maximum number of paired rows to compare. If None, compare until
            either file reaches EOF.

    Returns
    -------
        None.

    Raises
    ------
        ValueError
            If the bin sizes differ, if lines are malformed, or if positions
            are not parseable as integers.

    Notes
    -----
        Refactored out of 'compute_signal_ratio.py' for modularization.
    """
    if max_records is not None and max_records < 0:
        raise ValueError("'max_records' must be >= 0 or None.")

    #  Compare bins from each file to ensure they use the same bin size (if
    #  not, fast fail)
    with open_in(fil_A) as opn_A, open_in(fil_B) as opn_B:
        seen = 0
        while max_records is None or seen < max_records:
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


def _parse_coords_bdg(
    line: str,
    label: str
) -> tuple[str, int, int]:
    """
    Parse chromosome, start, and end from one BED or bedGraph data row.
    """
    fields = line.split()
    if len(fields) < 3:
        raise ValueError(
            f"Malformed BED or bedGraph line in {label}: {line!r}"
        )

    try:
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
    except ValueError as e:
        raise ValueError(
            f"Non-integer start/end encountered in {label}: {line!r}"
        ) from e

    if start < 0:
        raise ValueError(f"Negative start coordinate in {label}: {line!r}")
    if end <= start:
        raise ValueError(f"Non-positive interval width in {label}: {line!r}")

    return chrom, start, end


def check_grid_bin(
    fil_A: str,
    fil_B: str,
    skp_pfx: tuple[str, ...] = DEF_SKP_PFX
) -> None:
    """
    Verify that two BED or bedGraph files have the same ordered bin grid.

    Parameters
    ----------
        fil_A : str
            Path to first BED or bedGraph input file.
        fil_B : str
            Path to second BED or bedGraph input file.
        skp_pfx : tuple[str, ...]
            Prefixes to skip as BED or bedGraph header/meta lines.

    Returns
    -------
        None.

    Raises
    ------
        ValueError
            If either file has malformed coordinate rows, the ordered
            '(chrom, start, end)' grids differ, or one file has extra rows.
    """
    with open_in(fil_A) as opn_A, open_in(fil_B) as opn_B:
        idx = 0
        while True:
            lin_A = read_data_line(fh=opn_A, skp_pfx=skp_pfx)
            lin_B = read_data_line(fh=opn_B, skp_pfx=skp_pfx)

            if not lin_A and not lin_B:
                return
            if not lin_A or not lin_B:
                raise ValueError(
                    "Mismatched bedGraph grids: one file has additional data "
                    f"rows after row {idx}."
                )

            row_A = _parse_coords_bdg(lin_A, "first file (file A)")
            row_B = _parse_coords_bdg(lin_B, "second file (file B)")

            if row_A != row_B:
                raise ValueError(
                    "Mismatched bedGraph grids at paired data row "
                    f"{idx + 1}:\n"
                    f"  First file (file A):  {row_A}\n"
                    f"  Second file (file B): {row_B}"
                )

            idx += 1


def load_chr_sizes(path: str) -> dict[str, int]:
    """
    Parse a UCSC-style chromosome-size TSV.

    Parameters
    ----------
        path : str
            Path to a chromosome-size file with at least two whitespace-
            separated columns: chromosome name and positive integer size.

    Returns
    -------
        dict[str, int]
            Mapping from chromosome name to size in base pairs.

    Raises
    ------
        ValueError
            If the file is empty, malformed, contains duplicate chromosome
            names, or contains non-positive sizes.
        OSError
            If the file cannot be read.
    """
    chrom_sizes: dict[str, int] = {}

    with open_in(path) as fh:
        for idx, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split()
            if len(fields) < 2:
                raise ValueError(
                    f"Malformed chr.sizes line {idx}: expected at least two "
                    f"columns, got {raw.rstrip()!r}."
                )

            chrom, size_raw = fields[:2]
            if chrom in chrom_sizes:
                raise ValueError(
                    f"Duplicate chromosome in chr.sizes line {idx}: {chrom!r}."
                )

            try:
                size = int(size_raw)
            except ValueError as e:
                raise ValueError(
                    f"Non-integer chromosome size in chr.sizes line {idx}: "
                    f"{raw.rstrip()!r}."
                ) from e

            if size <= 0:
                raise ValueError(
                    f"Chromosome size must be > 0 in chr.sizes line {idx}: "
                    f"{raw.rstrip()!r}."
                )

            chrom_sizes[chrom] = size

    if not chrom_sizes:
        raise ValueError(f"No chromosome sizes found in {path!r}.")

    return chrom_sizes


def validate_bounds_bdg(
    path: str,
    chrom_sizes: dict[str, int],
    skp_pfx: tuple[str, ...] = DEF_SKP_PFX,
    label: str = "bedGraph"
) -> None:
    """
    Validate bedGraph coordinates against chromosome sizes.

    Parameters
    ----------
        path : str
            Path to a bedGraph file.
        chrom_sizes : dict[str, int]
            Chromosome-size mapping.
        skp_pfx : tuple[str, ...]
            Prefixes to skip as BED or bedGraph header/meta lines.
        label : str
            Human-readable file label used in error messages.

    Returns
    -------
        None.

    Raises
    ------
        ValueError
            If a data row is malformed, uses an unknown chromosome, has invalid
            coordinates, or extends beyond the declared chromosome size.
    """
    with open_in(path) as fh:
        for idx, raw in enumerate(fh, start=1):
            if is_header(raw, skp_pfx):
                continue

            chrom, start, end = _parse_coords_bdg(
                raw.rstrip("\n"), f"{label} line {idx}"
            )

            chrom_size = chrom_sizes.get(chrom)
            if chrom_size is None:
                raise ValueError(
                    f"Chromosome {chrom!r} in {label} line {idx} is absent "
                    "from the chr.sizes file."
                )
            if start >= chrom_size:
                raise ValueError(
                    f"Start coordinate in {label} line {idx} is outside "
                    f"chromosome {chrom!r} (size {chrom_size}): "
                    f"{raw.rstrip()!r}"
                )
            if end > chrom_size:
                raise ValueError(
                    f"End coordinate in {label} line {idx} extends beyond "
                    f"chromosome {chrom!r} (size {chrom_size}): "
                    f"{raw.rstrip()!r}"
                )


def key_bin(chrom: str, start: int) -> tuple[tuple[int, int, str], int]:
    """Return the sortable chromosome/start key for a bedGraph bin."""

    #  Use Roman-order first via sort_chrom(), then start coordinate
    #  (refactored out of 'compute_signal_ratio.py' for modularization)
    return (sort_chrom(chrom), start)


def write_bdg(cvg, fil_out, siz_bin, dp, chrom_sizes=None):
    """
    Write binned signal data to a bedGraph file.

    Parameters
    ----------
        cvg : dict
            Binned signal data, where keys are (chrom, bin_start) and float
            values.
        fil_out: str
            Path to the output file: '.bedGraph', '.bedgraph', '.bdg', or
            '.bg', optionally with '.gz'.
        siz_bin: int
            Bin size in base pairs.
        dp : int
            Maximum number of decimal places retained for signal values.
        chrom_sizes : dict[str, int] | None
            Optional chromosome sizes. If provided, output interval ends are
            clamped to chromosome ends and bins outside declared chromosome
            bounds are rejected.

    Returns
    -------
        None. Writes interval records to the output file.

    Raises
    ------
        ValueError
            If 'siz_bin' <= 0 or 'dp' < 0.
        OSError
            On filesystem/permission issues while writing.

    Notes
    -----
        - Signal values are rounded to at most 'dp' decimal places.
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
    if dp < 0:
        raise ValueError("'dp' must be >= 0.")

    #  (Updated approach: More memory-efficient stream-write)
    with open_out(fil_out) as fh:
        #  Iterate over the signal data, sorted by chromosome (using
        #  Roman-to-Arabic numeral conversion) and bin start position
        for (chrom, bin_start), value in sorted(
            cvg.items(),
            key=lambda kv: key_bin(kv[0][0], kv[0][1])
        ):
            v = round(value, dp)
            if v == 0.0:
                v = 0.0  # Collapse “-0”

            #  Format with at most 'dp' decimal places, then strip only
            #  non-informative trailing zeros and any trailing decimal point
            s_val = f"{v:.{dp}f}"
            if "." in s_val:
                s_val = s_val.rstrip("0").rstrip(".")
            if s_val == "-0":
                s_val = "0"

            bin_end = bin_start + siz_bin
            if chrom_sizes is not None:
                chrom_size = chrom_sizes.get(chrom)
                if chrom_size is None:
                    raise ValueError(
                        f"Missing chromosome size for bedGraph row: {chrom!r}."
                    )
                if bin_start < 0 or bin_start >= chrom_size:
                    raise ValueError(
                        "bedGraph bin start is outside chromosome bounds: "
                        f"{chrom}:{bin_start} (chromosome size {chrom_size})."
                    )
                bin_end = min(bin_end, chrom_size)

            #  Write chrom, start, end, and signal to the bedGraph file
            fh.write(f"{chrom}\t{bin_start}\t{bin_end}\t{s_val}\n")


def generate_name_track(fil_out: str) -> str:
    """
    Generate a track file name by inserting '.track' before the main extension.

    Parameters
    ----------
        fil_out : str
            Original output filename.

    Returns
    -------
        nam_trk : str
            Track filename with '.track' inserted before the main extension.

    Raises
    ------
        None.

    Notes
    -----
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
