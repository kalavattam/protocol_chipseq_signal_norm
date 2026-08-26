#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: utils_bdg.py
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
Utilities for bedGraph parsing, validation, sorting, and output.

Notes
-----
Helpers are importable and do not parse command-line arguments. BED and
bedGraph intervals are treated as 0-based, half-open coordinates.
"""

from __future__ import annotations

import math
import os
import sys
from collections.abc import Callable, Iterator
from dataclasses import dataclass
from typing import TextIO

from .utils_chrom import sort_chrom
from .utils_io import DEF_SKP_PFX, is_header, open_in, open_out, read_data_line

assert sys.version_info >= (3, 11), "Python >= 3.11 required."


def canon_nonfinite(token: str) -> str | None:
    """
    Return a canonical non-finite token or None.

    Parameters
    ----------
    token : str
        Raw string token.

    Returns
    -------
    token : str | None
        If the token case-insensitively matches any NaN/±inf spelling
        (including '+inf'); otherwise None.
    """

    normalized = token.strip().lower()

    if normalized in ("nan", "inf", "-inf", "+inf"):
        if normalized == "nan":
            return "nan"

        if normalized.startswith("-"):
            return "-inf"

        return "inf"

    return None


def try_float(token: str) -> float | None:
    """
    Parse a float token.

    Parameters
    ----------
    token : str
        Raw string token.

    Returns
    -------
    value : float | None
        On success, the parsed token as a float; on failure, None.

    Notes
    -----
    Non-finite tokens like 'nan'/'inf' are not canonicalized here; combine with
    'canon_nonfinite' to treat them specially.
    """

    try:
        return float(token)
    except ValueError:
        return None


def iter_rows_bdg(
    stream: TextIO,
    skip_predicate: Callable[[str], bool],
) -> Iterator[tuple[str, int, int, str, float | None]]:
    """
    Stream bedGraph-like rows.

    Parameters
    ----------
    stream : TextIO
        Text stream containing bedGraph records.
    skip_predicate : Callable[[str], bool]
        Predicate returning whether an input row should be skipped.

    Yields
    ------
    chromosome, start, end, value_text, numeric_value : tuple[
        str, int, int, str, float | None
    ]
        Chromosome, start, end, original value token, and parsed value. The
        parsed value is None for non-finite or invalid tokens.

    Notes
    -----
    - Skips lines where 'skip_predicate(line)' is True.
    - Skips lines with < 4 fields, non-integer coords, or non-positive width.
    - Coordinates are not clamped; caller may clamp to chrom bounds.
    """

    for raw_line in stream:
        if skip_predicate(raw_line):
            continue

        fields = raw_line.rstrip("\n").split()

        if len(fields) < 4:
            continue

        chromosome, start_text, end_text, value_text = fields[:4]

        try:
            start = int(start_text)
            end = int(end_text)
        except ValueError:
            continue

        if end <= start:
            continue

        nonfinite_name = canon_nonfinite(value_text)
        numeric_value = (
            None if nonfinite_name is not None else try_float(value_text)
        )

        yield chromosome, start, end, value_text, numeric_value


@dataclass(frozen=True)
class LibrarySizeBdg:
    """
    Represent a bedGraph library size and the bin width it rests on.

    The total is the column sum over fixed-width bins, which is edgeR's
    'lib.size' for the matrix the track represents. The bin width is carried
    alongside it because a caller that omitted the width still needs the one
    the track implied, and because the total means nothing without it.
    """

    total: float
    siz_bin: int


def sum_counts_bdg(
    path: str,
    siz_bin: int | None = None,
    skp_pfx: tuple[str, ...] | None = None,
) -> LibrarySizeBdg:
    """
    Sum a bedGraph over fixed-width bins to recover a library size.

    Parameters
    ----------
    path : str
        BedGraph-like path, or '-' for standard input.
    siz_bin : int | None
        Bin width in base pairs used to write the track. When None, infer it
        from the track. When supplied, cross-check it against the track and
        reject a value the track contradicts.
    skp_pfx : tuple[str, ...] | None
        Header prefixes to skip. The default is 'DEF_SKP_PFX'.

    Returns
    -------
    result : LibrarySizeBdg
        The library size and the bin width it was computed on.

    Raises
    ------
    ValueError
        If 'siz_bin' is supplied and is not a positive integer, if the track
        holds no usable interval, or if a supplied 'siz_bin' disagrees with the
        width the track implies.

    Notes
    -----
    A bedGraph run-length encodes equal neighbours, so an interval spanning
    'k' bins contributes 'k' times its value rather than once. Summing rows
    directly understates the library size by exactly the compression factor.

    Terminal intervals round up: a chromosome whose length is not a multiple of
    the bin width ends in a partial bin that still holds a real count. Only a
    chromosome's final interval can be partial, so the total folds in every
    earlier one and defers the terminal ones.

    The bin width is the greatest common divisor of the interval starts and of
    the non-terminal widths, both multiples of it. That divisor is always a
    multiple of the true width, and equals it unless every run shares a factor:
    a 10 bp track whose runs all span an even number of bins infers 20 and
    halves the library size, silently. Pass 'siz_bin' where the track cannot be
    trusted to resolve it.

    This is the column sum of the bin matrix, not the alignment count. Under
    the overlap counting in 'countReadsPerBin.py:705' one fragment increments
    every bin it touches, so the two differ by roughly '(F / siz_bin) + 1'.
    Passing the alignment count where this value is expected rescales every
    pseudocount by that factor.
    """

    if siz_bin is not None and (not isinstance(siz_bin, int) or siz_bin <= 0):
        raise ValueError("'siz_bin' must be a positive integer.")

    if skp_pfx is None:
        skp_pfx = DEF_SKP_PFX

    def skip_predicate(line: str) -> bool:
        """
        Return whether a raw input row is non-data content.
        """

        return is_header(line, skp_pfx)

    grain = 0
    total = 0.0
    order: list[str] = []
    final: dict[str, tuple[float, int]] = {}

    with open_in(path) as stream:
        for chrom, start, end, _text, value in iter_rows_bdg(
            stream,
            skip_predicate,
        ):
            if value is None or not math.isfinite(value):
                continue

            width = end - start

            if width <= 0:
                continue

            grain = math.gcd(grain, start)

            # Fold in the previous interval for this chromosome now that a
            # later one proves it was not terminal, so only genuinely final
            # intervals stay behind for the rounding-up step below.
            if chrom in final:
                value_prior, width_prior = final[chrom]
                total += value_prior * width_prior
                grain = math.gcd(grain, width_prior)
            else:
                order.append(chrom)

            final[chrom] = (value, width)

    if not order:
        raise ValueError(f"no usable bedGraph interval in '{path}'.")

    if grain <= 0:
        # A single-interval track starting at zero resolves nothing.
        grain = 0 if siz_bin is None else siz_bin

    if siz_bin is None:
        if grain <= 0:
            raise ValueError(
                f"cannot infer a bin width from '{path}'; pass '--siz_bin'.",
            )

        siz_bin = grain
    elif grain > 0 and grain % siz_bin:
        raise ValueError(
            f"'siz_bin' {siz_bin} disagrees with '{path}', whose intervals "
            f"imply a bin width of {grain}.",
        )

    total /= siz_bin

    for chrom in order:
        value, width = final[chrom]
        n_bin = width // siz_bin

        if width % siz_bin:
            n_bin += 1

        total += value * n_bin

    return LibrarySizeBdg(total=total, siz_bin=siz_bin)


def check_size_bin(
    fil_a: str,
    fil_b: str,
    skp_pfx: tuple[str, ...] = DEF_SKP_PFX,
    max_records: int | None = 5,
) -> None:
    """
    Verify that the bin sizes in two BED or bedGraph files are identical.

    Parameters
    ----------
    fil_a : str
        Path to first BED or bedGraph input file (file A; e.g., IP).
    fil_b : str
        Path to second BED or bedGraph input file (file B; e.g., input).
    skp_pfx : tuple[str, ...]
        Prefixes to skip as BED or bedGraph header/meta lines.
    max_records : int | None
        Maximum number of paired rows to compare. If None, compare until either
        file reaches EOF.

    Raises
    ------
    ValueError
        If the bin sizes differ, if lines are malformed, or if positions are
        not parseable as integers.

    Notes
    -----
    Refactored out of 'compute_signal_ratio.py' for modularization.
    """

    if max_records is not None and max_records < 0:
        raise ValueError("'max_records' must be >= 0 or None.")

    with open_in(fil_a) as opn_a, open_in(fil_b) as opn_b:
        seen = 0

        while max_records is None or seen < max_records:
            lin_a = read_data_line(handle=opn_a, skp_pfx=skp_pfx)
            lin_b = read_data_line(handle=opn_b, skp_pfx=skp_pfx)
            if not lin_a or not lin_b:
                break

            fields_a = lin_a.split()
            fields_b = lin_b.split()

            if len(fields_a) < 3:
                raise ValueError(
                    "Malformed BED or bedGraph line in first file (file A; "
                    f"e.g., IP) during bin-size check: {lin_a!r}",
                )

            if len(fields_b) < 3:
                raise ValueError(
                    "Malformed BED or bedGraph line in second file (file B; "
                    f"e.g., input) during bin-size check: {lin_b!r}",
                )

            try:
                bin_a = int(fields_a[2]) - int(fields_a[1])
                bin_b = int(fields_b[2]) - int(fields_b[1])
            except ValueError as error:
                raise ValueError(
                    "Non-integer start/end encountered during bin-size "
                    "check:\n"
                    f"  First file (file A; e.g., IP):     {lin_a!r}\n"
                    f"  Second file (file B; e.g., input): {lin_b!r}",
                ) from error

            if bin_a != bin_b:
                raise ValueError(
                    f"Error: Mismatched bin sizes detected.\n"
                    f"  - {fil_a}: {bin_a}-bp bins\n"
                    f"  - {fil_b}: {bin_b}-bp bins\n\n"
                    f"Ensure both files are binned identically.",
                )

            seen += 1


def _parse_coords_bdg(line: str, label: str) -> tuple[str, int, int]:
    """
    Parse chromosome, start, and end from one BED or bedGraph data row.

    Parameters
    ----------
    line : str
        Complete nonempty BED or bedGraph data row.
    label : str
        Input identity used in validation diagnostics.

    Returns
    -------
    chrom, start, end : tuple[str, int, int]
        Chromosome, nonnegative start, and strictly greater end coordinate.

    Raises
    ------
    ValueError
        If fields are missing, nonnumeric, negative, or non-increasing.
    """

    fields = line.split()
    if len(fields) < 3:
        raise ValueError(
            f"Malformed BED or bedGraph line in {label}: {line!r}",
        )

    try:
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
    except ValueError as e:
        raise ValueError(
            f"Non-integer start/end encountered in {label}: {line!r}",
        ) from e

    if start < 0:
        raise ValueError(f"Negative start coordinate in {label}: {line!r}")

    if end <= start:
        raise ValueError(f"Non-positive interval width in {label}: {line!r}")

    return chrom, start, end


def check_grid_bin(
    fil_a: str,
    fil_b: str,
    skp_pfx: tuple[str, ...] = DEF_SKP_PFX,
) -> None:
    """
    Verify that two BED or bedGraph files have the same ordered bin grid.

    Parameters
    ----------
    fil_a : str
        Path to first BED or bedGraph input file.
    fil_b : str
        Path to second BED or bedGraph input file.
    skp_pfx : tuple[str, ...]
        Prefixes to skip as BED or bedGraph header/meta lines.

    Raises
    ------
    ValueError
        If either file has malformed coordinate rows, the ordered
        '(chrom, start, end)' grids differ, or one file has extra rows.
    """

    with open_in(fil_a) as opn_a, open_in(fil_b) as opn_b:
        idx = 0

        while True:
            lin_a = read_data_line(handle=opn_a, skp_pfx=skp_pfx)
            lin_b = read_data_line(handle=opn_b, skp_pfx=skp_pfx)

            if not lin_a and not lin_b:
                return

            if not lin_a or not lin_b:
                raise ValueError(
                    "Mismatched bedGraph grids: one file has additional data "
                    f"rows after row {idx}.",
                )

            row_a = _parse_coords_bdg(lin_a, "first file (file A)")
            row_b = _parse_coords_bdg(lin_b, "second file (file B)")

            if row_a != row_b:
                raise ValueError(
                    "Mismatched bedGraph grids at paired data row "
                    f"{idx + 1}:\n"
                    f"  First file (file A):  {row_a}\n"
                    f"  Second file (file B): {row_b}",
                )

            idx += 1


def load_chromosome_sizes(path: str) -> dict[str, int]:
    """
    Parse a UCSC-style chromosome-size TSV.

    Parameters
    ----------
    path : str
        Path to a chromosome-size file with at least two whitespace-separated
        columns: chromosome name and positive integer size.

    Returns
    -------
    chromosome_sizes : dict[str, int]
        Mapping from chromosome name to size in base pairs.

    Raises
    ------
    ValueError
        If the file is empty, malformed, contains duplicate chromosome names,
        or contains non-positive sizes.
    OSError
        If the file cannot be read.
    """

    chromosome_sizes: dict[str, int] = {}

    with open_in(path) as handle:
        for row_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()

            if not line or line.startswith("#"):
                continue

            fields = line.split()

            if len(fields) < 2:
                raise ValueError(
                    f"Malformed chr.sizes line {row_number}: expected at "
                    f"least two columns, got {raw_line.rstrip()!r}.",
                )

            chrom, size_raw = fields[:2]

            if chrom in chromosome_sizes:
                raise ValueError(
                    f"Duplicate chromosome in chr.sizes line {row_number}: "
                    f"{chrom!r}.",
                )

            try:
                size = int(size_raw)
            except ValueError as error:
                raise ValueError(
                    f"Non-integer chromosome size in chr.sizes line "
                    f"{row_number}: {raw_line.rstrip()!r}.",
                ) from error

            if size <= 0:
                raise ValueError(
                    f"Chromosome size must be > 0 in chr.sizes line "
                    f"{row_number}: {raw_line.rstrip()!r}.",
                )

            chromosome_sizes[chrom] = size

    if not chromosome_sizes:
        raise ValueError(f"No chromosome sizes found in {path!r}.")

    return chromosome_sizes


def validate_bounds_bdg(
    path: str,
    chrom_sizes: dict[str, int],
    skp_pfx: tuple[str, ...] = DEF_SKP_PFX,
    label: str = "bedGraph",
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

    Raises
    ------
    ValueError
        If a data row is malformed, uses an unknown chromosome, has invalid
        coordinates, or extends beyond the declared chromosome size.
    """

    with open_in(path) as handle:
        for idx, raw in enumerate(handle, start=1):
            if is_header(raw, skp_pfx):
                continue

            chrom, start, end = _parse_coords_bdg(
                raw.rstrip("\n"),
                f"{label} line {idx}",
            )

            chrom_size = chrom_sizes.get(chrom)
            if chrom_size is None:
                raise ValueError(
                    f"Chromosome {chrom!r} in {label} line {idx} is absent "
                    "from the chr.sizes file.",
                )

            if start >= chrom_size:
                raise ValueError(
                    f"Start coordinate in {label} line {idx} is outside "
                    f"chromosome {chrom!r} (size {chrom_size}): "
                    f"{raw.rstrip()!r}",
                )

            if end > chrom_size:
                raise ValueError(
                    f"End coordinate in {label} line {idx} extends beyond "
                    f"chromosome {chrom!r} (size {chrom_size}): "
                    f"{raw.rstrip()!r}",
                )


def key_bin(chrom: str, start: int) -> tuple[tuple[int, int, str], int]:
    """
    Return the sortable chromosome/start key for a bedGraph bin.

    Parameters
    ----------
    chrom : str
        Chromosome identifier for the selected signal calculation.
    start : int
        Zero-based start coordinate for the selected bin.

    Returns
    -------
    key : tuple[tuple[int, int, str], int]
        The chromosome, start, and end key for one bedGraph row.
    """

    return (sort_chrom(chrom), start)


def write_bdg(
    coverage: dict[tuple[str, int], float],
    fil_out: str,
    siz_bin: int,
    decimal_places: int,
    chrom_sizes: dict[str, int] | None = None,
) -> None:
    """
    Write binned signal data to a bedGraph file.

    Parameters
    ----------
    coverage : dict[tuple[str, int], float]
        Binned signal data, where keys are (chrom, bin_start) and float values.
    fil_out : str
        Path to the output file: '.bedGraph', '.bedgraph', '.bdg', or '.bg',
        optionally with '.gz'.
    siz_bin : int
        Bin size in base pairs.
    decimal_places : int
        Maximum number of decimal places retained for signal values.
    chrom_sizes : dict[str, int] | None
        Optional chromosome sizes. If provided, output interval ends are
        clamped to chromosome ends and bins outside declared chromosome bounds
        are rejected.

    Raises
    ------
    ValueError
        If 'siz_bin' <= 0 or 'decimal_places' < 0.
    OSError
        On filesystem/permission issues while writing.

    Notes
    -----
    The function writes interval records to 'fil_out'.

    - Signal values are rounded to at most 'decimal_places' decimal places.
    - After rounding, non-informative trailing zeros are stripped from finite
      decimal representations, and any trailing decimal point is removed.
    - Negative zero is emitted as '0'.
    - Output is automatically gzip-compressed if 'fil_out' ends with '.gz'.
    - bedGraph format: 'chrom<tab>start<tab>end<tab>signal<newline>'.
    - Refactored out of 'compute_signal.py' for modularization.
    """

    if siz_bin <= 0:
        raise ValueError("'siz_bin' must be > 0.")

    if decimal_places < 0:
        raise ValueError("'dp' must be >= 0.")

    with open_out(fil_out) as handle:
        for (chrom, bin_start), value in sorted(
            coverage.items(),
            key=lambda kv: key_bin(kv[0][0], kv[0][1]),
        ):
            rounded_value = round(value, decimal_places)

            if rounded_value == 0.0:
                rounded_value = 0.0

            value_text = f"{rounded_value:.{decimal_places}f}"

            if "." in value_text:
                value_text = value_text.rstrip("0").rstrip(".")

            if value_text == "-0":
                value_text = "0"

            bin_end = bin_start + siz_bin

            if chrom_sizes is not None:
                chrom_size = chrom_sizes.get(chrom)
                if chrom_size is None:
                    raise ValueError(
                        "Missing chromosome size for bedGraph row: "
                        f"{chrom!r}.",
                    )

                if bin_start < 0 or bin_start >= chrom_size:
                    raise ValueError(
                        "bedGraph bin start is outside chromosome bounds: "
                        f"{chrom}:{bin_start} (chromosome size {chrom_size}).",
                    )

                bin_end = min(bin_end, chrom_size)

            handle.write(
                f"{chrom}\t{bin_start}\t{bin_end}\t{value_text}\n",
            )


def generate_name_track(fil_out: str) -> str:
    """
    Generate a track file name by inserting '.track' before the main extension.

    Parameters
    ----------
    fil_out : str
        Original output filename.

    Returns
    -------
    track_path : str
        Track filename with '.track' inserted before the main extension.

    Notes
    -----
    - Preserves the original main extension spelling (for example, '.bedGraph'
      remains '.track.bedGraph').
    - Preserves any trailing '.gz'.
    - Refactored out of 'compute_signal_ratio.py' for modularization.
    """

    is_gzip = fil_out.endswith(".gz")
    base_path = fil_out[:-3] if is_gzip else fil_out
    base_name, extension = os.path.splitext(base_path)
    track_path = f"{base_name}.track{extension}"

    if is_gzip:
        track_path += ".gz"

    return track_path
