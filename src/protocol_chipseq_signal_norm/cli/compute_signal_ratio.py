#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: compute_signal_ratio.py
#
# Copyright 2025-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
# GPT-5.6) were used in development and documentation.
#
# Distributed under the MIT license.


"""
Compute per-bin ratios between two bedGraph tracks.

The CLI accepts input bedGraphs, output path, ratio method, scaling,
pseudocount, denominator floor, filtering, and formatting options. It writes a
bedGraph-like ratio track and, when requested, a finite-valued '.track' sidecar
suitable for genome browsers.

Examples
--------
python -m protocol_chipseq_signal_norm.cli.compute_signal_ratio [options] \\
    --fil_A <file> --fil_B <file> --fil_out <file>
"""

from __future__ import annotations

import argparse
import math
import signal
import sys
from contextlib import nullcontext, redirect_stdout, suppress
from typing import TextIO

from protocol_chipseq_signal_norm.utilities.utils_bdg import (
    check_grid_bin,
    check_size_bin,
    generate_name_track,
    key_bin,
    load_chromosome_sizes,
    validate_bounds_bdg,
)
from protocol_chipseq_signal_norm.utilities.utils_check import (
    check_exists,
    check_writable,
    validate_comparison,
    validate_output_path,
)
from protocol_chipseq_signal_norm.utilities.utils_cli import (
    CapArgumentParser,
    add_help_cap,
)
from protocol_chipseq_signal_norm.utilities.utils_io import (
    DEF_SKP_PFX,
    open_in,
    open_out,
    parse_skp_pfx,
    read_data_line,
)

with suppress(AttributeError, ValueError):
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

assert sys.version_info >= (3, 11), "Python >= 3.11 required."

# TODO: Continue strengthening bin-grid validation beyond the quick upfront
# checks, while keeping default runtime practical for large bedGraph pairs.

# Restrict output to bedGraph filename extensions.
BEDGRAPH_FORMATS = ("bedGraph", "bedgraph", "bdg", "bg")

# Map accepted `--method` values to canonical internal names.
METHOD_CANON = {
    # Compute the simple unadjusted ratio as A / B.
    "r": "unadj",
    "raw": "unadj",
    "u": "unadj",
    "unadj": "unadj",
    "unadjusted": "unadj",
    "s": "unadj",
    "smp": "unadj",
    "simple": "unadj",
    # Compute the log2 ratio as log2(A / B).
    "2": "log2",
    "l2": "log2",
    "lg2": "log2",
    "log2": "log2",
    # Compute the reciprocal simple ratio as B / A.
    "rr": "unadj_r",
    "raw_r": "unadj_r",
    "ur": "unadj_r",
    "unadj_r": "unadj_r",
    "unadjusted_r": "unadj_r",
    "sr": "unadj_r",
    "smp_r": "unadj_r",
    "simple_r": "unadj_r",
    # Compute the reciprocal log2 ratio as log2(B / A).
    "2r": "log2_r",
    "l2r": "log2_r",
    "l2_r": "log2_r",
    "lg2_r": "log2_r",
    "log2_r": "log2_r",
}
METHOD_CHOICES = tuple(METHOD_CANON.keys())


def parse_pair(val: str, def_sec: float) -> tuple[float, float]:
    """
    Parse 'A' or 'A:B' into two floats.

    A missing second value uses the caller-supplied default.

    Parameters
    ----------
    val : str
        String, e.g., '2.0', '2.0:1.5', etc.
    def_sec : float
        Default for the second value when only 'A' is provided.

    Returns
    -------
    first, second : tuple[float, float]
        Tuple (A, B).

    Raises
    ------
    argparse.ArgumentTypeError
        If the string is not of the form 'A' or 'A:B' with numeric parts.
    """

    parts = val.split(":")

    if len(parts) == 1:
        try:
            return float(parts[0]), def_sec
        except ValueError as e:
            raise argparse.ArgumentTypeError("Expected number for 'A'.") from e

    if len(parts) == 2:
        try:
            return float(parts[0]), float(parts[1])
        except ValueError as e:
            raise argparse.ArgumentTypeError("Expected numeric 'A:B'.") from e

    raise argparse.ArgumentTypeError("Expected 'A' or 'A:B'.")


def calc_rat_bin(
    sig_a: float,
    sig_b: float,
    scl_a: float | None,
    scl_b: float | None,
    psc_a: float | None,
    psc_b: float | None,
    dep_min: float | None,
    log2: bool,
    recip: bool,
    skip_00: str | None,
    eps: float = 0.0,
) -> float | None:
    """
    Compute one optionally transformed 'A / B' ratio.

    Per-file scaling, pseudocounts, a denominator floor, logarithmic
    transformation, and reciprocal output are applied as requested.

    Parameters
    ----------
    sig_a : float
        First file (e.g., IP) signal for a bin (A).
    sig_b : float
        Second file (e.g., input) signal for a bin (B).
    scl_a : float | None
        Per-file multiplicative scale factor for A. If None or 1.0, treated
        as neutral.
    scl_b : float | None
        Per-file multiplicative scale factor for B. If None or 1.0, treated
        as neutral.
    psc_a : float | None
        Pseudocount added to A (post-scaling). If None or 0.0, skipped.
    psc_b : float | None
        Pseudocount added to B (post-scaling). If None or 0.0, skipped.
    dep_min : float | None
        Minimum allowed denominator after any optional scaling and/or
        pseudocount addition. If provided, B := max(B, dep_min) to avoid
        extreme and undefined (e.g., 'n / 0') divisions.
    log2 : bool
        If True, return 'log2(A / B)'. If False, return linear 'A / B'.
    recip : bool
        If True, return the reciprocal of the computed ratio:
            - linear: '1 / (A / B)'
            - log2: '-log2(A / B)' (since 'log2(1 / x) = -log2(x)')
    skip_00 : str | None
        Optional zero-bin ('0 / 0' or 'ε / ε') drop stage. One of
        "pre_scale", "post_scale", or None.
            - "pre_scale": Test on raw values (before optional scaling
                           and/or pseudocounts addition).
            - "post_scale": Test after optional scaling, before optional
                            pseudocount addition.
            - None: Do not drop '0 / 0' bins.
    eps : float, default 0.0
        Tolerance (epsilon value, ε) for treating values as zero in the
        pre-pseudocount zero-zero check. Use 0.0 for exact-zero behavior,
        similar to what is done in deepTools. A tiny value (e.g., say,
        '1e-12') can be used to ignore "float noise."

    Returns
    -------
    value : float | None
        The computed value (ratio or log2 ratio, possibly reciprocated).

        May also return the following:
            - None  When 'skip_00' is specified and the bin is '0 / 0' or
                    'ε / ε'.
            - -inf  When 'log2' is requested and 'ratio == 0'.
            -  inf  When 'reciprocal' is requested on a zero linear ratio.
            -  nan  When the computation is undefined (e.g., negative ratio
                    for 'log2', or 'B == 0' even after clamping).

        Caller skips writing the bin when a return value is None.

    Notes
    -----
    This function is defensive: domain issues are mapped to sentinel
    return values instead of exceptions. Specifically,
            - '0 / 0' (or 'ε / ε') bins may return None (i.e., when
              '--skip_00' applies).
            - denominator underflow after optional clamping yields
              'float('nan')'.
            - with 'log2=True', 'ratio == 0' yields '-inf', and 'ratio < 0'
              yields 'nan'.
            - with 'log2=False' and 'recip=True', 'ratio == 0' yields
              'inf'.

        Callers should validate option domains (e.g., 'scl_fct > 0',
        'eps >= 0') prior to calling. A TypeError may still propagate if
        non-numeric inputs are passed.

    Order of operations:
        1. Optionally skip zero-zero bins: '0 / 0' or 'ε / ε' (deepTools-
           like). ‡
        2. Optionally scale each file.
        3. Optionally skip scaled zero-zero bins:
           '[(sf_A × 0) / (sf_B × 0)]' or
           '[(sf_A × ε) / (sf_B × ε)]'. ‡
        4. Optionally add pseudocounts.
        5. Optionally clamp denominator by 'dep_min'.
        6. If the denominator is <= ε, treat the bin as undefined.
        7. Divide 'A / B'.
        8. Optionally perform log2 transformation.
        9. Optionally compute reciprocal of #7 (linear) or #8 (log2).

        ‡ Either of optional #1 or optional #3, not both.
    """

    # 1. Optionally skip an empty pair before scaling.
    if skip_00 == "pre_scale" and abs(sig_a) <= eps and abs(sig_b) <= eps:
        return None

    # 2. Apply non-neutral scale factors.
    num = sig_a if (scl_a is None or scl_a == 1.0) else (scl_a * sig_a)
    den = sig_b if (scl_b is None or scl_b == 1.0) else (scl_b * sig_b)

    # 3. Optionally skip an empty scaled pair before adding pseudocounts.
    if skip_00 == "post_scale" and abs(num) <= eps and abs(den) <= eps:
        return None

    # 4. Add non-neutral pseudocounts.
    if psc_a not in (None, 0.0):
        num += psc_a

    if psc_b not in (None, 0.0):
        den += psc_b

    # 5. Clamp the denominator when a minimum was specified.
    if dep_min is not None and den < dep_min:
        den = dep_min

    # 6. Treat an effectively zero denominator as undefined.
    if abs(den) <= eps:
        return float("nan")

    # 7. Compute the ratio.
    ratio = num / den

    # 8. Apply the optional log2 transformation.
    if log2:
        if ratio > 0.0:
            transformed = math.log2(ratio)
        elif ratio == 0.0:
            transformed = float("-inf")
        else:
            transformed = float("nan")

        # 9. Apply the optional reciprocal in log space.
        return -transformed if recip else transformed

    # 9. Apply the optional reciprocal in linear space.
    if recip:
        return (1.0 / ratio) if ratio != 0.0 else float("inf")

    return ratio


def _format_ratio_value(value: float, decimal_places: int) -> str:
    """
    Format one finite or non-finite bedGraph ratio value.

    Parameters
    ----------
    value : float
        Ratio value to render.
    decimal_places : int
        Maximum finite decimal places.

    Returns
    -------
    text : str
        Canonical bedGraph value spelling.
    """

    if not math.isfinite(value):
        if math.isnan(value):
            return "nan"

        return "-inf" if value < 0 else "inf"

    rounded = round(value, decimal_places)

    if rounded == 0.0:
        rounded = 0.0

    rendered = f"{rounded:.{decimal_places}f}"

    if "." in rendered:
        rendered = rendered.rstrip("0").rstrip(".")

    return "0" if rendered == "-0" else rendered


def _write_ratio_row(
    main_output: TextIO,
    track_output: TextIO | None,
    chrom: str,
    start: int,
    end: int,
    value: float | None,
    *,
    decimal_places: int,
    drop_nonfinite: bool,
) -> None:
    """
    Write one main ratio row and its optional finite track row.

    Parameters
    ----------
    main_output : TextIO
        Writable main bedGraph handle.
    track_output : TextIO | None
        Optional writable finite-only track handle.
    chrom : str
        Chromosome or reference-sequence name.
    start : int
        Zero-based interval start.
    end : int
        Half-open interval end.
    value : float | None
        Computed ratio, or None for an omitted zero-zero bin.
    decimal_places : int
        Maximum finite decimal places.
    drop_nonfinite : bool
        Whether the main output omits non-finite values.
    """

    if value is None:
        return

    if drop_nonfinite and not math.isfinite(value):
        return

    rendered = _format_ratio_value(value, decimal_places)
    main_output.write(f"{chrom}\t{start}\t{end}\t{rendered}\n")

    if track_output is not None and math.isfinite(value):
        track_output.write(f"{chrom}\t{start}\t{end}\t{rendered}\n")


def _parse_ratio_row(
    line: str | None,
    *,
    ordinal: str,
    example: str,
) -> tuple[str, int, int, float] | None:
    """
    Parse and validate one ratio-input bedGraph row.

    Parameters
    ----------
    line : str | None
        Current data line, or None at end of stream.
    ordinal : str
        Public file label such as "first" or "second".
    example : str
        Public example label such as "IP" or "input".

    Returns
    -------
    row : tuple[str, int, int, float] | None
        Parsed chromosome, start, end, and value, or None at end of stream.

    Raises
    ------
    ValueError
        If the row is incomplete, nonnumeric, or non-positive width.
    """

    if not line:
        return None

    fields = line.split()
    if len(fields) < 4:
        raise ValueError(
            f"Malformed bedGraph line in {ordinal} file (e.g., {example}): "
            f"{line!r}",
        )

    try:
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        value = float(fields[3])
    except ValueError as exc:
        raise ValueError(
            (
                f"Non-numeric start/end/signal in {ordinal}-file "
                f"(e.g., {example}) line: {line!r}"
            ),
        ) from exc

    if start >= end:
        raise ValueError(
            f"Non-positive width in {ordinal}-file (e.g., {example}) line: "
            f"{line!r}",
        )

    return chrom, start, end, value


def comp_sig_rat(
    fil_a: str,
    fil_b: str,
    fil_out: str,
    scl_a: float,
    scl_b: float,
    psc_a: float,
    psc_b: float,
    dep_min: float | None,
    decimal_places: int,
    log2: bool,
    recip: bool,
    skip_00: str | None,
    eps: float,
    track: bool,
    drp_nan: bool,
    skp_pfx: tuple[str, ...] = DEF_SKP_PFX,
    strict_bins: bool = False,
    chrom_sizes: dict[str, int] | None = None,
) -> None:
    """
    Stream ratios between two bedGraph tracks.

    Scaling, pseudocounts, denominator clamping, logarithmic transformation,
    and reciprocal output are applied as requested. The merge matches bins by
    chromosome and start coordinate.

    See module docstring for more details.

    Parameters
    ----------
    fil_a : str
        Path to first (e.g., IP) bedGraph file (can be gzipped; A).
    fil_b : str
        Path to second (e.g., input) bedGraph file (can be gzipped; B).
    fil_out : str
        Output file path (if '.gz' extension, gzip compression).
    scl_a : float
        Scale factor for A.
    scl_b : float
        Scale factor for B.
    psc_a : float
        Pseudocount added to A after scaling (0.0 means no addition).
    psc_b : float
        Pseudocount added to B after scaling (0.0 means no addition).
    dep_min : float | None
        Denominator clamp ("minimum input depth") to avoid extreme and/or
        erroneous division.
    decimal_places : int
        Maximum number of decimal places retained for finite emitted
        values.
    log2 : bool
        If True, return 'log2(A / B)'. If False, return linear 'A / B'.
    recip : bool
        If True, return reciprocal of the final value:
            - linear path: '1 / (A / B)'
            - log2 path: '-log2(A / B)'
    skip_00 : str | None
        Optional zero-bin ('0 / 0' or 'ε / ε') drop stage. One of
        "pre_scale", "post_scale", or None.
            - "pre_scale": Test on raw values (before optional scaling
                           and/or pseudocounts addition).
            - "post_scale": Test after optional scaling, before optional
                            pseudocount addition.
            - None: Do not drop '0 / 0' bins.
    eps : float
        Epsilon for zero tests used in the pre-optional-pseudocount '0 / 0'
        drop and for post-clamp denominator zero checks.
    track : bool
        If True, write a '.track' sidecar omitting non-finite values
        ('inf', '-inf', and 'nan').
    drp_nan : bool
        If True, omit non-finite values from the main output as well.
    skp_pfx : tuple[str, ...]
        Prefixes to skip as bedGraph header/meta lines.
    strict_bins : bool
        Whether to require identical chromosome/start/end grids.
    chrom_sizes : dict[str, int] | None
        Optional chromosome sizes used to validate input bounds.

    Raises
    ------
    FileNotFoundError
        If an input file or the output directory is missing.
    PermissionError
        If the output directory is not writable.
    ValueError
        On malformed bedGraph lines, inconsistent bin sizes, invalid
        numeric arguments, or non-finite formatting issues.
    OSError
        On I/O errors opening, reading, and/or writing files.

    Notes
    -----
    The standard bedGraph contains all computed ratios unless '--drp_nan'
    omits rows with non-finite values. Finite values are rounded to at most
    'decimal_places' places and stripped of non-informative trailing zeros.
    When 'track' is true, a '.track' sidecar excludes all non-finite values.
    """

    if fil_a == "-" or fil_b == "-" or fil_out == "-":
        raise ValueError(
            "Dash input/output is no longer supported; provide real file "
            "paths for '--fil_A', '--fil_B', and '--fil_out'.",
        )

    if chrom_sizes is not None:
        validate_bounds_bdg(
            fil_a,
            chrom_sizes,
            skp_pfx=skp_pfx,
            label="first file (file A)",
        )
        validate_bounds_bdg(
            fil_b,
            chrom_sizes,
            skp_pfx=skp_pfx,
            label="second file (file B)",
        )

    if strict_bins:
        check_grid_bin(fil_a=fil_a, fil_b=fil_b, skp_pfx=skp_pfx)
    else:
        check_size_bin(fil_a=fil_a, fil_b=fil_b, skp_pfx=skp_pfx)

    track_path = generate_name_track(fil_out) if track else None

    def calculate(sig_a: float, sig_b: float) -> float | None:
        """
        Calculate one row ratio with the configured transformations.
        """

        return calc_rat_bin(
            sig_a=sig_a,
            sig_b=sig_b,
            scl_a=scl_a,
            scl_b=scl_b,
            psc_a=psc_a,
            psc_b=psc_b,
            dep_min=dep_min,
            log2=log2,
            recip=recip,
            skip_00=skip_00,
            eps=eps,
        )

    with (
        open_in(fil_a) as opn_a,
        open_in(fil_b) as opn_b,
        open_out(fil_out) as f_out,
        open_out(track_path) if track else nullcontext() as track_handle,
    ):
        lin_a = read_data_line(handle=opn_a, skp_pfx=skp_pfx)
        lin_b = read_data_line(handle=opn_b, skp_pfx=skp_pfx)

        while lin_a or lin_b:
            row_a = _parse_ratio_row(
                lin_a,
                ordinal="first",
                example="IP",
            )
            row_b = _parse_ratio_row(
                lin_b,
                ordinal="second",
                example="input",
            )
            advance_a = False
            advance_b = False

            if row_a is not None and row_b is not None:
                chrom_a, start_a, end_a, signal_a = row_a
                chrom_b, start_b, end_b, signal_b = row_b
                key_a = key_bin(chrom_a, start_a)
                key_b = key_bin(chrom_b, start_b)

                if key_a == key_b:
                    output = (
                        chrom_a,
                        start_a,
                        end_a,
                        calculate(signal_a, signal_b),
                    )
                    advance_a = True
                    advance_b = True
                elif key_a < key_b:
                    output = (
                        chrom_a,
                        start_a,
                        end_a,
                        calculate(signal_a, 0.0),
                    )
                    advance_a = True
                else:
                    output = (
                        chrom_b,
                        start_b,
                        end_b,
                        calculate(0.0, signal_b),
                    )
                    advance_b = True
            elif row_a is not None:
                chrom_a, start_a, end_a, signal_a = row_a
                output = (
                    chrom_a,
                    start_a,
                    end_a,
                    calculate(signal_a, 0.0),
                )
                advance_a = True
            else:
                assert row_b is not None

                chrom_b, start_b, end_b, signal_b = row_b

                output = (
                    chrom_b,
                    start_b,
                    end_b,
                    calculate(0.0, signal_b),
                )
                advance_b = True

            _write_ratio_row(
                f_out,
                track_handle,
                *output,
                decimal_places=decimal_places,
                drop_nonfinite=drp_nan,
            )

            if advance_a:
                lin_a = read_data_line(handle=opn_a, skp_pfx=skp_pfx)

            if advance_b:
                lin_b = read_data_line(handle=opn_b, skp_pfx=skp_pfx)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse command-line arguments.

    Parameters
    ----------
    argv : list[str] | None
        Explicit arguments, or None to read the process arguments.

    Returns
    -------
    arguments : argparse.Namespace
        Parsed numerator, denominator, scaling, and output options.

    Raises
    ------
    SystemExit
        If argument parsing fails or help is requested.
    """

    parser = CapArgumentParser(
        description=(
            "Compute per-bin ratios between two bedGraphs [A (numerator) and "
            "B (denominator)] with optional per-file multiplicative scaling, "
            "pseudocount addition, 'dep_min' “clamping” (denominator "
            "thresholding, i.e., a divisor “floor”), log2 transformation, and "
            "reciprocal computation.\n"
            "\n"
            "Can assign an error tolerance (ε [an “epsilon value”]) for "
            "treating values (e.g., “float noise”) as zero. deepTools "
            "bamCompare-like behavior is 'ε = 0.0' (the default).\n"
            "\n"
            "Order of operations is generally deepTools bamCompare-like:\n"
            "    1. Optionally skip zero-zero bins: '0 / 0' or 'ε / ε' "
            "(deepTools-like). ‡\n"
            "    2. Optionally scale each file.\n"
            "    3. Optionally skip scaled zero-zero bins: '[(sf_A × 0) / "
            "(sf_B × 0)]' or '[(sf_A × ε) / (sf_B × ε)]'. ‡\n"
            "    4. Optionally add pseudocounts.\n"
            "    5. Optionally clamp divisor (denominator) by 'dep_min' (see "
            "below).\n"
            "    6. If the divisor (denominator) is <= ε, treat the bin as "
            "undefined.\n"
            "    7. Divide 'A / B'.\n"
            "    8. Optionally perform log2 transformation.\n"
            "    9. Optionally compute reciprocal of #7 (linear) or #8 "
            "(log2).\n"
            "\n"
            "‡ Either of optional #1 or optional #3 can be applied, not "
            "both.\n"
            "\n"
            "(See module docstring for more details.)"
        ),
    )
    add_help_cap(parser)
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Run script in verbose mode.\n\n",
    )

    parser.add_argument(
        "-fA",
        "--fil_A",
        dest="fil_A",
        required=True,
        type=str,
        help=(
            "First bedGraph input file, file A (e.g., IP). Supports plain "
            "text and '.gz'.\n\n"
        ),
    )
    parser.add_argument(
        "-fB",
        "--fil_B",
        dest="fil_B",
        required=True,
        type=str,
        help=(
            "Second bedGraph input file, file B (e.g., input). Supports plain "
            "text and '.gz'.\n\n"
        ),
    )

    parser.add_argument(
        "-fo",
        "--fil_out",
        dest="fil_out",
        required=True,
        type=str,
        help=(
            "Output file path. Path to the output bedGraph file. When a real "
            "path is provided, accepted extensions are '.bedGraph', "
            "'.bedgraph', '.bdg', and '.bg', each optionally followed by "
            "'.gz'.\n\n"
        ),
    )
    parser.add_argument(
        "-cs",
        "--chr_sizes",
        dest="chr_sizes",
        default=None,
        help=(
            "Chromosome sizes file in UCSC-style TSV format used to validate "
            "input bedGraph interval bounds. This validation is independent "
            "of '--strict_bins'.\n\n"
        ),
    )

    parser.add_argument(
        "-me",
        "--method",
        dest="method",
        choices=METHOD_CHOICES,
        default="unadj",
        help=(
            "Workflow method. Ratio-computation subtype (default: "
            "'%(default)s').\n"
            "  - Unadjusted aliases: 'r', 'raw', 'u', 'unadj', 'unadjusted', "
            "'s', 'smp', 'simple'. Internally standardized to 'unadj'.\n"
            "  - Log2 aliases: '2', 'l2', 'lg2', 'log2'. Internally "
            "standardized to 'log2'.\n"
            "  - Reciprocal-unadjusted aliases: 'rr', 'raw_r', 'ur', "
            "'unadj_r', 'unadjusted_r', 'sr', 'smp_r', 'simple_r'. Internally "
            "standardized to 'unadj_r'.\n"
            "  - Reciprocal-log2 aliases: '2r', 'l2r', 'l2_r', 'lg2_r', "
            "'log2_r'. Internally standardized to 'log2_r'.\n\n"
        ),
    )
    parser.add_argument(
        "-sf",
        "--scl_fct",
        "--scl-fct",
        dest="scl_fct",
        type=str,
        default=None,
        help=(
            "Scaling factor. Per-file scale factor(s) 'A[:B]'. If only 'A' is "
            "given, 'B' defaults to 1.0.\n\n"
        ),
    )
    parser.add_argument(
        "-ps",
        "--pseudo",
        dest="pseudo",
        type=str,
        default="0:0",
        help=(
            "Per-file pseudocount spec 'A[:B]' added after scaling (default: "
            "%(default)s).\n"
            "\n"
            "Note: This is primarily useful for log2-ratio methods, where it "
            "helps avoid undefined values such as 'log2(A / 0)'; for linear "
            "ratios, '--dep_min' may be preferable when the main concern is "
            "bounding low-depth extremes.\n"
            "\n"
            "Using '--pseudo' together with '--dep_min' is allowed, but "
            "usually makes low-depth stabilization harder to interpret.\n\n"
        ),
    )
    parser.add_argument(
        "-dm",
        "--dep_min",
        "--dep-min",
        dest="dep_min",
        type=float,
        default=None,
        help=(
            "Minimum allowed denominator threshold or “clamp” (i.e., “minimum "
            "input depth” as described in PMID 40364978). This is a “floor” "
            "ensuring denominators do not fall below this value.\n"
            "\n"
            "The 'dep_min' value is applied after any optional scaling and/or "
            "pseudocount addition ('B := max(B, dep_min)').\n"
            "\n"
            "With or without clamping, a zero-guard with ε ('--eps') is "
            "applied after this optional step in the order of operations; if "
            "'--dep_min <flt>' is specified, choose 'dep_min > ε'.\n"
            "\n"
            "Using '--dep_min' together with '--pseudo' is allowed, but "
            "usually makes low-depth stabilization harder to interpret.\n\n"
        ),
    )

    parser.add_argument(
        "-e",
        "--eps",
        dest="eps",
        type=float,
        default=0.0,
        help=(
            "Zero tolerance epsilon for zero checks: Values satisfying "
            "'|value| <= ε' are treated as zero (default: %(default)s).\n"
            "\n"
            "Applied in two places:\n"
            "  - '--skip_00':\n"
            "    + If '--skip_00 pre_scale', check raw values (before any "
            "optional scaling and/or pseudocount addition).\n"
            "    + If '--skip_00 post_scale', check scaled values (before any "
            "optional pseudocount addition).\n"
            "  - Denominator guard:\n"
            "    + If '|B| <= ε', then the bin is emitted as 'nan'.\n"
            "    + If optional denominator clamping is specified, the "
            "denominator guard takes place after that.\n"
            "\n"
            "Notes:\n"
            "  - With 'ε = 0.0', '--skip_00 pre_scale' and '--skip_00 "
            "post_scale' are equivalent.\n"
            "  - With 'ε > 0' and 'scaling != 1', '--skip_00 post_scale' "
            "interprets ε in scaled units (i.e., the same units that are "
            "divided).\n\n"
        ),
    )
    parser.add_argument(
        "-s0",
        "--skip_00",
        "--skip-00",
        dest="skip_00",
        choices=["pre_scale", "post_scale"],
        default=None,
        help=(
            "Skip rows where both compared values are zero. Values satisfying "
            "'|value| <= ε' are treated as zero for this check.\n"
            "  - '--skip_00 pre_scale':  Test on raw values (deepTools "
            "bamCompare-like behavior).\n"
            "  - '--skip_00 post_scale': Test after scaling.\n"
            "  - Omit entirely to disable the '0 / 0' (or 'ε / ε') drop.\n"
            "\n"
            "Notes:\n"
            "  - If not omitted, then 'pre_scale' and 'post_scale' modes are "
            "equivalent if 'ε = 0' or if no scaling is applied (i.e., 'scl_A "
            "= scl_B = 1').\n"
            "  - Pseudocounts do not affect this check in either mode.\n\n"
        ),
    )

    parser.add_argument(
        "-dn",
        "--drp_nan",
        "--drp-nan",
        dest="drp_nan",
        action="store_true",
        default=False,
        help=(
            "Drop non-finite values from main output. Values 'nan', 'inf', "
            "and '-inf' are skipped.\n"
            "\n"
            "Note: '--track' output will still drop 'inf', '-inf', and 'nan'."
            "\n\n"
        ),
    )
    parser.add_argument(
        "-tr",
        "--track",
        dest="track",
        action="store_true",
        default=False,
        help=(
            "Write a companion track file where rows with 'inf', '-inf', and "
            "'nan' values are excluded. The new file will have '.track' "
            "before the extension.\n"
            "\n"
            "Notes:\n"
            "  - This file is safe for use in genome browsers such as IGV.\n"
        ),
    )

    parser.add_argument(
        "-dp",
        "--dp",
        dest="dp",
        type=int,
        default=24,
        help=(
            "Maximum number of decimal places retained for finite emitted "
            "values (default: %(default)s). After rounding, non-informative "
            "trailing zeros are stripped.\n\n"
        ),
    )
    parser.add_argument(
        "-sp",
        "--skp_pfx",
        dest="skp_pfx",
        type=str,
        default=",".join(DEF_SKP_PFX),
        help=(
            "Comma-separated list of header prefixes to skip in bedGraph "
            "headers/metadata; to disable skipping, pass an empty string "
            "(default: %(default)s).\n\n"
        ),
    )
    parser.add_argument(
        "-sb",
        "--strict_bins",
        "--strict-bins",
        dest="strict_bins",
        action="store_true",
        default=False,
        help=(
            "Require strict bin compatibility: both input bedGraphs must have "
            "the same ordered '(chrom, start, end)' grid across all data "
            "rows. Without this flag, only the first few paired rows are "
            "checked for equal bin width.\n\n"
        ),
    )

    argv_parse = sys.argv[1:] if argv is None else argv

    if not argv_parse:
        parser.print_help(sys.stderr)
        raise SystemExit(0)

    return parser.parse_args(argv_parse)


def main(argv: list[str] | None = None) -> int:
    """
    Execute the primary control flow for the script.

    Parameters
    ----------
    argv : list[str] | None
        Optional list of command-line arguments. When None (the default),
        'sys.argv[1:]' is used (the typical CLI entry point).

    Returns
    -------
    status : int
        Zero on success after writing the main bedGraph and optional '.track'
        sidecar.

    Raises
    ------
    SystemExit
        For help, invalid paths or numeric arguments, malformed inputs,
        incompatible bins, or read and write failures.

    Notes
    -----
    Verbose argument banners and human-readable failure diagnostics are
    written to stderr.
    """

    args = parse_args(argv)

    if args.fil_A == "-" or args.fil_B == "-":
        raise SystemExit(
            "Dash input is no longer supported; provide file paths for "
            "'--fil_A' and '--fil_B'.",
        )

    if args.fil_out == "-":
        raise SystemExit(
            "Dash output is no longer supported; provide a file path for "
            "'--fil_out'.",
        )

    try:
        check_exists(args.fil_A, "file", "First file (A)")
    except FileNotFoundError as e:
        raise SystemExit(str(e)) from None

    try:
        check_exists(args.fil_B, "file", "Second file (B)")
    except FileNotFoundError as e:
        raise SystemExit(str(e)) from None

    if args.chr_sizes is not None:
        try:
            check_exists(args.chr_sizes, "file", "chr.sizes file")
            chrom_sizes = load_chromosome_sizes(args.chr_sizes)
        except (FileNotFoundError, ValueError, OSError) as e:
            raise SystemExit(str(e)) from None
    else:
        chrom_sizes = None

    try:
        fil_out, _, _ = validate_output_path(args.fil_out, BEDGRAPH_FORMATS)
        check_writable(fil_out, kind="file")
    except (
        ValueError,
        FileNotFoundError,
        PermissionError,
        IsADirectoryError,
    ) as e:
        raise SystemExit(str(e)) from None

    skp_pfx = parse_skp_pfx(args.skp_pfx, default=DEF_SKP_PFX)

    try:
        if args.scl_fct is not None:
            scl_a, scl_b = parse_pair(args.scl_fct, 1.0)
        else:
            scl_a, scl_b = 1.0, 1.0

        psc_a, psc_b = parse_pair(args.pseudo, 0.0)
    except argparse.ArgumentTypeError as e:
        raise SystemExit(str(e)) from None

    if args.skip_00 not in {None, "pre_scale", "post_scale"}:
        raise SystemExit(
            (
                "Invalid '--skip_00'; expected None, 'pre_scale', or "
                "'post_scale'."
            ),
        )

    try:
        for label, value in (("scl_fct:A", scl_a), ("scl_fct:B", scl_b)):
            validate_comparison(value, "gt", 0.0, label, allow_none=False)

        for label, value in (("pseudo:A", psc_a), ("pseudo:B", psc_b)):
            validate_comparison(value, "ge", 0.0, label, allow_none=False)

        validate_comparison(
            args.dep_min,
            "gt",
            0.0,
            "dep_min",
            allow_none=True,
        )
        validate_comparison(args.eps, "ge", 0.0, "eps", allow_none=False)
        validate_comparison(args.dp, "ge", 0, "dp", allow_none=False)

    except ValueError as e:
        raise SystemExit(str(e)) from None

    requested_method = args.method
    args.method = METHOD_CANON[args.method]

    log2 = args.method in {"log2", "log2_r"}
    recip = args.method in {"unadj_r", "log2_r"}

    if args.verbose and args.dep_min is not None and args.dep_min <= args.eps:
        # A clamp at or below epsilon does not prevent undefined ratios.
        print(
            f"Warning: 'dep_min' ({args.dep_min}) <= 'ε' ({args.eps}), so "
            "denominator guard will emit 'nan' whenever '|B| <= ε'. Consider "
            "choosing 'dep_min > ε'.",
            file=sys.stderr,
        )

    if (
        args.verbose
        and args.dep_min is not None
        and (psc_a > 0.0 or psc_b > 0.0)
    ):
        print(
            "Warning: Both '--dep_min' and '--pseudo' were supplied. This is "
            "allowed, but interpretability may be reduced because both "
            "arguments stabilize low-depth ratio behavior in different ways.",
            file=sys.stderr,
        )

    if (
        args.verbose
        and args.skip_00 is not None
        and (psc_a > 0.0 or psc_b > 0.0)
        and args.eps > 0.0
    ):
        print(
            "Warning: '--skip_00', '--eps', and non-zero '--pseudo' were all "
            "supplied. This is allowed, but note that the zero-zero skip is "
            "applied before pseudocount addition, so pseudocounts do not "
            "rescue bins dropped by '--skip_00'.",
            file=sys.stderr,
        )

    if args.verbose:
        with redirect_stdout(sys.stderr):
            print("##########################################")
            print("## Arguments for 'compute_signal_ratio' ##")
            print("##########################################")
            print("")

            print(f"--fil_A    {args.fil_A}")
            print(f"--fil_B    {args.fil_B}")
            print(f"--chr_sizes {args.chr_sizes}")
            print(f"--fil_out  {fil_out}")

            if requested_method != args.method:
                method_display = (
                    f"{requested_method}  (standardized internally to "
                    f"{args.method})"
                )

                print(f"--method   {method_display}")
            else:
                print(f"--method   {args.method}")

            if args.scl_fct is not None:
                print(f"--scl_fct {scl_a}:{scl_b}")

            print(f"--pseudo   {psc_a}:{psc_b}")

            if args.dep_min is not None:
                print(f"--dep_min  {args.dep_min}")

            if args.skip_00:
                print(f"--skip_00  {args.skip_00}")

            print(f"--eps      {args.eps}")

            if args.drp_nan:
                print("--drp_nan")

            print(f"--track    {args.track}")

            print(f"--dp       {args.dp}")
            print(f"--skp_pfx  {skp_pfx}")
            print(f"--strict_bins {args.strict_bins}")

            print("")
            print("")

    try:
        comp_sig_rat(
            fil_a=args.fil_A,
            fil_b=args.fil_B,
            fil_out=fil_out,
            scl_a=scl_a,
            scl_b=scl_b,
            psc_a=psc_a,
            psc_b=psc_b,
            dep_min=args.dep_min,
            decimal_places=args.dp,
            log2=log2,
            recip=recip,
            skip_00=args.skip_00,
            eps=args.eps,
            track=args.track,
            drp_nan=args.drp_nan,
            skp_pfx=skp_pfx,
            strict_bins=args.strict_bins,
            chrom_sizes=chrom_sizes,
        )

    except (ValueError, FileNotFoundError, PermissionError, OSError) as e:
        raise SystemExit(str(e)) from None

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        with suppress(Exception):
            sys.stdout.close()

        with suppress(Exception):
            sys.stderr.close()

        raise SystemExit(0) from None
