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

Usage
-----
python -m protocol_chipseq_signal_norm.cli.compute_signal_ratio [options] \\
    --fil_A <file> --fil_B <file> --fil_out <file>

Parameters
----------
Input bedGraphs, output path, ratio method, scaling, pseudocount, denominator
floor, filtering, and formatting options are parsed by 'parse_args()'.

Returns
-------
Writes a bedGraph-like ratio track and, when requested, a finite-valued
'.track' sidecar suitable for genome browsers.

See Also
--------
docs/dev/signal_stabilization.md
    Maintainer notes on ratio operation order, pseudocounts, denominator
    floors, and performance.
"""

from __future__ import annotations

import argparse
import math
import signal
import sys
from contextlib import nullcontext, redirect_stdout, suppress

from protocol_chipseq_signal_norm.utilities.utils_bdg import (
    check_grid_bin,
    check_size_bin,
    generate_name_track,
    key_bin,
    load_chr_sizes,
    validate_bounds_bdg,
)
from protocol_chipseq_signal_norm.utilities.utils_check import (
    check_cmp,
    check_exists,
    check_parse_fil_out,
    check_writable,
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

#  Set allowed fil_out extensions
ONLY_BDG = ("bedGraph", "bedgraph", "bdg", "bg")

#  Accepted '--method' values and their canonical internal names
METHOD_CANON = {
    #  Simple unadjusted ratio: A / B
    "r": "unadj",
    "raw": "unadj",
    "u": "unadj",
    "unadj": "unadj",
    "unadjusted": "unadj",
    "s": "unadj",
    "smp": "unadj",
    "simple": "unadj",

    #  Log2 ratio: log2(A / B)
    "2": "log2",
    "l2": "log2",
    "lg2": "log2",
    "log2": "log2",

    #  Reciprocal of simple ratio: B / A
    "rr": "unadj_r",
    "raw_r": "unadj_r",
    "ur": "unadj_r",
    "unadj_r": "unadj_r",
    "unadjusted_r": "unadj_r",
    "sr": "unadj_r",
    "smp_r": "unadj_r",
    "simple_r": "unadj_r",

    #  Reciprocal of log2 ratio: log2(B / A) = -log2(A / B)
    "2r": "log2_r",
    "l2r": "log2_r",
    "l2_r": "log2_r",
    "lg2_r": "log2_r",
    "log2_r": "log2_r",
}
METHOD_CHOICES = tuple(METHOD_CANON.keys())


def parse_pair(val: str, def_sec: float) -> tuple[float, float]:
    """
    Parse 'A' or 'A:B' into two floats. If only 'A' is provided, the second
    value defaults to a user-assigned default second value, 'def_sec'.

    Parameters
    ----------
        val : str
            String, e.g., '2.0', '2.0:1.5', etc.
        def_sec : float
            Default for the second value when only 'A' is provided.

    Returns
    -------
        (float, float)
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
    sig_A: float,
    sig_B: float,
    scl_A: float | None,
    scl_B: float | None,
    psc_A: float | None,
    psc_B: float | None,
    dep_min: float | None,
    log2: bool,
    recip: bool,
    skip_00: str | None,
    eps: float = 0.0
) -> float | None:
    """
    Compute 'A / B' (optionally log2-transformed and/or the reciprocal) with
    optional per-file multiplicative scaling and optional pseudocount addition,
    as well as an optional "denominator clamp" (e.g., a "minimum input depth
    value" as described in PMID 40364978).

    See module docstring for more details.

    Parameters
    ----------
        sig_A : float
            First file (e.g., IP) signal for a bin (A).
        sig_B : float
            Second file (e.g., input) signal for a bin (B).
        scl_A : float | None
            Per-file multiplicative scale factor for A. If None or 1.0, treated
            as neutral.
        scl_B : float | None
            Per-file multiplicative scale factor for B. If None or 1.0, treated
            as neutral.
        psc_A : float | None
            Pseudocount added to A (post-scaling). If None or 0.0, skipped.
        psc_B : float | None
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
        ratio | xfrm | -xfrm : float | None
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

    Notes
    -----
        Order of operations:
            1. Optionally skip zero-zero bins: '0 / 0' or 'ε / ε' (deepTools-
               like). ‡
            2. Optionally scale each file.
            3. Optionally skip scaled zero-zero bins:
               '[(sf_A × 0) / (sf_B × 0)]' or '[(sf_A × ε) / (sf_B × ε)]'. ‡
            4. Optionally add pseudocounts.
            5. Optionally clamp denominator by 'dep_min'.
            6. If the denominator is <= ε, treat the bin as undefined.
            7. Divide 'A / B'.
            8. Optionally perform log2 transformation.
            9. Optionally compute reciprocal of #7 (linear) or #8 (log2).

            ‡ Either of optional #1 or optional #3, not both.
    """
    #  1. Optionally skip on '0 / 0' or 'ε / ε' prior to optional scaling
    if (
        skip_00 == "pre_scale"
        and abs(sig_A) <= eps
        and abs(sig_B) <= eps
    ):
        return None

    #  2. Optionally scale each file, skipping neutral factors
    num = sig_A if (scl_A is None or scl_A == 1.0) else (scl_A * sig_A)
    den = sig_B if (scl_B is None or scl_B == 1.0) else (scl_B * sig_B)

    #  3. Optionally skip on '0 / 0' or 'ε / ε' after optional scaling, prior
    #     to optional pseudocount addition
    if skip_00 == "post_scale" and abs(num) <= eps and abs(den) <= eps:
        return None

    #  4. Optionally add pseudocounts, skipping neutral summands
    if psc_A not in (None, 0.0):
        num += psc_A
    if psc_B not in (None, 0.0):
        den += psc_B

    #  5. Optionally clamp the denominator (if a minimum denominator threshold
    #     is specified)
    if dep_min is not None and den < dep_min:
        den = dep_min

    #  6. Perform a zero guard: If the denominator (optionally clamped or not)
    #     remains effectively zero ('|den| <= ε'), treat the bin as undefined
    if abs(den) <= eps:  # Advice: pick 'dep_min > ε' so this never triggers
        return float("nan")

    #  7. Compute ratio
    ratio = num / den

    #  8. Optionally compute a log2 transformation of the ratio
    if log2:
        if ratio > 0.0:
            xfrm = math.log2(ratio)
        elif ratio == 0.0:
            xfrm = float("-inf")
        else:
            xfrm = float("nan")

        #  9. Optionally compute a reciprocal on the log2 value (log space)
        return -xfrm if recip else xfrm

    #  9. Optionally compute a reciprocal on the ratio (linear space)
    if recip:
        return (1.0 / ratio) if ratio != 0.0 else float("inf")

    return ratio


def comp_sig_rat(
    fil_A: str,
    fil_B: str,
    fil_out: str,
    scl_A: float,
    scl_B: float,
    psc_A: float,
    psc_B: float,
    dep_min: float | None,
    dp: int,
    log2: bool,
    recip: bool,
    skip_00: str | None,
    eps: float,
    track: bool,
    drp_nan: bool,
    skp_pfx: tuple[str, ...] = DEF_SKP_PFX,
    strict_bins: bool = False,
    chrom_sizes: dict[str, int] | None = None
) -> None:
    """
    Compute the ratio, log2 ratio, reciprocal ratio, or reciprocal log2 ratio
    of the first ('A'; e.g., IP) and second ('B'; e.g., input) bedGraph tracks,
    with optional per-file multiplicative scaling, addition of pseudocount(s),
    and denominator clamping (e.g., "minimum input depth"). Streams the merge.
    Handles cases where bins are missing in one of the files. The merge matches
    bins by chromosome and start coordinate and therefore assumes that matched
    bins also share the same end coordinate.

    See module docstring for more details.

    Parameters
    ----------
        fil_A : str
            Path to first (e.g., IP) bedGraph file (can be gzipped; A).
        fil_B : str
            Path to second (e.g., input) bedGraph file (can be gzipped; B).
        fil_out : str
            Output file path (if '.gz' extension, gzip compression).
        scl_A : float
            Scale factor for A.
        scl_B : float
            Scale factor for B.
        psc_A : float
            Pseudocount added to A after scaling (0.0 means no addition).
        psc_B : float
            Pseudocount added to B after scaling (0.0 means no addition).
        dep_min : float | None
            Denominator clamp ("minimum input depth") to avoid extreme and/or
            erroneous division.
        dp : int
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

    Returns
    -------
        None. Writes one or two bedGraphs to disk.
            - Standard bedGraph file: Contains all computed ratios unless
              '--drp_nan' is invoked, in which case rows with 'inf', '-inf',
              and 'nan' values are omitted. Finite values are rounded to at
              most 'dp' decimal places and then have non-informative trailing
              zeros stripped.
            - If track is 'True', a second file with '.track' before the
              extension is created, excluding rows with 'inf', '-inf', and
              'nan' values.

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
    """
    if fil_A == "-" or fil_B == "-" or fil_out == "-":
        raise ValueError(
            "Dash input/output is no longer supported; provide real file "
            "paths for '--fil_A', '--fil_B', and '--fil_out'."
        )

    if chrom_sizes is not None:
        validate_bounds_bdg(
            fil_A,
            chrom_sizes,
            skp_pfx=skp_pfx,
            label="first file (file A)"
        )
        validate_bounds_bdg(
            fil_B,
            chrom_sizes,
            skp_pfx=skp_pfx,
            label="second file (file B)"
        )

    #  Check bin size/grid consistency before proceeding
    if strict_bins:
        check_grid_bin(fil_A=fil_A, fil_B=fil_B, skp_pfx=skp_pfx)
    else:
        check_size_bin(fil_A=fil_A, fil_B=fil_B, skp_pfx=skp_pfx)

    #  Generate file name for optional track (insert '.track' before the ext)
    fil_trk = generate_name_track(fil_out) if track else None

    #  Open input/output safely; ensure sidecar (if any) is closed on error
    with (
        open_in(fil_A) as opn_A,
        open_in(fil_B) as opn_B,
        open_out(fil_out) as f_out,
        (open_out(fil_trk) if track else nullcontext()) as f_trk
    ):
        def write_line(chrom, start, end, val):
            """
            Helper function: Write one bedGraph row. Optionally drop non-finite
            values ('inf', '-inf', 'nan') in main output; if enabled, always
            drop non-finite values in '.track' sidecar output.
            """
            if val is None:
                return  # Skip '0 / 0' bins (assessed prior to pseudo addition)

            #  Optionally drop non-finite values in main output:
            if drp_nan and not math.isfinite(val):
                return

            #  Safely format value for main output
            if math.isfinite(val):
                v = round(val, dp)
                if v == 0.0:
                    v = 0.0  # Collapse "-0"

                #  Format with at most 'dp' decimal places, then strip only
                #  non-informative trailing zeros and any trailing decimal
                #  point
                s_val = f"{v:.{dp}f}"
                if "." in s_val:
                    s_val = s_val.rstrip("0").rstrip(".")
                if s_val == "-0":
                    s_val = "0"
            else:
                if math.isnan(val):
                    s_val = "nan"
                else:
                    s_val = "-inf" if val < 0 else "inf"

            f_out.write(f"{chrom}\t{start}\t{end}\t{s_val}\n")

            #  Write track 'sidecar': drop all non-finite values
            if (f_trk is not None) and math.isfinite(val):
                val_trk = round(val, dp)
                if val_trk == 0.0:
                    val_trk = 0.0

                s_trk = f"{val_trk:.{dp}f}"
                if "." in s_trk:
                    s_trk = s_trk.rstrip("0").rstrip(".")
                if s_trk == "-0":
                    s_trk = "0"

                f_trk.write(f"{chrom}\t{start}\t{end}\t{s_trk}\n")

        #  Prime the data streams
        lin_A = read_data_line(fh=opn_A, skp_pfx=skp_pfx)
        lin_B = read_data_line(fh=opn_B, skp_pfx=skp_pfx)

        #  Merge the two bedGraph streams in order; if a bin exists in only one
        #  file, we treat the missing partner as 0.0 for that bin, making it so
        #  that we can still compute 'A / B'
        while lin_A or lin_B:  # Continue until both files are exhausted
            fld_A = lin_A.split() if lin_A else None
            fld_B = lin_B.split() if lin_B else None

            #  Fail fast on malformed bedGraph rows
            if fld_A and len(fld_A) < 4:
                raise ValueError(
                    "Malformed bedGraph line in first file (e.g., IP): "
                    f"{lin_A!r}"
                )
            if fld_B and len(fld_B) < 4:
                raise ValueError(
                    "Malformed bedGraph line in second file (e.g., input): "
                    f"{lin_B!r}"
                )

            #  Parse bedGraph format, handling missing bins
            if fld_A:
                try:
                    chr_A = fld_A[0]
                    start_A = int(fld_A[1])
                    end_A = int(fld_A[2])
                    sig_A = float(fld_A[3])
                except ValueError as e:
                    raise ValueError(
                        "Non-numeric start/end/signal in first-file (e.g., "
                        f"IP) line: {lin_A!r}"
                    ) from e
            else:
                chr_A, start_A, end_A, sig_A = None, None, None, 0.0

            if fld_B:
                try:
                    chr_B = fld_B[0]
                    start_B = int(fld_B[1])
                    end_B = int(fld_B[2])
                    sig_B = float(fld_B[3])
                except ValueError as e:
                    raise ValueError(
                        "Non-numeric start/end/signal in second-file (e.g., "
                        f"input) line: {lin_B!r}"
                    ) from e
            else:
                chr_B, start_B, end_B, sig_B = None, None, None, 0.0

            #  Post-parsing sanity check for widths
            if fld_A and start_A >= end_A:
                raise ValueError(
                    "Non-positive width in first-file (e.g., IP) line: "
                    f"{lin_A!r}"
                )
            if fld_B and start_B >= end_B:
                raise ValueError(
                    "Non-positive width in second-file (e.g., input) line: "
                    f"{lin_B!r}"
                )

            #  Choose which bin to emit next
            if (fld_A is not None) and (fld_B is not None):
                # TODO: The merge key currently uses chromosome and start
                #       coordinate only (via 'key_bin'). For bins treated as
                #       matches ('key_A == key_B'), the end coordinate should
                #       also be checked. When stricter bin validation is added,
                #       enforce matching '(chrom, start, end)' across more (or
                #       all) rows rather than assuming that equal
                #       '(chrom, start)' implies a true bin match.
                key_A = key_bin(chr_A, start_A)
                key_B = key_bin(chr_B, start_B)

                if key_A == key_B:
                    #  Same bin in both files: compute ratio with A and B as-is
                    val = calc_rat_bin(
                        sig_A=sig_A,
                        sig_B=sig_B,
                        scl_A=scl_A,
                        scl_B=scl_B,
                        psc_A=psc_A,
                        psc_B=psc_B,
                        dep_min=dep_min,
                        log2=log2,
                        recip=recip,
                        skip_00=skip_00,
                        eps=eps
                    )
                    write_line(chr_A, start_A, end_A, val)
                    lin_A = read_data_line(fh=opn_A, skp_pfx=skp_pfx)
                    lin_B = read_data_line(fh=opn_B, skp_pfx=skp_pfx)

                elif key_A < key_B:
                    #  File A is "ahead", meaning this bin exists only in A at
                    #  the moment; thus, treat B as 0.0 for this bin and
                    #  compute 'A / 0' subject to guards
                    val = calc_rat_bin(
                        sig_A=sig_A,
                        sig_B=0.0,
                        scl_A=scl_A,
                        scl_B=scl_B,
                        psc_A=psc_A,
                        psc_B=psc_B,
                        dep_min=dep_min,
                        log2=log2,
                        recip=recip,
                        skip_00=skip_00,
                        eps=eps
                    )
                    write_line(chr_A, start_A, end_A, val)
                    lin_A = read_data_line(fh=opn_A, skp_pfx=skp_pfx)

                else:
                    #  File B is "ahead", meaning this bin exists only in B
                    #  right now; thus, treat A as 0.0 for this bin and compute
                    #  '0 / B' subject to guards
                    val = calc_rat_bin(
                        sig_A=0.0,
                        sig_B=sig_B,
                        scl_A=scl_A,
                        scl_B=scl_B,
                        psc_A=psc_A,
                        psc_B=psc_B,
                        dep_min=dep_min,
                        log2=log2,
                        recip=recip,
                        skip_00=skip_00,
                        eps=eps
                    )
                    write_line(chr_B, start_B, end_B, val)
                    lin_B = read_data_line(fh=opn_B, skp_pfx=skp_pfx)

            elif fld_A is not None:
                #  Only A has bins left; B is exhausted, so treat B as 0.0 for
                #  remaining bins
                val = calc_rat_bin(
                    sig_A=sig_A,
                    sig_B=0.0,
                    scl_A=scl_A,
                    scl_B=scl_B,
                    psc_A=psc_A,
                    psc_B=psc_B,
                    dep_min=dep_min,
                    log2=log2,
                    recip=recip,
                    skip_00=skip_00,
                    eps=eps
                )
                write_line(chr_A, start_A, end_A, val)
                lin_A = read_data_line(fh=opn_A, skp_pfx=skp_pfx)

            else:
                #  Only B has bins left; A is exhausted, so treat A as 0.0 for
                #  remaining bins
                val = calc_rat_bin(
                    sig_A=0.0,
                    sig_B=sig_B,
                    scl_A=scl_A,
                    scl_B=scl_B,
                    psc_A=psc_A,
                    psc_B=psc_B,
                    dep_min=dep_min,
                    log2=log2,
                    recip=recip,
                    skip_00=skip_00,
                    eps=eps
                )
                write_line(chr_B, start_B, end_B, val)
                lin_B = read_data_line(fh=opn_B, skp_pfx=skp_pfx)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse command-line arguments.
    """
    parser = CapArgumentParser(
        description=(
            "Compute per-bin ratios between two bedGraphs [A (numerator) and "
            "B (denominator)] with optional per-file multiplicative scaling, "
            "pseudocount addition, 'dep_min' “clamping” (denominator "
            "thresholding, i.e., a divisor “floor”), log2 transformation, "
            "and reciprocal computation.\n"
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
        )
    )
    add_help_cap(parser)
    parser.add_argument(
        "-v", "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Run script in verbose mode.\n\n"
    )
    parser.add_argument(
        "-fA", "--fil_A",
        dest="fil_A",
        required=True,
        type=str,
        help=(
            "First bedGraph input file, file A (e.g., IP). "
            "Supports plain text and '.gz'.\n\n"
        )
    )
    parser.add_argument(
        "-fB", "--fil_B",
        dest="fil_B",
        required=True,
        type=str,
        help=(
            "Second bedGraph input file, file B (e.g., input). "
            "Supports plain text and '.gz'.\n\n"
        )
    )
    parser.add_argument(
        "-fo", "--fil_out",
        dest="fil_out",
        required=True,
        type=str,
        help=(
            "Output file path. Path to the output bedGraph file. "
            "When a real path is provided, accepted extensions are "
            "'.bedGraph', '.bedgraph', '.bdg', and '.bg', each optionally "
            "followed by '.gz'.\n\n"
        )
    )
    parser.add_argument(
        "-cs", "--chr_sizes",
        dest="chr_sizes",
        default=None,
        help=(
            "Chromosome sizes file in UCSC-style TSV format used to validate "
            "input bedGraph interval bounds. This validation is independent "
            "of '--strict_bins'.\n\n"
        )
    )
    parser.add_argument(
        "-me", "--method",
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
        )
    )
    parser.add_argument(
        "-sf", "--scl_fct", "--scl-fct",
        dest="scl_fct",
        type=str,
        default=None,
        help=(
            "Scaling factor. Per-file scale factor(s) 'A[:B]'. If only 'A' is "
            "given, 'B' defaults to 1.0.\n\n"
        )
    )
    parser.add_argument(
        "-ps", "--pseudo",
        dest="pseudo",
        type=str,
        default="0:0",
        help=(
            "Per-file pseudocount spec 'A[:B]' added after scaling (default: "
            "%(default)s)\n"
            "\n"
            "Note: This is primarily useful for log2-ratio methods, where it "
            "helps avoid undefined values such as 'log2(A / 0)'; for linear "
            "ratios, '--dep_min' may be preferable when the main concern is "
            "bounding low-depth extremes.\n"
            "\n"
            "Using '--pseudo' together with '--dep_min' is allowed, but "
            "usually makes low-depth stabilization harder to interpret.\n\n"
        )
    )
    parser.add_argument(
        "-dm", "--dep_min", "--dep-min",
        dest="dep_min",
        type=float,
        default=None,
        help=(
            "Minimum allowed denominator threshold or “clamp” (i.e., “minimum "
            "input depth” as described in PMID 40364978). This is a “floor” "
            "(threshold) ensuring denominators do not fall below this value.\n"
            "\n"
            "'dep_min' is applied after any optional scaling and/or "
            "pseudocount addition ('B := max(B, dep_min)').\n"
            "\n"
            "With or without clamping, a zero-guard with ε ('--eps') is "
            "applied after this optional step in the order of operations; if "
            "'--dep_min <flt>' is specified, choose 'dep_min > ε'.\n"
            "\n"
            "Using '--dep_min' together with '--pseudo' is allowed, but "
            "usually makes low-depth stabilization harder to interpret.\n\n"
        )
    )
    parser.add_argument(
        "-e", "--eps",
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
        )
    )
    parser.add_argument(
        "-s0", "--skip_00", "--skip-00",
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
        )
    )
    parser.add_argument(
        "-dn", "--drp_nan", "--drp-nan",
        dest="drp_nan",
        action="store_true",
        default=False,
        help=(
            "Drop non-finite values from main output. Values 'nan', 'inf', "
            "and '-inf' are skipped.\n"
            "\n"
            "Note: '--track' output will still drop 'inf', '-inf', and "
            "'nan'.\n\n"
        )
    )
    parser.add_argument(
        "-tr", "--track",
        dest="track",
        action="store_true",
        default=False,
        help=(
            "Write a companion track file where rows with 'inf', "
            "'-inf', and 'nan' values are excluded. The new file will have "
            "'.track' before the extension.\n"
            "\n"
            "Notes:\n"
            "  - This file is safe for use in genome browsers such as IGV.\n"
        )
    )
    parser.add_argument(
        "-dp", "--dp",
        dest="dp",
        type=int,
        default=24,
        help=(
            "Maximum number of decimal places retained for finite emitted "
            "values (default: %(default)s). After rounding, "
            "non-informative trailing zeros are stripped.\n\n"
        )
    )
    parser.add_argument(
        "-sp", "--skp_pfx",
        dest="skp_pfx",
        type=str,
        default=",".join(DEF_SKP_PFX),
        help=(
            "Comma-separated list of header prefixes to skip in bedGraph headers/"
            "metadata; to disable skipping, pass an empty string (default: "
            "%(default)s).\n\n"
        )
    )
    parser.add_argument(
        "-sb", "--strict_bins", "--strict-bins",
        dest="strict_bins",
        action="store_true",
        default=False,
        help=(
            "Require strict bin compatibility: both input bedGraphs must have the same ordered "
            "'(chrom, start, end)' grid across all data rows. Without this "
            "flag, only the first few paired rows are checked for equal bin "
            "width.\n\n"
        )
    )

    #  If no arguments are provided, display help and exit
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
        int.
            On success, returns 0 and writes the main bedGraph, along with an
            optional '.track' sidecar when requested.

    Side effects:
        - When '--verbose' is set, prints a human-readable banner of the
          resolved arguments to stderr.
        - Prints error messages to stderr on validation/I/O failures.
        - Exits the interpreter (sys.exit) with status codes below.

    Exits:
        - 0 on success or when showing help with no arguments.
        - 1 on validation or computation errors, including (non-exhaustive):
            + missing inputs or unwritable output directory,
            + malformed bedGraph rows or inconsistent bin sizes,
            + invalid numeric arguments (e.g., 'scl_fct <= 0', 'dep_min <= 0',
              'eps < 0', 'dp < 0'), and
            + read/write I/O errors.
    """
    #  Parse CLI arguments
    args = parse_args(argv)

    if args.fil_A == "-" or args.fil_B == "-":
        raise SystemExit(
            "Dash input is no longer supported; provide file paths for "
            "'--fil_A' and '--fil_B'."
        )
    if args.fil_out == "-":
        raise SystemExit(
            "Dash output is no longer supported; provide a file path for "
            "'--fil_out'."
        )

    #  Validate fil_in existence
    try:
        check_exists(args.fil_A, "file", "First file (A)")
    except FileNotFoundError as e:
        #  Print a one-line message and exit cleanly
        raise SystemExit(str(e)) from None

    try:
        check_exists(args.fil_B, "file", "Second file (B)")
    except FileNotFoundError as e:
        #  Print a one-line message and exit cleanly
        raise SystemExit(str(e)) from None

    if args.chr_sizes is not None:
        try:
            check_exists(args.chr_sizes, "file", "chr.sizes file")
            chrom_sizes = load_chr_sizes(args.chr_sizes)
        except (FileNotFoundError, ValueError, OSError) as e:
            raise SystemExit(str(e)) from None
    else:
        chrom_sizes = None

    #  Check and standardize fil_out (bedGraph only in this script); validate
    #  writability too
    try:
        fil_out, _, _ = check_parse_fil_out(args.fil_out, ONLY_BDG)
        check_writable(fil_out, kind="file")
    except (
        ValueError, FileNotFoundError, PermissionError, IsADirectoryError
    ) as e:
        raise SystemExit(str(e)) from None

    #  Parse header-skip prefixes from CLI (empty string: no skipping)
    skp_pfx = parse_skp_pfx(args.skp_pfx, default=DEF_SKP_PFX)

    #  Parse per-file scale factors and pseudocounts
    try:
        if args.scl_fct is not None:
            scl_A, scl_B = parse_pair(args.scl_fct, 1.0)
        else:
            #  If not specified, default to '1.0:1.0'
            scl_A, scl_B = 1.0, 1.0

        psc_A, psc_B = parse_pair(args.pseudo, 0.0)
    except argparse.ArgumentTypeError as e:
        #  Surface malformed "A" or "A:B" string as a one-line error
        raise SystemExit(str(e)) from None

    #  Validate 'skip_00'
    if args.skip_00 not in {None, "pre_scale", "post_scale"}:
        raise SystemExit(
            "Invalid '--skip_00'; expected None, 'pre_scale', or 'post_scale'."
        )

    #  Validate numeric arguments
    try:
        for lbl, v in (("scl_fct:A", scl_A), ("scl_fct:B", scl_B)):
            check_cmp(v, "gt", 0.0, lbl, allow_none=False)              # > 0
        for lbl, v in (("pseudo:A", psc_A), ("pseudo:B", psc_B)):
            check_cmp(v, "ge", 0.0, lbl, allow_none=False)              # >= 0
        check_cmp(args.dep_min, "gt", 0.0, "dep_min", allow_none=True)  # > 0
        check_cmp(args.eps, "ge", 0.0, "eps", allow_none=False)         # >= 0
        check_cmp(args.dp, "ge", 0, "dp", allow_none=False)           # >= 0

    except ValueError as e:
        raise SystemExit(str(e)) from None

    #  Standardize '--method' to canonical internal name, then derive the
    #  internal computation flags from it
    mthd_in = args.method
    args.method = METHOD_CANON[args.method]

    log2 = args.method in {"log2", "log2_r"}
    recip = args.method in {"unadj_r", "log2_r"}

    #  Warn when the clamp threshold is at/below the zero tolerance
    if (
        args.verbose and args.dep_min is not None and
        args.dep_min <= args.eps
    ):
        #  If 'dep_min <= ε', the post-clamp denominator guard ('|B| <= ε')
        #  will still emit NaNs, effectively negating the clamp
        print(
            f"Warning: 'dep_min' ({args.dep_min}) <= 'ε' ({args.eps}), so "
            "denominator guard will emit 'nan' whenever '|B| <= ε'. Consider "
            "choosing 'dep_min > ε'.",
            file=sys.stderr
        )

    #  Warn when both stabilization strategies are supplied together
    if (
        args.verbose and args.dep_min is not None
        and (psc_A > 0.0 or psc_B > 0.0)
    ):
        print(
            "Warning: Both '--dep_min' and '--pseudo' were supplied. This is "
            "allowed, but interpretability may be reduced because both "
            "arguments stabilize low-depth ratio behavior in different ways.",
            file=sys.stderr
        )

    #  Warn when zero-zero skipping is combined with non-zero pseudocounts
    #  and a positive epsilon value
    if (
        args.verbose
        and args.skip_00 is not None
        and (psc_A > 0.0 or psc_B > 0.0)
        and args.eps > 0.0
    ):
        print(
            "Warning: '--skip_00', '--eps', and non-zero '--pseudo' were all "
            "supplied. This is allowed, but note that the zero-zero skip is "
            "applied before pseudocount addition, so pseudocounts do not "
            "rescue bins dropped by '--skip_00'.",
            file=sys.stderr
        )

    #  Print verbose output
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
            if mthd_in != args.method:
                print(
                    f"--method   {mthd_in}  (standardized internally to "
                    f"{args.method})"
                )
            else:
                print(f"--method   {args.method}")
            if args.scl_fct is not None:
                print(f"--scl_fct {scl_A}:{scl_B}")
            print(f"--pseudo   {psc_A}:{psc_B}")
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

    #  Call function for bin-wise ratio computations
    try:
        comp_sig_rat(
            fil_A=args.fil_A,
            fil_B=args.fil_B,
            fil_out=fil_out,
            scl_A=scl_A,
            scl_B=scl_B,
            psc_A=psc_A,
            psc_B=psc_B,
            dep_min=args.dep_min,
            dp=args.dp,
            log2=log2,
            recip=recip,
            skip_00=args.skip_00,
            eps=args.eps,
            track=args.track,
            drp_nan=args.drp_nan,
            skp_pfx=skp_pfx,
            strict_bins=args.strict_bins,
            chrom_sizes=chrom_sizes
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
