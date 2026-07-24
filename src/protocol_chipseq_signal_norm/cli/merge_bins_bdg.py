#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: merge_bins_bdg.py
#
# Copyright 2025-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5-series models; most recent: GPT-5.6) were
# used in development and documentation.
#
# Distributed under the MIT license.


"""
Merge adjacent same-valued bedGraph bins.

Usage
-----
python -m protocol_chipseq_signal_norm.cli.merge_bins_bdg \\
    --fil_in <file> --fil_out <file> [--dp <int>] [--eps <flt>]

Parameters
----------
Input/output paths, rounding, numeric tolerance, and header-prefix options are
parsed by 'parse_args()'.

Returns
-------
Writes a merged bedGraph stream or file. Returns 0 on success and 1 for
validated CLI errors.

Notes
-----
Merging is streaming and only combines contiguous intervals on the same
chromosome when values match under the selected comparison mode.
"""

from __future__ import annotations

import argparse
import math
import signal
import sys
from contextlib import redirect_stdout, suppress
from typing import TextIO

from protocol_chipseq_signal_norm.utilities.utils_bdg import (
    canon_nonfinite,
    try_float,
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
    is_header,
    open_in,
    open_out,
    parse_skp_pfx,
)

with suppress(AttributeError, ValueError):
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

assert sys.version_info >= (3, 11), "Python >= 3.11 required."


# TODO: Add focused tests for exact ties, near-ties, non-finite tokens,
#       combined comparison modes, and mid-stream header lines that restart a
#       merge run.
def format_value(
    tok: str, dp: int | None, cache: float | None
) -> str:
    """
    Decide what to print for a run’s value.

    Parameters
    ----------
        tok : str
            Original token from column 4 for the first bin in the run.
        dp : int | None
            If not None, number of decimals to print when a finite numeric
             cache is available.
        cache : float | None
            Parsed float for the first token in the run (if parseable).

    Returns
    -------
        tok | nf : str
            String to print in column 4 for the merged interval.

    Notes
    -----
        - If 'dp' is set and 'cache' is finite, print formatted to 'dp'.
        - Else if the original token is a recognized non-finite, return the
          canonical token.
        - Else return the original token unchanged.
    """
    nf = canon_nonfinite(tok)
    if (dp is not None) and (cache is not None) and math.isfinite(cache):
        v = round(cache, dp)
        if v == 0.0:  # Avoid printing “-0.0”
            v = 0.0
        return f"{v:.{dp}f}"

    if nf is not None:
        return nf

    return tok


def validate_merge_options(dp: int | None, eps: float | None) -> None:
    """
    Validate merge-mode options for CLI and programmatic callers.
    """
    check_cmp(dp, "ge", 0, "dp", allow_none=True)
    check_cmp(eps, "ge", 0, "eps", allow_none=True)


def merge_bins(
    fil_in: str,
    fil_out: str,
    *,
    dp: int | None = None,
    eps: float | None = None,
    skp_pfx: tuple[str, ...]
) -> None:
    """
    Stream through a (sorted) bedGraph, merging adjacent bins when values
    match.

    Parameters
    ----------
        fil_in : str
            Input bedGraph path ('.gz' ok) or '-' for stdin.
        fil_out : str
            Output bedGraph path ('.gz' ok) or '-' for stdout.
        dp : int | None = None
            If not None, compare by rounded-to-N-decimal strings and print to N
            decimals.
        eps : float | None = None
            If not None and 'dp' is None, compare numerically within
            |Δ| <= eps.
        skp_pfx : tuple[str, ...]
            Header/metadata prefixes to skip (after left-stripping).

    Returns
    -------
        None (writes merged bedGraph to 'fil_out').

    Raises
    ------
        ValueError
            If '--dp' < 0, if '--eps' < 0, or on malformed data lines
            (non-integer coords, non-positive width, or < 4 fields). Headers
            pass through.

    Notes
    -----
        - Matching modes (priority):
            + if dp is not None: compare on rounded string to 'dp' places.
            + elif eps is not None: compare numeric values within |Δ| <= eps.
            + else: compare exact text of the 4th field (after canonicalizing
              nan/inf).
        - Non-finite values ('nan', 'inf', '-inf') never merge with anything
          else.
    """
    validate_merge_options(dp, eps)

    #  Set the “rolling state” for the current merge “run;” a run extends while
    #  (a) chrom stays the same, (b) bins are contiguous, and (c) the value
    #  matches under the chosen rule (rounded/eps/text)
    last_chrom: str | None = None
    last_start: int | None = None
    last_end: int | None = None
    last_tok: str | None = None    # Original token for printing if needed
    last_num: float | None = None  # Numeric cache (float) if parseable
    last_key: tuple | None = None  # ("rounded", str) | ("eps", float) | ("text", str) | ("nonfin", "nan"/"inf"/"-inf")
    last_nonfin: bool = False      # True if last value is non-finite

    def flush(fo: TextIO) -> None:
        """
        Emit the current run (if any) to 'fo' and reset the rolling state.
        """
        nonlocal last_chrom, last_start, last_end
        nonlocal last_tok, last_num, last_key, last_nonfin
        if last_chrom is None:
            return
        out_val = format_value(last_tok, dp, last_num)
        fo.write(f"{last_chrom}\t{last_start}\t{last_end}\t{out_val}\n")
        last_chrom = last_start = last_end = last_tok = None
        last_num = None
        last_key = None
        last_nonfin = False

    #  Stream input to output; headers are forwarded verbatim and terminate
    #  runs
    with open_in(fil_in) as fi, open_out(fil_out) as fo:
        for raw in fi:
            #  Headers pass through and break runs
            if is_header(raw, skp_pfx):
                flush(fo)
                fo.write(raw)
                continue

            #  Blank/whitespace-only lines also pass through and break runs
            if not raw.strip():
                flush(fo)
                fo.write(raw)
                continue

            #  Validate data lines (need at least 4 fields)
            parts = raw.rstrip("\n").split()
            if len(parts) < 4:
                flush(fo)
                raise ValueError(
                    f"Malformed bedGraph line; needs >=4 fields: {raw!r}"
                )

            #  Parse and check coordinates: reject non-integers and
            #  non-positive spans, as bedGraph requires start < end
            chrom, s_str, e_str, v_str = parts[0], parts[1], parts[2], parts[3]
            try:
                s = int(s_str)
                e = int(e_str)
            except ValueError as error:
                flush(fo)
                raise ValueError(
                    f"Non-integer start/end: {raw!r}"
                ) from error
            if e <= s:
                flush(fo)
                raise ValueError(f"Non-positive interval width: {raw!r}")

            #  Build comparison key
            nf = canon_nonfinite(v_str)
            v_num = None if nf is not None else try_float(v_str)

            # MAYBE: extract the “build comparison key” block into a small
            #        helper function
            if nf is not None:
                #  Never merges across lines
                cmp_key = ("nonfin",)
            elif (
                (dp is not None) and
                (v_num is not None) and
                math.isfinite(v_num)
            ):
                #  Round once, using the rounded value for both comparisons and
                #  printing
                v_dp = round(v_num, dp)
                if eps is not None:
                    #  Both '--dp' and '--eps': compare on rounded floats
                    #  within eps
                    cmp_key = ("dp_eps", v_dp)
                else:
                    # '--dp' only: compare on rounded-string equality
                    cmp_key = ("rounded", f"{v_dp:.{dp}f}")
            else:
                if (
                    (eps is not None) and
                    (v_num is not None) and
                    math.isfinite(v_num)
                ):
                    #  '--eps' only: Compare on raw numeric values
                    cmp_key = ("eps", v_num)
                else:
                    #  Default: Compare on exact token text (non-finites
                    #  already handled)
                    cmp_key = ("text", v_str)

            #  Start a new run if (a) chrom changes, (b) bins are not
            #  contiguous, or (c) the mode-dependent comparison key no longer
            #  matches
            start_new = (
                last_chrom is None or
                chrom != last_chrom or
                s != last_end or
                nf is not None or last_nonfin or
                last_key is None or
                cmp_key[0] != last_key[0] or
                (
                    cmp_key[0] == "rounded" and
                    last_key[1] != cmp_key[1]
                ) or
                (
                    cmp_key[0] == "text" and
                    last_key[1] != cmp_key[1]
                ) or
                (
                    cmp_key[0] == "eps" and
                    not math.isclose(
                        last_key[1],
                        cmp_key[1],
                        rel_tol=0.0,
                        abs_tol=(eps or 0.0) + sys.float_info.epsilon,
                    )
                ) or
                (
                    cmp_key[0] == "dp_eps" and
                    not math.isclose(
                        last_key[1],
                        cmp_key[1],
                        rel_tol=0.0,
                        abs_tol=(eps or 0.0) + sys.float_info.epsilon,
                    )
                )
            )

            if start_new:
                flush(fo)
                last_chrom = chrom
                last_start = s
                last_end = e
                last_tok = v_str
                last_num = v_num
                last_nonfin = (nf is not None)
                last_key = cmp_key
            else:
                #  Extend current run
                last_end = e

        #  “Flush” the tail
        flush(fo)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse command-line arguments.
    """
    parser = CapArgumentParser(
        description="Merge adjacent same-valued bedGraph bins (streaming)."
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
        "-fi", "--fil_in",
        dest="fil_in",
        required=True,
        help=(
            "Input file path. Path to the input bedGraph file (.gz is"
            "handled), or '-' for stdin.\n\n"
        )
    )
    parser.add_argument(
        "-fo", "--fil_out",
        dest="fil_out",
        required=True,
        help=(
            "Output file path. Path to the output bedGraph file (.gz is"
            "handled), or '-' for stdout.\n\n"
        )
    )
    parser.add_argument(
        "-dp", "--dp",
        dest="dp",
        type=int,
        default=None,
        help=(
            "Maximum number of decimal places retained for finite emitted"
            "values. Number of decimals in the printed total.\n"
            "\n"
            "Controls the value comparisons by N-decimal rounded strings. "
            "With '--eps', rounding is applied first, then |Δ| (absolute "
            "difference) <= eps is tested on the rounded numbers.\n\n"
        )
    )
    parser.add_argument(
        "-e", "--eps",
        dest="eps",
        type=float,
        default=None,
        help=(
            "Zero tolerance epsilon for matching (merge if |Δ| <= eps).\n"
            "  - Alone: compare raw numeric values with |Δ| <= eps.\n"
            "  - With '--dp N': compare on rounded values (to N decimals) "
            "with |Δ| <= eps.\n\n"
        )
    )
    parser.add_argument(
        "-sp", "--skp_pfx",
        dest="skp_pfx",
        type=str,
        default=",".join(DEF_SKP_PFX),
        help=(
            "Comma-separated list of header prefixes to skip in bedGraph"
            "headers/metadata; to disable skipping, pass an empty string"
            "(default: %(default)s).\n\n"
        )
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
        argv : list[str] | None = None
            Optional argument vector for testing; defaults to sys.argv[1:].

    Returns
    -------
        Process exit code (int).

    Raises
    ------
        ValueError
            Propagated from 'merge_bins' when input data lines are malformed.
    """
    args = parse_args(argv)

    #  Perform argument checks
    try:
        if args.fil_in != "-":
            check_exists(args.fil_in, "file", "bedGraph")

        if args.fil_out != "-":
            fil_out, _, _ = check_parse_fil_out(
                args.fil_out, ("bedgraph", "bdg", "bg")
            )
            check_writable(fil_out, "file")
        else:
            fil_out = "-"

    except ValueError as e:
        raise SystemExit(str(e)) from None

    #  Standardize header-skip prefixes
    skp_pfx = parse_skp_pfx(args.skp_pfx, default=DEF_SKP_PFX)

    #  Print verbose output
    if args.verbose:
        with redirect_stdout(sys.stderr):
            print("####################################")
            print("## Arguments for 'merge_bins_bdg' ##")
            print("####################################")
            print("")
            print(f"--verbose {args.verbose}")
            print(f"--fil_in  {args.fil_in}")
            print(f"--fil_out {fil_out}")
            if args.dp is not None:
                print(f"--dp     {args.dp}")
            if args.eps is not None:
                print(f"--eps     {args.eps}")
            if args.skp_pfx:
                print(f"--skp_pfx {skp_pfx}")
            print("")
            print("")

    try:
        merge_bins(
            fil_in=args.fil_in,
            fil_out=args.fil_out,
            dp=args.dp,
            eps=args.eps,
            skp_pfx=skp_pfx
        )
    except (ValueError, FileNotFoundError, PermissionError, OSError) as e:
        print(str(e), file=sys.stderr)
        raise SystemExit(1) from None

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
