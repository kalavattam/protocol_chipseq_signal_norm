#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: sum_bdg.py
#
# Copyright 2025-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5-series models; most recent: GPT-5.6) were
# used in development and documentation.
#
# Distributed under the MIT license.


"""
Sum values from one or more bedGraph files.

Usage
-----
python -m protocol_chipseq_signal_norm.cli.sum_bdg \\
    [--weight] [--dp <int>] <file> [track_2 ...]

Returns
-------
Prints one tab-delimited 'path' and rounded sum line per input. With
'--weight', each value is multiplied by interval width before summing.

Notes
-----
Input is streamed and may be plain text, gzip-compressed, or stdin.
"""

from __future__ import annotations

import argparse
import signal
import sys
from contextlib import redirect_stdout, suppress

from protocol_chipseq_signal_norm.utilities.utils_bdg import iter_rows_bdg
from protocol_chipseq_signal_norm.utilities.utils_check import (
    check_cmp,
    check_exists,
)
from protocol_chipseq_signal_norm.utilities.utils_cli import (
    CapArgumentParser,
    add_help_cap,
)
from protocol_chipseq_signal_norm.utilities.utils_io import (
    DEF_SKP_PFX,
    ensure_single_stdin,
    is_header,
    open_in,
    parse_skp_pfx,
)

with suppress(AttributeError, ValueError):
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

assert sys.version_info >= (3, 11), "Python >= 3.11 required."


# TODO: Return a warning and exit if a user supplies interactive-mode flags

def sum_bdg(
    path: str, *, weight: bool, skp_pfx: tuple[str, ...]
) -> float:
    """
    Sum (integrate) a bedGraph track.

    Parameters
    ----------
        path : str
            Path to bedGraph file ('.gz' is handled), or '-' for stdin.
        weight : bool
            If True, accumulate 'val * (end - start)'; otherwise just 'val'.
        skp_pfx : tuple[str, ...]
            Header/metadata prefixes to skip (after left-stripping).

    Returns
    -------
        total : float
            The summed total.

    Notes
    -----
        - Lines beginning with any of 'skp_pfx' or blank lines are ignored.
        - Lines with fewer <4 fields, non-integer coordinates, non-positive
          widths, or non-finite values (NaN/±inf) are ignored.
        - Sorting is not required for summation.
    """
    total = 0.0

    #  Build the skip predicate once for this file
    def skp_prd(line: str) -> bool:
        return is_header(line, skp_pfx)

    with open_in(path) as fh:
        for _, s, e, _, v_num in iter_rows_bdg(fh, skp_prd):
            #  Headers, blanks, <4 fields, non-int coordinates, non-positive
            #  spans are filtered by 'iter_rows_bdg'
            if v_num is None:
                #  Skip non-finite or unparsable value
                continue
            total += v_num * (e - s) if weight else v_num

    return total


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse command-line arguments.
    """
    parser = CapArgumentParser(
        description=(
            "Sum (integrate) bedGraph file(s). Use '--weight' when values "
            "are per-base averages (e.g., deepTools CPM/BPM/RPKM/RPGC). For "
            "per-bin totals (e.g., 'compute_signal.py'), '--weight' is not "
            "needed."
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
        "-w", "--weight",
        dest="weight",
        action="store_true",
        default=False,
        help=(
            "Weight values by width (end - start) prior to integration.\n\n"
        )
    )
    parser.add_argument(
        "-dp", "--dp",
        dest="dp",
        type=int,
        default=24,
        help=(
            "Maximum number of decimal places retained for finite emitted "
            "values.(default: %(default)s).\n\n"
        )
    )
    parser.add_argument(
        "-sp", "--skp_pfx",
        dest="skp_pfx",
        type=str,
        default=",".join(DEF_SKP_PFX),
        help=(
            "Comma-separated list of header prefixes to skip in bedGraph "
            "headers/metadata; to disable skipping, pass an empty string "
            "(default: %(default)s).\n\n"
        )
    )
    parser.add_argument(
        "paths",
        nargs="+",
        help=(
            "bedGraph paths ('.gz' is handled), or '-' for stdin (at most, "
            "one '-' is allowed)."
        )
    )

    #  If no arguments are provided, use 'argv' to display help and exit
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
            Optional argument vector for testing. When None (default), uses
            sys.argv[1:].

    Returns
    -------
        int
            Process exit code. Returns 0 on success after printing
            "<path><tab><sum_rounded>" for each input.

            Nonzero codes are produced only by early exits elsewhere (e.g.,
            argument validation before calling this function or SystemExit
            raised by the program entry point).

    Raises
    ------
        SystemExit
            When invoked as a script, the module’s entry point calls
            raise SystemExit(main()). Argparse help (no arguments) also exits.
        ValueError
            Propagated only if lower-level helpers are extended to raise on
            malformed bedGraph lines. In the current implementation, malformed
            lines are skipped during summation.

    Notes
    -----
        - Validates that '--dp' >= 0 and at most one '-' (stdin) path is
          given.
        - Header prefixes are parsed once and reused for all inputs.
        - Each input path is processed independently; results are printed in
          the same order as provided.
        - With '--weight' (and equivalent flags), values are multiplied by
          (end - start) before summing; otherwise, the value column is summed
           as-is.
        - Printed totals are rounded to '--dp' decimal places.
    """
    args = parse_args(argv)

    try:
        #  Perform argument checks
        ensure_single_stdin(args.paths)

        #  Perform existence checks for files (but skip '-')
        for p in args.paths:
            if p != "-":
                check_exists(p, "file", "bedGraph")

        #  Handle scalars for 'dp': must be > 0
        check_cmp(args.dp, "ge", 0, "dp", allow_none=False)
    except ValueError as e:
        raise SystemExit(str(e)) from None

    #  Standardize header-skip prefixes
    skp_pfx = parse_skp_pfx(args.skp_pfx, default=DEF_SKP_PFX)

    #  If '--verbose', print banner
    if args.verbose:
        with redirect_stdout(sys.stderr):
            print("########################")
            print("## Arguments: sum_bdg ##")
            print("########################")
            print("")
            print("--verbose")
            if args.weight:
                print("--weight")
            print(f"--dp     {args.dp}")
            print(f"--skp_pfx {skp_pfx}")
            for p in args.paths:
                print(f"  {p}")
            print("")
            print("")

    #  Process each input independently
    rc = 0
    for p in args.paths:
        try:
            total = sum_bdg(p, weight=args.weight, skp_pfx=skp_pfx)
            v = round(total, args.dp)
            if v == 0.0:  # Avoid printing “-0.0”
                v = 0.0
            print(f"{p}\t{v:.{args.dp}f}")
        except SystemExit:
            #  Re-raise explicit exits from helpers if any
            raise
        except Exception as e:
            print(f"Error processing '{p}': {e}", file=sys.stderr)
            rc = 1

    if rc:
        raise SystemExit(rc)

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
