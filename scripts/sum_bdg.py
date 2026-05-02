#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2025 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Distributed under the MIT license.

"""
Script
------
sum_bdg.py


Description
-----------
Sum (integrate) values from one or more bedGraph files.

By default, this integrates the fourth column as-is, which is appropriate when
the bedGraph stores per-bin totals as in, e.g., output from 'compute_signal.py'
in the Bio-protocol workflow.

If the bedGraph stores per-base averages/densities as in, e.g., deepTools
CPM/BPM/RPKM/RPGC do, use '--weight' to weight each value by its width
(end - start) before summing.

Supports plain or gzip-compressed bedGraphs as well as stdin (via '-').


Usage
-----
python -m scripts.sum_bdg \
    [--verbose] [--weight] [--rnd <int>] <track_1.bdg[.gz]|-> [track_2 ...]


Output
------
For each input path, prints one tab-delimited line to stdout:
    <path> <tab> <sum_rounded>

- <path> is the literal argument provided (or '-' for stdin).
- <sum_rounded> is rounded to '--rnd' decimal places (default: 24).
- With '--weight', each value is multiplied by its interval width (end-start)
  before summing (appropriate for per-base densities like deepTools CPM/BPM/
  {F,R}PKM/{F,R}PGC). Without '--weight', the values are summed as-is
  (appropriate for per-bin totals as in 'compute_signal.py').

Examples:
    ip.totals.bdg.gz  123456.000000000000000000000000
    -                 0.123456789012345678901234


Examples
--------
1. Sum per-bin totals (no width weighting)
'''bash
python -m scripts.sum_bdg ip.totals.bdg.gz in.totals.bdg.gz
'''

2. Sum per-base densities (weight by width)
'''bash
python -m scripts.sum_bdg --weight ip.cpm.bdg.gz in.cpm.bdg.gz
'''

3. Pipe from stdin with 8 decimal places
'''bash
zcat ip.rpkm.bdg.gz | python -m scripts.sum_bdg --weight --rnd 8 -
'''

4. Pipe multiple files found via 'find' (each printed on its own line)
'''bash
find tracks/ -name "*.bdg.gz" -print0 \
    | xargs -0 -n 1 python -m scripts.sum_bdg
'''

5. Use process substitution to sum two bedGraph files in one command
'''bash
python -m scripts.sum_bdg <(zcat A.bdg.gz) <(zcat B.bdg.gz)
'''


General notes
-------------
- Sorting is not required for summation.
- If multiple files are provided, each is processed independently and reported
  on its own line.
- At most one '-' (stdin) path is allowed; enforced by 'ensure_single_stdin'.
- '--rnd' must be >= 0; this is validated via 'check_cmp'.
- Lines beginning with any of the header prefixes (defaults from 'DEF_SKP_PFX':
  '#', 'track', 'browser') are skipped; blank lines are also skipped.
- Malformed data lines are skipped by the iterator 'iter_bdg_rows': fewer than
  4 fields, non-integer coordinates, non-positive widths, or
  non-finite/unparsable value fields (yield 'v_num is None').


Performance notes
-----------------
- Streaming, single-pass implementation; memory usage is O(1) with respect to
  file size.
- Throughput is I/O-bound; due to decompression, gzip-compressed inputs are
  slower to process than plain text.
- To improve wall time, colocate data on fast storage.
"""

from __future__ import annotations

from contextlib import redirect_stdout

from scripts.functions.utils_bdg import iter_bdg_rows
from scripts.functions.utils_check import (
    check_cmp,
    check_exists
)
from scripts.functions.utils_cli import (
    add_help_cap,
    CapArgumentParser
)
from scripts.functions.utils_io import (
    DEF_SKP_PFX,
    ensure_single_stdin,
    is_header,
    open_in,
    parse_skp_pfx
)

import argparse
import signal
import sys

try:
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
except (AttributeError, ValueError):
    pass

assert sys.version_info >= (3, 10), "Python >= 3.10 required."


# TODO: Return a warning and exit if a user supplies interactive-mode flags

def sum_bdg(
    path: str, *, weight: bool, skp_pfx: tuple[str, ...]
) -> float:
    """
    Sum (integrate) a bedGraph track.

    Args:
        path : str
            Path to bedGraph file ('.gz' is handled), or '-' for stdin.
        weight : bool
            If True, accumulate 'val * (end - start)'; otherwise just 'val'.
        skp_pfx : tuple[str, ...]
            Header/metadata prefixes to skip (after left-stripping).

    Returns:
        total : float
            The summed total.

    Notes:
        - Lines beginning with any of 'skp_pfx' or blank lines are ignored.
        - Lines with fewer <4 fields, non-integer coordinates, non-positive
          widths, or non-finite values (NaN/±inf) are ignored.
        - Sorting is not required for summation.
    """
    total = 0.0

    #  Build the skip predicate once for this file
    skp_prd = (lambda line: is_header(line, skp_pfx))

    with open_in(path) as fh:
        for _, s, e, _, v_num in iter_bdg_rows(fh, skp_prd):
            #  Headers, blanks, <4 fields, non-int coordinates, non-positive
            #  spans are filtered by 'iter_bdg_rows'
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
        action="store_true",
        default=False,
        help="Increase output verbosity.\n\n"
    )
    parser.add_argument(
        "-w", "--weight", "--width",
        "--by-weight", "--by_weight",
        "--by-width",  "--by_width",
        action="store_true",
        default=False,
        help=(
            "Weight values by width (end - start) prior to integration.\n\n"
        )
    )
    parser.add_argument(
        "-r", "--rnd",
        type=int,
        default=24,
        help=(
            "Number of decimals in the printed total (default: "
            "%(default)s).\n\n"
        )
    )
    parser.add_argument(
        "-sp", "--skp_pfx",
        type=str,
        default=",".join(DEF_SKP_PFX),
        help=(
            "Comma-separated list of bedGraph prefixes to skip as headers/"
            "metadata; to disable skipping, pass an empty string (default: "
            "%(default)s).\n\n"
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
        sys.exit(0)

    return parser.parse_args(argv_parse)


def main(argv: list[str] | None = None) -> int:
    """
    Execute the primary control flow for the script.

    Args:
        argv : list[str] | None
            Optional argument vector for testing. When None (default), uses
            sys.argv[1:].

    Returns:
        int
            Process exit code. Returns 0 on success after printing
            "<path><tab><sum_rounded>" for each input.

            Nonzero codes are produced only by early exits elsewhere (e.g.,
            argument validation before calling this function or SystemExit
            raised by the program entry point).

    Raises:
        SystemExit
            When invoked as a script, the module’s entry point calls
            raise SystemExit(main()). Argparse help (no arguments) also exits.
        ValueError
            Propagated only if lower-level helpers are extended to raise on
            malformed bedGraph lines. In the current implementation, malformed
            lines are skipped during summation.

    Notes:
        - Validates that '--rnd' >= 0 and at most one '-' (stdin) path is
          given.
        - Header prefixes are parsed once and reused for all inputs.
        - Each input path is processed independently; results are printed in
          the same order as provided.
        - With '--weight' (and equivalent flags), values are multiplied by
          (end - start) before summing; otherwise, the value column is summed
           as-is.
        - Printed totals are rounded to '--rnd' decimal places.
    """
    args = parse_args(argv)

    try:
        #  Perform argument checks
        ensure_single_stdin(args.paths)

        #  Perform existence checks for files (but skip '-')
        for p in args.paths:
            if p != "-":
                check_exists(p, "file", "bedGraph")

        #  Handle scalars for 'rnd': must be > 0
        check_cmp(args.rnd, "ge", 0, "rnd", allow_none=False)
    except ValueError as e:
        raise SystemExit(str(e))

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
            print(f"--rnd     {args.rnd}")
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
            v = round(total, args.rnd)
            if v == 0.0:  # Avoid printing “-0.0”
                v = 0.0
            print(f"{p}\t{v:.{args.rnd}f}")
        except SystemExit:
            #  Re-raise explicit exits from helpers if any
            raise
        except Exception as e:
            print(f"Error processing '{p}': {e}", file=sys.stderr)
            rc = 1

    if rc:
        raise SystemExit(rc)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except BrokenPipeError:
        try:
            sys.stdout.close()
        except Exception:
            pass
        try:
            sys.stderr.close()
        except Exception:
            pass
        sys.exit(0)
