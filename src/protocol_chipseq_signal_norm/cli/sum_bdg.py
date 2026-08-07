#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: sum_bdg.py
#
# Copyright 2025-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5-series models; most recent: GPT-5.6) were
# used in design, development, and documentation, with all output reviewed,
# edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Sum values from one or more bedGraph files.

The CLI prints one tab-delimited 'path' and rounded sum line per input. With
'--weight', each value is multiplied by interval width before summing.

Notes
-----
Input is streamed and may be plain text, gzip-compressed, or stdin.

Examples
--------
python -m protocol_chipseq_signal_norm.cli.sum_bdg \\
    [--weight] [--dp <int>] <file> [track_2 ...]
"""

from __future__ import annotations

import argparse
import signal
import sys
from contextlib import redirect_stdout, suppress

from protocol_chipseq_signal_norm.utilities.utils_bdg import iter_rows_bdg
from protocol_chipseq_signal_norm.utilities.utils_check import (
    check_exists,
    validate_comparison,
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


# TODO: Warn and exit when a user supplies interactive-mode flags.


def sum_bdg(path: str, *, weight: bool, skp_pfx: tuple[str, ...]) -> float:
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

    def skip_predicate(line: str) -> bool:
        """
        Return whether one input line is header or metadata content.
        """

        return is_header(line, skp_pfx)

    with open_in(path) as handle:
        for (
            _chromosome,
            start,
            end,
            _token,
            numeric_value,
        ) in iter_rows_bdg(handle, skip_predicate):
            if numeric_value is None:
                continue

            total += numeric_value * (end - start) if weight else numeric_value

    return total


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
        Parsed bedGraph inputs and summation options.

    Raises
    ------
    SystemExit
        If argument parsing fails or help is requested.
    """

    parser = CapArgumentParser(
        description=(
            "Sum (integrate) bedGraph file(s). Use '--weight' when values are "
            "per-base averages (e.g., deepTools CPM/BPM/RPKM/RPGC). For "
            "per-bin totals (e.g., 'compute_signal.py'), '--weight' is not "
            "needed."
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
        "-w",
        "--weight",
        dest="weight",
        action="store_true",
        default=False,
        help=(
            "Weight values by width (end - start) prior to integration.\n\n"
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
            "values.(default: %(default)s).\n\n"
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
        "paths",
        nargs="+",
        help=(
            "bedGraph paths ('.gz' is handled), or '-' for stdin (at most, "
            "one '-' is allowed)."
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
        Optional argument vector for testing. When None (default), uses
        sys.argv[1:].

    Returns
    -------
    status : int
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
        ensure_single_stdin(args.paths)

        for p in args.paths:
            if p != "-":
                check_exists(p, "file", "bedGraph")

        validate_comparison(args.dp, "ge", 0, "dp", allow_none=False)
    except ValueError as e:
        raise SystemExit(str(e)) from None

    skp_pfx = parse_skp_pfx(args.skp_pfx, default=DEF_SKP_PFX)

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

            for path in args.paths:
                print(f"  {path}")

            print("")
            print("")

    return_code = 0

    for path in args.paths:
        try:
            total = sum_bdg(
                path,
                weight=args.weight,
                skp_pfx=skp_pfx,
            )
            rounded_value = round(total, args.dp)

            if rounded_value == 0.0:  # This avoids printing “-0.0”.
                rounded_value = 0.0

            print(f"{path}\t{rounded_value:.{args.dp}f}")
        except SystemExit:
            raise
        except Exception as error:
            print(f"Error processing '{path}': {error}", file=sys.stderr)
            return_code = 1

    if return_code:
        raise SystemExit(return_code)

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
