#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: merge_bins_bdg.py
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
Merge adjacent same-valued bedGraph bins.

The CLI accepts input and output paths, rounding, numeric tolerance, and
header-prefix options. It writes a merged bedGraph stream or file and exits
with status 0 on success or 1 for validated CLI errors.

Notes
-----
Merging is streaming and only combines contiguous intervals on the same
chromosome when values match under the selected comparison mode.

Examples
--------
python -m protocol_chipseq_signal_norm.cli.merge_bins_bdg \\
    --fil_in <file> --fil_out <file> [--dp <int>] [--eps <flt>]
"""

from __future__ import annotations

import argparse
import math
import signal
import sys
from contextlib import redirect_stdout, suppress
from dataclasses import dataclass
from typing import TextIO

from protocol_chipseq_signal_norm.utilities.utils_bdg import (
    canon_nonfinite,
    try_float,
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
    is_header,
    open_in,
    open_out,
    parse_skp_pfx,
)

with suppress(AttributeError, ValueError):
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

assert sys.version_info >= (3, 11), "Python >= 3.11 required."


# TODO: Add focused tests for exact ties, near-ties, non-finite tokens,
# combined comparison modes, and mid-stream header lines that restart a
# merge run.
def format_value(
    token: str,
    decimal_places: int | None,
    cache: float | None,
) -> str:
    """
    Decide what to print for a run’s value.

    Parameters
    ----------
    token : str
        Original token from column 4 for the first bin in the run.
    decimal_places : int | None
        If not None, number of decimals to print when a finite numeric cache is
        available.
    cache : float | None
        Parsed float for the first token in the run (if parseable).

    Returns
    -------
    text : str
        String to print in column 4 for the merged interval.

    Notes
    -----
    - If 'decimal_places' is set and 'cache' is finite, use that precision.
    - Else if the original token is a recognized non-finite, return the
      canonical token.
    - Else return the original token unchanged.
    """

    nonfinite_name = canon_nonfinite(token)

    if (
        decimal_places is not None
        and cache is not None
        and math.isfinite(cache)
    ):
        rounded_value = round(cache, decimal_places)

        # Avoid rendering negative zero.
        if rounded_value == 0.0:
            rounded_value = 0.0

        return f"{rounded_value:.{decimal_places}f}"

    if nonfinite_name is not None:
        return nonfinite_name

    return token


def validate_merge_options(
    decimal_places: int | None,
    eps: float | None,
) -> None:
    """
    Validate merge-mode options for CLI and programmatic callers.
    """

    validate_comparison(decimal_places, "ge", 0, "dp", allow_none=True)
    validate_comparison(eps, "ge", 0, "eps", allow_none=True)


@dataclass(frozen=True)
class _ComparedValue:
    """
    Hold one parsed bedGraph value and its merge-comparison key.
    """

    source: str
    numeric: float | None
    nonfinite: bool
    key: tuple[str | float, ...]


def _comparison_value(
    source: str,
    decimal_places: int | None,
    eps: float | None,
) -> _ComparedValue:
    """
    Parse one value and derive its mode-specific comparison key.
    """

    nonfinite_name = canon_nonfinite(source)
    numeric = None if nonfinite_name is not None else try_float(source)

    if nonfinite_name is not None:
        key: tuple[str | float, ...] = ("nonfin",)
    elif (
        decimal_places is not None
        and numeric is not None
        and math.isfinite(numeric)
    ):
        rounded = round(numeric, decimal_places)
        key = (
            ("dp_eps", rounded)
            if eps is not None
            else ("rounded", f"{rounded:.{decimal_places}f}")
        )
    elif eps is not None and numeric is not None and math.isfinite(numeric):
        key = ("eps", numeric)
    else:
        key = ("text", source)

    return _ComparedValue(
        source=source,
        numeric=numeric,
        nonfinite=nonfinite_name is not None,
        key=key,
    )


def _comparison_keys_match(
    current: tuple[str | float, ...],
    previous: tuple[str | float, ...],
    eps: float | None,
) -> bool:
    """
    Return whether two finite merge keys belong to the same run.
    """

    mode = current[0]
    if mode != previous[0]:
        return False

    if mode in {"rounded", "text"}:
        return current[1] == previous[1]

    if mode in {"eps", "dp_eps"}:
        current_value = current[1]
        previous_value = previous[1]
        numeric_values = isinstance(
            current_value, (int, float)
        ) and isinstance(previous_value, (int, float))

        if not numeric_values:
            return False

        return math.isclose(
            previous_value,
            current_value,
            rel_tol=0.0,
            abs_tol=(eps or 0.0) + sys.float_info.epsilon,
        )

    return False


def merge_bins(
    fil_in: str,
    fil_out: str,
    *,
    decimal_places: int | None = None,
    eps: float | None = None,
    skp_pfx: tuple[str, ...],
) -> None:
    """
    Merge adjacent equal-valued bins in a sorted bedGraph.

    Parameters
    ----------
    fil_in : str
        Input bedGraph path ('.gz' ok) or '-' for stdin.
    fil_out : str
        Output bedGraph path ('.gz' ok) or '-' for stdout.
    decimal_places : int | None = None
        If not None, compare by rounded-to-N-decimal strings and print to N
        decimals.
    eps : float | None = None
        If not None and 'dp' is None, compare numerically within |Δ| <= eps.
    skp_pfx : tuple[str, ...]
        Header/metadata prefixes to skip (after left-stripping).

    Raises
    ------
    ValueError
        If '--dp' < 0, if '--eps' < 0, or on malformed data lines (non-integer
        coords, non-positive width, or < 4 fields). Headers pass through.

    Notes
    -----
    The function writes the merged bedGraph to 'fil_out'.

    - Matching modes (priority):
        + if decimal_places is not None: compare rounded decimal strings.
        + elif eps is not None: compare numeric values within |Δ| <= eps.
        + else: compare exact text of the 4th field (after canonicalizing
          nan/inf).
    - Non-finite values ('nan', 'inf', '-inf') never merge with anything else.
    """

    validate_merge_options(decimal_places, eps)

    # A run extends across contiguous, equal-valued bins on one chromosome.
    last_chrom: str | None = None
    last_start: int | None = None
    last_end: int | None = None
    last_token: str | None = None
    last_number: float | None = None
    last_key: tuple[str | float, ...] | None = None
    last_is_nonfinite = False

    def flush(output_handle: TextIO) -> None:
        """
        Emit the current run, if present, and reset the rolling state.
        """

        nonlocal last_chrom, last_start, last_end
        nonlocal last_token, last_number, last_key, last_is_nonfinite

        if last_chrom is None:
            return

        assert last_start is not None
        assert last_end is not None
        assert last_token is not None

        output_value = format_value(
            last_token,
            decimal_places,
            last_number,
        )
        output_handle.write(
            f"{last_chrom}\t{last_start}\t{last_end}\t{output_value}\n",
        )
        last_chrom = last_start = last_end = last_token = None
        last_number = None
        last_key = None
        last_is_nonfinite = False

    with (
        open_in(fil_in) as input_handle,
        open_out(
            fil_out,
        ) as output_handle,
    ):
        for raw in input_handle:
            if is_header(raw, skp_pfx):
                flush(output_handle)
                output_handle.write(raw)

                continue

            if not raw.strip():
                flush(output_handle)
                output_handle.write(raw)

                continue

            parts = raw.rstrip("\n").split()

            if len(parts) < 4:
                flush(output_handle)

                raise ValueError(
                    f"Malformed bedGraph line; needs >=4 fields: {raw!r}",
                )

            chromosome, start_text, end_text, value_text = parts[:4]

            try:
                start = int(start_text)
                end = int(end_text)
            except ValueError as error:
                flush(output_handle)

                raise ValueError(f"Non-integer start/end: {raw!r}") from error

            if end <= start:
                flush(output_handle)

                raise ValueError(f"Non-positive interval width: {raw!r}")

            compared = _comparison_value(
                value_text,
                decimal_places,
                eps,
            )
            comparison_key = compared.key

            same_chromosome = chromosome == last_chrom
            contiguous = start == last_end
            finite_pair = not compared.nonfinite and not last_is_nonfinite
            comparable_values = last_key is not None and (
                _comparison_keys_match(comparison_key, last_key, eps)
            )
            continues_current_run = (
                last_chrom is not None
                and same_chromosome
                and contiguous
                and finite_pair
                and comparable_values
            )

            if not continues_current_run:
                flush(output_handle)
                last_chrom = chromosome
                last_start = start
                last_end = end
                last_token = compared.source
                last_number = compared.numeric
                last_is_nonfinite = compared.nonfinite
                last_key = comparison_key
            else:
                last_end = end

        flush(output_handle)


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
        Parsed bedGraph input, output, and merge options.

    Raises
    ------
    SystemExit
        If argument parsing fails or help is requested.
    """

    parser = CapArgumentParser(
        description="Merge adjacent same-valued bedGraph bins (streaming).",
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
        "-fi",
        "--fil_in",
        dest="fil_in",
        required=True,
        help=(
            "Input file path. Path to the input bedGraph file (.gz is "
            "handled), or '-' for stdin.\n\n"
        ),
    )
    parser.add_argument(
        "-fo",
        "--fil_out",
        dest="fil_out",
        required=True,
        help=(
            "Output file path. Path to the output bedGraph file (.gz is "
            "handled), or '-' for stdout.\n\n"
        ),
    )
    parser.add_argument(
        "-dp",
        "--dp",
        dest="dp",
        type=int,
        default=None,
        help=(
            "Maximum number of decimal places retained for finite emitted "
            "values.\n"
            "\n"
            "Controls the value comparisons by N-decimal rounded strings. "
            "With '--eps', rounding is applied first, then |Δ| (absolute "
            "difference) <= eps is tested on the rounded numbers.\n\n"
        ),
    )
    parser.add_argument(
        "-e",
        "--eps",
        dest="eps",
        type=float,
        default=None,
        help=(
            "Zero tolerance epsilon for matching (merge if |Δ| <= eps).\n"
            "  - Alone: compare raw numeric values with |Δ| <= eps.\n"
            "  - With '--dp N': compare on rounded values (to N decimals) "
            "with |Δ| <= eps.\n\n"
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
    status : int
        Process exit code.

    Raises
    ------
    ValueError
        Propagated from 'merge_bins' when input data lines are malformed.
    """

    args = parse_args(argv)

    try:
        if args.fil_in != "-":
            check_exists(args.fil_in, "file", "bedGraph")

        if args.fil_out != "-":
            fil_out, _, _ = validate_output_path(
                args.fil_out,
                ("bedgraph", "bdg", "bg"),
            )
            check_writable(fil_out, "file")
        else:
            fil_out = "-"

    except ValueError as e:
        raise SystemExit(str(e)) from None

    skp_pfx = parse_skp_pfx(args.skp_pfx, default=DEF_SKP_PFX)

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
            decimal_places=args.dp,
            eps=args.eps,
            skp_pfx=skp_pfx,
        )
    except (
        ValueError,
        FileNotFoundError,
        PermissionError,
        OSError,
    ) as error:
        print(str(error), file=sys.stderr)

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
