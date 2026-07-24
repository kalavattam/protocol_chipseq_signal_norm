#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: calculate_scaling_factor_spike.py
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
# GPT-5.6) were used in development and documentation.
#
# Distributed under the MIT license.


"""
Calculate spike-in-based ChIP-seq scaling coefficients.

Usage
-----
python -m protocol_chipseq_signal_norm.cli.calculate_scaling_factor_spike [options]

Parameters
----------
Coefficient, output-format, count, and formatting options are parsed by
'parse_args()'.

Returns
-------
Prints one coefficient in plain format, or one or more coefficients as TSV or
JSON when requested.

See Also
--------
docs/dev/scaling_factor_calculators.md
    Maintainer notes on supported coefficient families and counting
    assumptions.
"""

from __future__ import annotations

import argparse
import json
import re
import signal
import sys
from contextlib import redirect_stdout, suppress

from protocol_chipseq_signal_norm.utilities.utils_check import check_cmp
from protocol_chipseq_signal_norm.utilities.utils_cli import (
    CapArgumentParser,
    add_help_cap,
)
from protocol_chipseq_signal_norm.utilities.utils_format import format_value

with suppress(AttributeError, ValueError):
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

assert sys.version_info >= (3, 11), "Python >= 3.11 required."


# TODO: Add complete method references for the spike-in coefficient families.
# Revisit spike-in fraction warnings after empirical expected-range bounds are
# derived for the supported counting definitions.
#  Map user-facing '--coef' aliases to canonical coefficient names
COEF_ALIAS_CANON = {
    "fractional": "fractional",
    "bioprotocol": "fractional",
    "bio_protocol": "fractional",
    "alavattam": "fractional",                   # Not exposed in docs
    "tsukiyama": "fractional",                   # Not exposed in docs
    "s": "fractional",                           # For b/c; not exposed in docs

    "chiprx_alpha_ip": "chiprx_alpha_ip",
    "alpha_chiprx_ip": "chiprx_alpha_ip",
    "alpha_ip": "chiprx_alpha_ip",               # For b/c; not exposed in docs
    "chiprx_ip": "chiprx_alpha_ip",
    "orlando_ip": "chiprx_alpha_ip",             # Not exposed in docs

    "chiprx_alpha_in": "chiprx_alpha_in",
    "alpha_chiprx_in": "chiprx_alpha_in",
    "alpha_in": "chiprx_alpha_in",               # For b/c; not exposed in docs
    "chiprx_in": "chiprx_alpha_in",
    "orlando_in": "chiprx_alpha_in",             # Not exposed in docs

    "chiprx_alpha_ratio": "chiprx_alpha_ratio",
    "alpha_chiprx_ratio": "chiprx_alpha_ratio",
    "alpha_ratio": "chiprx_alpha_ratio",         # For consistency; not exposed
    "chiprx_ratio": "chiprx_alpha_ratio",
    "orlando_ratio": "chiprx_alpha_ratio",       # Not exposed in docs
    "r": "chiprx_alpha_ratio",                   # Not exposed in docs

    "rxinput_alpha": "rxinput_alpha",
    "alpha_rxinput": "rxinput_alpha",
    "rxi_alpha": "rxinput_alpha",
    "alpha_rxi": "rxinput_alpha",
    "rxinput": "rxinput_alpha",
    "rxi": "rxinput_alpha",
    "bressan": "rxinput_alpha",                  # Not exposed in docs
    "ma": "rxinput_alpha",                       # Not exposed in docs
    "niu": "rxinput_alpha",                      # Not exposed in docs
    "fursova": "rxinput_alpha",                  # Not exposed in docs

    "all": "all",
}

#  Canonical output order used when '--coef' all
COEF_ORDER = (
    "fractional", "chiprx_alpha_ip", "chiprx_alpha_in", "chiprx_alpha_ratio",
    "rxinput_alpha"
)


def normalize_coef(raw: str) -> str:
    """
    Normalize a user-supplied coefficient name or alias to its canonical name.
    """
    key = raw.strip().lower()

    #  Treat hyphens like underscores (and collapse runs)
    key = re.sub(r"[-_]+", "_", key)

    if key in COEF_ALIAS_CANON:
        return COEF_ALIAS_CANON[key]

    msg = (
        f"Invalid --coef '{raw}'.\n\n"
        "Accepted values and aliases (case-insensitive):\n"
        "    - fractional | bioprotocol | bio_protocol\n"
        "    - chiprx_alpha_ip | alpha_chiprx_ip | chiprx_ip\n"
        "    - chiprx_alpha_in | alpha_chiprx_in | chiprx_in\n"
        "    - chiprx_alpha_ratio | alpha_chiprx_ratio | chiprx_ratio\n"
        "    - rxinput_alpha | alpha_rxinput | rxi_alpha | alpha_rxi | "
        "rxinput | rxi\n"
        "    - all\n"
    )
    raise ValueError(msg)


def round_value(x: float, dp: int) -> float:
    """
    Round a value to the requested number of decimal places, avoiding -0.0.
    """
    y = round(x, dp)
    if y == 0.0:
        return 0.0
    return y


def emit_output(
    vals: dict[str, float],
    coef: str,
    fmt: str,
    dp: int
) -> None:
    """
    Emit one or more computed coefficients in plain, TSV, or JSON format.
    """
    if fmt == "plain" and coef == "all":
        raise ValueError(
            "Format 'plain' requires a single coefficient; use "
            "'--format tsv|json' or a specific '--coef' value."
        )

    if coef != "all" and coef not in vals:
        raise ValueError(
            f"Requested coefficient '{coef}' was not computed. "
            f"Computed keys: {', '.join(sorted(vals))}."
        )

    if coef != "all":
        key = coef
        v = round_value(vals[key], dp)
        if fmt == "plain":
            print(format_value(v, dp))
            return
        if fmt == "tsv":
            print("coef\tvalue")
            print(f"{key}\t{format_value(v, dp)}")
            return
        if fmt == "json":
            print(json.dumps({key: v}, sort_keys=False))
            return

    # coef == all
    if fmt == "tsv":
        print("coef\tvalue")
        for k in COEF_ORDER:
            v = round_value(vals[k], dp)
            print(f"{k}\t{format_value(v, dp)}")
        return

    if fmt == "json":
        obj = {k: round_value(vals[k], dp) for k in COEF_ORDER}
        print(json.dumps(obj, sort_keys=False))
        return

    raise ValueError(f"Unknown format '{fmt}'.")


def validate_counts(args: argparse.Namespace) -> None:
    """
    Ensure counts are valid and totals are non-zero.

    Parameters
    ----------
        args : argparse.Namespace
            main_ip, spike_ip, main_in, spike_in

    Raises
    ------
        ValueError
            If any count < 0 or if dp < 0.
        ZeroDivisionError
            If a per-sample total is zero.
    """
    for name in ("main_ip", "spike_ip", "main_in", "spike_in"):
        check_cmp(getattr(args, name), "ge", 0, name, allow_none=False)

    check_cmp(args.dp, "ge", 0, "dp", allow_none=False)

    #  Preserve raw user input before canonicalizing
    args.coef_raw = args.coef

    #  Canonicalize coef (aliases to canonical)
    args.coef = normalize_coef(args.coef)

    #  (Can’t be triggered, but keep it anyway)
    if args.fmt not in ("plain", "tsv", "json"):
        raise ValueError(f"Invalid --format '{args.fmt}'.")

    if args.fmt == "plain" and args.coef == "all":
        raise ValueError(
            "Format 'plain' requires a single coefficient; use "
            "'--format tsv|json' or a specific '--coef' value."
        )

    if (args.main_ip + args.spike_ip) == 0:
        raise ZeroDivisionError(
            "IP totals are zero ('main_ip' + 'spike_ip' = 0)."
        )
    if (args.main_in + args.spike_in) == 0:
        raise ZeroDivisionError(
            "Input totals are zero ('main_in' + 'spike_in' = 0)."
        )


def calculate_scaling_factors(
    main_ip: int, spike_ip: int,
    main_in: int, spike_in: int,
    required: tuple[str, ...] = COEF_ORDER
) -> dict[str, float]:
    """
    Compute spike-count algebra coefficients from main/spike counts in
    IP/input.

    Parameters
    ----------
        main_ip, spike_ip, main_in, spike_in
            Non-negative integer counts.
        required
            Tuple of canonical coefficient names to compute. Canonical names:
            fractional, chiprx_alpha_ip, chiprx_alpha_in, chiprx_alpha_ratio,
            rxinput_alpha.

    Returns
    -------
        dict[str, float]
            Keys are exactly those listed in 'required' (in any order).
    """
    #  Validate that all inputs are integers and are non-negative
    for name, count in (
        ("main_ip", main_ip), ("spike_ip", spike_ip),
        ("main_in", main_in), ("spike_in", spike_in)
    ):
        if not isinstance(count, int):
            raise TypeError(
                f"Expected type 'int' for '{name}', but got "
                f"'{type(count).__name__}'."
            )
        if count < 0:
            raise ValueError(
                f"Count for '{name}' must be >= 0, but got '{count}'."
            )

    t_ip = main_ip + spike_ip
    t_in = main_in + spike_in

    if t_ip == 0:
        raise ZeroDivisionError(
            "IP totals are zero ('main_ip' + 'spike_ip' = 0)."
        )
    if t_in == 0:
        raise ZeroDivisionError(
            "Input totals are zero ('main_in' + 'spike_in' = 0)."
        )

    #  Validate 'required' against the allowed canonical coefficient names
    req = set(required)
    unknown = req - set(COEF_ORDER)
    if unknown:
        raise ValueError(
            "Invalid coefficient name(s) in 'required': "
            f"{', '.join(sorted(unknown))}. Allowed: "
            f"{', '.join(COEF_ORDER)}."
        )

    #  Denominator checks only for requested coefficients
    if req & {
        "fractional", "chiprx_alpha_ip", "chiprx_alpha_ratio",
        "rxinput_alpha"
    } and spike_ip == 0:
        raise ZeroDivisionError(
            "'spike_ip' is 0; cannot compute requested coefficient(s) "
            "that require division by N_s^IP."
        )

    if req & {"chiprx_alpha_in", "chiprx_alpha_ratio"} and spike_in == 0:
        raise ZeroDivisionError(
            "'spike_in' is 0; cannot compute requested coefficient(s) "
            "that require division by N_s^in."
        )

    vals: dict[str, float] = {}

    #  Bio-protocol coefficient (Alavattam et al.)
    if "fractional" in req:
        #  fractional = (spike_in / t_in) / (spike_ip / t_ip)
        vals["fractional"] = (spike_in / t_in) / (spike_ip / t_ip)

    #  ChIP-Rx alpha for IP (Orlando et al.)
    if "chiprx_alpha_ip" in req:
        vals["chiprx_alpha_ip"] = 1e6 / spike_ip

    #  ChIP-Rx alpha for input (Orlando et al.)
    if "chiprx_alpha_in" in req:
        vals["chiprx_alpha_in"] = 1e6 / spike_in

    #  Ratio of ChIP-Rx alpha coefficients
    if "chiprx_alpha_ratio" in req:
        #  chiprx_alpha_ratio = (10^6 / spike_ip) / (10^6 / spike_in)
        #                     = spike_in / spike_ip
        vals["chiprx_alpha_ratio"] = spike_in / spike_ip

    #  Rx-input alpha (Bressan et al.)
    if "rxinput_alpha" in req:
        #  rxinput_alpha = (1e6 × spike_in) / (spike_ip × t_in)
        vals["rxinput_alpha"] = (1e6 * spike_in) / (spike_ip * t_in)

    return vals


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse command line arguments.
    """
    parser = CapArgumentParser(description=(
        "Calculate one or more spike-in-derived uniform scaling coefficients "
        "for a ChIP-seq IP/input pair."
    ))
    add_help_cap(parser)
    parser.add_argument(
        "-v", "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Run script in verbose mode.\n\n"
    )
    parser.add_argument(
        "-c", "--coef", "--coefficient",
        dest="coef",
        type=str,
        required=False,
        default="chiprx_alpha_ratio",
        help=(
            "Coefficient to output (case-insensitive; default: "
            "%(default)s).\n\n"
            "Aliases:\n"
            "    - fractional | bioprotocol | bio_protocol\n"
            "        (N_s^{in} / T^{in}) / (N_s^{IP} / T^{IP})\n"
            "    - chiprx_alpha_ip | alpha_chiprx_ip | chiprx_ip\n"
            "        10^6 / N_s^{IP}\n"
            "    - chiprx_alpha_in | alpha_chiprx_in | chiprx_in\n"
            "        10^6 / N_s^{in}\n"
            "    - chiprx_alpha_ratio | alpha_chiprx_ratio | chiprx_ratio\n"
            "        N_s^{in} / N_s^{IP}\n"
            "    - rxinput_alpha | alpha_rxinput | rxi_alpha |\n"
            "      alpha_rxi | rxinput | rxi\n"
            "        (10^6 * N_s^{in}) / (N_s^{IP} * T^{in})\n"
            "    - all\n"
            "\n"
            "Associated publications:\n"
            "    - fractional\tAlavattam et al.\n"
            "    - chiprx_alpha_ip\tOrlando et al.\n"
            "    - chiprx_alpha_in\tOrlando et al.\n"
            "    - rxinput_alpha\tNiu et al., Ma et al., Fursova et al., "
            "Bressan et al.\n\n"
        )
    )
    parser.add_argument(
        "-ft", "--fmt", "--format",
        dest="fmt",
        type=str,
        choices=("plain", "tsv", "json"),
        default="plain",
        help=(
            "Output format (default: %(default)s). 'plain' requires a single "
            "coefficient (i.e., not 'all').\n\n"
        )
    )
    parser.add_argument(
        "-mp", "-mip", "--mip", "--main_ip",
        dest="main_ip",
        type=int,
        required=True,
        help=(
            "Number of “main” alignments (reads or read-inferred fragments) "
            "for the ChIP-seq IP data.\n\n"
        )
    )
    parser.add_argument(
        "-sp", "-sip", "--sip", "--spike_ip",
        dest="spike_ip",
        type=int,
        required=True,
        help="Number of spike-in alignments for the ChIP-seq IP data.\n\n"
    )
    parser.add_argument(
        "-mn", "-min", "--min", "--main_in",
        dest="main_in",
        type=int,
        required=True,
        help=(
            "Number of “main” alignments for the corresponding ChIP-seq input "
            "data.\n\n"
        )
    )
    parser.add_argument(
        "-sn", "-sin", "--sin", "--spike_in",
        dest="spike_in",
        type=int,
        required=True,
        help=(
            "Number of spike-in alignments for the corresponding ChIP-seq "
            "input data.\n\n"
        )
    )
    parser.add_argument(
        "-dp", "--dp",
        dest="dp",
        type=int,
        default=24,
        required=False,
        help=(
            "Maximum number of decimal places retained for finite emitted values. Applies to emitted "
            "coefficient values. Plain and TSV output strip "
            "non-informative trailing zeros (default: %(default)s).\n\n"
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

    Parse arguments, validate counts, compute one or more spike-in scaling
    coefficients, and print the result to stdout.

    Parameters
    ----------
        argv: list[str] | None
            Optional command-line argument list. If None, arguments are read
            from 'sys.argv[1:]'.

    Returns
    -------
        int
            Exit status code. On success, prints the requested coefficient
            output and returns 0.

    Side effects:
        Emits warnings (e.g., unusually small or large IP or input ratios) to
        stderr. Prints human-readable error messages to stderr on failure.

    Exits:
        0 on success or when showing help with no arguments, 1 on validation or
        computation errors (e.g., negative counts, zero per-sample totals, or a
        requested division by zero).
    """
    #  Parse CLI arguments
    args = parse_args(argv)

    #  Compute requested coefficient(s), handling exceptions as necessary
    try:
        #  Check inputs
        validate_counts(args)

        #  Return “verbose banner” to stderr
        if args.verbose:
            with redirect_stdout(sys.stderr):
                print("####################################################")
                print("## Arguments for 'calculate_scaling_factor_spike' ##")
                print("####################################################")
                print("")
                print("--verbose")
                print(f"--coef     {args.coef_raw}  (canon: {args.coef})")
                print(f"--format   {args.fmt}")
                print(f"--main_ip  {args.main_ip}")
                print(f"--spike_ip {args.spike_ip}")
                print(f"--main_in  {args.main_in}")
                print(f"--spike_in {args.spike_in}")
                print(f"--dp      {args.dp}")
                print("")
                print("")

        #  Compute fractions for optional warnings and verbose echo
        ratio_ip = args.spike_ip / (args.main_ip + args.spike_ip)
        ratio_in = args.spike_in / (args.main_in + args.spike_in)

        #  Additional “implausible composition” warnings (non-fatal); these are
        #  mathematically valid but usually indicate a counting or main/spike-
        #  in partitioning problem
        for label, main, spike, frac in (
            ("IP", args.main_ip, args.spike_ip, ratio_ip),
            ("input", args.main_in, args.spike_in, ratio_in),
        ):
            if main == 0 and spike > 0:
                print(
                    f"Warning: '{label}' has 'main == 0' and 'spike > 0' "
                    f"(spike-in fraction = {frac:.3f}). This is unusual and "
                    "may indicate a counting issue, main/spike-in "
                    "partitioning problem, etc.",
                    file=sys.stderr,
                )
            if main > 0 and spike == 0:
                print(
                    f"Warning: '{label}' has 'spike == 0' and 'main > 0' "
                    "(no spike-in signal). This is unusual for spike-in "
                    "experiments and may indicate a counting issue, "
                    "main/spike-in partitioning problem, etc.",
                    file=sys.stderr,
                )

        #  Warn if spike-in fraction is vanishingly small or oddly large
        for label, frac in (("IP", ratio_ip), ("input", ratio_in)):
            if frac < 1e-6:
                print(
                    f"Warning: '{label}' spike-in fraction may be unusually "
                    f"low ({frac:.3e}).",
                    file=sys.stderr
                )
            elif frac > 0.5:
                print(
                    f"Warning: '{label}' spike-in fraction may be unusually "
                    f"high ({frac:.3f}).",
                    file=sys.stderr
                )

        required = COEF_ORDER if args.coef == "all" else (args.coef,)
        vals = calculate_scaling_factors(
            args.main_ip, args.spike_ip, args.main_in, args.spike_in,
            required=required
        )

        if args.verbose:
            with redirect_stdout(sys.stderr):
                print("###############################################")
                print("## Computed spike fractions and coefficients ##")
                print("###############################################")
                print("")

                #  Keep one shared “value column” start for everything in this
                #  block; ratio labels are “IP spike-in fraction” and ”input
                #  spike-in fraction”
                base_lbls = ("IP spike-in fraction", "input spike-in fraction")
                coef_lbls = COEF_ORDER if args.coef == "all" else (args.coef,)

                all_lbls = base_lbls + tuple(coef_lbls)
                width_lbl = max(len(s) for s in all_lbls)

                def print_dotted(lbl: str, value: str) -> None:
                    dots = "." * (width_lbl - len(lbl) + 8)
                    print(f"{lbl} {dots} {value}")

                print_dotted(
                    "IP spike-in fraction", format_value(ratio_ip, args.dp)
                )
                print_dotted(
                    "input spike-in fraction",
                    format_value(ratio_in, args.dp)
                )

                for k in coef_lbls:
                    print_dotted(k, format_value(vals[k], args.dp))

                print("")
                print("")

        emit_output(vals, args.coef, args.fmt, args.dp)
        return 0

    except (ValueError, TypeError, ZeroDivisionError) as e:
        raise SystemExit(str(e)) from None


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        with suppress(Exception):
            sys.stdout.close()
        with suppress(Exception):
            sys.stderr.close()
        raise SystemExit(0) from None
