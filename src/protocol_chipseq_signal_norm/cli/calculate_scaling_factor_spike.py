#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: calculate_scaling_factor_spike.py
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
# GPT-5.6) were used in design, development, and documentation, with all output
# reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Calculate spike-in-based ChIP-seq scaling coefficients.

The CLI accepts coefficient, output-format, count, and formatting options. It
prints one coefficient in plain format, or one or more coefficients as TSV or
JSON when requested.

Examples
--------
python -m protocol_chipseq_signal_norm.cli.calculate_scaling_factor_spike \\
    [options]
"""

from __future__ import annotations

import argparse
import json
import re
import signal
import sys
from contextlib import redirect_stdout, suppress

from protocol_chipseq_signal_norm.utilities.utils_check import (
    validate_comparison,
)
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
# Map user-facing '--coef' aliases to canonical coefficient names.
# Compatibility aliases remain accepted even when help omits them.
COEF_ALIAS_CANON = {
    "fractional": "fractional",
    "bioprotocol": "fractional",
    "bio_protocol": "fractional",
    "alavattam": "fractional",
    "tsukiyama": "fractional",
    "s": "fractional",
    "chiprx_alpha_ip": "chiprx_alpha_ip",
    "alpha_chiprx_ip": "chiprx_alpha_ip",
    "alpha_ip": "chiprx_alpha_ip",
    "chiprx_ip": "chiprx_alpha_ip",
    "orlando_ip": "chiprx_alpha_ip",
    "chiprx_alpha_in": "chiprx_alpha_in",
    "alpha_chiprx_in": "chiprx_alpha_in",
    "alpha_in": "chiprx_alpha_in",
    "chiprx_in": "chiprx_alpha_in",
    "orlando_in": "chiprx_alpha_in",
    "chiprx_alpha_ratio": "chiprx_alpha_ratio",
    "alpha_chiprx_ratio": "chiprx_alpha_ratio",
    "alpha_ratio": "chiprx_alpha_ratio",
    "chiprx_ratio": "chiprx_alpha_ratio",
    "orlando_ratio": "chiprx_alpha_ratio",
    "r": "chiprx_alpha_ratio",
    "rxinput_alpha": "rxinput_alpha",
    "alpha_rxinput": "rxinput_alpha",
    "rxi_alpha": "rxinput_alpha",
    "alpha_rxi": "rxinput_alpha",
    "rxinput": "rxinput_alpha",
    "rxi": "rxinput_alpha",
    "bressan": "rxinput_alpha",
    "ma": "rxinput_alpha",
    "niu": "rxinput_alpha",
    "fursova": "rxinput_alpha",
    "all": "all",
}

# Preserve the canonical output order when '--coef all' is selected.
COEF_ORDER = (
    "fractional",
    "chiprx_alpha_ip",
    "chiprx_alpha_in",
    "chiprx_alpha_ratio",
    "rxinput_alpha",
)


def normalize_coef(raw: str) -> str:
    """
    Normalize a user-supplied coefficient name or alias to its canonical name.
    """

    key = raw.strip().lower()

    # Treat hyphens like underscores and collapse repeated separators.
    key = re.sub(r"[-_]+", "_", key)

    if key in COEF_ALIAS_CANON:
        return COEF_ALIAS_CANON[key]

    message = (
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

    raise ValueError(message)


def round_value(value: float, decimal_places: int) -> float:
    """
    Round a value to the requested number of decimal places, avoiding -0.0.
    """

    rounded_value = round(value, decimal_places)

    if rounded_value == 0.0:
        return 0.0

    return rounded_value


def emit_output(
    factors: dict[str, float],
    coef: str,
    output_format: str,
    decimal_places: int,
) -> None:
    """
    Emit one or more computed coefficients in plain, TSV, or JSON format.

    Parameters
    ----------
    factors : dict[str, float]
        Computed coefficients keyed by canonical name.
    coef : str
        Requested coefficient name or 'all'.
    output_format : str
        Output format: 'plain', 'tsv', or 'json'.
    decimal_places : int
        Decimal precision for rendered values.

    Raises
    ------
    ValueError
        If the format or coefficient selection is invalid.
    """

    if output_format == "plain" and coef == "all":
        raise ValueError(
            "Format 'plain' requires a single coefficient; use "
            "'--format tsv|json' or a specific '--coef' value.",
        )

    if coef != "all" and coef not in factors:
        raise ValueError(
            f"Requested coefficient '{coef}' was not computed. "
            f"Computed keys: {', '.join(sorted(factors))}.",
        )

    if coef != "all":
        key = coef
        rendered_value = round_value(factors[key], decimal_places)

        if output_format == "plain":
            print(format_value(rendered_value, decimal_places))

            return

        if output_format == "tsv":
            print("coef\tvalue")
            print(f"{key}\t{format_value(rendered_value, decimal_places)}")

            return

        if output_format == "json":
            print(json.dumps({key: rendered_value}, sort_keys=False))

            return

    if output_format == "tsv":
        print("coef\tvalue")

        for key in COEF_ORDER:
            rendered_value = round_value(factors[key], decimal_places)
            print(f"{key}\t{format_value(rendered_value, decimal_places)}")

        return

    if output_format == "json":
        payload = {
            key: round_value(factors[key], decimal_places)
            for key in COEF_ORDER
        }
        print(json.dumps(payload, sort_keys=False))

        return

    raise ValueError(f"Unknown format '{output_format}'.")


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
        validate_comparison(
            getattr(args, name),
            "ge",
            0,
            name,
            allow_none=False,
        )

    validate_comparison(args.dp, "ge", 0, "dp", allow_none=False)

    args.coef_raw = args.coef

    args.coef = normalize_coef(args.coef)

    if args.fmt not in ("plain", "tsv", "json"):
        raise ValueError(f"Invalid --format '{args.fmt}'.")

    if args.fmt == "plain" and args.coef == "all":
        raise ValueError(
            "Format 'plain' requires a single coefficient; use "
            "'--format tsv|json' or a specific '--coef' value.",
        )

    if (args.main_ip + args.spike_ip) == 0:
        raise ZeroDivisionError(
            "IP totals are zero ('main_ip' + 'spike_ip' = 0).",
        )

    if (args.main_in + args.spike_in) == 0:
        raise ZeroDivisionError(
            "Input totals are zero ('main_in' + 'spike_in' = 0).",
        )


def calculate_scaling_factors(
    main_ip: int,
    spike_ip: int,
    main_in: int,
    spike_in: int,
    required: tuple[str, ...] = COEF_ORDER,
) -> dict[str, float]:
    """
    Compute spike-count algebra coefficients.

    The calculation uses main and spike counts for the IP and input samples.

    Parameters
    ----------
    main_ip, spike_ip, main_in, spike_in : int
        Non-negative integer counts.
    required : tuple[str, ...]
        Tuple of canonical coefficient names to compute. Canonical names:
        fractional, chiprx_alpha_ip, chiprx_alpha_in, chiprx_alpha_ratio,
        rxinput_alpha.

    Returns
    -------
    factors : dict[str, float]
        Keys are exactly those listed in 'required' (in any order).

    Raises
    ------
    TypeError
        If a count is not an integer.
    ValueError
        If a count is negative or a requested coefficient is unknown.
    ZeroDivisionError
        If a requested coefficient has a zero denominator.
    """

    for name, count in (
        ("main_ip", main_ip),
        ("spike_ip", spike_ip),
        ("main_in", main_in),
        ("spike_in", spike_in),
    ):
        if not isinstance(count, int):
            raise TypeError(
                f"Expected type 'int' for '{name}', but got "
                f"'{type(count).__name__}'.",
            )

        if count < 0:
            raise ValueError(
                f"Count for '{name}' must be >= 0, but got '{count}'.",
            )

    total_ip = main_ip + spike_ip
    total_input = main_in + spike_in

    if total_ip == 0:
        raise ZeroDivisionError(
            "IP totals are zero ('main_ip' + 'spike_ip' = 0).",
        )

    if total_input == 0:
        raise ZeroDivisionError(
            "Input totals are zero ('main_in' + 'spike_in' = 0).",
        )

    requested_names = set(required)
    unknown = requested_names - set(COEF_ORDER)

    if unknown:
        unknown_names = ", ".join(sorted(unknown))
        allowed_names = ", ".join(COEF_ORDER)

        raise ValueError(
            f"Invalid coefficient name(s) in 'required': {unknown_names}. "
            f"Allowed: {allowed_names}.",
        )

    spike_ip_coefficients = {
        "fractional",
        "chiprx_alpha_ip",
        "chiprx_alpha_ratio",
        "rxinput_alpha",
    }
    requires_spike_ip = bool(requested_names & spike_ip_coefficients)

    if requires_spike_ip and spike_ip == 0:
        raise ZeroDivisionError(
            "'spike_ip' is 0; cannot compute requested coefficient(s) that "
            "require division by N_s^IP.",
        )

    spike_input_coefficients = {
        "chiprx_alpha_in",
        "chiprx_alpha_ratio",
    }
    requires_spike_input = bool(requested_names & spike_input_coefficients)

    if requires_spike_input and spike_in == 0:
        raise ZeroDivisionError(
            "'spike_in' is 0; cannot compute requested coefficient(s) that "
            "require division by N_s^in.",
        )

    factors: dict[str, float] = {}

    # The bio-protocol coefficient follows Alavattam et al.
    if "fractional" in requested_names:
        factors["fractional"] = (spike_in / total_input) / (
            spike_ip / total_ip
        )

    # The ChIP-Rx IP coefficient follows Orlando et al.
    if "chiprx_alpha_ip" in requested_names:
        factors["chiprx_alpha_ip"] = 1e6 / spike_ip

    # The ChIP-Rx input coefficient follows Orlando et al.
    if "chiprx_alpha_in" in requested_names:
        factors["chiprx_alpha_in"] = 1e6 / spike_in

    if "chiprx_alpha_ratio" in requested_names:
        factors["chiprx_alpha_ratio"] = spike_in / spike_ip

    # The Rx-input coefficient follows Bressan et al.
    if "rxinput_alpha" in requested_names:
        factors["rxinput_alpha"] = (1e6 * spike_in) / (spike_ip * total_input)

    return factors


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse command line arguments.

    Parameters
    ----------
    argv : list[str] | None
        Explicit arguments, or None to read the process arguments.

    Returns
    -------
    arguments : argparse.Namespace
        Parsed spike-in counts, methods, and rendering options.

    Raises
    ------
    SystemExit
        If argument parsing fails or help is requested.
    """

    parser = CapArgumentParser(
        description=(
            "Calculate one or more spike-in-derived uniform scaling "
            "coefficients for a ChIP-seq IP/input pair."
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
        "-c",
        "--coef",
        "--coefficient",
        dest="coef",
        type=str,
        required=False,
        default="chiprx_alpha_ratio",
        help=(
            "Coefficient to output (case-insensitive; default: %(default)s)."
            "\n\n"
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
        ),
    )
    parser.add_argument(
        "-ft",
        "--fmt",
        "--format",
        dest="fmt",
        type=str,
        choices=("plain", "tsv", "json"),
        default="plain",
        help=(
            "Output format (default: %(default)s). 'plain' requires a single "
            "coefficient (i.e., not 'all').\n\n"
        ),
    )
    parser.add_argument(
        "-mp",
        "-mip",
        "--mip",
        "--main_ip",
        dest="main_ip",
        type=int,
        required=True,
        help=(
            "Number of “main” alignments (reads or read-inferred fragments) "
            "for the ChIP-seq IP data.\n\n"
        ),
    )
    parser.add_argument(
        "-sp",
        "-sip",
        "--sip",
        "--spike_ip",
        dest="spike_ip",
        type=int,
        required=True,
        help="Number of spike-in alignments for the ChIP-seq IP data.\n\n",
    )
    parser.add_argument(
        "-mn",
        "-min",
        "--min",
        "--main_in",
        dest="main_in",
        type=int,
        required=True,
        help=(
            "Number of “main” alignments for the corresponding ChIP-seq input "
            "data.\n\n"
        ),
    )
    parser.add_argument(
        "-sn",
        "-sin",
        "--sin",
        "--spike_in",
        dest="spike_in",
        type=int,
        required=True,
        help=(
            "Number of spike-in alignments for the corresponding ChIP-seq "
            "input data.\n\n"
        ),
    )
    parser.add_argument(
        "-dp",
        "--dp",
        dest="dp",
        type=int,
        default=24,
        required=False,
        help=(
            "Maximum number of decimal places retained for finite emitted "
            "values. Applies to emitted coefficient values. Plain and TSV "
            "output strip non-informative trailing zeros (default: "
            "%(default)s).\n\n"
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

    Parse arguments, validate counts, compute one or more spike-in scaling
    coefficients, and print the result to stdout.

    Parameters
    ----------
    argv : list[str] | None
        Optional command-line argument list. If None, arguments are read
        from 'sys.argv[1:]'.

    Returns
    -------
    status : int
        Zero on success after printing the requested coefficient output.

    Raises
    ------
    SystemExit
        For help, validation failures, or computation errors.

    Notes
    -----
    Ratio warnings and human-readable failure diagnostics are written to
    stderr.
    """

    args = parse_args(argv)

    try:
        validate_counts(args)

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

        ratio_ip = args.spike_ip / (args.main_ip + args.spike_ip)
        ratio_in = args.spike_in / (args.main_in + args.spike_in)

        # Additional composition warnings are non-fatal. These values are
        # mathematically valid but often indicate a counting or main/spike-in
        # partitioning problem.
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

        for label, frac in (("IP", ratio_ip), ("input", ratio_in)):
            if frac < 1e-6:
                print(
                    f"Warning: '{label}' spike-in fraction may be unusually "
                    f"low ({frac:.3e}).",
                    file=sys.stderr,
                )
            elif frac > 0.5:
                print(
                    f"Warning: '{label}' spike-in fraction may be unusually "
                    f"high ({frac:.3f}).",
                    file=sys.stderr,
                )

        required = COEF_ORDER if args.coef == "all" else (args.coef,)
        factors = calculate_scaling_factors(
            args.main_ip,
            args.spike_ip,
            args.main_in,
            args.spike_in,
            required=required,
        )

        if args.verbose:
            with redirect_stdout(sys.stderr):
                print("###############################################")
                print("## Computed spike fractions and coefficients ##")
                print("###############################################")
                print("")

                # Keep one shared value-column start for this block. Ratio
                # labels are “IP spike-in fraction” and “input spike-in
                # fraction”.
                base_labels = (
                    "IP spike-in fraction",
                    "input spike-in fraction",
                )
                coefficient_labels = (
                    COEF_ORDER if args.coef == "all" else (args.coef,)
                )

                all_labels = base_labels + tuple(coefficient_labels)
                label_width = max(len(label) for label in all_labels)

                def print_dotted(label: str, value: str) -> None:
                    dots = "." * (label_width - len(label) + 8)
                    print(f"{label} {dots} {value}")

                print_dotted(
                    "IP spike-in fraction",
                    format_value(ratio_ip, args.dp),
                )
                print_dotted(
                    "input spike-in fraction",
                    format_value(ratio_in, args.dp),
                )

                for key in coefficient_labels:
                    print_dotted(
                        key,
                        format_value(factors[key], args.dp),
                    )

                print("")
                print("")

        emit_output(factors, args.coef, args.fmt, args.dp)
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
