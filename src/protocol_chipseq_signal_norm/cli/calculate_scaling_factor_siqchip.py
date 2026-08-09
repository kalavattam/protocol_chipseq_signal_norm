#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: calculate_scaling_factor_siqchip.py
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
Calculate siQ-ChIP alpha scaling factors.

The CLI accepts equation choice, mass, volume, depth, fragment length, library
volume, and formatting options. It prints a positive floating-point alpha
scaling factor to stdout.

References
----------
- https://pubmed.ncbi.nlm.nih.gov/32994221
- https://pubmed.ncbi.nlm.nih.gov/37160995

Examples
--------
python -m protocol_chipseq_signal_norm.cli.calculate_scaling_factor_siqchip \\
    [options]
"""

from __future__ import annotations

import argparse
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


def check_val_pos(args: argparse.Namespace) -> None:
    """
    Ensure required numeric values are positive; check depths only for.

    '--eqn 5' or '--eqn 6'.

    Parameters
    ----------
    args : argparse.Namespace
        mass_ip, mass_in, vol_all, vol_in, len_ip, len_in, dp, eqn, and
        optionally dep_ip, dep_in.

    Raises
    ------
    ValueError
        If any required numeric argument (mass_ip, mass_in, vol_all, vol_in,
        len_ip, len_in) is <= 0.
    ValueError
        If 'vol_all <= vol_in'.
    ValueError
        If 'dp < 0'.
    ValueError
        If eqn ∈ {'5', '6'} and 'dep_ip' or 'dep_in' is missing or <= 0.
    """

    for name in (
        "mass_ip",
        "mass_in",
        "vol_all",
        "vol_in",
        "len_ip",
        "len_in",
    ):
        validate_comparison(
            getattr(args, name),
            "gt",
            0.0,
            name,
            allow_none=False,
        )

    # Input is taken from the whole, so the input fraction cannot exceed 1.
    if args.vol_all <= args.vol_in:
        raise ValueError(
            "Invalid volumes: 'vol_all' must be greater than 'vol_in', but "
            f"got 'vol_all={args.vol_all}' and 'vol_in={args.vol_in}'.",
        )

    validate_comparison(args.dp, "ge", 0, "dp", allow_none=False)

    # Library-loading volumes are optional only as a pair.
    #
    # Unequal normalized library loading is a benchwork reality (e.g., see
    # PRJNA857063), so reject one-sided metadata instead of silently treating
    # a missing side as 1.
    lib_vol_ip = getattr(args, "lib_vol_ip", None)
    lib_vol_in = getattr(args, "lib_vol_in", None)

    if (lib_vol_ip is None) != (lib_vol_in is None):
        raise ValueError(
            "Library-loading volume correction requires both '--lib_vol_ip' "
            "and '--lib_vol_in', or neither.",
        )

    if lib_vol_ip is not None and lib_vol_in is not None:
        validate_comparison(
            lib_vol_ip,
            "gt",
            0.0,
            "lib_vol_ip",
            allow_none=False,
        )
        validate_comparison(
            lib_vol_in,
            "gt",
            0.0,
            "lib_vol_in",
            allow_none=False,
        )

    if args.eqn in {"5", "6"}:
        dep_ip = getattr(args, "dep_ip", None)
        dep_in = getattr(args, "dep_in", None)

        if dep_ip is None or dep_in is None:
            raise ValueError(
                "Equations '5' and '6' include explicit sequencing-depth "
                "terms and require both '--dep_ip' and '--dep_in'. If signal "
                "tracks are already normalized by depth (e.g., as is the case "
                "for normalized coverage), use '5nd' or '6nd' instead, which "
                "omit depth terms.",
            )

        validate_comparison(dep_ip, "gt", 0, "dep_ip", allow_none=False)
        validate_comparison(dep_in, "gt", 0, "dep_in", allow_none=False)


def calculate_alpha(
    eqn: str,
    mass_ip: float,
    mass_in: float,
    vol_all: float,
    vol_in: float,
    dep_ip: int | None,
    dep_in: int | None,
    len_ip: float,
    len_in: float,
    lib_vol_ip: float | None = None,
    lib_vol_in: float | None = None,
) -> float:
    """
    Calculate a siQ-ChIP α (alpha) scaling factor using the provided values.

    This function computes the scaling factor α for siQ-ChIP experiments based
    on the selected equation. It supports equations 5 and 6 (PMID: 37160995)
    with or without sequencing depth terms (denoted as "nd" for "no depth").

    Parameters
    ----------
    eqn : str
        Alpha equation to compute. Options:
            - '5':   Equation 5 (for use with fragment length-adjusted raw
                     signal).
            - '5nd': Equation 5 without depth terms (for use with normalized
                     coverage).
            - '6':   Equation 6 (for use with fragment length-adjusted raw
                     signal).
            - '6nd': Equation 6 without depth terms (for use with normalized
                     coverage).
    mass_ip : float
        Mass of the IP sample (e.g., immunoprecipitated DNA; ng).
    mass_in : float
        Mass of the input sample (ng).
    vol_all : float
        Volume of sample before removal of input (µL).
    vol_in : float
        Volume of the input sample (µL).
    dep_ip : int | None
        Sequencing depth of the IP sample.
    dep_in : int | None
        Sequencing depth of the input sample.
    len_ip : float
        Summary fragment length of the IP sample (bp).
    len_in : float
        Summary fragment length of the input sample (bp).
    lib_vol_ip : float | None
        Volume of normalized IP library loaded into the sequencer (µL).
    lib_vol_in : float | None
        Volume of normalized input library loaded into the sequencer (µL).

    Returns
    -------
    alpha : float
        The calculated alpha scaling factor.

    Raises
    ------
    ValueError
        If an unsupported equation is provided, or if constraints required by
        the selected equation are violated (e.g., 'vol_all <= vol_in' for
        '--eqn 6' or '--eqn 6nd').
    """

    if eqn == "5":
        # Equation 5 uses mass, volume, depth, and fragment-length ratios.
        alpha = (
            (mass_ip / mass_in)
            * (vol_in / vol_all)
            * (dep_in / dep_ip)
            * (len_in / len_ip)
        )
    elif eqn == "5nd":
        # Equation 5 without depth excludes the sequencing-depth terms.
        alpha = (mass_ip / mass_in) * (vol_in / vol_all) * (len_in / len_ip)
    elif eqn in {"6", "6nd"}:
        # Equations 6 and 6nd require a positive IP volume.
        if vol_all <= vol_in:
            raise ValueError(
                f"For 'eqn={eqn}', 'vol_all' must be greater than 'vol_in'. "
                f"Received: 'vol_all={vol_all}' and 'vol_in={vol_in}'.",
            )

        if eqn == "6":
            # Equation 6 uses the concentration ratio with
            # fragment-length-adjusted raw-signal ratios.
            ip_concentration = (
                mass_ip / (660 * len_ip * (vol_all - vol_in))
            ) * (1 / dep_ip)
            input_concentration = (mass_in / (660 * len_in * vol_in)) * (
                1 / dep_in
            )
            alpha = ip_concentration / input_concentration
        elif eqn == "6nd":
            # Equation 6 without depth uses the concentration ratio with
            # normalized-coverage ratios.
            ip_concentration = mass_ip / (660 * len_ip * (vol_all - vol_in))
            input_concentration = mass_in / (660 * len_in * vol_in)
            alpha = ip_concentration / input_concentration
    else:
        # Raise an error for unsupported equations.
        raise ValueError(f"Unsupported equation specified: '{eqn}'")

    if lib_vol_ip is not None and lib_vol_in is not None:
        alpha *= lib_vol_in / lib_vol_ip

    return alpha


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
        Parsed siQ-ChIP input measurements and output options.

    Raises
    ------
    SystemExit
        If argument parsing fails or help is requested.
    """

    parser = CapArgumentParser(
        description=(
            "Calculate a siQ-ChIP alpha scaling factor for a ChIP-seq sample "
            "with IP and input data."
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
        "-eq",
        "--eqn",
        "--equation",
        dest="eqn",
        type=str,
        choices=["5", "5nd", "6", "6nd"],
        default="6nd",
        help=(
            "siQ-ChIP alpha equation to compute (PMID: 37160995; default: "
            "%(default)s).\n"
            "\n"
            "Equations '5' and '6' assume fragment length-adjusted raw "
            "signal, for which each fragment (not each read) contributes to "
            "the coverage signal. (Extending reads and including only first "
            "mates, or otherwise ensuring one fragment = one count, "
            "approximates that state).\n"
            "\n"
            "Variants '5nd' and '6nd' omit the sequencing-depth term for "
            "workflows where signal tracks are already normalized by depth "
            "(e.g., “normalized coverage”).\n\n"
        ),
    )

    parser.add_argument(
        "-mp",
        "--mass_ip",
        dest="mass_ip",
        type=float,
        required=True,
        help="Mass of IP sample (ng).\n\n",
    )
    parser.add_argument(
        "-mn",
        "--mass_in",
        dest="mass_in",
        type=float,
        required=True,
        help="Mass of input sample (ng).\n\n",
    )

    parser.add_argument(
        "-va",
        "--vol_all",
        dest="vol_all",
        type=float,
        required=True,
        help="Volume of sample before removal of input (µL).\n\n",
    )
    parser.add_argument(
        "-vn",
        "--vol_in",
        dest="vol_in",
        type=float,
        required=True,
        help="Volume of input sample (µL).\n\n",
    )

    parser.add_argument(
        "-di",
        "--dep_ip",
        dest="dep_ip",
        type=int,
        required=False,
        help=(
            "Sequencing depth of IP sample (alignments or alignment-inferred "
            "fragments; required for '--eqn 5' or '--eqn 6').\n\n"
        ),
    )
    parser.add_argument(
        "-dn",
        "--dep_in",
        dest="dep_in",
        type=int,
        required=False,
        help=(
            "Sequencing depth of input sample (alignments or alignment-"
            "inferred fragments; required for '--eqn 5' or '--eqn 6').\n\n"
        ),
    )

    parser.add_argument(
        "-lp",
        "--len_ip",
        dest="len_ip",
        type=float,
        required=True,
        help="Summary fragment length of IP sample (bp).\n\n",
    )
    parser.add_argument(
        "-ln",
        "--len_in",
        dest="len_in",
        type=float,
        required=True,
        help="Summary fragment length of input sample (bp).\n\n",
    )

    parser.add_argument(
        "-lvp",
        "--lib_vol_ip",
        "--lib-vol-ip",
        dest="lib_vol_ip",
        type=float,
        required=False,
        default=None,
        help=(
            "Volume of normalized IP library loaded into the sequencer (µL). "
            "Optional; if supplied, '--lib_vol_in' must also be supplied. "
            "When both library-loading volumes are supplied, alpha is "
            "multiplied by 'lib_vol_in / lib_vol_ip'.\n\n"
        ),
    )
    parser.add_argument(
        "-lvn",
        "--lib_vol_in",
        "--lib-vol-in",
        dest="lib_vol_in",
        type=float,
        required=False,
        default=None,
        help=(
            "Volume of normalized input library loaded into the sequencer "
            "(µL). Optional; if supplied, '--lib_vol_ip' must also be "
            "supplied. When both library-loading volumes are supplied, alpha "
            "is multiplied by 'lib_vol_in / lib_vol_ip'.\n\n"
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
            "values. Applies to alpha. Output strips non-informative trailing "
            "zeros (default: %(default)s).\n\n"
        ),
    )

    argv_parse = sys.argv[1:] if argv is None else argv

    if not argv_parse:
        parser.print_help(sys.stderr)
        raise SystemExit(0)

    return parser.parse_args(argv_parse)


def main(argv: list[str] | None = None) -> int:
    """
    Calculate and print the siQ-ChIP alpha scaling factor.

    The control flow parses arguments, validates inputs, and performs the
    selected equation.

    Parameters
    ----------
    argv : list[str] | None
        Arguments to parse. The process arguments are used by default.

    Returns
    -------
    status : int
        Zero on success after printing the alpha value.

    Raises
    ------
    SystemExit
        For help, validation failures, or computation errors.

    Notes
    -----
    Warnings and human-readable failure diagnostics are written to stderr.
    """

    args = parse_args(argv)

    if args.eqn in {"5", "5nd"}:
        ratio = args.vol_in / args.vol_all

        if ratio > 0.5:
            print(
                f"Warning: vol_in/vol_all = {ratio:.3f} (>0.5), which is "
                "unusual but not necessarily invalid.",
                file=sys.stderr,
            )

    if args.verbose:
        with redirect_stdout(sys.stderr):
            print("")
            print("#######################################################")
            print("## Arguments for 'calculate_scaling_factor_siqchip' ##")
            print("#######################################################")
            print("")
            print("--verbose")
            print(f"--eqn     {args.eqn}")
            print(f"--mass_ip {args.mass_ip}")
            print(f"--mass_in {args.mass_in}")
            print(f"--vol_all {args.vol_all}")
            print(f"--vol_in  {args.vol_in}")

            dep_ip = getattr(args, "dep_ip", None)
            dep_in = getattr(args, "dep_in", None)

            if dep_ip is not None or dep_in is not None:
                print(f"--dep_ip  {dep_ip}")
                print(f"--dep_in  {dep_in}")

            print(f"--len_ip  {args.len_ip}")
            print(f"--len_in  {args.len_in}")
            print(f"--lib_vol_ip {args.lib_vol_ip}")
            print(f"--lib_vol_in {args.lib_vol_in}")
            print(f"--dp     {args.dp}")
            print("")
            print("")

    try:
        check_val_pos(args)

        dep_ip = getattr(args, "dep_ip", None)
        dep_in = getattr(args, "dep_in", None)

        alpha = calculate_alpha(
            args.eqn,
            args.mass_ip,
            args.mass_in,
            args.vol_all,
            args.vol_in,
            dep_ip,
            dep_in,
            args.len_ip,
            args.len_in,
            args.lib_vol_ip,
            args.lib_vol_in,
        )

        print(format_value(alpha, args.dp))
    except (ValueError, TypeError, ZeroDivisionError) as e:
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
