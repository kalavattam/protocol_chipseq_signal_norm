#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2024-2025 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Distributed under the MIT license.

"""
Script
------
calculate_scaling_factor_siq_chip.py


Description
------------
Calculates the siQ-ChIP alpha scaling factor for a ChIP-seq sample. Requires
the following experimental values for both IP and input samples from ChIP-seq
benchwork: mass, volume, sequencing depth (if using '--eqn 5' or '--eqn 6'),
and fragment length.

The formulas used are based on details provided in the following publications:
    - https://pubmed.ncbi.nlm.nih.gov/32994221
    - https://pubmed.ncbi.nlm.nih.gov/37160995

See also the documentation in the following GitHub repository:
    - https://github.com/BradleyDickson/siQ-ChIP


Usage
-----
python -m scripts.calculate_scaling_factor_siq_chip \
    [--help] [--verbose] \
    --eqn <str> \
    --mass_ip <float_ng> --mass_in <float_ng> \
    --vol_all <float_µL> --vol_in <float_µL> \
    [--dep_ip <int_fragments>] [--dep_in <int_fragments>] \
    --len_ip <float_bp> --len_in <float_bp> \
    [--rnd <int>]


Arguments
---------
 -v, --verbose
    Increase output verbosity.

-eq, --eqn
    Equation to compute: '5', '5nd', '6', or '6nd' (default: 6nd).

-mp, --mass_ip
    Mass of the IP sample (ng).

-mn, --mass_in
    Mass of the input sample (ng).

-va, --vol_all
    Volume of sample before removal of input (µL).

-vn, --vol_in
    Volume of the input sample (µL).

-di, --dep_ip
    Sequencing depth (in the form of, ideally, alignment-inferred fragments) of
    the IP sample. Required for '--eqn 5' or '--eqn 6'; otherwise, ignored.

-dn, --dep_in
    Sequencing depth (in the form of, ideally, alignment-inferred fragments) of
    the input sample. Required for '--eqn 5' or '--eqn 6'; otherwise, ignored.

-lp, --len_ip
    Summary fragment length of the IP sample (bp).

-ln, --len_in
    Summary fragment length of the input sample (bp).

-dp, --dp, --rnd, --round, --decimals, --digits
    Number of decimal places for rounding alpha (default: 24).


Output
------
- A positive floating-point number representing the alpha scaling factor.
- Non-zero exit on error:
    + invalid/negative inputs
    + missing '--dep_*' for '--eqn 5' or '--eqn 6'
    + 'vol_all <= vol_in' for '--eqn 6' or '--eqn 6nd'
    + unsupported equation


Examples
--------
1. Workflow to compute alpha for normalized coverage (uses equation '6nd',
   which omits depth terms)
'''bash
python -m scripts.calculate_scaling_factor_siq_chip \
    --eqn '6nd' \
    --mass_ip 10.5 \
    --mass_in 8.0 \
    --vol_all 300 \
    --vol_in 12.4 \
    --len_ip 200.0 \
    --len_in 180.5 \
    --rnd 12
'''

2. Workflow to compute alpha for fragment length-adjusted raw signal (uses
   equation '6', which includes depth terms)
'''bash
python -m scripts.calculate_scaling_factor_siq_chip \
    --eqn '6' \
    --mass_ip 10.5 \
    --mass_in 8.0 \
    --vol_all 300 \
    --vol_in 12.4 \
    --dep_ip 5000000 \
    --dep_in 4500000 \
    --len_ip 200.0 \
    --len_in 180.5 \
    --rnd 12
'''


General notes
-------------
- Mass in ng.
- Volume in µL.
- Fragment length in bp.
- Sequencing depth in alignments or, preferably, alignment-inferred fragments.
- The constant 660 is the approximate molecular weight per base pair
  (approximately 660 g × mol⁻¹ × bp⁻¹).
- Alpha is dimensionless. Use consistent units across IP and input (e.g., ng
  for mass and µL for volume); the “660” factor cancels in the ratio.


Performance notes
-----------------
- O(1) arithmetic.
- Runtime is dominated by argument parsing, arithmetic, and printing, as
  there’s no file I/O. Input sizes do not affect runtime.
- Floating-point is double precision.
- Formatting controls final string rounding.
- Given the same inputs, output is deterministic.
"""

from __future__ import annotations

import argparse
import signal
import sys

from contextlib import redirect_stdout

from scripts.functions.utils_check import check_cmp
from scripts.functions.utils_cli import (
    add_help_cap,
    CapArgumentParser
)
from scripts.functions.utils_interactive import (
    echo_block,
    get_args_interactive
)

try:
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
except (AttributeError, ValueError):
    pass

assert sys.version_info >= (3, 10), "Python >= 3.10 required."

#  Run script in “interactive mode” (true) or “CLI mode” (false)
interactive = False


def set_interactive(
    echo: bool = False, ensure: bool = False
) -> argparse.Namespace:
    """
    Set paths and parameters for interactive mode.

    Notes:
        This module does not create or check any filesystem paths, so 'ensure'
        does nothing here. Passing 'ensure=True' is allowed but has no effect
        beyond emitting a warning.
    """
    if ensure:
        print(
            "Warning: 'set_interactive(ensure=True)' has no effect in this "
            "module, as no paths are created or checked.",
            file=sys.stderr
        )

    #  Set argument values
    verbose = True
    eqn = '5'          # '5nd' '6' '6nd'
    mass_ip = 3.575
    mass_in = 44.72
    vol_all = 300
    vol_in = 20
    dep_ip = 8491709   # ...
    dep_in = 11543808  # ...
    len_ip = 192.419   # 380
    len_in = 169.744   # 348
    rnd = 24

    #  Wrap the arguments in argparse.Namespace
    ns = argparse.Namespace(
        verbose=verbose,
        eqn=eqn,
        mass_ip=mass_ip,
        mass_in=mass_in,
        vol_all=vol_all,
        vol_in=vol_in,
        dep_ip=dep_ip,
        dep_in=dep_in,
        len_ip=len_ip,
        len_in=len_in,
        rnd=rnd
    )

    #  “Echo” the interactive argument assignments (etc.) if specified
    if echo:
        echo_block("args", vars(ns))

    #  Return the argparse.Namespace-wrapped arguments
    return ns


def check_val_pos(args: argparse.Namespace) -> None:
    """
    Ensure required numeric values are positive; check depths only for
    '--eqn 5' or '--eqn 6'.

    Args:
        args : argparse.Namespace
            mass_ip, mass_in, vol_all, vol_in, len_ip, len_in, rnd, eqn, and
            optionally dep_ip, dep_in.

    Raises:
        ValueError
            If any required numeric argument (mass_ip, mass_in, vol_all,
            vol_in, len_ip, len_in) is <= 0.
        ValueError
            If 'vol_all <= vol_in'.
        ValueError
            If 'rnd < 0'.
        ValueError
            If eqn ∈ {'5','6'} and 'dep_ip' or 'dep_in' is missing or <= 0.
    """
    #  Check values that should always be positive
    for name in (
        'mass_ip', 'mass_in', 'vol_all', 'vol_in', 'len_ip', 'len_in'
    ):
        check_cmp(getattr(args, name), "gt", 0.0, name, allow_none=False)

    #  Check physical constraint for all equations, as input is taken from the
    #  whole
    if args.vol_all <= args.vol_in:
        raise ValueError(
            "Invalid volumes: 'vol_all' must be greater than 'vol_in', but "
            f"got 'vol_all={args.vol_all}' and 'vol_in={args.vol_in}'."
        )

    #  Check that '--rnd' >= 0
    check_cmp(args.rnd, "ge", 0, "rnd", allow_none=False)

    #  Check depth only for '--eqn 5' or '--eqn 6' (>0)
    if args.eqn in {"5", "6"}:
        dep_ip = getattr(args, "dep_ip", None)
        dep_in = getattr(args, "dep_in", None)

        if dep_ip is None or dep_in is None:
            raise ValueError(
                "Equations '5' and '6' include explicit sequencing-depth "
                "terms and require both '--dep_ip' and '--dep_in'. If signal "
                "tracks are already normalized by depth (e.g., as is the case "
                "for normalized coverage), use '5nd' or '6nd' instead, which "
                "omit depth terms."
            )

        #  Now that we know they are present, enforce positivity
        check_cmp(dep_ip, "gt", 0, "dep_ip", allow_none=False)
        check_cmp(dep_in, "gt", 0, "dep_in", allow_none=False)


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
) -> float:
    """
    Calculate a siQ-ChIP α (alpha) scaling factor using the provided values.

    This function computes the scaling factor α for siQ-ChIP experiments based
    on the selected equation. It supports equations 5 and 6 (PMID: 37160995)
    with or without sequencing depth terms (denoted as "nd" for "no depth").

    Args:
        eqn : str
            Alpha equation to compute. Options:
                - '5':   Equation 5 (for use with fragment length-adjusted raw
                         signal).
                - '5nd': Equation 5 without depth terms (for use with
                         normalized coverage).
                - '6':   Equation 6 (for use with fragment length-adjusted raw
                         signal).
                - '6nd': Equation 6 without depth terms (for use with
                         normalized coverage).
        mass_ip : float
            Mass of the IP sample (e.g., immunoprecipitated DNA; ng).
        mass_in : float
            Mass of the input sample (ng).
        vol_all : float
            Volume of sample before removal of input (µL).
        vol_in : float
            Volume of the input sample (µL).
        dep_ip : int
            Sequencing depth of the IP sample.
        dep_in : int
            Sequencing depth of the input sample.
        len_ip : float
            Summary fragment length of the IP sample (bp).
        len_in : float
            Summary fragment length of the input sample (bp).

    Returns:
        alpha : float
            The calculated alpha scaling factor.

    Raises:
        ValueError
            If an unsupported equation is provided, or if constraints required
            by the selected equation are violated (e.g., 'vol_all <= vol_in'
            for '--eqn 6' or '--eqn 6nd').
    """
    if eqn == '5':
        #  Equation 5: Alpha is proportional to mass ratios, volume ratios,
        #  depth ratios, and fragment length ratios
        alpha = (
            (mass_ip / mass_in) *
            (vol_in / vol_all) *
            (dep_in / dep_ip) *
            (len_in / len_ip)
        )
    elif eqn == '5nd':
        #  Equation 5 (no depth): Excludes sequencing depth terms from
        #  calculation
        alpha = (
            (mass_ip / mass_in) *
            (vol_in / vol_all) *
            (len_in / len_ip)
        )
    elif eqn in {'6', '6nd'}:
        #  To avoid division by 0 or a negative integer, check that
        #  'vol_all > vol_in' for equations '6' and '6nd'
        if vol_all <= vol_in:
            raise ValueError(
                f"For 'eqn={eqn}', 'vol_all' must be greater than 'vol_in'. "
                f"Received: 'vol_all={vol_all}' and 'vol_in={vol_in}'."
            )
        if eqn == '6':
            #  Equation 6: Compute concentrations 'c_IP' and 'c_in' and use
            #  their ratio to calculate alpha (to be used with ratios of
            #  fragment length-adjusted raw signal)
            c_ip = (
                mass_ip / (660 * len_ip * (vol_all - vol_in))
            ) * (1 / dep_ip)
            c_in = (
                mass_in / (660 * len_in * vol_in)
            ) * (1 / dep_in)
            alpha = c_ip / c_in
        elif eqn == '6nd':
            #  Equation 6 (no depth): Compute concentrations 'c_IP' and 'c_in'
            #  without depth terms (to be used with ratios of normalized
            #  coverage)
            c_ip = (mass_ip / (660 * len_ip * (vol_all - vol_in)))
            c_in = (mass_in / (660 * len_in * vol_in))
            alpha = c_ip / c_in
    else:
        #  Raise an error for unsupported equations.
        raise ValueError(f"Unsupported equation specified: '{eqn}'")

    return alpha


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse command line arguments.
    """
    parser = CapArgumentParser(description=(
        "Calculate a siQ-ChIP alpha scaling factor for a ChIP-seq sample with "
        "IP and input data."
    ))
    add_help_cap(parser)
    parser.add_argument(
        "-v", "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Increase output verbosity.\n\n"
    )
    parser.add_argument(
        "-eq", "--eqn", "--equation",
        dest="eqn",
        type=str,
        choices=["5", "5nd", "6", "6nd"],
        default="6nd",
        help=(
            "Equation to compute the siQ-ChIP alpha scaling factor (PMID: "
            "37160995; default: %(default)s).\n"
            "\n"
            "Equations '5' and '6' assume fragment length-adjusted raw "
            "signal, for which each fragment (not each read) contributes to "
            "the coverage signal. (Extending reads and including only first "
            "mates, or otherwise ensuring one fragment = one count, "
            "approximates that state.)\n"
            "\n"
            "Variants '5nd' and '6nd' omit the sequencing-depth term for "
            "workflows where signal tracks are already normalized by depth "
            "(e.g., “normalized coverage”).\n\n"
        )
    )
    parser.add_argument(
        "-mp", "--mass_ip",
        dest="mass_ip",
        type=float,
        required=True,
        help="Mass of IP sample (ng).\n\n"
    )
    parser.add_argument(
        "-mn", "--mass_in",
        dest="mass_in",
        type=float,
        required=True,
        help="Mass of input sample (ng).\n\n"
    )
    parser.add_argument(
        "-va", "--vol_all",
        dest="vol_all",
        type=float,
        required=True,
        help="Volume of sample before removal of input (µL).\n\n"
    )
    parser.add_argument(
        "-vn", "--vol_in",
        dest="vol_in",
        type=float,
        required=True,
        help="Volume of input sample (µL).\n\n"
    )
    parser.add_argument(
        "-di", "--dep_ip",
        dest="dep_ip",
        type=int,
        required=False,
        help=(
            "Sequencing depth of IP sample (alignments or alignment-inferred "
            "fragments; required for '--eqn 5' or '--eqn 6').\n\n"
        )
    )
    parser.add_argument(
        "-dn", "--dep_in",
        dest="dep_in",
        type=int,
        required=False,
        help=(
            "Sequencing depth of input sample (alignments or alignment-"
            "inferred fragments; required for '--eqn 5' or '--eqn 6').\n\n"
        )
    )
    parser.add_argument(
        "-lp", "--len_ip",
        dest="len_ip",
        type=float,
        required=True,
        help="Summary fragment length of IP sample (bp).\n\n"
    )
    parser.add_argument(
        "-ln", "--len_in",
        dest="len_in",
        type=float,
        required=True,
        help="Summary fragment length of input sample (bp).\n\n"
    )
    parser.add_argument(
        "-dp", "--dp", "--rnd", "--round", "--decimals", "--digits",
        dest="rnd",
        type=int,
        default=24,
        required=False,
        help=(
            "Number of decimal places for rounding alpha (default: "
            "%(default)s).\n\n"
        )
    )

    #  If no arguments are provided, use 'argv' to display help and exit
    argv_parse = sys.argv[1:] if argv is None else argv
    if not argv_parse:
        parser.print_help(sys.stderr)
        sys.exit(0)

    return parser.parse_args(argv_parse)


def main(argv: list[str] | None = None) -> None:
    """
    Execute the primary control flow for the script: parses command-line
    arguments, validates inputs, computes the siQ-ChIP alpha scaling factor,
    and prints the result to stdout.

    Args:
        None.

    Returns:
        None. Prints the siQ-ChIP alpha value to stdout on success.

    Side effects:
        May emit warnings (e.g., unusual volume ratios) to stderr via logging.
        Prints human-readable error messages to stderr on failure.

    Exits:
        0 on success or when showing help with no arguments, 1 on validation or
        computation errors (e.g., invalid/negative inputs, missing depths for
        eqn 5/6, vol_all <= vol_in for eqn 6/6nd, unsupported equation).
    """
    #  Handle interactive or CLI arguments
    args, early_exit = get_args_interactive(
        argv, interactive, set_interactive, parse_args
    )
    if early_exit:
        return 0

    #  Warn if input fraction looks unexpectedly/unusually large
    if args.eqn in {'5', '5nd'}:
        ratio = args.vol_in / args.vol_all
        if ratio > 0.5:
            print(
                f"Warning: vol_in/vol_all = {ratio:.3f} (>0.5), which is "
                "unusual but not necessarily invalid.",
                file=sys.stderr
            )

    if args.verbose:
        with redirect_stdout(sys.stderr):
            print("")
            print("#######################################################")
            print("## Arguments for 'calculate_scaling_factor_siq_chip' ##")
            print("#######################################################")
            print("")
            print("--verbose")
            print(f"--eqn     {args.eqn}")
            print(f"--mass_ip {args.mass_ip}")
            print(f"--mass_in {args.mass_in}")
            print(f"--vol_all {args.vol_all}")
            print(f"--vol_in  {args.vol_in}")

            #  Safely handle optional/possibly-absent attributes
            dep_ip = getattr(args, "dep_ip", None)
            dep_in = getattr(args, "dep_in", None)
            if dep_ip is not None or dep_in is not None:
                print(f"--dep_ip  {dep_ip}")
                print(f"--dep_in  {dep_in}")

            print(f"--len_ip  {args.len_ip}")
            print(f"--len_in  {args.len_in}")
            print(f"--rnd     {args.rnd}")
            print("")
            print("")

    #  Calculate the siQ-ChIP alpha scaling factor
    try:
        #  Validate input values to ensure none are zero or negative
        check_val_pos(args)

        #  Safely handle depth attributes
        dep_ip = getattr(args, "dep_ip", None)
        dep_in = getattr(args, "dep_in", None)

        #  Compute alpha with the provided equation
        alpha = calculate_alpha(
            args.eqn,
            args.mass_ip, args.mass_in,
            args.vol_all, args.vol_in,
            dep_ip, dep_in,
            args.len_ip, args.len_in
        )

        out = round(alpha, args.rnd)
        if out == 0.0:  # Avoid printing “-0.0”
            out = 0.0
        print(f"{out:.{args.rnd}f}")
    except (ValueError, TypeError, ZeroDivisionError) as e:
        raise SystemExit(str(e))


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
