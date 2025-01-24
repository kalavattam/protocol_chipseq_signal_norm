#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2024 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Script: calculate_scaling_factor_alpha.py
#
# Description:
#     This script calculates the siQ-ChIP alpha scaling factor for a ChIP-seq
#     sample. It requires the following experimental values for both IP and
#     input samples from ChIP-seq benchwork: mass, volume, sequencing depth (if
#     using '--eqn 5' or '--eqn 6'), and fragment length. The formulae used are
#     based on details provided in the following publications:
#         - pubmed.ncbi.nlm.nih.gov/32994221
#         - pubmed.ncbi.nlm.nih.gov/37160995
#
#     See also the documentation in the following GitHub repository:
#     github.com/BradleyDickson/siQ-ChIP.
#
# Usage:
#     python calculate_scaling_factor_alpha.py \
#         --eqn <str> \
#         --mass_ip <float> \
#         --mass_in <float> \
#         --vol_all <float> \
#         --vol_in <float> \
#         --dep_ip <int> \
#         --dep_in <int> \
#         --len_ip <float> \
#         --len_in <float> \
#         --rnd <int>
#
# Arguments:
#     -eq, --eqn       (str): Equation to compute. Options: '5', '5nd', '6', or
#                             '6nd'.
#     -mp, --mass_ip (float): Mass of the IP sample.
#     -mn, --mass_in (float): Mass of the input sample.
#     -va, --vol_all (float): Volume of sample before removal of input.
#     -vn, --vol_in  (float): Volume of the input sample.
#     -dp, --dep_ip    (int): Sequencing depth of the IP sample.
#     -dn, --dep_in    (int): Sequencing depth of the input sample.
#     -lp, --len_ip  (float): Mean fragment length of the IP sample.
#     -ln, --len_in  (float): Mean fragment length of the input sample.
#      -r, --rnd       (int): Number of decimal places for rounding alpha.
#
# Examples:
#     ```bash
#     python calculate_scaling_factor_alpha.py \
#         --eqn '6nd' \
#         --mass_ip 10.5 \
#         --mass_in 8.0 \
#         --vol_all 300 \
#         --vol_in 12.4 \
#         --len_ip 200.0 \
#         --len_in 180.5 \
#         --rnd 12
#     ```
#
#     ```bash
#     python calculate_scaling_factor_alpha.py \
#         --eqn '6' \
#         --mass_ip 10.5 \
#         --mass_in 8.0 \
#         --vol_all 300 \
#         --vol_in 12.4 \
#         --dep_ip 5000000 \
#         --dep_in 4500000 \
#         --len_ip 200.0 \
#         --len_in 180.5 \
#         --rnd 12
#     ```
#
# Output:
#     A positive floating-point number representing the alpha scaling factor.
#
# License:
#     Distributed under the MIT License.

import argparse
import sys


#  Run script in interactive/test mode (True) or command-line mode (False)
interactive = False


def set_interactive():
    """Set up paths and parameters for interactive mode."""
    #  Set values
    eqn = '5'
    mass_ip = 10.5
    mass_in = 8.0
    vol_all = 300
    vol_in = 12.4
    dep_ip = 5000000
    dep_in = 4500000
    len_ip = 200.0
    len_in = 180.5
    rnd = 24

    #  Return the arguments wrapped in argparse.Namespace
    return argparse.Namespace(
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


def validate_positive_values(args):
    """Ensure all provided values are positive and non-zero."""
    for name, value in vars(args).items():
        #  Skip non-numeric attributes (e.g., strings)
        if isinstance(value, (int, float)):
            if value <= 0:
                raise ValueError(
                    f"{name} must be greater than zero, but got {value}."
                )


def calculate_alpha(
    eqn, mass_ip, mass_in, vol_all, vol_in, dep_ip, dep_in, len_ip, len_in
):
    """
    Calculate a siQ-ChIP 'alpha' scaling factor using the provided values.
    
    This function computes the scaling factor (alpha) for siQ-ChIP experiments 
    based on the selected equation. It supports equations 5 and 6 (PMID:
    37160995) with or without sequencing depth terms (denoted as "nd" for "no
    depth").
    
    Args:
        eqn       (str): Alpha equation to compute. Options:
                         - '5':   Equation 5 (for use with raw coverage).
                         - '5nd': Equation 5 without depth terms (for use with
                                  normalized coverage).
                         - '6':   Equation 6 (for use with raw coverage).
                         - '6nd': Equation 6 without depth terms (for use with
                                  normalized coverage).
        mass_ip (float): Mass of the IP sample (e.g., immunoprecipitated DNA).
        mass_in (float): Mass of the input sample.
        vol_all (float): Volume of sample before removal of input.
        vol_in  (float): Volume of the input sample.
        dep_ip    (int): Sequencing depth of the IP sample.
        dep_in    (int): Sequencing depth of the input sample.
        len_ip  (float): Mean fragment length of the IP sample.
        len_in  (float): Mean fragment length of the input sample.

    Returns:
        float: The calculated alpha scaling factor.
    
    Raises:
        ValueError: If an unsupported equation is provided.
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
        #  vol_all > vol_in for equations '6' and '6nd'
        if vol_all <= vol_in:
            raise ValueError(
                f"For 'eqn={eqn}', 'vol_all' must be greater than 'vol_in'. "
                f"Received: 'vol_all={vol_all}', 'vol_in={vol_in}'."
            )
        if eqn == '6':
            #  Equation 6: Compute concentrations c_IP and c_in and use their
            #  ratio to calculate alpha
            c_ip = (
                mass_ip / (660 * len_ip * (vol_all - vol_in))
            ) * (1 / dep_ip)
            c_in = (
                mass_in / (660 * len_in * vol_in)
            ) * (1 / dep_in)
            alpha = c_ip / c_in
        elif eqn == '6nd':
            #  Equation 6 (no depth): Excludes sequencing depth terms from
            #  calculation; compute concentrations c_IP and c_in adjusting for
            #  sequencing depth
            c_ip = (mass_ip / (660 * len_ip * (vol_all - vol_in)))
            c_in = (mass_in / (660 * len_in * vol_in))
            alpha = c_ip / c_in
    else:
        #  Raise an error for unsupported equations.
        raise ValueError(f"Unsupported equation specified: '{eqn}'")

    return alpha


def parse_args():
    """
    Parse command line arguments.

    Args:
        -eq, --eqn       (str): Equation to compute. Options: '5', '5nd', '6',
                                or '6nd'.
        -mp, --mass_ip (float): Mass of the IP sample.
        -mn, --mass_in (float): Mass of the input sample.
        -va, --vol_all (float): Volume of sample before removal of input.
        -vn, --vol_in  (float): Volume of the input sample.
        -dp, --dep_ip    (int): Sequencing depth of the IP sample (required for
                                --eqn 5 or --eqn 6; otherwise ignored).
        -dn, --dep_in    (int): Sequencing depth of the input sample (required
                                for --eqn 5 or --eqn 6; otherwise ignored).
        -lp, --len_ip  (float): Mean fragment length of the IP sample.
        -ln, --len_in  (float): Mean fragment length of the input sample.
         -r, --rnd       (int): Number of decimal places for rounding alpha.
    """
    parser = argparse.ArgumentParser(description=(
        'Calculate a siQ-ChIP alpha scaling factor for a ChIP-seq sample with '
        'IP and input data.'
    ))
    parser.add_argument(
        '-eq', '--eqn',
        type=str,
        required=True,
        choices=['5', '5nd', '6', '6nd'],
        default='6nd',
        help=(
            "Equation to compute the alpha scaling factor (PMID: 37160995; "
            "default: %(default)s). Options: '5' applies Equation 5 for use "
            "with fragment length-normalized coverage, '5nd' uses Equation 5 "
            "without depth terms for use with 'normalized' coverage, '6' "
            "applies Equation 6 for use with fragment length-normalized "
            "coverage, and '6nd' uses Equation 6 without depth terms for use "
            "with 'normalized' coverage."
        )
    )
    parser.add_argument(
        '-mp', '--mass_ip',
        type=float,
        required=True,
        help='Mass of IP sample.'
    )
    parser.add_argument(
        '-mn', '--mass_in',
        type=float,
        required=True,
        help='Mass of input sample.'
    )
    parser.add_argument(
        '-va', '--vol_all',
        type=float,
        required=True,
        help='Volume of sample before removal of input.'
    )
    parser.add_argument(
        '-vn', '--vol_in',
        type=float,
        required=True,
        help='Volume of input sample.'
    )
    parser.add_argument(
        '-dp', '--dep_ip',
        type=int,
        required=False,
        help=(
            'Sequencing depth of IP sample (required for --eqn \'5\' or --eqn'
            '\'6\').'
        )
    )
    parser.add_argument(
        '-dn', '--dep_in',
        type=int,
        required=False,
        help=(
            'Sequencing depth of input sample (required for --eqn \'5\' or '
            '--eqn \'6\').'
        )
    )
    parser.add_argument(
        '-lp', '--len_ip',
        type=float,
        required=True,
        help='Mean fragment length of IP sample.'
    )
    parser.add_argument(
        '-ln', '--len_in',
        type=float,
        required=True,
        help='Mean fragment length of input sample.'
    )
    parser.add_argument(
        '-r', '--rnd',
        type=int,
        default=24,
        required=False,
        help='Number of decimal places for rounding alpha.'
    )

    #  Display help and exit if no arguments were provided
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    return parser.parse_args()


def main():
    """
    Execute the primary control flow for the script.

    main() facilitates the calculation of the siQ-ChIP alpha scaling factor for
    ChIP-seq datasets by parsing command line arguments for the required IP
    and input experimental values. Then, the function calculates the siQ-ChIP
    alpha scaling factor based on these inputs and prints the result.

    Args:
        None.

    Returns:
        None. Outputs the siQ-ChIP alpha scaling factor.

    Raises:
        - ValueError: If invalid input values (e.g., zero or negative) are
                      provided.
        - TypeError: If incorrect data types are passed to the function.
        - ZeroDivisionError: If a division by zero occurs during calculation.
    """
    #  Use command line arguments or interactive setup based on 'interactive'
    args = set_interactive() if interactive else parse_args()

    #  Validate depth arguments for equations requiring depth
    if args.eqn in {'5', '6'}:
        if args.dep_ip is None or args.dep_in is None:
            print(
                "Error: Equations '5' and '6' require sequencing depth "
                "values for IP ('--dep_ip') and input ('--dep_in') samples.",
                file=sys.stderr
            )
            sys.exit(1)

    #  Calculate the siQ-ChIP alpha scaling factor
    try:
        #  Validate input values to ensure none are zero or negative
        validate_positive_values(args)

        #  Compute alpha with the provided equation
        alpha = round(calculate_alpha(
            args.eqn,
            args.mass_ip, args.mass_in,
            args.vol_all, args.vol_in,
            args.dep_ip or 1, args.dep_in or 1,
            args.len_ip, args.len_in
        ), args.rnd)
        print(f"{alpha}")
    except (ValueError, TypeError, ZeroDivisionError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
