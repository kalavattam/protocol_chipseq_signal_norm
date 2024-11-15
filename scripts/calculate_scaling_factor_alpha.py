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
#     input samples from ChIP-seq benchwork: mass, volume, sequencing depth,
#     and fragment length. The formula used is based on details provided in the
#     following publications:
#     - pubmed.ncbi.nlm.nih.gov/32994221
#     - pubmed.ncbi.nlm.nih.gov/37160995
#
#     See also the documentation in the following GitHub repository:
#     github.com/BradleyDickson/siQ-ChIP
#
# Usage:
#     python calculate_scaling_factor_alpha.py \
#         --mass_ip <float> \
#         --mass_in <float> \
#         --volume_ip <float> \
#         --volume_in <float> \
#         --depth_ip <int> \
#         --depth_in <int> \
#         --length_ip <float> \
#         --length_in <float>
#
# Arguments:
#     -mp, --mass_ip   (float): Mass of the IP sample.
#     -mn, --mass_in   (float): Mass of the input sample.
#     -vp, --volume_ip (float): Volume of the IP sample.
#     -vn, --volume_in (float): Volume of the input sample.
#     -dp, --depth_ip    (int): Sequencing depth of the IP sample.
#     -dn, --depth_in    (int): Sequencing depth of the input sample.
#     -lp, --length_ip (float): Mean fragment length of the IP sample.
#     -ln, --length_in (float): Mean fragment length of the input sample.
#
# Example:
#     python calculate_scaling_factor_alpha.py \
#         --mass_ip 10.5 \
#         --mass_in 8.0 \
#         --volume_ip 15.0 \
#         --volume_in 12.4 \
#         --depth_ip 5000000 \
#         --depth_in 4500000 \
#         --length_ip 200.0 \
#         --length_in 180.5
#
# Output:
#     The calculated alpha scaling factor as a positive floating point number. 
#
# License:
#     Distributed under terms of the MIT license.

import argparse
import sys


def calculate_alpha(
    mass_ip, mass_in,
    volume_ip, volume_in,
    depth_ip, depth_in,
    length_ip, length_in
):
    """
    Calculate a siQ-ChIP 'alpha' scaling factor using the provided values.
    
    Args:
        mass_ip   (float): Mass of IP sample.
        mass_in   (float): Mass of input sample.
        volume_ip (float): Volume of IP sample.
        volume_in (float): Volume of input sample.
        depth_ip    (int): Sequencing depth of IP sample.
        depth_in    (int): Sequencing depth of input sample.
        length_ip (float): Mean fragment length of IP sample.
        length_in (float): Mean fragment length of input sample.

    Returns:
        float: The calculated alpha scaling factor.
    """
    alpha = (
        (mass_ip / mass_in) *
        (volume_in / volume_ip) *
        (depth_in / depth_ip) *
        (length_in / length_ip)
    )
    return alpha


def main():
    """
    Execute the primary control flow for the script.

    main facilitates the calculation of the siQ-ChIP alpha scaling factor for
    ChIP-seq datasets by parsing command-line arguments for the required IP
    and input experimental values. Then, the function calculates the siQ-ChIP
    alpha scaling factor based on these inputs and prints the result.

    Command-line Arguments:
        -mp, --mass_ip   (float): Mass of the IP sample.
        -mn, --mass_in   (float): Mass of the input sample.
        -vp, --volume_ip (float): Volume of the IP sample.
        -vn, --volume_in (float): Volume of the input sample.
        -dp, --depth_ip    (int): Sequencing depth of the IP sample.
        -dn, --depth_in    (int): Sequencing depth of the input sample.
        -lp, --length_ip (float): Mean fragment length of the IP sample.
        -ln, --length_in (float): Mean fragment length of the input sample.
        -dp, --round   (int): Number of decimal places for rounding alpha.

    Returns:
        Outputs the siQ-ChIP alpha scaling factor.
    """
    # TODO Create a parse_args function
    parser = argparse.ArgumentParser(description=(
        'Calculate a siQ-ChIP alpha scaling factor for a ChIP-seq sample with '
        'IP and input data.'
    ))
    parser.add_argument(
        '-mp',
        '--mass_ip',
        type=float,
        required=True,
        help='Mass of IP sample'
    )
    parser.add_argument(
        '-mn',
        '--mass_in',
        type=float,
        required=True,
        help='Mass of input sample'
    )
    parser.add_argument(
        '-vp',
        '--volume_ip',
        type=float,
        required=True,
        help='Volume of IP sample.'
    )
    parser.add_argument(
        '-vn',
        '--volume_in',
        type=float,
        required=True,
        help='Volume of input sample.'
    )
    parser.add_argument(
        '-dp',
        '--depth_ip',
        type=int,
        required=True,
        help='Sequencing depth of IP sample.'
    )
    parser.add_argument(
        '-dn',
        '--depth_in',
        type=int,
        required=True,
        help='Sequencing depth of input sample.'
    )
    parser.add_argument(
        '-lp',
        '--length_ip',
        type=float,
        required=True,
        help='Mean fragment length of IP sample.'
    )
    parser.add_argument(
        '-ln',
        '--length_in',
        type=float,
        required=True,
        help='Mean fragment length of input sample.'
    )
    parser.add_argument(
        '-rp',
        '--round',
        type=int,
        default=6,
        required=False,
        help='Number of decimal places for rounding alpha.'
    )

    #  Display help and exit if no arguments were provided
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    #  Parse arguments
    args = parser.parse_args()

    #  Validate input values to ensure none are zero or negative
    if (
        args.mass_in == 0
        or args.volume_ip == 0
        or args.depth_ip == 0
        or args.length_ip == 0
    ):
        raise ValueError(
            "Input values cannot be zero to avoid division by zero."
        )

    if any(
        value <= 0 for value in [
            args.mass_ip, args.mass_in,
            args.volume_ip, args.volume_in,
            args.depth_ip, args.depth_in,
            args.length_ip, args.length_in
        ]
    ):
        raise ValueError("All input values must be positive.")
    
    #  Calculate the siQ-ChIP alpha scaling factor
    try:
        alpha = round(calculate_alpha(
            args.mass_ip, args.mass_in, args.volume_ip, args.volume_in,
            args.depth_ip, args.depth_in, args.length_ip, args.length_in
        ), args.round)
        print(f"{alpha}")
    except (ValueError, TypeError, ZeroDivisionError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
