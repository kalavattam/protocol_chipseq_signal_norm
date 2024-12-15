#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2024 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Script: calculate_scaling_factor_spike.py
#
# Description:
#    This script calculates a spike-in-derived scaling factor for ChIP-seq
#    datasets. It requires tallies of "main" and spike-in alignments for
#    corresponding IP (immunoprecipitation) and input data. The script computes
#    the ratio of spike-in to main alignments for both IP and input data, and
#    derives the final scaling factor by dividing the ratio of spike-in to main
#    counts in the input sample by the ratio of spike-in to main counts in the
#    IP sample.
#
#    For more details, see the following URL: biostars.org/p/9572653/#9572655
#
# Usage:
#     python calculate_scaling_factor_spike.py \
#         --main_ip <int> \
#         --spike_ip <int> \
#         --main_in <int> \
#         --spike_in <int>
#
# Arguments:
#     -mp, --main_ip  (int): Tally of main alignments in IP sample.
#     -sp, --spike_ip (int): Tally of spike-in alignments in IP sample.
#     -mn, --main_in  (int): Tally of main alignments in input sample.
#     -sn, --spike_in (int): Tally of spike-in alignments in input sample.
#
# Output:
#    The scaling factor as a positive floating point number.
#
# Example:
#     python calculate_scaling_factor_spike.py \
#         --main_ip 100000 \
#         --spike_ip 5000 \
#         --main_in 90000 \
#         --spike_in 4500
#
# License:
#     Distributed under terms of the MIT license.

import argparse
import sys


#  Run script in interactive/test mode (True) or command-line mode (False)
interactive = False


def set_interactive():
    """Set parameters for interactive mode."""
    main_ip = ...  # TODO
    spike_ip = ...  # TODO
    main_in = ...  # TODO
    spike_in = ...  # TODO

    #  Return the arguments wrapped in argparse.Namespace
    return argparse.Namespace(
        main_ip=main_ip,
        spike_ip=spike_ip,
        main_in=main_in,
        spike_in=spike_in
    )


def calculate_scaling_factor(main_ip, spike_ip, main_in, spike_in):
    """
    Calculate a spike-in-derived scaling factor based on alignment counts.

    calculate_scaling_factor computes two ratios: the ratio of spike-in
    alignment counts to "main" alignment counts for IP data and the same for
    its corresponding input data. calculate_scaling_factor then calculates and
    returns the scaling factor by dividing the input ratio by the IP ratio.

    Args:
        main_ip  (int): The number of "main" alignments in the IP data.
        spike_ip (int): The number of spike-in alignments in the IP data.
        main_in  (int): The number of "main" alignments in the input data.
        spike_in (int): The number of spike-in alignments in the input data.

    Returns:
        float: The scaling factor calculated as the ratio of input sample
               ratio to IP sample ratio.
    
    Raises:
        ValueError: If any of the counts are negative.
        TypeError: If any of the counts are not integers.
        ZeroDivisionError: If any of the main alignment counts are zero.
    """
    #  Validate that all inputs are integers
    for name, count in zip(
        ['main_ip', 'spike_ip', 'main_in', 'spike_in'],
        [main_ip, spike_ip, main_in, spike_in]
    ):
        if not isinstance(count, int):
            raise TypeError(
                f'Expected integer for {name}, but got {type(count).__name__}.'
            )

    #  Validate that all inputs are non-negative
    for name, count in zip(
        ['main_ip', 'spike_ip', 'main_in', 'spike_in'],
        [main_ip, spike_ip, main_in, spike_in]
    ):
        if count < 0:
            raise ValueError(f'{name} alignment counts are negative.')

    #  Check for division by zero
    if main_ip == 0 or main_in == 0:
        raise ValueError(
            'Main alignment counts cannot be zero to avoid division by zero.'
        )

    #  Calculate ratios
    ratio_ip = spike_ip / main_ip
    ratio_in = spike_in / main_in

    # #  Check for very high or low ratios
    # if not (0.1 <= ratio_ip <= 10):
    #     print(
    #         f'Warning: The ratio of spike-in to "main" alignments for IP '
    #         f'data ({ratio_ip:.2f}) may be outside the expected range.'
    #     )
    #
    # if not (0.1 <= ratio_in <= 10):
    #     print(
    #         f'Warning: The ratio of spike-in to "main" alignments for input '
    #         f'data ({ratio_in:.2f}) may be outside the expected range.'
    #     )

    #  Calculate the scaling factor
    return ratio_in / ratio_ip


def parse_args():
    """
    Parse command line arguments.

    Args:
        -mp, --main_ip  (int): The number of "main" alignments for the ChIP-seq
                               IP sample.
        -sp, --spike_ip (int): The number of spike-in alignments for the
                               ChIP-seq IP sample.
        -mn, --main_in  (int): The number of "main" alignments for the
                               corresponding ChIP-seq input sample.
        -sn, --spike_in (int): The number of spike-in alignments for the
                               corresponding ChIP-seq input sample.
    """
    parser = argparse.ArgumentParser(description=(
        'Calculate a spike-in-derived scaling factor for a ChIP-seq sample '
        'with IP and input data.'
    ))
    parser.add_argument(
        '-mp',
        '--main_ip',
        type=int,
        required=True,
        help='Number of "main" alignments for the ChIP-seq IP data.'
    )
    parser.add_argument(
        '-sp',
        '--spike_ip',
        type=int,
        required=True,
        help='Number of spike-in alignments for the ChIP-seq IP data.'
    )
    parser.add_argument(
        '-mn',
        '--main_in',
        type=int,
        required=True,
        help=(
            'Number of "main" alignments for the corresponding ChIP-seq input '
            'data.'
        )
    )
    parser.add_argument(
        '-sn',
        '--spike_in',
        type=int,
        required=True,
        help=(
            'Number of spike-in alignments for the corresponding ChIP-seq '
            'input data.'
        )
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

    return parser.parse_args()


def main():
    """
    Execute the primary control flow for the script.

    main() facilitates the calculation of a spike-in-derived scaling factor for
    ChIP-seq samples. It parses command line arguments for the number of
    alignments from both IP and input samples, calculates ratios of spike-in to
    main alignments for both sample types, and then computes a final scaling
    factor by dividing the input sample ratio by the IP sample ratio.

    The function will terminate early and print an error message if any of the
    main alignment counts are zero, preventing division-by-zero errors.
    
    Args:
        ...

    Returns:
        Outputs the final scaling factor.
    """
    #  Use command line arguments or interactive setup based on `interactive`
    if interactive:
        args = set_interactive()
    else:
        args = parse_args()

    #  Calculate the scaling factor, handling exceptions as necessary
    try:
        result = round(calculate_scaling_factor(
            args.main_ip,
            args.spike_ip,
            args.main_in,
            args.spike_in
        ), args.round)
        print(f'{result}')
    except (ValueError, TypeError, ZeroDivisionError) as e:
        print(f'Error encountered during calculation: {e}', file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
