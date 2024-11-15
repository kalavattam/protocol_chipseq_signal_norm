#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2024 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Script: relativize_scaling_factors.py
#
# Description:
#     This script reads a TSV file containing ChIP-seq metrics and calculates
#     a 'scaled' column by dividing each 'alpha' or 'sf' value by the maximum
#     'IP_*' value in the dataset, optionally including 'in_*' samples if the
#     --input flag is set. The scaled values represent a relative percentage,
#     with the largest 'IP_*' value set to 1. This script may be useful for
#     normalizing data across ChIP-seq samples for purposes of comparison.
#
#     By default, input (i.e., 'in_*') samples are excluded from the
#     normalization, and only IP (i.e., 'IP_*') samples are scaled by the 
#     maximum IP sample value. However, when the --input flag is used, input
#     samples are also scaled, but they are scaled by the largest value found
#     among the IP samples. For more details, see the following Biostars
#     comment: biostars.org/p/9572653/#9572962.
#
#     While we allow users to do so, we note that it is not appropriate to
#     scale 'alpha' values in this way since they represent physical quantities
#     of chromatin, whereas 'sf' values are arbitrary units ('scaled
#     coverage').
#
# Usage:
#     python relativize_scaling_factors.py [--input] --infile <input.tsv>
#
# Arguments:
#      -i, --infile  Input TSV file with ChIP-seq metrics.
#     -in, --input   Include 'in_*' samples in the scaling process. When this
#                    flag is used, the input samples are also scaled, but they
#                    are scaled by the largest 'IP_*' value, not the largest
#                    value overall.
#     -rp, --round   Set number of decimal places for rounding scaled values
#                    (default: 6).
#
# Example:
#     python relativize_scaling_factors.py \
#         --input \
#         --infile metrics.tsv \
#         --round 6
#
# Output:
#     Outputs the scaled table to stdout.
#
# License:
#     Distributed under terms of the MIT license.

#  Import libraries
import argparse
import pandas as pd
import sys


#  Run script in interactive/test mode (True) or command-line mode (False)
interactive = False


#  Define functions
def load_tsv(file_path):
    """Load a TSV file into a pandas DataFrame."""
    return pd.read_csv(file_path, sep='\t')


def determine_scaling_column(df):
    """Determine if 'alpha' or 'sf' should be used for scaling."""
    if 'alpha' in df.columns:
        return 'alpha'
    elif 'sf' in df.columns:
        return 'sf'
    else:
        raise ValueError(
            "Neither 'alpha' nor 'sf' columns are present in the file."
        )


def relativize(df, scaling_col, include_input, round_digits):
    """
    Calculate the relativized values. If include_input is True, scale 'in_*'
    samples by the largest 'IP_*' value; otherwise, exclude 'in_*' samples from
    the scaling.
    """
    #  Find the maximum value among 'IP_*' samples (i.e., samples that do not
    #  start with 'in_')
    ip_samples = df[~df['sample'].str.startswith('in_')]
    max_value = ip_samples[scaling_col].max()

    #  Check for condition in which --input flag is set and  there are no
    #  'in_*' samples in the TSV table
    if include_input and df['sample'].str.startswith('in_').sum() == 0:
        sys.stderr.write(
            "Warning: No 'in_*' samples found in the dataset, "
            "scaling only IP samples.\n"
        )

    if include_input:
        #  Scale all samples by the largest 'IP_*' value, including 'in_*'
        #  samples
        df['scaled'] = df.apply(
            lambda row: round(row[scaling_col] / max_value, round_digits),
            axis=1
        )
    else:
        #  Scale only non-input samples by the largest 'IP_*' value
        df['scaled'] = df.apply(
            lambda row: (
                round(row[scaling_col] / max_value, round_digits)
                if not row['sample'].startswith('in_')
                else 1
            ),
            axis=1
        )

    #  Dynamically reorder columns to place 'scaled' after 'alpha' or 'sf'
    cols = df.columns.tolist()
    scaling_index = cols.index(scaling_col) + 1
    reordered_cols = cols[:scaling_index] + ['scaled'] + cols[scaling_index:-1]
    df = df[reordered_cols]

    return df


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "This script reads a TSV file containing ChIP-seq metrics and "
            "calculates a 'scaled' column by dividing each 'alpha' or 'sf' "
            "value by the maximum 'IP_*' value in the dataset. It allows "
            "optional inclusion of 'in_*' input samples in the scaling using "
            "the --input flag."
        )
    )
    parser.add_argument(
        "-i",
        "--infile",
        required=True,
        help=(
            "Path to the input TSV file containing ChIP-seq metrics from "
            "running, e.g., execute_calculate_scaling_factor_spike.sh or "
            "execute_calculate_scaling_factor_alpha.sh. The file must "
            "contain either an 'alpha' or 'sf' column, and the script will "
            "use one of these to calculate scaled values."
        )
    )
    parser.add_argument(
        "-in",
        "--input",
        action="store_true",
        help=(
            "Include 'in_*' input samples in the relativization process. When "
            "this flag is set, both IP samples (e.g., 'IP_*') and input "
            "samples ('in_*') are scaled relative to the largest 'IP_*' "
            "sample value. The input samples are not used to determine the "
            "maximum value, but they are scaled using the same scaling factor "
            "derived from the IP samples. If the infile contains 'in_*' "
            "samples and '--input' is not specified, those samples will not "
            "be scaled."
        )
    )
    parser.add_argument(
        "-rp",
        "--round",
        type=int,
        default=6,
        help=(
            "Number of decimal places for rounding scaled values. The default "
            "is 6 decimal places."
        )
    )
    return parser.parse_args() if not interactive else argparse.Namespace(
        infile="",  # TODO
        input=False,
        round=6
    )


def main():
    #  Parse command-line arguments
    args = parse_args()

    #  Load the input TSV file
    df = load_tsv(args.infile)

    #  Determine whether to scale by 'alpha' or 'sf'
    scaling_col = determine_scaling_column(df)

    #  Calculate scaled values
    df = relativize(df, scaling_col, args.input, args.round)

    #  Output the DataFrame with scaled values
    print(df.to_csv(sep='\t', index=False))


if __name__ == "__main__":
    main()
