#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: relativize_scaling_factors.py
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models) were used in
# development.
#
# Distributed under the MIT license.

# Description:
#     This script reads a TSV file containing ChIP-seq metrics and calculates
#     a 'scaled' column by dividing each scaling-factor value by the maximum
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
#     scale siQ-ChIP values in this way since they represent physical
#     quantities of chromatin, whereas spike-in values are arbitrary units.
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

import argparse
import csv
import sys


def format_scaled(value: float, round_digits: int) -> str:
    """Format a rounded scaled value without non-informative trailing zeros."""
    text = f"{round(value, round_digits):.{round_digits}f}"
    return text.rstrip("0").rstrip(".") if "." in text else text


def load_tsv(file_path):
    """Load a TSV file into a header list and row dictionaries."""
    with open(file_path, "rt", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        rows = [dict(row) for row in reader]

    if not header:
        raise ValueError(f"Input TSV has no header: '{file_path}'.")

    if len(set(header)) != len(header):
        raise ValueError(
            f"Input TSV has duplicate column names: '{file_path}'."
        )

    return header, rows


def determine_scaling_column(header):
    """Determine which supported scaling-factor column should be used."""
    for column in ("siq", "spike"):
        if column in header:
            return column

    raise ValueError(
        "No supported scaling-factor column is present in the file: "
        "'siq' or 'spike'."
    )


def relativize(header, rows, scaling_col, include_input, round_digits):
    """
    Calculate the relativized values. If include_input is True, scale 'in_*'
    samples by the largest 'IP_*' value; otherwise, exclude 'in_*' samples from
    the scaling.
    """
    if "sample" not in header:
        raise ValueError("Input TSV is missing required column 'sample'.")

    vals_ip = []
    has_input = False

    #  Find the maximum value among 'IP_*' samples (i.e., samples that do not
    #  start with 'in_')
    for idx, row in enumerate(rows, start=2):
        sample = row.get("sample", "")
        if sample.startswith("in_"):
            has_input = True
            continue

        raw = row.get(scaling_col, "")
        try:
            vals_ip.append(float(raw))
        except ValueError as e:
            raise ValueError(
                f"Non-numeric '{scaling_col}' value on row {idx}: '{raw}'."
            ) from e

    if not vals_ip:
        raise ValueError("Input TSV contains no IP samples to define scaling.")

    max_value = max(vals_ip)
    if max_value == 0:
        raise ValueError(
            "Maximum IP scaling value is zero; cannot relativize."
        )

    #  Check for condition in which --input flag is set and  there are no
    #  'in_*' samples in the TSV table
    if include_input and not has_input:
        sys.stderr.write(
            "Warning: No 'in_*' samples found in the dataset, "
            "scaling only IP samples.\n"
        )

    for idx, row in enumerate(rows, start=2):
        sample = row.get("sample", "")
        if sample.startswith("in_") and not include_input:
            row["scaled"] = "1"
            continue

        raw = row.get(scaling_col, "")
        try:
            scaled = float(raw) / max_value
        except ValueError as e:
            raise ValueError(
                f"Non-numeric '{scaling_col}' value on row {idx}: '{raw}'."
            ) from e

        row["scaled"] = format_scaled(scaled, round_digits)

    #  Dynamically reorder columns to place 'scaled' after the scaling factor
    scaling_index = header.index(scaling_col) + 1
    header = header[:scaling_index] + ["scaled"] + header[scaling_index:]

    return header, rows


def write_tsv(header, rows) -> None:
    """Write row dictionaries as TSV to stdout."""
    writer = csv.DictWriter(
        sys.stdout,
        fieldnames=header,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore",
    )
    writer.writeheader()
    writer.writerows(rows)


def parse_args():
    """
    Parse command line arguments.

    Args:
        ...
    """
    parser = argparse.ArgumentParser(
        description=(
            "This script reads a TSV file containing ChIP-seq metrics and "
            "calculates a 'scaled' column by dividing each scaling-factor "
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
            "running, e.g., execute_calculate_scaling_factor.sh. The file "
            "must contain a supported scaling-factor column: 'siq' or "
            "'spike'."
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

    #  Display help and exit if no arguments were provided
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    return parser.parse_args()


def main():
    try:
        #  Parse CLI arguments
        args = parse_args()

        #  Load the input TSV file
        header, rows = load_tsv(args.infile)

        #  Determine which scaling-factor column to use
        scaling_col = determine_scaling_column(header)

        #  Calculate scaled values
        header, rows = relativize(
            header, rows, scaling_col, args.input, args.round
        )

        #  Output the table with scaled values
        write_tsv(header, rows)
    except (OSError, ValueError) as e:
        raise SystemExit(str(e))


if __name__ == "__main__":
    main()
