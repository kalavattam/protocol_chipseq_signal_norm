#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: relativize_scaling_factors.py
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
#   GPT-5.6);
# - Anthropic Claude Code (Opus 5).
#
# Distributed under the MIT license.


"""
Relativize scaling-factor table columns to a reference row.

IP samples define the maximum reference value. The optional input-sample mode
includes input rows in the transformed output without changing that reference.
The CLI writes a transformed scaling-factor table to stdout or the requested
output file.

Notes
-----
This relative transformation is intended for arbitrary-unit scaling factors.
Applying it to physical-quantity siQ-ChIP values changes their interpretation.

Examples
--------
python -m protocol_chipseq_signal_norm.cli.relativize_scaling_factors \\
    --fil_in <file> [options]
"""

import argparse
import csv
import sys

from protocol_chipseq_signal_norm.utilities.utils_cli import (
    CapArgumentParser,
    add_help_cap,
)


def format_scaled(value: float, round_digits: int) -> str:
    """
    Format a rounded scaled value without non-informative trailing zeros.
    """

    text = f"{round(value, round_digits):.{round_digits}f}"
    return text.rstrip("0").rstrip(".") if "." in text else text


def load_tsv(
    file_path: str,
) -> tuple[list[str], list[dict[str, str]]]:
    """
    Load a TSV file into a header list and row dictionaries.
    """

    with open(file_path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames or []
        rows = [dict(row) for row in reader]

    if not header:
        raise ValueError(f"Input TSV has no header: '{file_path}'.")

    if len(set(header)) != len(header):
        raise ValueError(
            f"Input TSV has duplicate column names: '{file_path}'.",
        )

    return header, rows


def determine_scaling_column(header: list[str]) -> str:
    """
    Determine which supported scaling-factor column should be used.
    """

    for column in ("siq", "spike"):
        if column in header:
            return column

    raise ValueError(
        "No supported scaling-factor column is present in the file: "
        "'siq' or 'spike'.",
    )


def relativize(
    header: list[str],
    rows: list[dict[str, str]],
    scaling_column: str,
    include_input: bool,
    round_digits: int,
) -> tuple[list[str], list[dict[str, str]]]:
    """
    Calculate the relativized values.

    Input samples are scaled only when requested.

    When 'include_input' is true, 'in_*' samples use the largest 'IP_*' value;
    otherwise, they are excluded.

    Parameters
    ----------
    header : list[str]
        Input TSV column names.
    rows : list[dict[str, str]]
        Mutable row dictionaries keyed by the input column names.
    scaling_column : str
        Column containing the scaling value to relativize.
    include_input : bool
        Whether to relativize input samples as well as IP samples.
    round_digits : int
        Maximum decimal precision for rendered scaled values.

    Returns
    -------
    header, rows : tuple[list[str], list[dict[str, str]]]
        Updated header and rows containing the inserted 'scaled' column.

    Raises
    ------
    ValueError
        If required values are absent, nonnumeric, or cannot define a scale.
    """

    if "sample" not in header:
        raise ValueError("Input TSV is missing required column 'sample'.")

    ip_values = []
    has_input = False

    for row_number, row in enumerate(rows, start=2):
        sample = row.get("sample", "")

        if sample.startswith("in_"):
            has_input = True
            continue

        raw_value = row.get(scaling_column, "")

        try:
            ip_values.append(float(raw_value))
        except ValueError as error:
            raise ValueError(
                f"Non-numeric '{scaling_column}' value on row {row_number}: "
                f"'{raw_value}'.",
            ) from error

    if not ip_values:
        raise ValueError("Input TSV contains no IP samples to define scaling.")

    max_value = max(ip_values)
    if max_value == 0:
        raise ValueError(
            "Maximum IP scaling value is zero; cannot relativize.",
        )

    if include_input and not has_input:
        sys.stderr.write(
            "Warning: No 'in_*' samples found in the dataset, scaling only IP "
            "samples.\n",
        )

    for row_number, row in enumerate(rows, start=2):
        sample = row.get("sample", "")

        if sample.startswith("in_") and not include_input:
            row["scaled"] = "1"
            continue

        raw_value = row.get(scaling_column, "")

        try:
            scaled_value = float(raw_value) / max_value
        except ValueError as error:
            raise ValueError(
                f"Non-numeric '{scaling_column}' value on row {row_number}: "
                f"'{raw_value}'.",
            ) from error

        row["scaled"] = format_scaled(scaled_value, round_digits)

    scaling_index = header.index(scaling_column) + 1
    header = [*header[:scaling_index], "scaled", *header[scaling_index:]]

    return header, rows


def write_tsv(
    header: list[str],
    rows: list[dict[str, str]],
) -> None:
    """
    Write row dictionaries as TSV to stdout.
    """

    writer = csv.DictWriter(
        sys.stdout,
        fieldnames=header,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore",
    )
    writer.writeheader()
    writer.writerows(rows)


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
        Parsed table, sample-selection, and output options.

    Raises
    ------
    SystemExit
        If argument parsing fails or help is requested.
    """

    parser = CapArgumentParser(
        description=(
            "This script reads a TSV file containing ChIP-seq metrics and "
            "calculates a 'scaled' column by dividing each scaling-factor "
            "value by the maximum 'IP_*' value in the dataset. It allows "
            "optional inclusion of 'in_*' input samples in the scaling using "
            "the --input flag."
        ),
    )
    add_help_cap(parser)
    parser.add_argument(
        "-fi",
        "--fil_in",
        dest="fil_in",
        required=True,
        help=(
            "Input file path. Path to the input TSV file containing ChIP-seq "
            "metrics from running, e.g., execute_calculate_scaling_factor.sh. "
            "The file must contain a supported scaling-factor column: 'siq' "
            "or 'spike'."
        ),
    )
    parser.add_argument(
        "-in",
        "--input",
        dest="input",
        action="store_true",
        default=False,
        help=(
            "Include 'in_*' input samples in the relativization process. When "
            "this flag is set, both IP samples (e.g., 'IP_*') and input "
            "samples ('in_*') are scaled relative to the largest 'IP_*' "
            "sample value. The input samples are not used to determine the "
            "maximum value, but they are scaled using the same scaling factor "
            "derived from the IP samples. If the fil_in contains 'in_*' "
            "samples and '--input' is not specified, those samples will not "
            "be scaled."
        ),
    )
    parser.add_argument(
        "-dp",
        "--dp",
        dest="dp",
        type=int,
        default=6,
        help=(
            "Maximum number of decimal places retained for finite emitted "
            "values (default: %(default)s)."
        ),
    )

    argv_parse = sys.argv[1:] if argv is None else argv

    if not argv_parse:
        parser.print_help(sys.stderr)
        raise SystemExit(0)

    return parser.parse_args(argv_parse)


def main(argv: list[str] | None = None) -> int:
    """
    Run scaling-factor relativization and return an exit status.

    Parameters
    ----------
    argv : list[str] | None
        Explicit arguments, or None to read the process arguments.

    Returns
    -------
    status : int
        Zero on success and a stable nonzero status for anticipated failures.

    Raises
    ------
    SystemExit
        If command-line validation terminates before computation.
    """

    try:
        args = parse_args(argv)
        header, rows = load_tsv(args.fil_in)
        scaling_column = determine_scaling_column(header)

        header, rows = relativize(
            header,
            rows,
            scaling_column,
            args.input,
            args.dp,
        )

        write_tsv(header, rows)
    except (OSError, ValueError) as e:
        raise SystemExit(str(e)) from None

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
