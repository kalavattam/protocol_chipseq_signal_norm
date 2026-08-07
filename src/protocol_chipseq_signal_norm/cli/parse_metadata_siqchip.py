#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: parse_metadata_siqchip.py
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
Parse one alignment filename against siQ-ChIP metadata.

The CLI prints configured metadata-derived calculator inputs for the uniquely
matched metadata row.

See Also
--------
docs/design/parse_metadata_siqchip.md
    Maintainer notes on filename parsing, YAML configuration, and metadata-row
    matching.
"""

from __future__ import annotations

import argparse
import csv
import io
import os
import re
import shlex
import sys
from contextlib import redirect_stdout

from protocol_chipseq_signal_norm.utilities.utils_check import check_exists
from protocol_chipseq_signal_norm.utilities.utils_cli import (
    CapArgumentParser,
    add_help_cap,
)
from protocol_chipseq_signal_norm.utilities.utils_io import (
    is_header,
    open_in,
    parse_skp_pfx,
)

try:
    import yaml
except ImportError as e:
    raise ImportError(
        "PyYAML is required for parse_metadata_siqchip.py. Install it in the "
        "project environment, e.g., with 'mamba install pyyaml'.",
    ) from e


assert sys.version_info >= (3, 11), "Python >= 3.11 required."


IDENTIFIER_PATTERN = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
PAIRED_OPTIONS = (
    ("dep_ip", "dep_in"),
    ("len_ip", "len_in"),
    ("lib_vol_ip", "lib_vol_in"),
)


def validate_id(name: str, label: str) -> None:
    """
    Validate that a configured name can be exported as a shell identifier.
    """

    if not isinstance(name, str) or not IDENTIFIER_PATTERN.fullmatch(name):
        raise ValueError(f"'{label}' must be a shell-safe identifier.")


def require_mapping(obj: object, label: str) -> dict[object, object]:
    """
    Return a configuration section after confirming it is a mapping.
    """

    if not isinstance(obj, dict):
        raise ValueError(f"Configuration section '{label}' must be a mapping.")

    return obj


def require_list_string(
    obj: object,
    label: str,
    allow_empty: bool = False,
) -> list[str]:
    """
    Return a YAML list after confirming every item is a non-empty string.

    Parameters
    ----------
    obj : object
        YAML-derived value expected to contain a list.
    label : str
        Configuration key used in validation diagnostics.
    allow_empty : bool
        Whether an empty list satisfies this field's contract.

    Returns
    -------
    values : list[str]
        Validated string values in configuration order.

    Raises
    ------
    ValueError
        If the value is not a permitted list of nonempty strings.
    """

    if not isinstance(obj, list):
        raise ValueError(f"Configuration key '{label}' must be a YAML list.")

    if not allow_empty and not obj:
        raise ValueError(
            f"Configuration key '{label}' must be a non-empty list.",
        )

    values: list[str] = []

    for index, value in enumerate(obj, start=1):
        if not isinstance(value, str) or not value:
            item_description = f"Configuration key '{label}' item {index}"

            raise ValueError(
                f"{item_description} must be a non-empty string.",
            )

        values.append(value)

    return values


def normalize_input_map(spec: object, label: str) -> dict[str, str]:
    """
    Normalize calculator input declarations to output-name -> metadata-column.

    Lists mean identity mappings. Mappings allow explicit exported names.

    Parameters
    ----------
    spec : object
        YAML-derived list, mapping, or null input declaration.
    label : str
        Configuration key used in validation diagnostics.

    Returns
    -------
    inputs : dict[str, str]
        Exported input names mapped to metadata column names.

    Raises
    ------
    ValueError
        If the declaration shape, identifier, or column name is invalid.
    """

    inputs: dict[str, str] = {}

    if spec is None:
        return inputs

    if isinstance(spec, list):
        for key in require_list_string(spec, label, allow_empty=True):
            validate_id(key, f"{label}.{key}")

            inputs[key] = key

        return inputs

    if isinstance(spec, dict):
        for key, column in spec.items():
            validate_id(str(key), f"{label}.{key}")

            if not isinstance(column, str) or not column:
                raise ValueError(
                    f"Configuration key '{label}.{key}' must name one "
                    "metadata column.",
                )

            inputs[str(key)] = column

        return inputs

    raise ValueError(
        f"Configuration key '{label}' must be a YAML list or mapping.",
    )


def load_config(path: str, verbose: bool = False) -> dict:
    """
    Load and validate the siQ-ChIP metadata parser configuration.

    Parameters
    ----------
    path : str
        YAML configuration path.
    verbose : bool
        Whether to report the loaded configuration.

    Returns
    -------
    configuration : dict
        Normalized and validated configuration.

    Raises
    ------
    ValueError
        If required sections, fields, mappings, or identifiers are invalid.
    """

    with open(path, encoding="utf-8") as handle:
        configuration = yaml.safe_load(handle) or {}

    if not isinstance(configuration, dict):
        raise ValueError(
            f"Configuration file '{path}' must contain a YAML mapping.",
        )

    for section in ("filename", "matching", "calculator_inputs"):
        if section not in configuration:
            raise ValueError(
                "Configuration file is missing required top-level section "
                f"'{section}'.",
            )

    filename_config = require_mapping(
        configuration["filename"],
        "filename",
    )
    matching_config = require_mapping(
        configuration["matching"],
        "matching",
    )
    calculator_config = require_mapping(
        configuration["calculator_inputs"],
        "calculator_inputs",
    )
    siq_config = require_mapping(
        calculator_config.get("siqchip"),
        "calculator_inputs.siqchip",
    )

    delimiter = filename_config.get("delimiter", "_")

    if not isinstance(delimiter, str) or delimiter == "":
        raise ValueError(
            "Configuration key 'filename.delimiter' must be a string.",
        )

    filename_config["delimiter"] = delimiter

    filename_config["strip_extensions"] = require_list_string(
        filename_config.get("strip_extensions", []),
        "filename.strip_extensions",
        allow_empty=True,
    )
    filename_config["strip_suffixes"] = require_list_string(
        filename_config.get("strip_suffixes", []),
        "filename.strip_suffixes",
        allow_empty=True,
    )
    filename_fields = require_list_string(
        filename_config.get("fields"),
        "filename.fields",
    )
    seen: set[str] = set()

    for field in filename_fields:
        validate_id(field, f"filename.fields.{field}")

        if field in seen:
            raise ValueError(
                f"Configuration key 'filename.fields' repeats '{field}'.",
            )

        seen.add(field)

    filename_config["fields"] = filename_fields

    matching_fields = require_list_string(
        matching_config.get("fields"),
        "matching.fields",
    )

    for field in matching_fields:
        validate_id(field, f"matching.fields.{field}")

        if field not in filename_fields:
            raise ValueError(
                f"Configuration key 'matching.fields' references '{field}', "
                "which is not declared in 'filename.fields'.",
            )

    matching_config["fields"] = matching_fields

    field_columns = configuration.get("field_to_column", {})

    if field_columns is None:
        field_columns = {}

    field_columns = require_mapping(field_columns, "field_to_column")

    for field, column in field_columns.items():
        validate_id(str(field), f"field_to_column.{field}")

        if field not in filename_fields:
            raise ValueError(
                f"Configuration key 'field_to_column.{field}' references a "
                "field not declared in 'filename.fields'.",
            )

        if not isinstance(column, str) or not column:
            raise ValueError(
                f"Configuration key 'field_to_column.{field}' must name one "
                "column.",
            )

    configuration["field_to_column"] = field_columns

    required_inputs = normalize_input_map(
        siq_config.get("required", {}),
        "calculator_inputs.siqchip.required",
    )
    optional_inputs = normalize_input_map(
        siq_config.get("optional", {}),
        "calculator_inputs.siqchip.optional",
    )
    overlap = sorted(set(required_inputs) & set(optional_inputs))

    if overlap:
        raise ValueError(
            "Calculator input(s) listed as both required and optional: "
            + ", ".join(overlap),
        )

    siq_config["required"] = required_inputs
    siq_config["optional"] = optional_inputs

    table_config = configuration.get("table", {})

    if table_config is None:
        table_config = {}

    table_config = require_mapping(table_config, "table")
    table_config["skip_prefixes"] = require_list_string(
        table_config.get("skip_prefixes", ["#"]),
        "table.skip_prefixes",
        allow_empty=True,
    )
    configuration["table"] = table_config

    if verbose:
        with redirect_stdout(sys.stderr):
            print(f"'load_cfg' path: {path}")
            print(f"'load_cfg' cfg: {configuration}")

    return configuration


def strip_filename(path: str, configuration: dict) -> str:
    """
    Remove configured extensions and suffixes from an alignment filename.
    """

    stem = os.path.basename(path)

    for extension in configuration["filename"]["strip_extensions"]:
        if stem.endswith(extension):
            stem = stem[: -len(extension)]
            break

    for suffix in configuration["filename"]["strip_suffixes"]:
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break

    return stem


def parse_filename(
    path: str,
    configuration: dict,
    verbose: bool = False,
) -> dict[str, str]:
    """
    Parse an alignment filename into configured metadata fields.
    """

    stem = strip_filename(path, configuration)
    delimiter = configuration["filename"]["delimiter"]
    filename_fields = configuration["filename"]["fields"]
    tokens = stem.split(delimiter)

    if len(tokens) != len(filename_fields):
        raise ValueError(
            f"Filename stem '{stem}' split into {len(tokens)} token(s), but "
            f"configuration declares {len(filename_fields)} field(s): "
            f"{', '.join(filename_fields)}.",
        )

    parsed = dict(zip(filename_fields, tokens, strict=True))

    if verbose:
        with redirect_stdout(sys.stderr):
            print(f"'parse_filename' stem={stem}")
            print(f"'parse_filename' parsed={parsed}")

    return parsed


def load_table(
    path: str,
    skp_pfx: tuple[str, ...],
    verbose: bool = False,
) -> list[dict]:
    """
    Load a metadata table as dictionaries keyed by header names.

    Parameters
    ----------
    path : str
        Metadata TSV path.
    skp_pfx : tuple[str, ...]
        Prefixes identifying non-data rows.
    verbose : bool
        Whether to report table-loading details.

    Returns
    -------
    rows : list[dict]
        Metadata rows keyed by unique header names.

    Raises
    ------
    ValueError
        If the header is missing or duplicated, or no data rows remain.
    """

    with open_in(path) as handle:
        text = handle.read()

    lines = [
        line
        for line in text.splitlines()
        if line.strip() and not (skp_pfx and is_header(line, skp_pfx))
    ]

    if not lines:
        raise ValueError("No non-header rows found in metadata table.")

    first_line = lines[0]
    delimiter = "\t" if "\t" in first_line else ","
    reader = csv.DictReader(
        io.StringIO("\n".join(lines)),
        delimiter=delimiter,
    )
    header = reader.fieldnames or []

    if not header:
        raise ValueError("Metadata table is missing a header row.")

    if len(set(header)) != len(header):
        raise ValueError("Metadata table contains duplicate header names.")

    rows = [
        {
            key: value.strip() if isinstance(value, str) else value
            for key, value in row.items()
        }
        for row in reader
    ]

    if verbose:
        with redirect_stdout(sys.stderr):
            print(f"'load_table' header={header}")
            print(f"'load_table' rows loaded={len(rows)}")

    return rows


def find_row_matching(
    rows: list[dict],
    parsed: dict[str, str],
    configuration: dict,
) -> dict:
    """
    Return the single metadata row matching parsed filename fields.

    Parameters
    ----------
    rows : list[dict]
        Loaded metadata rows.
    parsed : dict[str, str]
        Field values parsed from the input filename.
    configuration : dict
        Validated matching configuration.

    Returns
    -------
    row : dict
        The unique matching metadata row.

    Raises
    ------
    ValueError
        If columns are missing or the match count is not exactly one.
    """

    field_columns = configuration.get("field_to_column", {})
    criteria: dict[str, str] = {}

    for field in configuration["matching"]["fields"]:
        column = field_columns.get(field, field)
        criteria[column] = parsed[field]

    missing_columns = sorted(
        {column for column in criteria if not rows or column not in rows[0]},
    )

    if missing_columns:
        raise ValueError(
            "Metadata table is missing matching column(s): "
            + ", ".join(missing_columns),
        )

    matches = [
        row
        for row in rows
        if all(
            str(row.get(column, "")) == value
            for column, value in criteria.items()
        )
    ]

    description = ", ".join(
        f"{column}={value!r}" for column, value in criteria.items()
    )

    if len(matches) == 1:
        return matches[0]

    if not matches:
        raise ValueError(
            f"No metadata row matched criteria: {description}.",
        )

    raise ValueError(
        f"Multiple metadata rows matched criteria: {description}.",
    )


def is_missing(value: object) -> bool:
    """
    Return True when a metadata value is empty or an NA sentinel.
    """

    return value is None or str(value).strip() in {"", "NA", "N/A"}


def collect_outputs(row: dict, configuration: dict) -> dict[str, str]:
    """
    Collect configured calculator inputs from one matched metadata row.

    Parameters
    ----------
    row : dict
        Matched metadata row.
    configuration : dict
        Validated calculator-input configuration.

    Returns
    -------
    outputs : dict[str, str]
        Required and present optional calculator inputs.

    Raises
    ------
    ValueError
        If required values are absent or paired positive values are invalid.
    """

    siq_config = configuration["calculator_inputs"]["siqchip"]
    required_inputs: dict[str, str] = siq_config["required"]
    optional_inputs: dict[str, str] = siq_config["optional"]
    outputs: dict[str, str] = {}

    for key, column in required_inputs.items():
        if column not in row or is_missing(row.get(column)):
            raise ValueError(
                f"Matched metadata row is missing required calculator input "
                f"'{key}' from column '{column}'.",
            )

        outputs[key] = str(row[column])

    for key, column in optional_inputs.items():
        if column in row and not is_missing(row.get(column)):
            outputs[key] = str(row[column])
        else:
            outputs[key] = "NA"

    for key_a, key_b in PAIRED_OPTIONS:
        if key_a in outputs or key_b in outputs:
            missing_a = is_missing(outputs.get(key_a))
            missing_b = is_missing(outputs.get(key_b))

            if missing_a != missing_b:
                raise ValueError(
                    f"Optional metadata fields '{key_a}' and '{key_b}' must "
                    "be provided together or omitted together.",
                )

            if not missing_a:
                for key in (key_a, key_b):
                    try:
                        number = float(outputs[key])
                    except ValueError as error:
                        raise ValueError(
                            f"Optional metadata field '{key}' must be a "
                            "positive number when supplied, but got "
                            f"'{outputs[key]}'.",
                        ) from error

                    if number <= 0:
                        raise ValueError(
                            f"Optional metadata field '{key}' must be > 0 "
                            f"when supplied, but got '{outputs[key]}'.",
                        )

    return outputs


def output_shell(values: dict[str, str]) -> None:
    """
    Print shell export statements for configured metadata values.
    """

    for key in sorted(values):
        validate_id(key, key)

        print(f"export {key}={shlex.quote(str(values[key]))}")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse command-line arguments for metadata row lookup.

    Parameters
    ----------
    argv : list[str] | None
        Explicit arguments, or None to read the process arguments.

    Returns
    -------
    arguments : argparse.Namespace
        Parsed configuration, alignment, and output-selection options.

    Raises
    ------
    SystemExit
        If argument parsing fails or help is requested.
    """

    parser = CapArgumentParser(
        description=(
            "Deterministically map an alignment filename to one siQ-ChIP "
            "metadata row and export configured metadata values."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    add_help_cap(parser)
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Run script in verbose mode. Print parser diagnostics to stderr.",
    )
    parser.add_argument(
        "-a",
        "--alignment",
        dest="alignment",
        help="Alignment file whose basename should be parsed.",
    )
    parser.add_argument(
        "-tb",
        "--tbl_met",
        dest="tbl_met",
        help="siQ-ChIP metadata table containing rows to match.",
    )
    parser.add_argument(
        "-c",
        "--cfg",
        "--configure",
        dest="cfg",
        required=True,
        help="YAML parser configuration file.",
    )
    parser.add_argument(
        "-vc",
        "--validate_cfg",
        dest="validate_cfg",
        action="store_true",
        default=False,
        help="Validate the configuration file and exit.",
    )
    parser.add_argument(
        "-sh",
        "--shell",
        dest="shell",
        action="store_true",
        default=False,
        help="Emit shell export statements instead of key-value lines.",
    )
    parser.add_argument(
        "-sp",
        "--skp_pfx",
        dest="skp_pfx",
        default=None,
        help=(
            "Comma-separated list of header prefixes to skip in metadata "
            "lines."
        ),
    )

    args = parser.parse_args(argv)

    if not args.validate_cfg:
        if not args.alignment:
            parser.error(
                "'--alignment' is required unless '--validate_cfg' is "
                "supplied.",
            )

        if not args.tbl_met:
            parser.error(
                "'--tbl_met' is required unless '--validate_cfg' is supplied.",
            )

    return args


def main(argv: list[str] | None = None) -> int:
    """
    Run metadata lookup and return a process-style exit status.

    Parameters
    ----------
    argv : list[str] | None
        Explicit arguments, or None to read the process arguments.

    Returns
    -------
    status : int
        Zero on success and a stable nonzero status for anticipated failures.
    """

    args = parse_args(argv)

    try:
        check_exists(args.cfg, kind="file", label="--cfg")
        configuration = load_config(args.cfg, verbose=args.verbose)
    except (OSError, FileNotFoundError, ValueError, yaml.YAMLError) as error:
        print(str(error), file=sys.stderr)

        return 1

    if args.validate_cfg:
        print(f"Configuration OK: {args.cfg}")

        return 0

    if args.alignment == "-":
        print("Error: '--alignment -' is not supported.", file=sys.stderr)

        return 1

    try:
        check_exists(args.alignment, kind="file", label="--alignment")
        check_exists(args.tbl_met, kind="file", label="--tbl_met")
    except (FileNotFoundError, ValueError) as error:
        print(str(error), file=sys.stderr)

        return 1

    skp_pfx = parse_skp_pfx(
        args.skp_pfx,
        default=tuple(
            configuration.get("table", {}).get("skip_prefixes", []),
        ),
    )

    try:
        parsed = parse_filename(
            args.alignment,
            configuration,
            args.verbose,
        )
        rows = load_table(args.tbl_met, skp_pfx=skp_pfx, verbose=args.verbose)
        row = find_row_matching(rows, parsed, configuration)

        outputs = collect_outputs(row, configuration)
    except (OSError, ValueError) as error:
        print(str(error), file=sys.stderr)

        return 1

    if args.verbose:
        with redirect_stdout(sys.stderr):
            print(f"'main' matched row={row}")
            print(f"'main' exported values={outputs}")

    if args.shell:
        output_shell(outputs)
    else:
        for key in sorted(outputs):
            print(f"{key}: {outputs[key]}")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        raise SystemExit(0) from None
