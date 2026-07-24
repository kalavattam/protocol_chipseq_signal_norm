#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: parse_metadata_siqchip.py
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
# GPT-5.6) were used in development and documentation.
#
# Distributed under the MIT license.


"""
Parse one alignment filename against siQ-ChIP metadata.

Returns
-------
Prints configured metadata-derived calculator inputs for the uniquely matched
metadata row.

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
        "project environment, e.g., with 'mamba install pyyaml'."
    ) from e


assert sys.version_info >= (3, 11), "Python >= 3.11 required."


ID_RGX = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
OPT_PAIRED = (
    ("dep_ip", "dep_in"),
    ("len_ip", "len_in"),
    ("lib_vol_ip", "lib_vol_in"),
)


def validate_id(name: str, label: str) -> None:
    """
    Validate that a configured name can be exported as a shell identifier.
    """
    if not isinstance(name, str) or not ID_RGX.fullmatch(name):
        raise ValueError(f"'{label}' must be a shell-safe identifier.")


def require_mapping(obj: object, label: str) -> dict:
    """Return a configuration section after confirming it is a mapping."""
    if not isinstance(obj, dict):
        raise ValueError(f"Configuration section '{label}' must be a mapping.")
    return obj


def require_list_string(
    obj: object, label: str, allow_empty: bool = False
) -> list[str]:
    """Return a YAML list after confirming every item is a non-empty string."""
    if not isinstance(obj, list):
        raise ValueError(f"Configuration key '{label}' must be a YAML list.")
    if not allow_empty and not obj:
        raise ValueError(
            f"Configuration key '{label}' must be a non-empty list."
        )
    out: list[str] = []
    for i, val in enumerate(obj, start=1):
        if not isinstance(val, str) or not val:
            raise ValueError(
                f"Configuration key '{label}' item {i} must be a non-empty "
                "string."
            )
        out.append(val)
    return out


def normalize_input_map(spec: object, label: str) -> dict[str, str]:
    """
    Normalize calculator input declarations to output-name -> metadata-column.
    Lists mean identity mappings. Mappings allow explicit exported names.
    """
    out: dict[str, str] = {}
    if spec is None:
        return out

    if isinstance(spec, list):
        for key in require_list_string(spec, label, allow_empty=True):
            validate_id(key, f"{label}.{key}")
            out[key] = key
        return out

    if isinstance(spec, dict):
        for key, col in spec.items():
            validate_id(str(key), f"{label}.{key}")
            if not isinstance(col, str) or not col:
                raise ValueError(
                    f"Configuration key '{label}.{key}' must name one "
                    "metadata column."
                )
            out[str(key)] = col
        return out

    raise ValueError(
        f"Configuration key '{label}' must be a YAML list or mapping."
    )


def load_cfg(path: str, verbose: bool = False) -> dict:
    """
    Load and validate the siQ-ChIP metadata parser configuration.
    """
    with open(path, encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh) or {}

    if not isinstance(cfg, dict):
        raise ValueError(
            f"Configuration file '{path}' must contain a YAML mapping."
        )

    for section in ("filename", "matching", "calculator_inputs"):
        if section not in cfg:
            raise ValueError(
                "Configuration file is missing required top-level section "
                f"'{section}'."
            )

    cfg_fil = require_mapping(
        cfg["filename"], "filename"
    )
    cfg_mch = require_mapping(
        cfg["matching"], "matching"
    )
    cfg_clc = require_mapping(
        cfg["calculator_inputs"], "calculator_inputs"
    )
    cfg_siq = require_mapping(
        cfg_clc.get("siqchip"), "calculator_inputs.siqchip"
    )

    delim = cfg_fil.get("delimiter", "_")
    if not isinstance(delim, str) or delim == "":
        raise ValueError(
            "Configuration key 'filename.delimiter' must be a string."
        )
    cfg_fil["delimiter"] = delim

    cfg_fil["strip_extensions"] = require_list_string(
        cfg_fil.get("strip_extensions", []),
        "filename.strip_extensions",
        allow_empty=True,
    )
    cfg_fil["strip_suffixes"] = require_list_string(
        cfg_fil.get("strip_suffixes", []),
        "filename.strip_suffixes",
        allow_empty=True,
    )
    flds = require_list_string(cfg_fil.get("fields"), "filename.fields")
    seen: set[str] = set()
    for fld in flds:
        validate_id(fld, f"filename.fields.{fld}")
        if fld in seen:
            raise ValueError(
                f"Configuration key 'filename.fields' repeats '{fld}'."
            )
        seen.add(fld)
    cfg_fil["fields"] = flds

    fld_mch = require_list_string(
        cfg_mch.get("fields"), "matching.fields"
    )
    for fld in fld_mch:
        validate_id(fld, f"matching.fields.{fld}")
        if fld not in flds:
            raise ValueError(
                f"Configuration key 'matching.fields' references '{fld}', "
                "which is not declared in 'filename.fields'."
            )
    cfg_mch["fields"] = fld_mch

    fld_col = cfg.get("field_to_column", {})
    if fld_col is None:
        fld_col = {}
    fld_col = require_mapping(fld_col, "field_to_column")
    for fld, col in fld_col.items():
        validate_id(str(fld), f"field_to_column.{fld}")
        if fld not in flds:
            raise ValueError(
                f"Configuration key 'field_to_column.{fld}' references a "
                "field not declared in 'filename.fields'."
            )
        if not isinstance(col, str) or not col:
            raise ValueError(
                f"Configuration key 'field_to_column.{fld}' must name one "
                "column."
            )
    cfg["field_to_column"] = fld_col

    req = normalize_input_map(
        cfg_siq.get("required", {}),
        "calculator_inputs.siqchip.required",
    )
    opt = normalize_input_map(
        cfg_siq.get("optional", {}),
        "calculator_inputs.siqchip.optional",
    )
    overlap = sorted(set(req) & set(opt))
    if overlap:
        raise ValueError(
            "Calculator input(s) listed as both required and optional: "
            + ", ".join(overlap)
        )
    cfg_siq["required"] = req
    cfg_siq["optional"] = opt

    cfg_tbl = cfg.get("table", {})
    if cfg_tbl is None:
        cfg_tbl = {}
    cfg_tbl = require_mapping(cfg_tbl, "table")
    cfg_tbl["skip_prefixes"] = require_list_string(
        cfg_tbl.get("skip_prefixes", ["#"]),
        "table.skip_prefixes",
        allow_empty=True,
    )
    cfg["table"] = cfg_tbl

    if verbose:
        with redirect_stdout(sys.stderr):
            print(f"'load_cfg' path: {path}")
            print(f"'load_cfg' cfg: {cfg}")

    return cfg


def strip_filename(path: str, cfg: dict) -> str:
    """
    Remove configured extensions and suffixes from an alignment filename.
    """
    stem = os.path.basename(path)
    for ext in cfg["filename"]["strip_extensions"]:
        if stem.endswith(ext):
            stem = stem[: -len(ext)]
            break
    for sfx in cfg["filename"]["strip_suffixes"]:
        if stem.endswith(sfx):
            stem = stem[: -len(sfx)]
            break
    return stem


def parse_filename(
    path: str, cfg: dict, verbose: bool = False
) -> dict[str, str]:
    """
    Parse an alignment filename into configured metadata fields.
    """
    stem = strip_filename(path, cfg)
    delim = cfg["filename"]["delimiter"]
    flds = cfg["filename"]["fields"]
    toks = stem.split(delim)

    if len(toks) != len(flds):
        raise ValueError(
            f"Filename stem '{stem}' split into {len(toks)} token(s), but "
            f"configuration declares {len(flds)} field(s): "
            f"{', '.join(flds)}."
        )

    parsed = dict(zip(flds, toks, strict=True))

    if verbose:
        with redirect_stdout(sys.stderr):
            print(f"'parse_filename' stem={stem}")
            print(f"'parse_filename' parsed={parsed}")

    return parsed


def load_table(
    path: str, skp_pfx: tuple[str, ...], verbose: bool = False
) -> list[dict]:
    """
    Load a metadata table as dictionaries keyed by header names.
    """
    with open_in(path) as fh:
        buf = fh.read()

    lines = [
        line for line in buf.splitlines()
        if line.strip() and not (skp_pfx and is_header(line, skp_pfx))
    ]
    if not lines:
        raise ValueError("No non-header rows found in metadata table.")

    head = lines[0]
    delim = "\t" if "\t" in head else ","
    rdr = csv.DictReader(io.StringIO("\n".join(lines)), delimiter=delim)
    hdr = rdr.fieldnames or []
    if not hdr:
        raise ValueError("Metadata table is missing a header row.")
    if len(set(hdr)) != len(hdr):
        raise ValueError("Metadata table contains duplicate header names.")

    rows = [
        {k: (v.strip() if isinstance(v, str) else v) for k, v in row.items()}
        for row in rdr
    ]

    if verbose:
        with redirect_stdout(sys.stderr):
            print(f"'load_table' header={hdr}")
            print(f"'load_table' rows loaded={len(rows)}")

    return rows


def find_row_matching(
    rows: list[dict], parsed: dict[str, str], cfg: dict
) -> dict:
    """
    Return the single metadata row matching parsed filename fields.
    """
    fld_col = cfg.get("field_to_column", {})
    criteria: dict[str, str] = {}

    for fld in cfg["matching"]["fields"]:
        col = fld_col.get(fld, fld)
        criteria[col] = parsed[fld]

    col_missing = sorted({
        col for col in criteria
        if not rows or col not in rows[0]
    })
    if col_missing:
        raise ValueError(
            "Metadata table is missing matching column(s): "
            + ", ".join(col_missing)
        )

    matches = [
        row for row in rows
        if all(str(row.get(col, "")) == val for col, val in criteria.items())
    ]

    desc = ", ".join(f"{col}={val!r}" for col, val in criteria.items())
    if len(matches) == 1:
        return matches[0]
    if not matches:
        raise ValueError(f"No metadata row matched criteria: {desc}.")
    raise ValueError(f"Multiple metadata rows matched criteria: {desc}.")


def is_missing(value: object) -> bool:
    """Return True when a metadata value is empty or an NA sentinel."""
    return value is None or str(value).strip() in {"", "NA", "N/A"}


def collect_outputs(row: dict, cfg: dict) -> dict[str, str]:
    """Collect configured calculator inputs from one matched metadata row."""
    siq = cfg["calculator_inputs"]["siqchip"]
    req: dict[str, str] = siq["required"]
    opt: dict[str, str] = siq["optional"]
    out: dict[str, str] = {}

    for key, col in req.items():
        if col not in row or is_missing(row.get(col)):
            raise ValueError(
                f"Matched metadata row is missing required calculator input "
                f"'{key}' from column '{col}'."
            )
        out[key] = str(row[col])

    for key, col in opt.items():
        if col in row and not is_missing(row.get(col)):
            out[key] = str(row[col])
        else:
            out[key] = "NA"

    for key_a, key_b in OPT_PAIRED:
        if key_a in out or key_b in out:
            miss_a = is_missing(out.get(key_a))
            miss_b = is_missing(out.get(key_b))
            if miss_a != miss_b:
                raise ValueError(
                    f"Optional metadata fields '{key_a}' and '{key_b}' must "
                    "be provided together or omitted together."
                )
            if not miss_a:
                for key in (key_a, key_b):
                    try:
                        num = float(out[key])
                    except ValueError as e:
                        raise ValueError(
                            f"Optional metadata field '{key}' must be a "
                            "positive number when supplied, but got "
                            f"'{out[key]}'."
                        ) from e
                    if num <= 0:
                        raise ValueError(
                            f"Optional metadata field '{key}' must be > 0 "
                            f"when supplied, but got '{out[key]}'."
                        )

    return out


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
        "-v", "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Run script in verbose mode. Print parser diagnostics to stderr.",
    )
    parser.add_argument(
        "-a", "--alignment",
        dest="alignment",
        help="Alignment file whose basename should be parsed.",
    )
    parser.add_argument(
        "-tb", "--tbl_met",
        dest="tbl_met",
        help="siQ-ChIP metadata table containing rows to match.",
    )
    parser.add_argument(
        "-c", "--cfg", "--configure",
        dest="cfg",
        required=True,
        help="YAML parser configuration file.",
    )
    parser.add_argument(
        "-vc", "--validate_cfg",
        dest="validate_cfg",
        action="store_true",
        default=False,
        help="Validate the configuration file and exit.",
    )
    parser.add_argument(
        "-sh", "--shell",
        dest="shell",
        action="store_true",
        default=False,
        help="Emit shell export statements instead of key-value lines.",
    )
    parser.add_argument(
        "-sp", "--skp_pfx",
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
                "supplied."
            )
        if not args.tbl_met:
            parser.error(
                "'--tbl_met' is required unless '--validate_cfg' is supplied."
            )
    return args


def main(argv: list[str] | None = None) -> int:
    """Run metadata lookup and return a process-style exit status."""
    args = parse_args(argv)

    try:
        check_exists(args.cfg, kind="file", label="--cfg")
        cfg = load_cfg(args.cfg, verbose=args.verbose)
    except (OSError, FileNotFoundError, ValueError, yaml.YAMLError) as e:
        print(str(e), file=sys.stderr)
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
    except (FileNotFoundError, ValueError) as e:
        print(str(e), file=sys.stderr)
        return 1

    skp_pfx = parse_skp_pfx(
        args.skp_pfx,
        default=tuple(cfg.get("table", {}).get("skip_prefixes", [])),
    )

    try:
        parsed = parse_filename(args.alignment, cfg, args.verbose)
        rows = load_table(args.tbl_met, skp_pfx=skp_pfx, verbose=args.verbose)
        row = find_row_matching(rows, parsed, cfg)
        out = collect_outputs(row, cfg)
    except (OSError, ValueError) as e:
        print(str(e), file=sys.stderr)
        return 1

    if args.verbose:
        with redirect_stdout(sys.stderr):
            print(f"'main' matched row={row}")
            print(f"'main' exported values={out}")

    if args.shell:
        output_shell(out)
    else:
        for key in sorted(out):
            print(f"{key}: {out[key]}")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        raise SystemExit(0) from None
