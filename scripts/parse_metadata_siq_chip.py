#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: parse_metadata_siq_chip.py
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
#
# Distributed under the MIT license.

"""
Script
-----
parse_metadata_siq_chip.py


Description
-----------
Parses a metadata table (CSV/TSV; raw or gzip-compressed) to extract siQ-ChIP
measurement values (e.g., volume, mass, concentration, fragment length, and
sequencing depth) by parsing identifying attributes from a BAM filename (e.g.,
assay, genotype, state, treatment, factor, and strain/replicate) and matching
them to a metadata table row.


Usage
-----
python parse_metadata_siq_chip.py \\
    [--verbose] \\
    --bam sample.bam \\
    --tbl_met table.tsv[.gz] \\
    --cfg parse_metadata_siq_chip.yaml \\
    [--eqn {5,5nd,6,6nd}] \\
    [--shell] \\
    [--skp_pfx "#,//"]


Arguments
---------
 -v, --verbose
    Enables verbose output for debugging or information.

 -b, --bam
    Input BAM file used to identify metadata row (not stdin).

-tb, --tbl_met
    CSV or TSV siQ-ChIP metadata infile.

 -c, --cfg, --configure
    YAML configuration file defining filename parsing, column aliases, and
    metadata matching.

-eq, --eqn
    Equation to compute: '5', '5nd', '6', or '6nd' (default: 6nd).

-sh, --shell
    Outputs values in a shell-parseable format (for export).

-sp, --skp_pfx
    Comma-separated prefixes to skip as headers/comments in the metadata.
    Overrides configured/default prefixes if supplied. (Example: '#,//'.)


Output
------
When '--shell' is specified, prints key-value pairs in a shell-parseable
format to stdout. Otherwise, prints human-readable measurement values to
stdout.


Examples
--------
1. Export shell-ready values
'''bash
python -m parse_metadata_siq_chip \
    --tbl_met "${HOME}/path/to/table.tsv" \
    --bam "${HOME}/path/to/sample.bam" \
    --shell
'''

2. Return non-shell output with an explicit equation
'''bash
python -m parse_metadata_siq_chip \
    --tbl_met table.csv \
    --bam IP_WT_G1_Isw1-flag_S288C.bam \
    --eqn 5
'''


General notes
-------------
- Filename parsing is controlled by '--cfg'.
- If 'filename.pattern' is non-null, the script uses regex mode with the
  configured named-group pattern.
- Otherwise, the script uses layout mode and tries the user-defined layouts in
  'filename.layouts' in the order given.
- Layout mode supports configurable assay/state/genotype/factor detection
  without requiring a full regular expression.
- 'state' tokens are defined in the YAML configuration under 'filename.states'.
  Matching is case-insensitive.
- 'factor' may be a histone mark (e.g., H3K4me3) or a *-flag protein.
- Column aliases in the table are standardized to canonical names via '--cfg'.


Performance notes
-----------------
- The metadata file is read into memory to allow header standardization and
  filtering.


#TODO
-----
- Continue to expand YAML alias coverage for metadata column names used across
  siQ-ChIP tables.
- Consider documenting the YAML schema more formally for lab users.
- Consider documenting example YAML patterns for alternative filename orders
  used in the lab.
"""

from __future__ import annotations

from contextlib import redirect_stdout

import argparse
import csv
import io
import os
import re
import shlex
import sys

from scripts.functions.utils_check import check_exists
from scripts.functions.utils_cli import add_help_cap, CapArgumentParser
from scripts.functions.utils_io import is_header, open_in, parse_skp_pfx

try:
    import yaml
except ImportError as e:
    raise ImportError(
        "PyYAML is required for 'parse_metadata_siq_chip.py' configuration "
        "support. Add it to the environment or install with, e.g., 'mamba "
        "install conda-forge::pyyaml'."
    ) from e


assert sys.version_info >= (3, 10), "Python >= 3.10 required."


def load_cfg(
    path: str,
    verbose: bool = False
) -> dict:
    """
    Load YAML configuration from 'path' and validate required sections/keys.
    """
    with open(path, "rt", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh) or {}

    if not isinstance(cfg, dict):
        raise ValueError(
            f"Configuration file '{path}' must contain a YAML mapping."
        )

    req_top = (
        "filename",
        "matching",
        "table",
        "columns",
        "required_keys",
        "optional_keys",
    )
    missing = [key for key in req_top if key not in cfg]
    if missing:
        raise ValueError(
            "Configuration file is missing required top-level section(s): "
            + ", ".join(missing)
        )

    req_filename = (
        "delimiter",
        "strip_extensions",
        "assay_tokens",
        "states",
        "pattern_group_defaults",
        "layouts",
        "genotype",
        "factor",
    )
    missing = [
        key for key in req_filename if key not in cfg["filename"]
    ]
    if missing:
        raise ValueError(
            "Configuration file section 'filename' is missing required "
            "key(s): " + ", ".join(missing)
        )

    req_genotype = ("exact", "suffixes", "exclude_suffixes")
    missing = [
        key for key in req_genotype if key not in cfg["filename"]["genotype"]
    ]
    if missing:
        raise ValueError(
            "Configuration file section 'filename.genotype' is missing "
            "required key(s): " + ", ".join(missing)
        )

    req_factor = ("regex", "whitelist")
    missing = [
        key for key in req_factor if key not in cfg["filename"]["factor"]
    ]
    if missing:
        raise ValueError(
            "Configuration file section 'filename.factor' is missing "
            "required key(s): " + ", ".join(missing)
        )

    req_matching = (
        "factor_uses_abbrev",
        "match_treatment",
        "wildcard_values",
    )
    missing = [
        key for key in req_matching if key not in cfg["matching"]
    ]
    if missing:
        raise ValueError(
            "Configuration file section 'matching' is missing required "
            "key(s): " + ", ".join(missing)
        )

    req_table = ("skip_prefixes",)
    missing = [
        key for key in req_table if key not in cfg["table"]
    ]
    if missing:
        raise ValueError(
            "Configuration file section 'table' is missing required key(s): "
            + ", ".join(missing)
        )

    req_pat_dfl = ("genotype", "state", "treatment")
    missing = [
        key for key in req_pat_dfl
        if key not in cfg["filename"]["pattern_group_defaults"]
    ]
    if missing:
        raise ValueError(
            "Configuration file section 'filename.pattern_group_defaults' "
            "is missing required key(s): " + ", ".join(missing)
        )

    req_cols = (
        "factor",
        "factor_abbrev",
        "strain",
        "genotype",
        "state",
        "treatment",
    )
    missing = [key for key in req_cols if key not in cfg["columns"]]
    if missing:
        raise ValueError(
            "Configuration file section 'columns' is missing required key(s): "
            + ", ".join(missing)
        )

    for key in cfg["required_keys"]:
        if key not in cfg["columns"]:
            raise ValueError(
                f"Required key '{key}' listed in 'required_keys' is not "
                "defined in configuration section 'columns'."
            )

    for key in cfg["optional_keys"]:
        if key not in cfg["columns"]:
            raise ValueError(
                f"Optional key '{key}' listed in 'optional_keys' is not "
                "defined in configuration section 'columns'."
            )

    layouts = cfg["filename"].get("layouts")
    if not isinstance(layouts, list):
        raise ValueError(
            "Configuration section 'filename.layouts' must be a YAML list."
        )

    pattern = cfg["filename"].get("pattern")
    if not pattern and not layouts:
        raise ValueError(
            "Configuration must provide either 'filename.pattern' or at least "
            "one entry in 'filename.layouts'."
        )

    valid_fields = {
        "assay",
        "genotype",
        "state",
        "treatment",
        "factor",
        "strain",
    }
    valid_positions = {"first", "last"}

    for i, layout in enumerate(layouts, start=1):
        if not isinstance(layout, dict):
            raise ValueError(
                f"Layout {i} in 'filename.layouts' must be a mapping."
            )

        if "fields" not in layout:
            raise ValueError(
                f"Layout {i} in 'filename.layouts' is missing required key "
                "'fields'."
            )

        fields = layout["fields"]
        if not isinstance(fields, list) or not fields:
            raise ValueError(
                f"Layout {i} in 'filename.layouts.fields' must be a non-empty "
                "list."
            )

        seen_names: set[str] = set()
        n_first = 0
        n_last = 0

        for j, field in enumerate(fields, start=1):
            if not isinstance(field, dict):
                raise ValueError(
                    f"Layout {i}, field {j} in 'filename.layouts.fields' must "
                    "be a mapping."
                )

            name = field.get("name")
            if name not in valid_fields:
                raise ValueError(
                    f"Layout {i}, field {j} has invalid name '{name}'. Valid "
                    "field names are: assay, genotype, state, treatment, "
                    "factor, strain."
                )

            if name in seen_names:
                raise ValueError(
                    f"Layout {i} repeats field name '{name}'. Each field may "
                    "appear at most once per layout."
                )
            seen_names.add(name)

            if "required" in field and not isinstance(field["required"], bool):
                raise ValueError(
                    f"Layout {i}, field '{name}' has non-boolean 'required'."
                )

            pos = field.get("position")
            if pos is not None:
                if pos not in valid_positions:
                    raise ValueError(
                        f"Layout {i}, field '{name}' has invalid position "
                        f"'{pos}'. Valid positions are 'first' and 'last'."
                    )
                if pos == "first":
                    n_first += 1
                elif pos == "last":
                    n_last += 1

        if n_first != 1:
            raise ValueError(
                f"Layout {i} must define exactly one field with "
                "'position: first'."
            )
        if n_last != 1:
            raise ValueError(
                f"Layout {i} must define exactly one field with "
                "'position: last'."
            )

    pattern = cfg["filename"].get("pattern")
    if pattern:
        try:
            rgx = re.compile(pattern)
        except re.error as e:
            raise ValueError(
                "Configuration value 'filename.pattern' is not a valid "
                f"regular expression: {e}"
            ) from e

        req_groups = {"assay", "factor", "strain"}
        missing = sorted(req_groups - set(rgx.groupindex))
        if missing:
            raise ValueError(
                "Configuration regex 'filename.pattern' is missing required "
                "named capture group(s): " + ", ".join(missing)
            )

    if verbose:
        with redirect_stdout(sys.stderr):
            print(f"'load_cfg' path: {path}")
            print(f"'load_cfg' cfg: {cfg}")

    return cfg


def is_state_tok(tok: str, states: tuple[str, ...]) -> bool:
    return tok.lower() in {s.lower() for s in states}


def is_genotype_tok(tok: str, gen_cfg: dict) -> bool:
    tl = tok.lower()
    exact_gen = {s.lower() for s in gen_cfg.get("exact", [])}
    geno_suffixes = tuple(s.lower() for s in gen_cfg.get("suffixes", []))
    geno_exclude = tuple(
        s.lower() for s in gen_cfg.get("exclude_suffixes", [])
    )

    if tl in exact_gen:
        return True
    if any(tl.endswith(sfx) for sfx in geno_exclude):
        return False
    if any(tl.endswith(sfx) for sfx in geno_suffixes):
        return True
    return False


def is_factor_tok(tok: str, fac_cfg: dict) -> bool:
    fac_regex = fac_cfg.get("regex", [])
    fac_whitelist = {s.lower() for s in fac_cfg.get("whitelist", [])}

    if tok.lower() in fac_whitelist:
        return True

    return any(
        re.fullmatch(pat, tok, flags=re.IGNORECASE)
        for pat in fac_regex
    )


def is_assay_tok(tok: str, assay_tokens: tuple[str, ...]) -> bool:
    return tok in assay_tokens


def token_matches_field(
    tok: str,
    field_name: str,
    assay_tokens: tuple[str, ...],
    states: tuple[str, ...],
    gen_cfg: dict,
    fac_cfg: dict,
) -> bool:
    if field_name == "assay":
        return is_assay_tok(tok, assay_tokens)
    if field_name == "genotype":
        return is_genotype_tok(tok, gen_cfg)
    if field_name == "state":
        return is_state_tok(tok, states)
    if field_name == "factor":
        return is_factor_tok(tok, fac_cfg)
    if field_name == "strain":
        return True
    if field_name == "treatment":
        #  Treatment is handled specially in ordered layout parsing because no
        #  lexical recognizer is currently defined for it in the YAML
        return False
    raise ValueError(f"Unrecognized field name '{field_name}'.")


def parse_layout_middle_fields(
    mid_toks: list[str],
    floating_fields: list[dict],
    assay_tokens: tuple[str, ...],
    states: tuple[str, ...],
    gen_cfg: dict,
    fac_cfg: dict,
    layout_name: str,
    stem: str,
) -> dict[str, str]:
    """
    Parse the non-anchored middle tokens of a filename according to the ordered
    floating fields in a layout.

    Optional fields may be skipped. Required fields must consume a compatible
    token. Fields such as 'treatment' may be handled specially by ordered
    backtracking when no lexical recognizer is defined.
    """
    def rec(
        tok_idx: int,
        fld_idx: int,
        out_mid: dict[str, str],
    ) -> dict[str, str] | None:
        if fld_idx == len(floating_fields):
            if tok_idx == len(mid_toks):
                return out_mid
            return None

        field = floating_fields[fld_idx]
        name = field["name"]
        required = field.get("required", False)

        #  Option 1: consume the next token if it matches this field
        if tok_idx < len(mid_toks):
            tok = mid_toks[tok_idx]

            can_consume = token_matches_field(
                tok, name, assay_tokens, states, gen_cfg, fac_cfg
            )
            if name == "treatment":
                can_consume = True

            if can_consume:
                out_next = out_mid.copy()
                out_next[name] = tok
                hit = rec(tok_idx + 1, fld_idx + 1, out_next)
                if hit is not None:
                    return hit

        #  Option 2: skip this field if it is optional
        if not required:
            out_next = out_mid.copy()
            out_next[name] = "NA"
            hit = rec(tok_idx, fld_idx + 1, out_next)
            if hit is not None:
                return hit

        return None

    parsed = rec(0, 0, {})
    if parsed is None:
        raise ValueError(
            f"Layout '{layout_name}' could not parse middle tokens "
            f"{mid_toks} for filename stem '{stem}'."
        )

    return parsed


def parse_bam_filename_layout(
    stem: str,
    layout: dict,
    cfg: dict,
) -> dict[str, str]:
    fil_cfg = cfg["filename"]
    delim = fil_cfg.get("delimiter", "_")
    toks = stem.split(delim)

    if len(toks) < 2:
        raise ValueError(f"Filename stem '{stem}' is too short to parse.")

    assay_tokens = tuple(fil_cfg.get("assay_tokens", ["IP", "in"]))
    states = tuple(fil_cfg.get("states", []))
    gen_cfg = fil_cfg.get("genotype", {})
    fac_cfg = fil_cfg.get("factor", {})

    fields = layout["fields"]
    out = {
        "assay": "NA",
        "genotype": "NA",
        "state": "NA",
        "treatment": "NA",
        "factor": "NA",
        "strain": "NA",
    }

    layout_name = layout.get("name", "unnamed")
    floating_fields: list[dict] = []

    first_name = None
    last_name = None

    for field in fields:
        name = field["name"]
        pos = field.get("position")

        if pos == "first":
            first_name = name
        elif pos == "last":
            last_name = name
        else:
            floating_fields.append(field)

    if first_name is None or last_name is None:
        raise ValueError(
            f"Layout '{layout_name}' must define both a first-position and "
            "last-position field."
        )

    first_tok = toks[0]
    last_tok = toks[-1]

    if not token_matches_field(
        first_tok, first_name, assay_tokens, states, gen_cfg, fac_cfg
    ):
        raise ValueError(
            f"Layout '{layout_name}' expected field '{first_name}' in first "
            f"position, but token '{first_tok}' did not match."
        )

    if not token_matches_field(
        last_tok, last_name, assay_tokens, states, gen_cfg, fac_cfg
    ):
        raise ValueError(
            f"Layout '{layout_name}' expected field '{last_name}' in last "
            f"position, but token '{last_tok}' did not match."
        )

    out[first_name] = first_tok
    out[last_name] = last_tok

    mid_toks = toks[1:-1]
    parsed_mid = parse_layout_middle_fields(
        mid_toks,
        floating_fields=floating_fields,
        assay_tokens=assay_tokens,
        states=states,
        gen_cfg=gen_cfg,
        fac_cfg=fac_cfg,
        layout_name=layout_name,
        stem=stem,
    )
    out.update(parsed_mid)

    return out


def parse_bam_filename(
    filename: str,
    cfg: dict,
    verbose: bool = False
) -> dict[str, str]:
    """
    Parse a BAM filename into its components.

    Supports two modes:
        1. Regex mode, if 'cfg["filename"]["pattern"]' is non-null.
        2. Otherwise, layout mode, which tries the configured layouts in
           'cfg["filename"]["layouts"]' in order.

    Returns:
        dict[str, str]:
            {'assay', 'genotype', 'state', 'treatment', 'factor', 'strain'}
    """
    fil_cfg = cfg["filename"]
    nam_bas = os.path.basename(filename)

    stem = nam_bas
    for ext in fil_cfg.get("strip_extensions", []):
        if stem.endswith(ext):
            stem = stem[: -len(ext)]
            break

    strip_tag_re = fil_cfg.get("strip_trailing_tag_regex")
    if strip_tag_re:
        stem = re.sub(strip_tag_re, "", stem)

    pattern = fil_cfg.get("pattern")
    if pattern:
        m = re.fullmatch(pattern, stem)
        if not m:
            raise ValueError(
                f"Filename '{filename}' does not match configured regex "
                f"pattern."
            )

        grp = m.groupdict()
        dflt = fil_cfg.get("pattern_group_defaults", {})

        assay = grp.get("assay")
        factor = grp.get("factor")
        strain = grp.get("strain")

        if not assay:
            raise ValueError(
                f"Configured regex did not capture required group 'assay' "
                f"from filename '{filename}'."
            )
        if not factor:
            raise ValueError(
                f"Configured regex did not capture required group 'factor' "
                f"from filename '{filename}'."
            )
        if not strain:
            raise ValueError(
                f"Configured regex did not capture required group 'strain' "
                f"from filename '{filename}'."
            )

        assay_tokens = tuple(fil_cfg.get("assay_tokens", ["IP", "in"]))
        if assay not in assay_tokens:
            raise ValueError(
                f"Assay must be one of {assay_tokens}, but got '{assay}' "
                f"in '{filename}'."
            )

        out = {
            "assay": assay,
            "genotype": grp.get("genotype") or dflt.get("genotype", "NA"),
            "state": grp.get("state") or dflt.get("state", "NA"),
            "treatment": grp.get("treatment") or dflt.get("treatment", "NA"),
            "factor": factor,
            "strain": strain,
        }

        if verbose:
            with redirect_stdout(sys.stderr):
                print(f"'parse_bam_filename' regex stem={stem}")
                print(f"'parse_bam_filename' regex groups={grp}")
                print(f"'parse_bam_filename' resolved={out}")

        return out

    layouts = fil_cfg.get("layouts", [])
    errs: list[str] = []

    for layout in layouts:
        try:
            out = parse_bam_filename_layout(
                stem,
                layout=layout,
                cfg=cfg,
            )

            if verbose:
                with redirect_stdout(sys.stderr):
                    print(
                        "'parse_bam_filename' layout mode: "
                        f"name={layout.get('name', 'unnamed')}"
                    )
                    print(f"'parse_bam_filename' layout stem={stem}")
                    print(f"'parse_bam_filename' resolved={out}")

            return out

        except ValueError as e:
            errs.append(str(e))

    raise ValueError(
        f"Could not parse filename '{filename}' using configured layout mode. "
        f"Tried {len(layouts)} layout(s). Last error: {errs[-1]}"
    )


def standardize_header(
    head_raw: list[str],
    cfg: dict
) -> dict[str, str]:
    """
    Build a mapping from raw header names to canonical names using the
    configured column-alias mapping.
    """
    col_map = cfg["columns"]

    lookup: dict[str, str] = {}
    for canon, alts in col_map.items():
        for alias in (alts + [canon]):
            lookup[alias.lower()] = canon

    std_map: dict[str, str] = {}
    for col in head_raw:
        key = col.strip().lower()
        std_map[col] = lookup.get(key, col)

    return std_map


def load_table(
    path: str,
    skp_pfx: tuple[str, ...],
    cfg: dict,
    verbose: bool = False
) -> list[dict]:
    """
    Load CSV/TSV file as a list of dictionaries, standardizing column names
    and rejecting canonical-name collisions after alias normalization.
    """
    with open_in(path) as fh:
        buf = fh.read()

    lines: list[str] = []
    for ln in buf.splitlines():
        if not ln.strip():
            continue
        if skp_pfx and is_header(ln, skp_pfx):
            continue
        lines.append(ln)

    if not lines:
        raise ValueError("No non-header rows found in table.")

    head = lines[0]
    delim = "\t" if "\t" in head else ","

    rdr = csv.DictReader(io.StringIO("\n".join(lines)), delimiter=delim)
    head_raw = rdr.fieldnames or []
    std_map = standardize_header(head_raw, cfg)

    canon_to_raw: dict[str, list[str]] = {}
    for raw in head_raw:
        canon = std_map.get(raw, raw)
        canon_to_raw.setdefault(canon, []).append(raw)

    collisions = {
        canon: raws
        for canon, raws in canon_to_raw.items()
        if len(raws) > 1
    }
    if collisions:
        msg_parts = []
        for canon, raws in sorted(collisions.items()):
            msg_parts.append(
                f"{canon}: {', '.join(repr(r) for r in raws)}"
            )
        raise ValueError(
            "Metadata table contains multiple columns that normalize to the "
            "same canonical name: " + "; ".join(msg_parts)
        )

    if verbose:
        with redirect_stdout(sys.stderr):
            std_nam = [std_map.get(h, h) for h in head_raw]
            print(f"'load_table' header (raw): {head_raw}")
            print(f"'load_table' header (std): {std_nam}")
            print(f"'load_table' std_map: {std_map}")

    out: list[dict] = []
    for row in rdr:
        out.append({
            std_map.get(k, k): (v.strip() if isinstance(v, str) else v)
            for k, v in row.items()
        })

    if verbose:
        with redirect_stdout(sys.stderr):
            print(f"'load_table' rows loaded: {len(out)}")

    return out


def find_matching_row(
    data: list[dict[str, str]],
    factor: str,
    strain: str,
    cfg: dict,
    genotype: str = "NA",
    state: str = "NA",
    treatment: str = "NA",
    verbose: bool = False
) -> dict[str, str]:
    """
    Find the first row matching 'factor' (or 'factor_abbrev'), 'strain', and
    optional 'genotype', 'state', and/or 'treatment' fields. Matching is
    case-insensitive; configured wildcard values act as wildcards, and
    'treatment' matching is applied only when enabled in the YAML
    configuration.
    """
    match_cfg = cfg.get("matching", {})
    use_abbrev = bool(match_cfg.get("factor_uses_abbrev", True))
    match_treatment = bool(match_cfg.get("match_treatment", False))

    def standardize(s: str | None) -> str:
        return (s or "").strip().lower()

    wildcards = {
        standardize(v)
        for v in match_cfg.get("wildcard_values", ["", "NA", "N/A"])
    }

    sf, ss, sg, sst, strt = map(
        standardize, (factor, strain, genotype, state, treatment)
    )
    wild_sg = sg in wildcards
    wild_sst = sst in wildcards
    wild_strt = strt in wildcards

    if verbose:
        with redirect_stdout(sys.stderr):
            print(
                f"'find_matching_row' match search: factor={sf}, strain={ss}, "
                f"genotype={sg or 'n/a'}, state={sst or 'n/a'}, "
                f"treatment={strt or 'n/a'}"
            )

    def row_matches(row: dict[str, str], fac: str, st: str) -> bool:
        rfac = standardize(row.get("factor"))
        rfab = standardize(row.get("factor_abbrev"))
        rstr = standardize(row.get("strain"))
        rgen = standardize(row.get("genotype"))
        rstate = standardize(row.get("state"))
        rtrt = standardize(row.get("treatment"))

        fac_ok = (fac == rfac) or (use_abbrev and fac == rfab)
        gen_ok = wild_sg or (sg == rgen)
        st_ok = wild_sst or (st == rstate)
        trt_ok = (not match_treatment) or wild_strt or (strt == rtrt)

        if verbose and (rstr == ss):
            with redirect_stdout(sys.stderr):
                print(
                    f"'find_matching_row' matches: row strain={rstr}, "
                    f"fac=({rfac}|{rfab}), gen={rgen}, state={rstate}, "
                    f"treatment={rtrt}"
                )
                print(
                    f"'find_matching_row' status: fac_ok={fac_ok}, "
                    f"gen_ok={gen_ok}, st_ok={st_ok}, trt_ok={trt_ok}"
                )

        return fac_ok and (ss == rstr) and gen_ok and st_ok and trt_ok

    for row in data:
        if row_matches(row, sf, sst):
            if verbose:
                with redirect_stdout(sys.stderr):
                    print("'find_matching_row' status: Pass 1 matched")
            return row

    for row in data:
        if row_matches(row, sst, sf):
            if verbose:
                with redirect_stdout(sys.stderr):
                    print(
                        "'find_matching_row' status: Pass 2 matched with "
                        "swapped factor/state"
                    )
            return row

    raise ValueError(
        f"No matching row found for state '{state or 'NA'}', "
        f"factor '{factor}', strain '{strain}', genotype "
        f"'{genotype or 'NA'}', treatment '{treatment or 'NA'}'."
    )


def output_for_shell(**kwargs):
    """
    Print key/value pairs in a shell-parseable format, one per line:
        export key=value
    """
    for key, value in kwargs.items():
        print(f"export {key}={shlex.quote(str(value))}")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse command line arguments.
    """
    parser = CapArgumentParser(
        description=(
            "Parse a BAM filename and retrieve a matching row from a TSV/CSV "
            "table of siQ-ChIP metadata."
        )
    )
    add_help_cap(parser)
    parser.add_argument(
        "-v", "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Increase output verbosity.\n\n"
    )
    parser.add_argument(
        "-b", "--bam",
        dest="bam",
        help="Path to input BAM file (not stdin).\n\n",
        required=True
    )
    parser.add_argument(
        "-tb", "--tbl-met", "--tbl_met",
        dest="tbl_met",
        required=True,
        help="TSV or CSV table of siQ-ChIP metadata; .gz is handled.\n\n"
    )
    parser.add_argument(
        "-c", "--cfg", "--configure",
        dest="cfg",
        type=str,
        required=True,
        help=(
            "YAML configuration file defining filename parsing, column "
            "aliases, and metadata matching. If 'filename.pattern' is set, "
            "the script uses regex mode; otherwise it uses layout mode and "
            "tries the configured layouts in order.\n\n"
        )
    )
    parser.add_argument(
        '-eq', '--eqn',
        dest="eqn",
        type=str,
        choices=['5', '5nd', '6', '6nd'],
        default='6nd',
        help=(
            "Equation to compute the alpha scaling factor (PMID: 37160995; "
            "default: %(default)s). Options:\n"
            "\n"
            "  - '5' applies Equation 5 for use with fragment length-"
            "normalized coverage,\n"
            "  - '5nd' uses Equation 5 without depth terms for use with "
            "“normalized coverage”,\n"
            "  - '6' applies Equation 6 for use with fragment length-"
            "normalized coverage, and\n"
            "  - '6nd' uses Equation 6 without depth terms for use with "
            "“normalized coverage”.\n\n"
        )
    )
    parser.add_argument(
        "-sh", "--shell",
        dest="shell",
        action="store_true",
        help="Output values in a shell-parseable format.\n\n"
    )
    parser.add_argument(
        "-sp", "--skp-pfx", "--skp_pfx",
        dest="skp_pfx",
        type=str,
        default=None,
        help=(
            "Comma-separated prefixes to skip as headers/comments in the "
            "metadata. Overrides configured/default prefixes if supplied. "
            "(Example: '#,//'.)\n\n"
        )
    )

    argv_parse = sys.argv[1:] if argv is None else argv
    if not argv_parse:
        parser.print_help(sys.stderr)
        sys.exit(0)

    return parser.parse_args(argv_parse)


def main(argv: list[str] | None = None) -> int:
    """
    Execute the primary control flow for the script.
    """
    #  Parse CLI arguments
    args = parse_args(argv)

    #  Disallow BAM from stdin: its filename is explicitly parsed
    if args.bam == "-":
        print(
            "Error: '--bam -' not supported; provide a real BAM path.",
            file=sys.stderr
        )
        return 1

    #  Perform BAM, metadata table, configuration YAML existence checks
    check_exists(args.bam, "--bam")
    check_exists(args.tbl_met, "--tbl_met")
    check_exists(args.cfg, "--cfg")

    #  Load and validate YAML configuration
    try:
        cfg = load_cfg(args.cfg, verbose=args.verbose)
    except (OSError, ValueError, yaml.YAMLError) as e:
        print(str(e), file=sys.stderr)
        return 1

    #  Standardize skip prefixes once
    skp_pfx = parse_skp_pfx(
        args.skp_pfx,
        default=tuple(cfg.get("table", {}).get("skip_prefixes", []))
    )

    #  Return “verbose banner” to stderr
    if args.verbose:
        with redirect_stdout(sys.stderr):
            print("#############################################")
            print("## Arguments for 'parse_metadata_siq_chip' ##")
            print("#############################################")
            print("")
            print("--verbose")
            print(f"--bam      {args.bam}")
            print(f"--tbl_met  {args.tbl_met}")
            print(f"--cfg      {args.cfg}")
            print(f"--eqn      {args.eqn}")
            if args.shell:
                print("--shell")
            print(f"--skp_pfx  {skp_pfx}")
            print("")
            print("")

    try:
        comp = parse_bam_filename(
            args.bam,
            cfg=cfg,
            verbose=args.verbose
        )
    except ValueError as e:
        print(str(e), file=sys.stderr)
        return 1

    assay = comp['assay']
    genotype = comp.get('genotype', 'NA')
    state = comp.get('state', 'NA')
    treatment = comp.get('treatment', 'NA')
    factor = comp['factor']
    strain = comp['strain']

    if args.verbose:
        with redirect_stdout(sys.stderr):
            print("####################################")
            print("## Parsed BAM filename components ##")
            print("####################################")
            print("")
            print(f"eqn={args.eqn}\n")
            print(f"assay={assay}\n")
            print(f"genotype={genotype}\n")
            print(f"state={state}\n")
            print(f"treatment={treatment}\n")
            print(f"factor={factor}\n")
            print(f"strain/replicate={strain}")
            print("")
            print("")

    #  Load the TSV/CSV table
    try:
        data = load_table(
            args.tbl_met,
            skp_pfx=skp_pfx,
            cfg=cfg,
            verbose=args.verbose
        )
    except (OSError, ValueError) as e:
        print(str(e), file=sys.stderr)
        return 1

    #  Find the matching row
    try:
        row = find_matching_row(
            data,
            factor,
            strain,
            cfg=cfg,
            genotype=genotype,
            state=state,
            treatment=treatment,
            verbose=args.verbose
        )
        if args.verbose:
            with redirect_stdout(sys.stderr):
                print(f"Matching row found:\n{row}")
    except ValueError as e:
        print(e, file=sys.stderr)
        return 1

    #  Check that required fields are present
    req_keys = tuple(cfg.get("required_keys", []))
    opt_keys = tuple(cfg.get("optional_keys", []))

    missing = [
        k for k in req_keys
        if (row.get(k) is None or str(row.get(k)).strip() == "")
    ]
    if missing:
        print(
            "Missing required column(s) in matched row: "
            + ", ".join(missing),
            file=sys.stderr
        )
        return 1

    #  Fill optional holes for stable downstream behavior
    for k in opt_keys:
        row.setdefault(k, "NA")

    #  Output based on '--shell' flag
    if args.shell:
        output_for_shell(
            eqn=args.eqn,
            vol_in=row['vol_in'],
            vol_all=row['vol_all'],
            mass_in=row['mass_in'],
            mass_ip=row['mass_ip'],
            conc_in=row.get('conc_in', 'NA'),
            conc_ip=row.get('conc_ip', 'NA'),
            len_in=row['len_in'],
            len_ip=row['len_ip'],
            dep_in=row.get('dep_in', 'NA'),
            dep_ip=row.get('dep_ip', 'NA')
        )
    else:
        print(f"vol_in: {row['vol_in']}")
        print(f"vol_all: {row['vol_all']}")
        print(f"mass_in: {row['mass_in']}")
        print(f"mass_ip: {row['mass_ip']}")
        print(f"conc_in: {row.get('conc_in', 'NA')}")
        print(f"conc_ip: {row.get('conc_ip', 'NA')}")
        print(f"len_in: {row['len_in']}")
        print(f"len_ip: {row['len_ip']}")
        print(f"dep_in: {row.get('dep_in', 'NA')}")
        print(f"dep_ip: {row.get('dep_ip', 'NA')}")

    return 0


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
