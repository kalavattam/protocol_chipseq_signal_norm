#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: help_combine_parts_scaling_factor.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6);
# - Anthropic Claude Code (Opus 5).
#
# Distributed under the MIT license.


function help_combine_parts_scaling_factor() {
    cat << EOM
Usage
-----
  combine_parts_scaling_factor.sh
    [--help] [--dry_run]
    --mode <mode> --csv_fil_in <csv> --fil_out <file>
    [--force] [--no_parts]

  Combine generated calculate_scaling_factor part files into one data-only tab-delimited output table.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  -dr, --dry, --dry_run : flag
    Run script in dry-run mode. Validate inputs and print the ordered combination plan without writing output.

  -md, --mode : {'siq', 'spike'}
    Workflow mode for scaling-factor part files. Required.

  -ci, --csv_fil_in : list of file
    Comma-separated list of input file paths for generated scaling-factor part files. Required.

  -fo, --fil_out : file
    Output file path. Final tab-delimited output table. Required.

  -f, --force : flag
    Replace an existing output file only when this flag is supplied.

  -np, --no_parts : flag
    Retain validated part files by default; remove them after successfully writing the final table when this flag is supplied.

Notes
-----
  Runtime requirements:
    - awk
    - basename
    - bash >= 4.4
    - cat
    - dirname
    - mktemp
    - mv
    - rm
    - sed
    - sort

  - Input basenames must end in '.part.<digits>'.
  - Each input must contain exactly one non-empty core scaling-factor row.
  - Core rows contain 12 or 14 fields for 'siq' and 10 fields for 'spike'.
  - Rows are ordered by numeric part-file index, not by input-list order.
  - This script does not write a header; use 'write_header.sh' separately or let 'execute_calculate_scaling_factor.sh' orchestrate header insertion.
  - Duplicate paths, duplicate indexes, and header rows are rejected.

Examples
--------
  1. Combine spike-in scaling-factor parts.
    '''bash
    bash combine_parts_scaling_factor.sh \\
        --mode spike \\
        --csv_fil_in results.spike.tsv.part.000001,results.spike.tsv.part.000000 \\
        --fil_out results.tsv
    '''

  2. Preview siQ-ChIP part combination.
    '''bash
    bash combine_parts_scaling_factor.sh \\
        --dry_run \\
        --mode siq \\
        --csv_fil_in results.siq.tsv.part.000001,results.siq.tsv.part.000000 \\
        --fil_out results.tsv
    '''
EOM
}
