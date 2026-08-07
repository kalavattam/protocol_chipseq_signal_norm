#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: help_write_header.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in design, development,
# and documentation, with all output reviewed, edited, and approved by the
# author.
#
# Distributed under the MIT license.


function help_write_header() {
    # The production owner initializes interpolated defaults before invocation.
    # shellcheck disable=SC2154
    cat << EOM
Usage
-----
  write_header.sh
    [--help] [--verbose] [--dry_run] [--mode <mode>] [--fil_in <file>] [--fil_out <file>] [--in_place]

  Create a header-only scaling-factor table, write a headered copy of a data table, or add a header to a data table in place.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  -v, --verbose : flag
    Run script in verbose mode. Print the header before writing.

  -dr, --dry, --dry_run : flag
    Run script in dry-run mode. Print the header and planned file action without creating or modifying a file.

  -md, --mode : {'siq', 'spike'}
    Workflow mode. Type of header to write: 'siq' (siQ-ChIP normalization) or 'spike' (normalization with a spike-in coefficient) (default: '${mode}').

  -fi, --fil_in : file
    Input file path. Input data table to header. If omitted, '--fil_out' creates a header-only utility table.

  -fo, --fil_out : file
    Output file path. Output file to create. With '--fil_in', writes a headered copy. Without '--fil_in', writes a header-only table.

  --in_place : flag
    Modify '--fil_in' in place instead of writing a separate output file.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - cat
    - head
    - mv

  - The script writes the selected core scaling-factor header.
  - 'execute_calculate_scaling_factor.sh' calls this script with '--fil_in <final_table> --in_place' after successful part-file combination unless '--no_header' is supplied there.
  - If '--fil_in' already begins with the expected header, the copy or in-place operation does not duplicate the header.
  - If '--fil_in' is data-only, the selected header is prepended.
  - If '--fil_in' is omitted, '--fil_out' creates a header-only utility table.
  - If '--dry_run' is enabled, the script validates arguments and reports the planned action without creating or modifying files.

Examples
--------
  1. Create a header-only siQ-ChIP-mode utility table
    '''bash
    bash write_header.sh \\
        --mode siq \\
        --fil_out results/header.siq.tsv
    '''

  2. Add a spike-mode header to a data-only table in place
    '''bash
    bash write_header.sh \\
        --mode spike \\
        --fil_in results/ChIP_samples_spike.tsv \\
        --in_place
    '''

  3. Preview writing a headered siQ-ChIP copy
    '''bash
    bash write_header.sh \\
        --dry_run \\
        --mode siq \\
        --fil_in results/ChIP_samples_siq_6nd.data.tsv \\
        --fil_out results/ChIP_samples_siq_6nd.tsv
    '''
EOM
}
