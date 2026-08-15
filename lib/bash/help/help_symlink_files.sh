#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: help_symlink_files.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.4, GPT-5.5, GPT-5.6);
# - Anthropic Claude Code (Opus 5).
#
# Distributed under the MIT license.


function help_symlink_files() {
    cat >&2 << EOM
Usage
-----
  symlink_files.sh
    [--help] [--dry_run] --csv_fil_in <csv> (--csv_fil_out <csv> | --dir_out <dir>) [--no_force] [--quiet]

  Create symbolic links for a comma-separated list of input files.

Parameters
----------
  -h, --help : flag
    Print this help message and exit.

  -dr, --dry, --dry_run : flag
    Run script in dry-run mode. Validate inputs, print planned link commands, and exit.

  -ci, --csv_fil_in : list of file
    Comma-separated list of input file paths as a serialized string.

  -co, --csv_fil_out : list of file
    Comma-separated list of output file paths as a serialized string. Must have the same number of entries as '--csv_fil_in'.

  -do, --dir_out : dir
    Output directory. If used, output links are named using the basenames of '--csv_fil_in'.

  -nf, --no_force : flag
    Refuse to replace existing destination symlinks.

  -q, --quiet : flag
    Suppress the final stderr summary line.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - ln

  - Exactly one of '--dir_out' or '--csv_fil_out' must be provided.
  - All input files are validated before any links are created.
  - Resolved output paths must be unique; e.g., under '--dir_out', basename collisions are rejected.
  - Existing non-symlink paths at destination are never overwritten.
  - By default, existing destination symlinks are replaced with 'ln -sf'.
  - Paths containing commas are not supported.

Examples
--------
  1. Preview links whose destination names come from input basenames.
    '''bash
    bash bin/symlink_files.sh \\
        --dry_run \\
        --csv_fil_in data/IP_rep1.bam,data/IP_rep2.bam \\
        --dir_out links
    '''

  2. Preview explicit destination paths without replacing existing symlinks.
    '''bash
    bash bin/symlink_files.sh \\
        --dry_run \\
        --csv_fil_in data/IP_rep1.bam,data/input_rep1.bam \\
        --csv_fil_out links/chip.bam,links/input.bam \\
        --no_force
    '''
EOM
}
