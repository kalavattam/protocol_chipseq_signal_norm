#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: help_symlink_files.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.4) was used in development.
#
# Distributed under the MIT license.


function help_symlink_files() {
    cat << EOM
Usage:
  symlink_files.sh
    [--help] [--dry_run] --csv_infile <str,str,...> (--csv_outfile <str,str,...> | --dir_out <str>) [--no_force] [--quiet]

Description:
  Create symbolic links for a comma-separated list of input files.

Keyword arguments:
  -h, --hlp, --help  <flg>
    Print this help message and exit.

  -dr, --dry, --dry_run  <flg>
    Validate inputs, print planned link commands, and exit.

  -i, --csv_infile  <str>
    Comma-separated list (serialized string) of input files.

    Compatibility aliases include '--infile', '--infiles', '--fil_in', and '--csv_infiles'.

  -o, --csv_outfile  <str>
    Comma-separated list (serialized string) of output paths. Must have the same number of entries as '--csv_infile'.

    Compatibility aliases include '--outfile', '--outfiles', '--fil_out', and '--csv_outfiles'.

  -do, --dir_out  <str>
    Output directory. If used, output links are named using the basenames of '--csv_infile'.

  -nf, --no_force  <flg>
    Refuse to replace existing destination symlinks.

  -q, --quiet  <flg>
    Suppress the final stderr summary line.

Notes:
  - Exactly one of '--dir_out' or '--csv_outfile' must be provided.
  - All input files are validated before any links are created.
  - Resolved output paths must be unique; e.g., under '--dir_out', basename collisions are rejected.
  - Existing non-symlink paths at destination are never overwritten.
  - By default, existing destination symlinks are replaced with 'ln -sf'.
  - Paths containing commas are not supported.

#TODO:
  - Write examples.
EOM
}
