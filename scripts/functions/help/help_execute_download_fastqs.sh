#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: help_execute_download_fastqs.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models) were used in
# development.
#
# Distributed under the MIT license.


# shellcheck disable=SC2154
function help_execute_download_fastqs() {
    cat << EOM
Usage:
  execute_download_fastqs.sh
    [--help] [--verbose] [--dry_run]
    [--threads <int>]
    --infile <file> --dir_out <dir> --dir_sym <dir>
    [--nam_job <str>] [--dir_eo <dir>] [--slurm] [--time <time>]

Description:
  Download FASTQ files listed in a TSV file and create symbolic links with custom names specified in the file. Supports single- and paired-end short reads from FTP or HTTPS sources.

  Jobs run serially by default. When '--threads' is greater than 1, jobs run with GNU Parallel, optionally inside one Slurm job.

Keyword arguments:
  -h, --help  <flag>
    Display this help message and exit.

  -v, --verbose  <flag>
    Print detailed script state and commands.

  -dr, --dry, --dry_run  <flag>
    Print commands without running them.

  -t, --thr, --threads  <int>
    Number of download jobs to run concurrently (default: ${threads}).

  -i, -fi, --infile, --fil_in  <file>
    Input TSV file containing 'run_accession', 'custom_name', and 'fastq_ftp' or 'fastq_https' columns.

  -do, --dir_out  <dir>
    Output directory for accession-named downloaded FASTQ files.

  -dy, --dir_sym, --dir_symlink  <dir>
    Output directory for custom-named FASTQ symlinks.

  -nj, --nam_job  <str>
    Job-name prefix used for logs and scheduler submissions (default: '${nam_job}').

  -eo, --dir_eo  <dir>
    Directory for stderr/stdout logs and GNU Parallel config files (default: '\${dir_out}/logs').

  -sl, --slurm  <flag>
    Submit a GNU Parallel download job through Slurm.

  -tm, --time  <time>
    Slurm walltime in h:mm:ss format; used only with '--slurm' (default: '${time}').

Dependencies:
  Recommended environment:
    - env_protocol

  External programs:
    - Bash >= 4.4
    - cut
    - GNU Parallel, when '--threads' is greater than 1
    - ln
    - Slurm, when '--slurm' is set
    - wget

Notes:
  - The TSV header must include 'run_accession', 'custom_name', and either 'fastq_ftp' or 'fastq_https'.
  - Paired-end FASTQ URLs must be stored in one URL field separated by a semicolon.
  - Downloaded files are written to 'dir_out' using run-accession names. Symlinks are written to 'dir_sym' using 'custom_name' values.
  - Repeated 'run_accession' rows with identical URLs are downloaded once and linked to every requested 'custom_name'.
  - Repeated 'run_accession' rows with conflicting URLs and duplicate 'custom_name' values are rejected.
  - Slurm mode requires '--threads' greater than 1 because serial Slurm submission is not implemented.

Examples:
  1. Download FASTQs and create custom symlinks.
    '''bash
    bash "\${dir_scr}/execute_download_fastqs.sh"
      --threads "\${threads}"
      --infile "\${pth_tsv}"
      --dir_out "\${dir_raw}"
      --dir_sym "\${dir_sym}"
      --dir_eo "\${dir_raw}/logs"
      --slurm
    '''
EOM
}
