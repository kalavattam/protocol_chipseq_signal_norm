#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: help_submit_download_fastqs.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6);
# - Anthropic Claude Code (Opus 5, Opus 5.5).
#
# Distributed under the MIT license.


# TODO: do we need '--dir_scr' in the examples? Only if using 'sbatch'.
function help_submit_download_fastqs() {
    cat >&2 << EOM
Usage
-----
  submit_download_fastqs.sh
    [--help] [--dir_scr <dir>]
    --srr <str> --url_1 <str> [--url_2 <str>]
    --dir_out <dir> --dir_sym <dir> --nam_cus <str> --dir_eo <dir> --nam_job <str>

  Download one single-end or paired-end FASTQ entry and create custom symlink(s).

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  -ds, --dir_scr : dir
    Maintained entrypoint directory, the repository 'bin'; shared functions are read from the adjacent 'lib/bash'.

    Passed by the 'execute_*.sh' wrappers, and needed when this script runs from a copy, as 'sbatch <script>' does.

  -sr, --srr : str
    NCBI SRA database run accession code.

  -u1, --url_1 : str
    URL (FTP or HTTPS) for FASTQ file.

  -u2, --url_2 : str
    Second FASTQ URL for paired-end data (default: 'NA', for single-end data).

  -do, --dir_out : dir
    Output directory for FASTQ files.

  -dy, --dir_sym, --dir_symlink : dir
    Directory for symlink(s) to FASTQ file(s).

  -nc, --nam_cus : str
    Custom name for symlink(s).

  -deo, --dir_eo : dir
    Directory for stderr and stdout log files.

  -nj, --nam_job : str
    Job name.

Notes
-----
  Runtime requirements:
    - basename
    - bash >= 4.4
    - cat
    - dirname
    - ln
    - Network access
    - wget

  - Omit '--url_2', or pass 'NA', for single-end data.

Examples
--------
  1. Download one single-end FASTQ and create a custom symlink.
    '''bash
    bash "\${dir_scr}/submit_download_fastqs.sh" \\
        --dir_scr "\${dir_scr}" \\
        --srr SRR_SINGLE \\
        --url_1 "\${url_single}" \\
        --dir_out "\${dir_out}" \\
        --dir_sym "\${dir_sym}" \\
        --nam_cus sample_single \\
        --dir_eo "\${dir_eo}" \\
        --nam_job download_fastqs
    '''

  2. Download one paired-end FASTQ pair and create custom symlinks.
    '''bash
    bash "\${dir_scr}/submit_download_fastqs.sh" \\
        --dir_scr "\${dir_scr}" \\
        --srr SRR_PAIRED \\
        --url_1 "\${url_r1}" \\
        --url_2 "\${url_r2}" \\
        --dir_out "\${dir_out}" \\
        --dir_sym "\${dir_sym}" \\
        --nam_cus sample_paired \\
        --dir_eo "\${dir_eo}" \\
        --nam_job download_fastqs
    '''
EOM
}
