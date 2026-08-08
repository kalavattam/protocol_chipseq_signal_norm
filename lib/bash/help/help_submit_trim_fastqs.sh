#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: help_submit_trim_fastqs.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in design, development,
# and documentation, with all output reviewed, edited, and approved by the
# author.
#
# Distributed under the MIT license.


function help_submit_trim_fastqs() {
    # The submit owner initializes interpolated defaults before invocation.
    # shellcheck disable=SC2154
    cat << EOM >&2
Usage
-----
  submit_trim_fastqs.sh
    [--help]
    [--env_nam <str>] --dir_scr <dir>
    [--threads <int>]
    --csv_fil_in <csv> --dir_out <dir>
    --sfx_se <str> --sfx_pe <str>
    --dir_eo <dir> [--nam_job <str>]

  Submit or execute one or more FASTQ-trimming jobs by calling downstream program 'atria'.

  This wrapper
    - parses a semicolon-delimited list of FASTQ input entries,
    - derives sample names,
    - activates the requested Conda environment, and then
    - runs read trimming either under Slurm array execution or by serial/GNU-Parallel-style iteration, depending on how the script is invoked.

  For each input entry, this script writes log files to:

    \${dir_eo}/\${nam_job}.\${samp}.stdout.txt
    \${dir_eo}/\${nam_job}.\${samp}.stderr.txt

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  -en, --env, --env_nam : str
    Conda environment to activate.

  -ds, --dir_scr : dir
    Directory containing scripts and functions.

  -t, --thr, --threads : int
    Number of threads to use.

  -ci, --csv_fil_in : list of structured string
    Comma-separated list of input file paths. Semicolon-delimited serialized string of FASTQ input entries.

    For single-end data, each entry is one FASTQ file. For paired-end data, each entry contains a comma-delimited FASTQ pair.

  -do, --dir_out : dir
    Output directory. Directory for trimmed FASTQ output files.

  -sxs, --sfx_se, --suffix_se : str
    Suffix to strip from single-end FASTQ filenames.

  -sxp, --sfx_pe, --suffix_pe : str
    Suffix to strip from paired-end FASTQ read-1 filenames.

  -deo, --dir_eo : dir
    Directory for stderr and stdout log files.

  -nj, --nam_job : str
    Job name.

Notes
-----
  Runtime requirements:
    - A compatible Conda environment providing the listed tools
    - atria
    - bash >= 4.4
    - conda (when the requested environment is not active)
    - pbzip2
    - pigz
    - Slurm allocation (when run as a Slurm array task)

  - All arguments are required with the following notes and exceptions:
    + '--env_nam' defaults to 'env_nam=${env_nam}' if not specified.
    + '--threads' defaults to 'threads=${threads}' if not specified.
    + '--nam_job' defaults to 'nam_job=${nam_job}' if not specified.

Examples
--------
  1. Trim one single-end FASTQ file during local execution.
    '''bash
    bash "\${dir_scr}/submit_trim_fastqs.sh" \\
        --env_nam "env_protocol" \\
        --dir_scr "\${dir_scr}" \\
        --threads 4 \\
        --csv_fil_in "\${dir_fastq}/sample.fastq.gz" \\
        --dir_out "\${dir_out}" \\
        --sfx_se ".fastq.gz" \\
        --sfx_pe "_R1.fastq.gz" \\
        --dir_eo "\${dir_eo}" \\
        --nam_job "trim_se"
    '''

  2. Trim two paired-end FASTQ entries during local iteration.
    '''bash
    bash "\${dir_scr}/submit_trim_fastqs.sh" \\
        --env_nam "env_protocol" \\
        --dir_scr "\${dir_scr}" \\
        --threads 4 \\
        --csv_fil_in "\${dir_fastq}/sample_1_R1.fastq.gz,\${dir_fastq}/sample_1_R2.fastq.gz;\${dir_fastq}/sample_2_R1.fastq.gz,\${dir_fastq}/sample_2_R2.fastq.gz" \\
        --dir_out "\${dir_out}" \\
        --sfx_se ".fastq.gz" \\
        --sfx_pe "_R1.fastq.gz" \\
        --dir_eo "\${dir_eo}" \\
        --nam_job "trim_pe"
    '''
EOM
}
