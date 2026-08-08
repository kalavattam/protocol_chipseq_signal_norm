#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: help_submit_qsort_bam.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in design, development,
# and documentation, with all output reviewed, edited, and approved by the
# author.
#
# Distributed under the MIT license.


function help_submit_qsort_bam() {
    # The submit owner initializes interpolated defaults before invocation.
    # shellcheck disable=SC2154
    cat << EOM
Usage
-----
  submit_qsort_bam.sh
    [--help] [--env_nam <str>] --dir_scr <dir> [--threads <int>]
    --csv_fil_in <csv> [--ref_fa <file>]
    --dir_out <dir> --dir_eo <dir> [--nam_job <str>]

  Sort BAM or CRAM input files by query name.

  When run inside a Slurm array task, this script processes the indexed input selected by 'SLURM_ARRAY_TASK_ID'. Otherwise, it processes every input file serially in the current shell.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  -en, --env, --env_nam : str
    Conda environment to activate (default: '${env_nam}').

  -ds, --dir_scr : dir
    Directory containing scripts and functions.

  -t, --thr, --threads : int
    Number of threads to use per task (default: ${threads}).

  -ci, --csv_fil_in : list of file
    Comma-separated list of input file paths for BAM or CRAM files.

  -rf, --ref_fa : file
    Reference FASTA file. Reference FASTA required when any input file is CRAM.

  -do, --dir_out : dir
    Output directory. Directory to save queryname-sorted alignment files.

  -deo, --dir_eo : dir
    Directory for stderr and stdout log files.

  -nj, --nam_job : str
    Job name (default: '${nam_job}').

Notes
-----
  Runtime requirements:
    - A compatible Conda environment providing the listed tools
    - bash >= 4.4
    - conda (when the requested environment is not active)
    - Reference FASTA and required index (when processing CRAM)
    - samtools

  - BAM inputs are written as '*.qnam.bam'.
  - CRAM inputs are written as '*.qnam.cram' and require '--ref_fa'.
  - This is a submit-wrapper script: it supports serial/local iteration and Slurm-array task execution, but it does not submit Slurm jobs itself.

Examples
--------
  1. Queryname-sort two BAM files during local serial execution.
    '''bash
    bash "\${dir_scr}/submit_qsort_bam.sh" \\
        --env_nam "env_protocol" \\
        --dir_scr "\${dir_scr}" \\
        --threads 4 \\
        --csv_fil_in "\${dir_aln}/sample_1.bam,\${dir_aln}/sample_2.bam" \\
        --dir_out "\${dir_qsort}" \\
        --dir_eo "\${dir_eo}" \\
        --nam_job "qsort_bam"
    '''

  2. Queryname-sort CRAM input with its reference FASTA.
    '''bash
    bash "\${dir_scr}/submit_qsort_bam.sh" \\
        --env_nam "env_protocol" \\
        --dir_scr "\${dir_scr}" \\
        --threads 4 \\
        --csv_fil_in "\${dir_aln}/sample.cram" \\
        --ref_fa "\${dir_ref}/genome.fa" \\
        --dir_out "\${dir_qsort}" \\
        --dir_eo "\${dir_eo}" \\
        --nam_job "qsort_cram"
    '''
EOM
}
