#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: help_submit_convert_bam_bed.sh
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


function help_submit_convert_bam_bed() {
    # The submit owner initializes interpolated defaults before invocation.
    # shellcheck disable=SC2154
    cat << EOM
Usage
-----
  submit_convert_bam_bed.sh
    [--help] [--env_nam <str>] --dir_scr <dir> [--threads <int>]
    --csv_fil_in <csv> [--ref_fa <file>]
    [(--pth_scr_py <file> | --use_awk)]
    --dir_out <dir> --dir_eo <dir> [--nam_job <str>]

  Convert BAM or CRAM input files to BED output files.

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

  -pp, --pth_scr_py : file
    Python converter script (default: '\${dir_scr}/compute_signal.py'). Mutually exclusive with '--use_awk'.

  -awk, --use_awk : flag
    Run AWK processing code rather than the Python script. Mutually exclusive with explicit '--pth_scr_py'. Do not use with single-end data.

  -do, --dir_out : dir
    Output directory. Directory to save BED output files.

  -deo, --dir_eo : dir
    Directory for stderr and stdout log files.

  -nj, --nam_job : str
    Job name (default: '${nam_job}').

Notes
-----
  Runtime requirements:
    - A compatible Conda environment providing the listed tools
    - awk (when '--use_awk' is specified)
    - bash >= 4.4
    - conda (when the requested environment is not active)
    - gzip
    - python >= 3.11 (when '--use_awk' is not specified)
    - Reference FASTA and required index (when processing CRAM)
    - samtools
    - sort (when '--use_awk' is specified)

  - The AWK and Python branches are not equivalent:
    + The AWK branch assumes QNAME-sorted paired-end records in adjacent pairs and writes paired-fragment intervals.
    + The Python branch uses 'compute_signal.py', which does not require QNAME sorting, handles single-end data, and emits fragments according to its own 'parse_bam()' policy.
  - Python conversion is the default. Use '--pth_scr_py' only to supply a custom converter, or use '--use_awk' to select the AWK branch.
  - CRAM inputs require '--ref_fa' in both AWK and Python branches.
  - This is a submit-wrapper script: it supports serial/local iteration and Slurm-array task execution, but it does not submit Slurm jobs itself.

Examples
--------
  1. Convert coordinate-sorted BAM and CRAM files with the Python converter.
    '''bash
    bash "\${dir_scr}/submit_convert_bam_bed.sh" \\
        --env_nam "env_protocol" \\
        --dir_scr "\${dir_scr}" \\
        --threads 4 \\
        --csv_fil_in "\${dir_aln}/sample.bam,\${dir_aln}/sample.cram" \\
        --ref_fa "\${dir_ref}/genome.fa" \\
        --pth_scr_py "\${dir_scr}/compute_signal.py" \\
        --dir_out "\${dir_bed}" \\
        --dir_eo "\${dir_eo}" \\
        --nam_job "convert_python"
    '''

  2. Convert a queryname-sorted paired-end BAM with the AWK branch.
    '''bash
    bash "\${dir_scr}/submit_convert_bam_bed.sh" \\
        --env_nam "env_protocol" \\
        --dir_scr "\${dir_scr}" \\
        --threads 4 \\
        --csv_fil_in "\${dir_aln}/sample.qnam.bam" \\
        --use_awk \\
        --dir_out "\${dir_bed}" \\
        --dir_eo "\${dir_eo}" \\
        --nam_job "convert_awk"
    '''
EOM
}
