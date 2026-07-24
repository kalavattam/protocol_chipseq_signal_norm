#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: help_submit_filter_alignments.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


function help_submit_filter_alignments() {
    # The submit owner initializes interpolated defaults before invocation.
    # shellcheck disable=SC2154
    cat << EOM >&2
Usage
-----
  submit_filter_alignments.sh
    [--help] [--env_nam <str>] --dir_scr <dir> [--threads <int>]
    --csv_fil_in <csv> [--ref_fa <file>]
    --dir_out <dir> [--out_ext <format>]
    [--retain <choice>] [--mito] [--tg] [--mtr] [--chk_chr]
    --dir_eo <dir> [--nam_job <str>]

  Submit or execute one or more alignment-filtering jobs by calling downstream functions 'filter_alignment_sc' or 'filter_alignment_sp'.

  This wrapper
    - parses a comma-delimited list of BAM or CRAM input files,
    - determines which downstream filtering function to run based on '--retain',
    - activates the requested Conda environment, and then
    - runs filtering either under Slurm array execution or by serial or GNU-Parallel iteration, depending on how the script is invoked.

  For each input BAM or CRAM file, this script writes BAM or CRAM output and log files to:

    \${dir_eo}/\${nam_job}.\${samp}.stdout.txt
    \${dir_eo}/\${nam_job}.\${samp}.stderr.txt

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  -en, --env, --env_nam : str
    Conda environment to activate. (default: '${env_nam}').

  -ds, --dir_scr : dir
    Directory containing scripts and functions.

  -t, --thr, --threads : int
    Number of threads to use (default: ${threads}).

  -ci, --csv_fil_in : list of file
    Comma-separated list of input file paths. Comma-delimited list of input BAM or CRAM files.

  -do, --dir_out : dir
    Output directory. Directory in which filtered alignment files will be written.

  -ox, --out_ext : {'bam', 'cram'}
    Final output extension. Filtered output extension: 'bam' or 'cram' (default: '${out_ext}').

  -rt, --retain : {'sc', 'sp'}
    Species chromosomes to retain: 'sc' or 'sp' (default: '${retain}').

  -rf, --ref_fa : file
    Reference FASTA file. Reference FASTA required when any input file is CRAM or '--out_ext cram'.

  -m, --mito : flag
    Retain mitochondrial chromosome.

  -tg, --tg : flag
    Retain SP_II_TG chromosome.

  -mr, --mtr : flag
    Retain SP_MTR chromosome.

  -cc, --chk_chr : flag
    Check chromosomes in output alignment files.

  -deo, --dir_eo : dir
    Directory for stderr and stdout log files.

  -nj, --nam_job : str
    Job name used in log-file naming (default: '${nam_job}').

Notes
-----
  Runtime requirements:
    - A compatible Conda environment providing the listed tools
    - awk
    - bash >= 4.4
    - conda (when the requested environment is not active)
    - grep
    - mv
    - Reference FASTA and required index (when processing CRAM)
    - rm
    - samtools

  - BAM/CRAM input files must be coordinate-sorted.
  - CRAM inputs require '--ref_fa'.
  - Input and output paths supplied to this wrapper interface must not contain spaces, commas, or semicolons.
  - This wrapper does not support '-' for stdin/stdout. Use the underlying Python scripts directly for streaming input/output workflows.
  - Use consistent file ordering in input and output lists.
  - To run in debug mode, set hardcoded variable 'debug=true'.
  - Flags '--tg' and '--mtr' are only meaningful with '--retain sp'; if supplied with '--retain sc', they are ignored with a warning.

Examples
--------
  1. Retain canonical S. cerevisiae chromosomes in BAM output.
    '''bash
    bash "\${dir_scr}/submit_filter_alignments.sh" \\
        --env_nam "env_protocol" \\
        --dir_scr "\${dir_scr}" \\
        --threads 4 \\
        --csv_fil_in "\${dir_aln}/sample.bam" \\
        --dir_out "\${dir_out}" \\
        --out_ext "bam" \\
        --retain "sc" \\
        --mito \\
        --chk_chr \\
        --dir_eo "\${dir_eo}" \\
        --nam_job "filter_sc"
    '''

  2. Retain S. pombe chromosomes from CRAM input and write CRAM output.
    '''bash
    bash "\${dir_scr}/submit_filter_alignments.sh" \\
        --env_nam "env_protocol" \\
        --dir_scr "\${dir_scr}" \\
        --threads 4 \\
        --csv_fil_in "\${dir_aln}/sample.cram" \\
        --ref_fa "\${dir_ref}/genome.fa" \\
        --dir_out "\${dir_out}" \\
        --out_ext "cram" \\
        --retain "sp" \\
        --tg \\
        --mtr \\
        --chk_chr \\
        --dir_eo "\${dir_eo}" \\
        --nam_job "filter_sp"
    '''
EOM
}
