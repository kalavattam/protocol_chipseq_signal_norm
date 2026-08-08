#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: help_execute_filter_alignments.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
# GPT-5.6) were used in design, development, and documentation, with all output
# reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


# shellcheck disable=SC2154
function help_execute_filter_alignments() {
    cat << EOM
Usage
-----
  execute_filter_alignments.sh
    [--help] [--verbose] [--dry_run]
    [--env_nam <str>] [--threads <int>]
    --csv_fil_in <csv> [--ref_fa <file>]
    --dir_out <dir> [--out_ext <format>]
    [--retain <choice>] [--mito] [--tg] [--mtr] [--chk_chr]
    [--dir_eo <dir>] [--nam_job <str>] [--max_job <int>] [--slurm] [--time <time>]

  Filter BAM or CRAM input files retaining species-specific chromosomes for S. cerevisiae ("main" or "experimental" alignments) or S. pombe ("spike-in" alignments). Output defaults to BAM and can be set to CRAM with '--out_ext cram'.

  Optional features include retaining mitochondrial (S. cerevisiae and S. pombe: '--mito') and additional chromosomes (S. pombe: '--tg', '--mtr'), and performing checks on chromosomes in filtered output files ('--chk_chr'). CRAM input or output requires an explicit reference FASTA via '--ref_fa'.

  The script supports parallel execution via Slurm or GNU Parallel, or can run serially.

Parameters
----------
  -h, --help : flag
    Print this help message and exit.

  -v, --verbose : flag
    Run script in verbose mode.

  -dr, --dry, --dry_run : flag
    Run script in dry-run mode. Run the command in "check" mode.

  -en, --env, --env_nam : str
    Conda environment to activate (default: '${env_nam}').

  -t, --thr, --threads : int
    Number of threads to use (default: '${threads}').

  -ci, --csv_fil_in : list of file
    Comma-separated list of input file paths. Comma-delimited serialized string of coordinate-sorted BAM or CRAM input files.

  -rf, --ref_fa : file
    Reference FASTA file. Reference FASTA required when '--csv_fil_in' contains CRAM input or '--out_ext cram'.

  -do, --dir_out : dir
    Output directory. The directory to store species-filtered and -reheadered output files.

  -ox, --out_ext : {'bam', 'cram'}
    Final output extension. Filtered output extension: 'bam' or 'cram' (default: '${out_ext}').

  -rt, --retain : {'sc', 'sp'}
    Species chromosomes to retain: S. cerevisiae, "sc"; S. pombe, "sp" (default: '${retain}').

  -m, --mito : flag
    Retain mitochondrial chromosome.

  -tg, --tg : flag
    Retain SP_II_TG chromosome (sp only).

  -mr, --mtr : flag
    Retain SP_MTR chromosome (sp only).

  -cc, --chk_chr : flag
    Check chromosomes in output alignment files.

  -deo, --dir_eo : dir
    Directory for stderr and stdout log files. The directory to store stderr and stdout TXT output files (default: '\${dir_out}/logs').

  -nj, --nam_job : str
    Job name. The name of the job, which is used when writing stderr and stdout TXT files (default: '${nam_job}').

  -mj, --max_job : int
    Maximum number of jobs to run concurrently (default: '${max_job}').
      - If '--slurm' is specified, controls Slurm array tasks.
      - If '--slurm' is not specified:
        + If 'max_job' is greater than 1, jobs run in parallel via GNU Parallel.
        + If 'max_job' is 1, jobs run sequentially (serial mode).

  -sl, --slurm : flag
    Submit jobs to the Slurm scheduler; otherwise, run them in serial.

  -tm, --time : time
    Slurm job time limit. The length of time, in 'h:mm:ss' format, for the Slurm job (required if '--slurm' is specified, ignored if not; default: '${time}').

Notes
-----
  Runtime requirements:
    - A compatible Conda environment providing the listed tools
    - awk
    - bash >= 4.4
    - conda (when the requested environment is not active)
    - grep
    - mv
    - parallel (when '--slurm' is not specified and multiple jobs are run)
    - rm
    - samtools
    - sbatch (when '--slurm' is specified)

  - When the '--slurm' flag is used, jobs are parallelized via Slurm array tasks; otherwise, if multiple jobs are to be run, they are parallelized locally via GNU Parallel; if only one job is to be run, execution is serial.
  - BAM and CRAM input files must be coordinate-sorted.
  - '--out_ext' defaults to 'bam'.
  - CRAM input or CRAM output requires '--ref_fa'.
  - Flag '--mito' applies to either S. cerevisiae or S. pombe data.
  - Flags '--tg' and '--mtr' apply only to S. pombe data; if supplied with '--retain sc', they are ignored with a warning.

Examples
--------
  1. Use Slurm to filter S. cerevisiae alignments.
    '''bash
    retain="sc"
    bash "\${dir_scr}/execute_filter_alignments.sh"
      --verbose
      --threads "\${threads}"
      --csv_fil_in "\${csv_fil_in}"
      --dir_out "\${dir_out}/\${retain}"
      --dir_eo "\${dir_out}/\${retain}/logs"
      --retain "\${retain}"
      --slurm
    '''

  2. Use GNU Parallel to filter S. pombe alignments.
    '''bash
    retain="sp"
    bash "\${dir_scr}/execute_filter_alignments.sh"
      --verbose
      --threads "\${threads}"
      --csv_fil_in "\${csv_fil_in}"
      --dir_out "\${dir_out}/\${retain}"
      --dir_eo "\${dir_out}/\${retain}/logs"
      --retain "\${retain}"
    '''
EOM
}
