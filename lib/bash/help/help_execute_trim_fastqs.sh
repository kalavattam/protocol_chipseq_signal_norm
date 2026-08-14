#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: help_execute_trim_fastqs.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
#   GPT-5.6);
# - Anthropic Claude Code (Opus 5).
#
# Distributed under the MIT license.


# shellcheck disable=SC2154
function help_execute_trim_fastqs() {
    cat << EOM
Usage
-----
  execute_trim_fastqs.sh
    [--help] [--verbose] [--dry_run]
    [--env_nam <str>] [--threads <int>] --csv_fil_in <csv> --dir_out <dir>
    [--sfx_se <str>] [--sfx_pe <str>] [--dir_eo <dir>]
    [--nam_job <str>] [--max_job <int>] [--slurm] [--time <time>]

  Perform read-adapter and quality trimming with the program Atria, working with either single- or paired-end FASTQ files.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  -v, --verbose : flag
    Run script in verbose mode.

  -dr, --dry, --dry_run : flag
    Run script in dry-run mode. Run the command in check mode.

  -en, --env, --env_nam : str
    Conda environment to activate (default: '${env_nam}').

  -t, --thr, --threads : int
    Number of threads to use (default: '${threads}').

  -ci, --csv_fil_in : list of structured string
    Comma-separated list of input file paths. Semicolon-delimited serialized string of FASTQ input entries. For single-end data, each entry is one FASTQ file. For paired-end data, each entry contains a comma-delimited FASTQ pair, e.g., "${HOME}/path/samp_1.fastq.gz;${HOME}/path/samp_2_R1.fastq.gz,${HOME}/path/samp_2_R2.fastq.gz;${HOME}/path/samp_3.fastq.gz".

  -do, --dir_out : dir
    Output directory. Directory for Atria-trimmed FASTQ output files.

  -sxs, --sfx_se, --suffix_se : str
    Suffix to strip from single-end FASTQ filenames (default: '${sfx_se}').

  -sxp, --sfx_pe, --suffix_pe : str
    Suffix to strip from paired-end FASTQ read-1 filenames (default: '${sfx_pe}').

  -deo, --dir_eo : dir
    Directory for stderr and stdout log files. The directory to store stderr and stdout TXT output files (default: '\${dir_out}/logs').

  -nj, --nam_job : str
    Job name. The name of the job, which is used when writing stderr and stdout (default: '${nam_job}').

  -mj, --max_job : int
    Maximum number of jobs to run concurrently (default: '${max_job}').
      - If '--slurm' is specified, controls Slurm array tasks.
      - If '--slurm' is not specified:
        + If 'max_job' is greater than 1, jobs run in parallel via GNU Parallel.
        + If 'max_job' is 1, jobs run sequentially (serial mode).

  -sl, --slurm : flag
    Submit jobs to the Slurm scheduler. Without this flag, jobs run in serial.

  -tm, --time : time
    Slurm job time limit. The length of time, in 'h:mm:ss' format, for the Slurm job (required if '--slurm' is specified, ignored if not; default: '${time}').

Notes
-----
  Runtime requirements:
    - A compatible Conda environment providing the listed tools
    - atria
    - bash >= 4.4
    - conda (when the requested environment is not active)
    - parallel (when '--slurm' is not specified and multiple jobs are run)
    - pbzip2
    - pigz
    - sbatch (when '--slurm' is specified)

  - When the '--slurm' flag is used, jobs are parallelized via SLURM array tasks; otherwise, if multiple jobs are to be run, they are parallelized locally via GNU Parallel; if only one job is to be run, execution is serial.
  - Atria is set to not allow read lengths less than 35 bp. It is also set to search for and trim known adapters, among other things. For more details, see the Atria documentation: github.com/cihga39871/Atria/blob/master/docs/2.Atria_trimming_methods_and_usages.md

Examples
--------
  1. Preview mixed single-end and paired-end trimming on the local host.
    '''bash
    bash "\${HOME}/bin/execute_trim_fastqs.sh" \\
        --verbose \\
        --dry_run \\
        --threads 4 \\
        --csv_fil_in "\${HOME}/path/samp_1.fastq.gz;\${HOME}/path/samp_2_R1.fastq.gz,\${HOME}/path/samp_2_R2.fastq.gz" \\
        --dir_out "\${HOME}/path/output/dir" \\
        --sfx_se ".fastq.gz" \\
        --sfx_pe "_R1.fastq.gz"
    '''

  2. Submit paired-end trimming as a Slurm array job.
    '''bash
    bash "\${HOME}/bin/execute_trim_fastqs.sh" \\
        --threads 8 \\
        --csv_fil_in "\${HOME}/path/samp_1_R1.fastq.gz,\${HOME}/path/samp_1_R2.fastq.gz;\${HOME}/path/samp_2_R1.fastq.gz,\${HOME}/path/samp_2_R2.fastq.gz" \\
        --dir_out "\${HOME}/path/output/dir" \\
        --sfx_se ".fastq.gz" \\
        --sfx_pe "_R1.fastq.gz" \\
        --dir_eo "\${HOME}/path/output/logs" \\
        --max_job 2 \\
        --slurm \\
        --time "1:00:00"
    '''
EOM
}
