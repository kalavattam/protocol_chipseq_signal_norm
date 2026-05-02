#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: help_execute_trim_fastqs.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
#
# Distributed under the MIT license.


# shellcheck disable=SC2154
function help_execute_trim_fastqs() {
    cat << EOM
Usage:
  execute_trim_fastqs.sh
    [--help] [--verbose] [--dry_run]
    --threads <int> --csv_infile <str,str,...> --dir_out <str>
    --sfx_se <str> --sfx_pe <str> --err_out <str> --nam_job <str>
    --max_job <int> [--slurm] [--time <str>]

Description:
  Perform read-adapter and quality trimming with the program Atria, working with either single- or paired-end FASTQ files.

Arguments:
  -h, --help  <flg>
    Display this help message and exit.

  -v, --verbose  <flg>
    Run script in 'verbose' mode .

  -dr, --dry_run  <flg>
    Run the command in check mode.

  -t, --thr, --threads  <int>
    Number of threads to use (default: '${threads}').

  -i, --csv_infile
    Semicolon-delimited serialized string of FASTQ input entries. For single-end data, each entry is one FASTQ file. For paired-end data, each entry contains a comma-delimited FASTQ pair, e.g., "${HOME}/path/samp_1.fastq.gz;${HOME}/path/samp_2_R1.fastq.gz,${HOME}/path/samp_2_R2.fastq.gz;${HOME}/path/samp_3.fastq.gz".

    Compatibility aliases include '--infile', '--infiles', '--fil_in', and '--csv_infiles'.

  -do, --dir_out
    Directory for Atria-trimmed FASTQ outfile(s).

  -sxs, --sfx_se, --suffix_se
    Suffix to strip from single-end sequenced FASTQ files (default: '${sfx_se}').

  -sxp, --sfx_pe, --suffix_pe
    Suffix to strip from the first of two paired-end sequenced FASTQ files (default: '${sfx_pe}').

  -eo, --err_out
    The directory to store stderr and stdout TXT outfiles (default: '\${dir_out}/err_out').

  -nj, --nam_job
    The name of the job, which is used when writing stderr and stdout (default: '${nam_job}').

  -mj, --max_job
    Maximum number of jobs to run concurrently (default: '${max_job}').
      - If '--slurm' is specified, controls Slurm array tasks.
      - If '--slurm' is not specified:
        + If 'max_job' is greater than 1, jobs run in parallel via GNU Parallel.
        + If 'max_job' is 1, jobs run sequentially (serial mode).

  -sl, --slurm  <flg>
    Submit jobs to the Slurm scheduler; otherwise, run them in serial (optional).

  -tm, --time
    The length of time, in 'h:mm:ss' format, for the Slurm job (required if '--slurm' is specified, ignored if not; default: '${time}').

Dependencies:
  External programs:
    - Atria
    - Bash
    - GNU Parallel (when '--slurm' is not specified and multiple jobs are run)
    - pbzip2
    - pigz
    - Slurm (when '--slurm' is specified)

  Sourced function scripts:
    - check_args.sh
        + check_arg_supplied
    - check_env.sh
        + check_env_installed
        + check_pgrm_path
    - check_inputs.sh
        + check_file_dir_exists
        + debug_arr_contents (used by print_parallel_info)
    - check_numbers.sh
        + check_format_time
        + check_int_pos
    - format_outputs.sh
        + echo_error
        + echo_warning
    - handle_env.sh
        + handle_env
    - help/help_execute_trim_fastqs.sh
        + help_execute_trim_fastqs
    - manage_parallel.sh
        + print_parallel_info
        + reset_max_job
        + set_params_parallel
    - process_sequences.sh
        + check_string_fastqs

Notes:
  - When the '--slurm' flag is used, jobs are parallelized via SLURM array tasks; otherwise, if multiple jobs are to be run, they are parallelized locally via GNU Parallel; if only one job is to be run, execution is serial.
  - Atria is set to not allow read lengths less than 35 bp. It is also set to search for and trim known adapters, among other things. For more details, see the Atria documentation: github.com/cihga39871/Atria/blob/master/docs/2.Atria_trimming_methods_and_usages.md

Example:
  '''bash
  bash execute_trim_fastqs.sh
      --verbose
      --dry_run
      --threads 4
      --csv_infile "\${HOME}/path/samp_1.fastq.gz;\${HOME}/path/samp_2_R1.fastq.gz,\${HOME}/path/samp_2_R2.fastq.gz"
      --dir_out "\${HOME}/path/output/dir"
      --sfx_se ".fastq.gz"
      --sfx_pe "_R1.fastq.gz"
      --slurm
  '''
EOM
}
