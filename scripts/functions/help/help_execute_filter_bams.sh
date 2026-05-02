#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: help_execute_filter_bams.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
#
# Distributed under the MIT license.


# shellcheck disable=SC2154
function help_execute_filter_bams() {
    cat << EOM
Usage:
  execute_filter_bams.sh
    [--help] [--verbose] [--dry_run]
    --threads <int> --csv_infile <str> --dir_out <str>
    --retain <str> [--mito] [--tg] [--mtr] [--chk_chr] --err_out <str>
    --nam_job <str> --max_job <int> [--slurm] [--time <str>]

Description:
  Filter BAM infiles to retain species-specific chromosomes for S. cerevisiae ("main" or "experimental" alignments) or S. pombe ("spike-in" alignments).

  Optional features include retaining mitochondrial (S. cerevisiae and S. pombe: '--mito') and additional chromosomes (S. pombe: '--tg', '--mtr'), and performing checks on chromosomes in filtered BAM outfiles ('--chk_chr').

  The script supports parallel execution via Slurm or GNU Parallel, or can run serially.

Arguments:
   -h, --help
    Print this help message and exit.

   -v, --verbose
    Run script in "verbose" mode.

  -dr, --dry_run
    Run the command in "check" mode.

   -t, --threads
    Number of threads to use (default: '${threads}').

   -i, --csv_infile
    Comma-delimited serialized string of coordinate-sorted BAM input files.

    Compatibility aliases include '--infile', '--infiles', '--fil_in', and '--csv_infiles'.

  -do, --dir_out
    The directory to store species-filtered and -reheadered BAM output files.

  -rt, --retain
    Specify species chromosomes to retain: S. cerevisiae, "sc"; S. pombe, "sp" (default: '${retain}').

   -m, --mito
    Retain mitochondrial chromosome.

  -tg, --tg
    Retain SP_II_TG chromosome (sp only).

  -mr, --mtr
    Retain SP_MTR chromosome (sp only).

  -cc, --chk_chr
    Check chromosomes in filtered BAM outfile (optional)

  -eo, --err_out
    The directory to store stderr and stdout TXT outfiles (default: '\${dir_out}/err_out').

  -nj, --nam_job
    The name of the job, which is used when writing stderr and stdout TXT files (default: '${nam_job}').

  -mj, --max_job
    Maximum number of jobs to run concurrently (default: '${max_job}').
      - If '--slurm' is specified, controls Slurm array tasks.
      - If '--slurm' is not specified:
        + If 'max_job' is greater than 1, jobs run in parallel via GNU Parallel.
        + If 'max_job' is 1, jobs run sequentially (serial mode).

  -sl, --slurm
    Submit jobs to the Slurm scheduler; otherwise, run them in serial.

  -tm, --time
    The length of time, in 'h:mm:ss' format, for the Slurm job (required if '--slurm' is specified, ignored if not; default: '${time}').

Dependencies:
  External programs:
    - AWK
    - Bash
    - GNU Parallel (when '--slurm' is not specified and multiple jobs are run)
    - grep
    - Samtools
    - Slurm (when '--slurm' is specified)

  Sourced function scripts:
    - check_args.sh
        + check_arg_supplied
        + check_str_delim
    - check_env.sh
        + check_env_installed
        + check_pgrm_path
    - check_inputs.sh
        + check_file_dir_exists
        + debug_arr_contents (used by print_parallel_info)
    - check_numbers.sh
        + check_format_time
        + check_int_pos
    - filter_bam.sh
        + filter_bam_sc
        + filter_bam_sp
    - format_outputs.sh
        + echo_error
        + echo_warning
    - handle_env.sh
        + handle_env
    - help/help_execute_filter_bams.sh
        + help_execute_filter_bams
    - manage_parallel.sh
        + print_parallel_info
        + reset_max_job
        + set_params_parallel

Notes:
  - When the '--slurm' flag is used, jobs are parallelized via Slurm array tasks; otherwise, if multiple jobs are to be run, they are parallelized locally via GNU Parallel; if only one job is to be run, execution is serial.
  - BAM infiles must be coordinate-sorted.
  - Flag '--mito' applies to either S. cerevisiae or S. pombe data.
  - Flags '--tg' and '--mtr' apply only to S. pombe data; if supplied with '--retain sc', they are ignored with a warning.

Examples:
  1. Use Slurm to filter BAM files for S. cerevisiae ("sc") chromosomes (i.e., "main" alignments)
  '''bash
  retain="sc"
  bash "\${dir_scr}/execute_filter_bams.sh"
      --verbose
      --threads "\${threads}"
      --csv_infile "\${csv_infile}"
      --dir_out "\${dir_out}/\${retain}"
      --err_out "\${dir_out}/\${retain}/logs"
      --retain  "\${retain}"
      --slurm
  '''

  2. Use GNU Parallel to filter BAM files for S. pombe ("sp") chromosomes (i.e., "spike-in" alignments)
  '''bash
  retain="sp"
  bash "\${dir_scr}/execute_filter_bams.sh"
      --verbose
      --threads "\${threads}"
      --csv_infile "\${csv_infile}"
      --dir_out "\${dir_out}/\${retain}"
      --err_out "\${dir_out}/\${retain}/logs"
      --retain  "\${retain}"
  '''
EOM
}
