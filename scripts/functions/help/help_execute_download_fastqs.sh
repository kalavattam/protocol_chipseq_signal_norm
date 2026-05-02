#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: help_execute_download_fastqs.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
#
# Distributed under the MIT license.


# shellcheck disable=SC2154
function help_execute_download_fastqs() {
    cat << EOM
Usage:
  execute_download_fastqs.sh
    [--help] [--verbose] [--dry_run]
    --threads <int> --infile <str> --dir_out <str> --dir_sym <str>
    --nam_job <str> --err_out <str> --slurm --time <str>

Description:
  Downloads FASTQ files listed in a TSV file and creates symbolic links with custom names specified in the file. It handles single- and paired-end reads and supports downloads from FTP and HTTPS sources.

  The script can execute jobs in parallel using GNU Parallel, optionally with Slurm, or run them serially.

Arguments:
   -h, --help  <flg>
       Display this help message and exit.

   -v, --verbose  <flg>
       Run script in 'verbose' mode.

  -dr, --dry, --dry_run, --dry-run  <flg>
       Print commands without running them.

   -t, --thr, --threads  <int>
       Number of threads to use (default: ${threads}; see 'Notes' below).

   -i, -fi, --infile, --fil_in, --fil-in
       Input TSV file containing download metadata.

  -do, --dir_out, --dir-out  <str>
       Output directory for downloaded FASTQ files

  -dy, --dir_sym, --dir-sym, --dir_symlink, --dir-symlink  <str>
       Output directory for symlinked FASTQ files.

  -nj, --nam_job, --nam-job  <str>
       The name of the job (default: ${nam_job}).

  -eo, --err_out, --err-out  <str>
       The directory to store stderr and stdout TXT outfiles (default: \${dir_out}/logs).

  -sl, --slurm  <flg>
       Submit jobs to the Slurm scheduler (optional; see 'Notes' below).

  -tm, --time  <str>
       The length of time (e.g., h:mm:ss) for the Slurm job (required if '--slurm' is specified, ignored if not; default: ${time}).

Dependencies:
  - Programs
    + Bash >= 4.4
    + cut
    + GNU Parallel (if 'threads > 1')
    + ln
    + Slurm (if '--slurm' is specified)
    + wget
  - Sourced function scripts
    + check_args.sh
    + check_env.sh
    + check_inputs.sh
    + check_numbers.sh
    + format_outputs.sh
    + handle_env.sh
    + wrap_cmd.sh

Notes:
  - The script requires a properly formatted TSV (tab-separated value) file with a header and columns for run accession numbers, custom file names, and URLs (FTP or HTTPS). For paired-end files, URLs in the TSV should be separated by semicolons. See TSV files in 'protocol_chipseq_signal_norm/data/raw/docs' for examples.
  - Symbolic links are created in 'dir_sym' with names specified by the 'custom_name' column in the TSV file.
  - If 'threads' is a positive integer greater than 1, the job submission script will be executed using GNU Parallel with the '--jobs' option set to \${threads}.
  - If the '--slurm' flag is specified and 'threads' is greater than 1, the job submission script will run via GNU Parallel within a Slurm job submission.
  - If '--slurm' is specified and 'threads' is set to 1, the execution script will terminate with an error, as serial job submissions to Slurm are not permitted (and array job submission has not been implemented).

Example:
  '''bash
  bash "\${dir_scr}/execute_download_fastqs.sh"
      --threads "\${threads}"
      --infile "\${pth_tsv}"
      --dir_out "\${dir_raw}"
      --dir_sym "\${dir_sym}"
      --err_out "\${dir_raw}/logs"
      --slurm
  '''
EOM
}
