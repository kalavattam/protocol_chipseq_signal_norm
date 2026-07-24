#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: help_execute_download_fastqs.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
# GPT-5.6) were used in development and documentation.
#
# Distributed under the MIT license.


function help_execute_download_fastqs() {
    #  The owning execute wrapper initializes these documented defaults before
    #+ rendering help; standalone sourcing is outside this helper's contract
    # shellcheck disable=SC2154
    cat << EOM
Usage
-----
  execute_download_fastqs.sh
    [--help] [--verbose] [--dry_run]
    [--env_nam <str>] [--threads <int>]
    --fil_in <file> --dir_out <dir> --dir_sym <dir>
    [--nam_job <str>] [--dir_eo <dir>] [--slurm] [--time <time>]

  Download FASTQ files listed in a TSV file and create symbolic links with custom names specified in the file. Supports single- and paired-end short reads from FTP or HTTPS sources.

  Jobs run serially by default. When '--threads' is greater than 1, jobs run with GNU Parallel, optionally inside one Slurm job.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  -v, --verbose : flag
    Run script in verbose mode. Print detailed script state and commands.

  -dr, --dry, --dry_run : flag
    Run script in dry-run mode. Print commands without running them.

  -en, --env, --env_nam : str
    Conda environment to activate. (default: '${env_nam}').

  -t, --thr, --threads : int
    Number of threads to use. Number of download jobs to run concurrently (default: ${threads}).

  -fi, --fil_in : file
    Input file path. Input TSV file containing 'run_accession', 'custom_name', and 'fastq_ftp' or 'fastq_https' columns.

  -do, --dir_out : dir
    Output directory for accession-named downloaded FASTQ files.

  -dy, --dir_sym, --dir_symlink : dir
    Output directory for custom-named FASTQ symlinks.

  -nj, --nam_job : str
    Job name. Job-name prefix used for logs and scheduler submissions (default: '${nam_job}').

  -deo, --dir_eo : dir
    Directory for stderr and stdout log files plus GNU Parallel config files (default: '\${dir_out}/logs').

  -sl, --slurm : flag
    Submit jobs to the Slurm scheduler. Submit a GNU Parallel download job through Slurm.

  -tm, --time : time
    Slurm job time limit. Slurm walltime in h:mm:ss format; used only with '--slurm' (default: '${time}').

Notes
-----
  Runtime requirements:
    - A compatible Conda environment providing the listed tools
    - awk
    - basename
    - bash >= 4.4
    - cat
    - conda
    - cut
    - dirname
    - grep (when resolving or activating the requested environment)
    - ln
    - Network access (when '--dry_run' is not specified)
    - parallel (when '--threads' is greater than 1)
    - rm
    - sbatch (when '--slurm' is specified)
    - wget

  - The TSV header must include 'run_accession', 'custom_name', and either 'fastq_ftp' or 'fastq_https'.
  - Paired-end FASTQ URLs must be stored in one URL field separated by a semicolon.
  - Downloaded files are written to 'dir_out' using run-accession names. Symlinks are written to 'dir_sym' using 'custom_name' values.
  - Repeated 'run_accession' rows with identical URLs are downloaded once and linked to every requested 'custom_name' during local execution. Slurm mode rejects repeated accessions.
  - Repeated 'run_accession' rows with conflicting URLs and duplicate 'custom_name' values are rejected.
  - Slurm mode requires '--threads' greater than 1 because serial Slurm submission is not implemented.

Examples
--------
  1. Download FASTQs serially on the local host.
    '''bash
    bash "\${dir_scr}/execute_download_fastqs.sh" \\
        --fil_in "\${pth_tsv}" \\
        --dir_out "\${dir_raw}" \\
        --dir_sym "\${dir_sym}"
    '''

  2. Submit concurrent download jobs through Slurm.
    '''bash
    bash "\${dir_scr}/execute_download_fastqs.sh" \\
        --threads 8 \\
        --fil_in "\${pth_tsv}" \\
        --dir_out "\${dir_raw}" \\
        --dir_sym "\${dir_sym}" \\
        --dir_eo "\${dir_raw}/logs" \\
        --slurm
    '''
EOM
}
