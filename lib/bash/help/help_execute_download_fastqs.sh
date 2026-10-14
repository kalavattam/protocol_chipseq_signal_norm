#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: help_execute_download_fastqs.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
#   GPT-5.6);
# - Anthropic Claude Code (Opus 5, Opus 5.5).
#
# Distributed under the MIT license.


function help_execute_download_fastqs() {
    # The owning execute wrapper initializes these documented defaults before
    # rendering help; standalone sourcing is outside this helper's contract.
    # shellcheck disable=SC2154
    cat >&2 << EOM
Usage
-----
  execute_download_fastqs.sh
    [--help] [--verbose] [--dry_run]
    [--env_nam <str>] [--threads <int>]
    --fil_in <file> --dir_out <dir> --dir_sym <dir>
    [--use_existing] [--keep_names]
    [--nam_job <str>] [--dir_eo <dir>] [--slurm] [--time <time>]

  Download FASTQ files listed in a TSV file and create symbolic links with custom names specified in the file. Supports single- and paired-end short reads from FTP or HTTPS sources.

  Jobs run serially by default. When '--threads' is greater than 1, jobs run with GNU Parallel, optionally inside one Slurm job.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  -v, --verbose : flag
    Run script in verbose mode: print detailed script state and commands.

  -dr, --dry, --dry_run : flag
    Run script in dry-run mode: print commands without running them.

  -en, --env, --env_nam : str
    Conda environment to activate (default: '${env_nam}').

  -t, --thr, --threads : int
    Number of download jobs to run concurrently, each a single-threaded process; with '--slurm', also the number of CPUs requested (default: ${threads}).

  -fi, --fil_in : file
    Input TSV file containing 'run_accession', 'custom_name', and 'fastq_ftp' or 'fastq_https' columns.

  -do, --dir_out : dir
    Output directory for accession-named downloaded FASTQ files.

  -dy, --dir_sym, --dir_symlink : dir
    Output directory for custom-named FASTQ symlinks.

  -ue, --use_existing : flag
    Check and keep a FASTQ file already in '--dir_out' instead of downloading it again.

  -kn, --keep_names : flag
    Keep the read names as the archive wrote them; do not restore the Illumina names.

  -nj, --nam_job : str
    Job-name prefix used for logs and scheduler submissions (default: '${nam_job}').

  -deo, --dir_eo : dir
    Directory for stderr and stdout log files plus GNU Parallel config files (default: '\${dir_out}/logs').

  -sl, --slurm : flag
    Submit a GNU Parallel download job through Slurm.

  -tm, --time : time
    Slurm job time limit, in 'h:mm:ss' format (default: '${time}'); used only with '--slurm', ignored with a warning otherwise.

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
    - gzip
    - head
    - ln
    - md5sum (when the TSV has 'fastq_md5' values)
    - mv
    - Network access (when '--dry_run' is not specified)
    - parallel (when '--threads' is greater than 1)
    - rm
    - sbatch (when '--slurm' is specified)
    - wc
    - wget

  - The TSV header must include 'run_accession', 'custom_name', and either 'fastq_ftp' or 'fastq_https'.
  - Paired-end FASTQ URLs must be stored in one URL field separated by a semicolon.
  - URLs without a scheme, as in ENA's 'fastq_ftp' values, are downloaded over 'https://'.
  - Optional 'fastq_md5' and 'fastq_bytes' columns, separated by a semicolon like the URLs, give each file's expected MD5 and size, and each download is checked against them.
    + They are recommended, not required: a download can have the right size but corrupt content, which only an MD5 check catches.
    + ENA's filereport returns both when asked, for example 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=<project or run>&result=read_run&format=tsv&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes'; add the 'custom_name' column yourself.
    + Files without an MD5 are checked only with 'gzip -t', and the run warns once.
  - A file that fails its check is downloaded again, up to 3 attempts in all; then it is removed and its entry fails.
  - Read names: ENA moves each read's Illumina name into the header comment ('@<run>.<n> <Illumina name>/<mate>'), and aligners keep only '<run>.<n>', so optical duplicates cannot be told apart. After its check, each entry's headers are rewritten to '@<Illumina name> <run>.<n>'; sequence and quality lines are unchanged. '--keep_names' turns this off.
    + Every header of both mates must have the archive's form, with spot numbers in order, mate suffixes matching the file, and equal names in both mates; otherwise the files stay as downloaded and one warning gives the count.
    + Names the archive did not store ('@<run>.<n> <n>/1') are left as is, with one note.
    + Names that are not '<run>.<n>', as in FASTQ files straight from a sequencer, are already what is wanted; they are left as is, with one note, and only the first header is read.
    + The rewrite takes about 5 minutes per pair of 9 million reads and needs free space about the size of the entry while it runs.
  - With '--use_existing', a file already in '--dir_out' is checked rather than downloaded again. The check uses the file's MD5 and size when they are given, and 'gzip -t' when they are not.
    + A file whose names are already restored no longer matches the archive's MD5, so it gets 'gzip -t' only, and a rerun notes that its names are already restored.
    + A file that fails its check fails its entry and is not downloaded again.
    + A missing file is downloaded as usual.
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
