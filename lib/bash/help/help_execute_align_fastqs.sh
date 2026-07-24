#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: help_execute_align_fastqs.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
# GPT-5.6) were used in development and documentation.
#
# Distributed under the MIT license.


# shellcheck disable=SC2154
function help_execute_align_fastqs() {
    cat << EOM
Usage
-----
  execute_align_fastqs.sh
    [--help] [--verbose] [--dry_run]
    [--env_nam <str>] [--threads <int>]
    [--aligner <spec>] [--bt2_mode <spec>] [--bwa_alg <spec>] [--mapq <int>] [--req_flg]
    --index <path> --csv_fil_in <csv> [--ref_fa <file>]
    --dir_out <dir> [--out_ext <format>] [--qname] [--sfx_se <str>] [--sfx_pe <str>]
    [--dir_eo <dir>] [--nam_job <str>]
    [--max_job <int>] [--slurm] [--time <time>]

  Align single- or paired-end short-read FASTQ files using Bowtie 2, BWA, or BWA-MEM2, followed by post-alignment processing with Samtools, including MAPQ-based filtering, sorting, mate fixing (for paired-end alignments), duplicate marking, and indexing.

  Jobs may be run through Slurm ('--slurm'), GNU Parallel, or serial execution, depending on user arguments and the resolved number of jobs.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  -v, --verbose : flag
    Run script in verbose mode.

  -dr, --dry, --dry_run : flag
    Run script in dry-run mode. Do not execute commands.

  -en, --env, --env_nam : str
    Conda environment to activate. (default: '${env_nam}').

  -t, --thr, --threads : int
    Number of threads to use (default: ${threads}).

  -a, --aln, --aligner : {'bowtie2', 'bwa', 'bwa-mem2'}
    Alignment program to use: 'bowtie2', 'bwa', or 'bwa-mem2' (default: '${aligner}').

  -2m, --bt2_mode : {'local', 'global', 'end-to-end'}
    Bowtie 2 alignment type when '--aligner bowtie2': 'local', 'global', or 'end-to-end' (default: ${bt2_mode}; ignored otherwise).

  -ba, --bwa_alg : {'mem', 'aln'}
    BWA algorithm when '--aligner bwa': 'mem' or 'aln' (default: ${bwa_alg}; ignored otherwise, and must remain 'mem' when '--aligner bwa-mem2').

  -mq, --mapq : int
    MAPQ threshold for filtering alignment output files (default: ${mapq}). To disable MAPQ-based filtering, specify 0.

  -rq, --req_flg : flag
    Require SAM flag bit 2 for properly paired alignments. Ignored for single-end data.

  -ix, --index : path
    Path to the aligner index/reference. If using Bowtie 2, the path should end with the index stem, e.g., "\${HOME}/path/stem". If using BWA or BWA-MEM2, the path should be the indexed reference FASTA path, e.g., "\${HOME}/path/stem.fa".

  -ci, --csv_fil_in : list of structured string
    Comma-separated list of input file paths. Semicolon-delimited serialized string of FASTQ input entries. For single-end data, each entry is one FASTQ file. For paired-end data, each entry contains a comma-delimited FASTQ pair, e.g., "\${HOME}/path/samp_1.fastq.gz;\${HOME}/path/samp_2_R1.fastq.gz,\${HOME}/path/samp_2_R2.fastq.gz;\${HOME}/path/samp_3.fastq.gz".

  -rf, --ref_fa : file
    Reference FASTA file. Reference FASTA path required when '--out_ext cram' (ignored otherwise).

  -do, --dir_out : dir
    Output directory. Directory in which to write alignment output files.

  -ox, --out_ext : {'bam', 'cram'}
    Final output extension for alignment files: 'bam' or 'cram' (default: '${out_ext}').

  -qn, --qname : flag
    Retain queryname-sorted intermediate alignment files.

  -sxs, --sfx_se, --suffix_se : str
    Suffix to strip from single-end FASTQ filenames (default: '${sfx_se}').

  -sxp, --sfx_pe, --suffix_pe : str
    Suffix to strip from paired-end FASTQ read-1 filenames (default: '${sfx_pe}').

  -deo, --dir_eo : dir
    Directory for stderr and stdout log files in TXT format (default: \${dir_out}/logs).

  -nj, --nam_job : str
    Job name used when writing stderr and stdout TXT files (default: '${nam_job}').

  -mj, --max_job : int
    Maximum number of jobs to run concurrently (default: ${max_job}).
      - If '--slurm' is specified, controls Slurm array-task concurrency.
      - If '--slurm' is not specified:
        + If 'max_job' is greater than 1, jobs run in parallel via GNU Parallel.
        + If 'max_job' is 1, jobs run in serial.

  -sl, --slurm : flag
    Submit jobs to the Slurm scheduler; otherwise, run them through GNU Parallel or in serial, as resolved locally.

  -tm, --time : time
    Slurm job time limit. Length of time, in 'h:mm:ss' format, for the Slurm job (default: ${time}; required if '--slurm' is specified, ignored otherwise).

Notes
-----
  Runtime requirements:
    - A compatible Conda environment providing the listed tools
    - bash >= 4.4
    - bowtie2 (when '--aligner bowtie2' is specified)
    - bwa (when '--aligner bwa_aln' or 'bwa_mem' is specified)
    - bwa-mem2 (when '--aligner bwa_mem2' is specified)
    - conda (when the requested environment is not active)
    - parallel (when '--slurm' is not specified and multiple jobs are run)
    - samtools
    - sbatch (when '--slurm' is specified)

  - When the '--slurm' flag is used, jobs are parallelized via Slurm array tasks. Otherwise, if multiple jobs are to be run, they are parallelized locally via GNU Parallel; if only one job is to be run, execution proceeds in serial.
  - If using Bowtie 2, ensure the path to index files ends with the index stem, e.g., "\${HOME}/path/stem" or "\${HOME}/path/sc_sp_proc". If using BWA or BWA-MEM2, the path should include the reference FASTA filename, e.g., "\${HOME}/path/stem.fa" or "\${HOME}/path/sc_sp_proc.fa".
  - Calling the script with '--qname' retains an intermediate queryname-sorted alignment file used during mate fixing for paired-end alignments.
  - Retained queryname-sorted output files will share the same path and stem as the final alignment output file, but with '.qnam' inserted before the final extension (for example, '.qnam.bam' or '.qnam.cram').
  - When '--out_ext cram' is used, '--ref_fa' must also be supplied. Although intermediate work files are processed in BAM format, the final output, and any retained queryname-sorted output, are written as CRAM, which requires a reference FASTA.

Examples
--------
  1. Run local Bowtie 2 alignment and write BAM files.
    '''bash
    bash "\${dir_scr}/execute_align_fastqs.sh" \\
        --threads "\${threads}" \\
        --aligner bowtie2 \\
        --bt2_mode global \\
        --index "\${idx_bt2}" \\
        --csv_fil_in "\${csv_fil_in}" \\
        --dir_out "\${dir_out}/bam" \\
        --out_ext bam \\
        --dir_eo "\${dir_out}/logs"
    '''

  2. Run local BWA-MEM2 alignment and write CRAM files.
    '''bash
    bash "\${dir_scr}/execute_align_fastqs.sh" \\
        --threads "\${threads}" \\
        --aligner bwa-mem2 \\
        --index "\${idx_bwa_mem2}" \\
        --csv_fil_in "\${csv_fil_in}" \\
        --ref_fa "\${ref_fa}" \\
        --dir_out "\${dir_out}/cram" \\
        --out_ext cram \\
        --dir_eo "\${dir_out}/logs"
    '''

  3. Submit a Slurm array job.
    '''bash
    bash "\${dir_scr}/execute_align_fastqs.sh" \\
        --threads "\${threads}" \\
        --aligner bowtie2 \\
        --index "\${idx_bt2}" \\
        --csv_fil_in "\${csv_fil_in}" \\
        --dir_out "\${dir_out}/bam" \\
        --dir_eo "\${dir_out}/logs" \\
        --max_job 12 \\
        --slurm \\
        --time 1:00:00
    '''
EOM
}
