#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: help_execute_align_fastqs.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
#
# Distributed under the MIT license.


# shellcheck disable=SC2154
function help_execute_align_fastqs() {
    cat << EOM
Usage:
  execute_align_fastqs.sh
    [--help] [--verbose] [--dry_run]
    [--threads <int>]
    [--aligner <spec>] [--bt2_aln <spec>] [--bwa_alg <spec>] [--mapq <int>] [--req_flg]
    --index <str> --csv_infile <spec> [--ref <file>]
    --dir_out <dir> [--out_ext <spec>] [--qname] [--sfx_se <str>] [--sfx_pe <str>]
    [--err_out <dir>] [--nam_job <str>]
    [--max_job <int>] [--slurm] [--time <time>]

Description:
  Align single- or paired-end short-read FASTQ files using Bowtie 2, BWA, or BWA-MEM2, followed by post-alignment processing with Samtools, including MAPQ-based filtering, sorting, mate fixing (for paired-end alignments), duplicate marking, and indexing.

  Jobs may be run through Slurm ('--slurm'), GNU Parallel, or serial execution, depending on user arguments and the resolved number of jobs.

Keyword arguments:
  -h, --hlp, --help  <flag>
    Display this help message and exit.

  -v, --verbose  <flag>
    Run script in verbose mode.

  -dr, --dry, --dry_run  <flag>
    Perform a dry run without executing commands.

  -t, --thr, --threads  <int>
    Number of threads to use (default: ${threads}).

  -a, --aln, --aligner  <spec>
    Alignment program to use: 'bowtie2', 'bwa', or 'bwa-mem2' (default: ${aligner}).

  -2a, -bn, --bt2_aln  <spec>
    Bowtie 2 alignment type when '--aligner bowtie2': 'local', 'global', or 'end-to-end' (default: ${bt2_aln}; ignored otherwise).

  -ba, --bwa_alg  <spec>
    BWA algorithm when '--aligner bwa': 'mem' or 'aln' (default: ${bwa_alg}; ignored otherwise, and must remain 'mem' when '--aligner bwa-mem2').

  -mq, --mapq  <int>
    MAPQ threshold for filtering alignment outfiles (default: ${mapq}). To disable MAPQ-based filtering, specify 0.

  -rf, --req_flg  <flag>
    Require flag bit 2, signifying that paired-end alignments are properly paired, when filtering alignment outfiles (optional; ignored for single-end data).

  -ix, --index  <str>
    Path to the aligner index/reference. If using Bowtie 2, the path should end with the index stem, e.g., "\${HOME}/path/stem". If using BWA or BWA-MEM2, the path should be the indexed reference FASTA path, e.g., "\${HOME}/path/stem.fa".

  -i, -fi, -ci, --infile, --infiles, --fil_in, --csv_infile, --csv_infiles  <spec>
    Semicolon-delimited serialized string of FASTQ input entries. For single-end data, each entry is one FASTQ file. For paired-end data, each entry contains a comma-delimited FASTQ pair, e.g., "\${HOME}/path/samp_1.fastq.gz;\${HOME}/path/samp_2_R1.fastq.gz,\${HOME}/path/samp_2_R2.fastq.gz;\${HOME}/path/samp_3.fastq.gz".

  -r, --ref  <file>
    Reference FASTA path required when '--out_ext cram' (ignored otherwise).

  -do, --dir_out  <dir>
    Directory in which to write alignment outfiles.

  -ox, --out_ext  <spec>
    Final alignment output extension: 'bam' or 'cram' (default: ${out_ext}).

  -qn, --qname  <flag>
    Retain queryname-sorted intermediate alignment files.

  -sxs, --sfx_se, --suffix_se  <str>
    Suffix to strip from single-end FASTQ files (default: ${sfx_se}).

  -sxp, --sfx_pe, --suffix_pe  <str>
    Suffix to strip from paired-end FASTQ files (default: ${sfx_pe}).

  -eo, --err_out  <dir>
    Directory in which to store stderr and stdout TXT outfiles (default: \${dir_out}/logs).

  -nj, --nam_job  <str>
    Job name used when writing stderr and stdout TXT files (default: ${nam_job}).

  -mj, --max_job  <int>
    Maximum number of jobs to run concurrently (default: ${max_job}).
      - If '--slurm' is specified, controls Slurm array-task concurrency.
      - If '--slurm' is not specified:
        + If 'max_job' is greater than 1, jobs run in parallel via GNU Parallel.
        + If 'max_job' is 1, jobs run in serial.

  -sl, --slurm  <flag>
    Submit jobs to the Slurm scheduler; otherwise, run them through GNU Parallel or in serial, as resolved locally.

  -tm, --time  <time>
    Length of time, in 'h:mm:ss' format, for the Slurm job (default: ${time}; required if '--slurm' is specified, ignored otherwise).

Dependencies:
  Recommended environment:
    - env_protocol

  External programs:
    - Bash >= 4.4
    - Bowtie 2, BWA, or BWA-MEM2, depending on '--aligner'
    - GNU Parallel, when '--slurm' is not specified and multiple jobs are run
    - Samtools
    - Slurm, when '--slurm' is specified

  Shell scripts:
    - submit_align_fastqs.sh

Notes:
  - When the '--slurm' flag is used, jobs are parallelized via Slurm array tasks. Otherwise, if multiple jobs are to be run, they are parallelized locally via GNU Parallel; if only one job is to be run, execution proceeds in serial.
  - If using Bowtie 2, ensure the path to index files ends with the index stem, e.g., "\${HOME}/path/stem" or "\${HOME}/path/sc_sp_proc". If using BWA or BWA-MEM2, the path should include the reference FASTA filename, e.g., "\${HOME}/path/stem.fa" or "\${HOME}/path/sc_sp_proc.fa".
  - Calling the script with '--qname' retains an intermediate queryname-sorted alignment file used during mate fixing for paired-end alignments.
  - Retained queryname-sorted outfiles will share the same path and stem as the final alignment outfile, but with '.qnam' inserted before the final extension (for example, '.qnam.bam' or '.qnam.cram').
  - When '--out_ext cram' is used, '--ref' must also be supplied. Although intermediate work files are processed in BAM format, the final output, and any retained queryname-sorted output, are written as CRAM, which requires a reference FASTA.

Examples:
  1. Run local Bowtie 2 alignment and write BAM files.
    '''bash
    bash "\${dir_scr}/execute_align_fastqs.sh"
      --threads "\${threads}"
      --aligner bowtie2
      --bt2_aln global
      --index "\${idx_bt2}"
      --csv_infile "\${csv_infile}"
      --dir_out "\${dir_out}/bam"
      --out_ext bam
      --err_out "\${dir_out}/logs"
    '''

  2. Run local BWA-MEM2 alignment and write CRAM files.
    '''bash
    bash "\${dir_scr}/execute_align_fastqs.sh"
      --threads "\${threads}"
      --aligner bwa-mem2
      --index "\${idx_bwa_mem2}"
      --csv_infile "\${csv_infile}"
      --ref "\${ref_fa}"
      --dir_out "\${dir_out}/cram"
      --out_ext cram
      --err_out "\${dir_out}/logs"
    '''

  3. Submit a Slurm array job.
    '''bash
    bash "\${dir_scr}/execute_align_fastqs.sh"
      --threads "\${threads}"
      --aligner bowtie2
      --index "\${idx_bt2}"
      --csv_infile "\${csv_infile}"
      --dir_out "\${dir_out}/bam"
      --err_out "\${dir_out}/logs"
      --max_job 12
      --slurm
      --time 1:00:00
    '''
EOM
}
