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
    [--threads <int>] [--aligner <str>] [--bt2_aln <str>] [--bwa_alg <str>] [--ref <str>] [--mapq <int>] [--req_flg] --index <str>
    --csv_infile <str> --dir_out <str> [--out_ext <str>] [--qname] [--sfx_se <str>] [--sfx_pe <str>]
    [--err_out <str>] [--nam_job <str>] [--max_job <int>] [--slurm] [--time <str>]

Description:
  Align single- or paired-end short-read FASTQ files using Bowtie 2, BWA, or BWA-MEM2, followed by post-alignment processing with Samtools, including MAPQ-based filtering, sorting, mate fixing (for paired-end alignments), duplicate marking, and indexing.

  Jobs may be run through Slurm ('--slurm'), GNU Parallel, or serial execution, depending on user arguments and the resolved number of jobs.

Keyword arguments, flags:
  -h, --hlp, --help
    Display this help message and exit (0).

  -v, --verbose
    Run script in verbose mode.

  -dr, --dry, --dry_run
    Perform a dry run without executing commands.

  -t, --thr, --threads
    Number of threads to use (default: ${threads}).

  -a, --aligner
    Alignment program to use: 'bowtie2', 'bwa', or 'bwa-mem2' (default: ${aligner}).

  -2a, --bt2_aln
    Bowtie 2 alignment type when '--aligner bowtie2': 'local', 'global', or 'end-to-end' (default: ${bt2_aln}; ignored otherwise).

  -ba, --bwa_alg
    BWA algorithm when '--aligner bwa': 'mem' or 'aln' (default: ${bwa_alg}; ignored otherwise, and must remain 'mem' when '--aligner bwa-mem2').

  -r, --ref
    Reference FASTA path required when '--out_ext cram' (ignored otherwise).

  -mq, --mapq
    MAPQ threshold for filtering alignment outfiles (default: ${mapq}). To disable MAPQ-based filtering, specify 0.

  -rf, --req_flg
    Require flag bit 2, signifying that paired-end alignments are properly paired, when filtering alignment outfiles (optional; ignored for single-end data).

  -ix, --index
    Path to the aligner index/reference. If using Bowtie 2, the path should end with the index stem, e.g., "\${HOME}/path/stem". If using BWA or BWA-MEM2, the path should be the indexed reference FASTA path, e.g., "\${HOME}/path/stem.fa".

  -i, --csv_infile
    Semicolon-delimited serialized string of FASTQ input entries. For single-end data, each entry is one FASTQ file. For paired-end data, each entry contains a comma-delimited FASTQ pair, e.g., "${HOME}/path/samp_1.fastq.gz;${HOME}/path/samp_2_R1.fastq.gz,${HOME}/path/samp_2_R2.fastq.gz;${HOME}/path/samp_3.fastq.gz".

    Compatibility aliases include '--infile', '--infiles', '--fil_in', and '--csv_infiles'.

  -do, --dir_out
    Directory in which to write alignment outfiles.

  -ox, --out_ext
    Final alignment output extension: 'bam' or 'cram' (default: ${out_ext}).

  -qn, --qname
    Retain queryname-sorted intermediate alignment files.

  -sxs, --sfx_se, --sfx-se, --suffix_se, --suffix-se
    Suffix to strip from single-end FASTQ files (default: ${sfx_se}).

  -sxp, --sfx_pe, --sfx-pe, --suffix_pe, --suffix-pe
    Suffix to strip from paired-end FASTQ files (default: ${sfx_pe}).

  -eo, --err_out
    Directory in which to store stderr and stdout TXT outfiles (default: \${dir_out}/logs).

  -nj, --nam_job
    Job name used when writing stderr and stdout TXT files (default: ${nam_job}).

  -mj, --max_job
    Maximum number of jobs to run concurrently (default: ${max_job}).
      - If '--slurm' is specified, controls Slurm array-task concurrency.
      - If '--slurm' is not specified:
        + If 'max_job' is greater than 1, jobs run in parallel via GNU Parallel.
        + If 'max_job' is 1, jobs run in serial.

  -sl, --slurm
    Submit jobs to the Slurm scheduler; otherwise, run them through GNU Parallel or in serial, as resolved locally.

  -tm, --time
    Length of time, in 'h:mm:ss' format, for the Slurm job (default: ${time}; required if '--slurm' is specified, ignored otherwise).


Dependencies:
  External programs:
    - AWK >= 5
    - Bash >= 4.4
    - Bowtie 2, BWA, or BWA-MEM2, depending on '--aligner'
    - GNU Parallel, when '--slurm' is not specified and multiple jobs are run
    - Samtools
    - Slurm, when '--slurm' is specified

  Sourced function scripts:
    - check_args.sh
        + require_optarg
    - check_env.sh
        + check_env_installed
        + check_pgrm_path
    - check_inputs.sh
        + validate_var
        + validate_var_dir
        + validate_var_file
    - check_numbers.sh
        + check_format_time
        + check_int_nonneg
        + check_int_pos
    - format_outputs.sh
        + echo_err
        + print_banner_pretty
    - handle_env.sh
        + handle_env
    - help/help_execute_align_fastqs.sh
        + help_execute_align_fastqs
    - manage_parallel.sh
        + print_parallel_info
        + reset_max_job
        + set_params_parallel
    - process_sequences.sh
        + check_string_fastqs
    - wrap_cmd.sh
        + get_submit_logs
        + print_built_cmd

Notes:
  - When the '--slurm' flag is used, jobs are parallelized via Slurm array tasks. Otherwise, if multiple jobs are to be run, they are parallelized locally via GNU Parallel; if only one job is to be run, execution proceeds in serial.
  - If using Bowtie 2, ensure the path to index files ends with the index stem, e.g., "\${HOME}/path/stem" or "\${HOME}/path/sc_sp_proc". If using BWA or BWA-MEM2, the path should include the reference FASTA filename, e.g., "\${HOME}/path/stem.fa" or "\${HOME}/path/sc_sp_proc.fa".
  - Calling the script with '--qname' retains an intermediate queryname-sorted alignment file used during mate fixing for paired-end alignments.
  - Retained queryname-sorted outfiles will share the same path and stem as the final alignment outfile, but with '.qnam' inserted before the final extension (for example, '.qnam.bam' or '.qnam.cram').
  - When '--out_ext cram' is used, '--ref' must also be supplied. Although intermediate work files are processed in BAM format, the final output, and any retained queryname-sorted output, are written as CRAM, which requires a reference FASTA.

Example:
  '''bash
  bash "\${dir_scr}/execute_align_fastqs.sh"
      --verbose
      --threads "\${threads}"
      --aligner "\${aligner}"
      --bt2_aln "\${bt2_aln}"
      --bwa_alg "\${bwa_alg}"
      --ref "\${ref_fa}"
      --mapq "\${mapq}"
      --req_flg
      --index "\${pth_idx}"
      --csv_infile "\${csv_infile}"
      --dir_out "\${dir_out}/init"
      --out_ext "cram"
      --err_out "\${dir_out}/init/logs"
      --slurm
  '''

#TODO:
  - More examples.
EOM
}
