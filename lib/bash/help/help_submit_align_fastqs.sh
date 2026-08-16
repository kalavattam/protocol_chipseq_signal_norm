#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: help_submit_align_fastqs.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6);
# - Anthropic Claude Code (Opus 5).
#
# Distributed under the MIT license.


# TODO: do we need '--dir_scr' in the examples? Only if using 'sbatch'.
# TODO: include Conda, aligner, and Samtools versions.
# shellcheck disable=SC2154
function help_submit_align_fastqs() {
    cat >&2 << EOM
Usage
-----
  submit_align_fastqs.sh
    [--help]
    [--env_nam <str>] [--dir_scr <dir>] [--threads <int>]
    [--aligner <spec>] [--bt2_mode <spec>] [--bwa_alg <spec>] [--mapq <int>] [--req_flg]
    --index <path> --csv_fil_in <csv> [--ref_fa <file>]
    --dir_out <dir> [--out_ext <format>]
    [--qname] --sfx_se <str> --sfx_pe <str>
    --dir_eo <dir> [--nam_job <str>]

  Submit or execute one or more FASTQ-alignment jobs by calling the downstream
  function 'align_fastqs'.

  This wrapper
    - parses a semicolon-delimited list of FASTQ input entries,
    - derives sample names and alignment fil_out paths,
    - activates the requested Conda environment, and then
    - runs alignment either under Slurm array execution or by serial/local iteration, depending on how the script is invoked.

  For each input entry, this script writes log files to:

    \${dir_eo}/\${nam_job}.\${samp}.stdout.txt
    \${dir_eo}/\${nam_job}.\${samp}.stderr.txt

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  -en, --env, --env_nam : str
    Conda environment to activate (default: '${env_nam}').

  -ds, --dir_scr : dir
    Directory containing scripts and functions. Passed by the 'execute_*.sh' wrappers, and needed when this script is run from a copy, as 'sbatch <script>' does, rather than from its real path.

  -t, --thr, --threads : int
    Number of threads to use (default: ${threads}).

  -a, --aln, --aligner : {'bowtie2', 'bwa', 'bwa-mem2'}
    Alignment program to use: 'bowtie2', 'bwa', or 'bwa-mem2' (default: '${aligner}').

  -2m, --bt2_mode : {'local', 'global', 'end-to-end'}
    Bowtie 2 alignment type: 'local', 'global', or 'end-to-end' (default: '${bt2_mode}').

  -ba, --bwa_alg : {'mem', 'aln'}
    BWA algorithm when '--aligner bwa': 'mem' or 'aln' (default: '${bwa_alg}').

  -mq, --mapq : int
    MAPQ threshold for filtering alignment output files (default: ${mapq}).

  -rq, --req_flg : flag
    Require SAM flag bit 2 for properly paired alignments.

  -ix, --index : path
    Path to the aligner index/reference.

  -ci, --csv_fil_in : list of structured string
    Comma-separated list of input file paths. Semicolon-delimited serialized string of FASTQ input entries. For single-end data, each entry is one FASTQ file. For paired-end data, each entry contains a comma-delimited FASTQ pair.

  -rf, --ref_fa : file
    Reference FASTA file. Reference FASTA path required when '--out_ext cram'.

  -do, --dir_out : dir
    Output directory. Directory to write alignment output files.

  -ox, --out_ext : {'bam', 'cram'}
    Final output extension for alignment files: 'bam' or 'cram' (default: '${out_ext}').

  -qn, --qnam, --qname : flag
    Retain queryname-sorted intermediate alignment files.

  -sxs, --sfx_se, --suffix_se : str
    Suffix to strip from single-end FASTQ filenames.

  -sxp, --sfx_pe, --suffix_pe : str
    Suffix to strip from paired-end FASTQ read-1 filenames.

  -deo, --dir_eo : dir
    Directory for stderr and stdout log files.

  -nj, --nam_job : str
    Job name (default: '${nam_job}').

Notes
-----
  Runtime requirements:
    - A compatible Conda environment providing the listed tools
    - bash >= 4.4
    - bowtie2 (when '--aligner bowtie2' is specified)
    - bwa (when '--aligner bwa_aln' or 'bwa_mem' is specified)
    - bwa-mem2 (when '--aligner bwa_mem2' is specified)
    - conda (when the requested environment is not active)
    - samtools
    - Slurm allocation (when run as a Slurm array task)

  - '--ref_fa' is required when '--out_ext cram'.
  - '--req_flg' and '--qname' are optional flags.

Examples
--------
  1. Align one single-end library with Bowtie 2 and write BAM output.
    '''bash
    bash "\${dir_scr}/submit_align_fastqs.sh" \\
        --env_nam "env_protocol" \\
        --dir_scr "\${dir_scr}" \\
        --threads 4 \\
        --aligner "bowtie2" \\
        --bt2_mode "global" \\
        --mapq 1 \\
        --index "\${dir_ref}/bowtie2/genome" \\
        --csv_fil_in "\${dir_fastq}/sample.fastq.gz" \\
        --dir_out "\${dir_out}" \\
        --out_ext "bam" \\
        --sfx_se ".fastq.gz" \\
        --sfx_pe "_R1.fastq.gz" \\
        --dir_eo "\${dir_eo}" \\
        --nam_job "align_bowtie2"
    '''

  2. Align paired-end reads with BWA-MEM2 and write CRAM output.
    '''bash
    bash "\${dir_scr}/submit_align_fastqs.sh" \\
        --env_nam "env_protocol" \\
        --dir_scr "\${dir_scr}" \\
        --threads 8 \\
        --aligner "bwa-mem2" \\
        --mapq 1 \\
        --index "\${dir_ref}/bwa-mem2/genome.fa" \\
        --csv_fil_in "\${dir_fastq}/sample_R1.fastq.gz,\${dir_fastq}/sample_R2.fastq.gz" \\
        --ref_fa "\${dir_ref}/genome.fa" \\
        --dir_out "\${dir_out}" \\
        --out_ext "cram" \\
        --req_flg \\
        --sfx_se ".fastq.gz" \\
        --sfx_pe "_R1.fastq.gz" \\
        --dir_eo "\${dir_eo}" \\
        --nam_job "align_bwa_mem2"
    '''
EOM
}
