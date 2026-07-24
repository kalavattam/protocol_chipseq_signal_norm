#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: help_submit_calculate_scaling_factor.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


# TODO: Add more direct-submit examples for siQ-ChIP and spike-in modes.
function help_submit_calculate_scaling_factor() {
    # The submit owner initializes interpolated defaults before invocation.
    # shellcheck disable=SC2154
    cat << EOM >&2
Usage
-----
  submit_calculate_scaling_factor.sh
    [--help] [--env_nam <str>] --dir_scr <dir> [--threads <int>]
    [--mode <mode>] [--method <method>]
    --csv_mip <csv> --csv_min <csv> [--csv_sip <csv>] [--csv_sin <csv>] [--aln_typ <layout>] [--ref_fa <file>] --fil_out <file> [--idx_out <int>]
    [--tbl_met <file>] [--cfg_met <file>] [--eqn <equation>]
    [--len_def <int>] [--len_mip <csv>] [--len_min <csv>] [--dep_mip <csv>] [--dep_min <csv>] [--dep_sip <csv>] [--dep_sin <csv>]
    [--dp <int>] --dir_eo <dir> [--nam_job <str>]


  Compute one per-sample siQ-ChIP or spike-in scaling-factor row from comma-separated BAM/CRAM lists and write it to a deterministic part file.


Parameters
----------
  -h, --help : flag
    Print this help message and exit.

  -en, --env, --env_nam : str
    Conda environment to activate (default: '${env_nam}').

  -ds, --dir_scr : dir
    Directory containing scripts and functions.

  -t, --thr, --threads : int
    Number of threads to use for alignment-processing steps (default: ${threads}).

  -md, --mode : {'siq', 'spike'}
    Workflow mode. Scaling-factor framework: 'siq' or 'spike' (default: '${mode}').

  -me, --method : {'fractional', 'chiprx_alpha_ratio', 'chiprx_alpha_ip', 'chiprx_alpha_in', 'rxinput_alpha'}
    Workflow method. Spike-in coefficient to compute: 'fractional', 'chiprx_alpha_ratio', 'chiprx_alpha_ip', 'chiprx_alpha_in', 'rxinput_alpha', or aliases ('--mode spike'; default: 'chiprx_alpha_ratio').

  -cmip, --csv_mip : list of file
    Comma-separated list of main IP alignment files (BAM or CRAM).

  -cmin, --csv_min : list of file
    Comma-separated list of main input alignment files (BAM or CRAM).

  -csip, --csv_sip : list of file
    Comma-separated list of spike-in IP alignment files (BAM or CRAM; '--mode spike').

  -csin, --csv_sin : list of file
    Comma-separated list of spike-in input alignment files (BAM or CRAM; '--mode spike').

  -at, --aln_typ, --align_typ : {'pe', 'se', 'auto'}
    Alignment layout type for input alignment files: 'pe', 'se', or 'auto' (default: '${aln_typ}').

  -rf, --ref_fa : file
    Reference FASTA file. Reference FASTA required when any input alignment file is CRAM.

  -fo, --fil_out : file
    Output file path. Base output TSV path used to derive '<fil_out>.part.<idx>'.

  -io, --idx_out : int
    Optional output-part index override for execute-layer local dispatch.

  -tb, --tbl_met : file
    siQ-ChIP metadata table. Used in 'siq' mode.

  -cm, --cfg_met : file
    YAML configuration file for metadata parsing. Used in 'siq' mode.

  -eq, --eqn : {'5', '5nd', '6', '6nd'}
    siQ-ChIP alpha equation to compute when '--mode siq' is active: '5', '5nd', '6', or '6nd' (default: '${eqn}').

  -ld, --len_def : int
    Default fragment length for single-end libraries when no per-sample override is provided.

  -lmp, --len_mip : list of int
    Fragment length value(s) for main IP alignment files. Optional comma-separated list of precomputed fragment lengths.

  -lmn, --len_min : list of int
    Fragment length value(s) for main input alignment files. Optional comma-separated list of precomputed fragment lengths.

  -dmp, --dep_mip : list of int
    Sequencing/alignment depth value(s) for main IP alignment files. Optional comma-separated list of precomputed alignment counts.

  -dmn, --dep_min : list of int
    Optional comma-separated list of precomputed alignment counts for main input alignment files.

  -dsp, --dep_sip : list of int
    Sequencing/alignment depth value(s) for spike-in IP alignment files. Optional comma-separated list of precomputed alignment counts.

  -dsn, --dep_sin : list of int
    Sequencing/alignment depth value(s) for spike-in input alignment files. Optional comma-separated list of precomputed alignment counts.

  -dp, --dp : int
    Maximum number of decimal places retained for finite emitted values (default: ${dp}).

  -deo, --dir_eo : dir
    Directory for stderr and stdout log files.

  -nj, --nam_job : str
    Job name. Job-name prefix (default depends on resolved mode/method).


Notes
-----
  Runtime requirements:
    - A compatible Conda environment providing the listed tools
    - bash >= 4.4
    - bc
    - conda (when the requested environment is not active)
    - gawk
    - python >= 3.11
    - Reference FASTA and required index (when processing CRAM)
    - samtools

  - Input alignment paths supplied to this wrapper should not contain spaces, commas, or semicolons.
  - This wrapper expects coordinate-sorted BAM or CRAM files.
  - CRAM inputs require '--ref_fa'.
  - Optional override vectors may be omitted, may contain a single broadcast value, or may contain one value per sample.
  - For '--mode siq' with SE data, '--len_def' or both '--len_mip' and '--len_min' must be supplied.
  - This worker writes one '<fil_out>.part.<idx>' row. The execute-layer wrapper combines successful part files into the final TSV.
  - '--idx_out' is an internal execute-layer coordination option. Omit it for direct multi-sample submit calls and Slurm array tasks.
  - Compute downstream denominator floors separately with 'compute_input_floor'.
  - To run in "debug mode", set hardcoded variable 'debug=true'                   [debug=${debug:-UNSET}]
  - To run in "parse-only mode", set hardcoded variable 'p_only=true'             [p_only=${p_only:-UNSET}]
  - To run in "parse-and-check-only mode", set hardcoded variable 'pc_only=true'  [pc_only=${pc_only:-UNSET}]


Examples
--------
  1. Compute a spike-in ChIP-Rx alpha ratio from paired-end BAM files.
    '''bash
    bash \${HOME}/repos/protocol_chipseq_signal_norm/bin/submit_calculate_scaling_factor.sh \\
        --dir_scr \${HOME}/repos/protocol_chipseq_signal_norm/scripts \\
        --env_nam env_protocol \\
        --threads 4 \\
        --mode spike \\
        --method chiprx_alpha_ratio \\
        --csv_mip \${HOME}/repos/protocol_chipseq_signal_norm/data/processed/align_reads/bowtie2_global_flag-2_mapq-1/sc/IP_WT_Q_Hmo1_7751.sc.bam \\
        --csv_min \${HOME}/repos/protocol_chipseq_signal_norm/data/processed/align_reads/bowtie2_global_flag-2_mapq-1/sc/in_WT_Q_Hmo1_7751.sc.bam \\
        --csv_sip \${HOME}/repos/protocol_chipseq_signal_norm/data/processed/align_reads/bowtie2_global_flag-2_mapq-1/sp/IP_WT_Q_Hmo1_7751.sp.bam \\
        --csv_sin \${HOME}/repos/protocol_chipseq_signal_norm/data/processed/align_reads/bowtie2_global_flag-2_mapq-1/sp/in_WT_Q_Hmo1_7751.sp.bam \\
        --fil_out \${HOME}/repos/protocol_chipseq_signal_norm/data/processed/compute_signal/bowtie2_global_flag-2_mapq-1/tables/test_spike.tsv \\
        --dp 24 \\
        --dir_eo \${HOME}/repos/protocol_chipseq_signal_norm/data/processed/compute_signal/bowtie2_global_flag-2_mapq-1/tables/logs \\
        --nam_job calc_sf_spike
    '''

  2. Compute an Equation 5 siQ-ChIP scaling factor from paired-end BAM files.
    '''bash
    bash \${HOME}/repos/protocol_chipseq_signal_norm/bin/submit_calculate_scaling_factor.sh \\
        --dir_scr \${HOME}/repos/protocol_chipseq_signal_norm/scripts \\
        --env_nam env_protocol \\
        --threads 4 \\
        --mode siq \\
        --csv_mip \${HOME}/project/alignments/IP_sample.sc.bam \\
        --csv_min \${HOME}/project/alignments/in_sample.sc.bam \\
        --aln_typ pe \\
        --fil_out \${HOME}/project/tables/siq_scaling.tsv \\
        --tbl_met \${HOME}/project/metadata/measurements.tsv \\
        --cfg_met \${HOME}/repos/protocol_chipseq_signal_norm/data/raw/docs/parse_metadata_siqchip.yml \\
        --eqn 5 \\
        --dir_eo \${HOME}/project/tables/logs \\
        --nam_job calc_sf_siq
    '''

EOM
}
