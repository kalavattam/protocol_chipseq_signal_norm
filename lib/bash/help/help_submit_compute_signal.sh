#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: help_submit_compute_signal.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


function help_submit_compute_signal() {
    # The submit owner initializes interpolated defaults before invocation.
    # shellcheck disable=SC2154
    cat << EOM >&2
Usage
-----
  submit_compute_signal.sh
    [--help]
    [--env_nam <str>] --dir_scr <dir> [--threads <int>]
    [--mode <mode>] [--method <method>]
    (--csv_fil_in <csv> [--ref_fa <file>] [--chr_sizes <file>] | --csv_fil_A <csv> --csv_fil_B <csv> [--chr_sizes <file>])
    --csv_fil_out <csv>
    [--siz_bin <int>] [--chunk_size <int>] [--engine <engine>] [--csv_scl_fct <csv>] [--csv_usr_frg <csv>]
    [--csv_dep_min <csv>] [--csv_pseudo <csv>] [--eps <flt>] [--skip_00 <choice>] [--strict_bins] [--drp_nan]
    [--skp_pfx <csv>] [--track] [--dp <int>]
    --dir_eo <dir> [--nam_job <str>]

  Submit per-sample signal, ratio, or fragment-coordinate jobs from comma-separated file lists to 'compute_signal.py' or 'compute_signal_ratio.py'.

Parameters
----------
  -h, --help : flag
    Print this help message and exit.

  -en, --env, --env_nam : str
    Conda environment to activate. Environment to activate (default: '${env_nam}').

  -ds, --dir_scr : dir
    Directory containing scripts and functions.

  -t, --thr, --threads : int
    Number of threads to use per job (default: ${threads}).

  -md, --mode : {'signal', 'ratio', 'coord'}
    Workflow mode: 'signal', 'ratio', or 'coord' (default: '${mode}').

  -me, --method : {'unadj', 'frag', 'norm', 'log2', 'unadj_r', 'log2_r'}
    Workflow method. Computation subtype.

    If '--mode signal', defaults to 'norm'. If '--mode ratio', defaults to 'unadj'. If '--mode coord', this argument is ignored.

  -ci, --csv_fil_in : list of file
    Comma-separated list of input file paths for BAM or CRAM files.

    Used with '--mode signal' or '--mode coord'.

  -rf, --ref_fa : file
    Reference FASTA file for CRAM input files.

    Required when '--csv_fil_in' contains CRAM input.

  -cs, --chr_sizes, --chrom_sizes : file
    Chromosome sizes file in optional UCSC-style TSV format.

    Used with '--mode signal' or '--mode coord' to supplement BAM/CRAM header sizes, and with '--mode ratio' to validate bedGraph interval bounds.

  -cA, --csv_fil_A : list of file
    Comma-separated list of file A paths for numerator bedGraph files.

    Used with '--mode ratio'.

  -cB, --csv_fil_B : list of file
    Comma-separated list of file B paths for denominator bedGraph files.

    Used with '--mode ratio'.

  -co, --csv_fil_out : list of file
    Comma-separated list of output file paths.

  -tr, --track : flag
    Write a companion track file. If '--mode ratio', write a companion bedGraph without non-finite rows.

  -sb, --siz_bin : int
    Bin size in base pairs for signal computation (default: ${siz_bin}).

    Used with '--mode signal'.

  -ck, --chunk_size : int
    Number of records to process per chunk. Reserved compatibility option for compute_signal.py (default: ${chunk_size}).

    Used with '--mode signal'. Ignored by the public 'chrom' and 'window' engines.

  -eg, --engine : {'chrom', 'window'}
    Processing engine for signal computation (default: '${engine}').

    Used with '--mode signal'.

    Engine selection:
      - 'chrom': default; best general choice and current best CRAM choice.
      - 'window': recommended to try for large BAM inputs.

  -csf, --csv_scl_fct : list of structured string
    Comma-separated list of scaling factors or sentinels.

    Used with '--mode signal' or '--mode ratio'.

  -cuf, --csv_usr_frg : list of int
    Comma-separated list of fixed fragment-length values. Optional comma-separated list of fragment lengths or sentinels.

    Used with '--mode signal' or '--mode coord'.

  -cdm, --csv_dep_min : list of number
    Comma-separated list of minimum-depth values. Optional comma-separated list of minimum input depth values or sentinels.

    Used with '--mode ratio'.

  -cps, --csv_pseudo : list of structured string
    Comma-separated list of pseudocount values. Optional comma-separated list of per-sample pseudocount specs 'A[:B]'.

    Used with '--mode ratio'.

  -e, --eps : float
    Zero tolerance epsilon or sentinel used for ratio-mode zero checks.

  -s0, --skip_00 : {'pre_scale', 'post_scale'}
    Skip rows where both compared values are zero. Shared zero-zero skip mode or sentinel for ratio computation.

    Accepted zero-zero skip modes are 'pre_scale' and 'post_scale'.

  --strict_bins : flag
    Require strict bin compatibility. If '--mode ratio', require both input bedGraphs to have the same ordered '(chrom, start, end)' grid across all data rows.

  -dn, --drp_nan : flag
    Drop non-finite values from main output. Applies to ratio mode.

  -sp, --skp_pfx : list of str
    Comma-separated list of header prefixes to skip. Shared comma-separated bedGraph header prefixes or sentinel to skip.

  -dp, --dp : int
    Maximum number of decimal places retained for finite emitted values (default: ${dp}).

  -deo, --dir_eo : dir
    Directory for stderr and stdout log files.

  -nj, --nam_job : str
    Job name. Prefix for job names (default: 'compute_\${mode}_\${method}' for '--mode signal' and '--mode ratio'; 'compute_\${mode}' for '--mode coord').

Notes
-----
  Runtime requirements:
    - A compatible Conda environment providing the listed tools
    - bash >= 4.4
    - conda (when the requested environment is not active)
    - python >= 3.11
    - Reference FASTA and required index (when processing CRAM)

  - BAM/CRAM and bedGraph input files must be coordinate-sorted.
  - CRAM inputs in '--mode signal' or '--mode coord' require '--ref_fa'.
  - Input and output paths supplied to this wrapper interface must not contain spaces, commas, or semicolons.
  - Dash input/output ('-') is not supported by this wrapper or the underlying Python scripts.
  - Use consistent file ordering in input and output lists.
  - For ratio mode, keep numerator and denominator files in corresponding order.
  - To run in debug mode, set hardcoded variable 'debug=true'.
  - To run in dry-run mode, set hardcoded variable 'dry_run=true'.
  - To run in parse-only mode, set hardcoded variable 'p_only=true'.
  - To run in parse-and-check-only mode, set hardcoded variable 'pc_only=true'.

Examples
--------
  1. Compute normalized signal from two coordinate-sorted BAM files.
    '''bash
    bash "\${dir_scr}/submit_compute_signal.sh" \\
        --env_nam "env_protocol" \\
        --dir_scr "\${dir_scr}" \\
        --threads 4 \\
        --mode "signal" \\
        --method "norm" \\
        --csv_fil_in "\${dir_bam}/sample_1.bam,\${dir_bam}/sample_2.bam" \\
        --csv_fil_out "\${dir_out}/sample_1.bdg.gz,\${dir_out}/sample_2.bdg.gz" \\
        --siz_bin 10 \\
        --engine "chrom" \\
        --csv_scl_fct "NA,NA" \\
        --csv_usr_frg "NA,NA" \\
        --dir_eo "\${dir_eo}" \\
        --nam_job "compute_signal_norm"
    '''

  2. Compute log2 ratios from paired numerator and denominator bedGraphs.
    '''bash
    bash "\${dir_scr}/submit_compute_signal.sh" \\
        --env_nam "env_protocol" \\
        --dir_scr "\${dir_scr}" \\
        --threads 1 \\
        --mode "ratio" \\
        --method "log2" \\
        --csv_fil_A "\${dir_bdg}/IP_1.bdg.gz,\${dir_bdg}/IP_2.bdg.gz" \\
        --csv_fil_B "\${dir_bdg}/in_1.bdg.gz,\${dir_bdg}/in_2.bdg.gz" \\
        --csv_fil_out "\${dir_out}/ratio_1.bdg.gz,\${dir_out}/ratio_2.bdg.gz" \\
        --csv_scl_fct "NA,NA" \\
        --csv_dep_min "NA,NA" \\
        --csv_pseudo "NA,NA" \\
        --eps 0 \\
        --skip_00 "pre_scale" \\
        --track \\
        --dp 6 \\
        --dir_eo "\${dir_eo}" \\
        --nam_job "compute_ratio_log2"
    '''
EOM
}
