#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: help_submit_compute_signal.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6);
# - Anthropic Claude Code (Opus 5, Fable 5).
#
# Distributed under the MIT license.


# TODO: do we need '--dir_scr' in the examples? Only if using 'sbatch'.
function help_submit_compute_signal() {
    # The submit owner initializes interpolated defaults before invocation.
    # shellcheck disable=SC2154
    cat >&2 << EOM
Usage
-----
  submit_compute_signal.sh
    [--help]
    [--env_nam <str>] [--dir_scr <dir>] [--threads <int>]
    [--mode <mode>] [--method <method>]
    (--csv_fil_in <csv> [--ref_fa <file>] | --csv_fil_A <csv> --csv_fil_B <csv> [--chr_siz <file>])
    --csv_fil_out <csv>
    [--siz_bin <int>] [--engine <engine>] [--siz_win <int>] [--csv_scl_fct <csv>] [--csv_usr_frg <csv>]
    [--csv_dep_min <csv>] [--csv_pseudo <csv>] [--eps <flt>] [--skip_00 <choice>] [--strict_bins] [--drp_nan] [--skp_pfx <csv>]
    [--csv_report_n_frg <csv>] [--csv_report_n_ovlp <csv>]
    [--track] [--dp <int>]
    --dir_eo <dir> [--nam_job <str>]

  Submit per-sample signal, ratio, or fragment-coordinate jobs from comma-separated file lists to 'compute_signal.py' or 'compute_signal_ratio.py'.

Parameters
----------
  -h, --help : flag
    Print this help message and exit.

  -en, --env, --env_nam : str
    Conda environment to activate (default: '${env_nam}').

  -ds, --dir_scr : dir
    Directory containing scripts and functions. Passed by the 'execute_*.sh' wrappers, and needed when this script is run from a copy, as 'sbatch <script>' does, rather than from its real path.

  -t, --thr, --threads : int
    Number of threads to use per job (default: ${threads}).

  -md, --mode : {'signal', 'ratio', 'coord'}
    Workflow mode: 'signal', 'ratio', or 'coord' (default: '${mode}').

  -me, --method : {'unadj', 'frag', 'norm', 'count', 'cpm', 'linear', 'log2', 'linear_r', 'log2_r'}
    Workflow method. Signal or ratio computation subtype.

    With '--mode signal':
      - 'unadj' gives base-pair overlap with no adjustment;
      - 'frag' gives base-pair overlap divided by fragment length;
      - 'norm' (alias 'nc') gives base-pair overlap divided by fragment length and total fragments, so genome-wide coverage sums to 1;
      - 'count' deposits one whole count per bin a fragment touches, so the track sums to the overlap count 'L'; and
      - 'cpm' rescales those whole counts to sum to one million.

    With '--mode ratio':
      - 'linear' gives 'fil_A / fil_B' (e.g., IP/input) and 'log2' (alias 'l2') gives its log2, which is symmetric about zero.
      - The 'linear_r' and 'log2_r' (alias 'l2_r') variants invert the comparison to 'fil_B / fil_A' ("r" stands for "reciprocal").

    If '--mode signal', defaults to 'norm'. If '--mode ratio', defaults to 'linear'. If '--mode coord', this argument is ignored.

  -ci, --csv_fil_in : list of file
    Comma-separated list of input file paths for BAM or CRAM files.

    Used with '--mode signal' or '--mode coord'.

  -rf, --ref_fa : file
    Reference FASTA file for CRAM input files.

    Required when '--csv_fil_in' contains CRAM input.

  -cs, --chr_siz : file
    Chromosome sizes file in UCSC-style TSV format.

    Used only with '--mode ratio' to validate bedGraph interval bounds. Signal and coordinate modes take chromosome sizes from the BAM/CRAM header.

  -cA, --csv_fil_A : list of file
    Comma-separated list of file A paths for numerator bedGraph files.

    Used with '--mode ratio'.

  -cB, --csv_fil_B : list of file
    Comma-separated list of file B paths for denominator bedGraph files.

    Used with '--mode ratio'.

  -co, --csv_fil_out : list of file
    Comma-separated list of output file paths.

  -sb, --siz_bin : int
    Bin size in base pairs for signal computation (default: ${siz_bin}).

    Used with '--mode signal'.

  -eg, --engine : {'chrom', 'window'}
    Processing engine for signal computation (default: '${engine}').

    Used with '--mode signal'.

    Both engines dispatch indexed fetch tasks and produce the same signal; they differ only in how fetch work is divided among threads.
      - 'chrom': one fetch task per chromosome. Task size tracks chromosome size, so the longest chromosomes dominate wall time.
      - 'window': each chromosome is split into fixed-size coordinate windows, with one fetch task per window. Task sizes are uniform, giving finer load balance across threads, but at the cost of more fetch calls. Window size is set by '--siz_win'.

    Recommended: keep 'chrom' as the general choice and the current best choice for CRAM input; try 'window' for large BAM inputs.

  -sw, --siz_win : int
    Window size in base pairs for the 'window' engine's indexed fetch tasks (default: ${siz_win}). Ignored by the 'chrom' engine. Used only with '--mode signal'.

  -csf, --csv_scl_fct : list of structured string
    Comma-separated list of scaling factors or sentinels.

    Used with '--mode signal' or '--mode ratio'.

    A factor multiplies the finished track verbatim, whatever '--method' produced it.

    A siQ-ChIP alpha is derived against a ratio of signal types: equations 5 and 6 for 'frag', 5nd and 6nd for 'norm'.

    A spike-in alpha is likewise derived against a ratio of signal types: either 'count' or 'unadj', since the per-bin deposition difference largely cancels in the division. Earlier spike-in work used 'count', so prefer it where continuity with those results matters.

  -cuf, --csv_usr_frg : list of int
    Comma-separated list of fixed fragment-length values or sentinels.

    Used with '--mode signal' or '--mode coord'.

  -cdm, --csv_dep_min : list of number
    Comma-separated list of minimum-depth values or sentinels; here 'min' abbreviates minimum.

    Used with '--mode ratio'.

  -cps, --csv_pseudo : list of structured string
    Comma-separated list of per-sample pseudocount specs 'A[:B]'.

    Used with '--mode ratio'.

  -e, --eps : float
    Zero tolerance epsilon or sentinel used for ratio-mode zero checks.

  -s0, --skip_00 : {'pre_scale', 'post_scale'}
    Skip rows where both compared values are zero. Shared zero-zero skip mode or sentinel for ratio computation.

    Accepted zero-zero skip modes are 'pre_scale' and 'post_scale'.

  -stn, --strict_bins : flag
    Require strict bin compatibility. If '--mode ratio', require both input bedGraphs to have the same ordered '(chrom, start, end)' grid across all data rows.

  -dn, --drp_nan : flag
    Drop non-finite values from main output. Applies to ratio mode.

  -sp, --skp_pfx : list of str
    Comma-separated list of header prefixes to skip. Shared comma-separated bedGraph header prefixes or sentinel to skip.

  -crf, --csv_report_n_frg : list of file
    Comma-separated list of paths for per-sample fragment-count reports. Supply one path per '--csv_fil_in' element. Used only with '--mode signal'.

  -cro, --csv_report_n_ovlp : list of file
    Comma-separated list of paths for per-sample overlap-count reports. 'L' counts fragment-bin overlaps, so it depends on '--siz_bin'. Used only with '--mode signal'.

  -tr, --track : flag
    Write a companion track file. If '--mode ratio', write a companion bedGraph without non-finite rows.

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
        --csv_fil_out "\${dir_out}/sample_1.bedGraph.gz,\${dir_out}/sample_2.bedGraph.gz" \\
        --siz_bin 10 \\
        --engine "chrom" \\
        --csv_scl_fct "0.8734,1.1290" \\
        --csv_usr_frg "150,150" \\
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
        --csv_fil_A "\${dir_bdg}/IP_1.bedGraph.gz,\${dir_bdg}/IP_2.bedGraph.gz" \\
        --csv_fil_B "\${dir_bdg}/in_1.bedGraph.gz,\${dir_bdg}/in_2.bedGraph.gz" \\
        --csv_fil_out "\${dir_out}/ratio_1.bedGraph.gz,\${dir_out}/ratio_2.bedGraph.gz" \\
        --csv_scl_fct "0.8734,1.1290" \\
        --csv_pseudo "0.5,0.4" \\
        --eps 0 \\
        --skip_00 "pre_scale" \\
        --track \\
        --dp 6 \\
        --dir_eo "\${dir_eo}" \\
        --nam_job "compute_ratio_log2"
    '''

  3. Write per-sample counts beside each signal track.
    '''bash
    bash "\${dir_scr}/submit_compute_signal.sh" \\
        --env_nam "env_protocol" \\
        --dir_scr "\${dir_scr}" \\
        --threads 4 \\
        --mode "signal" \\
        --method "norm" \\
        --csv_fil_in "\${dir_bam}/sample_1.bam,\${dir_bam}/sample_2.bam" \\
        --csv_fil_out "\${dir_out}/sample_1.bedGraph.gz,\${dir_out}/sample_2.bedGraph.gz" \\
        --csv_report_n_frg "\${dir_out}/sample_1.n_frg.txt,\${dir_out}/sample_2.n_frg.txt" \\
        --csv_report_n_ovlp "\${dir_out}/sample_1.n_ovlp.txt,\${dir_out}/sample_2.n_ovlp.txt" \\
        --siz_bin 10 \\
        --dir_eo "\${dir_eo}" \\
        --nam_job "compute_signal_norm"
    '''

    Supply one path per '--csv_fil_in' element. Counting happens before the output branch, so the values match those of a run that writes no track.

  4. Count fragments and fragment-bin overlaps without writing a track.
    '''bash
    bash "\${dir_scr}/submit_compute_signal.sh" \\
        --env_nam "env_protocol" \\
        --dir_scr "\${dir_scr}" \\
        --threads 4 \\
        --mode "signal" \\
        --csv_fil_in "\${dir_bam}/sample_1.bam,\${dir_bam}/sample_2.bam" \\
        --csv_report_n_frg "\${dir_out}/sample_1.n_frg.txt,\${dir_out}/sample_2.n_frg.txt" \\
        --csv_report_n_ovlp "\${dir_out}/sample_1.n_ovlp.txt,\${dir_out}/sample_2.n_ovlp.txt" \\
        --siz_bin 10 \\
        --dir_eo "\${dir_eo}" \\
        --nam_job "compute_signal_counts"
    '''

    Omitting '--csv_fil_out' is permitted only because a report list is given; 'L' counts fragment-bin overlaps, so it still depends on '--siz_bin'.

  5. Compute counts-per-million signal from whole-count deposition.
    '''bash
    bash "\${dir_scr}/submit_compute_signal.sh" \\
        --env_nam "env_protocol" \\
        --dir_scr "\${dir_scr}" \\
        --threads 4 \\
        --mode "signal" \\
        --method "cpm" \\
        --csv_fil_in "\${dir_bam}/sample_1.bam,\${dir_bam}/sample_2.bam" \\
        --csv_fil_out "\${dir_out}/sample_1.bedGraph.gz,\${dir_out}/sample_2.bedGraph.gz" \\
        --csv_report_n_ovlp "\${dir_out}/sample_1.n_ovlp.txt,\${dir_out}/sample_2.n_ovlp.txt" \\
        --siz_bin 10 \\
        --dir_eo "\${dir_eo}" \\
        --nam_job "compute_signal_cpm"
    '''

    Each fragment deposits one whole count per touched bin, and the track is rescaled to sum to one million. '--csv_report_n_ovlp' writes the divisor 'L' beside each track.
EOM
}
