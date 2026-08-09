#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: help_execute_compute_signal.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
#   GPT-5.6);
# - Anthropic Claude Code (Opus 5).
#
# Distributed under the MIT license.


# TODO: Add more examples and PMIDs for all studies mentioned in detailed help.
usage=$(cat << EOM
Usage
-----
  execute_compute_signal.sh
    [--help] [--details] [--all_help]
    [--verbose] [--dry_run]
    [--env_nam <str>] [--threads <int>]
    [--mode <mode>] [--method <method>]
    (--csv_fil_in <csv> [--ref_fa <file>] [--chr_sizes <file>] | --csv_fil_A <csv> --csv_fil_B <csv> [--chr_sizes <file>])
    --dir_out <dir> [--typ_out <format>] [--prefix <str>]
    [--siz_bin <int>] [--engine <engine>] [--chunk_size <int>] [--csv_usr_frg <csv>] [--csv_scl_fct <csv>]
    [--csv_dep_min <csv>] [--csv_pseudo <csv>] [--eps <num>] [--skip_00 <choice>] [--strict_bins] [--drp_nan] [--skp_pfx <csv>]
    [--track] [--dp <int>]
    [--dir_eo <dir>] [--nam_job <str>] [--max_job <int>] [--slurm] [--time <time>]

EOM
)

# shellcheck disable=SC2154
function help_execute_compute_signal() {
cat << EOM
${usage}


  Coordinate and automate the computation of signal tracks, ratio tracks, or fragment-coordinate files from BAM/CRAM or bedGraph input files. Supports multiple normalization strategies and runs computations in serial or parallel via GNU Parallel or Slurm.

  For more details on what this script can do, including notes and usage examples, run either of the following:
    '''bash
    bash path/to/execute_compute_signal.sh --details
    bash path/to/execute_compute_signal.sh --all_help
    '''


Parameters
----------
  -h, --help : flag
    Print this short help message and exit.

  -d, --details : flag
    Print full documentation with notes and examples, then exit.

  -ah, --all_help : flag
    Print both the short and full documentation, then exit.

  -v, --verbose : flag
    Run script in verbose mode.

  -dr, --dry, --dry_run : flag
    Run script in dry-run mode.

  -en, --env, --env_nam : str
    Conda environment to activate (default: '${env_nam}').

  -t, --thr, --threads : int
    Number of threads to use (default: ${threads}).

  -md, --mode : {'signal', 'ratio', 'coord'}
    Workflow mode: 'signal', 'ratio', or 'coord' (default: '${mode}').

    See '--details' for full synonym lists.

  -me, --method : {'unadj', 'frag', 'norm', 'log2', 'unadj_r', 'log2_r'}
    Workflow method. Signal or ratio computation subtype (used only with '--mode signal' or '--mode ratio'; default if '--mode signal': norm; default if '--mode ratio': unadj).

    For '--mode signal', the main choices are
      - 'unadj' (per-bin totals),
      - 'frag' (fragment-length adjusted), and
      - 'norm' (fragment- and library-size normalized so that coverage integrates to 1).

    For '--mode ratio', the main choices are
      - 'unadj' ('file A / file B'),
      - 'log2' ['log2(file A / file B)'],
      - 'unadj_r' ('file B / file A'), and
      - 'log2_r' ['log2(file B / file A)' = '-log2(file A / file B)'].

    See '--details' for full synonym lists and references.

  -ci, --csv_fil_in : list of file
    Comma-separated list of input file paths for coordinate-sorted BAM/CRAM files (used only with '--mode signal' or '--mode coord').

  -rf, --ref_fa : file
    Reference FASTA file for CRAM input files (required if any '--csv_fil_in' element ends in '.cram'; used only with '--mode signal' or '--mode coord').

  -cs, --chr_sizes, --chrom_sizes : file
    Chromosome sizes file in optional UCSC-style TSV format.

    Used with '--mode signal' or '--mode coord' to supplement BAM/CRAM header sizes, and with '--mode ratio' to validate bedGraph interval bounds.

  -cA, --csv_fil_A : list of file
    Comma-separated list of file A paths for coordinate-sorted numerator bedGraph files (e.g., IP signal; used only with '--mode ratio').

  -cB, --csv_fil_B : list of file
    Comma-separated list of file B paths for coordinate-sorted denominator bedGraph files (e.g., input signal; used only with '--mode ratio').

  -do, --dir_out : dir
    Output directory for generated files.

  -to, --typ_out : {'bedGraph', 'bedGraph.gz', 'bedgraph', 'bedgraph.gz', 'bdg', 'bdg.gz', 'bg', 'bg.gz', 'bed', 'bed.gz'}
    Output file format (default: '${typ_out}').

    For '--mode signal' or '--mode ratio', the typical choice is a bedGraph-style track (e.g., 'bedGraph.gz' or 'bdg.gz').

    For '--mode coord', the typical choice is a BED-like coordinate file (e.g., 'bed.gz').

    If an incompatible combination is requested, the value is coerced to a sensible default with a warning; see '--details' for the full list and mode-specific behavior.

  -px, --pfx, --prefix : str
    Custom prefix to prepend to output filenames.

  -sb, --siz_bin : int
    Bin size in base pairs for signal computation (used only with '--mode signal'; default: 10).

  -eg, --engine : {'chrom', 'window'}
    Processing engine for signal computation (used only with '--mode signal'; default: 'chrom').

    Engine selection:
      - 'chrom': default; best general choice and current best CRAM choice.
      - 'window': recommended to try for large BAM inputs.

  -ck, --chunk_size : int
    Number of records to process per chunk. Reserved compatibility option for compute_signal.py (used only with '--mode signal'; default: 100000).

    Ignored by the public 'chrom' and 'window' engines.

  -csf, --csv_scl_fct : list of structured string
    Comma-separated list of scaling factors or sentinels. Used only with '--mode signal' or '--mode ratio'.

    For '--mode signal', each element must be 'NA' or a positive scalar float.

    For '--mode ratio', each element may be 'NA', a positive scalar float, or a positive 'A:B' spec, where A scales '--csv_fil_A' and B scales '--csv_fil_B' before ratio calculation.

  -cuf, --csv_usr_frg : list of int
    Comma-separated list of fixed fragment-length values or sentinels. Used with '--mode signal' or '--mode coord'.

  -cdm, --csv_dep_min : list of number
    Comma-separated list of minimum-depth values or sentinels; here 'min' abbreviates minimum. Used only with '--mode ratio'.

  -cps, --csv_pseudo : list of structured string
    Comma-separated list of pseudocount values as per-sample specs 'A[:B]' or sentinels. Used only with '--mode ratio'.

  -e, --eps : number
    Zero tolerance epsilon or sentinel used for ratio-mode zero checks. Used only with '--mode ratio'.

  -s0, --skp_00, --skip_00 : {'pre_scale', 'post_scale'}
    Skip rows where both compared values are zero. Shared zero-zero skip mode or sentinel for ratio computation: 'pre_scale' or 'post_scale'. Used only with '--mode ratio'.

  --strict_bins : flag
    Require strict bin compatibility. If '--mode ratio', require both input bedGraph files to have the same ordered '(chrom, start, end)' grid across all data rows.

  -dn, --drp_nan, --drop_nan : flag
    Drop non-finite ratio rows ('inf', '-inf', and 'nan') from the main ratio output. Used only with '--mode ratio'.

  -sp, --skp_pfx : list of str
    Comma-separated list of header prefixes to skip. Shared comma-separated list of bedGraph header prefixes or sentinel to skip while parsing ratio-mode input files. Used only with '--mode ratio'.

  -tr, --trk, --track : flag
    Write a companion track file. If '--mode ratio', also write a companion bedGraph with all non-finite rows ('inf', '-inf', and 'nan') removed.

  -dp, --dp : int
    Maximum number of decimal places retained for finite emitted values (default: ${dp}).

    After rounding, non-informative trailing zeros are stripped.

  -deo, --dir_eo : dir
    Directory for stderr and stdout log files in TXT format (default: '\${dir_out}/logs').

  -nj, --nam_job : str
    Job name. Prefix for job names (default depends on resolved '--mode' and '--method'; e.g., 'compute_signal_norm', 'compute_ratio_unadj', or 'compute_coord').

  -mj, --max_job : int
    Maximum number of jobs to run concurrently (default: ${max_job}).

  -sl, --slurm : flag
    Submit jobs to the Slurm scheduler.

  -tm, --time : time
    Slurm job time limit in 'h:mm:ss' format (required if '--slurm'; default: '${time}').


Notes
-----
  Runtime requirements:
    - A compatible Conda environment providing the listed tools
    - basename
    - bash >= 4.4
    - dirname
    - parallel (when '--slurm' is not specified and multiple jobs are run)
    - python >= 3.11
    - rm (when '--slurm' is not specified and multiple jobs are run)
    - sbatch (when '--slurm' is specified)
    - tr

  - BAM/CRAM and bedGraph input files must be coordinate-sorted.
  - CRAM inputs require '--ref_fa'.
  - Input and output paths supplied to this wrapper interface must not contain spaces, commas, or semicolons.
  - See 'execute_compute_signal.sh --details' for more notes.

Examples
--------
  1. Print the concise compute-signal wrapper help to stdout.
    '''bash
    help_execute_compute_signal
    '''

  2. Save the concise help for offline review.
    '''bash
    help_execute_compute_signal > compute_signal.help.txt
    '''
EOM
}


function detail_execute_compute_signal() {
    local mode="${1:-}"

    #  Only print the top-level 'Usage' block if '--no-usage' is not invoked
    if [[ ! "${mode}" =~ ^--no[_-]usage$ ]]; then
cat << EOM
${usage}


EOM
    fi

#TODO: detailed descriptions of new arguments:
#+     - --chr_sizes <file>
#+     - --engine <engine>
#+     - --chunk_size <int>
cat << EOM
  Driver script automating the computations of bedGraph signal or ratio tracks, or BED-like fragment coordinate files, from BAM/CRAM (for signal tracks or fragment coordinate files) or bedGraph (for ratio tracks) input files.

  Supports multiple signal normalization strategies, including
    - unadjusted (raw) signal (i.e., per-bin totals with no fragment-length and/or library-size adjustments)
    - fragment-length adjusted signal (Dickson et al., JBC 2020; Dickson et al., Sci Rep 2023)
    - normalized coverage (Dickson et al., Sci Rep 2023)
    - siQ-ChIP IP efficiency (input-normalized ratio-based analyses; Dickson et al., JBC 2020; Dickson et al., Sci Rep 2023)
    - log2-transformed input-normalized signal ratios (e.g., as described in "Data Analysis G" of Alavattam et al., Bio-protocol 2025)
    - spike-in-normalized signal (i.e., unmodified ChIP-Rx; Orlando et al., Cell Rep 2014)
    - input- and spike-in-normalized signal ratios (modified ChIP-Rx as described in "Data Analysis I" of Alavattam et al., Bio-protocol 2025)

  Supports serial job execution and parallel job execution via Slurm or GNU Parallel.


Parameters
----------
  -h, --help : flag
    Print a short help message and exit.

  -d, --details : flag
    Print this full documentation with notes and examples, then exit.

  -ah, --all_help : flag
    Print both the short and full documentation, then exit.

  -v, --verbose : flag
    Run script in verbose mode.

  -dr, --dry, --dry_run : flag
    Run script in dry-run mode.

  -en, --env, --env_nam : str
    Conda environment to activate (default: '${env_nam}').

  -t, --thr, --threads : int
    Number of threads to use (default: ${threads}).

  -md, --mode : {'signal', 'ratio', 'coord'}
    Workflow mode: 'signal', 'ratio', or 'coord' (default: '${mode}'). Available options:
      - 's', 'sig', 'signal':
        + Compute signal tracks directly from BAM/CRAM input files.
        + Supports unadjusted signal, fragment-length-adjusted signal, or normalized coverage.
        + See '--method' for calculation styles.
        + Use this option with appropriately computed scaling factors ('--csv_scl_fct') to compute unmodified ChIP-Rx (Orlando et al., Cell Rep 2014) spike-in-normalized signal.
        + If 's' or 'sig' are supplied, variable 'mode' is set to "signal".

      - 'r', 'rat', 'ratio':
        + Compute IP/input signal ratios from bedGraph files.
        + Supports linear ratios, log2 ratios, and reciprocal variants of both.
        + Untransformed ratios are used in the modified ChIP-Rx spike-in and siQ-ChIP normalizations described in Alavattam et al., Bio-protocol 2025.
        + Use this option with '--method log2' to compute log2(IP/input) ratios.
        + If 'r' or 'rat' are supplied, variable 'mode' is set to "ratio".

      - 'c', 'coord', 'coordinates':
        + Instead of computing signal or ratio tracks, output fragment coordinates in BED-like format from BAM/CRAM input files.
        + Use this option to prepare input files for the original siQ-ChIP implementation (Dickson et al., JBC 2020; Dickson et al., Sci Rep 2023).
        + For more details, see github.com/BradleyDickson/siQ-ChIP or github.com/kalavattam/siQ-ChIP.
        + This mode disables '--siz_bin' and '--method', and sets '--typ_out' to 'bed.gz' by default (or 'bed' if '--typ_out bed' is specified).
        + If 'c' or 'coordinates' are supplied, variable 'mode' is set to "coord".

  -me, --method : {'unadj', 'frag', 'norm', 'log2', 'unadj_r', 'log2_r'}
    Workflow method. Signal or ratio computation subtype, used only with '--mode signal' or '--mode ratio' (default if '--mode signal': norm; default if '--mode ratio': unadj).
      - If '--mode signal', then the available options are
        + 'u', 'unadj', 'unadjusted', 's', 'smp', 'simple', 'r', 'raw':
          - Compute unadjusted signal (per-bin totals with no fragment-length and/or library-size adjustments).
          - Internally, all of these values are standardized to 'method=unadj'.

        + 'f', 'frg', 'frag', 'frg_len', 'frag_len', 'l', 'len', 'len_frg', 'len_frag':
          - Adjust signal by fragment length.
          - Internally, all of these values are standardized to 'method=frag'.
          - For example, use this option to compute siQ-ChIP-scaled signal with the initial equation described in Dickson et al., JBC 2020 or Equation 5 in Dickson et al., Sci Rep 2023.

        + 'n', 'nrm', 'norm', 'normalized':
          - Compute normalized coverage per Dickson et al., Sci Rep 2023, adjusting by both fragment length and the total number of fragments so that the genome-wide coverage sums to 1.
            + That is, the coverage integrates to unity and can be interpreted as a probability distribution over the genome.
          - Internally, all of these values are standardized to 'method=norm'.

      - If '--mode ratio', then the available options are
        + 'u', 'unadj', 'unadjusted', 's', 'smp', 'simple', 'r', 'raw':
          - Compute simple, unadjusted (non-log2) fil_A/fil_B ratio (e.g., IP/input): 'ratio = fil_A / fil_B'.
          - Internally, all of these values are standardized to 'method=unadj'.

        + '2', 'l2', 'lg2', 'log2':
          - Compute log2(fil_A/fil_B) ratio [e.g., log2(IP/input)]: 'ratio = log2(fil_A / fil_B)'.
          - Internally, all of these values are standardized to 'method=log2'.

        + 'ur', 'unadj_r', 'unadjusted_r', 'sr', 'smp_r', 'simple_r', 'rr', 'raw_r':
          - Compute the reciprocal of the simple, unadjusted (non-log2) ratio: 'ratio = fil_B / fil_A = 1 / (fil_A / fil_B)'.
          - Internally, all of these values are standardized to 'method=unadj_r'.

        + '2r', 'l2r', 'l2_r', 'lg2_r', 'log2_r':
          - Compute the reciprocal of the log2(fil_A/fil_B) ratio: 'ratio = log2(fil_B / fil_A) = -log2(fil_A / fil_B)'.
          - Internally, all of these values are standardized to 'method=log2_r'.

  -ci, --csv_fil_in : list of file
    Comma-separated list of input file paths for coordinate-sorted BAM/CRAM files.

    Required when '--mode signal' or '--mode coord'. Ignored for '--mode ratio'.

  -rf, --ref_fa : file
    Reference FASTA file for CRAM input files.

    Required when any '--csv_fil_in' element ends in '.cram'. Used only with '--mode signal' or '--mode coord'; ignored for '--mode ratio'.

  -cA, --csv_fil_A : list of file
    Comma-separated list of file A paths for coordinate-sorted numerator bedGraph files (e.g., ChIP IP signal tracks).

    Use with '--mode ratio'.

    The list order must match that of '--csv_fil_B' files.

  -cB, --csv_fil_B : list of file
    Comma-separated list of file B paths for coordinate-sorted denominator bedGraph files (e.g., input signal tracks).

    Use with '--mode ratio'.

    The list order must match that of '--csv_fil_A' files.

  -do, --dir_out : dir
    Output directory for generated files:
      - Signal tracks if '--mode signal'.
      - Ratio tracks if '--mode ratio'.
      - BED-like files of fragment coordinates if '--mode coord'.

  -to, --typ_out : {'bedGraph', 'bedGraph.gz', 'bedgraph', 'bedgraph.gz', 'bdg', 'bdg.gz', 'bg', 'bg.gz', 'bed', 'bed.gz'}
    Output file format for signal track output files (default: '${typ_out}'). Available options:
      - 'bedGraph', 'bedgraph', 'bdg', 'bg':
        + Signal/ratio in bedGraph format.
        + Intended for '--mode signal' or '--mode ratio'; with '--mode coord' these values are accepted but coerced to 'bed.gz' (see Notes).

      - 'bedGraph.gz', 'bedgraph.gz', 'bdg.gz', 'bg.gz':
        + Signal/ratio in gzip-compressed bedGraph format.
        + Intended for '--mode signal' or '--mode ratio'; with '--mode coord' these values are accepted but coerced to 'bed.gz' (see Notes).

      - 'bed', 'bed.gz':
        + BED-like format for fragment coordinates instead of signal.
        + Intended for '--mode coord'; with '--mode signal' or '--mode ratio' these values are accepted but coerced to 'bdg.gz' (see Notes).

  -px, --pfx, --prefix : str
    Custom prefix to prepend to output filenames.

    When '--mode signal' or '--mode coord':
      - If not specified, no prefix is added.
      - If specified, the prefix is prepended to the base filename, and any leading 'IP_' or 'in_' in that base filename is stripped before applying it.
      - Note: If you use the same '--prefix' for both IP and input BAMs/CRAMs in the same output directory, stripping 'IP_'/ 'in_' may cause filename collisions.
        + In that case, prefer distinct prefixes (e.g., 'IP', 'in').

    When '--mode ratio':
      - If not specified, a default prefix is automatically constructed based on '--method' and '--csv_scl_fct'; for example:
        + 'rat' (default)
        + 'log2_rat' (if '--method log2')
        + 'recip_rat' (if '--method unadj_r')
        + 'log2_recip_rat' (if '--method log2_r')
        + 'scl_rat' (if '--method unadj' and '--csv_scl_fct' is supplied)
        + 'scl_log2_rat' (if '--method log2' and '--csv_scl_fct' is supplied)
        + 'scl_recip_rat' (if '--method unadj_r' and '--csv_scl_fct' is supplied)
        + 'scl_log2_recip_rat' (if '--method log2_r' and '--csv_scl_fct' is supplied)
      - If specified, the custom prefix replaces the default.
      - Whether specified or not, any leading 'IP_' string in the base name is stripped before the prefix.

  -sb, --siz_bin : int
    Bin size in base pairs for signal computation.

    Used only with '--mode signal' (default: 10); ignored otherwise.

  -csf, --csv_scl_fct : list of structured string
    Comma-separated list of scaling factors or sentinels to apply to signal or ratio values.

    Used with either '--mode signal' or '--mode ratio'; ignored otherwise.

    List size must match the number of input files via '--csv_fil_in' or '--csv_fil_A'/'--csv_fil_B'.

    For '--mode signal', each non-sentinel element must be a positive scalar float.

    For '--mode ratio', each non-sentinel element may be either:
      - 'A'    Scale file A by A and file B by 1.0.
      - 'A:B'  Scale file A by A and file B by B.

  -cuf, --csv_usr_frg : list of int
    Comma-separated list of fixed fragment-length values or sentinels to use instead of read lengths (single-end alignments) or template lengths (paired-end alignments).

    Used with either '--mode signal' or '--mode coord'; ignored otherwise.

    List size must match the number of input files via '--csv_fil_in'.

  -cdm, --csv_dep_min : list of number
    Comma-separated list of minimum-depth values or sentinels used to avoid extreme division operations; here 'min' abbreviates minimum.

    Used only with '--mode ratio'; ignored otherwise.

    List size must match the number of input files via '--csv_fil_A'/'--csv_fil_B'.

    Although allowed, using '--csv_dep_min' together with '--csv_pseudo' is usually harder to interpret, since both stabilize low-depth ratio behavior in different ways.

  -cps, --csv_pseudo : list of structured string
    Comma-separated list of pseudocount values or sentinels used during ratio computation.

    Used only with '--mode ratio'; ignored otherwise.

    List size must match the number of input files via '--csv_fil_A'/'--csv_fil_B'.

    Each non-sentinel element may be either:
      - 'A'    Add pseudocount A symmetrically.
      - 'A:B'  Add pseudocount A to file A and B to file B.

    Although allowed, using '--csv_pseudo' together with '--csv_dep_min' is usually harder to interpret, since both stabilize low-depth ratio behavior in different ways.

  -e, --eps : number
    Zero tolerance epsilon or sentinel used for ratio-mode zero checks.

    Used only with '--mode ratio'; ignored otherwise.

    Non-sentinel values must be non-negative floats.

  -s0, --skp_00, --skip_00 : {'pre_scale', 'post_scale'}
    Skip rows where both compared values are zero. Shared zero-zero skip mode or sentinel for ratio computation.

    Used only with '--mode ratio'; ignored otherwise.

    Non-sentinel values must be one of 'pre_scale' or 'post_scale'.

  -dn, --drp_nan, --drop_nan : flag
    Drop non-finite values from the main ratio output.

    Used only with '--mode ratio'; ignored otherwise.

    If set, rows yielding 'inf', '-inf', or 'nan' are omitted from the main ratio output.

  -sp, --skp_pfx : list of str
    Comma-separated list of header prefixes to skip. Shared comma-separated list of bedGraph header prefixes or sentinel to skip.

    Used only with '--mode ratio'; ignored otherwise.

    Passed through to 'submit_compute_signal.sh' and then to 'compute_signal_ratio.py'.

  -tr, --trk, --track : flag
    Write a companion track file. If '--mode ratio', also write a companion bedGraph with all non-finite rows ('inf', '-inf', and 'nan') removed.

    The new file will include '.track' before the extension.

    This cleaned version is ideal for visualization in genome browsers such as IGV, avoiding issues caused by 'inf', '-inf' or 'nan' values.

  -dp, --dp : int
    Maximum number of decimal places retained for finite emitted values (default: ${dp}).

    After rounding, non-informative trailing zeros are stripped.

  -deo, --dir_eo : dir
    Directory for stderr and stdout log files in TXT format (default: '\${dir_out}/logs').

  -nj, --nam_job : str
    Job name. Prefix for job names (default depends on resolved '--mode' and '--method'; e.g., 'compute_signal_norm', 'compute_ratio_log2', or 'compute_coord').

  -mj, --max_job : int
    Maximum number of jobs to run concurrently (default: ${max_job}).

    With '--slurm': maximum number of Slurm array tasks allowed to run concurrently.

    Without '--slurm':
      - If the resolved number of parallel jobs is greater than 1, jobs run in parallel via GNU Parallel.
      - If the resolved number of parallel jobs is 1, jobs run serially.

  -sl, --slurm : flag
    Submit jobs to the Slurm scheduler.

    If '--slurm' is not specified, this script uses '--threads', '--max_job', and the detected CPU core count (via 'set_params_parallel') to decide whether to run jobs in parallel with GNU Parallel or in serial.
      - In this non-Slurm path, '--threads' is treated as a total CPU/thread budget for the local machine, and 'set_params_parallel' converts that into threads per job plus number of parallel jobs.
      - If the resolved number of parallel jobs is greater than 1, jobs are run with GNU Parallel.
      - If the resolved number of parallel jobs is 1, jobs are run serially.

  -tm, --time : time
    Slurm job time limit. The length of time, in 'h:mm:ss' format, for the Slurm job (required if '--slurm' is specified, ignored if not; default: '${time}').


Notes
-----
  Runtime requirements:
    - A compatible Conda environment providing the listed tools
    - basename
    - bash >= 4.4
    - conda (when the requested environment is not active)
    - dirname
    - parallel (when '--slurm' is not specified and multiple jobs are run)
    - python >= 3.11
    - rm (when '--slurm' is not specified and multiple jobs are run)
    - sbatch (when '--slurm' is specified)
    - tr

  - BAM/CRAM and bedGraph input files must be coordinate-sorted.
  - CRAM inputs require '--ref_fa'.
  - Input and output paths supplied to this wrapper interface must not contain spaces, commas, or semicolons.
    + Commas are used internally as list delimiters.
    + Semicolons are also considered unsafe in this wrapper workflow.
    + Spaces in paths are not supported by the current list-serialization and reconstruction logic.
  - If applicable, use consistent file ordering between IP and input files.
  - '--typ_out' must be compatible with the selected '--mode'.
    + With '--mode signal' or '--mode ratio', bedGraph-style values ('bedGraph', 'bedgraph', 'bdg', and 'bg', and their '.gz' variants) are allowed; 'bed'/'bed.gz' are accepted but are automatically converted to 'bdg.gz' with a warning.
    + With '--mode coord', 'bed'/'bed.gz' are allowed; bedGraph-style values are accepted but are automatically converted to 'bed.gz' with a warning.
  - Output file path. Output filenames are derived from BAM/CRAM or bedGraph input files and the value associated with '--typ_out'.
  - For bedGraph-style output, '--dp' sets the maximum number of decimal places retained for finite emitted values; after rounding, non-informative trailing zeros and any trailing decimal point are stripped.
  - BED-like files of fragment coordinates are, e.g., used as input to the original siQ-ChIP implementation (Dickson et al., JBC 2020; Dickson et al., Sci Rep 2023).
  - Job execution mode (serial, GNU Parallel, or Slurm array) is chosen automatically from '--slurm', '--threads', and '--max_job':
    + If '--slurm' is specified:
      - Jobs are submitted as a Slurm array.
      - '--max_job' (after adjustment) sets the maximum number of array tasks running concurrently.
      - '--threads' controls '--cpus-per-task' for each array element.
    + If '--slurm' is not specified:
      - Helper function 'set_params_parallel' uses '--threads', '--max_job', and the detected CPU core count to determine a safe combination of threads per job and the number of parallel jobs.
      - If the resulting number of parallel jobs is greater than 1, commands are written to a configuration file and executed with GNU Parallel.
      - If the resulting number of parallel jobs is 1, all jobs are run serially in a single Bash process (neither GNU Parallel nor Slurm).


Examples
--------
  1. Compute normalized coverage using GNU Parallel.
    '''bash
    bash "\${HOME}/bin/execute_compute_signal.sh" \\
        --threads 8 \\
        --mode "signal" \\
        --method "norm" \\
        --csv_fil_in "\${HOME}/project/samples/sample_1.bam,\${HOME}/project/samples/sample_2.bam" \\
        --dir_out "\${HOME}/project/tracks" \\
        --typ_out "bdg.gz" \\
        --siz_bin 50 \\
        --dir_eo "\${HOME}/project/logs" \\
        --nam_job "norm_sig"
    '''

  2. Compute log2 IP/input ratios from bedGraph files in serial.
    '''bash
    bash "\${HOME}/bin/execute_compute_signal.sh" \\
        --threads 1 \\
        --mode "ratio" \\
        --method "log2" \\
        --csv_fil_A "\${HOME}/project/norm/IP_1.bdg,\${HOME}/project/norm/IP_2.bdg" \\
        --csv_fil_B "\${HOME}/project/norm/in_1.bdg,\${HOME}/project/norm/in_2.bdg" \\
        --dir_out "\${HOME}/project/ratios" \\
        --typ_out "bg"
    '''

EOM
}
