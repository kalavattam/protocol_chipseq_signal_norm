#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: help_execute_compute_signal.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
#
# Distributed under the MIT license.


usage=$(cat << EOM
Usage:
  execute_compute_signal.sh
    [--help] [--details] [--all_hlp] [--verbose] [--dry_run]
    [--threads <int>]
    [--mode <enum:signal,ratio,coord>] [--method <enum:unadj,frag,norm,log2,unadj_r,log2_r>]
    (--csv_infile <csv:file> | --csv_fil_A <csv:file> --csv_fil_B <csv:file>)
    --dir_out <dir> [--typ_out <enum:bedGraph|bedgraph|bdg|bg|bed[.gz]>] [--prefix <str>]
    [--siz_bin <int>] [--csv_usr_frg <csv:int>] [--csv_scl_fct <csv:num>]
    [--csv_dep_min <csv:num>] [--csv_pseudo <csv:spec>] [--eps <num>] [--skip_00 <enum:pre_scale|post_scale>] [--drp_nan] [--skp_pfx <csv:str>] [--track]
    [--dp <int>]
    [--err_out <dir>] [--nam_job <str>] [--max_job <int>] [--slurm] [--time <time>]

EOM
)

# shellcheck disable=SC2154
function help_execute_compute_signal() {
cat << EOM
${usage}


Description:
  Coordinate and automate the computation of signal tracks, ratio tracks, or fragment-coordinate files from BAM or bedGraph input files. Supports multiple normalization strategies and runs computations in serial or parallel via GNU Parallel or Slurm.

  For more details on what this script can do, including notes and usage examples, run either of the following:
    '''bash
    bash path/to/execute_compute_signal.sh --details
    bash path/to/execute_compute_signal.sh --all_hlp
    '''


Arguments:
  -h, --hlp, --help  <flag>
    Print this short help message and exit.

  -d, --det, --dets, --detail, --details  <flag>
    Print full documentation with notes and examples, then exit.

  -ah, --allhlp, --allhelp, --all_hlp, --all_help  <flag>
    Print both the short and full documentation, then exit.

  -v, --verbose  <flag>
    Run script in "verbose mode" (optional).

  -dr, --dry, --dry_run  <flag>
    Run script in "dry-run mode" (optional).

  -t, --thr, --threads  <int>
    Number of threads to use (default: ${threads}).

  -md, --mode  <enum:signal,ratio,coord>
    Type of computation to perform: 'signal', 'ratio', or 'coord' (default: '${mode}').

    See '--details' for full synonym lists.

  -me, --method  <enum:unadj,frag,norm,log2,unadj_r,log2_r>
    Signal or ratio computation subtype (used only with '--mode signal' or '--mode ratio'; default if '--mode signal': norm; default if '--mode ratio': unadj).

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

  -i, -fi, -ci, --infile, --infiles, --fil_in, --csv_infile, --csv_infiles  <csv:file>
    Comma-separated list of coordinate-sorted BAM files (used only with '--mode signal' or '--mode coord').

  -fA, -f1, -cA, -c1, --fil_A, --fil_1, --csv_A, --csv_1, --csv_fil_A  <csv:file>
    Comma-separated list of coordinate-sorted numerator bedGraph files (i.e. "file A"; e.g., for IP signal; used only with '--mode ratio').

  -fB, -f2, -cB, -c2, --fil_B, --fil_2, --csv_B, --csv_2, --csv_fil_B  <csv:file>
    Comma-separated list of coordinate-sorted denominator bedGraph files (i.e. "file B"; e.g., for input signal; used only with '--mode ratio').

  -do, --dir_out  <dir>
    Output directory for generated files.

  -to, --typ_out  <enum:bedGraph|bedgraph|bdg|bg|bed[.gz]>
    Format of output files (default: '${typ_out}').

    For '--mode signal' or '--mode ratio', the typical choice is a bedGraph-style track (e.g., 'bedGraph.gz' or 'bdg.gz').

    For '--mode coord', the typical choice is a BED-like coordinate file (e.g., 'bed.gz').

    If an incompatible combination is requested, the value is coerced to a sensible default with a warning; see '--details' for the full list and mode-specific behavior.

  -px, -pr, --pfx, --prfx, --prefix  <str>
    Custom prefix to prepend to output filenames.

  -sb, --siz_bin  <int>
    Bin size in base pairs for signal computation (used only with '--mode signal'; default: 10).

  -sf, --scale, --scl_fct, --csv_scl_fct  <csv:num>
    Comma-separated list of scaling factors or sentinels (optional; used only with '--mode signal' or '--mode ratio').

  -uf, --usr_frg, --csv_usr_frg  <csv:int>
    Comma-separated list of fragment lengths or sentinels (optional; used with '--mode signal' or '--mode coord').

  -dm, --dep_min, --depth_min, --csv_dep_min  <csv:num>
    Comma-separated list of minimum input depth values or sentinels (optional; used only with '--mode ratio').

  -ps, --pseudo, --pseudocount, --csv_pseudo  <csv:spec>
    Comma-separated list of per-sample pseudocount specs 'A[:B]' or sentinels (optional; used only with '--mode ratio').

  -e, --eps  <num>
    Shared epsilon or sentinel used for ratio-mode zero checks (optional; used only with '--mode ratio').

  -s0, --skp_00, --skip_00  <enum:pre_scale|post_scale>
    Shared zero-zero skip mode or sentinel for ratio computation: 'pre_scale' or 'post_scale' (optional; used only with '--mode ratio').

  -dn, --drp_nan, --drop_nan  <flag>
    Drop non-finite ratio rows ('inf', '-inf', and 'nan') from the main ratio output (optional; used only with '--mode ratio').

  -sk, --skp_pfx, --skip_pfx, --skip_prefix  <csv:str>
    Shared comma-separated list of bedGraph header prefixes or sentinel to skip while parsing ratio-mode input files (optional; used only with '--mode ratio').

  -tr, --trk, --track  <flag>
    If '--mode ratio', also write a companion bedGraph with all non-finite rows ('inf', '-inf', and 'nan') removed (optional).

  -dp, --dp, --rnd, --round, --decimals, --digits  <int>
    Maximum number of decimal places retained for finite emitted signal or ratio values (default: ${rnd}).

    After rounding, non-informative trailing zeros are stripped.

  -eo, --err_out  <dir>
    Directory for stderr and stdout TXT output files (default: '\${dir_out}/logs').

  -nj, --nam_job, --name_job  <str>
    Prefix for job names (default depends on resolved '--mode' and '--method'; e.g., 'compute_signal_norm', 'compute_ratio_unadj', or 'compute_coord').

  -mj, --max_job  <int>
    Maximum concurrent jobs to run (default: ${max_job}).

  -sl, --slurm  <flag>
    Submit jobs to the Slurm scheduler (optional).

  -tm, --time  <time>
    Slurm job time in 'h:mm:ss' format (required if '--slurm'; default: '${time}').


Notes:
  - BAM and bedGraph input files must be coordinate-sorted.
  - Input and output paths supplied to this wrapper interface must not contain spaces, commas, or semicolons.
  - See 'execute_compute_signal.sh --details' for more notes.
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

cat << EOM
Description:
  Driver script automating the computations of bedGraph signal or ratio tracks, or BED-like fragment coordinate files, from BAM (for signal tracks or fragment coordinate files) or bedGraph (for ratio tracks) input files.

  Supports multiple signal normalization strategies, including
    - unadjusted (raw) signal (i.e., per-bin totals with no fragment-length and/or library-size adjustments)
    - fragment-length adjusted signal (Dickson et al., JBC 2020; Dickson et al., Sci Rep 2023)
    - normalized coverage (Dickson et al., Sci Rep 2023)
    - siQ-ChIP IP efficiency (input-normalized ratio-based analyses; Dickson et al., JBC 2020; Dickson et al., Sci Rep 2023)
    - log2-transformed input-normalized signal ratios (e.g., as described in "Data Analysis G" of Alavattam et al., Bio-protocol 2025)
    - spike-in-normalized signal (i.e., unmodified ChIP-Rx; Orlando et al., Cell Rep 2014)
    - input- and spike-in-normalized signal ratios (modified ChIP-Rx as described in "Data Analysis I" of Alavattam et al., Bio-protocol 2025)

  Supports serial job execution and parallel job execution via Slurm or GNU Parallel.


Arguments:
  -h, --hlp, --help  <flag>
    Print a short help message and exit.

  -d, --det, --dets, --detail, --details  <flag>
    Print this full documentation with notes and examples, then exit.

  -ah, --allhlp, --allhelp, --all_hlp, --all_help  <flag>
    Print both the short and full documentation, then exit.

  -v, --verbose  <flag>
    Run script in 'verbose mode' (optional).

  -dr, --dry, --dry_run  <flag>
    Run script in 'dry-run mode' (optional).

  -t, --thr, --threads  <int>
    Number of threads to use (default: ${threads}).

  -md, --mode  <enum:signal,ratio,coord>
    Type of computation to perform: 'signal', 'ratio', or 'coord' (default: '${mode}'). Available options:
      - 's', 'sig', 'signal':
        + Compute signal tracks directly from BAM input files.
        + Supports unadjusted signal, fragment-length-adjusted signal, or normalized coverage.
        + See '--method' for calculation styles.
        + Use this option with appropriately computed scaling factors ('--scl_fct') to compute unmodified ChIP-Rx (Orlando et al., Cell Rep 2014) spike-in-normalized signal.
        + If 's' or 'sig' are supplied, variable 'mode' is set to "signal".

      - 'r', 'rat', 'ratio':
        + Compute IP/input signal ratios from bedGraph files.
        + Supports linear ratios, log2 ratios, and reciprocal variants of both.
        + Untransformed ratios are used in the modified ChIP-Rx spike-in and siQ-ChIP normalizations described in Alavattam et al., Bio-protocol 2025.
        + Use this option with '--method log2' to compute log2(IP/input) ratios.
        + If 'r' or 'rat' are supplied, variable 'mode' is set to "ratio".

      - 'c', 'coord', 'coordinates':
        + Instead of computing signal or ratio tracks, output fragment coordinates in BED-like format from BAM input files.
        + Use this option to prepare input files for the original siQ-ChIP implementation (Dickson et al., JBC 2020; Dickson et al., Sci Rep 2023).
        + For more details, see github.com/BradleyDickson/siQ-ChIP or github.com/kalavattam/siQ-ChIP.
        + This mode disables '--siz_bin' and '--method', and sets '--typ_out' to 'bed.gz' by default (or 'bed' if '--typ_out bed' is specified).
        + If 'c' or 'coordinates' are supplied, variable 'mode' is set to "coord".

  -me, --method  <enum:unadj,frag,norm,log2,unadj_r,log2_r>
    Signal or ratio computation subtype, used only with '--mode signal' or '--mode ratio' (default if '--mode signal': norm; default if '--mode ratio': unadj).
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

  -i, -fi, -ci, --infile, --infiles, --fil_in, --csv_infile, --csv_infiles  <csv:file>
    Comma-separated list of coordinate-sorted BAM infiles.

    Required when '--mode signal' or '--mode coord'. Ignored for '--mode ratio'.

  -fA, -f1, -cA, -c1, --fil_A, --fil_1, --csv_A, --csv_1, --csv_fil_A  <csv:file>
    Comma-separated list of coordinate-sorted numerator bedGraph files (i.e. "file A"; e.g., representing ChIP IP signal tracks).

    Required when '--mode ratio'; ignored otherwise.

    The list order must match that of '--fil_B' files.

  -fB, -f2, -cB, -c2, --fil_B, --fil_2, --csv_B, --csv_2, --csv_fil_B  <csv:file>
    Comma-separated list of coordinate-sorted denominator bedGraph files (i.e. "file B"; e.g., representing input signal tracks).

    Required when '--mode ratio'; ignored otherwise.

    The list order must match that of '--fil_A' files.

  -do, --dir_out  <dir>
    Output directory for generated files:
      - Signal tracks if '--mode signal'.
      - Ratio tracks if '--mode ratio'.
      - BED-like files of fragment coordinates if '--mode coord'.

  -to, --typ_out  <enum:bedGraph|bedgraph|bdg|bg|bed[.gz]>
    Format of signal track output files (default: '${typ_out}'). Available options:
      - 'bedGraph', 'bedgraph', 'bdg', 'bg':
        + Signal/ratio in bedGraph format.
        + Intended for '--mode signal' or '--mode ratio'; with '--mode coord' these values are accepted but coerced to 'bed.gz' (see Notes).

      - 'bedGraph.gz', 'bedgraph.gz', 'bdg.gz', 'bg.gz':
        + Signal/ratio in gzip-compressed bedGraph format.
        + Intended for '--mode signal' or '--mode ratio'; with '--mode coord' these values are accepted but coerced to 'bed.gz' (see Notes).

      - 'bed', 'bed.gz':
        + BED-like format for fragment coordinates instead of signal.
        + Intended for '--mode coord'; with '--mode signal' or '--mode ratio' these values are accepted but coerced to 'bdg.gz' (see Notes).

  -px, -pr, --pfx, --prfx, --prefix  <str>
    Custom prefix to prepend to output filenames.

    When '--mode signal' or '--mode coord':
      - If not specified, no prefix is added.
      - If specified, the prefix is prepended to the base filename, and any leading 'IP_' or 'in_' in that base filename is stripped before applying it.
      - Note: If you use the same '--prefix' for both IP and input BAMs in the same output directory, stripping 'IP_'/ 'in_' may cause filename collisions.
        + In that case, prefer distinct prefixes (e.g., 'IP', 'in').

    When '--mode ratio':
      - If not specified, a default prefix is automatically constructed based on '--method' and '--scl_fct'; for example:
        + 'rat' (default)
        + 'log2_rat' (if '--method log2')
        + 'recip_rat' (if '--method unadj_r')
        + 'log2_recip_rat' (if '--method log2_r')
        + 'scl_rat' (if '--method unadj' and '--scl_fct' is supplied)
        + 'scl_log2_rat' (if '--method log2' and '--scl_fct' is supplied)
        + 'scl_recip_rat' (if '--method unadj_r' and '--scl_fct' is supplied)
        + 'scl_log2_recip_rat' (if '--method log2_r' and '--scl_fct' is supplied)
      - If specified, the custom prefix replaces the default.
      - Whether specified or not, any leading 'IP_' string in the base name is stripped before the prefix.

  -sb, --siz_bin  <int>
    Bin size in base pairs for signal computation.

    Used only with '--mode signal' (default: 10); ignored otherwise.

  -sf, --scale, --scl_fct, --csv_scl_fct  <csv:num>
    Comma-separated list of scaling factors or sentinels to apply to signal or ratio values.

    Used with either '--mode signal' or '--mode ratio'; ignored otherwise.

    List size must match the number of input files via '--csv_infile' or '--csv_fil_A'/'--csv_fil_B'.

  -uf, --usr_frg, --csv_usr_frg  <csv:int>
    Comma-separated list of fragment lengths or sentinels to use instead of read lengths (single-end alignments) or template lengths (paired-end alignments; optional).

    Used with either '--mode signal' or '--mode coord'; ignored otherwise.

    List size must match the number of input files via '--csv_infile'.

  -dm, --dep_min, --depth_min, --csv_dep_min  <csv:num>
    Comma-separated list of minimum input depth values or sentinels used to avoid extreme division operations (optional).

    Used only with '--mode ratio'; ignored otherwise.

    List size must match the number of input files via '--csv_fil_A'/'--csv_fil_B'.

    Although allowed, using '--dep_min' together with '--pseudo' is usually harder to interpret, since both stabilize low-depth ratio behavior in different ways.

  -ps, --pseudo, --pseudocount, --csv_pseudo  <csv:spec>
    Comma-separated list of pseudocount specifications or sentinels used during ratio computation (optional).

    Used only with '--mode ratio'; ignored otherwise.

    List size must match the number of input files via '--csv_fil_A'/'--csv_fil_B'.

    Each non-sentinel element may be either:
      - 'A'    Add pseudocount A symmetrically.
      - 'A:B'  Add pseudocount A to file A and B to file B.

    Although allowed, using '--pseudo' together with '--dep_min' is usually harder to interpret, since both stabilize low-depth ratio behavior in different ways.

  -e, --eps  <num>
    Shared epsilon or sentinel used for ratio-mode zero checks (optional).

    Used only with '--mode ratio'; ignored otherwise.

    Non-sentinel values must be non-negative floats.

  -s0, --skp_00, --skip_00  <enum:pre_scale|post_scale>
    Shared zero-zero skip mode or sentinel for ratio computation (optional).

    Used only with '--mode ratio'; ignored otherwise.

    Non-sentinel values must be one of 'pre_scale' or 'post_scale'.

  -dn, --drp_nan, --drop_nan  <flag>
    Drop non-finite values from the main ratio output (optional).

    Used only with '--mode ratio'; ignored otherwise.

    If set, rows yielding 'inf', '-inf', or 'nan' are omitted from the main ratio output.

  -sk, --skp_pfx, --skip_pfx, --skip_prefix  <csv:str>
    Shared comma-separated list of bedGraph header prefixes or sentinel to skip (optional).

    Used only with '--mode ratio'; ignored otherwise.

    Passed through to 'submit_compute_signal.sh' and then to 'compute_signal_ratio.py'.

  -tr, --trk, --track  <flag>
    If '--mode ratio', also write a companion bedGraph with all non-finite rows ('inf', '-inf', and 'nan') removed (optional).

    The new file will include '.track' before the extension.

    This cleaned version is ideal for visualization in genome browsers such as IGV, avoiding issues caused by 'inf', '-inf' or 'nan' values.

  -dp, --dp, --rnd, --round, --decimals, --digits  <int>
    Maximum number of decimal places retained for finite emitted signal or ratio values (default: ${rnd}).

    After rounding, non-informative trailing zeros are stripped.

  -eo, --err_out  <dir>
    Directory for stderr and stdout TXT output files (default: '\${dir_out}/logs').

  -nj, --nam_job, --name_job  <str>
    Prefix for job names (default depends on resolved '--mode' and '--method'; e.g., 'compute_signal_norm', 'compute_ratio_log2', or 'compute_coord').

  -mj, --max_job  <int>
    Maximum concurrent jobs to run (default: ${max_job}).

    With '--slurm': maximum number of Slurm array tasks allowed to run concurrently.

    Without '--slurm':
      - If the resolved number of parallel jobs is greater than 1, jobs run in parallel via GNU Parallel.
      - If the resolved number of parallel jobs is 1, jobs run serially.

  -sl, --slurm  <flag>
    Submit jobs to the Slurm scheduler (optional).

    If '--slurm' is not specified, this script uses '--threads', '--max_job', and the detected CPU core count (via 'set_params_parallel') to decide whether to run jobs in parallel with GNU Parallel or in serial.
      - In this non-Slurm path, '--threads' is treated as a total CPU/thread budget for the local machine, and 'set_params_parallel' converts that into threads per job plus number of parallel jobs.
      - If the resolved number of parallel jobs is greater than 1, jobs are run with GNU Parallel.
      - If the resolved number of parallel jobs is 1, jobs are run serially.

  -tm, --time  <time>
    The length of time, in 'h:mm:ss' format, for the Slurm job (required if '--slurm' is specified, ignored if not; default: '${time}').


Dependencies:
  External programs:
    - Bash >= 4.4
    - basename
    - dirname
    - GNU Parallel, when '--slurm' is not specified and multiple jobs are run
    - python
    - rm, when '--slurm' is not specified and multiple jobs are run
    - sbatch, when '--slurm'
    - tr

  Shell scripts:
    - submit_compute_signal.sh

  Sourced function scripts:
    - source_helpers.sh
      + source_helpers
    - check_args.sh
      + require_optarg
      + check_str_delim
    - check_env.sh
      + check_env_installed
      + check_pgrm_path
    - check_inputs.sh
      + validate_var
      + validate_var_dir
      + validate_var_file
      + check_arr_nonempty
      + check_arr_lengths
    - check_numbers.sh
      + check_flt_nonneg
      + check_flt_pos
      + check_format_time
      + check_int_pos
      + check_int_nonneg
    - format_outputs.sh
      + echo_err
      + echo_warn
      + print_banner_pretty
      + summarize_sig_norm
    - handle_env.sh
      + handle_env
    - help/help_execute_compute_signal.sh
      + help_execute_compute_signal
      + detail_execute_compute_signal
    - manage_parallel.sh
      + print_parallel_info
      + reset_max_job
      + set_params_parallel
    - populate_array_empty.sh
      + populate_array_empty
    - wrap_cmd.sh
      + get_submit_logs
      + print_built_cmd


Notes:
  - BAM and bedGraph input files must be coordinate-sorted.
  - Input and output paths supplied to this wrapper interface must not contain spaces, commas, or semicolons.
    + Commas are used internally as list delimiters.
    + Semicolons are also considered unsafe in this wrapper workflow.
    + Spaces in paths are not supported by the current list-serialization and reconstruction logic.
  - If applicable, use consistent file ordering between IP and input files.
  - '--typ_out' must be compatible with the selected '--mode'.
    + With '--mode signal' or '--mode ratio', bedGraph-style values ('bedGraph', 'bedgraph', 'bdg', and 'bg', and their '.gz' variants) are allowed; 'bed'/'bed.gz' are accepted but are automatically converted to 'bdg.gz' with a warning.
    + With '--mode coord', 'bed'/'bed.gz' are allowed; bedGraph-style values are accepted but are automatically converted to 'bed.gz' with a warning.
  - Output filenames are derived from BAM or bedGraph input files and the value associated with '--typ_out'.
  - For bedGraph-style output, '--rnd' sets the maximum number of decimal places retained for finite emitted values; after rounding, non-informative trailing zeros and any trailing decimal point are stripped.
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


Examples:
  1. Compute normalized coverage using GNU Parallel.
    '''bash
    bash "\${HOME}/scripts/execute_compute_signal.sh"
      --threads 8
      --mode "signal"
      --method "norm"
      --csv_infile "\${HOME}/project/samples/sample_1.bam,\${HOME}/project/samples/sample_2.bam"
      --dir_out "\${HOME}/project/tracks"
      --typ_out "bdg.gz"
      --siz_bin 50
      --err_out "\${HOME}/project/logs"
      --nam_job "norm_sig"
    '''

  2. Compute log2 IP/input ratios from bedGraph files in serial.
    '''bash
    bash "\${HOME}/scripts/execute_compute_signal.sh"
      --threads 1
      --mode "ratio"
      --method "log2"
      --csv_fil_A "\${HOME}/project/norm/IP_1.bdg,\${HOME}/project/norm/IP_2.bdg"
      --csv_fil_B "\${HOME}/project/norm/in_1.bdg,\${HOME}/project/norm/in_2.bdg"
      --dir_out "\${HOME}/project/ratios"
      --typ_out "bg"
    '''


#TODO:
  - More examples.
  - Add PMIDs for all studies mentioned.
EOM
}
