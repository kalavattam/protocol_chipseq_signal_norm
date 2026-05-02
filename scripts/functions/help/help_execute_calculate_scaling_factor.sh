#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: help_execute_calculate_scaling_factor.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
#
# Distributed under the MIT license.


# shellcheck disable=SC2154
function help_execute_calculate_scaling_factor() {
    cat << EOM
Usage:
  execute_calculate_scaling_factor.sh
    [--help] [--verbose] [--dry_run]
    [--threads <int>]
    [--mode <enum:siq,spike>] [--method <enum:fractional,chiprx_alpha_ratio,chiprx_alpha_ip,chiprx_alpha_in,rxinput_alpha>]
    --csv_mip <csv:file> --csv_min <csv:file> [--csv_sip <csv:file>] [--csv_sin <csv:file>] [--aln_typ <enum:pe,se,auto>]
    --fil_out <file>
    [--tbl_met <file>] [--cfg_met <file>] [--eqn <enum:5,5nd,6,6nd>]
    [--len_def <int>] [--len_mip <csv:num>] [--len_min <csv:num>] [--dep_mip <csv:int>] [--dep_min <csv:int>] [--dep_sip <csv:int>] [--dep_sin <csv:int>]
    [--dp <int>]
    [--err_out <dir>] [--nam_job <str>] [--max_job <int>] [--slurm] [--time <time>]


Description:
  Coordinate calculation of siQ-ChIP or spike-in scaling factors for ChIP-seq data across one or more samples.

  In 'siq' mode, the script uses main-organism IP and input alignment files, together with a metadata table and YAML configuration file, to calculate siQ-ChIP alpha scaling factors via the downstream 'submit_calculate_scaling_factor.sh' wrapper and the Python helper scripts 'calculate_scaling_factor_siq_chip.py' and 'parse_metadata_siq_chip.py'. This workflow uses 'parse_metadata_siq_chip.yaml', which must be configured appropriately if input filenames do not follow the filename conventions described in the Tsukiyama Lab Bio-protocol manuscript.

  In 'spike' mode, the script uses main-organism and spike-in-organism IP and input alignment files to calculate spike-in scaling factors via the downstream 'submit_calculate_scaling_factor.sh' wrapper and the associated Python/helper-script workflow, including 'calculate_scaling_factor_spike.py' and supporting shell/Python utilities for obtaining fragment-length and alignment-depth values when needed.

  Jobs may be run through Slurm, GNU Parallel, or serial execution, depending on user arguments and the resolved number of jobs.


Arguments:
  -h, --hlp, --help  <flag>
    Display this help message and exit.

  -v, --verbose  <flag>
    Run in verbose mode.

  -dr, --dry, --dry_run  <flag>
    Print commands that would be executed without running them.

  -t, --thr, --threads  <int>
    Number of threads to use (default: ${threads}).

  -md, --mode  <enum:siq,spike>
    Scaling-factor mode to run: 'siq' or 'spike' (default: '${mode}').

  -me, --method  <enum:fractional,chiprx_alpha_ratio,chiprx_alpha_ip,chiprx_alpha_in,rxinput_alpha>
    Spike-in scaling method to compute when '--mode spike' is active (default if '--mode spike': 'fractional'; no default if '--mode siq').

    List of accepted canonical method names (first), aliases (subsequent), and calculations:
      - fractional | bioprotocol | bio_protocol
        (N_s^{in} / T^{in}) / (N_s^{IP} / T^{IP})
      - chiprx_alpha_ratio | alpha_chiprx_ratio | chiprx_ratio
        N_s^{in} / N_s^{IP}
      - chiprx_alpha_ip | alpha_chiprx_ip | chiprx_ip
        10^6 / N_s^{IP}
      - chiprx_alpha_in | alpha_chiprx_in | chiprx_in
        10^6 / N_s^{in}
      - rxinput_alpha | alpha_rxinput | rxi_alpha | alpha_rxi | rxinput | rxi
        (10^6 * N_s^{in}) / (N_s^{IP} * T^{in})

    Supported aliases are normalized internally.

  -at, --aln_typ, --align_typ  <enum:pe,paired,se,single,auto>
    Alignment layout type for input alignment files: 'pe' / 'paired', 'se' / 'single', or 'auto' (default: '${aln_typ}').

  -mp, --csv_mip  <csv:file>
    Comma-separated list of main-organism IP alignment files.

  -mn, --csv_min  <csv:file>
    Comma-separated list of main-organism input alignment files.

  -sp, --csv_sip  <csv:file>
    Comma-separated list of spike-in-organism IP alignment files. Required if '--mode spike'; ignored otherwise.

  -sn, --csv_sin  <csv:file>
    Comma-separated list of spike-in-organism input alignment files. Required if '--mode spike'; ignored otherwise.

  -tb, --tbl_met  <file>
    Tab-delimited siQ-ChIP metadata table. Required if '--mode siq'; ignored otherwise.

  -cm, --cfg_met  <file>
    YAML configuration file used to parse metadata in 'siq' mode.

  -eq, --eqn  <enum:5,5nd,6,6nd>
    Alpha equation to compute when '--mode siq' is active: '5', '5nd', '6', or '6nd' (default: '${eqn}'; ignored if '--mode spike').

    For descriptions of these equations, see Dickson et al., Sci Rep 2023 (PMID: 37160995). '5' corresponds to Equation 5 in the paper, and '6' corresponds to Equation 6.

    The 'nd' suffix denotes versions of those equations without depth terms (i.e., "no depth"), meaning forms that omit terms containing \hat{R} and/or \hat{R}_\mathrm{in}. Use the 'nd' versions when applying them to ratios of normalized coverage.

  -ld, --len_def  <int>
    Default fragment length to use for SE libraries when a per-file fragment length is not otherwise available.

  -lmp, --len_mip  <csv:num>
    Override fragment length(s) for main-organism IP files. May be a single broadcast value or a comma-separated list aligned to samples.

  -lmn, --len_min  <csv:num>
    Override fragment length(s) for main-organism input files. May be a single broadcast value or a comma-separated list aligned to samples.

  -dmp, --dep_mip  <csv:int>
    Override sequencing/alignment depth value(s) for main-organism IP files. May be a single broadcast value or a comma-separated list aligned to samples.

  -dmn, --dep_min  <csv:int>
    Override sequencing/alignment depth value(s) for main-organism input files. May be a single broadcast value or a comma-separated list aligned to samples.

  -dsp, --dep_sip  <csv:int>
    Override sequencing/alignment depth value(s) for spike-in-organism IP files. May be a single broadcast value or a comma-separated list aligned to samples. Used only in 'spike' mode.

  -dsn, --dep_sin  <csv:int>
    Override sequencing/alignment depth value(s) for spike-in-organism input files. May be a single broadcast value or a comma-separated list aligned to samples. Used only in 'spike' mode.

  -fo, --fil_out  <file>
    Tab-delimited output file to which scaling factors and related values are written.

  -dp, --dp, --rnd, --round, --decimals, --digits  <int>
    Number of decimal places used when rounding output values (default: ${dp}).

  -eo, --err_out  <dir>
    Directory in which wrapper-level stdout/stderr logs are written for local GNU Parallel or serial execution (default: '\$(dirname "\${fil_out}")/logs').

  -nj, --nam_job  <str>
    Job name prefix. If omitted, a mode-specific default is constructed:
      - "calc_sf_siq_\${eqn}" for 'siq' mode
      - "calc_sf_spike_\${method}" for 'spike' mode

  -mj, --max_job  <int>
    Maximum number of jobs to run concurrently (default: ${max_job}).

  -sl, --slurm  <flag>
    Submit jobs through Slurm.

  -tm, --time  <time>
    Slurm wall-clock time in 'h:mm:ss' format (default: '${time}'; used only if '--slurm' is active).


Dependencies:
  External programs:
    - awk
    - Bash >= 4.4
    - basename
    - dirname
    - GNU Parallel, when '--slurm' is not specified and multiple jobs are run
    - python
    - rm, when '--slurm' is not specified and multiple jobs are run
    - samtools
    - sbatch, when '--slurm'

  Shell scripts:
    - submit_calculate_scaling_factor.sh

  Python scripts:
    - calculate_scaling_factor_siq_chip.py, via submit_calculate_scaling_factor.sh
    - parse_metadata_siq_chip.py, via submit_calculate_scaling_factor.sh
    - calculate_scaling_factor_spike.py, via submit_calculate_scaling_factor.sh

  Configuration files:
    - parse_metadata_siq_chip.yaml, when '--mode siq'

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
      + check_arr_len_bcst
    - check_numbers.sh
      + check_arr_int_pos
      + check_arr_num_pos
      + check_format_time
      + check_int_nonneg
      + check_int_pos
    - format_outputs.sh
      + echo_err
      + print_banner_pretty
    - handle_env.sh
      + handle_env
    - help/help_execute_calculate_scaling_factor.sh
      + help_execute_calculate_scaling_factor
    - manage_parallel.sh
      + print_parallel_info
      + reset_max_job
      + set_params_parallel
    - wrap_cmd.sh
      + get_submit_logs
      + print_built_cmd


Notes:
  - In 'spike' mode, '--method' is normalized internally to one of the supported spike-in calculation modes.
  - In 'siq' mode, '--method' is not used.
  - For required per-sample input vectors, reconstructed arrays must be non-empty and of matching length.
  - Optional override vectors such as '--len_mip' or '--dep_min' may contain either:
    1. no value,
    2. one broadcast value, or
    3. one value per sample.
  - When '--slurm' is used, execution is parallelized via Slurm array tasks.
  - When '--slurm' is not used:
    + if the resolved local job count is greater than 1, execution uses GNU Parallel;
    + otherwise, execution is serial.
  - If '--dry_run' is enabled, commands are printed but not executed.


Examples:
  1. Compute spike-in scaling factors locally.
    '''bash
    execute_calculate_scaling_factor.sh
      --mode spike
      --method fractional
      --csv_mip IP1.sc.bam,IP2.sc.bam
      --csv_min in1.sc.bam,in2.sc.bam
      --csv_sip IP1.sp.bam,IP2.sp.bam
      --csv_sin in1.sp.bam,in2.sp.bam
      --fil_out scaling_factors.spike.tsv
    '''

  2. Compute siQ-ChIP alpha scaling factors.
    '''bash
    execute_calculate_scaling_factor.sh
      --mode siq
      --csv_mip IP1.sc.bam,IP2.sc.bam
      --csv_min in1.sc.bam,in2.sc.bam
      --tbl_met metadata.tsv
      --cfg_met parse_metadata_siq_chip.yaml
      --eqn 6nd
      --fil_out scaling_factors.siq.tsv
    '''

#TODO:
  - Add examples making use of optional arguments, including value overrides.
  - Hardcoded variable 'scr_hdr' is defined and validated, but not directly used in the top-level driver script:
    + remove 'scr_hdr' from this wrapper if unused or
    + pass it downstream if the submit script expects it.
EOM
}
