#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: help_execute_calculate_scaling_factor.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
# GPT-5.6) were used in development and documentation.
#
# Distributed under the MIT license.


# TODO: Add examples that use optional arguments, including value overrides.
# shellcheck disable=SC2154
function help_execute_calculate_scaling_factor() {
    cat << EOM
Usage
-----
  execute_calculate_scaling_factor.sh
    [--help] [--verbose] [--dry_run]
    [--env_nam <str>] [--threads <int>]
    [--mode <mode>] [--method <method>]
    [--aln_typ <layout>] [--ref_fa <file>]
    --csv_mip <csv> --csv_min <csv> [--csv_sip <csv>] [--csv_sin <csv>]
    --fil_out <file> [--force] [--no_parts] [--no_header]
    [--tbl_met <file>] [--cfg_met <file>] [--eqn <equation>]
    [--len_def <int>] [--len_mip <csv>] [--len_min <csv>] [--dep_mip <csv>] [--dep_min <csv>] [--dep_sip <csv>] [--dep_sin <csv>]
    [--dp <int>]
    [--dir_eo <dir>] [--nam_job <str>]
    [--max_job <int>] [--slurm] [--time <time>]


  Coordinate calculation of siQ-ChIP or spike-in scaling factors for ChIP-seq data across one or more samples.

  In 'siq' mode, the script uses main-organism IP and input alignment files, together with a metadata table and YAML configuration file, to calculate siQ-ChIP alpha scaling factors via the downstream 'submit_calculate_scaling_factor.sh' wrapper and the Python helper scripts 'calculate_scaling_factor_siqchip.py' and 'parse_metadata_siqchip.py'. This workflow uses 'parse_metadata_siqchip.yml', which must be configured appropriately if input filenames do not follow the filename conventions described in the Tsukiyama Lab Bio-protocol manuscript.

  In 'spike' mode, the script uses main-organism and spike-in-organism IP and input alignment files to calculate spike-in scaling factors via the downstream 'submit_calculate_scaling_factor.sh' wrapper and the associated Python/helper-script workflow, including 'calculate_scaling_factor_spike.py' and supporting shell/Python utilities for obtaining fragment-length and alignment-depth values when needed.

  Jobs may be run through Slurm, GNU Parallel, or serial execution, depending on user arguments and the resolved number of jobs. After successful worker completion, the script combines deterministic per-sample part files into the requested final TSV, then writes the mode-specific header unless '--no_header' is supplied.


Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  -v, --verbose : flag
    Run script in verbose mode.

  -dr, --dry, --dry_run : flag
    Run script in dry-run mode. Print commands that would be executed without running them.

  -en, --env, --env_nam : str
    Conda environment to activate. (default: '${env_nam}').

  -t, --thr, --threads : int
    Number of threads to use (default: ${threads}).

  -md, --mode : {'siq', 'spike'}
    Workflow mode. Scaling-factor mode to run: 'siq' or 'spike' (default: '${mode}').

  -me, --method : {'fractional', 'chiprx_alpha_ratio', 'chiprx_alpha_ip', 'chiprx_alpha_in', 'rxinput_alpha'}
    Workflow method. Spike-in scaling method to compute when '--mode spike' is active (default if '--mode spike': 'chiprx_alpha_ratio'; no default if '--mode siq').

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

  -at, --align_typ, --aln_typ : {'pe', 'paired', 'se', 'single', 'auto'}
    Alignment layout type for input alignment files: 'pe', 'paired', 'se', 'single', or 'auto' (default: '${aln_typ}').

  -rf, --ref_fa : file
    Reference FASTA file. Reference FASTA required when any input alignment file is CRAM.

  -cmip, --csv_mip : list of file
    Comma-separated list of main IP alignment files.

  -cmin, --csv_min : list of file
    Comma-separated list of main input alignment files.

  -csip, --csv_sip : list of file
    Comma-separated list of spike-in IP alignment files. Required if '--mode spike'; ignored otherwise.

  -csin, --csv_sin : list of file
    Comma-separated list of spike-in input alignment files. Required if '--mode spike'; ignored otherwise.

  -fo, --fil_out : file
    Output file path. Tab-delimited output file to which scaling factors and related values are written.

  -f, --force : flag
    Replace an existing final TSV after successful part-file assembly.

  -np, --no_parts : flag
    Remove per-sample part files after successful final-table assembly.

  -nh, --no_header : flag
    Leave the final TSV data-only instead of prepending the mode-specific header.

  -tb, --tbl_met : file
    siQ-ChIP metadata table. Required if '--mode siq'; ignored otherwise.

  -cm, --cfg_met : file
    YAML configuration file for metadata parsing. Used in 'siq' mode.

  -eq, --eqn : {'5', '5nd', '6', '6nd'}
    siQ-ChIP alpha equation to compute when '--mode siq' is active: '5', '5nd', '6', or '6nd' (default: '${eqn}'; ignored if '--mode spike').

    For descriptions of these equations, see Dickson et al., Sci Rep 2023 (PMID: 37160995). '5' corresponds to Equation 5 in the paper, and '6' corresponds to Equation 6.

    The 'nd' suffix denotes versions of those equations without depth terms (i.e., 'no depth'), meaning forms that omit terms containing \hat{R} and/or \hat{R}_\mathrm{in}. Use the 'nd' versions when applying them to ratios of normalized coverage.

  -ld, --len_def : int
    Default fragment length for single-end libraries when a per-file fragment length is not otherwise available.

  -lmp, --len_mip : list of number
    Fragment length value(s) for main IP alignment files. May be a single broadcast value or a comma-separated list aligned to samples.

  -lmn, --len_min : list of number
    Fragment length value(s) for main input alignment files. May be a single broadcast value or a comma-separated list aligned to samples.

  -dmp, --dep_mip : list of int
    Sequencing/alignment depth value(s) for main IP alignment files. May be a single broadcast value or a comma-separated list aligned to samples.

  -dmn, --dep_min : list of int
    Override sequencing/alignment depth value(s) for main-organism input files. May be a single broadcast value or a comma-separated list aligned to samples.

  -dsp, --dep_sip : list of int
    Sequencing/alignment depth value(s) for spike-in IP alignment files. May be a single broadcast value or a comma-separated list aligned to samples. Used only in 'spike' mode.

  -dsn, --dep_sin : list of int
    Sequencing/alignment depth value(s) for spike-in input alignment files. May be a single broadcast value or a comma-separated list aligned to samples. Used only in 'spike' mode.

  -dp, --dp : int
    Maximum number of decimal places retained for finite emitted values (default: ${dp}).

  -deo, --dir_eo : dir
    Directory for stderr and stdout log files written for local GNU Parallel or serial execution (default: '\$(dirname "\${fil_out}")/logs').

  -nj, --nam_job : str
    Job name prefix. If omitted, a mode-specific default is constructed:
      - 'calc_sf_siq_\${eqn}' for 'siq' mode
      - 'calc_sf_spike_\${method}' for 'spike' mode

  -mj, --max_job : int
    Maximum number of jobs to run concurrently (default: ${max_job}).

  -sl, --slurm : flag
    Submit jobs to the Slurm scheduler. Submit jobs through Slurm.

  -tm, --time : time
    Slurm job time limit. Slurm wall-clock time in 'h:mm:ss' format (default: '${time}'; used only if '--slurm' is active).


Notes
-----
  Runtime requirements:
    - A compatible Conda environment providing the listed tools
    - awk
    - basename
    - bash >= 4.4
    - conda (when the requested environment is not active)
    - dirname
    - parallel (when '--slurm' is not specified and multiple jobs are run)
    - python >= 3.11
    - rm (when '--slurm' is not specified and multiple jobs are run)
    - samtools
    - sbatch (when '--slurm' is specified)

  - In 'spike' mode, '--method' is normalized internally to one of the supported spike-in calculation modes.
  - In 'siq' mode, '--method' is not used.
  - CRAM inputs require '--ref_fa'.
  - For required per-sample input vectors, reconstructed arrays must be non-empty and of matching length.
  - Optional override vectors such as '--len_mip' or '--dep_min' may contain either:
    1. no value,
    2. one broadcast value, or
    3. one value per sample.
  - When '--slurm' is used, execution is parallelized via Slurm array tasks.
    + A dependent combiner job is submitted with 'afterok:<array_job_id>'.
    + Unless '--no_header' is supplied, a dependent header job is submitted with 'afterok:<combiner_job_id>'.
  - When '--slurm' is not used:
    + if the resolved local job count is greater than 1, execution uses GNU Parallel;
    + otherwise, execution is serial.
  - The final TSV is assembled automatically after all workers finish.
  - The final TSV receives a mode-specific header unless '--no_header' is supplied.
  - Existing final TSVs are replaced only when '--force' is supplied.
  - Per-sample part files are retained unless '--no_parts' is supplied.
  - Compute downstream denominator floors separately with 'compute_input_floor'.
  - If '--dry_run' is enabled, commands are printed but not executed.


Examples
--------
  1. Compute spike-in scaling factors locally.
    '''bash
    execute_calculate_scaling_factor.sh \\
        --mode spike \\
        --method fractional \\
        --csv_mip IP1.sc.bam,IP2.sc.bam \\
        --csv_min in1.sc.bam,in2.sc.bam \\
        --csv_sip IP1.sp.bam,IP2.sp.bam \\
        --csv_sin in1.sp.bam,in2.sp.bam \\
        --fil_out scaling_factors.spike.tsv
    '''

  2. Compute siQ-ChIP alpha scaling factors.
    '''bash
    execute_calculate_scaling_factor.sh \\
        --mode siq \\
        --csv_mip IP1.sc.bam,IP2.sc.bam \\
        --csv_min in1.sc.bam,in2.sc.bam \\
        --tbl_met metadata.tsv \\
        --cfg_met parse_metadata_siqchip.yml \\
        --eqn 6nd \\
        --fil_out scaling_factors.siq.tsv
    '''
EOM
}
