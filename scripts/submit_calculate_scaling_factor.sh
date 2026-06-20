#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: submit_calculate_scaling_factor.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models) were used in
# development.
#
# Distributed under the MIT license.


#  Require Bash >= 4.4 before doing any work
if [[ -z "${BASH_VERSION:-}" ]]; then
    echo "error(shell):" \
        "this script must be run under Bash >= 4.4." >&2
    exit 1
elif ((
    BASH_VERSINFO[0] < 4 || ( BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] < 4 )
)); then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "this script requires Bash >= 4.4; current version is" \
        "'${BASH_VERSION}'." >&2
    exit 1
fi

#  Run in safe mode, exiting on errors, unset variables, and pipe failures
set -euo pipefail


#  Define script-specific functions
function derive_samp_sf() {
    local fil_aln="${1:-}"
    local samp=""
    local show_help

    show_help=$(cat << EOM
Usage:
  derive_samp_sf [--help] fil_aln

Description:
  Derive a sample name from a BAM or CRAM path for logging.

Positional argument:
  1  fil_aln  <str>  Alignment path from which to derive the sample name.

Returns:
  Prints the derived sample name to stdout; returns 1 if argument parsing fails.
EOM
    )

    if [[ "${fil_aln}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${fil_aln}" ]]; then
        echo "error(derive_samp_sf): requires an alignment path." >&2
        return 1
    fi

    samp="$(basename "${fil_aln}")"
    samp="${samp%.bam}"
    samp="${samp%.cram}"
    samp="${samp#IP_}"
    samp="${samp#in_}"
    samp="${samp//./_}"

    echo "${samp}"
}


#  Resolve '--dir_scr' before sourced parser helpers are available
function resolve_dir_scr() {
    local scr="${1:-}"
    shift

    local -a args=( "$@" )
    local i=0

    if [[ -z "${scr}" ]]; then
        scr="unknown_script"
    fi

    for (( i = 0; i < ${#args[@]}; i++ )); do
        case "${args[i]}" in
            -ds|--dir[_-]scr)
                if (( i + 1 >= ${#args[@]} )) \
                    || [[ -z "${args[i + 1]:-}" || "${args[i + 1]}" == -* ]]
                then
                    echo "error(${scr}):" \
                        "option '${args[i]}' requires a value." >&2
                    echo >&2
                    show_help_main
                    return 1
                fi

                printf "%s\n" "${args[i + 1]}"
                return 0
                ;;
        esac
    done

    echo "error(${scr}):" \
        "required option '--dir_scr' was not supplied." >&2
    echo >&2
    show_help_main
    return 1
}


#  Source 'source_helpers.sh' and requested helper scripts from 'dir_scr'
function source_helpers_submit() {
    local scr="${1:-}"
    local dir_scr_arg="${2:-}"
    local fnc_src

    if (( $# < 2 )); then
        echo "error(${scr:-unknown_script}):" \
            "expected at least two arguments: 'script' and 'dir_scr_arg'." >&2
        return 1
    fi

    shift 2

    if [[ -z "${scr}" ]]; then
        scr="unknown_script"
    fi

    if [[ -z "${dir_scr_arg}" ]]; then
        echo "error(${scr}):" \
            "positional argument 2, 'dir_scr_arg', is missing." >&2
        return 1
    elif [[ ! -d "${dir_scr_arg}" ]]; then
        echo "error(${scr}):" \
            "script directory not found: '${dir_scr_arg}'." >&2
        return 1
    elif (( $# < 1 )); then
        echo "error(${scr}):" \
            "at least one helper script name must be supplied." >&2
        return 1
    fi

    dir_fnc="${dir_scr_arg}/functions"
    fnc_src="${dir_fnc}/source_helpers.sh"

    if [[ ! -f "${fnc_src}" ]]; then
        echo "error(${scr}):" \
            "script not found: '${fnc_src}'." >&2
        return 1
    fi

    # shellcheck disable=SC1090
    source "${fnc_src}" || {
        echo "error(${scr}):" \
            "failed to source '${fnc_src}'." >&2
        return 1
    }

    source_helpers "${dir_fnc}" "$@" || {
        echo "error(${scr}):" \
            "failed to source required helper scripts." >&2
        return 1
    }
}


#  Parse keyword arguments after helper scripts have been sourced
function parse_args() {
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
            -en|--env|--env[_-]nam)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                env_nam="${2}"
                shift 2
                ;;

            -ds|--dir[_-]scr)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                dir_scr="${2}"
                shift 2
                ;;

            -t|--thr|--threads)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                threads="${2}"
                shift 2
                ;;

            -md|--mode)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                mode="$(printf '%s\n' "${2}" | tr '[:upper:]' '[:lower:]')"
                shift 2
                ;;

            -me|--method)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                method="$(printf '%s\n' "${2}" | tr '[:upper:]' '[:lower:]')"
                shift 2
                ;;

            -mp|--csv[_-]mip)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                csv_mip="${2}"
                shift 2
                ;;

            -mn|--csv[_-]min)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                csv_min="${2}"
                shift 2
                ;;

            -sp|--csv[_-]sip)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                csv_sip="${2}"
                shift 2
                ;;

            -sn|--csv[_-]sin)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                csv_sin="${2}"
                shift 2
                ;;

            -at|--aln[_-]typ|--align[_-]typ)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                aln_typ="$(printf '%s\n' "${2}" | tr '[:upper:]' '[:lower:]')"
                shift 2
                ;;

            -r|--ref|--ref[_-]fa|--reference)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                ref_fa="${2}"
                shift 2
                ;;

            -fo|--fil[_-]out)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                fil_out="${2}"
                shift 2
                ;;

            -io|--idx[_-]out)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                idx_out="${2}"
                shift 2
                ;;

            -tb|--tbl[_-]met)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                tbl_met="${2}"
                shift 2
                ;;

            -cm|--cfg[_-]met)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                cfg_met="${2}"
                shift 2
                ;;

            -eq|--eqn)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                eqn="$(printf '%s\n' "${2}" | tr '[:upper:]' '[:lower:]')"
                shift 2
                ;;

            -ld|--len[_-]def)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                len_def="${2}"
                shift 2
                ;;

            -lmp|--len[_-]mip)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                len_mip="${2}"
                shift 2
                ;;

            -lmn|--len[_-]min)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                len_min="${2}"
                shift 2
                ;;

            -dmp|--dep[_-]mip)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                dep_mip="${2}"
                shift 2
                ;;

            -dmn|--dep[_-]min)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                dep_min="${2}"
                shift 2
                ;;

            -dsp|--dep[_-]sip)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                dep_sip="${2}"
                shift 2
                ;;

            -dsn|--dep[_-]sin)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                dep_sin="${2}"
                shift 2
                ;;

            -dp|--dp|--rnd|--round|--decimals|--digits)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                rnd="${2}"
                shift 2
                ;;

            -eo|--err[_-]out)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                err_out="${2}"
                shift 2
                ;;

            -nj|--nam[_-]job)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                nam_job="${2}"
                shift 2
                ;;

            *)
                echo_err "unknown parameter passed: '${1}'."
                echo >&2
                show_help_main
                return 1
                ;;
        esac
    done
}


#  Canonicalize mode, method, alignment type, coefficient name, and job name
function canonicalize_args() {
    case "${mode}" in
        siq)
            mode="siq"

            if [[ -n "${method}" ]]; then
                echo_err "'--method' is valid only when '--mode spike'."
                return 1
            fi

            method=""
            ;;

        spike|spk)
            mode="spike"

            if [[ -z "${method}" ]]; then
                method="chiprx_alpha_ratio"
            fi

            case "${method}" in
                fractional|bioprotocol|bio_protocol|s)
                    method="fractional"
                    ;;
                chiprx_alpha_ratio|alpha_chiprx_ratio|chiprx_ratio|r)
                    method="chiprx_alpha_ratio"
                    ;;
                chiprx_alpha_ip|alpha_chiprx_ip|chiprx_ip)
                    method="chiprx_alpha_ip"
                    ;;
                chiprx_alpha_in|alpha_chiprx_in|chiprx_in)
                    method="chiprx_alpha_in"
                    ;;
                rxinput_alpha|alpha_rxinput|rxi_alpha|alpha_rxi|rxinput|rxi)
                    method="rxinput_alpha"
                    ;;
                *)
                    echo_err "invalid '--method' value: '${method}'."
                    return 1
                    ;;
            esac
            ;;

        *)
            echo_err "'--mode' must be 'siq' or 'spike', but got '${mode}'."
            return 1
            ;;
    esac

    case "${aln_typ}" in
        pe|paired) aln_typ="pe" ;;
        se|single) aln_typ="se" ;;
        auto)      : ;;
        *)
            echo_err \
                "'--aln_typ' must be 'pe', 'se', or 'auto', but got" \
                "'${aln_typ}'."
            return 1
            ;;
    esac

    case "${eqn}" in
        5|5nd|6|6nd) : ;;
        *)
            echo_err \
                "'--eqn' must be '5', '5nd', '6', or '6nd', but got" \
                "'${eqn}'."
            return 1
            ;;
    esac

    if [[ -z "${nam_job}" ]]; then
        if [[ "${mode}" == "siq" ]]; then
            nam_job="calc_sf_siq_${eqn}"
        else
            nam_job="calc_sf_spike_${method}"
        fi
    fi

    if [[ "${mode}" == "spike" ]]; then
        coef_spk="${method}"
    else
        coef_spk=""
    fi
}


#  Validate required arguments and simple numeric values
function validate_args() {
    validate_var     "env_nam" "${env_nam}"         || return 1
    validate_var_dir "dir_scr" "${dir_scr}" 0 false || return 1
    validate_var     "threads" "${threads}"         || return 1
    check_int_pos    "${threads}" "threads"         || return 1
    check_int_pos    "${rnd}"     "rnd"             || return 1

    validate_var "csv_mip" "${csv_mip}" || return 1
    validate_var "csv_min" "${csv_min}" || return 1
    validate_var "fil_out" "${fil_out}" || return 1

    if [[ "${fil_out}" == "-" ]]; then
        echo_err "'-' is not allowed for '--fil_out' in this wrapper."
        return 1
    fi

    validate_var_dir "fil_out parent directory" "$(dirname "${fil_out}")" \
        || return 1

    if [[ -n "${idx_out}" ]]; then
        check_int_nonneg "${idx_out}" "idx_out" || return 1
    fi

    if [[ "${mode}" == "spike" ]]; then
        validate_var "csv_sip" "${csv_sip}" || return 1
        validate_var "csv_sin" "${csv_sin}" || return 1
    fi

    if [[ "${mode}" == "siq" ]]; then
        validate_var_file "tbl_met" "${tbl_met}" || return 1
        validate_var_file "cfg_met" "${cfg_met}" || return 1

        if [[ "${aln_typ}" == "se" ]]; then
            if [[
                     -z "${len_def}" \
                && ( -z "${len_mip}" || -z "${len_min}" )
            ]]; then
                echo_err \
                    "for '--mode siq' with SE data, supply '--len_def' or" \
                    "both '--len_mip' and '--len_min'."
                return 1
            fi
        fi
    fi

    if [[ -n "${len_def}" ]]; then
        check_int_pos "${len_def}" "len_def" || return 1
    fi

    validate_var "nam_job" "${nam_job}" || return 1

    if [[ -z "${err_out}" ]]; then
        echo_err "'--err_out' is required."
        return 1
    fi

    validate_var_dir "err_out" "${err_out}" || return 1
}


#  Resolve script paths after 'dir_scr' has been validated and helpers sourced
function resolve_script_paths() {
    validate_var_dir "dir_scr" "${dir_scr}" 0 false || return 1

    dir_rep="$(cd "${dir_scr}/.." && pwd)"

    scr_met="${dir_scr}/parse_metadata_siqchip.py"
    scr_siq="${dir_scr}/calculate_scaling_factor_siqchip.py"
    scr_spk="${dir_scr}/calculate_scaling_factor_spike.py"

    validate_var_file "scr_met" "${scr_met}" || return 1
    validate_var_file "scr_siq" "${scr_siq}" || return 1
    validate_var_file "scr_spk" "${scr_spk}" || return 1
}


#  Initialize hardcoded argument variables
function init_args_hardcoded() {
    #  If true, print verbose/debug Bash-level logging
    debug=true

    #  If true, parse arguments and exit before validation or execution
    p_only=false

    #  If true, parse and check arguments, then exit before execution
    pc_only=false
}


#  Initialize argument variables, assigning default values where applicable
function init_arg_defs() {
    env_nam="env_protocol"
    dir_scr=""
    threads=1

    mode="spike"    # siq | spike
    method=""       # chiprx_alpha_ratio | fractional | rxinput_alpha | ...

    csv_mip=""
    csv_min=""
    csv_sip=""
    csv_sin=""
    aln_typ="auto"  # paired | pe | single | se | auto
    ref_fa=""
    fil_out=""
    idx_out=""

    tbl_met=""
    cfg_met=""
    eqn="6nd"

    len_def=""      # SE fallback length (bp)
    len_mip=""      # Precomputed lengths: main IP
    len_min=""      # Precomputed lengths: main input
    dep_mip=""      # Precomputed counts: main IP
    dep_min=""      # Precomputed counts: main input
    dep_sip=""      # Precomputed counts: spike IP
    dep_sin=""      # Precomputed counts: spike input

    rnd=24
    err_out=""
    nam_job=""

    unset arr_mip arr_min arr_sip arr_sin
    unset arr_len_mip arr_len_min
    unset arr_dep_mip arr_dep_min arr_dep_sip arr_dep_sin

    declare -ga arr_mip arr_min arr_sip arr_sin
    declare -ga arr_len_mip arr_len_min
    declare -ga arr_dep_mip arr_dep_min arr_dep_sip arr_dep_sin
}


#  Initialize hardcoded arguments and user-facing argument defaults
function init_defs() {
    init_args_hardcoded
    init_arg_defs
}

#  Define the help message
function show_help_main() {
    cat << EOM >&2
Usage:
  submit_calculate_scaling_factor.sh
    [--help] [--env_nam <str>] --dir_scr <str> [--threads <int>]
    [--mode {siq,spike}] [--method {fractional,chiprx_alpha_ratio,chiprx_alpha_ip,chiprx_alpha_in,rxinput_alpha}]
    --csv_mip <csv:file> --csv_min <csv:file> [--csv_sip <csv:file>] [--csv_sin <csv:file>] [--aln_typ {pe,se,auto}] [--ref_fa <file>] --fil_out <str> [--idx_out <int>]
    [--tbl_met <str>] [--cfg_met <str>] [--eqn {5,5nd,6,6nd}]
    [--len_def <int>] [--len_mip <csv:int>] [--len_min <csv:int>] [--dep_mip <csv:int>] [--dep_min <csv:int>] [--dep_sip <csv:int>] [--dep_sin <csv:int>]
    [--dp <int>] --err_out <str> [--nam_job <str>]


Description:
  Compute one per-sample siQ-ChIP or spike-in scaling-factor row from comma-separated BAM/CRAM lists and write it to a deterministic part file.


Keyword arguments:
  -h, --hlp, --help  <flag>
    Print this help message and exit.

  -en, --env, --env_nam  <str>
    Mamba environment to activate (default: ${env_nam}).

  -ds, --dir_scr  <str>
    Directory containing scripts and functions.

  -t, --thr, --threads  <int>
    Number of threads for alignment-processing steps (default: ${threads}).

  -md, --mode  <str>
    Scaling-factor framework: 'siq' or 'spike' (default: ${mode}).

  -me, --method  <str>
    Spike-in coefficient to compute: 'fractional', 'chiprx_alpha_ratio', 'chiprx_alpha_ip', 'chiprx_alpha_in', 'rxinput_alpha', or aliases ('--mode spike'; default: chiprx_alpha_ratio).

  -mp, --csv_mip  <csv:file>
    Comma-separated list of main IP BAM or CRAM files.

  -mn, --csv_min  <csv:file>
    Comma-separated list of main input BAM or CRAM files.

  -sp, --csv_sip  <csv:file>
    Comma-separated list of spike-in IP BAM or CRAM files ('--mode spike').

  -sn, --csv_sin  <csv:file>
    Comma-separated list of spike-in input BAM or CRAM files ('--mode spike').

  -at, --aln_typ, --align_typ  <str>
    Library type override: 'pe', 'se', or 'auto' (default: ${aln_typ}).

  -r, --ref, --ref_fa, --reference  <file>
    Reference FASTA required when any input alignment file is CRAM.

  -fo, --fil_out  <str>
    Base output TSV path used to derive '<fil_out>.part.<idx>'.

  -io, --idx_out  <int>
    Optional output-part index override for execute-layer local dispatch.

  -tb, --tbl_met  <str>
    siQ-ChIP metadata table ('--mode siq').

  -cm, --cfg_met  <str>
    YAML configuration for 'parse_metadata_siqchip.py' ('--mode siq').

  -eq, --eqn  <str>
    siQ-ChIP equation: '5', '5nd', '6', or '6nd' ('--mode siq'; default: ${eqn}).

  -ld, --len_def  <int>
    Default fragment length for SE libraries when no per-sample override is provided.

  -lmp, --len_mip  <csv:int>
    Optional comma-separated list of precomputed fragment lengths for main IP alignment files.

  -lmn, --len_min  <csv:int>
    Optional comma-separated list of precomputed fragment lengths for main input alignment files.

  -dmp, --dep_mip  <csv:int>
    Optional comma-separated list of precomputed alignment counts for main IP alignment files.

  -dmn, --dep_min  <csv:int>
    Optional comma-separated list of precomputed alignment counts for main input alignment files.

  -dsp, --dep_sip  <csv:int>
    Optional comma-separated list of precomputed alignment counts for spike IP alignment files.

  -dsn, --dep_sin  <csv:int>
    Optional comma-separated list of precomputed alignment counts for spike input alignment files.

  -dp, --dp, --rnd, --round, --decimals, --digits  <int>
    Maximum number of decimal places retained for computed values (default: ${rnd}).

  -eo, --err_out  <str>
    Directory for stderr/stdout log files.

  -nj, --nam_job  <str>
    Job-name prefix (default depends on resolved mode/method).


Dependencies:
  External programs:
    - Bash >= 4.4
    - bc
    - GNU AWK
    - Python
    - Samtools
    - Slurm

  Python scripts:
    - calculate_scaling_factor_siqchip.py
    - calculate_scaling_factor_spike.py
    - parse_metadata_siqchip.py

  Configuration file:
    - parse_metadata_siqchip.yml

  Sourced function scripts:
    - calculate_scaling_factor.sh
    - check_args.sh
    - check_inputs.sh
    - check_numbers.sh
    - format_outputs.sh
    - handle_env.sh
    - manage_slurm.sh


Notes:
  - Input alignment paths supplied to this wrapper should not contain spaces, commas, or semicolons.
  - This wrapper expects coordinate-sorted BAM or CRAM files.
  - CRAM inputs require '--ref_fa'.
  - Optional override vectors may be omitted, may contain a single broadcast value, or may contain one value per sample.
  - For '--mode siq' with SE data, '--len_def' or both '--len_mip' and '--len_min' must be supplied.
  - This worker writes one '<fil_out>.part.<idx>' row. The execute-layer wrapper combines successful part files into the final TSV.
  - '--idx_out' is an internal execute-layer coordination option. Omit it for direct multi-sample submit calls and Slurm array tasks.
  - Compute downstream denominator floors separately with 'python -m scripts.compute_input_floor'.
  - To run in "debug mode", set hardcoded variable 'debug=true'                   [debug=${debug:-UNSET}]
  - To run in "parse-only mode", set hardcoded variable 'p_only=true'             [p_only=${p_only:-UNSET}]
  - To run in "parse-and-check-only mode", set hardcoded variable 'pc_only=true'  [pc_only=${pc_only:-UNSET}]


Example:
  '''bash
  bash \${HOME}/repos/protocol_chipseq_signal_norm/scripts/submit_calculate_scaling_factor.sh
      -ds \${HOME}/repos/protocol_chipseq_signal_norm/scripts
      -en env_protocol
      -t 4
      -md spike
      -me chiprx_alpha_ratio
      -mp \${HOME}/repos/protocol_chipseq_signal_norm/data/processed/align_reads/bowtie2_global_flag-2_mapq-1/sc/IP_WT_Q_Hmo1_7751.sc.bam
      -mn \${HOME}/repos/protocol_chipseq_signal_norm/data/processed/align_reads/bowtie2_global_flag-2_mapq-1/sc/in_WT_Q_Hmo1_7751.sc.bam
      -sp \${HOME}/repos/protocol_chipseq_signal_norm/data/processed/align_reads/bowtie2_global_flag-2_mapq-1/sp/IP_WT_Q_Hmo1_7751.sp.bam
      -sn \${HOME}/repos/protocol_chipseq_signal_norm/data/processed/align_reads/bowtie2_global_flag-2_mapq-1/sp/in_WT_Q_Hmo1_7751.sp.bam
      -fo \${HOME}/repos/protocol_chipseq_signal_norm/data/processed/compute_signal/bowtie2_global_flag-2_mapq-1/tables/test_spike.tsv
      -dp 24
      -eo \${HOME}/repos/protocol_chipseq_signal_norm/data/processed/compute_signal/bowtie2_global_flag-2_mapq-1/tables/logs
      -nj calc_sf_spike
  '''


#TODO:
  - More examples.
EOM
}


#  Print parsed scalar arguments when debugging is enabled
function print_state_debug() {
    if [[ "${debug}" == "true" ]]; then
        debug_var \
            "env_nam=${env_nam}" \
            "dir_scr=${dir_scr}" \
            "threads=${threads}" \
            "mode=${mode}" \
            "method=${method:-UNSET}"

        if [[ "${mode}" == "spike" ]]; then
            debug_var "coef_spk=${coef_spk}"
        fi

        debug_var \
            "csv_mip=${csv_mip}" \
            "csv_min=${csv_min}"

        if [[ "${mode}" == "spike" ]]; then
            debug_var \
                "csv_sip=${csv_sip}" \
                "csv_sin=${csv_sin}"
        fi

        debug_var \
            "aln_typ=${aln_typ}" \
            "ref_fa=${ref_fa:-UNSET}" \
            "fil_out=${fil_out}" \
            "idx_out=${idx_out:-UNSET}"

        if [[ "${mode}" == "siq" ]]; then
            debug_var \
                "tbl_met=${tbl_met}" \
                "cfg_met=${cfg_met}" \
                "eqn=${eqn}"
        fi

        debug_var \
            "len_def=${len_def:-UNSET}" \
            "len_mip=${len_mip:-UNSET}" \
            "len_min=${len_min:-UNSET}" \
            "dep_mip=${dep_mip:-UNSET}" \
            "dep_min=${dep_min:-UNSET}" \
            "dep_sip=${dep_sip:-UNSET}" \
            "dep_sin=${dep_sin:-UNSET}" \
            "rnd=${rnd}" \
            "err_out=${err_out}" \
            "nam_job=${nam_job}"
    fi
}


#  Reconstruct required and optional arrays from serialized argument strings
function prepare_vecs() {
    IFS=',' read -r -a arr_mip <<< "${csv_mip}"
    IFS=',' read -r -a arr_min <<< "${csv_min}"

    if [[ "${mode}" == "spike" ]]; then
        IFS=',' read -r -a arr_sip <<< "${csv_sip}"
        IFS=',' read -r -a arr_sin <<< "${csv_sin}"
    fi

    if [[ -n "${len_mip}" ]]; then
        IFS=',' read -r -a arr_len_mip <<< "${len_mip}"
    else
        unset arr_len_mip && declare -ga arr_len_mip=()
    fi

    if [[ -n "${len_min}" ]]; then
        IFS=',' read -r -a arr_len_min <<< "${len_min}"
    else
        unset arr_len_min && declare -ga arr_len_min=()
    fi

    if [[ -n "${dep_mip}" ]]; then
        IFS=',' read -r -a arr_dep_mip <<< "${dep_mip}"
    else
        unset arr_dep_mip && declare -ga arr_dep_mip=()
    fi

    if [[ -n "${dep_min}" ]]; then
        IFS=',' read -r -a arr_dep_min <<< "${dep_min}"
    else
        unset arr_dep_min && declare -ga arr_dep_min=()
    fi

    if [[ -n "${dep_sip}" ]]; then
        IFS=',' read -r -a arr_dep_sip <<< "${dep_sip}"
    else
        unset arr_dep_sip && declare -ga arr_dep_sip=()
    fi

    if [[ -n "${dep_sin}" ]]; then
        IFS=',' read -r -a arr_dep_sin <<< "${dep_sin}"
    else
        unset arr_dep_sin && declare -ga arr_dep_sin=()
    fi
}


#  Validate vector lengths, sentinels, CRAM reference needs, and values
function validate_vecs() {
    local fil_aln need_ref=false n_samp
    local -a arr_aln

    check_arr_nonempty "arr_mip" "csv_mip" || return 1
    check_arr_nonempty "arr_min" "csv_min" || return 1

    if [[ "${mode}" == "spike" ]]; then
        check_arr_nonempty "arr_sip" "csv_sip" || return 1
        check_arr_nonempty "arr_sin" "csv_sin" || return 1
    fi

    arr_aln=( "${arr_mip[@]}" "${arr_min[@]}" )

    if [[ "${mode}" == "spike" ]]; then
        arr_aln+=( "${arr_sip[@]}" "${arr_sin[@]}" )
    fi

    for fil_aln in "${arr_aln[@]}"; do
        if [[ "${fil_aln,,}" == *.cram ]]; then
            need_ref=true
            break
        fi
    done

    if [[ "${need_ref}" == "true" && -z "${ref_fa}" ]]; then
        echo_err "'--ref_fa' is required when an input alignment file is CRAM."
        return 1
    fi

    if [[ -n "${ref_fa}" ]]; then
        validate_var_file "ref_fa" "${ref_fa}" || return 1
    fi

    if [[ "${mode}" == "siq" ]]; then
        check_arr_lengths "arr_mip" "arr_min" || return 1
        n_samp="${#arr_mip[@]}"
    else
        check_arr_lengths "arr_mip" "arr_min" || return 1
        check_arr_lengths "arr_mip" "arr_sip" || return 1
        check_arr_lengths "arr_mip" "arr_sin" || return 1
        n_samp="${#arr_mip[@]}"
    fi

    check_arr_len_bcst \
        "${n_samp}" \
        arr_len_mip arr_len_min \
        arr_dep_mip arr_dep_min \
        arr_dep_sip arr_dep_sin \
        || return 1

    if (( ${#arr_len_mip[@]} > 0 )); then
        check_arr_num_pos arr_len_mip len_mip || return 1
    fi

    if (( ${#arr_len_min[@]} > 0 )); then
        check_arr_num_pos arr_len_min len_min || return 1
    fi

    if (( ${#arr_dep_mip[@]} > 0 )); then
        check_arr_int_pos arr_dep_mip dep_mip || return 1
    fi

    if (( ${#arr_dep_min[@]} > 0 )); then
        check_arr_int_pos arr_dep_min dep_min || return 1
    fi

    if (( ${#arr_dep_sip[@]} > 0 )); then
        check_arr_int_pos arr_dep_sip dep_sip || return 1
    fi

    if (( ${#arr_dep_sin[@]} > 0 )); then
        check_arr_int_pos arr_dep_sin dep_sin || return 1
    fi
}


#  Print reconstructed array values when debugging is enabled
function print_vecs_debug() {
    if [[ "${debug}" == "true" ]]; then
        debug_var \
            "\${#arr_mip[@]}=${#arr_mip[@]}" \
            "\${#arr_min[@]}=${#arr_min[@]}"
        echo "arr_mip=( ${arr_mip[*]} )" >&2 && echo >&2
        echo "arr_min=( ${arr_min[*]} )" >&2 && echo >&2

        if [[ "${mode}" == "spike" ]]; then
            debug_var \
                "\${#arr_sip[@]}=${#arr_sip[@]}" \
                "\${#arr_sin[@]}=${#arr_sin[@]}"
            echo "arr_sip=( ${arr_sip[*]} )" >&2 && echo >&2
            echo "arr_sin=( ${arr_sin[*]} )" >&2 && echo >&2
        fi

        echo "arr_len_mip=( ${arr_len_mip[*]} )" >&2 && echo >&2
        echo "arr_len_min=( ${arr_len_min[*]} )" >&2 && echo >&2
        echo "arr_dep_mip=( ${arr_dep_mip[*]} )" >&2 && echo >&2
        echo "arr_dep_min=( ${arr_dep_min[*]} )" >&2 && echo >&2

        if [[ "${mode}" == "spike" ]]; then
            echo "arr_dep_sip=( ${arr_dep_sip[*]} )" >&2 && echo >&2
            echo "arr_dep_sin=( ${arr_dep_sin[*]} )" >&2 && echo >&2
        fi
    fi
}


#  Validate required input alignment files
function validate_input_files() {
    local idx

    for idx in "${!arr_mip[@]}"; do
        validate_var_file "arr_mip" "${arr_mip[${idx}]}" "${idx}" || return 1
        validate_var_file "arr_min" "${arr_min[${idx}]}" "${idx}" || return 1
    done

    if [[ "${mode}" == "spike" ]]; then
        for idx in "${!arr_sip[@]}"; do
            validate_var_file "arr_sip" "${arr_sip[${idx}]}" "${idx}" \
                || return 1
            validate_var_file "arr_sin" "${arr_sin[${idx}]}" "${idx}" \
                || return 1
        done
    fi
}


#  Activate the environment and ensure project imports are discoverable
function setup_env() {
    handle_env "${env_nam}" || return 1

    if [[ -z "${PYTHONPATH:-}" ]]; then
        export PYTHONPATH="${dir_rep}"
    else
        export PYTHONPATH="${dir_rep}:${PYTHONPATH}"
    fi
}


#  Dispatch Slurm-array or serial work
function run_jobs() {
    local err_dsc err_ini id_job id_tsk idx out_dsc out_ini samp

    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        id_job=${SLURM_ARRAY_JOB_ID}
        id_tsk=${SLURM_ARRAY_TASK_ID}

        if ! [[ "${id_tsk}" =~ ^[1-9][0-9]*$ ]]; then
            echo_err "Slurm task ID is invalid: '${id_tsk}'."
            return 1
        elif [[ "${id_tsk}" -gt "${#arr_mip[@]}" ]]; then
            echo_err \
                "Slurm task ID '${id_tsk}' exceeds number of samples:" \
                "'${#arr_mip[@]}'."
            return 1
        else
            idx=$(( id_tsk - 1 ))
        fi

        if [[ "${debug}" == "true" ]]; then
            debug_var "id_job=${id_job}" "id_tsk=${id_tsk}" "idx=${idx}"
        fi

        samp="$(derive_samp_sf "${arr_mip[idx]}")" || return 1

        if [[ "${debug}" == "true" ]]; then debug_var "samp=${samp}"; fi

        IFS=',' read -r err_ini out_ini err_dsc out_dsc < <(
            set_logs_slurm \
                "${id_job}" "${id_tsk}" "${samp}" "${err_out}" "${nam_job}"
        ) || return 1
        unset IFS

        if [[ "${debug}" == "true" ]]; then
            debug_var \
                "err_ini=${err_ini}" "out_ini=${out_ini}" \
                "err_dsc=${err_dsc}" "out_dsc=${out_dsc}"
        fi

        case "${mode}" in
            siq)   process_samp_siq   "${idx}" || return 1 ;;
            spike) process_samp_spike "${idx}" || return 1 ;;
        esac

        rm -f "${err_ini}" "${out_ini}" || {
            echo_warn \
                "failed to remove initial Slurm log file(s):" \
                "'${err_ini}' and/or '${out_ini}'."
        }
    else
        for idx in "${!arr_mip[@]}"; do
            case "${mode}" in
                siq)   process_samp_siq   "${idx}" || return 1 ;;
                spike) process_samp_spike "${idx}" || return 1 ;;
            esac
        done
    fi
}


#  Main script execution
function main() {
    init_defs

    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        show_help_main
        echo >&2
        return 0
    fi

    dir_scr="$(resolve_dir_scr "${0##*/}" "$@")" || return 1

    source_helpers_submit "${0##*/}" "${dir_scr}" \
        calculate_scaling_factor \
        check_args \
        check_inputs \
        check_numbers \
        format_outputs \
        handle_env \
        manage_slurm \
        || return 1

    parse_args "$@" || return 1

    if [[ "${p_only}" == "true" ]]; then
        if [[ "${debug}" == "true" ]]; then debug_var "p_only=true"; fi
        return 0
    fi

    canonicalize_args    || return 1
    validate_args        || return 1
    resolve_script_paths || return 1
    print_state_debug    || return 1
    prepare_vecs         || return 1
    validate_vecs        || return 1
    print_vecs_debug     || return 1
    validate_input_files || return 1

    if [[ "${pc_only}" == "true" ]]; then
        if [[ "${debug}" == "true" ]]; then debug_var "pc_only=true"; fi
        return 0
    fi

    setup_env || return 1
    run_jobs  || return 1
}


main "$@"
