#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: submit_calculate_scaling_factor.sh
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


# Require Bash >= 4.4 before doing any work.
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

# Run in safe mode, exiting on errors, unset variables, and pipe failures.
set -euo pipefail

# 'sbatch' copies this script, so '--dir_scr' must outrank 'BASH_SOURCE'.
function resolve_dir_scr() {
    local -a args=( "$@" )
    local i

    for (( i = 0; i < ${#args[@]}; i++ )); do
        case "${args[i]}" in
            -ds|--dir[_-]scr)
                if \
                       (( i + 1 < ${#args[@]} )) \
                    && [[ -n "${args[i + 1]}" ]] \
                    && [[ "${args[i + 1]}" != -* ]]
                then
                    printf '%s\n' "${args[i + 1]}"
                    return 0
                fi
                ;;
        esac
    done

    cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
}


dir_scr="$(resolve_dir_scr "$@")" || {
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "failed to resolve the script directory." >&2
    exit 1
}

fil_hlp="${dir_scr}/../lib/bash/help/help_submit_calculate_scaling_factor.sh"

# shellcheck source=lib/bash/help/help_submit_calculate_scaling_factor.sh
source "${fil_hlp}" || {
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "failed to source main help." >&2
    exit 1
}
unset fil_hlp


# Define script-specific functions.
function derive_samp_sf() {
    local fil_aln="${1:-}"
    local samp=""
    local show_help

    show_help=$(cat << EOM
Usage
-----
  derive_samp_sf
    [--help] fil_aln

  Derive a sample name from a BAM or CRAM path for logging.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  fil_aln : file
    Input BAM or CRAM alignment file from which to derive the sample name.

Returns
-------
  Prints the derived sample name to stdout; returns 1 if argument parsing fails.

Notes
-----
  Runtime requirements:
    - basename
    - bash >= 4.4

Examples
--------
  1. Strip an IP prefix and BAM extension from a main alignment path.
    '''bash
    derive_samp_sf results/IP_WT_G1.bam
    '''

  2. Capture a sample derived from an input-prefixed CRAM path.
    '''bash
    if samp="\$(derive_samp_sf results/in_WT_G1.cram)"; then
        printf '%s\n' "\${samp}"
    fi
    '''
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


# Source 'source_helpers.sh' and requested helper scripts from 'dir_scr'.
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

    dir_fnc="${dir_scr_arg}/../lib/bash"
    fnc_src="${dir_fnc}/core/source_helpers.sh"

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


# Parse keyword arguments after helper scripts have been sourced.
function parse_args() {
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
            -en|--env|--env[_-]nam)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_calculate_scaling_factor
                    return 1
                }
                env_nam="${2}"
                shift 2
                ;;

            -ds|--dir[_-]scr)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_calculate_scaling_factor
                    return 1
                }
                dir_scr="${2}"
                shift 2
                ;;

            -t|--thr|--threads)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_calculate_scaling_factor
                    return 1
                }
                threads="${2}"
                shift 2
                ;;

            -md|--mode)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_calculate_scaling_factor
                    return 1
                }
                mode="$(printf '%s\n' "${2}" | tr '[:upper:]' '[:lower:]')"
                shift 2
                ;;

            -me|--method)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_calculate_scaling_factor
                    return 1
                }
                method="$(printf '%s\n' "${2}" | tr '[:upper:]' '[:lower:]')"
                shift 2
                ;;

            -at|--aln[_-]typ|--align[_-]typ)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_calculate_scaling_factor
                    return 1
                }
                aln_typ="$(printf '%s\n' "${2}" | tr '[:upper:]' '[:lower:]')"
                shift 2
                ;;

            -rf|--ref[_-]fa)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_calculate_scaling_factor
                    return 1
                }
                ref_fa="${2}"
                shift 2
                ;;

            -cmip|--csv[_-]mip)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_calculate_scaling_factor
                    return 1
                }
                csv_mip="${2}"
                shift 2
                ;;

            -cmin|--csv[_-]min)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_calculate_scaling_factor
                    return 1
                }
                csv_min="${2}"
                shift 2
                ;;

            -csip|--csv[_-]sip)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_calculate_scaling_factor
                    return 1
                }
                csv_sip="${2}"
                shift 2
                ;;

            -csin|--csv[_-]sin)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_calculate_scaling_factor
                    return 1
                }
                csv_sin="${2}"
                shift 2
                ;;

            -fo|--fil[_-]out)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_calculate_scaling_factor
                    return 1
                }
                fil_out="${2}"
                shift 2
                ;;

            -io|--idx[_-]out)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_calculate_scaling_factor
                    return 1
                }
                idx_out="${2}"
                shift 2
                ;;

            -tb|--tbl[_-]met)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_calculate_scaling_factor
                    return 1
                }
                tbl_met="${2}"
                shift 2
                ;;

            -cm|--cfg[_-]met)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_calculate_scaling_factor
                    return 1
                }
                cfg_met="${2}"
                shift 2
                ;;

            -eq|--eqn)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_calculate_scaling_factor
                    return 1
                }
                eqn="$(printf '%s\n' "${2}" | tr '[:upper:]' '[:lower:]')"
                shift 2
                ;;

            -ld|--len[_-]def)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_calculate_scaling_factor
                    return 1
                }
                len_def="${2}"
                shift 2
                ;;

            -clmp|--csv[_-]len[_-]mip)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_calculate_scaling_factor
                    return 1
                }
                csv_len_mip="${2}"
                shift 2
                ;;

            -clmn|--csv[_-]len[_-]min)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_calculate_scaling_factor
                    return 1
                }
                csv_len_min="${2}"
                shift 2
                ;;

            -cdmp|--csv[_-]dep[_-]mip)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_calculate_scaling_factor
                    return 1
                }
                csv_dep_mip="${2}"
                shift 2
                ;;

            -cdmn|--csv[_-]dep[_-]min)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_calculate_scaling_factor
                    return 1
                }
                csv_dep_min="${2}"
                shift 2
                ;;

            -cdsp|--csv[_-]dep[_-]sip)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_calculate_scaling_factor
                    return 1
                }
                csv_dep_sip="${2}"
                shift 2
                ;;

            -cdsn|--csv[_-]dep[_-]sin)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_calculate_scaling_factor
                    return 1
                }
                csv_dep_sin="${2}"
                shift 2
                ;;

            -dp|--dp)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_calculate_scaling_factor
                    return 1
                }
                dp="${2}"
                shift 2
                ;;

            -deo|--dir[_-]eo)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_calculate_scaling_factor
                    return 1
                }
                dir_eo="${2}"
                shift 2
                ;;

            -nj|--nam[_-]job)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_calculate_scaling_factor
                    return 1
                }
                nam_job="${2}"
                shift 2
                ;;

            *)
                echo_err "unknown option/parameter passed: '${1}'."
                echo >&2
                help_submit_calculate_scaling_factor
                return 1
                ;;
        esac
    done
}


# Canonicalize mode, method, alignment type, coefficient name, and job name.
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


# Validate required arguments and simple numeric values.
function validate_args() {
    validate_var     "env_nam" "${env_nam}"         || return 1
    validate_var_dir "dir_scr" "${dir_scr}" 0 false || return 1
    validate_var     "threads" "${threads}"         || return 1
    check_int_pos    "${threads}" "threads"         || return 1
    check_int_pos    "${dp}"     "dp"             || return 1

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
                && ( -z "${csv_len_mip}" || -z "${csv_len_min}" )
            ]]; then
                echo_err \
                    "for '--mode siq' with SE data, supply '--len_def' or" \
                    "both '--csv_len_mip' and '--csv_len_min'."
                return 1
            fi
        fi
    fi

    if [[ -n "${len_def}" ]]; then
        check_int_pos "${len_def}" "len_def" || return 1
    fi

    validate_var "nam_job" "${nam_job}" || return 1

    if [[ -z "${dir_eo}" ]]; then
        echo_err "'--dir_eo' is required."
        return 1
    fi

    validate_var_dir "dir_eo" "${dir_eo}" || return 1
}


# Resolve installed Python command names after helper sourcing.
function resolve_script_paths() {
    validate_var_dir "dir_scr" "${dir_scr}" 0 false || return 1

    # 'run_python.sh' reads the caller-owned repository root dynamically.
    # shellcheck disable=SC2034
    dir_rep="$(cd "${dir_scr}/.." && pwd)"

    scr_met=parse_metadata_siqchip
    scr_siq=calculate_scaling_factor_siqchip
    scr_spk=calculate_scaling_factor_spike

    if [[ "${PY_INVOKE}" == "console" ]]; then
        check_pgrm_path "${scr_met}" || return 1
        check_pgrm_path "${scr_siq}" || return 1
        check_pgrm_path "${scr_spk}" || return 1
    fi
}


# Initialize hardcoded argument variables.
function init_args_hardcoded() {
    # If true, print verbose/debug Bash-level logging.
    debug=true

    # If true, parse arguments and exit before validation or execution.
    p_only=false

    # If true, parse and check arguments, then exit before execution.
    pc_only=false
}


# Initialize argument variables, assigning default values where applicable.
function init_arg_defs() {
    env_nam="env_protocol"
    threads=1

    mode="spike"    # siq | spike
    method=""       # chiprx_alpha_ratio | fractional | rxinput_alpha | ...

    aln_typ="auto"  # paired | pe | single | se | auto
    ref_fa=""

    csv_mip=""
    csv_min=""
    csv_sip=""
    csv_sin=""
    fil_out=""
    idx_out=""

    tbl_met=""
    cfg_met=""
    eqn="6nd"

    len_def=""      # SE fallback length (bp)
    csv_len_mip=""      # Precomputed lengths: main IP
    csv_len_min=""      # Precomputed lengths: main input
    csv_dep_mip=""      # Precomputed counts: main IP
    csv_dep_min=""      # Precomputed counts: main input
    csv_dep_sip=""      # Precomputed counts: spike IP
    csv_dep_sin=""      # Precomputed counts: spike input

    dp=24
    dir_eo=""
    nam_job=""

    unset arr_mip arr_min arr_sip arr_sin
    unset arr_len_mip arr_len_min
    unset arr_dep_mip arr_dep_min arr_dep_sip arr_dep_sin

    declare -ga arr_mip arr_min arr_sip arr_sin
    declare -ga arr_len_mip arr_len_min
    declare -ga arr_dep_mip arr_dep_min arr_dep_sip arr_dep_sin
}


# Initialize hardcoded arguments and user-facing argument defaults.
function init_defs() {
    init_args_hardcoded
    init_arg_defs
}


# Print parsed scalar arguments when debugging is enabled.
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
            "aln_typ=${aln_typ}" \
            "ref_fa=${ref_fa:-UNSET}" \
            "csv_mip=${csv_mip}" \
            "csv_min=${csv_min}"

        if [[ "${mode}" == "spike" ]]; then
            debug_var \
                "csv_sip=${csv_sip}" \
                "csv_sin=${csv_sin}"
        fi

        debug_var \
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
            "csv_len_mip=${csv_len_mip:-UNSET}" \
            "csv_len_min=${csv_len_min:-UNSET}" \
            "csv_dep_mip=${csv_dep_mip:-UNSET}" \
            "csv_dep_min=${csv_dep_min:-UNSET}" \
            "csv_dep_sip=${csv_dep_sip:-UNSET}" \
            "csv_dep_sin=${csv_dep_sin:-UNSET}" \
            "dp=${dp}" \
            "dir_eo=${dir_eo}" \
            "nam_job=${nam_job}"
    fi
}


# Reconstruct required and optional arrays from serialized argument strings.
function prepare_vecs() {
    IFS=',' read -r -a arr_mip <<< "${csv_mip}"
    IFS=',' read -r -a arr_min <<< "${csv_min}"

    if [[ "${mode}" == "spike" ]]; then
        IFS=',' read -r -a arr_sip <<< "${csv_sip}"
        IFS=',' read -r -a arr_sin <<< "${csv_sin}"
    fi

    if [[ -n "${csv_len_mip}" ]]; then
        IFS=',' read -r -a arr_len_mip <<< "${csv_len_mip}"
    else
        unset arr_len_mip && declare -ga arr_len_mip=()
    fi

    if [[ -n "${csv_len_min}" ]]; then
        IFS=',' read -r -a arr_len_min <<< "${csv_len_min}"
    else
        unset arr_len_min && declare -ga arr_len_min=()
    fi

    if [[ -n "${csv_dep_mip}" ]]; then
        IFS=',' read -r -a arr_dep_mip <<< "${csv_dep_mip}"
    else
        unset arr_dep_mip && declare -ga arr_dep_mip=()
    fi

    if [[ -n "${csv_dep_min}" ]]; then
        IFS=',' read -r -a arr_dep_min <<< "${csv_dep_min}"
    else
        unset arr_dep_min && declare -ga arr_dep_min=()
    fi

    if [[ -n "${csv_dep_sip}" ]]; then
        IFS=',' read -r -a arr_dep_sip <<< "${csv_dep_sip}"
    else
        unset arr_dep_sip && declare -ga arr_dep_sip=()
    fi

    if [[ -n "${csv_dep_sin}" ]]; then
        IFS=',' read -r -a arr_dep_sin <<< "${csv_dep_sin}"
    else
        unset arr_dep_sin && declare -ga arr_dep_sin=()
    fi
}


# Validate vector lengths, sentinels, CRAM reference needs, and values.
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
        check_arr_num_pos arr_len_mip csv_len_mip || return 1
    fi

    if (( ${#arr_len_min[@]} > 0 )); then
        check_arr_num_pos arr_len_min csv_len_min || return 1
    fi

    if (( ${#arr_dep_mip[@]} > 0 )); then
        check_arr_int_pos arr_dep_mip csv_dep_mip || return 1
    fi

    if (( ${#arr_dep_min[@]} > 0 )); then
        check_arr_int_pos arr_dep_min csv_dep_min || return 1
    fi

    if (( ${#arr_dep_sip[@]} > 0 )); then
        check_arr_int_pos arr_dep_sip csv_dep_sip || return 1
    fi

    if (( ${#arr_dep_sin[@]} > 0 )); then
        check_arr_int_pos arr_dep_sin csv_dep_sin || return 1
    fi
}


# Print reconstructed array values when debugging is enabled.
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


# Validate required input alignment files.
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


# Activate the environment containing the managed editable installation.
function setup_env() {
    handle_env "${env_nam}" || return 1
}


# Dispatch Slurm-array or serial work.
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
                "${id_job}" "${id_tsk}" "${samp}" "${dir_eo}" "${nam_job}"
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


# Main script execution.
function main() {
    init_defs

    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        help_submit_calculate_scaling_factor
        echo >&2
        return 0
    fi

    source_helpers_submit "${0##*/}" "${dir_scr}" \
        calculate_scaling_factor \
        check_args \
        check_env \
        check_inputs \
        check_numbers \
        format_outputs \
        handle_env \
        manage_slurm \
        help/help_submit_calculate_scaling_factor \
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
