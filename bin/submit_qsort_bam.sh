#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: submit_qsort_bam.sh
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

# 'sbatch' runs a copy, so 'BASH_SOURCE' can resolve outside the repo; if so,
# fail here rather than at a subsequent 'source'.
if [[ ! -d "${dir_scr}/../lib/bash" ]]; then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "cannot locate helper libraries from '${dir_scr}'; pass '--dir_scr'" \
        "when this script is run from a copy." >&2
    exit 1
fi

# shellcheck source=lib/bash/help/help_submit_qsort_bam.sh
source "${dir_scr}/../lib/bash/help/help_submit_qsort_bam.sh" || {
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "failed to source main help." >&2
    exit 1
}


# Initialize hardcoded argument variables.
function init_args_hardcoded() {
    # If true, run script in debug mode.
    debug=true
}


# Initialize argument variables, assigning default values where applicable.
function init_arg_defs() {
    env_nam="env_protocol"
    threads=1
    csv_fil_in=""
    dir_out=""
    dir_eo=""
    nam_job="qsort_bam"
    ref_fa=""

    dir_fnc=""

    unset arr_fil_ins
    declare -ga arr_fil_ins
}


# Initialize hardcoded arguments and user-facing argument defaults.
function init_defs() {
    init_args_hardcoded
    init_arg_defs
}


# Source 'source_helpers.sh' and requested helper scripts from 'dir_scr'.
function source_helpers_submit() {
    local scr="${1:-}"
    local dir_scr_arg="${2:-}"
    local fnc_src

    if (( $# < 2 )); then
        echo "error(${scr:-unknown_script}):" \
            "expected at least two arguments: 'scr' and 'dir_scr_arg'." >&2
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
        echo "error(${scr}): script not found: '${fnc_src}'." >&2
        return 1
    fi

    # shellcheck disable=SC1090
    source "${fnc_src}" || {
        echo "error(${scr}): failed to source '${fnc_src}'." >&2
        return 1
    }

    source_helpers "${dir_fnc}" "$@" || {
        echo "error(${scr}): failed to source required helper scripts." >&2
        return 1
    }
}


# Parse keyword arguments after helper scripts have been sourced.
function parse_args() {
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
            -h|--hlp|--help)
                help_submit_qsort_bam
                echo >&2
                return 0
                ;;

            -en|--env|--env[_-]nam)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_qsort_bam
                    return 1
                }
                env_nam="${2}"
                shift 2
                ;;

            -ds|--dir[_-]scr)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_qsort_bam
                    return 1
                }
                dir_scr="${2}"
                shift 2
                ;;

            -t|--thr|--threads)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_qsort_bam
                    return 1
                }
                threads="${2}"
                shift 2
                ;;

            -ci|--csv[_-]fil_in)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_qsort_bam
                    return 1
                }
                csv_fil_in="${2}"
                shift 2
                ;;

            -rf|--ref[_-]fa)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_qsort_bam
                    return 1
                }
                ref_fa="${2}"
                shift 2
                ;;

            -do|--dir[_-]out)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_qsort_bam
                    return 1
                }
                dir_out="${2}"
                shift 2
                ;;

            -deo|--dir[_-]eo)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_qsort_bam
                    return 1
                }
                dir_eo="${2}"
                shift 2
                ;;

            -nj|--nam[_-]job)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_qsort_bam
                    return 1
                }
                nam_job="${2}"
                shift 2
                ;;

            *)
                echo_err "unknown option/parameter passed: '${1}'."
                echo >&2
                help_submit_qsort_bam
                return 1
                ;;
        esac
    done
}


# Validate required arguments and scalar values.
function validate_args() {
    validate_var     "env_nam"    "${env_nam}"         || return 1
    validate_var_dir "dir_scr"    "${dir_scr}" 0 false || return 1
    validate_var     "threads"    "${threads}"         || return 1
    validate_var     "csv_fil_in" "${csv_fil_in}"      || return 1
    validate_var_dir "dir_out"    "${dir_out}"         || return 1
    validate_var_dir "dir_eo"     "${dir_eo}"          || return 1
    validate_var     "nam_job"    "${nam_job}"         || return 1

    check_int_pos "${threads}" "threads" || return 1
}


# Print parsed scalar arguments when debugging is enabled.
function print_state_debug() {
    if [[ "${debug}" != "true" ]]; then
        return 0
    fi

    debug_var \
        "env_nam=${env_nam}" \
        "dir_scr=${dir_scr}" \
        "threads=${threads}" \
        "csv_fil_in=${csv_fil_in}" \
        "dir_out=${dir_out}" \
        "dir_eo=${dir_eo}" \
        "nam_job=${nam_job}" \
        "ref_fa=${ref_fa:-UNSET}"
}


# Activate the requested Conda environment.
function setup_env() {
    handle_env "${env_nam}" || return 1
}


# Check runtime tool dependencies.
function check_tools() {
    check_pgrm_path samtools
}


# Reconstruct input alignment vector.
function prepare_vecs() {
    unset arr_fil_ins && declare -ga arr_fil_ins
    IFS=',' read -r -a arr_fil_ins <<< "${csv_fil_in}"
}


# Validate reconstructed input alignment vector.
function validate_vecs() {
    local idx fil_in need_ref=false

    check_arr_nonempty "arr_fil_ins" "csv_fil_in" || return 1

    for idx in "${!arr_fil_ins[@]}"; do
        fil_in="${arr_fil_ins[idx]}"
        validate_var_file "arr_fil_ins" "${fil_in}" "${idx}" || return 1

        case "${fil_in,,}" in
            *.bam) : ;;
            *.cram) need_ref=true ;;
            *)
                echo_err \
                    "input file must end in '.bam' or '.cram': '${fil_in}'."
                return 1
                ;;
        esac
    done

    if [[ "${need_ref}" == "true" && -z "${ref_fa}" ]]; then
        echo_err "'--ref_fa' is required when '--csv_fil_in' contains CRAM."
        return 1
    fi

    if [[ -n "${ref_fa}" ]]; then
        validate_var_file "ref_fa" "${ref_fa}" || return 1
    fi
}


# Print reconstructed array values when debugging is enabled.
function print_vecs_debug() {
    if [[ "${debug}" != "true" ]]; then
        return 0
    fi

    debug_var "\${#arr_fil_ins[@]}=${#arr_fil_ins[@]}"
    echo "arr_fil_ins=( ${arr_fil_ins[*]} )" >&2 && echo >&2
}


# Derive the qname-sorted output path for one input alignment file.
function derive_fil_out() {
    local fil_in="${1:-}"
    local base ext_out stem

    base="$(basename "${fil_in}")"

    case "${base,,}" in
        *.bam)
            stem="${base%.bam}"
            ext_out="bam"
            ;;

        *.cram)
            stem="${base%.cram}"
            ext_out="cram"
            ;;

        *)
            echo_err_func "${FUNCNAME[0]}" \
                "input file must end in '.bam' or '.cram': '${fil_in}'."
            return 1
            ;;
    esac

    printf '%s\n' "${dir_out}/${stem}.qnam.${ext_out}"
}


# Sort one alignment file by query name.
function run_job() {
    local idx="${1:-}"
    local fil_in fil_out samp log_out log_err

    if ! [[ "${idx}" =~ ^[0-9]+$ ]]; then
        echo_err_func "${FUNCNAME[0]}" "invalid task index: '${idx}'."
        return 1
    elif (( idx >= ${#arr_fil_ins[@]} )); then
        echo_err_func "${FUNCNAME[0]}" "task index is out of bounds: '${idx}'."
        return 1
    fi

    fil_in="${arr_fil_ins[idx]}"
    fil_out="$(derive_fil_out "${fil_in}")" || return 1
    samp="$(basename "${fil_out}")"
    log_out="${dir_eo}/${nam_job}.${samp}.stdout.txt"
    log_err="${dir_eo}/${nam_job}.${samp}.stderr.txt"

    if [[ "${debug}" == "true" ]]; then
        debug_var \
            "idx=${idx}" \
            "fil_in=${fil_in}" \
            "fil_out=${fil_out}" \
            "samp=${samp}"
    fi

    if ! \
        qsort_file_alignments \
            "${threads}" \
            "${fil_in}" \
            "${fil_out}" \
            "${ref_fa}" \
            "${log_out}" \
            "${log_err}"
    then
        echo_err "failed to sort file: '${fil_in}'. See log: '${log_err}'."
        return 1
    fi

    echo "Successfully sorted '${fil_in}' to '${fil_out}'." >&2
}


# Parse Slurm task state and run one array task.
function run_job_slurm() {
    local err_dsc err_ini id_job id_tsk idx out_dsc out_ini
    local fil_in fil_out samp

    id_job="${SLURM_ARRAY_JOB_ID:-}"
    id_tsk="${SLURM_ARRAY_TASK_ID:-}"

    if [[ -z "${id_job}" ]]; then
        echo_err "Slurm array job ID is missing."
        return 1
    elif ! [[ "${id_tsk}" =~ ^[1-9][0-9]*$ ]]; then
        echo_err "Slurm task ID is invalid: '${id_tsk}'."
        return 1
    elif (( id_tsk > ${#arr_fil_ins[@]} )); then
        echo_err \
            "Slurm task ID '${id_tsk}' exceeds number of input files:" \
            "'${#arr_fil_ins[@]}'."
        return 1
    fi

    idx=$(( id_tsk - 1 ))
    fil_in="${arr_fil_ins[idx]}"
    fil_out="$(derive_fil_out "${fil_in}")" || return 1
    samp="$(basename "${fil_out}")"

    # shellcheck disable=SC2034
    IFS=',' read -r err_ini out_ini err_dsc out_dsc < <(
        set_logs_slurm \
            "${id_job}" "${id_tsk}" "${samp}" "${dir_eo}" "${nam_job}"
    ) || return 1

    run_job "${idx}" || return 1

    rm -f "${err_ini}" "${out_ini}" || {
        echo_warn \
            "failed to remove initial Slurm log file(s):" \
            "'${err_ini}' and/or '${out_ini}'."
    }
}


# Dispatch one Slurm-array task or local loop over input files.
function run_jobs() {
    local idx

    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        run_job_slurm || return 1
    else
        for idx in "${!arr_fil_ins[@]}"; do
            run_job "${idx}" || return 1
        done
    fi
}


# Main script execution.
function main() {
    init_defs

    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        help_submit_qsort_bam
        echo >&2
        return 0
    fi

    source_helpers_submit "${0##*/}" "${dir_scr}" \
        check_args \
        check_env \
        check_inputs \
        check_numbers \
        format_outputs \
        handle_env \
        manage_slurm \
        process_alignments \
        help/help_submit_qsort_bam \
        || return 1

    parse_args "$@"   || return 1
    validate_args     || return 1
    print_state_debug || return 1
    setup_env         || return 1
    check_tools       || return 1
    prepare_vecs      || return 1
    validate_vecs     || return 1
    print_vecs_debug  || return 1
    run_jobs          || return 1
}


main "$@"
