#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: submit_trim_fastqs.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
# GPT-5.6) were used in development and documentation.
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

#  Source main help before bootstrap argument parsing
_dir_scr_help="$(
    cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
)"

# shellcheck source=lib/bash/help/help_submit_trim_fastqs.sh
source "${_dir_scr_help}/../lib/bash/help/help_submit_trim_fastqs.sh" || {
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "failed to source main help." >&2
    exit 1
}
unset _dir_scr_help


#  Resolve '--dir_scr' before sourced parser helpers are available
function resolve_dir_scr() {
    local script="${1:-}"
    shift

    local -a args=( "$@" )
    local i=0

    if [[ -z "${script}" ]]; then
        script="unknown_script"
    fi

    for (( i = 0; i < ${#args[@]}; i++ )); do
        case "${args[i]}" in
            -ds|--dir[_-]scr)
                if (( i + 1 >= ${#args[@]} )) \
                    || [[ -z "${args[i + 1]:-}" || "${args[i + 1]}" == -* ]]
                then
                    echo "error(${script}):" \
                        "option '${args[i]}' requires a value." >&2
                    echo >&2
                    help_submit_trim_fastqs
                    return 1
                fi

                printf "%s\n" "${args[i + 1]}"
                return 0
                ;;
        esac
    done

    echo "error(${script}):" \
        "required option '--dir_scr' was not supplied." >&2
    echo >&2
    help_submit_trim_fastqs
    return 1
}


#  Source 'source_helpers.sh' and requested helper scripts from 'dir_scr'
function source_helpers_submit() {
    local script="${1:-}"
    local dir_scr_arg="${2:-}"
    local fnc_src

    if (( $# < 2 )); then
        echo "error(${script:-unknown_script}):" \
            "expected at least two arguments: 'script' and 'dir_scr_arg'." >&2
        return 1
    fi

    shift 2

    if [[ -z "${script}" ]]; then
        script="unknown_script"
    fi

    if [[ -z "${dir_scr_arg}" ]]; then
        echo "error(${script}):" \
            "positional argument 2, 'dir_scr_arg', is missing." >&2
        return 1
    elif [[ ! -d "${dir_scr_arg}" ]]; then
        echo "error(${script}):" \
            "script directory not found: '${dir_scr_arg}'." >&2
        return 1
    elif (( $# < 1 )); then
        echo "error(${script}):" \
            "at least one helper script name must be supplied." >&2
        return 1
    fi

    dir_fnc="${dir_scr_arg}/../lib/bash"
    fnc_src="${dir_fnc}/core/source_helpers.sh"

    if [[ ! -f "${fnc_src}" ]]; then
        echo "error(${script}):" \
            "script not found: '${fnc_src}'." >&2
        return 1
    fi

    # shellcheck disable=SC1090
    source "${fnc_src}" || {
        echo "error(${script}):" \
            "failed to source '${fnc_src}'." >&2
        return 1
    }

    source_helpers "${dir_fnc}" "$@" || {
        echo "error(${script}):" \
            "failed to source required helper scripts." >&2
        return 1
    }
}


#  Initialize hardcoded argument variables
function init_args_hardcoded() {
    debug=true
}


#  Initialize argument variables, assigning default values where applicable
function init_arg_defs() {
    env_nam="env_protocol"
    dir_scr=""
    threads=4
    csv_fil_in=""
    dir_out=""
    sfx_se=""
    sfx_pe=""
    dir_eo=""
    nam_job="trim_fastqs"

    unset arr_fil_in
    declare -ga arr_fil_in
}


#  Initialize all script defaults
function init_defs() {
    init_args_hardcoded
    init_arg_defs
}


#  Parse keyword arguments after helper scripts have been sourced
function parse_args() {
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
            -en|--env|--env[_-]nam)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_trim_fastqs
                    return 1
                }
                env_nam="${2}"
                shift 2
                ;;

            -ds|--dir[_-]scr)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_trim_fastqs
                    return 1
                }
                dir_scr="${2}"
                shift 2
                ;;

            -t|--thr|--threads)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_trim_fastqs
                    return 1
                }
                threads="${2}"
                shift 2
                ;;

            -ci|--csv[_-]fil_in)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_trim_fastqs
                    return 1
                }
                csv_fil_in="${2}"
                shift 2
                ;;

            -do|--dir[_-]out)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_trim_fastqs
                    return 1
                }
                dir_out="${2}"
                shift 2
                ;;

            -sxs|--sfx[_-]se|--suffix[_-]se)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_trim_fastqs
                    return 1
                }
                sfx_se="${2}"
                shift 2
                ;;

            -sxp|--sfx[_-]pe|--suffix[_-]pe)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_trim_fastqs
                    return 1
                }
                sfx_pe="${2}"
                shift 2
                ;;

            -deo|--dir[_-]eo)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_trim_fastqs
                    return 1
                }
                dir_eo="${2}"
                shift 2
                ;;

            -nj|--nam[_-]job)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_trim_fastqs
                    return 1
                }
                nam_job="${2}"
                shift 2
                ;;

            *)
                echo_err "unknown option/parameter passed: '${1}'."
                echo >&2
                help_submit_trim_fastqs
                return 1
                ;;
        esac
    done
}


#  Validate required arguments and paths
function validate_args() {
    validate_var     "env_nam"    "${env_nam}"         || return 1
    validate_var_dir "dir_scr"    "${dir_scr}" 0 false || return 1
    validate_var     "threads"    "${threads}"         || return 1
    validate_var     "csv_fil_in" "${csv_fil_in}"      || return 1
    validate_var_dir "dir_out"    "${dir_out}"         || return 1
    validate_var     "sfx_se"     "${sfx_se}"          || return 1
    validate_var     "sfx_pe"     "${sfx_pe}"          || return 1
    validate_var_dir "dir_eo"    "${dir_eo}"         || return 1
    validate_var     "nam_job"    "${nam_job}"         || return 1

    check_int_pos "${threads}" "threads" || return 1
}


#  Print debug argument variable assignments
function print_debug_args() {
    if [[ "${debug}" == "true" ]]; then
        echo
        debug_var \
            "env_nam=${env_nam}" \
            "dir_scr=${dir_scr}" \
            "threads=${threads}" \
            "csv_fil_in=${csv_fil_in}" \
            "dir_out=${dir_out}" \
            "sfx_se=${sfx_se}" \
            "sfx_pe=${sfx_pe}" \
            "dir_eo=${dir_eo}" \
            "nam_job=${nam_job}"
    fi
}


#  Reconstruct the input FASTQ vector from serialized input
function prepare_vecs() {
    IFS=';' read -r -a arr_fil_in <<< "${csv_fil_in}"
    check_arr_nonempty "arr_fil_in" "csv_fil_in" || return 1
}


#  Activate the requested environment
function setup_env() {
    handle_env "${env_nam}" || return 1
}


#  Check tools required by the submit worker
function check_tools() {
    check_pgrm_path atria || return 1
}


#  Print parsed vector state in debug mode
function print_vecs_debug() {
    if [[ "${debug}" != "true" ]]; then
        return 0
    fi

    echo "\${#arr_fil_in[@]}=${#arr_fil_in[@]}" && echo
    echo "arr_fil_in=( ${arr_fil_in[*]} )"      && echo
}


#  Parse one FASTQ entry into 'fq_1', 'fq_2', and 'samp'
function parse_entry_trim_fastq() {
    local fil_in="${1:-}"  # Input FASTQ file(s)
    local sfx_se="${2:-}"  # Suffix for SE FASTQ files
    local sfx_pe="${3:-}"  # Suffix for PE FASTQ files (FASTQ #1)
    local fq_1             # FASTQ file #1
    local fq_2             # FASTQ file #2, or 'NA' for SE
    local samp             # Sample name
    local show_help        # Help message

    show_help=$(cat << EOM
Usage
-----
  parse_entry_trim_fastq
    [--help] fil_in sfx_se sfx_pe

  Parse one Input file path. FASTQ input entry into 'fq_1', 'fq_2', and 'samp'.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  fil_in : file
    Input file path. FASTQ input entry. For single-end data, this is one FASTQ file. For paired-end data, this is a comma-delimited FASTQ pair.

  2  sfx_se : str
    Suffix to strip from single-end FASTQ filenames.

  3  sfx_pe : str
    Suffix to strip from paired-end FASTQ read-1 filenames.

Returns
-------
  Prints a semicolon-delimited record to stdout:

    fq_1;fq_2;samp

  where 'fq_2' is 'NA' for single-end data. Returns 0 when parsing succeeds; otherwise 1.

Notes
-----
  Runtime requirements:
    - awk
    - basename
    - bash >= 4.4

Examples
--------
  1. Parse one single-end trimmed FASTQ entry.
    '''bash
    parse_entry_trim_fastq \\
        reads/sample.atria.fastq.gz \\
        .atria.fastq.gz \\
        _R1.atria.fastq.gz
    '''

  2. Capture the parsed fields for one paired-end FASTQ entry.
    '''bash
    if \\
        parsed="\$(
            parse_entry_trim_fastq \\
                reads/sample_R1.atria.fastq.gz,reads/sample_R2.atria.fastq.gz \\
                .atria.fastq.gz \\
                _R1.atria.fastq.gz
        )" \\
    then
        printf '%s\n' "\${parsed}"
    fi
    '''
EOM
    )

    if [[ "${fil_in}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${fil_in}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'fil_in', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    validate_var "fil_in" "${fil_in}" || return 1
    validate_var "sfx_se" "${sfx_se}" || return 1
    validate_var "sfx_pe" "${sfx_pe}" || return 1

    IFS=';' read -r fq_1 fq_2 samp < <(
        parse_fastq_entry "${fil_in}" "${sfx_se}" "${sfx_pe}"
    ) || return 1

    validate_var "fq_1" "${fq_1}" || return 1

    if [[ "${fq_2}" != "NA" ]]; then
        validate_var "fq_2" "${fq_2}" || return 1
    fi

    echo "${fq_1};${fq_2};${samp}"
}


#  Parse one input entry, print debug state, and run Atria
function run_one_entry() {
    local fil_in="${1:-}"
    local fq_1=""
    local fq_2=""
    local samp=""
    local log_out=""
    local log_err=""

    if [[ "${debug}" == "true" ]]; then debug_var "fil_in=${fil_in}"; fi

    IFS=';' read -r fq_1 fq_2 samp < <(
        parse_entry_trim_fastq "${fil_in}" "${sfx_se}" "${sfx_pe}"
    ) || return 1

    log_out="${dir_eo}/${nam_job}.${samp}.stdout.txt"
    log_err="${dir_eo}/${nam_job}.${samp}.stderr.txt"

    if [[ "${debug}" == "true" ]]; then
        debug_var \
            "fq_1=${fq_1}" \
            "fq_2=${fq_2}" \
            "samp=${samp}" \
            "log_out=${log_out}" \
            "log_err=${log_err}"
    fi

    if ! \
        trim_fastqs_atria \
            "${threads}" "${fq_1}" "${fq_2}" "${dir_out}" "${log_out}" \
            "${log_err}" "${samp}"
    then
        echo_err "failed to perform read trimming."
        return 1
    fi
}


#  Run the Slurm-array task selected by 'SLURM_ARRAY_TASK_ID'
function run_slurm_task() {
    local id_job="${SLURM_ARRAY_JOB_ID:-}"
    local id_tsk="${SLURM_ARRAY_TASK_ID:-}"
    local idx=""
    local fil_in=""
    local fq_1=""
    local fq_2=""
    local samp=""
    local err_ini=""
    local out_ini=""
    local err_dsc=""
    local out_dsc=""

    if [[ -z "${id_job}" ]]; then
        echo_err "Slurm array job ID is missing."
        return 1
    fi

    if ! [[ "${id_tsk}" =~ ^[1-9][0-9]*$ ]]; then
        echo_err "Slurm task ID is invalid: '${id_tsk}'."
        return 1
    elif (( id_tsk > ${#arr_fil_in[@]} )); then
        echo_err \
            "Slurm task ID '${id_tsk}' exceeds number of FASTQ entries:" \
            "'${#arr_fil_in[@]}'."
        return 1
    else
        idx=$(( id_tsk - 1 ))
    fi

    fil_in="${arr_fil_in[idx]}"

    if [[ "${debug}" == "true" ]]; then debug_var "fil_in=${fil_in}"; fi

    IFS=';' read -r fq_1 fq_2 samp < <(
        parse_entry_trim_fastq "${fil_in}" "${sfx_se}" "${sfx_pe}"
    ) || return 1

    IFS=',' read -r err_ini out_ini err_dsc out_dsc < <(
        set_logs_slurm \
            "${id_job}" "${id_tsk}" "${samp}" "${dir_eo}" "${nam_job}"
    ) || return 1

    if [[ "${debug}" == "true" ]]; then
        debug_var \
            "fq_1=${fq_1}" \
            "fq_2=${fq_2}" \
            "samp=${samp}" \
            "err_ini=${err_ini}" \
            "out_ini=${out_ini}" \
            "err_dsc=${err_dsc}" \
            "out_dsc=${out_dsc}"
    fi

    if ! \
        trim_fastqs_atria \
            "${threads}" "${fq_1}" "${fq_2}" "${dir_out}" "${out_dsc}" \
            "${err_dsc}" "${samp}"
    then
        echo_err "failed to perform read trimming."
        return 1
    fi

    rm "${err_ini}" "${out_ini}"
}


#  Run all input entries locally in serial order
function run_local_jobs() {
    local idx=""
    local fil_in=""

    for idx in "${!arr_fil_in[@]}"; do
        fil_in="${arr_fil_in[idx]}"
        run_one_entry "${fil_in}" || return 1
    done
}


#  Dispatch one Slurm task or all local jobs
function run_jobs() {
    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        run_slurm_task || return 1
    else
        run_local_jobs || return 1
    fi
}


#  Main script execution
function main() {
    init_defs

    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        help_submit_trim_fastqs
        echo >&2
        return 0
    fi

    dir_scr="$(resolve_dir_scr "${0##*/}" "$@")" || return 1

    source_helpers_submit "${0##*/}" "${dir_scr}" \
        check_args \
        check_env \
        check_inputs \
        check_numbers \
        format_outputs \
        handle_env \
        manage_slurm \
        process_sequences \
        help/help_submit_trim_fastqs \
        || return 1

    parse_args "$@"    || return 1
    validate_args      || return 1
    prepare_vecs       || return 1
    print_debug_args   || return 1
    print_vecs_debug   || return 1
    setup_env          || return 1
    check_tools        || return 1
    run_jobs           || return 1
}


main "$@"
