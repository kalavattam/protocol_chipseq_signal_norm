#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: submit_trim_fastqs.sh
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
                    show_help_main
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
    show_help_main
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

    dir_fnc="${dir_scr_arg}/functions"
    fnc_src="${dir_fnc}/source_helpers.sh"

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
    csv_infile=""
    dir_out=""
    sfx_se=""
    sfx_pe=""
    err_out=""
    nam_job="trim_fastqs"

    unset arr_infile
    declare -ga arr_infile
}


#  Initialize all script defaults
function init_defs() {
    init_args_hardcoded
    init_arg_defs
}


function show_help_main() {
    cat << EOM >&2
Usage:
  submit_trim_fastqs.sh
    [-h|--hlp|--help]
    [-en|--env_nam <str>] -ds|--dir_scr <str>
    [-t|--threads <int>]
    -i|--csv_infile <str> -do|--dir_out <str>
    -sxs|--sfx_se <str> -sxp|--sfx_pe <str>
    -eo|--err_out <str> [-nj|--nam_job <str>]

Description:
  Submit or execute one or more FASTQ-trimming jobs by calling downstream program 'atria'.

  This wrapper
    - parses a semicolon-delimited list of FASTQ input entries,
    - derives sample names,
    - activates the requested Conda/Mamba environment, and then
    - runs read trimming either under Slurm array execution or by serial/GNU-Parallel-style iteration, depending on how the script is invoked.

  For each input entry, this script writes log files to:

    \${err_out}/\${nam_job}.\${samp}.stdout.txt
    \${err_out}/\${nam_job}.\${samp}.stderr.txt

Keyword arguments:
  -en, --env, --env_nam  <str>
    Conda/Mamba environment to activate.

  -ds, --dir_scr  <str>
    Directory containing scripts and functions.

  -t, --thr, --threads  <int>
    Number of threads to use.

  -i, --csv_infile  <str>
    Semicolon-delimited serialized string of FASTQ input entries.

    For single-end data, each entry is one FASTQ file. For paired-end data, each entry contains a comma-delimited FASTQ pair.

  -do, --dir_out  <str>
    Directory for trimmed FASTQ output files.

  -sxs, --sfx_se, --suffix_se  <str>
    Suffix to strip from SE FASTQ files.

  -sxp, --sfx_pe, --suffix_pe  <str>
    Suffix to strip from PE FASTQ read-1 files.

  -eo, --err_out  <str>
    Directory to store stderr and stdout outfiles.

  -nj, --nam_job  <str>
    Name of job.

Notes:
  - All arguments are required with the following notes and exceptions:
    + '--env_nam' defaults to 'env_nam=${env_nam}' if not specified.
    + '--threads' defaults to 'threads=${threads}' if not specified.
    + '--nam_job' defaults to 'nam_job=${nam_job}' if not specified.
EOM
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

            -i|-fi|-ci|--infile|--infiles|--fil[_-]in|--csv[_-]infile|--csv[_-]infiles)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                csv_infile="${2}"
                shift 2
                ;;

            -do|--dir[_-]out)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                dir_out="${2}"
                shift 2
                ;;

            -sxs|--sfx[_-]se|--suffix[_-]se)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                sfx_se="${2}"
                shift 2
                ;;

            -sxp|--sfx[_-]pe|--suffix[_-]pe)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    show_help_main
                    return 1
                }
                sfx_pe="${2}"
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


#  Validate required arguments and paths
function validate_args() {
    validate_var     "env_nam"    "${env_nam}"         || return 1
    validate_var_dir "dir_scr"    "${dir_scr}" 0 false || return 1
    validate_var     "threads"    "${threads}"         || return 1
    validate_var     "csv_infile" "${csv_infile}"      || return 1
    validate_var_dir "dir_out"    "${dir_out}"         || return 1
    validate_var     "sfx_se"     "${sfx_se}"          || return 1
    validate_var     "sfx_pe"     "${sfx_pe}"          || return 1
    validate_var_dir "err_out"    "${err_out}"         || return 1
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
            "csv_infile=${csv_infile}" \
            "dir_out=${dir_out}" \
            "sfx_se=${sfx_se}" \
            "sfx_pe=${sfx_pe}" \
            "err_out=${err_out}" \
            "nam_job=${nam_job}"
    fi
}


#  Reconstruct the input FASTQ vector from serialized input
function prepare_vecs() {
    IFS=';' read -r -a arr_infile <<< "${csv_infile}"
    check_arr_nonempty "arr_infile" "csv_infile" || return 1
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

    echo "\${#arr_infile[@]}=${#arr_infile[@]}" && echo
    echo "arr_infile=( ${arr_infile[*]} )"      && echo
}


#  Parse one FASTQ entry into 'fq_1', 'fq_2', and 'samp'
function parse_entry_trim_fastq() {
    local infile="${1:-}"  # Input FASTQ file(s)
    local sfx_se="${2:-}"  # Suffix for SE FASTQ files
    local sfx_pe="${3:-}"  # Suffix for PE FASTQ files (FASTQ #1)
    local fq_1             # FASTQ file #1
    local fq_2             # FASTQ file #2, or 'NA' for SE
    local samp             # Sample name
    local show_help        # Help message

    show_help=$(cat << EOM
Usage:
  parse_entry_trim_fastq
    [-h|--hlp|--help] infile sfx_se sfx_pe

Description:
  Parse one FASTQ input entry into 'fq_1', 'fq_2', and 'samp'.

Positional arguments:
  1  infile  <str>
    FASTQ input entry. For single-end data, this is one FASTQ file. For paired-end data, this is a comma-delimited FASTQ pair.

  2  sfx_se  <str>
    Suffix to strip from single-end FASTQ filenames.

  3  sfx_pe  <str>
    Suffix to strip from paired-end FASTQ read-1 filenames.

Returns:
  Prints a semicolon-delimited record to stdout:

    fq_1;fq_2;samp

  where 'fq_2' is 'NA' for single-end data. Returns 0 when parsing succeeds; otherwise 1.
EOM
    )

    if [[ "${infile}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${infile}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'infile', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    validate_var "infile" "${infile}" || return 1
    validate_var "sfx_se" "${sfx_se}" || return 1
    validate_var "sfx_pe" "${sfx_pe}" || return 1

    IFS=';' read -r fq_1 fq_2 samp < <(
        parse_fastq_entry "${infile}" "${sfx_se}" "${sfx_pe}"
    ) || return 1

    validate_var "fq_1" "${fq_1}" || return 1

    if [[ "${fq_2}" != "NA" ]]; then
        validate_var "fq_2" "${fq_2}" || return 1
    fi

    echo "${fq_1};${fq_2};${samp}"
}


#  Parse one input entry, print debug state, and run Atria
function run_one_entry() {
    local infile="${1:-}"
    local fq_1=""
    local fq_2=""
    local samp=""
    local log_out=""
    local log_err=""

    if [[ "${debug}" == "true" ]]; then debug_var "infile=${infile}"; fi

    IFS=';' read -r fq_1 fq_2 samp < <(
        parse_entry_trim_fastq "${infile}" "${sfx_se}" "${sfx_pe}"
    ) || return 1

    log_out="${err_out}/${nam_job}.${samp}.stdout.txt"
    log_err="${err_out}/${nam_job}.${samp}.stderr.txt"

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
    local infile=""
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
    elif (( id_tsk > ${#arr_infile[@]} )); then
        echo_err \
            "Slurm task ID '${id_tsk}' exceeds number of FASTQ entries:" \
            "'${#arr_infile[@]}'."
        return 1
    else
        idx=$(( id_tsk - 1 ))
    fi

    infile="${arr_infile[idx]}"

    if [[ "${debug}" == "true" ]]; then debug_var "infile=${infile}"; fi

    IFS=';' read -r fq_1 fq_2 samp < <(
        parse_entry_trim_fastq "${infile}" "${sfx_se}" "${sfx_pe}"
    ) || return 1

    IFS=',' read -r err_ini out_ini err_dsc out_dsc < <(
        set_logs_slurm \
            "${id_job}" "${id_tsk}" "${samp}" "${err_out}" "${nam_job}"
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
    local infile=""

    for idx in "${!arr_infile[@]}"; do
        infile="${arr_infile[idx]}"
        run_one_entry "${infile}" || return 1
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
        show_help_main
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
