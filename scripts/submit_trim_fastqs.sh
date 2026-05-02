#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: submit_trim_fastqs.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
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

#  If true, run script in debug mode
debug=true


#  Define functions
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

  The input entry may represent either:
    - single-end data: one FASTQ path
    - paired-end data: two comma-delimited FASTQ paths

  This helper calls 'parse_fastq_entry' to parse the FASTQ input and derive
  the sample name.

Positional arguments:
  1  infile  <str>  FASTQ input entry. For SE data, this is one FASTQ file. For PE data, this is a comma-delimited FASTQ pair.
  2  sfx_se  <str>  Suffix to strip from SE FASTQ filenames when deriving the sample name.
  3  sfx_pe  <str>  Suffix to strip from PE FASTQ read-1 filenames when deriving the sample name.

Returns:
  Prints a semicolon-delimited record to stdout:

    fq_1;fq_2;samp

  where 'fq_2' is 'NA' for SE data.

Notes:
  - This helper validates required inputs with 'validate_var'.
  - Parsing of the FASTQ entry itself is delegated to 'parse_fastq_entry'.
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

    #MAYBE: these 'fq_(1|2)' checks are a little redundant
    validate_var "fq_1" "${fq_1}" || return 1

    if [[ "${fq_2}" != "NA" ]]; then
        validate_var "fq_2" "${fq_2}" || return 1
    fi

    echo "${fq_1};${fq_2};${samp}"
}
#MAYBE: 'parse_entry_trim_fastq' still uses 'validate_var' rather than file-
#       aware helpers such as 'validate_file' and 'validate_var_file'; revisit
#       if stronger file validation is later wanted here


#  Run Atria with explicitly passed argument values
function run_atria() {
    local threads="${1:-}"  # Number of threads to use
    local fq_1="${2:-}"     # First FASTQ file (required)
    local fq_2="${3:-}"     # Second FASTQ file (optional, for PE reads)
    local dir_out="${4:-}"  # Output directory for processed FASTQ files
    local err_out="${5:-}"  # Directory for logging stderr/stdout
    local nam_job="${6:-}"  # Job name for log file naming
    local samp="${7:-}"     # Sample name for log file naming
    local log_out log_err   # Atria stdout and stderr log files
    local show_help

    show_help=$(cat << EOM
Usage:
  run_atria
    [-h|--hlp|--help] threads fq_1 fq_2 dir_out err_out nam_job samp

Description:
  Run 'atria' with explicitly passed argument values and write stdout/stderr logs to sample-specific files.

  Log files are written to:

    \${err_out}/\${nam_job}.\${samp}.stdout.txt
    \${err_out}/\${nam_job}.\${samp}.stderr.txt

Positional arguments:
  1  threads  <int>  Number of threads to use.
  2  fq_1     <str>  FASTQ file #1.
  3  fq_2     <str>  FASTQ file #2, or 'NA' for SE data.
  4  dir_out  <str>  Output directory for processed FASTQ files.
  5  err_out  <str>  Directory for stderr/stdout log files.
  6  nam_job  <str>  Job name used in log-file naming.
  7  samp     <str>  Sample name used in log-file naming.

Notes:
  - This helper is a thin wrapper around 'atria'.
  - '-R <fq_2>' is passed only when 'fq_2 != NA'.
  - The current wrapper always passes '--length-range 35:500'.
EOM
    )

    if [[ "${threads}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${threads}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'threads', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  Define paths for log output files
    log_out="${err_out}/${nam_job}.${samp}.stdout.txt"
    log_err="${err_out}/${nam_job}.${samp}.stderr.txt"

    #  Run Atria with the specified arguments
    # shellcheck disable=SC2046
    if ! \
        atria \
            -t "${threads}" \
            -r "${fq_1}" \
            $(if [[ "${fq_2}" != "NA" ]]; then echo "-R ${fq_2}"; fi) \
            -o "${dir_out}" \
            --length-range 35:500 \
                 > "${log_out}" \
                2> "${log_err}"
    then
        echo_err_func "${FUNCNAME[0]}" \
            "read trimming failed for sample '${samp}'. See log: '${log_err}'."
        return 1
    fi
    #TODO: 'fq_2' handling is vulnerable to word splitting if paths contain
    #+     spaces
}
#MAYBE: spin the call to 'atria' out into its own separate function (consistent
#+      with other scripts, e.g., 'align_fastqs.sh::align_fastqs')?
#MAYBE: maybe 'align_fastqs.sh' should become, e.g., 'process_fastqs.sh' and
#+      'run_atria' should be included in there?


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
function source_submit_helpers() {
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


#  Initialize argument variables, assigning default values where applicable
env_nam="env_protocol"
dir_scr=""
threads=4
csv_infile=""
dir_out=""
sfx_se=""
sfx_pe=""
err_out=""
nam_job="trim_fastqs"


function show_help_main() {
    cat << EOM >&2
Usage:
  submit_trim_fastqs.sh
    [-h|--hlp|--help] [-en|--env_nam <str>] -ds|--dir_scr <str> [-t|--threads <int>]
    -i|--csv_infile <str> -do|--dir_out <str>
    -sxs|--sfx_se <str> -sxp|--sfx_pe <str> -eo|--err_out <str> [-nj|--nam_job <str>]

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


#  Main script execution
function main() {
    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        show_help_main
        echo >&2
        exit 0
    fi

    dir_scr="$(resolve_dir_scr "${0##*/}" "$@")" || exit 1

    source_submit_helpers "${0##*/}" "${dir_scr}" \
        check_args \
        check_inputs \
        check_numbers \
        format_outputs \
        handle_env \
        manage_slurm \
        process_sequences \
        || exit 1

    parse_args "$@"  || exit 1
    validate_args    || exit 1
    print_debug_args || exit 1

    #  Activate environment
    handle_env "${env_nam}" || exit 1

    IFS=';' read -r -a arr_infile <<< "${csv_infile}"
    check_arr_nonempty "arr_infile" "csv_infile" || exit 1

    if [[ "${debug}" == "true" ]]; then
        echo "\${#arr_infile[@]}=${#arr_infile[@]}" && echo
        echo "arr_infile=( ${arr_infile[*]} )"      && echo
    fi

    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        id_job="${SLURM_ARRAY_JOB_ID}"
        id_tsk="${SLURM_ARRAY_TASK_ID}"

        if ! [[ "${id_tsk}" =~ ^[1-9][0-9]*$ ]]; then
            echo_err "Slurm task ID is invalid: '${id_tsk}'."
            exit 1
        elif (( id_tsk > ${#arr_infile[@]} )); then
            echo_err \
                "Slurm task ID '${id_tsk}' exceeds number of FASTQ entries:" \
                "'${#arr_infile[@]}'."
            exit 1
        else
            idx=$(( id_tsk - 1 ))
        fi

        infile="${arr_infile[idx]}"

        if [[ "${debug}" == "true" ]]; then debug_var "infile=${infile}"; fi

        IFS=';' read -r fq_1 fq_2 samp < <(
            parse_entry_trim_fastq "${infile}" "${sfx_se}" "${sfx_pe}"
        ) || exit 1

        if [[ "${debug}" == "true" ]]; then
            debug_var \
                "fq_1=${fq_1}" \
                "fq_2=${fq_2}" \
                "samp=${samp}"
        fi

        IFS=',' read -r err_ini out_ini err_dsc out_dsc < <(
            set_logs_slurm \
                "${id_job}" "${id_tsk}" "${samp}" "${err_out}" "${nam_job}"
        ) || exit 1

        if [[ "${debug}" == "true" ]]; then
            debug_var \
                "err_ini=${err_ini}" \
                "out_ini=${out_ini}" \
                "err_dsc=${err_dsc}" \
                "out_dsc=${out_dsc}"
        fi

        if ! \
            run_atria \
                "${threads}" "${fq_1}" "${fq_2}" "${dir_out}" "${err_out}" \
                "${nam_job}" "${samp}"
        then
            echo_err "failed to perform read trimming."
            exit 1
        fi

        rm "${err_ini}" "${out_ini}"
    else
        for idx in "${!arr_infile[@]}"; do
            infile="${arr_infile[idx]}"

            if [[ "${debug}" == "true" ]]; then debug_var "infile=${infile}"; fi

            IFS=';' read -r fq_1 fq_2 samp < <(
                parse_entry_trim_fastq "${infile}" "${sfx_se}" "${sfx_pe}"
            ) || exit 1

            if [[ "${debug}" == "true" ]]; then
                debug_var \
                    "fq_1=${fq_1}" \
                    "fq_2=${fq_2}" \
                    "samp=${samp}"
            fi

            if ! \
                run_atria \
                    "${threads}" "${fq_1}" "${fq_2}" "${dir_out}" "${err_out}" \
                    "${nam_job}" "${samp}"
            then
                echo_err "failed to perform read trimming."
                exit 1
            fi
        done
    fi
}


main "$@"
