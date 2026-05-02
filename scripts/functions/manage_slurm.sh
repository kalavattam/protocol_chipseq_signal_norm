#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: manage_slurm.sh
#
# Copyright 2025-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5-series models) was used in development.
#
# Distributed under the MIT license.


# set_logs_slurm


#  Require Bash >= 4.4 before defining functions
if [[ -z "${BASH_VERSION:-}" ]]; then
    echo "error(shell):" \
        "this script must be sourced or run under Bash >= 4.4." >&2

    if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
        return 1
    else
        exit 1
    fi
elif ((
    BASH_VERSINFO[0] < 4 || ( BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] < 4 )
)); then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "this script requires Bash >= 4.4; current version is" \
        "'${BASH_VERSION}'." >&2

    if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
        return 1
    else
        exit 1
    fi
fi

#  Source required helper functions if needed
# shellcheck disable=SC1091
{
    _dir_src_slurm="$(
        cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
    )"

    source "${_dir_src_slurm}/source_helpers.sh" || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source '${_dir_src_slurm}/source_helpers.sh'." >&2

        if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
            return 1
        else
            exit 1
        fi
    }

    source_helpers "${_dir_src_slurm}" \
        check_inputs check_source format_outputs || {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to source required helper dependencies." >&2

            if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
                return 1
            else
                exit 1
            fi
        }

    unset _dir_src_slurm
}


#  Derive Slurm log-file paths and create descriptive hard-linked log names
function set_logs_slurm() {
    local id_job="${1:-}"    # Slurm job ID
    local id_tsk="${2:-}"    # Slurm task ID within the job array
    local samp="${3:-}"      # Sample name associated with log files
    local err_out="${4:-}"   # Directory for stderr and stdout log files
    local nam_job="${5:-}"   # Name of the job (used in log file naming)
    local err_ini out_ini    # Slurm initial stderr and stdout log files
    local err_dsc out_dsc    # Hard-linked stderr/stdout with descriptive names
    local show_help          # Help message

    show_help=$(cat << EOM
Usage:
  set_logs_slurm
    [-h|--hlp|--help] id_job id_tsk samp err_out nam_job

Description:
  Derive initial and descriptive Slurm log-file paths and create hard-linked descriptive log names.

Positional arguments:
  1  id_job   <int>  Slurm job ID.
  2  id_tsk   <int>  Slurm task ID within the job array.
  3  samp     <str>  Sample name associated with the log files.
  4  err_out  <str>  Directory for stderr/stdout log files.
  5  nam_job  <str>  Job name used in log-file naming.

Returns:
  Prints a comma-delimited record to stdout:

    err_ini,out_ini,err_dsc,out_dsc

Notes:
  - Initial log names follow the form:
      \${err_out}/\${nam_job}.\${id_job}-\${id_tsk}.stderr.txt
      \${err_out}/\${nam_job}.\${id_job}-\${id_tsk}.stdout.txt
  - Descriptive hard-linked names additionally include '\${samp}'.
EOM
    )

    if [[ "${id_job}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${id_job}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'id_job', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${id_tsk}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'id_tsk', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${samp}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 3, 'samp', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${err_out}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 4, 'err_out', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${nam_job}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 5, 'nam_job', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    validate_var_dir "err_out" "${err_out}" || return 1

    #  Set log paths
    err_ini="${err_out}/${nam_job}.${id_job}-${id_tsk}.stderr.txt"
    out_ini="${err_out}/${nam_job}.${id_job}-${id_tsk}.stdout.txt"
    err_dsc="${err_out}/${nam_job}.${samp}.${id_job}-${id_tsk}.stderr.txt"
    out_dsc="${err_out}/${nam_job}.${samp}.${id_job}-${id_tsk}.stdout.txt"

    #  Create hard-linked log files
    if ! \
        ln -f "${err_ini}" "${err_dsc}"
    then
        echo_err_func "${FUNCNAME[0]}" \
            "failed to create hard link: '${err_ini}' to '${err_dsc}'."
        return 1
    fi

    if ! \
        ln -f "${out_ini}" "${out_dsc}"
    then
        echo_err_func "${FUNCNAME[0]}" \
            "failed to create hard link: '${out_ini}' to '${out_dsc}'."
        return 1
    fi

    #  Return values
    echo "${err_ini},${out_ini},${err_dsc},${out_dsc}"
}


#  Print an error message when function script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    err_source_only "${BASH_SOURCE[0]}"
fi
