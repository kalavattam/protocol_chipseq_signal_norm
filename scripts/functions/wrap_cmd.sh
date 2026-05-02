#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: wrap_cmd.sh
#
# Copyright 2025-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5-series models) was used in development.
#
# Distributed under the MIT license.


# get_submit_logs
# print_built_cmd


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
    _dir_src_cmd="$(
        cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
    )"

    source "${_dir_src_cmd}/source_helpers.sh" || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source '${_dir_src_cmd}/source_helpers.sh'." >&2

        if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
            return 1
        else
            exit 1
        fi
    }

    source_helpers "${_dir_src_cmd}" \
        check_source format_outputs || {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to source required helper dependencies." >&2

            if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
                return 1
            else
                exit 1
            fi
        }

    unset _dir_src_cmd
}


# shellcheck disable=SC2154
function get_submit_logs() {
    local fil_smp="${1:-}"
    local base log_out log_err
    local show_help

    show_help=$(cat << EOM
Usage:
  get_submit_logs [-h|--hlp|--help] fil_smp

Description:
  Derive stdout/stderr log paths for one local per-sample call to 'submit_*.sh'.

Positional arguments:
  1  fil_smp  <str>  Input file path used to derive the log-file basename.

Expected globals:
  err_out  Existing directory in which wrapper-level stdout/stderr logs are written.
  nam_job  Job-name prefix used in derived log filenames.

Returns:
  0 if log paths are derived successfully; otherwise 1.

Output:
  Prints one semicolon-delimited line:

      log_out;log_err
EOM
    )

    if [[ "${fil_smp}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${fil_smp}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'fil_smp', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    if [[ -z "${err_out:-}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "global variable 'err_out' is empty or unset."
        return 1
    elif [[ ! -d "${err_out}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "global variable 'err_out' is not an existing directory:" \
            "'${err_out}'."
        return 1
    elif [[ -z "${nam_job:-}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "global variable 'nam_job' is empty or unset."
        return 1
    fi

    base="$(basename "${fil_smp}")"
    base="${base%.gz}"
    base="${base%.bam}"
    base="${base%.cram}"
    base="${base%.sam}"
    base="${base%.bedGraph}"
    base="${base%.bedgraph}"
    base="${base%.bdg}"
    base="${base%.bg}"
    base="${base%.bed}"

    log_out="${err_out}/${nam_job}.${base}.stdout.txt"
    log_err="${err_out}/${nam_job}.${base}.stderr.txt"

    printf '%s;%s\n' "${log_out}" "${log_err}"
}


#MAYBE: rename to 'print_cmd_built'?
# shellcheck disable=SC2154
function print_built_cmd() {
    local log_out="${1:-}"
    local log_err="${2:-}"
    local show_help

    show_help=$(cat << EOM
Usage:
  print_built_cmd [-h|--hlp|--help] [log_out] [log_err]

Description:
  Print the current command array 'cmd_bld' as one shell-escaped command line.

  If both 'log_out' and 'log_err' are supplied, append stdout and stderr redirections to the printed command.

Positional arguments:
  1  log_out  <str>  Optional stdout log path to append as '>> log_out'.
  2  log_err  <str>  Optional stderr log path to append as '2>> log_err'.

Expected global:
  cmd_bld  Global indexed array containing the command to print; must be set and non-empty.

Returns:
  0 if the command is printed successfully; otherwise 1.

Notes:
  - Redirections are appended only when both 'log_out' and 'log_err' are supplied.
  - Output is shell-escaped via 'printf %q'.
EOM
    )

    if [[ "${log_out}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -n "${log_err}" && -z "${log_out}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'log_err', was supplied without" \
            "positional argument 1, 'log_out'."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -n "${log_out}" && -z "${log_err}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'log_out', was supplied without" \
            "positional argument 2, 'log_err'."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    if ! declare -p cmd_bld > /dev/null 2>&1; then
        echo_err_func "${FUNCNAME[0]}" \
            "global indexed array 'cmd_bld' is unset."
        return 1
    fi

    if [[ "$(declare -p cmd_bld 2> /dev/null)" != declare\ -a* ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'cmd_bld' must be an indexed array."
        return 1
    fi

    if (( ${#cmd_bld[@]} == 0 )); then
        echo_err_func "${FUNCNAME[0]}" \
            "global indexed array 'cmd_bld' is empty."
        return 1
    fi

    printf '%q ' "${cmd_bld[@]}"

    if [[ -n "${log_out}" && -n "${log_err}" ]]; then
        printf '>> %q 2>> %q' "${log_out}" "${log_err}"
    fi

    echo
}


#  Print an error message when function script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    err_source_only "${BASH_SOURCE[0]}"
fi
