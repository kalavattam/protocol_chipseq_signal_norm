#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: check_env.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
#
# Distributed under the MIT license.


# check_env_installed
# check_pgrm_path


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
    _dir_src_env_chk="$(
        cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
    )"

    source "${_dir_src_env_chk}/source_helpers.sh" || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source '${_dir_src_env_chk}/source_helpers.sh'." >&2

        if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
            return 1
        else
            exit 1
        fi
    }

    source_helpers "${_dir_src_env_chk}" \
        check_source format_outputs || {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to source required helper dependencies." >&2

            if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
                return 1
            else
                exit 1
            fi
        }

    unset _dir_src_env_chk
}


function check_env_installed() {
    local env_nam="${1:-}"
    local quiet="${2:-false}"
    local quiet_lc
    local show_help

    show_help=$(cat << EOM
Usage:
  check_env_installed [-h|--hlp|--help] env_nam [quiet]

Description:
  Check that a specific Conda (or Mamba) environment is installed.

Positional arguments:
  1  env_nam  <str>  Name of Conda (or Mamba) environment to check.
  2  quiet    <bol>  If 'true' or 't', suppress expected "not installed" error messages and return status only (default: 'false').

Returns:
  0 if the environment is installed; otherwise 1 and, unless 'quiet' is true-like, an error message.

Dependencies:
  - Bash >= 5
  - Conda

Notes:
  - Quiet mode does not suppress the "'conda' is not available in PATH" message.
  - Quiet mode accepts 'true', 't', 'false', or 'f' in any letter case.

Examples:
  1. Check that environment "env_protocol" is installed
  '''bash
  check_env_installed "env_protocol"  # Returns 0
  '''

  2. Confirm that fictitious environment "env_nonexistent" errors
  '''bash
  check_env_installed "env_nonexistent"  # Returns 1 and error message
  '''

  3. Confirm that fictitious environment "env_nonexistent" errors quietly
  '''bash
  check_env_installed "env_nonexistent" "true"  # Returns 1
  '''
EOM
    )

    #  Parse and check function arguments
    if [[ "${env_nam}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${env_nam}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'env_nam', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ "${env_nam}" == -* ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "unknown option '${env_nam}'. This function accepts only" \
            "'-h', '--hlp', or '--help' as options; otherwise supply" \
            "'env_nam' as positional argument 1."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  Lowercase-convert 'quiet' assignment
    quiet_lc="${quiet,,}"

    case "${quiet_lc}" in
        t|true)  quiet=true  ;;
        f|false) quiet=false ;;
        *)
            echo_err_func "${FUNCNAME[0]}" \
                "positional argument 2, 'quiet', must be Boolean-like" \
                "'true', 't', 'false', or 'f': '${quiet}'."
            echo >&2
            echo "${show_help}" >&2
            return 1
            ;;
    esac

    #  Check availability of 'conda'
    if ! command -v conda >/dev/null 2>&1; then
        echo_err_func "${FUNCNAME[0]}" \
            "'conda' is not available in PATH."
        return 1
    fi

    #  Check that the specific Conda or Mamba environment is installed
    if ! \
        conda env list | awk -v env="${env_nam}" '
            $1 == env { found = 1 } END { exit !found }
        '
    then
        if [[ "${quiet}" == "false" ]]; then
            echo_err_func "${FUNCNAME[0]}" \
                "environment '${env_nam}' appears not to be installed."
        fi
        return 1
    fi
}


function check_pgrm_path() {
    local prog="${1:-}"
    local show_help

    show_help=$(cat << EOM
Usage:
  check_pgrm_path [-h|--hlp|--help] prog

Description:
  Checks that a given program is available in PATH.

Positional argument:
  1  prog  <str>  Name of program to check.

Returns:
  0 if program is found in PATH; otherwise, 1 with an error message.

Example:
  1. Confirm success when program 'samtools' is in PATH
  '''bash
  check_pgrm_path "samtools"  # Returns 0
  '''

  2. Confirm error when program 'samtools' is not in PATH
  '''bash
  check_pgrm_path "samtools"  # Returns 1 and error message
  '''
EOM
    )

    #  Parse and check function arguments
    if [[ "${prog}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${prog}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'prog', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ "${prog}" == -* ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "unknown option '${prog}'. This function accepts only '-h'," \
            "'--hlp', or '--help' as options; otherwise supply 'prog' as" \
            "positional argument 1."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    if ! command -v "${prog}" &> /dev/null; then
        echo_err_func "${FUNCNAME[0]}" \
            "'${prog}' is not in PATH. Please install '${prog}' and/or add" \
            "it to PATH."
        return 1
    fi
}


#  Print an error message when function script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    err_source_only "${BASH_SOURCE[0]}"
fi
