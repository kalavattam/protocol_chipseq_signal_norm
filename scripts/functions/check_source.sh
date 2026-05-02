#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: check_source.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.4) was used in development.
#
# Distributed under the MIT license.


# err_source_only


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

# _msg="this script must be run under Bash"
#
#
# function require_bash() {
#     if [[ -z "${BASH_VERSION:-}" ]]; then
#         echo "error(shell):" \
#             "${_msg}." >&2
#         return 1
#     fi
# }
#
#
# function require_bash_5() {
#     if [[ -z "${BASH_VERSION:-}" ]]; then
#         echo "error(shell):" \
#             "${_msg} >= 5." >&2
#         return 1
#     fi
#
#     if (( BASH_VERSINFO[0] < 5 )); then
#         echo "error(shell):" \
#             "${_msg} >= 5." >&2
#         return 1
#     fi
# }


function err_source_only() {
    local pth_scr="${1:-}"
    local nam_scr
    local show_help

    show_help=$(cat << EOM
Usage:
  err_source_only [-h|--hlp|--help] pth_scr

Description:
  Print a warning if the supplied script path is being executed directly rather than being sourced.

Positional argument:
  1  pth_scr  <str>  Script path, typically passed as "\${BASH_SOURCE[0]}".

Returns:
  0 in all cases. Prints a warning to stderr only when the script is being executed directly.

Example:
  '''bash
  err_source_only "\${BASH_SOURCE[0]}"
  '''
EOM
    )

    if [[ "${pth_scr}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${pth_scr}" ]]; then
        echo "error(check_source.sh::err_source_only):" \
            "positional argument 1, 'pth_scr' is missing." >&2
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    nam_scr="$(basename "${pth_scr}")"

    if [[ "${pth_scr}" == "${0}" ]]; then
        echo "error(${nam_scr}):" \
            "this script is intended to be sourced before use; do not run it" \
            "as, e.g., 'bash ${pth_scr}'." >&2
    fi

    return 0
}


#  Print an error message when function script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    err_source_only "${BASH_SOURCE[0]}"
fi
