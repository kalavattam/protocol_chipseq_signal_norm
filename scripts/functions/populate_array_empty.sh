#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: populate_array_empty.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
#
# Distributed under the MIT license.


# populate_array_empty


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
    _dir_src_pop="$(
        cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
    )"

    source "${_dir_src_pop}/source_helpers.sh" || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source '${_dir_src_pop}/source_helpers.sh'." >&2

        if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
            return 1
        else
            exit 1
        fi
    }

    source_helpers "${_dir_src_pop}" \
        check_source format_outputs || {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to source required helper dependencies." >&2

            if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
                return 1
            else
                exit 1
            fi
        }

    unset _dir_src_pop
}


# shellcheck disable=SC2034
function populate_array_empty() {
    local arr_nam="${1:-}"  # Name of target array variable (string)
    local len_tgt="${2:-}"  # Desired length if empty ('int >= 0')
    local fill="${3:-NA}"   # Value to insert
    local i                 # Loop index
    local len_cur           # Current length of target array
    local decl              # Output of 'declare -p' for target variable
    local show_help         # Help message

    show_help=$(cat << EOM
Usage:
  populate_array_empty [-h|--hlp|--help] arr_nam len_tgt [fill]

Description:
  Populate a named indexed array with '\${fill}' (default: "NA") only if it is currently empty.
    - If the named variable does not exist, it is initialized as an empty indexed array.
    - If the named variable already exists, it must be an indexed array.
    - Emptiness is determined from the array length.
    - Existing, non-empty arrays are left unchanged.

Positional arguments:
  1  arr_nam  <str>  Name of the target array variable (not a reference).
  2  len_tgt  <int>  Desired length to populate (must be an integer >= 0).
  3  fill     <str>  Value to write into each element (default: "NA").

Returns:
  0 on success, non-0 on validation or runtime error.

Notes:
  - Intended for Bash indexed arrays.
  - Preserves any pre-existing contents if the array is already non-empty.

Examples:
  '''bash
  #  If 'arr_foo' is unset or empty, make it 4 × "NA"
  populate_array_empty "arr_foo" 4

  #  If 'arr_bar' is unset or empty, make it 3 × "0"
  populate_array_empty "arr_bar" 3 "0"
  '''

#TODO:
  - What if the array has elements but not the correct number?
EOM
    )

    #  Parse and check function arguments
    if [[ "${arr_nam}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${arr_nam}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'arr_nam', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${len_tgt}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'len_tgt', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  Check that target length is int >= 0
    if ! [[ "${len_tgt}" =~ ^[0-9]+$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'len_tgt', must be a non-negative" \
            "integer: '${len_tgt}'."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  Check that 'arr_nam' is a valid Bash variable name
    if ! [[ "${arr_nam}" =~ ^[a-zA-Z_][a-zA-Z0-9_]*$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'arr_nam', must be a valid shell" \
            "variable name: '${arr_nam}'."
        return 1
    fi

    #  Ensure the target exists and is an indexed array
    if ! \
        decl="$(declare -p "${arr_nam}" 2> /dev/null)"
    then
        declare -g -a "${arr_nam}"
        len_cur=0
    elif [[ "${decl}" != declare\ -a* ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'arr_nam', must name an indexed array:" \
            "'${arr_nam}'."
        return 1
    elif [[ "${decl}" != *"["* ]]; then
        len_cur=0
    else
        eval "len_cur=\${#${arr_nam}[@]}"
    fi

    #  Only populate if the array is empty
    if (( len_cur == 0 )); then
        for (( i = 0; i < len_tgt; i++ )); do
            printf -v "${arr_nam}[${i}]" '%s' "${fill}"
        done
    fi
}


#  Print an error message when function script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    err_source_only "${BASH_SOURCE[0]}"
fi
