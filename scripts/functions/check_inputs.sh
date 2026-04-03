#!/bin/bash
# -*- coding: utf-8 -*-
# 
# Script: check_inputs.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
# 
# Distributed under the MIT license.


#  Check that all supplied file paths exist
function check_array_files() {
    local desc="${1:-}"
    local file

    #  Check that 'desc' input is not empty
    if [[ -z "${desc}" ]]; then
        echo "Error: Positional argument 1, 'desc', is required." >&2
        return 1
    fi

    shift

    #  Check that files were supplied
    if (( $# < 1 )); then
        echo "Error: No files supplied to validate for '${desc}'." >&2
        return 1
    fi
    
    #  Check each supplied file path
    for file in "$@"; do
        if [[ ! -f "${file}" ]]; then
            echo "Error: '${desc}' file does not exist: '${file}'." >&2
            return 1
        fi
    done

    return 0
}


#  Check that two indexed arrays have matching lengths (usable with Bash >=3.2)
function check_arrays_lengths() {
    local arr_nam_1="${1:-}"
    local arr_nam_2="${2:-}"
    local arr_siz_1 arr_siz_2

    #  Ensure array names are valid
    if [[
           ! "${arr_nam_1}" =~ ^[a-zA-Z_][a-zA-Z0-9_]*$ \
        || ! "${arr_nam_2}" =~ ^[a-zA-Z_][a-zA-Z0-9_]*$
    ]]; then
        echo "Error: Invalid array name '${arr_nam_1}' and/or '${arr_nam_2}'." >&2
        return 1
    fi

    #  Ensure both names refer to arrays
    if ! \
        eval 'declare -p '"${arr_nam_1}"' &>/dev/null'
    then
        echo "Error: '${arr_nam_1}' is unset." >&2
        return 1
    fi

    if ! \
        eval 'declare -p '"${arr_nam_2}"' &>/dev/null'
    then
        echo "Error: '${arr_nam_2}' is unset." >&2
        return 1
    fi

    if ! \
        eval 'declare -p '"${arr_nam_1}"' 2>/dev/null | grep -q "^declare \-a "'
    then
        echo "Error: '${arr_nam_1}' is not an indexed array." >&2
        return 1
    fi

    if ! \
        eval 'declare -p '"${arr_nam_2}"' 2>/dev/null | grep -q "^declare \-a "'
    then
        echo "Error: '${arr_nam_2}' is not an indexed array." >&2
        return 1
    fi

    #  Get array sizes using indirect references
    eval "arr_siz_1=\${#${arr_nam_1}[@]}"
    eval "arr_siz_2=\${#${arr_nam_2}[@]}"

    #  Check that array sizes match
    # shellcheck disable=SC2154
    if [[ "${arr_siz_1}" -ne "${arr_siz_2}" ]]; then
        echo \
            "Error: Array length mismatch. '${arr_nam_1}' has '${arr_siz_1}'" \
            "element(s), whereas '${arr_nam_2}' has '${arr_siz_2}'." >&2
        return 1
    fi

    return 0
}


function check_file_dir_exists() {
    local type="${1}"
    local item="${2}"
    local name="${3}"
    local item_type
    local check_flag
    local not_exist_msg
    local show_help

show_help=$(cat << EOM
Usage:
  check_file_dir_exists type item [name]

Description:
  Check the existence of a file or directory based on the provided type.

Positional arguments:
  1  type  <str>  The type to check for existence. Use 'f' for file or 'd' for directory (required).
  2  item  <str>  The file or directory, including its path, to check (required).
  3  name  <str>  The name to associate with the item for error messages (optional).

Returns:
  0 if the file or directory exists; otherwise, 1 and an error message.

Dependencies:
  - Bash

Examples:
  '''bash
  #  Check if a file exists
  check_file_dir_exists "f" "/path/to/file.txt"

  #  Check if a directory exists
  check_file_dir_exists "d" "/path/to/directory"

  #  Check if a file exists with an associated name for error message
  check_file_dir_exists "f" "/path/to/file.txt" "infile"

  #  Check if a directory exists with an associated name for error message
  check_file_dir_exists "d" "/path/to/directory" "outdir"
  '''
EOM
    )

    #  Parse and check function parameters
    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    #  Validate variable type
    if [[ "${type}" != "f" && "${type}" != "d" ]]; then
        echo \
            "Error: Positionl parameter 1, 'type', is invalid: '${type}'." \
            "Expected 'f' for file or 'd' for directory." >&2
        echo "" >&2
        return 1
    fi

    #  Check that variable item is defined and not empty
    if [[ -z "${item}" ]]; then
        echo \
            "Error: Positional parameter 2, 'item', is not defined or is" \
            "empty." >&2
        echo "" >&2
        return 1
    fi

    #  Check file or directory existence; construct and return error message if
    #+ applicable
    if [[ "${type}" == "f" ]]; then
        item_type="File"
        check_flag="-f"
        not_exist_msg="File does not exist"
    elif [[ "${type}" == "d" ]]; then
        item_type="Directory"
        check_flag="-d"
        not_exist_msg="Directory does not exist"
    fi

    if [[ -n "${name}" ]]; then
        not_exist_msg="${item_type} associated with '--${name}' does not exist"
    fi

    if ! eval "[[ ${check_flag} \"${item}\" ]]"; then
        echo "Error: ${not_exist_msg}: '${item}'." >&2
        return 1
    fi
}


#  Debug array contents (usable with Bash >=3.2)
function debug_array_contents() {
    local arr_nam
    local -a arr

    for arr_nam in "$@"; do
        #  Skip invalid variable names
        if ! [[ "${arr_nam}" =~ ^[a-zA-Z_][a-zA-Z0-9_]*$ ]]; then
            continue
        fi

        #  Check array exists (if unset/not an array, behavior can get messy)
        if ! eval 'declare -p '"${arr_nam}"' >/dev/null 2>&1'; then
            continue
        fi

        #  Skip non-indexed variables
        if ! eval '
            declare -p '"${arr_nam}"' 2>/dev/null | grep -q "^declare \-a "
        '; then
            continue
        fi

        #  Access the array indirectly using eval
        eval "arr=( \"\${${arr_nam}[@]}\" )"
        if [[ -n "${arr[*]}" ]]; then
            echo "  - ${arr_nam}=( ${arr[*]} )"
        fi
    done
}
