#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: check_inputs.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
# GPT-5.6) were used in development and documentation.
#
# Distributed under the MIT license.


# validate_var
# validate_file
# validate_dir
# validate_var_file
# validate_var_dir
# debug_var
# check_arr_files
# check_arr_lengths
# check_file_dir_exists
# debug_arr_contents
# check_arr_nonempty
# check_arr_len_bcst


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
{
    _dir_src_in="$(
        cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
    )"

    # shellcheck source=lib/bash/core/source_helpers.sh
    source "${_dir_src_in}/source_helpers.sh" || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source '${_dir_src_in}/source_helpers.sh'." >&2

        if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
            return 1
        else
            exit 1
        fi
    }

    source_helpers "${_dir_src_in}" \
        check_source format_outputs || {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to source required helper dependencies." >&2

            if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
                return 1
            else
                exit 1
            fi
        }

    unset _dir_src_in
}


function validate_var() {
    local var_nam="${1:-}"  # Variable name (for error messages)
    local var_val="${2:-}"  # Value to check for emptiness/unset state
    local idx="${3:-0}"     # Optional index (for arrays); defaults to '0'
    local show_help         # Help message

    show_help=$(cat << EOM
Usage
-----
  validate_var
    [--help] var_nam var_val [idx]

  Check that a variable value is non-empty.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  var_nam : str
    Variable name to use in error messages.

  2  var_val : str
    Value to check.

  3  idx : int
    Optional array index for error messages (default: 0).

Returns
-------
  0 if 'var_val' is non-empty; otherwise 1.

Notes
-----
  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. After sourcing this file, accept a nonempty scalar value.
    '''bash
    validate_var "mode" "signal"
    '''

  2. Reject an empty value at a reported array index.
    '''bash
    if ! validate_var "csv_fil_in" "" 2; then
      echo "csv_fil_in element 2 is empty" >&2
    fi
    '''
EOM
    )

    if [[ "${var_nam}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${var_nam}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'var_nam', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    if [[ -z "${var_val}" ]]; then
        if [[ "${idx}" -ne 0 ]]; then
            echo_err_func "${FUNCNAME[0]}" \
                "'${var_nam}' is empty or unset for array index '${idx}'."
        else
            echo_err_func "${FUNCNAME[0]}" \
                "'${var_nam}' is empty or unset."
        fi
        return 1
    fi
}


function validate_file() {
    local var_nam="${1:-}"    # Variable name (for error messages)
    local var_val="${2:-}"    # Value (file) to check for existence
    local idx="${3:-0}"       # Optional index (for arrays); defaults to '0'
    local chk_r="${4:-false}" # If true, require file to be readable
    local show_help           # Help message

    show_help=$(cat << EOM
Usage
-----
  validate_file
    [--help] var_nam var_val [idx] [chk_r]

  Check that a file path assigned to a variable exists as a regular file.
  Optionally check that the file is readable.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  var_nam : str
    Variable name to use in error messages.

  2  var_val : file
    File path to check.

  3  idx : int
    Optional array index for error messages (default: 0).

  4  chk_r : bool
    Whether to require readability: case-insensitive Boolean-like true|t|yes|y|1 or false|f|no|n|0 (default: false).

Returns
-------
  0 if 'var_val' exists as a regular file and, when 'chk_r=true', is readable; otherwise 1.

Notes
-----
  Runtime requirements:
    bash >= 4.4

  Accepted true-like values: true, t, yes, y, 1. Accepted false-like values: false, f, no, n, 0. Matching is case-insensitive; surrounding whitespace and empty values are invalid. Successful normalization emits 'true' or 'false'.

Examples
--------
  1. Require the sourced helper file to exist and be readable.
    '''bash
    validate_file "helper" "\${BASH_SOURCE[0]}" 0 "true"
    '''

  2. Reject a path that is not a regular file.
    '''bash
    if ! validate_file "fil_in" "/path/that/does/not/exist"; then
      echo "fil_in is unavailable" >&2
    fi
    '''
EOM
    )

    if [[ "${var_nam}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${var_nam}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'var_nam', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${var_val}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'var_val', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    chk_r="$(normalize_bool "${chk_r}" "chk_r")" || return 1

    if [[ ! -f "${var_val}" ]]; then
        if [[ "${idx}" -ne 0 ]]; then
            echo_err_func "${FUNCNAME[0]}" \
                "'${var_nam}' does not exist as a regular file for array" \
                "index '${idx}'."
        else
            echo_err_func "${FUNCNAME[0]}" \
                "'${var_nam}' does not exist as a regular file."
        fi
        return 1
    elif [[ "${chk_r}" == "true" && ! -r "${var_val}" ]]; then
        if [[ "${idx}" -ne 0 ]]; then
            echo_err_func "${FUNCNAME[0]}" \
                "'${var_nam}' is not readable for array index '${idx}'."
        else
            echo_err_func "${FUNCNAME[0]}" \
                "'${var_nam}' is not readable."
        fi
        return 1
    fi
}


function validate_dir() {
    local var_nam="${1:-}"    # Variable name (for error messages)
    local var_val="${2:-}"    # Value (directory) to check
    local idx="${3:-0}"       # Optional index (for arrays); defaults to '0'
    local chk_w="${4:-true}"  # If true, require directory to be writable
    local show_help           # Help message

    show_help=$(cat << EOM
Usage
-----
  validate_dir
    [--help] var_nam var_val [idx] [chk_w]

  Check that a directory path assigned to a variable exists and, by default, is writable.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  var_nam : str
    Variable name to use in error messages.

  2  var_val : dir
    Directory path to check.

  3  idx : int
    Optional array index for error messages (default: 0).

  4  chk_w : bool
    Whether to require writability: case-insensitive Boolean-like true|t|yes|y|1 or false|f|no|n|0 (default: true).

Returns
-------
  0 if 'var_val' exists as a directory and, when 'chk_w=true', is writable; otherwise 1.

Notes
-----
  Runtime requirements:
    bash >= 4.4

  Accepted true-like values: true, t, yes, y, 1. Accepted false-like values: false, f, no, n, 0. Matching is case-insensitive; surrounding whitespace and empty values are invalid. Successful normalization emits 'true' or 'false'.

Examples
--------
  1. Require the current directory to exist and be writable.
    '''bash
    validate_dir "dir_out" "\${PWD}" 0 "true"
    '''

  2. Reject a path that is not a directory.
    '''bash
    if ! validate_dir "dir_out" "/path/that/does/not/exist"; then
      echo "dir_out is unavailable" >&2
    fi
    '''
EOM
    )

    if [[ "${var_nam}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${var_nam}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'var_nam', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${var_val}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'var_val', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    chk_w="$(normalize_bool "${chk_w}" "chk_w")" || return 1

    if [[ ! -d "${var_val}" ]]; then
        if [[ "${idx}" -ne 0 ]]; then
            echo_err_func "${FUNCNAME[0]}" \
                "directory for '${var_nam}' does not exist for array index" \
                "'${idx}'."
        else
            echo_err_func "${FUNCNAME[0]}" \
                "directory for '${var_nam}' does not exist."
        fi
        return 1
    elif [[ "${chk_w}" == "true" && ! -w "${var_val}" ]]; then
        if [[ "${idx}" -ne 0 ]]; then
            echo_err_func "${FUNCNAME[0]}" \
                "directory for '${var_nam}' is not writable for array index" \
                "'${idx}'."
        else
            echo_err_func "${FUNCNAME[0]}" \
                "directory for '${var_nam}' is not writable."
        fi
        return 1
    fi
}


function validate_var_file() {
    local var_nam="${1:-}"    # Variable name to use in error messages
    local pth_fil="${2:-}"    # File path to validate
    local idx="${3:-0}"       # Optional index (for arrays); defaults to '0'
    local chk_r="${4:-false}" # If true, require file to be readable
    local show_help           # Help message

    show_help=$(cat << EOM
Usage
-----
  validate_var_file
    [--help] var_nam pth_fil [idx] [chk_r]

  Check that a variable value is set/non-empty and points to an existing
  regular file. Optionally check that the file is readable.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  var_nam : str
    Variable name to use in error messages.

  2  pth_fil : file
    File path to validate.

  3  idx : int
    Optional array index for error messages (default: 0).

  4  chk_r : bool
    Whether to require readability: case-insensitive Boolean-like true|t|yes|y|1 or false|f|no|n|0 (default: false).

Returns
-------
  0 if the value is non-empty and the file exists as a regular file and, when 'chk_r=true', is readable; otherwise 1.

Notes
-----
  Runtime requirements:
    bash >= 4.4

  Accepted true-like values: true, t, yes, y, 1. Accepted false-like values: false, f, no, n, 0. Matching is case-insensitive; surrounding whitespace and empty values are invalid. Successful normalization emits 'true' or 'false'.

Examples
--------
  1. Validate a nonempty readable helper-file path.
    '''bash
    validate_var_file "helper" "\${BASH_SOURCE[0]}" 0 "true"
    '''

  2. Reject an empty required file value.
    '''bash
    if ! validate_var_file "fil_in" ""; then
      echo "fil_in is required" >&2
    fi
    '''
EOM
    )

    if [[ "${var_nam}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${var_nam}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'var_nam', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${pth_fil}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'pth_fil', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    validate_var  "${var_nam}" "${pth_fil}" "${idx}" || return 1
    validate_file "${var_nam}" "${pth_fil}" "${idx}" "${chk_r}" || return 1
}


function validate_var_dir() {
    local var_nam="${1:-}"    # Variable name to use in error messages
    local pth_dir="${2:-}"    # Directory path to validate
    local idx="${3:-0}"       # Optional index (for arrays); defaults to '0'
    local chk_w="${4:-true}"  # If true, require directory to be writable
    local show_help           # Help message

    show_help=$(cat << EOM
Usage
-----
  validate_var_dir
    [--help] var_nam pth_dir [idx] [chk_w]

  Check that a variable value is set/non-empty and points to an existing directory. By default, also check that the directory is writable.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  var_nam : str
    Variable name to use in error messages.

  2  pth_dir : dir
    Directory path to validate.

  3  idx : int
    Optional array index for error messages (default: 0).

  4  chk_w : bool
    Whether to require writability: case-insensitive Boolean-like true|t|yes|y|1 or false|f|no|n|0 (default: true).

Returns
-------
  0 if the value is non-empty and the directory exists and, when 'chk_w=true', is writable; otherwise 1.

Notes
-----
  Runtime requirements:
    bash >= 4.4

  Accepted true-like values: true, t, yes, y, 1. Accepted false-like values: false, f, no, n, 0. Matching is case-insensitive; surrounding whitespace and empty values are invalid. Successful normalization emits 'true' or 'false'.

Examples
--------
  1. Validate a nonempty writable output directory.
    '''bash
    validate_var_dir "dir_out" "\${PWD}" 0 "true"
    '''

  2. Reject an empty required directory value.
    '''bash
    if ! validate_var_dir "dir_out" ""; then
      echo "dir_out is required" >&2
    fi
    '''
EOM
    )

    if [[ "${var_nam}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${var_nam}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'var_nam', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${pth_dir}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'pth_dir', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    validate_var "${var_nam}" "${pth_dir}" "${idx}" || return 1
    validate_dir "${var_nam}" "${pth_dir}" "${idx}" "${chk_w}" || return 1
}


#  Print one or more variable-assignment strings to stderr
function debug_var() { printf "%s\n\n" "$@" >&2; }


function check_arr_files() {
    local desc="${1:-}"  # Description of the file group being checked
    local file           # File path to validate
    local show_help      # Help message

    show_help=$(cat << EOM
Usage
-----
  check_arr_files
    [--help] desc file1 [file2 ...]

  Check that all supplied file paths exist as regular files.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1   desc : str
    Short description of the file group being validated; used in error messages.

  2   file1 : file
    First file path to validate.

  3+  file2 : file
    Additional file path(s) to validate.

Returns
-------
  0 if all supplied files exist; otherwise 1.

Notes
-----
  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Validate the sourced helper file as a one-file group.
    '''bash
    check_arr_files "helper" "\${BASH_SOURCE[0]}"
    '''

  2. Reject a group containing a missing file.
    '''bash
    if ! check_arr_files "input BAM" "\${BASH_SOURCE[0]}" "/missing.bam"; then
      echo "input BAM group is incomplete" >&2
    fi
    '''
EOM
    )

    if [[ "${desc}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${desc}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'desc', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    shift

    #  Check that files were supplied
    if (( $# < 1 )); then
        echo_err_func "${FUNCNAME[0]}" \
            "no files supplied to validate for '${desc}'."
        return 1
    fi

    #  Check each supplied file path
    for file in "$@"; do
        if [[ ! -f "${file}" ]]; then
            echo_err_func "${FUNCNAME[0]}" \
                "'${desc}' file does not exist as a regular file: '${file}'."
            return 1
        fi
    done

    return 0
}


function check_arr_lengths() {
    local arr_nam_1="${1:-}"  # Name of first indexed array
    local arr_nam_2="${2:-}"  # Name of second indexed array
    local arr_siz_1           # Size of first array
    local arr_siz_2           # Size of second array
    local decl
    local show_help           # Help message

    show_help=$(cat << EOM
Usage
-----
  check_arr_lengths
    [--help] arr_nam_1 arr_nam_2

  Check that two named indexed arrays exist and have matching lengths.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  arr_nam_1 : str
    Name of the first indexed array.

  2  arr_nam_2 : str
    Name of the second indexed array.

Returns
-------
  0 if both arrays exist, are indexed arrays, and have the same length; otherwise 1.

Notes
-----
  Runtime requirements:
    bash >= 4.4

  - Array names must be valid Bash variable names.
  - This helper uses namerefs and is intended for indexed arrays.

Examples
--------
  1. Accept two indexed arrays with matching lengths.
    '''bash
    arr_fil_in=( "a.bam" "b.bam" )
    arr_fil_out=( "a.bdg" "b.bdg" )
    check_arr_lengths "arr_fil_in" "arr_fil_out"
    '''

  2. Reject indexed arrays with different lengths.
    '''bash
    arr_fil_in=( "a.bam" "b.bam" )
    arr_fil_out=( "a.bdg" )
    if ! check_arr_lengths "arr_fil_in" "arr_fil_out"; then
      echo "input and output arrays differ in length" >&2
    fi
    '''
EOM
    )

    if [[ "${arr_nam_1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${arr_nam_1}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'arr_nam_1', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${arr_nam_2}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'arr_nam_2', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    if [[
        ! "${arr_nam_1}" =~ ^[a-zA-Z_][a-zA-Z0-9_]*$ \
        || ! "${arr_nam_2}" =~ ^[a-zA-Z_][a-zA-Z0-9_]*$
    ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "invalid array name '${arr_nam_1}' and/or '${arr_nam_2}'."
        return 1
    fi

    if ! decl="$(declare -p "${arr_nam_1}" 2> /dev/null)"; then
        echo_err_func "${FUNCNAME[0]}" \
            "'${arr_nam_1}' is unset."
        return 1
    elif [[ "${decl}" != declare\ -a* ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'${arr_nam_1}' is not an indexed array."
        return 1
    fi

    if ! decl="$(declare -p "${arr_nam_2}" 2> /dev/null)"; then
        echo_err_func "${FUNCNAME[0]}" \
            "'${arr_nam_2}' is unset."
        return 1
    elif [[ "${decl}" != declare\ -a* ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'${arr_nam_2}' is not an indexed array."
        return 1
    fi

    local -n arr_1="${arr_nam_1}"
    local -n arr_2="${arr_nam_2}"

    arr_siz_1="${#arr_1[@]}"
    arr_siz_2="${#arr_2[@]}"

    if (( arr_siz_1 != arr_siz_2 )); then
        echo_err_func "${FUNCNAME[0]}" \
            "array length mismatch. '${arr_nam_1}' has '${arr_siz_1}'" \
            "element(s), whereas '${arr_nam_2}' has '${arr_siz_2}'."
        return 1
    fi

    return 0
}


#TODO: retire this function once all call sites have been moved to the
#+     'validate_*' family
function check_file_dir_exists() {
    local type="${1:-}"   # Item type: 'f' for file or 'd' for directory
    local item="${2:-}"   # File or directory path to check
    local name="${3:-}"   # Optional argument name for error messages
    local typ_itm msg
    # local check_flag
    local show_help

    show_help=$(cat << EOM
Usage
-----
  check_file_dir_exists
    [--help] type item [name]

  Check the existence of a file or directory based on the provided type.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  type : str
    Type to check for existence: 'f' for file or 'd' for directory.

  2  item : path
    File or directory path to check.

  3  name : str
    Optional argument name to associate with the item in error messages.

Returns
-------
  0 if the file or directory exists; otherwise, 1 and an error message.

Notes
-----
  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Check an existing regular file.
    '''bash
    check_file_dir_exists "f" "\${BASH_SOURCE[0]}" "fil_in"
    '''

  2. Check an existing directory.
    '''bash
    check_file_dir_exists "d" "\${PWD}" "dir_out"
    '''
EOM
    )

    #  Parse and check function parameters
    if [[ "${type}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${type}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'type', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${item}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'item', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  Validate 'type'
    if [[ "${type}" != "f" && "${type}" != "d" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'type', is invalid: '${type}'. Expected" \
            "'f' for file or 'd' for directory."
        return 1
    fi

    #  Check file or directory existence; construct and return error message if
    #+ applicable
    if [[ "${type}" == "f" ]]; then
        typ_itm="file"
        msg="file does not exist"
    else
        typ_itm="directory"
        msg="directory does not exist"
    fi

    if [[ -n "${name}" ]]; then
        msg="${typ_itm} associated with '--${name}' does not exist"
    fi

    case "${type}" in
        f)
            if [[ ! -f "${item}" ]]; then
                echo_err_func "${FUNCNAME[0]}" \
                    "${msg}: '${item}'."
                return 1
            fi
            ;;
        d)
            if [[ ! -d "${item}" ]]; then
                echo_err_func "${FUNCNAME[0]}" \
                    "${msg}: '${item}'."
                return 1
            fi
            ;;
    esac
}


function debug_arr_contents() {
    local arr_nam    # Name of array variable to inspect
    local decl
    local show_help  # Help message

    show_help=$(cat << EOM
Usage
-----
  debug_arr_contents
    [--help] arr1 [arr2 ...]

  Print the contents of one or more named indexed arrays to stderr.

  For each supplied name, this helper
    - skips invalid variable names,
    - skips unset variables,
    - skips variables that are not indexed arrays, and
    - prints non-empty indexed arrays in a compact debug form.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1   arr1 : str
    Name of the first indexed array to inspect.

  2+  arr2 : str
    Name(s) of additional indexed arrays to inspect.

Returns
-------
  0; invalid, unset, or non-indexed variables are skipped silently.

Notes
-----
  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Print one nonempty indexed array.
    '''bash
    arr_fil_in=( "a.bam" "b.bam" )
    debug_arr_contents "arr_fil_in"
    '''

  2. Inspect several arrays while silently skipping an unset name.
    '''bash
    arr_fil_in=( "a.bam" )
    arr_fil_out=( "a.bdg" )
    debug_arr_contents "arr_fil_in" "arr_fil_out" "arr_unset"
    '''
EOM
    )

    if [[ "${1:-}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${1:-}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "at least one array name must be supplied."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    for arr_nam in "$@"; do
        if ! [[ "${arr_nam}" =~ ^[a-zA-Z_][a-zA-Z0-9_]*$ ]]; then
            continue
        fi

        if ! decl="$(declare -p "${arr_nam}" 2> /dev/null)"; then
            continue
        elif [[ "${decl}" != declare\ -a* ]]; then
            continue
        fi

        local -n arr_ref="${arr_nam}"
        if (( ${#arr_ref[@]} > 0 )); then
            echo "  - ${arr_nam}=( ${arr_ref[*]} )" >&2
        fi
    done
}


function check_arr_nonempty() {
    local arr_nam="${1:-}"
    local src_nam="${2:-}"
    local arr_len=0
    local decl
    local show_help

    show_help=$(cat << EOM
Usage
-----
  check_arr_nonempty
    [--help] arr_nam src_nam

  Check that reconstructed array 'arr_nam' has at least one element after being built from serialized source variable 'src_nam'.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  arr_nam : str
    Name of the reconstructed array variable to check.

  2  src_nam : str
    Name of the serialized source variable from which the array was reconstructed.

Returns
-------
  0 if the reconstructed array is non-empty; 1 if argument parsing fails or if the array is empty.

Notes
-----
  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Accept a reconstructed array with one element.
    '''bash
    arr_fil_in=( "sample.bam" )
    check_arr_nonempty "arr_fil_in" "csv_fil_in"
    '''

  2. Reject an empty reconstructed array.
    '''bash
    arr_fil_in=()
    if ! check_arr_nonempty "arr_fil_in" "csv_fil_in"; then
      echo "csv_fil_in produced no entries" >&2
    fi
    '''
EOM
    )

    if [[ "${arr_nam}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${arr_nam}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'arr_nam', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${src_nam}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'src_nam', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    if ! [[ "${arr_nam}" =~ ^[a-zA-Z_][a-zA-Z0-9_]*$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "invalid array name '${arr_nam}'."
        return 1
    fi

    if ! decl="$(declare -p "${arr_nam}" 2> /dev/null)"; then
        echo_err_func "${FUNCNAME[0]}" \
            "reconstructed array '${arr_nam}' is unset. Check input for" \
            "serialized variable '${src_nam}'."
        return 1
    elif [[ "${decl}" != declare\ -a* ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'${arr_nam}' is not an indexed array."
        return 1
    fi

    local -n arr_ref="${arr_nam}"
    arr_len="${#arr_ref[@]}"

    if (( arr_len == 0 )); then
        echo_err_func "${FUNCNAME[0]}" \
            "reconstructed array '${arr_nam}' is empty. Check input for" \
            "serialized variable '${src_nam}'."
        return 1
    fi

    return 0
}


function check_arr_len_bcst() {
    local n_req="${1:-}"  # Required full length for per-sample arrays
    local arr_nam         # Name of array being checked
    local arr_len         # Length of current array
    local decl
    local show_help       # Help message

    show_help=$(cat << EOM
Usage
-----
  check_arr_len_bcst
    [--help] n_req arr1 [arr2 ...]

  Check that each named array has length 0, 1, or 'n_req'.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1   n_req : int
    Required full length for per-sample arrays.

  2   arr1 : str
    Name of the first array variable to check.

  3+  arr2 : str
    Name(s) of additional array variable(s) to check.

Returns
-------
  0 if every named array has length 0, 1, or 'n_req'; 1 if argument parsing fails or if any array has an invalid length.

Notes
-----
  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Accept a single broadcast value for three samples.
    '''bash
    arr_scl_fct=( "1" )
    check_arr_len_bcst 3 "arr_scl_fct"
    '''

  2. Reject a two-element array for three samples.
    '''bash
    arr_scl_fct=( "1" "2" )
    if ! check_arr_len_bcst 3 "arr_scl_fct"; then
      echo "scaling-factor vector cannot broadcast" >&2
    fi
    '''
EOM
    )

    if [[ "${n_req}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${n_req}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'n_req', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ $# -lt 2 ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "at least one array name must be supplied after 'n_req'."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    if ! [[ "${n_req}" =~ ^[0-9]+$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'n_req', must be a non-negative integer."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    shift

    for arr_nam in "$@"; do
        if ! [[ "${arr_nam}" =~ ^[a-zA-Z_][a-zA-Z0-9_]*$ ]]; then
            echo_err_func "${FUNCNAME[0]}" \
                "invalid array name '${arr_nam}'."
            return 1
        fi

        if ! decl="$(declare -p "${arr_nam}" 2> /dev/null)"; then
            echo_err_func "${FUNCNAME[0]}" \
                "array '${arr_nam}' is unset."
            return 1
        elif [[ "${decl}" != declare\ -a* ]]; then
            echo_err_func "${FUNCNAME[0]}" \
                "'${arr_nam}' is not an indexed array."
            return 1
        fi

        local -n arr_ref="${arr_nam}"

        if [[ -v "arr_ref[@]" ]]; then
            arr_len="${#arr_ref[@]}"
        else
            arr_len=0
        fi

        if (( arr_len != 0 && arr_len != 1 && arr_len != n_req )); then
            echo_err_func "${FUNCNAME[0]}" \
                "array '${arr_nam}' has length '${arr_len}', but expected 0," \
                "1, or ${n_req}."
            return 1
        fi
    done

    return 0
}


#  Print an error message when function script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    err_source_only "${BASH_SOURCE[0]}"
fi
