#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: check_numbers.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
# GPT-5.6) were used in design, development, and documentation, with all output
# reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


# check_flt_nonneg
# check_flt_pos
# check_format_time
# check_int_nonneg
# check_int_pos
# check_arr_int_pos
# check_arr_num_pos
# check_scl_fct


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
    _dir_src_nos="$(
        cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
    )"

    # shellcheck source=lib/bash/core/source_helpers.sh
    source "${_dir_src_nos}/source_helpers.sh" || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source '${_dir_src_nos}/source_helpers.sh'." >&2

        if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
            return 1
        else
            exit 1
        fi
    }

    source_helpers "${_dir_src_nos}" \
        check_source format_outputs || {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to source required helper dependencies." >&2

            if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
                return 1
            else
                exit 1
            fi
        }

    unset _dir_src_nos
}


function check_flt_nonneg() {
    local val="${1:-}"
    local nam="${2:-}"
    local show_help

    show_help=$(cat << EOM
Usage
-----
  check_flt_nonneg
    [--help] val [nam]

  Check that a value is a non-negative integer or float.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  val : float
    The value to check.

  2  nam : str
    Name of argument or variable associated with the value. If empty, a generic label is used.

Returns
-------
  0 if the check passes; otherwise 1.

Notes
-----
  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. After sourcing this file, accept zero as non-negative.
    '''bash
    check_flt_nonneg 0 "scale"
    '''

  2. Reject a negative floating-point value.
    '''bash
    if ! check_flt_nonneg -0.5 "scale"; then
      echo "scale must be non-negative" >&2
    fi
    '''
EOM
    )

    #  Parse and check function arguments
    if [[ "${val}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${val}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'val', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    if [[ "${val}" =~ ^[+]?([0-9]+([.][0-9]*)?|[.]?[0-9]+)$ ]]; then
        return 0
    else
        if [[ -n "${nam}" ]]; then
            echo_err_func "${FUNCNAME[0]}" \
                "'--${nam}' was assigned '${val}' but must be a" \
                "non-negative integer or float."
        else
            echo_err_func "${FUNCNAME[0]}" \
                "'${val}' is not a non-negative integer or float."
        fi
        return 1
    fi
}


function check_flt_pos() {
    local val="${1:-}"
    local nam="${2:-}"
    local show_help

    show_help=$(cat << EOM
Usage
-----
  check_flt_pos
    [--help] val [nam]

  Check that a value is a positive integer or float.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  val : float
    The value to check.

  2  nam : str
    Name of argument or variable associated with the value. If empty, a generic label is used.

Returns
-------
  0 if the check passes; 1 if not.

Notes
-----
  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Accept a positive floating-point scaling factor.
    '''bash
    check_flt_pos 4.5 "scaling_factor"
    '''

  2. Reject zero because it is not positive.
    '''bash
    if ! check_flt_pos 0 "scaling_factor"; then
      echo "scaling_factor must be positive" >&2
    fi
    '''
EOM
    )

    #  Parse and check function arguments
    if [[ "${val}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${val}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'val', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  Check that the value is numeric and strictly greater than 0
    if \
        [[ "${val}" =~ ^[+]?([0-9]+([.][0-9]*)?|[.]?[0-9]+)$ ]] \
        && ! [[ "${val}" =~ ^[+]?(0+([.]0*)?|[.]0+)$ ]]
    then
        return 0
    else
        if [[ -n "${nam}" ]]; then
            echo_err_func "${FUNCNAME[0]}" \
                "'--${nam}' was assigned '${val}' but must be a" \
                "positive integer or float."
        else
            echo_err_func "${FUNCNAME[0]}" \
                "'${val}' is not a positive integer or float."
        fi
        return 1
    fi
}


function check_format_time() {
    local time="${1:-}"
    local show_help

    show_help=$(cat << EOM
Usage
-----
  check_format_time
    [--help] time

  Checks that a string is formatted as 'mm:ss', 'h:mm:ss', or 'hh:mm:ss', where minutes and seconds must be between 00 and 59.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  time : str
    Slurm job time limit.

Returns
-------
  0 if the time string is correctly formatted; otherwise 1.

Notes
-----
  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Accept an hours-minutes-seconds Slurm time limit.
    '''bash
    check_format_time "2:30:15"
    '''

  2. Reject a time with minutes outside 00 through 59.
    '''bash
    if ! check_format_time "2:61:44"; then
      echo "time is malformed" >&2
    fi
    '''
EOM
    )

    #  Parse and check function arguments
    if [[ "${time}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${time}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'time', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  Check for 'mm:ss', 'h:mm:ss', or 'hh:mm:ss'
    if [[ ! "${time}" =~ ^([0-9]{1,2}:)?[0-5][0-9]:[0-5][0-9]$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'${time}' is not a valid time format. Expected format is" \
            "'mm:ss', 'h:mm:ss', or 'hh:mm:ss'."
        return 1
    fi
}


function check_int_nonneg() {
    local val="${1:-}"
    local nam="${2:-}"
    local show_help

    show_help=$(cat << EOM
Usage
-----
  check_int_nonneg
    [--help] val [nam]

  Check that a value is an integer greater than or equal (gte) to 0, i.e., a non-negative integer.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  val : int
    The value to check.

  2  nam : str
    Name of argument or variable associated with the value. If empty, a generic label is used.

Returns
-------
  0 if the check passes, otherwise returns an error message and exit code 1.

Notes
-----
  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Accept zero as a non-negative integer.
    '''bash
    check_int_nonneg 0 "max_job"
    '''

  2. Reject a negative integer.
    '''bash
    if ! check_int_nonneg -2 "max_job"; then
      echo "max_job must be non-negative" >&2
    fi
    '''
EOM
    )

    #  Parse and check function arguments, printing help message as appropriate
    if [[ "${val}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${val}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'val', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  Perform the check and return an error message if it fails
    if ! [[ "${val}" =~ ^[0-9]+$ ]]; then
        if [[ -n "${nam}" ]]; then
            echo_err_func "${FUNCNAME[0]}" \
                "'--${nam}' was assigned '${val}' but must be a" \
                "non-negative integer."
        else
            echo_err_func "${FUNCNAME[0]}" \
                "'${val}' is not a non-negative integer."
        fi
        return 1
    fi
}


function check_int_pos() {
    local val="${1:-}"
    local nam="${2:-}"
    local show_help

    show_help=$(cat << EOM
Usage
-----
  check_int_pos
    [--help] val [nam]

  Check that a value is an integer greater than or equal to (gte) 1.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  val : int
    The value to check.

  2  nam : str
    Name of argument or variable associated with the value. If empty, a generic label is used.

Returns
-------
  0 if the check passes, otherwise returns an error message and exit code 1.

Notes
-----
  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Accept a positive thread count.
    '''bash
    check_int_pos 4 "threads"
    '''

  2. Reject zero as a thread count.
    '''bash
    if ! check_int_pos 0 "threads"; then
      echo "threads must be positive" >&2
    fi
    '''
EOM
    )

    #  Parse and check function arguments, printing help message as appropriate
    if [[ "${val}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${val}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'val', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  Perform the check and return an error message if it fails
    if ! [[ "${val}" =~ ^[1-9][0-9]*$ ]]; then
        if [[ -n "${nam}" ]]; then
            echo_err_func "${FUNCNAME[0]}" \
                "'--${nam}' was assigned '${val}' but must be an integer" \
                "greater than or equal to 1."
        else
            echo_err_func "${FUNCNAME[0]}" \
                "'${val}' is not an integer greater than or equal to 1."
        fi
        return 1
    fi
}


function check_arr_int_pos() {
    local arr_nam="${1:-}"
    local src_nam="${2:-}"
    local idx val decl
    local show_help

    show_help=$(cat << EOM
Usage
-----
  check_arr_int_pos
    [--help] arr_nam src_nam

  Check that every element in named indexed array 'arr_nam' is a positive integer.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  arr_nam : str
    Name of indexed array to validate.

  2  src_nam : str
    Source argument or variable name for error messages.

Returns
-------
  0 if all array elements are positive integers; otherwise 1.

Notes
-----
  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Accept a vector of positive integer fragment lengths.
    '''bash
    arr_usr_frg=( 150 200 )
    check_arr_int_pos "arr_usr_frg" "csv_usr_frg"
    '''

  2. Reject a vector containing zero.
    '''bash
    arr_usr_frg=( 150 0 )
    if ! check_arr_int_pos "arr_usr_frg" "csv_usr_frg"; then
      echo "csv_usr_frg contains a non-positive value" >&2
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

    for idx in "${!arr_ref[@]}"; do
        val="${arr_ref[${idx}]}"

        if ! check_int_pos "${val}" "${src_nam}"; then
            echo_err_func "${FUNCNAME[0]}" \
                "'${src_nam}' element at index '${idx}' failed positive" \
                "integer validation."
            return 1
        fi
    done
}


function check_arr_num_pos() {
    local arr_nam="${1:-}"
    local src_nam="${2:-}"
    local idx val decl
    local show_help

    show_help=$(cat << EOM
Usage
-----
  check_arr_num_pos
    [--help] arr_nam src_nam

  Check that every element in named indexed array 'arr_nam' is a positive integer or float.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  arr_nam : str
    Name of indexed array to validate.

  2  src_nam : str
    Source argument or variable name for error messages.

Returns
-------
  0 if all array elements are positive numbers; otherwise 1.

Notes
-----
  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Accept a vector of positive numeric scaling factors.
    '''bash
    arr_scl_fct=( 0.5 2 )
    check_arr_num_pos "arr_scl_fct" "csv_scl_fct"
    '''

  2. Reject a vector containing a negative value.
    '''bash
    arr_scl_fct=( 0.5 -2 )
    if ! check_arr_num_pos "arr_scl_fct" "csv_scl_fct"; then
      echo "csv_scl_fct contains a non-positive value" >&2
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

    for idx in "${!arr_ref[@]}"; do
        val="${arr_ref[${idx}]}"

        if ! check_flt_pos "${val}" "${src_nam}"; then
            echo_err_func "${FUNCNAME[0]}" \
                "'${src_nam}' element at index '${idx}' failed positive" \
                "number validation."
            return 1
        fi
    done
}


#TODO: audit current usage; keep this helper even if unused
function check_scl_fct() {
    local scl_fct="${1:-}"  # Scaling factor. Comma-separated scaling factors to validate
    local entries=()        # Array for individual comma-separated components
    local entry             # Individual scaling factor: 'num:den' or 'num'
    local num den           # Values for validation and formatting
    local val_fmt=()        # Array for validated and formatted entries
    local out_str           # Final formatted, comma-separated string output
    local show_help         # Help message/function documentation

    show_help=$(cat << EOM
Usage
-----
  check_scl_fct
    [--help] scl_fct

  Validate and format positive floats for precomputed scaling factors used with deepTools 'bamCompare' or 'bigwigCompare'.

  This helper ensures the following:
    - scaling factors are provided as a single comma-separated string,
    - each factor pair (numerator and denominator) is in the format 'num:den',
    - single values (numerators) are automatically paired with a denominator of 1,
    - all values are validated as positive floats, and
    - leading decimal points (e.g., ".5") are zero-padded (e.g., "0.5").

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  scl_fct : str
    Scaling factor. Comma-separated scaling factors to validate.

Returns
-------
  0 and a validated, comma-separated string of formatted scaling factors; example: '0.4:1,1.2:1,0.65:0.8'. Otherwise, returns 1 with an error message.

Notes
-----
  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Format numerator-only and numerator-denominator factors.
    '''bash
    check_scl_fct "0.5,1.2:3,.65"
    '''

  2. Reject a nonnumeric denominator.
    '''bash
    if ! check_scl_fct "0.33:den"; then
      echo "scaling-factor specification is invalid" >&2
    fi
    '''
EOM
    )

    #  Parse and check function parameter
    if [[ "${scl_fct}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${scl_fct}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'scl_fct', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  Split the input string by commas
    IFS=',' read -ra entries <<< "${scl_fct}"

    for entry in "${entries[@]}"; do
        #  Check that entry is a 'num:den' pair (containing a colon)
        if [[ "${entry}" == *:* ]]; then
            #  Split into num and den
            IFS=':' read -r num den <<< "${entry}"

            #  Validate that both 'num' and 'den' are positive floats
            if ! check_flt_pos "${num}" >/dev/null 2>&1; then
                echo_err_func "${FUNCNAME[0]}" \
                    "invalid 'num' in '${entry}' (must be a positive float)."
                return 1
            fi

            if ! check_flt_pos "${den}" >/dev/null 2>&1; then
                echo_err_func "${FUNCNAME[0]}" \
                    "invalid 'den' in '${entry}' (must be a positive float)."
                return 1
            fi

            #  0-pad values starting with a decimal point
            if [[ "${num}" == .* ]]; then num="0${num}"; fi
            if [[ "${den}" == .* ]]; then den="0${den}"; fi

            #  Add back to the formatted values as 'num:den'
            val_fmt+=( "${num}:${den}" )
        else
            #  Validate single numerator-only entry as a positive float
            if ! check_flt_pos "${entry}" >/dev/null 2>&1; then
                echo_err_func "${FUNCNAME[0]}" \
                    "invalid numerator-only scaling factor '${entry}' (must" \
                    "be a positive float)."
                return 1
            fi

            #  Zero-pad values starting with a decimal point
            if [[ "${entry}" == .* ]]; then entry="0${entry}"; fi

            #  Append as 'entry:1' (numerator-only case)
            val_fmt+=( "${entry}:1" )
        fi
    done

    #  Join the validated, formatted values back into a comma-separated string
    out_str=$(IFS=','; echo "${val_fmt[*]}")

    #  Output the final formatted string
    echo "${out_str}"
}


#  Print an error message when function script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    err_source_only "${BASH_SOURCE[0]}"
fi
