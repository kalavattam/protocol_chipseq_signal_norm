#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: check_numbers.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
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
# shellcheck disable=SC1091
{
    _dir_src_nos="$(
        cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
    )"

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
Usage:
  check_flt_nonneg [-h|--hlp|--help] val [nam]

Description:
  Check that a value is a non-negative integer or float.

Positional arguments:
  1  val  <flt>  The value to check.
  2  nam  <str>  Name of argument or variable associated with the value (optional).

Returns:
  0 if the check passes; otherwise 1.

Examples:
  '''bash
  check_flt_nonneg   0 "scale"  # Returns 0
  check_flt_nonneg 4.5 "scale"  # Returns 0
  check_flt_nonneg  -1 "scale"  # Returns 1 and prints an error message to stderr
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
Usage:
  check_flt_pos [-h|--hlp|--help] val [nam]

Description:
  Check that a value is a positive integer or float.

Positional arguments:
  1  val  <flt>  The value to check.
  2  nam  <str>  Name of argument or variable associated with the value (optional).

Returns:
  0 if the check passes; 1 if not.

Examples:
  Check that the value assigned to the argument '--scaling_factor' is a positive float:
  '''bash
  check_flt_pos 4.5 "scaling_factor"  # Returns 0
  '''

  '''bash
  check_flt_pos a "scaling_factor"  # Returns 1 and prints an error message to stderr
  '''

  '''bash
  check_flt_pos -2.5  # Returns 1 and prints an error message to stderr
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
Usage:
  check_format_time [-h|--hlp|--help] time

Description:
  Checks that a string is formatted as 'mm:ss', 'h:mm:ss', or 'hh:mm:ss', where minutes and seconds must be between 00 and 59.

Positional argument:
  1  time  <str>  The time string to check.

Returns:
  0 if the time string is correctly formatted; otherwise 1.

Examples:
  '''bash
  check_format_time "2:30:15"    # Returns 0
  check_format_time "2:61:44"    # Shows error message and returns 1
  check_format_time "25:165:80"  # Shows error message and returns 1
  check_format_time --help       # Shows help message and returns 0
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
Usage:
  check_int_nonneg [-h|--hlp|--help] val [nam]

Description:
  Check that a value is an integer greater than or equal (gte) to 0, i.e., a non-negative integer.

Positional arguments:
  1  val  <int>  The value to check.
  2  nam  <str>  Name of argument or variable associated with the value (optional).

Returns:
  0 if the check passes, otherwise returns an error message and exit code 1.

Examples:
  '''bash
  check_int_nonneg 4 "threads"  # Returns 0
  '''

  '''bash
  check_int_nonneg a "threads"  # Returns 1 and prints an error message to stderr
  '''

  '''bash
  check_int_nonneg -2  # Returns 1 and prints an error message to stderr
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
Usage:
  check_int_pos [-h|--hlp|--help] val [nam]

Description:
  Check that a value is an integer greater than or equal to (gte) 1.

Positional arguments:
  1  val  <int>  The value to check.
  2  nam  <str>  Name of argument or variable associated with the value (optional).

Returns:
  0 if the check passes, otherwise returns an error message and exit code 1.

Examples:
  '''bash
  check_int_pos 4 "threads"  # Returns 0
  '''

  '''bash
  check_int_pos a "threads"  # Returns 1 and prints an error message to stderr
  '''

  '''bash
  check_int_pos -2  # Returns 1 and prints an error message to stderr
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
Usage:
  check_arr_int_pos [-h|--hlp|--help] arr_nam src_nam

Description:
  Check that every element in named indexed array 'arr_nam' is a positive integer.

Positional arguments:
  1  arr_nam  <str>  Name of indexed array to validate.
  2  src_nam  <str>  Source argument or variable name for error messages.

Returns:
  0 if all array elements are positive integers; otherwise 1.
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
Usage:
  check_arr_num_pos [-h|--hlp|--help] arr_nam src_nam

Description:
  Check that every element in named indexed array 'arr_nam' is a positive integer or float.

Positional arguments:
  1  arr_nam  <str>  Name of indexed array to validate.
  2  src_nam  <str>  Source argument or variable name for error messages.

Returns:
  0 if all array elements are positive numbers; otherwise 1.
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
    local scl_fct="${1:-}"  # Comma-separated scaling factors to validate
    local entries=()        # Array for individual comma-separated components
    local entry             # Individual scaling factor: 'num:den' or 'num'
    local num den           # Values for validation and formatting
    local val_fmt=()        # Array for validated and formatted entries
    local out_str           # Final formatted, comma-separated string output
    local show_help         # Help message/function documentation

    show_help=$(cat << EOM
Usage:
  check_scl_fct [-h|--hlp|--help] scl_fct

Description:
  Validate and format positive floats for precomputed scaling factors used with deepTools 'bamCompare' or 'bigwigCompare'.

  This helper ensures the following:
    - scaling factors are provided as a single comma-separated string,
    - each factor pair (numerator and denominator) is in the format 'num:den',
    - single values (numerators) are automatically paired with a denominator of 1,
    - all values are validated as positive floats, and
    - leading decimal points (e.g., ".5") are zero-padded (e.g., "0.5").

Positional argument:
  1  scl_fct  <str>  Comma-separated scaling factors to validate.

Returns:
  0 and a validated, comma-separated string of formatted scaling factors; example: '0.4:1,1.2:1,0.65:0.8'. Otherwise, returns 1 with an error message.

Examples:
  '''bash
  scl_fct=0.5,1.2:3,.65
  check_scl_fct \${scl_fct}  # Returns 0.5:1,1.2:3,0.65:1
  '''

  '''bash
  scl_fct=0.00456,0.00789
  check_scl_fct \${scl_fct}  # Returns 0.00456:1,0.00789:1
  '''

  '''bash
  scl_fct="num:den"
  check_scl_fct \${scl_fct}  # Returns 1 and prints an error message to stderr
  '''

  '''bash
  scl_fct="0.33:den"
  check_scl_fct \${scl_fct}  # Returns 1 and prints an error message to stderr
  '''

  '''bash
  scl_fct=0.5,1.2a,.65
  check_scl_fct \${scl_fct}  # Returns 1 and prints an error message to stderr
  '''

  '''bash
  scl_fct=hotdog
  check_scl_fct \${scl_fct}  # Returns 1 and prints an error message to stderr
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
