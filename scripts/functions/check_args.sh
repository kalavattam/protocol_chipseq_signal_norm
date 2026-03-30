#!/bin/bash
# -*- coding: utf-8 -*-
# 
# Script: check_args.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
# 
# Distributed under the MIT license.


function check_arg_supplied() {
    local asm=""
    local nam=""
    local show_help

    show_help=$(cat << EOM
Usage:
  check_arg_supplied [<empty>|-h|--hlp|--help] --asm <str> --nam <str>

Description:
  Checks that a given argument has been assigned a value.

Keyword arguments:
  -a, --asm  <str>  The value (i.e., assignment) of the argument to check.
  -n, --nam  <str>  The name of the argument to check.

Returns:
  0 if the argument has been assigned a value; otherwise, returns 1 if the argument is empty or not set.

Example:
  Check that a required argument ('req_arg') is supplied:
  '''bash
  req_arg="some_value"
  
  check_arg_supplied -a "\${req_arg}" -n "req_arg"
  '''
EOM
    )

    #  Parse and check function arguments
    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
            -a|--asm)
                if [[ $# -lt 2 ]]; then
                    echo "Error: '--asm' requires a value." >&2
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                fi
                asm="${2}"
                shift 2
                ;;

            -n|--nam)
                if [[ $# -lt 2 ]]; then
                    echo "Error: '--nam' requires a value." >&2
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                fi
                nam="${2}"
                shift 2
                ;;

            *)
                echo "## Unknown argument passed: '${1}' ##" >&2
                echo >&2
                echo "${show_help}" >&2
                return 1
                ;;
        esac
    done

    #  Check that required arguments are supplied
    if [[ -z "${nam}" ]]; then
        echo "Error: '--nam' is required." >&2
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  Check that the argument has been supplied; return an informative error
    #+ message if not
    if [[ -z "${asm}" ]]; then
        echo "Error: '--${nam}' is required." >&2
        return 1
    fi
}


function check_args_mut_excl() {
    local nam_1="${1:-}"
    local val_1="${2:-}"
    local nam_2="${3:-}"
    local val_2="${4:-}"
    local show_help

    show_help=$(cat << EOM
Usage:
  check_args_mut_excl [<empty>|-h|--hlp|--help] nam_1 val_1 nam_2 val_2

Description:
  Checks that two mutually exclusive arguments are not specified at the same time. If both are specified or neither is specified, an error is returned.

Positional arguments:
  1  nam_1  <str>  Name of the first argument (e.g., "norm").
  2  val_1  <str>  Value of the first argument.
  3  nam_2  <str>  Name of the second argument (e.g., "raw").
  4  val_2  <str>  Value of the second argument.

Returns:
  0 if the arguments are mutually exclusive and valid, 1 and an error message if both arguments are specified or neither is specified.

Example:
  '''bash
  check_args_mut_excl
      "arg_1_name" "\${arg_1_value}"
      "arg_2_name" "\${arg_2_value}"
  '''
EOM
    )

    #  Print help message if no arguments are passed or if help is requested
    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    #  Validate mutually exclusive arguments
    if [[ -n "${val_1}" && -n "${val_2}" ]]; then
        echo \
            "Error: Only one of '--${nam_1}' or '--${nam_2}' can be" \
            "specified at a time." >&2
        return 1
    elif [[ -z "${val_1}" && -z "${val_2}" ]]; then
        echo \
            "Error: One of '--${nam_1}' or '--${nam_2}' must be specified." >&2
        return 1
    fi
}


function check_flags_mut_excl() {
    local flg_1="${1:-}"
    local nam_1="${2:-}"
    local flg_2="${3:-}"
    local nam_2="${4:-}"
    local show_help

    show_help=$(cat << EOM
Usage:
  check_flags_mut_excl [<empty>|-h|--hlp|--help] nam_1 val_1 nam_2 val_2

Description:
  Checks that two mutually exclusive flags are not specified at the same time. If both are specified or neither is specified, an error is returned.

Positional parameters:
  1  flg_1  <bol>  The first flag to check.
  2  nam_1  <str>  The name of the first flag.
  3  flg_2  <bol>  The second flag to check.
  4  nam_2  <str>  The name of the second flag.

Returns:
  0 if the arguments are mutually exclusive and valid, 1 and an error message if both arguments are specified or neither is specified.

Examples:
  '''bash
  #  Check that either '--flg_1' or '--flg_2' is specified, but not both.
  flg_1=true
  flg_2=false
  check_flags_mut_excl "\${flg_1}" "flg_1" "\${flg_2}" "flg_2"  # Returns 0

  #  What happens if '--flg_1' are '--flg_2' are both specified?
  flg_1=true
  flg_2=true
  check_flags_mut_excl "\${flg_1}" "flg_1" "\${flg_2}" "flg_2"
  # Error: Only one of '--flg_1' or '--flg_2' can be specified at a time.
  '''
EOM
    )

    #  Parse and check function parameters
    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    #  Perform the check for two-flag mutual exclusivity
    if [[ "${flg_1}" == "true" && "${flg_2}" == "true" ]]; then
        echo \
            "Error: Only one of '--${nam_1}' or '--${nam_2}' can be" \
            "specified at a time." >&2
        return 1
    fi
}


#  Function to validate matching values between file pairs, i.e.,
#+ first file/numerator and second file/denominator
function check_match() {
    local var_dsc="${1:-}"  # Description of variable being compared
    local val_num="${2:-}"  # Numerator value
    local val_den="${3:-}"  # Denominator value
    local out_nam="${4:-}"  # Name of output variable to assign if matched

    if [[ "${val_num}" == "${val_den}" ]]; then
        #  Assign the matched value to the output variable if specified
        if [[ -n "${out_nam}" ]]; then eval "${out_nam}=\"${val_num}\""; fi
    else
        echo \
            "Error: ${var_dsc} does not match between numerator" \
            "(${val_num}) and denominator (${val_den})." >&2
        return 1
    fi
}


function check_str_delim() {
    local nam="${1:-}"
    local val="${2:-}"
    local show_help

    show_help=$(cat << EOM
Usage:
  check_str_delim [-h|--hlp|--help] nam val

Description:
  Checks that a string value does not contain improperly formatted comma or semicolon delimiters, including consecutive, leading, or trailing delimiters, mixed delimiters (e.g., ",;" or ";,"), or spaces before/after delimiters. Additionally, ensures the string is not empty.

Positional arguments:
  1  nam  <str>  Name or description of string to be checked (e.g., "infiles", "scl_fct").
  2  val  <str>  String value to validate.

Returns:
  0 if the string passes all checks; 1 and an error message if the string fails any of the checks.

Notes:
  - Validates no consecutive, leading, or trailing commas (e.g., ',,', ',', etc. at the start/end) or semicolons (e.g., ';;', ';', etc. at the start/end).
  - Validates no mixed adjacent delimiters (e.g., ',;' or ';,').
  - Validates no spaces before or after commas or semicolons (e.g., " , ", " ;", etc.).
  - Checks that the string is not empty.

Dependencies:
  - Bash

Examples:
  '''bash
  check_str_delim "infiles" "file1.bam,file2.bam;file3.bam"  # Returns 0

  check_str_delim "infiles" ""  # Returns 1 and the below message
  # Error: Improperly formatted infiles value: '' (empty value).

  check_str_delim "infiles" "file1.bam,,file2.bam"  # Returns 1 and the below message
  # Error: Improperly formatted infiles value: 'file1.bam,,file2.bam' (invalid commas).

  check_str_delim "scl_fct" "0.5,;1.0"  # Returns 1 and the below message
  # Error: Improperly formatted scl_fct value: '0.5,;1.0' (mixed invalid delimiters).

  check_str_delim "usr_frg" "300, 200;400"  # Returns 1 and the below message
  # Error: Improperly formatted usr_frg value: '300, 200;400' (spaces around commas).
  '''
EOM
    )

    #  Parse and check function arguments
    if [[ "${nam}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif (( $# < 1 )); then
        echo "Error: Positional argument 1, 'nam', is missing." >&2
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif (( $# < 2 )); then
        echo "Error: Positional argument 2, 'val', is missing." >&2
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  Check for empty value
    if [[ -z "${val}" ]]; then
        echo \
            "Error: Improperly formatted '${nam}' value: '${val}' (empty" \
            "value)." >&2
        return 1
    fi

    #  Check for improper commas
    if [[ "${val}" =~ ,{2,} || "${val}" =~ ^, || "${val}" =~ ,$ ]]; then
        echo \
            "Error: Improperly formatted '${nam}' value: '${val}' (invalid" \
            "commas)." >&2
        return 1
    fi

    #  Check for improper semicolons
    if [[ "${val}" =~ \;{2,} || "${val}" =~ ^\; || "${val}" =~ \;$ ]]; then
        echo \
            "Error: Improperly formatted '${nam}' value: '${val}' (invalid" \
            "semicolons)." >&2
        return 1
    fi

    #  Check for mixed improper delimiters
    if [[ "${val}" =~ ,\; || "${val}" =~ \;, ]]; then
        echo \
            "Error: Improperly formatted '${nam}' value: '${val}' (mixed" \
            "invalid delimiters)." >&2
        return 1
    fi

    #  Check for spaces before/after commas
    if [[
           "${val}" =~ [[:space:]]+,[[:space:]]* \
        || "${val}" =~ [[:space:]]*,[[:space:]]+
    ]]; then
        echo \
            "Error: Improperly formatted '${nam}' value: '${val}' (spaces" \
            "around commas)." >&2
        return 1
    fi

    #  Check for spaces before/after semicolons
    if [[
           "${val}" =~ [[:space:]]+\;[[:space:]]* \
        || "${val}" =~ [[:space:]]*\;[[:space:]]+
    ]]; then
        echo \
            "Error: Improperly formatted '${nam}' value: '${val}' (spaces" \
            "around semicolons)." >&2
        return 1
    fi
}
