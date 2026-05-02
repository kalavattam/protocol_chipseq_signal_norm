#!/usr/bin/env bash
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


# require_optarg
# check_arg_supplied
# check_args_mut_excl
# check_flags_mut_excl
# check_match
# check_str_delim


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
    _dir_src_args="$(
        cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
    )"

    source "${_dir_src_args}/source_helpers.sh" || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source '${_dir_src_args}/source_helpers.sh'." >&2

        if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
            return 1
        else
            exit 1
        fi
    }

    source_helpers "${_dir_src_args}" \
        check_source format_outputs || {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to source required helper dependencies." >&2

            if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
                return 1
            else
                exit 1
            fi
        }

    unset _dir_src_args
}


#  Require that an option be followed by a non-empty value that is not another
#+ option-like token
function require_optarg() {
    local opt="${1:-}"
    local val="${2:-}"
    local func="${3:-${FUNCNAME[1]:-${FUNCNAME[0]}}}"
    local show_help

    show_help=$(cat << EOM
Usage:
  require_optarg [-h|--hlp|--help] opt val [func]

Description:
  Check that an option requiring an argument is followed by a usable value.

Positional arguments:
  1  opt   <str>  Option token being checked (for example, '--threads').
  2  val   <str>  Candidate value following the option.
  3  func  <str>  Optional function name to use in the error message.

Returns:
  0 if 'val' is present, non-empty, and not option-like; otherwise 1.
EOM
    )

    if [[ "${opt}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${opt}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'opt', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    if (( $# < 2 )) || [[ -z "${val}" || "${val}" == -* ]]; then
        echo_err_func "${func}" \
            "option '${opt}' requires a value."
        return 1
    fi
}


#MAYBE: 'check_arg_supplied' overlaps 'validate_var' in 'check_inputs.sh';
#+      may retire one after call sites settle
#NOTE 2026-04-20: has been removed from all repo scripts, etc.
function check_arg_supplied() {
    local asm=""
    local nam=""
    local show_help

    show_help=$(cat << EOM
Usage:
  check_arg_supplied [-h|--hlp|--help] --asm <str> --nam <str>

Description:
  Checks that a given argument has been assigned a value.

Keyword arguments:
  -as, --asm, --asgmt  <str>
    The value (i.e., assignment) of the argument to check.

  -n, --nam, --name  <str>
    The name of the argument to check.

Returns:
  0 if the argument has been assigned a non-empty value; otherwise 1.

Example:
  Check that a required argument ('req_arg') is supplied:
  '''bash
  req_arg="some_value"

  check_arg_supplied -as "\${req_arg}" -n "req_arg"
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
            -as|--asm|--asmgt)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                asm="${2}"
                shift 2
                ;;

            -n|--nam|--name)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
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
        echo_err_func "${FUNCNAME[0]}" \
            "'--name' is required."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  Check that the argument has been supplied; return an informative error
    #+ message if not
    if [[ -z "${asm}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'--${nam}' is required."
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
  check_args_mut_excl [-h|--hlp|--help] nam_1 val_1 nam_2 val_2

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
  check_args_mut_excl "arg_1_name" "\${arg_1_value}" "arg_2_name" "\${arg_2_value}"
  '''
EOM
    )

    #  Print help message if no arguments are passed or if help is requested
    if [[ "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${nam_1}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'nam_1', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${nam_2}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 3, 'nam_2', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  Validate mutually exclusive arguments
    if [[ -n "${val_1}" && -n "${val_2}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "only one of '--${nam_1}' or '--${nam_2}' can be specified at a" \
            "time."
        return 1
    elif [[ -z "${val_1}" && -z "${val_2}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "one of '--${nam_1}' or '--${nam_2}' must be specified."
        return 1
    fi
}


function check_flags_mut_excl() {
    local flg_1="${1:-}"
    local nam_1="${2:-}"
    local flg_2="${3:-}"
    local nam_2="${4:-}"
    local flg_1_lc flg_2_lc
    local show_help

    show_help=$(cat << EOM
Usage:
  check_flags_mut_excl [-h|--hlp|--help] flg_1 nam_1 flg_2 nam_2

Description:
  Checks that two mutually exclusive flags are not specified at the same time. If both are specified or neither is specified, an error is returned.

Positional arguments:
  1  flg_1  <bol>  The first flag to check. Boolean-like string required.
  2  nam_1  <str>  The name of the first flag.
  3  flg_2  <bol>  The second flag to check. Boolean-like string required.
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
  check_flags_mut_excl "\${flg_1}" "flg_1" "\${flg_2}" "flg_2"  # Returns 1 and error message
  '''
EOM
    )

    #  Parse and check function arguments
    if [[ "${flg_1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${flg_1}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'flg_1', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${nam_1}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'nam_1', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${flg_2}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 3, 'flg_2', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${nam_2}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 4, 'nam_2', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    flg_1_lc="${flg_1,,}"
    flg_2_lc="${flg_2,,}"

    case "${flg_1_lc}" in
        t|true)   flg_1=true  ;;
        f|false)  flg_1=false ;;
        *)
            echo_err_func "${FUNCNAME[0]}" \
                "positional argument 1, 'flg_1', must be 'true', 't'," \
                "'false', or 'f': '${flg_1}'."
            return 1
            ;;
    esac

    case "${flg_2_lc}" in
        t|true)   flg_2=true  ;;
        f|false)  flg_2=false ;;
        *)
            echo_err_func "${FUNCNAME[0]}" \
                "positional argument 3, 'flg_2', must be 'true', 't'," \
                "'false', or 'f': '${flg_2}'."
            return 1
            ;;
    esac

    #  Perform the checks
    if [[ "${flg_1}" == "true" && "${flg_2}" == "true" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "only one of '--${nam_1}' or '--${nam_2}' can be specified at a" \
            "time."
        return 1
    # elif [[ "${flg_1}" != "true" && "${flg_2}" != "true" ]]; then
    #     echo_err_func "${FUNCNAME[0]}" \
    #         "one of '--${nam_1}' or '--${nam_2}' must be specified."
    #     return 1
    fi
}


#TODO: audit current usage; keep this helper even if unused
function check_match() {
    local var_dsc="${1:-}"
    local val_num="${2:-}"
    local val_den="${3:-}"
    local out_nam="${4:-}"
    local show_help

    show_help=$(cat << EOM
Usage:
  check_match [-h|--hlp|--help] var_dsc val_num val_den [out_nam]

Description:
  Check that two values match, such as values associated with paired files or numerator/denominator inputs.

Positional arguments:
  1  var_dsc  <str>  Description of the value being compared.
  2  val_num  <str>  First value to compare.
  3  val_den  <str>  Second value to compare.
  4  out_nam  <str>  Optional output variable name to assign if values match.

Returns:
  0 if the values match; otherwise 1.

Notes:
  - If 'out_nam' is supplied and the values match, the matched value is assigned via 'printf -v'.

Example:
  '''bash
  check_match "bin size" "50" "50" "bin_size"
  '''
EOM
    )

    if [[ "${var_dsc}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${var_dsc}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'var_dsc', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${val_num}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'val_num', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${val_den}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 3, 'val_den', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    if [[ -n "${out_nam}" && ! "${out_nam}" =~ ^[a-zA-Z_][a-zA-Z0-9_]*$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 4, 'out_nam', must be a valid shell" \
            "variable name: '${out_nam}'."
        return 1
    fi

    if [[ "${val_num}" == "${val_den}" ]]; then
        if [[ -n "${out_nam}" ]]; then
            printf -v "${out_nam}" '%s' "${val_num}"
        fi
    else
        echo_err_func "${FUNCNAME[0]}" \
            "'${var_dsc}' does not match between numerator (${val_num}) and" \
            "denominator (${val_den})."
        return 1
    fi
}


#TODO: reassess long-term home for helper: e.g., perhaps 'check_inputs.sh';
#+     'check_args.sh' is OK for now because it validates argument-string
#+     formatting
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
  check_str_delim "infiles" ""                               # Returns 1 and an error message
  check_str_delim "infiles" "file1.bam,,file2.bam"           # Returns 1 and an error message
  check_str_delim "scl_fct" "0.5,;1.0"                       # Returns 1 and an error message
  check_str_delim "usr_frg" "300, 200;400"                   # Returns 1 and an error message
  '''
EOM
    )

    if [[ "${nam}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${nam}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'nam', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif (( $# < 2 )); then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'val', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  Check for empty value
    if [[ -z "${val}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "improperly formatted '${nam}' value: '${val}' (empty value)."
        return 1
    fi

    #  Check for improper commas
    if [[ "${val}" =~ ,{2,} || "${val}" =~ ^, || "${val}" =~ ,$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "improperly formatted '${nam}' value: '${val}' (invalid commas)."
        return 1
    fi

    #  Check for improper semicolons
    if [[ "${val}" =~ \;{2,} || "${val}" =~ ^\; || "${val}" =~ \;$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "improperly formatted '${nam}' value: '${val}' (invalid" \
            "semicolons)."
        return 1
    fi

    #  Check for mixed improper delimiters
    if [[ "${val}" =~ ,\; || "${val}" =~ \;, ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "improperly formatted '${nam}' value: '${val}' (mixed invalid" \
            "delimiters)."
        return 1
    fi

    #  Check for spaces before/after commas
    if [[
           "${val}" =~ [[:space:]]+,[[:space:]]* \
        || "${val}" =~ [[:space:]]*,[[:space:]]+
    ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "improperly formatted '${nam}' value: '${val}' (spaces around" \
            "commas)."
        return 1
    fi

    #  Check for spaces before/after semicolons
    if [[
           "${val}" =~ [[:space:]]+\;[[:space:]]* \
        || "${val}" =~ [[:space:]]*\;[[:space:]]+
    ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "improperly formatted '${nam}' value: '${val}' (spaces around" \
            "semicolons)."
        return 1
    fi
}


#  Print an error message when function script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    err_source_only "${BASH_SOURCE[0]}"
fi
