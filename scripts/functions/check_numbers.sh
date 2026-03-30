#!/bin/bash
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


function check_flt_nonneg() {
    local val="${1:-}"
    local nam="${2:-}"
    local show_help

    show_help=$(cat << EOM
Usage:
  check_flt_nonneg [-h|--hlp|--help] val [nam]

Description:
  Checks that a value is a non-negative integer or float.

Positional arguments:
  1  val  <flt>  The value to check.
  2  nam  <str>  Name of argument or variable associated with the value (optional).

Returns:
  0 if the check passes; 1 if it doesn’t.

#TODO:
  - Add usage example(s).
EOM
    )

    #  Parse and check function arguments
    if [[ "${val}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${val}" ]]; then
        echo "Error: Positional argument 1, 'val', is missing." >&2
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    if [[ "${val}" =~ ^[+]?([0-9]+([.][0-9]*)?|[.]?[0-9]+)$ ]]; then
        return 0
    else
        if [[ -n "${nam}" ]]; then
            echo \
                "Error: '--${nam}' was assigned '${val}' but must be a" \
                "non-negative integer or float." >&2
        else
            echo \
                "Error: '${val}' is not a non-negative integer or float." >&2
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
  Checks that a value is a positive integer or float.

Positional arguments:
  1  val  <flt>  The value to check.
  2  nam  <str>  Name of argument or variable associated with the value (optional).

Returns:
  0 if the check passes; 1 if not.

Examples:
  Check that the value assigned to the argument '--scaling_factor' is a positive float:
  '''bash
  check_flt_pos 4.5 "scaling_factor"  # Returns 0

  check_flt_pos a "scaling_factor"  # Returns 1 and the below message
  # Error: '--scaling_factor' was assigned 'a' but must be a positive integer or float.

  check_flt_pos -2.5  # Returns 1 and the below message
  # Error: '-2.5' is not a positive integer or float.
  '''
EOM
    )

    #  Parse and check function arguments
    if [[ "${val}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${val}" ]]; then
        echo "Error: Positional argument 1, 'val', is missing." >&2
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
            echo \
                "Error: '--${nam}' was assigned '${val}' but must be a" \
                "positive integer or float." >&2
        else
            echo \
                "Error: '${val}' is not a positive integer or float." >&2
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
  0 if the time string is correctly formatted; otherwise, 1 and an error message.

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
        echo "Error: Positional argument 1, 'time', is missing." >&2
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  Check for 'mm:ss', 'h:mm:ss', or 'hh:mm:ss'
    if [[ ! "${time}" =~ ^([0-9]{1,2}:)?[0-5][0-9]:[0-5][0-9]$ ]]; then
        echo \
            "Error: '${time}' is not a valid time format. Expected format is" \
            "'mm:ss', 'h:mm:ss', or 'hh:mm:ss'." >&2
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

Dependencies:
  - Bash or Zsh

Examples:
  #TODO: Brief description of this example
  '''bash
  ❯ check_int_nonneg 4 "threads"  # Returns 0
  '''

  '''txt
  0 ❯
  '''

  #TODO: Brief description of this example
  '''bash
  check_int_nonneg a "threads"  # Returns 1 and the below message
  '''

  '''txt
  1 ❯ Error: --threads was assigned 'a' but must be a non-negative integer.
  '''

  #TODO: Brief description of this example
  '''bash
  check_int_nonneg -2  # Returns 1 and the below message
  '''

  '''txt
  1 ❯ Error: -2 is not a non-negative integer.
  '''
EOM
    )

    #  Parse and check function arguments, printing help message as appropriate
    if [[ "${val}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${val}" ]]; then
        echo "Error: Positional argument 1, 'val', is missing." >&2
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  Perform the check and return an error message if it fails
    if ! [[ "${val}" =~ ^[0-9]+$ ]]; then
        if [[ -n "${nam}" ]]; then
            echo \
                "Error: '--${nam}' was assigned '${val}' but must be a" \
                "non-negative integer." >&2
        else
            echo \
                "Error: '${val}' is not a non-negative integer." >&2
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
  #TODO: Brief description of this example
  '''bash
  check_int_pos 4 "threads"  # Returns 0
  '''

  '''txt
  0 ❯
  '''

  #TODO: Brief description of this example
  '''bash
  check_int_pos a "threads"  # Returns 1 and the below message
  '''

  '''txt
  1 ❯ Error: '--threads' was assigned 'a' but must be a positive integer greater than or equal to 1.
  '''

  #TODO: Brief description of this example
  '''bash
  check_int_pos -2  # Returns 1 and the below message
  '''

  '''txt
  1 ❯ Error: '-2' is not a positive integer greater than or equal to 1.
  '''
EOM
    )

    #  Parse and check function arguments, printing help message as appropriate
    if [[ "${val}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${val}" ]]; then
        echo "Error: Positional argument 1, 'val', is missing." >&2
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  Perform the check and return an error message if it fails
    if ! [[ "${val}" =~ ^[1-9][0-9]*$ ]]; then
        if [[ -n "${nam}" ]]; then
            echo \
                "Error: '--${nam}' was assigned '${val}' but must be an" \
                "integer greater than or equal to 1." >&2
        else
            echo \
                "Error: '${val}' is not an integer greater than or equal to" \
                "1." >&2
        fi
        return 1
    fi
}


function check_scl_fct() {
    local scl_fct="${1}"  # Comma-separated scaling factors to validate
    local entries=()      # Array for individual comma-separated components
    local entry           # Individual scaling factor: 'num:den' or 'num'
    local num den         # Values for validation and formatting
    local val_fmt=()      # Array for validated and formatted entries
    local out_str         # Final formatted, comma-separated string output
    local show_help       # Help message/function documentation

    show_help=$(cat << EOM
-------------
check_scl_fct
-------------

Description:
  Function (subroutine) to validate and format positive floats for precomputed
  scaling factors used with deepTools 'bamCompare' or 'bigwigCompare'. This
  subroutine ensures the following:
    - Scaling factors are provided as a single comma-separated string.
    - Each factor pair (numerator and denominator) is in the format 'num:den'.
    - Single values (numerators) are automatically paired with a denominator
      of 1.
    - All values are validated as positive floats.
    - Leading decimal points (e.g., ".5") are zero-padded (e.g., "0.5").

Positional argument:
  1, scl_fct (str): Comma-separated scaling factors to validate.

Returns:
  0 and a validated, comma-separated string of formatted scaling factors;
  example: '0.4:1,1.2:1,0.65:0.8'. Otherwise, returns 1 with an error message.

Usage:
  check_scl_fct "\${scl_fct}"

Examples:
  '''txt
  ❯ scl_fct=0.5,1.2:3,.65
  ❯ check_scl_fct \${scl_fct}
  0.5:1,1.2:3,0.65:1

  ❯ scl_fct=0.00456,0.00789
  ❯ check_scl_fct \${scl_fct}
  0.00456:1,0.00789:1

  ❯ scl_fct="num:den"
  ❯ check_scl_fct \${scl_fct}
  Error: Invalid 'num' in 'num:den' (must be a positive float).

  ❯ scl_fct="0.33:den"
  ❯ check_scl_fct \${scl_fct}
  Error: Invalid 'den' in '0.33:den' (must be a positive float).

  ❯ scl_fct=0.5,1.2a,.65
  ❯ check_scl_fct \${scl_fct}
  Error: Invalid numerator-only scaling factor '1.2a' (must be a positive float).

  ❯ scl_fct=hotdog
  ❯ check_scl_fct \${scl_fct}
  Error: Invalid numerator-only scaling factor 'hotdog' (must be a positive float)
  '''
EOM
    )

    #  Parse and check function parameter
    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    #  Split the input string by commas
    IFS=',' read -ra entries <<< "${scl_fct}"

    for entry in "${entries[@]}"; do
        #  Check that entry is a 'num:den' pair (containing a colon)
        if [[ "${entry}" == *:* ]]; then
            #  Split into num and den
            IFS=':' read -r num den <<< "${entry}"

            #  Validate that both 'num' and 'den' are positive floats
            if [[ ! "${num}" =~ ^[0-9]*\.?[0-9]+$ ]]; then
                echo \
                    "Error: Invalid 'num' in '${entry}' (must be a positive" \
                    "float)." >&2
                return 1
            fi

            if [[ ! "${den}" =~ ^[0-9]*\.?[0-9]+$ ]]; then
                echo \
                    "Error: Invalid 'den' in '${entry}' (must be a positive" \
                    "float)." >&2
                return 1
            fi

            #  0-pad values starting with a decimal point
            if [[ "${num}" == .* ]]; then num="0${num}"; fi
            if [[ "${den}" == .* ]]; then den="0${den}"; fi

            #  Add back to the formatted values as 'num:den'
            val_fmt+=( "${num}:${den}" )
        else
            #  Validate single numerator-only entry as a positive float
            if [[ ! "${entry}" =~ ^[0-9]*\.?[0-9]+$ ]]; then
                echo \
                    "Error: Invalid numerator-only scaling factor '${entry}'" \
                    "(must be a positive float)." >&2
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
