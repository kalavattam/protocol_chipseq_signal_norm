#!/bin/bash
# -*- coding: utf-8 -*-
# 
# Script: check_str_delim.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
# 
# Distributed under the MIT license.

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
