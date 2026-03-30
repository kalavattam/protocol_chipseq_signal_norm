#!/bin/bash
# -*- coding: utf-8 -*-
# 
# Script: check_int_nonneg.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
# 
# Distributed under the MIT license.


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
