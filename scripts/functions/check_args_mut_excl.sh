#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: check_args_mut_excl.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
# 
# Distributed under the MIT license.


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
