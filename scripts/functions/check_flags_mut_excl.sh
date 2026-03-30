#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: check_flags_mut_excl.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
# 
# Distributed under the MIT license.


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
