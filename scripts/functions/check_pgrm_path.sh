#!/bin/bash
# -*- coding: utf-8 -*-
# 
# Script: check_pgrm_path.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
# 
# Distributed under the MIT license.

function check_pgrm_path() {
    local prog="${1:-}"
    local show_help

    show_help=$(cat << EOM
Usage:
  check_pgrm_path [-h|--hlp|--help] prog

Description:
  Checks that a given program is available in PATH.

Positional argument:
  1  prog  <str>  Name of program to check.

Returns:
  0 if program is found in PATH; otherwise, 1 with an error message.

Dependencies:
  - Bash

Example:
  1. Confirm success when program 'samtools' is in PATH
  '''bash
  check_pgrm_path "samtools"
  '''

  '''txt
  0 ❯
  '''

  2. Confirm error when program 'samtools' is not in PATH
  '''bash
  check_pgrm_path "samtools"
  '''

  '''txt
  1 ❯ Error: 'samtools' is not in PATH. Please install 'samtools' and/or add it to PATH.
  '''
EOM
    )

    #  Parse and check function arguments
    if [[ "${prog}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${prog}" ]]; then
        echo "Error: Positional argument 1, 'prog', is missing." >&2
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    if ! command -v "${prog}" &> /dev/null; then
        echo \
            "Error: '${prog}' is not in PATH. Please install '${prog}'" \
            "and/or add it to PATH." >&2
        return 1
    fi
}
