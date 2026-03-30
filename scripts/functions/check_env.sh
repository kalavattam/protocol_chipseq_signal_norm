#!/bin/bash
# -*- coding: utf-8 -*-
# 
# Script: check_env.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
# 
# Distributed under the MIT license.


function check_env_installed() {
    local env_nam="${1:-}"
    local quiet="${2:-false}"
    local quiet_lc
    local show_help

    show_help=$(cat << EOM
Usage:
  check_env_installed [-h|--hlp|--help] env_nam [quiet]

Description:
  Check that a specific Conda (or Mamba) environment is installed.

Positional arguments:
  1  env_nam  <str>  Name of Conda (or Mamba) environment to check.
  2  quiet    <bol>  If 'true', suppress error messages and return status only (default: 'false').

Returns:
  0 if the environment is installed; otherwise, 1 and, unless 'quiet=true', an error message.

Dependencies:
  - Bash or Zsh
  - Conda

Note:
  Quiet mode does not suppress the "'conda' is not available in PATH" message, as the purpose of the script is to suppress the expected "env not installed" message during Boolean testing.

Examples:
  1. Check that environment "env_protocol" is installed
  '''bash
  check_env_installed "env_protocol"
  '''

  '''txt
  0 ❯
  '''

  2. Confirm that fictitious environment "env_nonexistent" errors
  '''bash
  check_env_installed "env_nonexistent"
  '''

  '''txt
  1 ❯ Error: Environment 'env_nonexistent' appears not to be installed.
  '''

  3. Confirm that fictitious environment "env_nonexistent" errors quietly
  '''bash
  check_env_installed "env_nonexistent" "true"
  '''

  '''txt
  1 ❯
  '''
EOM
    )

    #  Lowercase-convert 'quiet' assignment
    quiet_lc=$(printf '%s' "${quiet}" | tr '[:upper:]' '[:lower:]')

    #  Parse and check function arguments
    if [[ "${env_nam}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${env_nam}" ]]; then
        if [[ ! "${quiet_lc}" =~ ^(t|true)$ ]]; then
            echo "Error: Positional argument 1, 'env_nam', is missing." >&2
            echo >&2
            echo "${show_help}" >&2
        fi
        return 1
    fi

    #  Check availability of 'conda'
    if ! command -v conda >/dev/null 2>&1; then
        echo "Error: 'conda' is not available in PATH." >&2
        return 1
    fi

    #  Check that the specific Conda or Mamba environment is installed
    if ! \
        conda env list | awk -v env="${env_nam}" '
            $1 == env { found = 1 } END { exit !found }
        '
    then
        if [[ ! "${quiet_lc}" =~ ^(t|true)$ ]]; then
            echo \
                "Error: Environment '${env_nam}' appears not to be" \
                "installed." >&2
        fi
        return 1
    fi
}


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
