#!/bin/bash
# -*- coding: utf-8 -*-
# 
# Script: check_flt_nonneg.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.4) was used in development.
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