#!/bin/bash
# -*- coding: utf-8 -*-
# 
# Script: check_format_time.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
# 
# Distributed under the MIT license.

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
