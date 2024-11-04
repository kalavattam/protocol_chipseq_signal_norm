#!/bin/bash

function check_format_time() {
    local time="${1}"
    local show_help

    show_help=$(cat << EOM
-----------------
check_format_time
-----------------

Description:
  This function checks if a variable is assigned a time formatted as 'h:mm:ss' or 'hh:mm:ss'.

Positional parameter:
  1, time (str): The time string to check (required).

Returns:
  0 if the time string is correctly formatted; otherwise, returns 1.

Dependencies:
  - Bash or Zsh

Example:
  \`\`\`
  check_format_time "2:30:15"  # Returns 0

  check_format_time "25:165:80"  # Returns 1 (invalid time format)
  \`\`\`
EOM
    )

    #  Parse and check function parameter
    if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
        echo "${show_help}"
        return 1
    fi

    #  Regular expression to check the time format 'h:mm:ss'
    if [[ "${time}" =~ ^[0-9]{1,2}:[0-5][0-9]:[0-5][0-9]$ ]]; then
        return 0
    else
        echo \
            "Error: '${time}' is not a valid time format. Expected format is" \
            "'h:mm:ss' or 'hh:mm:ss'." >&2
        return 1
    fi
}
