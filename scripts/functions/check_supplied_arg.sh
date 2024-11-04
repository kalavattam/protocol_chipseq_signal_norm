#!/bin/bash

function check_supplied_arg() {
    local asgmt
    local name
    local show_help

    show_help=$(cat << EOM
------------------
check_supplied_arg
------------------

Description:
  check_supplied_arg checks that a given argument has been assigned a value.
  If the argument is empty or not set, it prints an error message and return
  exit code 1.

Keyword parameters:
  -a, --asgmt  The value (i.e., assignment or "asgmt") of the argument to
               check (required).
  -n, --name   The name of the argument to check (required).

Returns:
  0 if the argument has been assigned a value; otherwise, returns exit code 1
  if the argument is empty or not set.

Dependencies:
  - Bash or Zsh

Example:
  \`\`\`
  #  Check that a required argument is supplied:
  required_arg="some_value"
  
  check_supplied_arg -a "\${required_arg}" -n "required_arg"
  \`\`\`
EOM
    )

    #  Parse and check function parameters
    if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
        echo "${show_help}"
        return 0
    fi

    #  Parse and check function parameters
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
            -a|--asgmt) asgmt="${2}"; shift 2 ;;
            -n|--name)  name="${2}";  shift 2 ;;
            *)
                echo "## Unknown parameter passed: ${1} ##" >&2
                echo "" >&2
                echo "${show_help}" >&2
                return 1
                ;;
        esac
    done

    #  Check that required parameters are supplied
    if [[ -z "${name}" ]]; then
        echo "Error: --name is required." >&2
        echo "${show_help}"
        echo "" >&2
        return 1
    fi

    #  Check that the argument has been supplied; return an informative error
    #+ message if not
    if [[ -z "${asgmt}" ]]; then
        echo "Error: --${name} is required." >&2
        return 1
    fi
}