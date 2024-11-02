#!/bin/bash

function check_exists_file_dir() {
    local type="${1}"
    local item="${2}"
    local name="${3}"
    local show_help
    local item_type
    local check_flag
    local not_exist_msg

show_help=$(cat << EOM
---------------------
check_exists_file_dir
---------------------

Description:
  check_exists_file_dir checks the existence of a file or directory based on
  the provided type (positional parameter 1: 'f' for file, 'd' for directory).
  If the specified file or directory (positional parameter 2) does not exist,
  an error message is returned along with exit code 1.

Positional parameters:
  1, type (str): The type to check for existence. Use 'f' for file or 'd' for
                 directory (required).
  2, item (str): The file or directory, including its path, to check
                 (required).
  3, name (str): An optional name to associate with the item for error
                 messages (optional).

Returns:
  0 if the file or directory exists; otherwise, prints an error message and
  returns exit code 1.

Dependencies:
  - Bash or Zsh

Examples:
  \`\`\`
  #  Check if a file exists
  check_exists_file_dir "f" "/path/to/file.txt"

  #  Check if a directory exists
  check_exists_file_dir "d" "/path/to/directory"

  #  Check if a file exists with an associated name for error message
  check_exists_file_dir "f" "/path/to/file.txt" "infile"

  #  Check if a directory exists with an associated name for error message
  check_exists_file_dir "d" "/path/to/directory" "outdir"
  \`\`\`
EOM
    )

    #  Parse and check function parameters
    if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
        echo "${show_help}"
        return 0
    fi

    #  Validate variable type
    if [[ "${type}" != "f" && "${type}" != "d" ]]; then
        echo \
            "Error: Positionl parameter 1, type, is invalid: '${type}'." \
            "Expected 'f' for file or 'd' for directory." >&2
        echo "" >&2
        return 1
    fi

    #  Check that variable item is defined and not empty
    if [[ -z "${item}" ]]; then
        echo \
            "Error: Positional parameter 2, item, is not defined or is" \
            "empty." >&2
        echo "" >&2
        return 1
    fi

    #  Check file or directory existence; construct and return error message if
    #+ applicable
    if [[ "${type}" == "f" ]]; then
        item_type="File"
        check_flag="-f"
        not_exist_msg="File does not exist"
    elif [[ "${type}" == "d" ]]; then
        item_type="Directory"
        check_flag="-d"
        not_exist_msg="Directory does not exist"
    fi

    if [[ -n "${name}" ]]; then
        not_exist_msg="${item_type} associated with --${name} does not exist"
    fi

    if ! eval "[[ ${check_flag} \"${item}\" ]]"; then
        echo "Error: ${not_exist_msg}: ${item}." >&2
        return 1
    fi
}
