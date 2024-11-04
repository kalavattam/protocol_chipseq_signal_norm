#!/bin/bash

function build_command_find() {
    local dir_fnd
    local pattern
    local follow=false  #TODO Initialize and assign local variable separately
    local depth
    local include
    local exclude
    local show_help
    local array=()
    local find_command
    
    show_help=$(cat << EOM
------------------
build_command_find
------------------

Description:
  build_command_find sonstructs a find command string to locate files in a
  specified directory based on one or more specified patterns. The function
  also allows for the inclusion and/or exclusion of files that match additional
  given patterns.

Keyword parameters:
  -df, --dir_fnd (str): The directory in which to search for files (required).
  -pa, --pattern (str): The pattern to match file names, including shell
                        wildcard characters (required).
  -fl, --follow (flag): Follow symbolic links during the search (optional).
  -de, --depth   (int): The maximum depth to search within the directory
                        (optional).
  -in, --include (str): A comma-separated list of patterns, including shell
                        wildcard characters, to include in the search
                        (optional).
  -ex, --exclude (str): A comma-separated list of patterns, including shell
                        wildcard characters, to exclude from the search
                        (optional).

Returns:
  str: A find command string that can be evaluated to perform the file search.

Dependencies:
  - Bash or Zsh
  - find

Note:
  At its most basic, running build_command_find results in a call to find like
  this:
  \`\`\`
  find "\${dir_fnd}" -type f -name "\${pattern}"
  \`\`\`

  build_command_find is not meant for extensive searches; instead, the function
  is meant to find a relatively small number of, e.g., FASTQ, BAM, or TXT files
  in a single directory with no or minimal subdirectories. Found files are then
  intended for the performance of various operations (in serial or, more
  likely, parallel).

  For extensive searches, consider using
  (a) locate,
  (b) fd (github.com/sharkdp/fd), or
  (c) ripgrep (github.com/BurntSushi/ripgrep; e.g., see
      github.com/BurntSushi/ripgrep/issues/193).

Example:
  \`\`\`
  find_command=\$(
      build_command_find
          --dir_fnd "/path/to/dir"
          --pattern "*.txt"
          --depth 1
          --include "*include_pattern_1*,*include_pattern_2*"
          --exclude "*exclude_pattern_1*,*exclude_pattern_2*"
  )

  eval "\${find_command}" | while IFS= read -r file; do
      echo "Found file: \${file}"
  done
  \`\`\`
EOM
)
    #  Parse and check function parameters
    if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
        echo "${show_help}"
        return 0
    fi

    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
            -df|--dir_fnd) dir_fnd="${2}"; shift 2 ;;
            -pa|--pattern) pattern="${2}"; shift 2 ;;
            -fl|--follow)  follow=true;    shift 1 ;;
            -de|--depth)   depth="${2}";   shift 2 ;;
            -in|--include) include="${2}"; shift 2 ;;
            -ex|--exclude) exclude="${2}"; shift 2 ;;
            *)
                echo "## Unknown parameter passed: ${1} ##" >&2
                echo "" >&2
                echo "${show_help}" >&2
                return 1
                ;;
        esac
    done

    if [[ -z "${dir_fnd}" ]]; then
        echo "Error: --dir_fnd is required." >&2
        echo "" >&2
        return 1
    fi

    if [[ ! -d "${dir_fnd}" ]]; then
        echo \
            "Error: Directory associated with --dir_fnd does not exist:" \
            "${dir_fnd}." >&2
        return 1
    fi  

    if [[ -z "${pattern}" ]]; then
        echo "Error: --pattern is required." >&2
        echo "" >&2
        return 1
    fi

    if [[ -n "${depth}" ]]; then
        if [[ ! "${depth}" =~ ^[1-9][0-9]*$ ]]; then
            echo \
                "Error: --depth was assigned '${depth}' but must be a" \
                "positive integer greater than or equal to 1." >&2
            echo "" >&2
            return 1
        fi
    fi

    #  Dynamically build the find command based on what and how arguments are
    #+ specified
    find_command="find"
    
    if ${follow}; then
        find_command+=" -L"
    fi
    
    find_command+=" \"${dir_fnd}\""

    if [[ -n "${depth}" ]]; then
        find_command+=" -maxdepth ${depth}"
    fi
    
    find_command+=" -type f -name \"${pattern}\""

    if [[ -n "${include}" ]]; then
        IFS=',' read -r -a array <<< "${include}"
        for string in "${array[@]}"; do
            find_command+=" -name \"${string}\""
        done
    fi

    if [[ -n "${exclude}" ]]; then
        IFS=',' read -r -a array <<< "${exclude}"
        for string in "${array[@]}"; do
            find_command+=" ! -name \"${string}\""
        done
    fi

    #  Return the find command as a string to be evaluated
    echo "${find_command}"
}
