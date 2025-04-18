#!/bin/bash

#  compress_remove_files.sh
#  KA


#  Run script in interactive mode (true) or command-line mode (false)
interactive=false

#  Exit on errors or pipe failures if not in "interactive mode"
if ! ${interactive}; then set -eo pipefail; fi

#  Set the path to the "scripts" directory
if ${interactive:-false}; then
    ## WARNING: If 'interactive=true', change path as needed ##
    dir_scr="${HOME}/repos/protocol_chipseq_signal_norm/scripts"
else
    dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi


#  Source and define functions ================================================
#  Set the path to the "functions" directory
dir_fnc="${dir_scr}/functions"

# shellcheck disable=SC1091
{
    source "${dir_fnc}/build_command_find.sh"
    source "${dir_fnc}/check_exists_file_dir.sh"
    source "${dir_fnc}/check_int_pos.sh"
    source "${dir_fnc}/check_mut_excl_flags.sh"
    source "${dir_fnc}/check_program_path.sh"
    source "${dir_fnc}/check_supplied_arg.sh"
    source "${dir_fnc}/echo_error.sh"
    source "${dir_fnc}/handle_env.sh"
}


#  Set paths, values, arguments, etc. for interactive mode
function set_interactive() {
    #  Set hardcoded paths, values, etc.
    ## WARNING: If interactive=true, change values as needed ##
    dir_bas="${HOME}/repos" 
    dir_rep="${dir_bas}/protocol_chipseq_signal_norm"
    nam_qc="04_qc"
    {
        #  Assign variables used in file paths
        aligner="bowtie2"
        a_type="global"
        flag=2
        mapq=1
        details="${aligner}/${a_type}/flag-${flag}_mapq-${mapq}"
        prog="preseq"  # "ssp"
        retain="sc"
    }
    dir_qc="${dir_rep}/${nam_qc}/${details}/${prog}/${retain}"

    #  Set hardcoded argument assignments, etc.
    dir_fnd="${dir_qc}/err_out"
    pattern="*.std???.txt"  # "*.txt"  # "*.bam"
    size=1
    depth=1
    include=""
    exclude="*.gz"
    chk_con=true
    chk_exc=false
}


#  Initialize argument variables, check and parse arguments, etc. =============
#  Initialize hardcoded argument variables
env_nam="env_protocol"

#  Initialize variables along with default assignments
threads=1
dir_fnd=""
pattern="*.std???.txt"
size=1  # Default minimum size to compress: 1 kilobyte
depth=""
include=""
exclude="*.gz"
chk_con=false
chk_exc=false

#  Help message
show_help=$(cat << EOM
Usage:
  compress_remove_files.sh
    --dir_fnd <dir> --pattern <str> --size <int> [--depth <int>]
    [--include <str>] [--exclude <str>] [--chk_con] [--chk_exc]

Description:
  Find files in a specified directory that match the given patterns, compress
  files larger than the specified size, and delete files that are 0 in size.
  The user can also include or exclude files based on name patterns and limit
  the search to a specified directory depth.

Arguments:
   -h, --help     Display this help message and exit (0).
   -t, --threads  Number of threads to use (default: ${threads}).
  -df, --dir_fnd  Directory to search with the program find.
  -pa, --pattern  File pattern, including shell wildcard characters, used in
                  the construction of the find command (default: ${pattern}).
  -sz, --size     Minimum size in kilobytes for compression (default: ${size}).
  -de, --depth    Maximum depth to search within the directory (optional).
  -in, --include  Comma-separated vector of patterns to include with respect to
                  --pattern, including shell wildcards (optional). '--include'
                  is subordinate to '--pattern'.
  -ex, --exclude  Comma-separated vector of patterns to exclude with respect to
                  --pattern, including shell wildcards (optional; default:
                  ${exclude}). '--exclude' is subordinate to '--pattern'.
  -cc, --chk_con  Check the construction of the find command and exit
                  (optional).
  -ce, --chk_exc  Check the construction and execution of the find command and
                  exit (optional).

Dependencies:
  - Programs
    + Bash or Zsh
    + find
    + GNU Parallel (if threads > 1)
    + sort
  - Functions
    + build_command_find
    + check_exists_file_dir
    + check_int_pos
    + check_mut_excl_flags
    + check_program_path
    + check_supplied_arg
    + echo_error
    + handle_env
    + handle_env_activate
    + handle_env_deactivate

Notes:
  - This script does not handle logical OR operations, just AND and AND NOT.
  - compress_remove_files.sh will exit with an error message if it is run from
    the target directory being searched, as doing so causes file-globbing
    errors.
  - The script will exit with an error message if the string assigned to
    '--pattern' matches anything in the current working directory.
  - If --threads is assigned a positive integer greater than 1, then the script
    will use GNU Parallel to parallelize file handling and processing.

Examples:
  \`\`\`
  bash process_files.sh
      --dir_fnd "\${HOME}/path/to/dir"
      --size 2
      --depth 2
      --include "*.log"

  bash process_files.sh
      --threads 2
      --dir_fnd "./err_out"
  \`\`\`
EOM
)

#  Parse arguments
if [[ -z "${1:-}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    echo "${show_help}"
    if ! ${interactive}; then exit 0; fi
fi

if ${interactive:-false}; then
    set_interactive
else
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
             -t|--threads) threads="${2}"; shift 2 ;;
            -df|--dir_fnd) dir_fnd="${2}"; shift 2 ;;
            -pa|--pattern) pattern="${2}"; shift 2 ;;
            -sz|--size)    size="${2}";    shift 2 ;;
            -de|--depth)   depth="${2}";   shift 2 ;;
            -in|--include) include="${2}"; shift 2 ;;
            -ex|--exclude) exclude="${2}"; shift 2 ;;
            -cc|--chk_con) chk_con=true;   shift 1 ;;
            -ce|--chk_exc) chk_exc=true;   shift 1 ;;
            *)
                echo "## Unknown parameter passed: '${1}' ##" >&2
                echo "" >&2
                echo "${show_help}" >&2
                exit 1
                ;;
        esac
    done
fi

#  Check arguments
check_supplied_arg -a "${env_nam}" -n "env_nam"

check_supplied_arg -a "${threads}" -n "threads"
check_int_pos "${threads}" "threads"

check_supplied_arg -a "${dir_fnd}" -n "dir_fnd"
check_exists_file_dir "d" "${dir_fnd}" "dir_fnd"

if [[ "$(realpath "${dir_fnd}")" == "$(realpath "${PWD}")" ]]; then
    echo_error \
        "compress_remove_files.sh cannot be run from the target directory" \
        "being searched, as this causes file-globbing errors. Please run" \
        "the script from a different directory. Currently," \
        "dir_fnd=\"${dir_fnd}\"."
    exit_1
fi

check_supplied_arg -a "${pattern}" -n "pattern"

if compgen -G "${pattern}" > /dev/null; then
    echo_error \
        "The specified pattern '${pattern}' matches files in the current" \
        "working directory. Running this script with such patterns in the" \
        "current directory will lead to downstream errors. Please ensure you" \
        "run the script from a different directory or adjust the pattern."
    exit_1
fi

check_supplied_arg -a "${size}" -n "size"
check_int_pos "${size}" "size"

if [[ -n "${depth}" ]]; then check_int_pos "${depth}" "depth"; fi

check_mut_excl_flags ${chk_con} "chk_con" ${chk_exc} "chk_exc"

#  Activate environment and check that dependencies are in PATH
handle_env "${env_nam}" > /dev/null

check_program_path find
check_program_path parallel
check_program_path sort


#  Do the main work ===========================================================
#  Build find command
# shellcheck disable=SC2046
cmd_find="$(
    build_command_find \
        --dir_fnd "${dir_fnd}" \
        --pattern "${pattern}" \
        $(if [[ -n "${depth}" ]]; then echo "--depth ${depth}"; fi) \
        $(if [[ -n "${include}" ]]; then echo "--include ${include}"; fi) \
        $(if [[ -n "${exclude}" ]]; then echo "--exclude ${exclude}"; fi)
)"

if ${chk_con} || ${chk_exc}; then
    echo "## Call to find for files larger than ${size}k ##"
    echo "${cmd_find} -size +${size}k | sort"
    echo ""

    echo "## Call to find for files with size 0 ##"
    echo "${cmd_find} -size 0 | sort"
    echo ""

    if [[ ${threads} -gt 1 ]]; then
cat << EOM
#  Use GNU Parallel to compress files larger than the size threshold
eval "${cmd_find} -size +${size}k" \\
    | sort \\
    | parallel -j ${threads} gzip

#  Use GNU Parallel to delete empty files
eval "${cmd_find} -size 0" \\
    | sort \\
    | parallel -j ${threads} rm

EOM
    else
cat << EOM
#  Serially compress files larger than the size threshold
eval "${cmd_find} -size +${size}k -exec gzip '{}' \;"

#  Serially delete empty files
eval "${cmd_find} -size 0 -delete"

EOM
    fi
fi

if ${chk_con}; then
    if ! ${interactive}; then exit 0; fi
fi

if ${chk_exc}; then
    echo "## Results of find command for files larger than ${size}k ##"
    eval "${cmd_find} -size +${size}k" | sort
    echo ""

    echo "## Results of find command for files with size 0 ##"
    eval "${cmd_find} -size 0" | sort
    echo ""
fi

if ${chk_exc}; then
    if ! ${interactive}; then exit 0; fi
fi

#  Execute the find command and process the files
if [[ ${threads} -gt 1 ]]; then
    #  Using GNU Parallel, compress files larger than the size threshold and
    #+ delete empty files
    eval "${cmd_find} -size +${size}k" \
        | sort \
        | parallel -j "${threads}" gzip

    eval "${cmd_find} -size 0" \
        | sort \
        | parallel -j "${threads}" rm
else
    #  Serially compress files larger than the size threshold and delete empty
    #+ files
    eval "${cmd_find} -size +${size}k -exec gzip '{}' \;"

    eval "${cmd_find} -size 0 -delete"
fi
