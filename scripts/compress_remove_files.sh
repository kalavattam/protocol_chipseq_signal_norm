#!/bin/bash

#  compress_remove_files.sh
#  KA


#  Run script in interactive/test mode (true) or command-line mode (false)
interactive=false

#  Exit on errors, unset variables, or pipe failures if not in "interactive
#+ mode"
if ! ${interactive}; then
    # set -e
    # set -euo pipefail
    set -eo pipefail
fi

#  Set the path to the "scripts" directory
if ${interactive}; then
    ## WARNING: Change path if you're not Kris and `interactive=true` ##
    dir_sc="${HOME}/tsukiyamalab/Kris/202X_protocol_ChIP/scripts"
else
    dir_sc="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi


#  Source and define functions ================================================
#  Set the path to the "functions" directory
dir_fn="${dir_sc}/functions"

# shellcheck disable=SC1091
{
    source "${dir_fn}/build_command_find.sh"
    source "${dir_fn}/check_exists_file_dir.sh"
    source "${dir_fn}/check_int_pos.sh"
    source "${dir_fn}/check_mut_excl_flags.sh"
    source "${dir_fn}/check_program_path.sh"
    source "${dir_fn}/check_supplied_arg.sh"
    source "${dir_fn}/handle_env.sh"
    source "${dir_fn}/handle_env_activate.sh"
    source "${dir_fn}/handle_env_deactivate.sh"
}


function set_interactive() {
    #  Set hardcoded paths, values, etc.
    ## WARNING: Change the values if you're not Kris and `interactive=true` ##
    dir_bas="${HOME}/tsukiyamalab/Kris" 
    dir_rep="${dir_bas}/202X_protocol_ChIP"
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
env_nam="env_analyze"

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
                  --pattern, including shell wildcards (optional).
  -ex, --exclude  Comma-separated vector of patterns to exclude with respect to
                  --pattern, including shell wildcards (optional; default:
                  ${exclude}).
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
    + handle_env
    + handle_env_activate
    + handle_env_deactivate

Notes:
  - If --threads is assigned a positive integer greater than 1, then the script
    will use GNU Parallel to parallelize file handling and processing.

Examples:
  \`\`\`
  bash process_files.sh
      --dir_fnd "/path/to/dir"
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
if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    echo "${show_help}"
    if ! ${interactive}; then exit 0; fi
fi

if ${interactive}; then
    set_interactive
else
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
            -df|--dir_fnd) dir_fnd="${2}"; shift 2 ;;
            -pa|--pattern) pattern="${2}"; shift 2 ;;
            -sz|--size)    size="${2}";    shift 2 ;;
            -de|--depth)   depth="${2}";   shift 2 ;;
            -in|--include) include="${2}"; shift 2 ;;
            -ex|--exclude) exclude="${2}"; shift 2 ;;
            -cc|--chk_con) chk_con=true;   shift 1 ;;
            -ce|--chk_exc) chk_exc=true;   shift 1 ;;
            *)
                echo "## Unknown parameter passed: ${1} ##" >&2
                echo "" >&2
                echo "${show_help}" >&2
                exit 1
                ;;
        esac
    done
fi

check_supplied_arg -a "${dir_fnd}" -n "dir_fnd"
check_exists_file_dir "d" "${dir_fnd}" "dir_fnd"

check_supplied_arg -a "${pattern}" -n "pattern"

check_supplied_arg -a "${size}" -n "size"
check_int_pos "${size}"

if [[ -n "${depth}" ]]; then check_int_pos "${depth}"; fi

check_mut_excl_flags ${chk_con} "chk_con" ${chk_exc} "chk_exc"

#  Activate environment and check that dependencies are in PATH
handle_env "${env_nam}"

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
