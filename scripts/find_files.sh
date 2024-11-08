#!/bin/bash

#  find_files.sh
#  KA


#  Run script in interactive/test mode (true) or command-line mode (false)
interactive=false

#  Exit on errors or pipe failures if not in "interactive mode"
if ! ${interactive}; then set -eo pipefail; fi

#  Set the path to the "scripts" directory
if ${interactive}; then
    ## WARNING: Change path if you're not Kris and `interactive=true` ##
    dir_scr="${HOME}/tsukiyamalab/Kris/202X_protocol_ChIP/scripts"
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
    source "${dir_fnc}/exit_1.sh"
    source "${dir_fnc}/pair_fastqs.sh"
}


set_interactive() {
    #  Hardcoded paths
    ## WARNING: Change the values if you're not Kris and `interactive=true` ##
    dir_rep="${HOME}/tsukiyamalab/Kris/202X_protocol_ChIP"
    dir_fil="${dir_rep}/data/symlinked"

    #  Hardcoded argument assignments
    dir_fnd="${dir_rep}/${dir_fil}"
    pattern="*.fastq.gz"  # "*.bam"
    follow=true  # false
    depth=1
    include=""
    exclude=""
    fastqs=true  # false
    chk_con=false
    chk_exc=true
}


#  Initialize argument variables, check and parse arguments, etc. =============
#  Initialize variables along with default assignments
dir_fnd=""
pattern=""
depth=""
follow=false
fastqs=false
include=""
exclude=""
chk_con=false
chk_exc=false

show_help=$(cat << EOM
Usage:
  find_files.sh
    --dir_fnd <str> --pattern <str> [--depth <int>] [--follow] [--fastqs]
    [--include <str>] [--exclude <str>] [--chk_con] [--chk_exc]

Description:
  Search for files in a specified directory using the *nix 'find' command.
  find_files.sh returns the files as a single comma-separated string to stdout.

Arguments:
   -h, --help     Display this help message and exit (0).
  -df, --dir_fnd  Directory to search with the program find.
  -pa, --pattern  Comma-separated serialized list of file patterns, including
                  shell wildcard characters, used in the construction of the
                  find command.
  -de, --depth    Maximum depth to search within the directory (optional).
  -fl, --follow   Follow symbolic links during the search (optional).
  -fq, --fastqs   Find FASTQ files, returning them in a semicolon-separated
                  string with paired-end sequenced files in comma-separated
                  substrings (optional).
  -in, --include  Comma-separated serialized list of patterns to include,
                  including shell wildcards (optional).
  -ex, --exclude  Comma-separated serialized list of patterns to exclude,
                  including shell wildcards (optional).
  -cc, --chk_con  Check the construction of the find command and exit
                  (optional).
  -ce, --chk_exc  Check the construction and execution of the find command and
                  exit (optional).

Dependencies:
  - Programs
    + Bash or Zsh
    + find
    + paste
    + sed
    + sort
  - Functions
    + build_command_find
    + check_exists_file_dir
    + check_int_pos
    + check_mut_excl_flags
    + check_program_path
    + check_supplied_arg
    + echo_error
    + exit_1
    + pair_fastqs

Notes:
  - This script doesn't handle logical OR operations, just AND and AND NOT.
  - find_files.sh will exit with an error message if it is run from the target
    directory being searched, as doing so causes file-globbing errors.

Examples:
  \`\`\`
  #  For BAM files
  bash /home/user/path/to/find_files.sh
      --dir_fnd "/path/to/directory"
      --pattern "*.bam"
      --depth 1
      --include "*Hho1*,*Q*"
      --exclude "*Hmo1*,*G2M*,*G1*"

  #  For FASTQ files
  bash /home/user/path/to/find_files.sh
      --dir_fnd "/path/to/another/directory"
      --pattern "*.fastq.gz"
      --follow
      --fastqs
  \`\`\`
EOM
)

#  Parse arguments
if [[ -z "${1:-}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
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
            -de|--depth)   depth="${2}";   shift 2 ;;
            -fl|--follow)  follow=true;    shift 1 ;;
            -fq|--fastqs)  fastqs=true;    shift 1 ;;
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

#  Check arguments
check_supplied_arg -a "${dir_fnd}" -n "dir_fnd"
check_exists_file_dir "d" "${dir_fnd}" "dir_fnd"

if [[ "$(realpath "${dir_fnd}")" == "$(realpath "${PWD}")" ]]; then
    echo_error \
        "find_files.sh cannot be run from the target directory" \
        "being searched, as this causes file-globbing errors. Please run" \
        "the script from a different directory. Currently," \
        "dir_fnd=\"${dir_fnd}\"."
    exit_1
fi

check_supplied_arg -a "${pattern}" -n "pattern"

if [[ -n "${depth}" ]]; then check_int_pos "${depth}" "depth"; fi

check_mut_excl_flags ${chk_con} "chk_con" ${chk_exc} "chk_exc"

#  Check programs
check_program_path find
check_program_path paste
check_program_path sed
check_program_path sort


#  Do the main work ===========================================================
#  Build find command
# shellcheck disable=SC2046
cmd_find="$(
    build_command_find \
        --dir_fnd "${dir_fnd}" \
        --pattern "${pattern}" \
        $(if ${follow}; then echo "--follow"; fi) \
        $(if [[ -n "${depth}" ]]; then echo "--depth ${depth}"; fi) \
        $(if [[ -n "${include}" ]]; then echo "--include ${include}"; fi) \
        $(if [[ -n "${exclude}" ]]; then echo "--exclude ${exclude}"; fi)
)"

if ${chk_con} || ${chk_exc}; then
    echo "## Call to find ##"

    if ${fastqs}; then
cat << EOM
eval "${cmd_find}" \\
    | sort \\
    | pair_fastqs \\
    | paste -s \\
    | sed \\
        -e 's:\t::g' \\
        -e 's/;$//'

EOM
    else
cat << EOM
eval "${cmd_find}" \\
    | sort \\
    | paste -sd "," -

EOM
    fi
fi

if ${chk_con}; then
    if ! ${interactive}; then exit 0; fi
fi

if ${chk_exc}; then
    echo "## Results of find command ##"
    eval "${cmd_find}" | sort
    echo ""

    if ${fastqs}; then
        echo \
            "## Results of find command as single semicolon- and" \
            "comma-separated string ##"
        eval "${cmd_find}" \
            | sort \
            | pair_fastqs \
            | paste -s \
            | sed \
                -e 's:\t::g' \
                -e 's/;$//'
    else
        echo "## Results of find command as single comma-separated string ##"
        eval "${cmd_find}" \
            | sort \
            | paste -sd "," -
        echo ""
    fi
fi

if ${chk_exc}; then
    if ! ${interactive}; then exit 0; fi
fi

if ${fastqs}; then
    eval "${cmd_find}" \
        | sort \
        | pair_fastqs \
        | paste -s \
        | sed \
            -e 's:\t::g' \
            -e 's/;$//'
else
    eval "${cmd_find}" \
        | sort \
        | paste -sd "," -
fi
