#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: compress_remove_files.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
#
# Distributed under the MIT license.


#  Require Bash >= 4.4 before doing any work
if [[ -z "${BASH_VERSION:-}" ]]; then
    echo "error(shell):" \
        "this script must be run under Bash >= 4.4." >&2
    exit 1
elif ((
    BASH_VERSINFO[0] < 4 || ( BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] < 4 )
)); then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "this script requires Bash >= 4.4; current version is" \
        "'${BASH_VERSION}'." >&2
    exit 1
fi

#  Run in safe mode, exiting on errors, unset variables, and pipe failures
set -euo pipefail

#  Set the path to the 'scripts' directory
dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"


#  Source and define functions ================================================
dir_fnc="${dir_scr}/functions"
fnc_src="${dir_fnc}/source_helpers.sh"

if [[ ! -f "${fnc_src}" ]]; then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "script not found: '${fnc_src}'." >&2
    exit 1
fi

# shellcheck disable=SC1090
source "${fnc_src}" || {
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "failed to source '${fnc_src}'." >&2
    exit 1
}

source_helpers "${dir_fnc}" \
    check_args \
    check_env \
    check_inputs \
    check_numbers \
    construct_find \
    format_outputs \
    handle_env \
    help/help_compress_remove_files \
    || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source required helper scripts." >&2
        exit 1
    }

unset fnc_src


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

#  Parse arguments
if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
    help_compress_remove_files >&2
    exit 0
fi

while [[ "$#" -gt 0 ]]; do
    case "${1}" in
        -t|--thr|--threads)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_compress_remove_files >&2
                exit 1
            }
            threads="${2}"
            shift 2
            ;;

        -df|--dir[_-]fnd)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_compress_remove_files >&2
                exit 1
            }
            dir_fnd="${2}"
            shift 2
            ;;

        -pa|--pttrn|--pattern)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_compress_remove_files >&2
                exit 1
            }
            pattern="${2}"
            shift 2
            ;;

        -sz|--size)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_compress_remove_files >&2
                exit 1
            }
            size="${2}"
            shift 2
            ;;

        -de|--dpth|--depth)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_compress_remove_files >&2
                exit 1
            }
            depth="${2}"
            shift 2
            ;;

        -in|--incld|--include)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_compress_remove_files >&2
                exit 1
            }
            include="${2}"
            shift 2
            ;;

        -ex|--excld|--exclude)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_compress_remove_files >&2
                exit 1
            }
            exclude="${2}"
            shift 2
            ;;

        -cn|--chk[_-]con)
            chk_con=true
            shift 1
            ;;

        -ce|-cu|--chk[_-]exc|--chk[_-]exu)
            chk_exc=true
            shift 1
            ;;

        *)
            echo_err "unknown option/parameter passed: '${1}'."
            echo >&2
            help_compress_remove_files >&2
            exit 1
            ;;
    esac
done

#  Check arguments
validate_var "env_nam" "${env_nam}"
check_env_installed "${env_nam}"

validate_var "threads" "${threads}"
check_int_pos "${threads}" "threads"

validate_var_dir "dir_fnd" "${dir_fnd}"

dir_fnd_real="$(realpath "${dir_fnd}")"
dir_pwd_real="$(realpath "${PWD}")"
if \
       [[ "${dir_pwd_real}" == "${dir_fnd_real}" ]] \
    || [[ "${dir_pwd_real}" == "${dir_fnd_real}"/* ]]
then
    echo_err \
        "'compress_remove_files.sh' cannot be run from the target directory" \
        "being searched, or from one of its subdirectories. This restriction" \
        "is kept as a conservative safety policy. Please run the script from" \
        "a different directory. Currently, 'dir_fnd=${dir_fnd}' and" \
        "'PWD=${PWD}'."
    exit 1
fi

validate_var "pattern" "${pattern}"
if compgen -G "${pattern}" > /dev/null; then
    echo_warn \
        "the specified pattern '${pattern}' matches one or more paths in the" \
        "current working directory. This is no longer treated as a hard" \
        "error, but it can still be a sign that the script is being run from" \
        "an inconvenient location. Current 'PWD=${PWD}'."
fi

validate_var "size" "${size}"
check_int_pos "${size}" "size"

if [[ -n "${depth}" ]]; then check_int_pos "${depth}" "depth"; fi

check_flags_mut_excl ${chk_con} "chk_con" ${chk_exc} "chk_exc"

#  Activate environment and check that dependencies are in PATH
env_msg=(
    "'handle_env' failed for 'env_nam=${env_nam}'. Check that Conda/Mamba are"
    "available and that the environment exists."
)

if ! handle_env "${env_nam}" > /dev/null 2>&1; then
    echo_err "${env_msg[*]}"
    exit 1
fi

check_pgrm_path find
check_pgrm_path gzip
if [[ "${threads}" -gt 1 ]]; then check_pgrm_path parallel; fi
check_pgrm_path sort


#  Do the main work ===========================================================
#  Build base 'find' command as an array
unset arr_cmd_find && declare -a arr_cmd_find

arr_build=(
    --arr_nam arr_cmd_find
    --dir_fnd "${dir_fnd}"
    --pattern "${pattern}"
)

if [[ -n "${depth}" ]]; then
    arr_build+=( --depth "${depth}" )
fi

if [[ -n "${include}" ]]; then
    arr_build+=( --include "${include}" )
fi

if [[ -n "${exclude}" ]]; then
    arr_build+=( --exclude "${exclude}" )
fi

construct_find "${arr_build[@]}"
unset arr_build

arr_cmd_big=( "${arr_cmd_find[@]}" -size "+${size}k" )
arr_cmd_zero=( "${arr_cmd_find[@]}" -size 0 )

if [[ "${chk_con}" == "true" || "${chk_exc}" == "true" ]]; then
    echo "## Call to find for files larger than ${size}k ##"
    print_cmd_array arr_cmd_big
    echo
    echo

    echo "## Call to find for files with size 0 ##"
    print_cmd_array arr_cmd_zero
    echo
    echo

    if [[ "${threads}" -gt 1 ]]; then
        cat << EOM
## Use GNU Parallel to compress files larger than the size threshold ##
"\${arr_cmd_big[@]}" | sort | parallel -j ${threads} gzip

    which is

"${arr_cmd_big[@]}" | sort | parallel -j ${threads} gzip


## Use GNU Parallel to delete empty files ##
"\${arr_cmd_zero[@]}" | sort | parallel -j ${threads} rm -- {}

    which is

"${arr_cmd_zero[@]}" | sort | parallel -j ${threads} rm -- {}


EOM
    else
        cat << EOM
## Serially compress files larger than the size threshold ##
"\${arr_cmd_big[@]}" -exec gzip '{}' \;

    which is

"${arr_cmd_big[@]}" -exec gzip '{}' \;


## Serially delete empty files ##
"\${arr_cmd_zero[@]}" -delete

    which is

"${arr_cmd_zero[@]}" -delete


EOM
    fi
fi

if [[ "${chk_con}" == "true" ]]; then exit 0; fi

if [[ "${chk_exc}" == "true" ]]; then
    echo "## Results of find command for files larger than ${size}k ##"
    "${arr_cmd_big[@]}" | sort
    echo
    echo

    echo "## Results of find command for files with size 0 ##"
    "${arr_cmd_zero[@]}" | sort
    echo
    echo
fi

if [[ "${chk_exc}" == "true" ]]; then exit 0; fi

if [[ "${threads}" -gt 1 ]]; then
    "${arr_cmd_big[@]}" \
        | sort \
        | parallel -j "${threads}" gzip

    "${arr_cmd_zero[@]}" \
        | sort \
        | parallel -j "${threads}" rm -- {}
else
    "${arr_cmd_big[@]}" -exec gzip '{}' \;
    "${arr_cmd_zero[@]}" -delete
fi
