#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: find_files.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
# GPT-5.6) were used in design, development, and documentation, with all output
# reviewed, edited, and approved by the author.
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

#  Set path to the 'scripts' directory
dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"


#  Source shared helpers
function source_helpers_script() {
    local fnc_src

    dir_fnc="${dir_scr}/../lib/bash"
    fnc_src="${dir_fnc}/core/source_helpers.sh"

    if [[ ! -f "${fnc_src}" ]]; then
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "script not found: '${fnc_src}'." >&2
        return 1
    fi

    # shellcheck disable=SC1090
    source "${fnc_src}" || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source '${fnc_src}'." >&2
        return 1
    }

    #TODO: source function to normalize Boolean strings (do not rely on just
    #+     "true", "false")
    source_helpers "${dir_fnc}" \
        check_args \
        check_env \
        check_inputs \
        check_numbers \
        construct_find \
        format_outputs \
        help/help_find_files \
        process_sequences \
        || {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to source required helper scripts." >&2
            return 1
        }
}


function init_arg_defs() {
    dir_fnd=""
    pattern=""
    depth=""
    follow=false
    fastqs=false
    include=""
    exclude=""
    chk_con=false
    chk_exc=false
    shw_str=false

    unset arr_cmd_find
    declare -ga arr_cmd_find
}


function parse_args() {
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
            -df|--dir[_-]fnd)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_find_files >&2
                    return 1
                }
                dir_fnd="${2}"
                shift 2
                ;;

            -pa|--pttrn|--pattern)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_find_files >&2
                    return 1
                }
                pattern="${2}"
                shift 2
                ;;

            -de|--dpth|--depth)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_find_files >&2
                    return 1
                }
                depth="${2}"
                shift 2
                ;;

            -fl|-sy|--fllw|--follow|--symlink)
                follow=true
                shift 1
                ;;

            -fq|--fqs|--fastqs)
                fastqs=true
                shift 1
                ;;

            -in|--incld|--include)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_find_files >&2
                    return 1
                }
                include="${2}"
                shift 2
                ;;

            -ex|--excld|--exclude)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_find_files >&2
                    return 1
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

            -ss|--shw[_-]str|--show[_-]str)
                shw_str=true
                shift 1
                ;;

            *)
                echo_err "unknown option/parameter passed: '${1}'."
                echo >&2
                help_find_files >&2
                return 1
                ;;
        esac
    done
}


function validate_args() {
    local dir_fnd_real
    local dir_pwd_real

    validate_var_dir "dir_fnd" "${dir_fnd}" 0 false

    dir_fnd_real="$(realpath "${dir_fnd}")"
    dir_pwd_real="$(realpath "${PWD}")"
    if \
           [[ "${dir_pwd_real}" == "${dir_fnd_real}" ]] \
        || [[ "${dir_pwd_real}" == "${dir_fnd_real}"/* ]]
    then
        echo_err \
            "'find_files.sh' cannot be run from the target directory being" \
            "searched, or from one of its subdirectories. This restriction" \
            "is kept as a conservative safety policy. Please run the script" \
            "from a different directory. Currently, 'dir_fnd=${dir_fnd}'" \
            "and 'PWD=${PWD}'."
        return 1
    fi

    validate_var "pattern" "${pattern}"
    if compgen -G "${pattern}" > /dev/null; then
        echo_warn \
            "the specified pattern '${pattern}' matches one or more paths" \
            "in the current working directory. This is no longer treated as" \
            "a hard error, but it can still be a sign that the script is" \
            "being run from an inconvenient location. Current 'PWD=${PWD}'."
    fi

    if [[ -n "${depth}" ]]; then
        check_int_pos "${depth}" "depth"
    fi

    check_flags_mut_excl ${chk_con} "chk_con" ${chk_exc} "chk_exc"
}


function check_tools() {
    check_pgrm_path find
    check_pgrm_path paste
    check_pgrm_path sed
    check_pgrm_path sort
    check_pgrm_path tr
}


function build_find_cmd() {
    local -a arr_build

    arr_build=(
        --arr_nam arr_cmd_find
        --dir_fnd "${dir_fnd}"
        --pattern "${pattern}"
    )

    if [[ "${follow}" == "true" ]]; then
        arr_build+=( --follow )
    fi

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
}


#TODO: make it a sourced external primitive?
function print_find_string() {
    if [[ "${fastqs}" == "true" ]]; then
        "${arr_cmd_find[@]}" \
            | sort \
            | pair_fastqs \
            | tr -d '\n' \
            | sed 's/;$//'
        echo
    else
        "${arr_cmd_find[@]}" \
            | sort \
            | tr '\n' ',' \
            | sed 's/,$//'
        echo
    fi
}


function run_find() {
    if [[ "${chk_con}" == "true" || "${chk_exc}" == "true" ]]; then
        print_banner_pretty "Call to 'find'"
        print_cmd_array arr_cmd_find
        echo
    fi

    if [[ "${chk_con}" == "true" ]]; then
        return 0
    fi

    if [[ "${chk_exc}" == "true" ]]; then
        print_banner_pretty "Results of 'find' command"
        "${arr_cmd_find[@]}" | sort
        echo

        if [[ "${shw_str}" == "true" ]]; then
            if [[ "${fastqs}" == "true" ]]; then
                print_banner_pretty \
                    "Results of 'find' command as a single semicolon- and" \
                    "comma-separated string"
            else
                print_banner_pretty \
                    "Results of 'find' command as a single comma-separated" \
                    "string"
            fi
            print_find_string
        fi

        return 0
    fi

    print_find_string
}


function main() {
    init_arg_defs
    source_helpers_script || return 1

    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        help_find_files >&2
        return 0
    fi

    parse_args "$@" || return 1
    validate_args   || return 1
    check_tools     || return 1
    build_find_cmd  || return 1
    run_find        || return 1
}


main "$@"
