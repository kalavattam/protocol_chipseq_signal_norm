#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: construct_find.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
#
# Distributed under the MIT license.


# construct_find


#  Require Bash >= 4.4 before defining functions
if [[ -z "${BASH_VERSION:-}" ]]; then
    echo "error(shell):" \
        "this script must be sourced or run under Bash >= 4.4." >&2

    if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
        return 1
    else
        exit 1
    fi
elif ((
    BASH_VERSINFO[0] < 4 || ( BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] < 4 )
)); then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "this script requires Bash >= 4.4; current version is" \
        "'${BASH_VERSION}'." >&2

    if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
        return 1
    else
        exit 1
    fi
fi

#  Source required helper functions if needed
# shellcheck disable=SC1091
{
    _dir_src_find="$(
        cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
    )"

    source "${_dir_src_find}/source_helpers.sh" || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source '${_dir_src_find}/source_helpers.sh'." >&2

        if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
            return 1
        else
            exit 1
        fi
    }

    source_helpers "${_dir_src_find}" \
        check_args check_source format_outputs || {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to source required helper dependencies." >&2

            if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
                return 1
            else
                exit 1
            fi
        }

    unset _dir_src_find
}


function construct_find() {
    local dir_fnd=""
    local pattern=""
    local follow=false
    local depth=""
    local include=""
    local exclude=""
    local arr_nam=""
    local show_help
    local patt
    local -a arr_inc
    local -a arr_exc

    show_help=$(cat << EOM
Usage:
  construct_find
    [--help] --arr_nam <str> --dir_fnd <str> --pattern <str> [--follow] [--depth <int>] [--include <str,str,...>] [--exclude <str,str,...>]

Description:
  Build a 'find' command as an argument array for safe downstream execution.

Keyword arguments:
  -an, --arr_nam, --arr-nam      <str>  Name of the output array variable to populate.
  -df, --dir_fnd, --dir-fnd      <str>  Directory in which to search for files (required).
  -pa, --pattern                 <str>  Primary filename pattern to match (required).
  -fl, -sy, --follow, --symlink  <flg>  Follow symbolic links during the search.
  -de, --dpth, --depth           <int>  Maximum search depth.
  -in, --incld, --include        <str>  Comma-separated list of additional filename patterns to include; each is applied as an additional '-name' predicate, i.e., logical AND.
  -ex, --excld, --exclude        <str>  Comma-separated list of filename patterns to exclude; each is applied with logical AND NOT.

Returns:
  Populates the named array variable with a 'find' argument array.

Notes:
  - This function builds an argument array and is intended to be consumed as "\${arr_cmd[@]}" rather than by reconstructing a command string.
  - Logical OR is not implemented; '--include' patterns are applied as repeated '-name' predicates, i.e., logical AND.
EOM
    )

    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    while (( $# > 0 )); do
        case "${1}" in
            -an|--arr[_-]nam)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                arr_nam="${2:-}"
                shift 2
                ;;

            -df|--dir[_-]fnd)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                dir_fnd="${2:-}"
                shift 2
                ;;

            -pa|--pattern)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                pattern="${2:-}"
                shift 2
                ;;

            -fl|-sy|--follow|--symlink)
                follow=true
                shift 1
                ;;

            -de|--dpth|--depth)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                depth="${2:-}"
                shift 2
                ;;

            -in|--incld|--include)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                include="${2:-}"
                shift 2
                ;;

            -ex|--excld|--exclude)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                exclude="${2:-}"
                shift 2
                ;;

            *)
                echo "## Unknown argument passed: '${1}' ##" >&2
                echo >&2
                echo "${show_help}" >&2
                return 1
                ;;
        esac
    done

    if [[ -z "${arr_nam}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'--arr_nam' is required."
        return 1
    elif [[ ! "${arr_nam}" =~ ^[a-zA-Z_][a-zA-Z0-9_]*$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'--arr_nam' must be a valid shell variable name: '${arr_nam}'."
        return 1
    fi

    if [[ -z "${dir_fnd}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'--dir_fnd' is required."
        return 1
    fi

    if [[ ! -d "${dir_fnd}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "directory associated with '--dir_fnd' does not exist:" \
            "'${dir_fnd}'."
        return 1
    fi

    if [[ ! -r "${dir_fnd}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "directory associated with '--dir_fnd' is not readable:" \
            "'${dir_fnd}'."
        return 1
    fi

    if [[ -z "${pattern}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'--pattern' is required."
        return 1
    fi

    if [[ -n "${depth}" && ! "${depth}" =~ ^[1-9][0-9]*$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'--depth' must be a positive integer: '${depth}'."
        return 1
    fi

    local -n arr_out="${arr_nam}"
    arr_out=()

    arr_out+=( find )

    if [[ "${follow}" == "true" ]]; then
        arr_out+=( -L )
    fi

    arr_out+=( "${dir_fnd}" )

    if [[ -n "${depth}" ]]; then
        arr_out+=( -maxdepth "${depth}" )
    fi

    arr_out+=( -type f -name "${pattern}" )

    if [[ -n "${include}" ]]; then
        IFS=',' read -r -a arr_inc <<< "${include}"
        for patt in "${arr_inc[@]}"; do
            [[ -n "${patt}" ]] || continue
            arr_out+=( -name "${patt}" )
        done
    fi

    if [[ -n "${exclude}" ]]; then
        IFS=',' read -r -a arr_exc <<< "${exclude}"
        for patt in "${arr_exc[@]}"; do
            [[ -n "${patt}" ]] || continue
            arr_out+=( ! -name "${patt}" )
        done
    fi
}


#  Print an error message when function script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    err_source_only "${BASH_SOURCE[0]}"
fi
