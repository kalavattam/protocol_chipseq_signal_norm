#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: symlink_files.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.4, GPT-5.5, GPT-5.6) were used in design,
# development, and documentation, with all output reviewed, edited, and
# approved by the author.
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

    source_helpers "${dir_fnc}" \
        check_args \
        check_inputs \
        format_outputs \
        help/help_symlink_files \
        || {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to source required helper scripts." >&2
            return 1
        }
}


function init_arg_defs() {
    csv_fil_in=""
    csv_fil_out=""
    dir_out=""
    dry_run=false
    no_force=false
    quiet=false

    n_in=0
    n_out=0
    n_planned=0
    n_linked=0

    unset arr_fil_in arr_fil_out arr_seen_out
    declare -ga arr_fil_in arr_fil_out arr_seen_out
}


function parse_args() {
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
            -dr|--dry|--dry[_-]run)
                dry_run=true
                shift 1
                ;;

            -ci|--csv[_-]fil_in)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_symlink_files >&2
                    return 1
                }
                csv_fil_in="${2}"
                shift 2
                ;;

            -co|--csv[_-]fil_out)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_symlink_files >&2
                    return 1
                }
                csv_fil_out="${2}"
                shift 2
                ;;

            -do|--dir[_-]out)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_symlink_files >&2
                    return 1
                }
                dir_out="${2}"
                shift 2
                ;;

            -nf|--no[_-]force)
                no_force=true
                shift 1
                ;;

            -q|--quiet)
                quiet=true
                shift 1
                ;;

            *)
                echo_err "unknown option/parameter passed: '${1}'."
                echo >&2
                help_symlink_files >&2
                return 1
                ;;
        esac
    done
}


function validate_args() {
    validate_var "csv_fil_in" "${csv_fil_in}"

    if [[ -n "${dir_out}" && -n "${csv_fil_out}" ]]; then
        echo_err "specify exactly one of '--dir_out' or '--csv_fil_out'."
        return 1
    fi

    if [[ -z "${dir_out}" && -z "${csv_fil_out}" ]]; then
        echo_err "one of '--dir_out' or '--csv_fil_out' must be specified."
        return 1
    fi
}


function prepare_vecs() {
    local base
    local fil_in

    IFS=',' read -r -a arr_fil_in <<< "${csv_fil_in}"
    n_in="${#arr_fil_in[@]}"

    if [[ -n "${dir_out}" ]]; then
        for fil_in in "${arr_fil_in[@]}"; do
            base="$(basename "${fil_in}")"
            arr_fil_out+=( "${dir_out}/${base}" )
        done
    else
        IFS=',' read -r -a arr_fil_out <<< "${csv_fil_out}"
        n_out="${#arr_fil_out[@]}"
    fi
}


function validate_vecs() {
    local i
    local fil_in
    local fil_out
    local parent
    local seen

    if [[ "${n_in}" -eq 0 ]]; then
        echo_err "'--csv_fil_in' resolved to an empty vector."
        return 1
    fi

    for fil_in in "${arr_fil_in[@]}"; do
        if [[ -z "${fil_in}" ]]; then
            echo_err "'--csv_fil_in' contains an empty element."
            return 1
        fi

        if [[ "${fil_in}" == *","* ]]; then
            echo_err "input paths must not contain commas: '${fil_in}'."
            return 1
        fi

        validate_var_file "fil_in" "${fil_in}"
    done

    if [[ -n "${dir_out}" ]]; then
        if [[ "${dir_out}" == *","* ]]; then
            echo_err "'--dir_out' must not contain commas: '${dir_out}'."
            return 1
        fi

        validate_var_dir "dir_out" "${dir_out}"
    else
        if [[ "${n_out}" -ne "${n_in}" ]]; then
            echo_err \
                "'--csv_fil_out' must contain the same number of elements" \
                "as '--csv_fil_in' (${n_out} vs ${n_in})."
            return 1
        fi

        for fil_out in "${arr_fil_out[@]}"; do
            if [[ -z "${fil_out}" ]]; then
                echo_err "'--csv_fil_out' contains an empty element."
                return 1
            fi

            if [[ "${fil_out}" == *","* ]]; then
                echo_err "output paths must not contain commas: '${fil_out}'."
                return 1
            fi

            parent="$(dirname "${fil_out}")"
            validate_var_dir "parent" "${parent}"
        done
    fi

    #  Reject duplicate resolved output paths, including basename collisions
    #+ under '--dir_out'
    for (( i = 0; i < n_in; i++ )); do
        fil_out="${arr_fil_out[${i}]}"

        for seen in "${arr_seen_out[@]}"; do
            if [[ "${fil_out}" == "${seen}" ]]; then
                echo_err \
                    "resolved output path appears more than once:" \
                    "'${fil_out}'."
                return 1
            fi
        done
        arr_seen_out+=( "${fil_out}" )

        if [[ -e "${fil_out}" && ! -L "${fil_out}" ]]; then
            echo_err \
                "output path already exists and is not a symlink:" \
                "'${fil_out}'."
            return 1
        fi

        if [[ -L "${fil_out}" && "${no_force}" == "true" ]]; then
            echo_err \
                "output symlink already exists and '--no_force' was" \
                "specified: '${fil_out}'."
            return 1
        fi
    done
}


function run_jobs() {
    local i
    local fil_in
    local fil_out

    for (( i = 0; i < n_in; i++ )); do
        fil_in="${arr_fil_in[${i}]}"
        fil_out="${arr_fil_out[${i}]}"

        if [[ "${dry_run}" == "true" ]]; then
            if [[ "${no_force}" == "true" ]]; then
                printf "Dry run: ln -s \\\n" >&2
            else
                printf "Dry run: ln -sf \\\n" >&2
            fi
            printf "    %q \\\n" "${fil_in}" >&2
            printf "    %q\n" "${fil_out}" >&2
            n_planned=$(( n_planned + 1 ))
            continue
        fi

        if [[ "${no_force}" == "true" ]]; then
            ln -s "${fil_in}" "${fil_out}"
        else
            ln -sf "${fil_in}" "${fil_out}"
        fi
        n_linked=$(( n_linked + 1 ))
    done

    if [[ "${quiet}" == "false" ]]; then
        printf "\nsummary\t"                  >&2
        printf "inputs=%d\t"   "${n_in}"      >&2
        printf "linked=%d\t"   "${n_linked}"  >&2
        printf "planned=%d\t"  "${n_planned}" >&2
        printf "dry_run=%s\t"  "${dry_run}"   >&2
        printf "no_force=%s\n" "${no_force}"  >&2
    fi
}


function main() {
    init_arg_defs
    source_helpers_script || return 1

    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        help_symlink_files >&2
        return 0
    fi

    parse_args "$@" || return 1
    validate_args   || return 1
    prepare_vecs    || return 1
    validate_vecs   || return 1
    run_jobs        || return 1
}


main "$@"
