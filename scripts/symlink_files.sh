#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: symlink_files.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.4, GPT-5.5) were used in development.
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

    dir_fnc="${dir_scr}/functions"
    fnc_src="${dir_fnc}/source_helpers.sh"

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
    csv_infile=""
    csv_outfile=""
    dir_out=""
    dry_run=false
    no_force=false
    quiet=false

    n_in=0
    n_out=0
    n_planned=0
    n_linked=0

    unset arr_infile arr_outfile arr_seen_out
    declare -ga arr_infile arr_outfile arr_seen_out
}


function parse_args() {
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
            -dr|--dry|--dry[_-]run)
                dry_run=true
                shift 1
                ;;

            -i|-fi|-ci|--infile|--infiles|--fil[_-]in|--csv[_-]infile|--csv[_-]infiles)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_symlink_files >&2
                    return 1
                }
                csv_infile="${2}"
                shift 2
                ;;

            -o|-co|-fo|--outfile|--outfiles|--fil[_-]out|--csv[_-]outfile|--csv[_-]outfiles)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_symlink_files >&2
                    return 1
                }
                csv_outfile="${2}"
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
    validate_var "csv_infile" "${csv_infile}"

    if [[ -n "${dir_out}" && -n "${csv_outfile}" ]]; then
        echo_err "specify exactly one of '--dir_out' or '--csv_outfile'."
        return 1
    fi

    if [[ -z "${dir_out}" && -z "${csv_outfile}" ]]; then
        echo_err "one of '--dir_out' or '--csv_outfile' must be specified."
        return 1
    fi
}


function prepare_vecs() {
    local base
    local infile

    IFS=',' read -r -a arr_infile <<< "${csv_infile}"
    n_in="${#arr_infile[@]}"

    if [[ -n "${dir_out}" ]]; then
        for infile in "${arr_infile[@]}"; do
            base="$(basename "${infile}")"
            arr_outfile+=( "${dir_out}/${base}" )
        done
    else
        IFS=',' read -r -a arr_outfile <<< "${csv_outfile}"
        n_out="${#arr_outfile[@]}"
    fi
}


function validate_vecs() {
    local i
    local infile
    local outfile
    local parent
    local seen

    if [[ "${n_in}" -eq 0 ]]; then
        echo_err "'--csv_infile' resolved to an empty vector."
        return 1
    fi

    for infile in "${arr_infile[@]}"; do
        if [[ -z "${infile}" ]]; then
            echo_err "'--csv_infile' contains an empty element."
            return 1
        fi

        if [[ "${infile}" == *","* ]]; then
            echo_err "input paths must not contain commas: '${infile}'."
            return 1
        fi

        validate_var_file "infile" "${infile}"
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
                "'--csv_outfile' must contain the same number of elements" \
                "as '--csv_infile' (${n_out} vs ${n_in})."
            return 1
        fi

        for outfile in "${arr_outfile[@]}"; do
            if [[ -z "${outfile}" ]]; then
                echo_err "'--csv_outfile' contains an empty element."
                return 1
            fi

            if [[ "${outfile}" == *","* ]]; then
                echo_err "output paths must not contain commas: '${outfile}'."
                return 1
            fi

            parent="$(dirname "${outfile}")"
            validate_var_dir "parent" "${parent}"
        done
    fi

    #  Reject duplicate resolved output paths, including basename collisions
    #+ under '--dir_out'
    for (( i = 0; i < n_in; i++ )); do
        outfile="${arr_outfile[${i}]}"

        for seen in "${arr_seen_out[@]}"; do
            if [[ "${outfile}" == "${seen}" ]]; then
                echo_err \
                    "resolved output path appears more than once:" \
                    "'${outfile}'."
                return 1
            fi
        done
        arr_seen_out+=( "${outfile}" )

        if [[ -e "${outfile}" && ! -L "${outfile}" ]]; then
            echo_err \
                "output path already exists and is not a symlink:" \
                "'${outfile}'."
            return 1
        fi

        if [[ -L "${outfile}" && "${no_force}" == "true" ]]; then
            echo_err \
                "output symlink already exists and '--no_force' was" \
                "specified: '${outfile}'."
            return 1
        fi
    done
}


function run_jobs() {
    local i
    local infile
    local outfile

    for (( i = 0; i < n_in; i++ )); do
        infile="${arr_infile[${i}]}"
        outfile="${arr_outfile[${i}]}"

        if [[ "${dry_run}" == "true" ]]; then
            if [[ "${no_force}" == "true" ]]; then
                printf "Dry run: ln -s \\\n" >&2
            else
                printf "Dry run: ln -sf \\\n" >&2
            fi
            printf "    %q \\\n" "${infile}" >&2
            printf "    %q\n" "${outfile}" >&2
            n_planned=$(( n_planned + 1 ))
            continue
        fi

        if [[ "${no_force}" == "true" ]]; then
            ln -s "${infile}" "${outfile}"
        else
            ln -sf "${infile}" "${outfile}"
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
