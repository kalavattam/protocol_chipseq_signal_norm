#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: submit_download_fastqs.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
#   GPT-5.6);
# - Anthropic Claude Code (Opus 5).
#
# Distributed under the MIT license.


# TODO: switch from positional to keyword option parameters, consistent with
# all other wrappers.
# Require bash >= 4.4 before doing any work.
if [[ -z "${BASH_VERSION:-}" ]]; then
    echo "error(shell):" \
        "this script must be run under bash >= 4.4." >&2
    exit 1
elif ((
    BASH_VERSINFO[0] < 4 || ( BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] < 4 )
)); then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "this script requires bash >= 4.4; current version is" \
        "'${BASH_VERSION}'." >&2
    exit 1
fi

# Run in safe mode, exiting on errors, unset variables, and pipe failures.
set -euo pipefail

# 'sbatch' copies this script, so '--dir_scr' must outrank 'BASH_SOURCE'.
function resolve_dir_scr() {
    local -a args=( "$@" )
    local i

    for (( i = 0; i < ${#args[@]}; i++ )); do
        case "${args[i]}" in
            -ds|--dir[_-]scr)
                if \
                       (( i + 1 < ${#args[@]} )) \
                    && [[ -n "${args[i + 1]}" ]] \
                    && [[ "${args[i + 1]}" != -* ]]
                then
                    printf '%s\n' "${args[i + 1]}"
                    return 0
                fi
                ;;
        esac
    done

    cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
}


dir_scr="$(resolve_dir_scr "$@")" || {
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "failed to resolve the script directory." >&2
    exit 1
}

fil_hlp="${dir_scr}/../lib/bash/help/help_submit_download_fastqs.sh"

# shellcheck source=lib/bash/help/help_submit_download_fastqs.sh
source "${fil_hlp}" || {
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "failed to source main help." >&2
    exit 1
}
unset fil_hlp


# Initialize argument variables.
function init_defs() {
    hlp_req="false"
    srr=""
    url_1=""
    url_2=""
    dir_out=""
    dir_sym=""
    nam_cus=""
    dir_eo=""
    nam_job=""
    dir_fnc=""
}


# Detect terminal help before resolving bootstrap arguments.
function scan_help_args() {
    local arg

    for arg in "$@"; do
        case "${arg}" in
            -h|--hlp|--help)
                hlp_req="true"
                return 0
                ;;
        esac
    done
}


# Source shared helpers after bootstrap argument parsing.
function source_helpers_submit() {
    local scr="${1:-}"
    local dir_scr_arg="${2:-}"
    local fnc_src

    if (( $# < 2 )); then
        echo "error(${scr:-unknown_script}):" \
            "expected at least two arguments: 'scr' and 'dir_scr_arg'." >&2
        return 1
    fi

    if [[ -z "${scr}" ]]; then
        scr="unknown_script"
    fi

    if [[ -z "${dir_scr_arg}" ]]; then
        echo "error(${scr}):" \
            "positional argument 2, 'dir_scr_arg', is missing." >&2
        return 1
    elif [[ ! -d "${dir_scr_arg}" ]]; then
        echo "error(${scr}):" \
            "script directory not found: '${dir_scr_arg}'." >&2
        return 1
    fi

    dir_scr="${dir_scr_arg}"
    dir_fnc="${dir_scr_arg}/../lib/bash"
    fnc_src="${dir_fnc}/core/source_helpers.sh"

    if [[ ! -f "${fnc_src}" ]]; then
        echo "error(${scr}): script not found: '${fnc_src}'." >&2
        return 1
    fi

    # shellcheck disable=SC1090
    source "${fnc_src}" || {
        echo "error(${scr}): failed to source '${fnc_src}'." >&2
        return 1
    }

    source_helpers "${dir_fnc}" \
        check_args \
        check_env \
        format_outputs \
        || {
            echo "error(${scr}): failed to source required helper scripts." >&2
            return 1
        }
}


# Parse the bootstrap option and fixed positional submit interface.
function parse_args() {
    local msg
    local -a args_pos=()

    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
            -h|--hlp|--help)
                hlp_req="true"
                return 0
                ;;

            -ds|--dir[_-]scr)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_download_fastqs >&2
                    return 1
                }
                dir_scr="${2}"
                shift 2
                ;;

            -*)
                echo_err "unknown option/parameter passed: '${1}'."
                echo >&2
                help_submit_download_fastqs >&2
                return 1
                ;;

            *)
                args_pos+=( "${1}" )
                shift 1
                ;;
        esac
    done

    if [[ ${#args_pos[@]} -ne 8 ]]; then
        msg="but ${#args_pos[@]} were supplied."

        if [[ ${#args_pos[@]} -eq 1 ]]; then
            msg="but only ${#args_pos[@]} was supplied."
        fi

        cat >&2 << EOM
error: '$(basename "${0}")' requires 8 positional arguments, ${msg}

$(help_submit_download_fastqs)
EOM
        return 1
    fi

    srr="${args_pos[0]}"
    url_1="${args_pos[1]}"
    url_2="${args_pos[2]}"
    dir_out="${args_pos[3]}"
    dir_sym="${args_pos[4]}"
    nam_cus="${args_pos[5]}"
    dir_eo="${args_pos[6]}"
    nam_job="${args_pos[7]}"
}


# Validate submit-layer directories and runtime tool availability.
function validate_args() {
    local dir

    for dir in dir_out dir_sym dir_eo; do
        local value="${!dir}"

        if [[ ! -d "${value}" ]]; then
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "directory does not exist: '${value}'" >&2
            return 1
        fi

        value="$(cd "${value}" > /dev/null 2>&1 && pwd -P)" || {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to resolve directory: '${value}'" >&2
            return 1
        }
        printf -v "${dir}" '%s' "${value}"
    done

    check_pgrm_path ln   || return 1
    check_pgrm_path wget || return 1
}


# Download one FASTQ URL to one output file.
function download_fastq() {
    local url="${1:-}"
    local fil_out="${2:-}"
    local log_out="${3:-}"
    local log_err="${4:-}"

    echo "Downloading ${srr} from ${url}." >&2
    wget --progress=dot:mega -O "${fil_out}" "${url}" \
        > "${log_out}" \
        2> "${log_err}" ||
        {
            echo_err "failed to download '${url}'."
            return 1
        }
}


# Create one custom-name symlink.
function link_fastq() {
    local src="${1:-}"
    local dst="${2:-}"

    ln -sf "${src}" "${dst}" ||
        {
            echo_err "failed to create symlink for '${src}'"
            return 1
        }
}


# Download the selected SE or PE entry and create custom symlink(s).
function run_downloads() {
    if [[ "${url_2}" != "NA" ]]; then
        download_fastq \
            "${url_1}" \
            "${dir_out}/${srr}_R1.fastq.gz" \
            "${dir_eo}/${nam_job}.${srr}_R1.stdout.txt" \
            "${dir_eo}/${nam_job}.${srr}_R1.stderr.txt" || return 1

        download_fastq \
            "${url_2}" \
            "${dir_out}/${srr}_R2.fastq.gz" \
            "${dir_eo}/${nam_job}.${srr}_R2.stdout.txt" \
            "${dir_eo}/${nam_job}.${srr}_R2.stderr.txt" || return 1
    else
        download_fastq \
            "${url_1}" \
            "${dir_out}/${srr}.fastq.gz" \
            "${dir_eo}/${nam_job}.${srr}.stdout.txt" \
            "${dir_eo}/${nam_job}.${srr}.stderr.txt" || return 1
    fi

    echo "Symlinking ${srr} to ${nam_cus}." >&2

    if [[ "${url_2}" != "NA" ]]; then
        link_fastq \
            "${dir_out}/${srr}_R1.fastq.gz" \
            "${dir_sym}/${nam_cus}_R1.fastq.gz" || return 1

        link_fastq \
            "${dir_out}/${srr}_R2.fastq.gz" \
            "${dir_sym}/${nam_cus}_R2.fastq.gz" || return 1
    else
        link_fastq \
            "${dir_out}/${srr}.fastq.gz" \
            "${dir_sym}/${nam_cus}.fastq.gz" || return 1
    fi
}


# Main script execution.
function main() {
    init_defs
    scan_help_args "$@"

    if [[ "${hlp_req}" == "true" ]]; then
        help_submit_download_fastqs
        return 0
    fi

    if [[ -z "${1:-}" ]]; then
        help_submit_download_fastqs
        return 0
    fi

    source_helpers_submit "${0##*/}" "${dir_scr}" || return 1
    parse_args "$@"                               || return 1

    if [[ "${hlp_req}" == "true" ]]; then
        help_submit_download_fastqs
        return 0
    fi

    validate_args                                 || return 1
    run_downloads                                 || return 1
}


main "$@"
