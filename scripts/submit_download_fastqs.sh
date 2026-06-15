#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: submit_download_fastqs.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models) were used in
# development.
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


#  Print main usage text
function show_help_main() {
    cat << EOM
Usage:
  submit_download_fastqs.sh
    srr url_1 url_2 dir_out dir_sym nam_cus dir_eo nam_job

Description:
  Download one single-end or paired-end FASTQ entry and create custom symlink(s).

Positional arguments:
  01  srr  <str>
    NCBI SRA database run accession code.

  02  url_1  <str>
    URL (FTP or HTTPS) for FASTQ file.

  03  url_2  <str>
    Second FASTQ URL for paired-end data ("NA" for single-end data).

  04  dir_out  <dir>
    Directory to save FASTQ file(s).

  05  dir_sym  <dir>
    Directory for symlink(s) to FASTQ file(s).

  06  nam_cus  <str>
    Custom name for symlink(s).

  07  dir_eo  <dir>
    Directory for stderr and stdout files.

  08  nam_job  <str>
    Job name.

Notes:
  - Use 'NA' for 'url_2' when processing single-end data.
EOM
}


#  Initialize argument variables
function init_defs() {
    srr=""
    url_1=""
    url_2=""
    dir_out=""
    dir_sym=""
    nam_cus=""
    dir_eo=""
    nam_job=""
    dir_scr=""
    dir_fnc=""
}


#  Source shared helpers after help-only handling
function source_submit_helpers() {
    local fnc_src

    dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
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
        check_env \
        check_inputs \
        format_outputs \
        || {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to source required helper scripts." >&2
            return 1
        }
}


#  Parse the fixed positional submit interface
function parse_args() {
    local msg

    if [[ $# -ne 8 ]]; then
        msg="but $# were supplied."
        [[ $# -eq 1 ]] && msg="but only $# was supplied."
        cat >&2 << EOM
error: '$(basename "${0}")' requires 8 positional arguments, ${msg}

The necessary positional arguments:
$(show_help_main)

EOM
        return 1
    fi

    srr="${1}"
    url_1="${2}"
    url_2="${3}"
    dir_out="${4}"
    dir_sym="${5}"
    nam_cus="${6}"
    dir_eo="${7}"
    nam_job="${8}"
}


#  Validate submit-layer directories and runtime tool availability
function validate_args() {
    local dir

    for dir in "${dir_out}" "${dir_sym}" "${dir_eo}"; do
        if [[ ! -d "${dir}" ]]; then
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "directory does not exist: '${dir}'" >&2
            return 1
        fi
    done

    check_pgrm_path wget || return 1
}


#  Download one FASTQ URL to one output file
function download_fastq() {
    local url="${1:-}"
    local outfile="${2:-}"
    local log_out="${3:-}"
    local log_err="${4:-}"

    echo "Downloading ${srr} from ${url}." >&2
    wget --progress=dot:mega -O "${outfile}" "${url}" \
         > "${log_out}" \
        2> "${log_err}" ||
        {
            echo_err "failed to download '${url}'."
            return 1
        }
}


#  Create one custom-name symlink
function link_fastq() {
    local src="${1:-}"
    local dst="${2:-}"

    ln -sf "${src}" "${dst}" ||
        {
            echo_err "failed to create symlink for '${src}'"
            return 1
        }
}


#  Download the selected SE or PE entry and create custom symlink(s)
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


#  Main script execution
function main() {
    init_defs

    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        cat >&2 << EOM
'$(basename "${0}")' requires 8 positional arguments:

The necessary positional arguments:
$(show_help_main)

EOM
        return 0
    fi

    source_submit_helpers || return 1
    parse_args "$@"       || return 1
    validate_args         || return 1
    run_downloads         || return 1
}


main "$@"
