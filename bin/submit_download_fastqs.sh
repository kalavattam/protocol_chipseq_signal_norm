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
# - Anthropic Claude Code (Opus 5, Opus 5.5).
#
# Distributed under the MIT license.


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

# 'sbatch' runs a copy, so 'BASH_SOURCE' can resolve outside the repo; if so,
# fail here rather than at a subsequent 'source'.
if [[ ! -d "${dir_scr}/../lib/bash" ]]; then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "cannot locate helper libraries from '${dir_scr}'; pass '--dir_scr'" \
        "when this script is run from a copy." >&2
    exit 1
fi

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
    url_2="NA"
    md5_1="NA"
    md5_2="NA"
    siz_1="NA"
    siz_2="NA"
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
        check_inputs \
        format_outputs \
        || {
            echo "error(${scr}): failed to source required helper scripts." >&2
            return 1
        }
}


# Parse the keyword submit interface.
function parse_args() {
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

            -sr|--srr)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_download_fastqs >&2
                    return 1
                }
                srr="${2}"
                shift 2
                ;;

            -u1|--url[_-]1)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_download_fastqs >&2
                    return 1
                }
                url_1="${2}"
                shift 2
                ;;

            -u2|--url[_-]2)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_download_fastqs >&2
                    return 1
                }
                url_2="${2}"
                shift 2
                ;;

            -m1|--md5[_-]1)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_download_fastqs >&2
                    return 1
                }
                md5_1="${2}"
                shift 2
                ;;

            -m2|--md5[_-]2)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_download_fastqs >&2
                    return 1
                }
                md5_2="${2}"
                shift 2
                ;;

            -s1|--siz[_-]1)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_download_fastqs >&2
                    return 1
                }
                siz_1="${2}"
                shift 2
                ;;

            -s2|--siz[_-]2)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_download_fastqs >&2
                    return 1
                }
                siz_2="${2}"
                shift 2
                ;;

            -do|--dir[_-]out)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_download_fastqs >&2
                    return 1
                }
                dir_out="${2}"
                shift 2
                ;;

            -dy|--dir[_-]sym|--dir[_-]symlink)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_download_fastqs >&2
                    return 1
                }
                dir_sym="${2}"
                shift 2
                ;;

            -nc|--nam[_-]cus)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_download_fastqs >&2
                    return 1
                }
                nam_cus="${2}"
                shift 2
                ;;

            -deo|--dir[_-]eo)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_download_fastqs >&2
                    return 1
                }
                dir_eo="${2}"
                shift 2
                ;;

            -nj|--nam[_-]job)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_submit_download_fastqs >&2
                    return 1
                }
                nam_job="${2}"
                shift 2
                ;;

            *)
                echo_err "unknown option/parameter passed: '${1}'."
                echo >&2
                help_submit_download_fastqs >&2
                return 1
                ;;
        esac
    done
}


# Validate submit-layer directories and runtime tool availability.
function validate_args() {
    local dir nam val

    validate_var "srr"     "${srr}"     || return 1
    validate_var "url_1"   "${url_1}"   || return 1
    validate_var "url_2"   "${url_2}"   || return 1
    validate_var "md5_1"   "${md5_1}"   || return 1
    validate_var "md5_2"   "${md5_2}"   || return 1
    validate_var "siz_1"   "${siz_1}"   || return 1
    validate_var "siz_2"   "${siz_2}"   || return 1
    validate_var "dir_out" "${dir_out}" || return 1
    validate_var "dir_sym" "${dir_sym}" || return 1
    validate_var "nam_cus" "${nam_cus}" || return 1
    validate_var "dir_eo"  "${dir_eo}"  || return 1
    validate_var "nam_job" "${nam_job}" || return 1

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

    for nam in md5_1 md5_2; do
        val="${!nam}"

        if [[ "${val}" != "NA" && ! "${val}" =~ ^[0-9A-Fa-f]{32}$ ]]; then
            echo_err \
                "'--${nam}' must be 32 hexadecimal characters or 'NA':" \
                "'${val}'."
            return 1
        fi
    done

    for nam in siz_1 siz_2; do
        val="${!nam}"

        if [[ "${val}" != "NA" && ! "${val}" =~ ^[1-9][0-9]*$ ]]; then
            echo_err \
                "'--${nam}' must be a positive integer number of bytes or" \
                "'NA': '${val}'."
            return 1
        fi
    done

    check_pgrm_path gzip || return 1
    check_pgrm_path ln   || return 1
    check_pgrm_path wget || return 1

    if [[ "${md5_1}" != "NA" || "${md5_2}" != "NA" ]]; then
        check_pgrm_path md5sum || return 1
    fi
}


# Check a downloaded file against its expected size and MD5; without an MD5,
# test its gzip integrity. Print the reason for a failure.
function verify_fastq() {
    local fil="${1:-}"
    local md5="${2:-}"
    local siz="${3:-}"
    local val

    if [[ "${siz}" != "NA" ]]; then
        val="$(wc -c < "${fil}")"
        val="${val//[[:space:]]/}"

        if [[ "${val}" != "${siz}" ]]; then
            echo "size is ${val} bytes, expected ${siz}"
            return 1
        fi
    fi

    if [[ "${md5}" != "NA" ]]; then
        val="$(md5sum < "${fil}")"
        val="${val%% *}"

        if [[ "${val,,}" != "${md5,,}" ]]; then
            echo "MD5 is ${val}, expected ${md5,,}"
            return 1
        fi
    elif ! \
        gzip -t < "${fil}" 2> /dev/null
    then
        echo "file fails 'gzip -t'"
        return 1
    fi
}


# Download one FASTQ URL to a partial file, verify it, and move it into place;
# a file that fails verification is downloaded again, up to 'n_try' times.
function download_fastq() {
    local url="${1:-}"
    local fil_out="${2:-}"
    local md5="${3:-}"
    local siz="${4:-}"
    local log_out="${5:-}"
    local log_err="${6:-}"
    local fil_prt="${fil_out}.part"
    local n_try=3
    local i_try reason
    local -a cmd_wget

    # ENA's 'fastq_ftp' values carry no scheme, and 'wget' would then use HTTP
    # on port 80.
    if [[ "${url}" != *://* ]]; then
        url="https://${url}"
    fi

    cmd_wget=(
        wget
            --progress=dot:mega
            --tries=5
            --retry-connrefused
            --waitretry=10
            -O "${fil_prt}"
            "${url}"
    )

    for (( i_try = 1; i_try <= n_try; i_try++ )); do
        echo \
            "Downloading ${srr} from ${url} (attempt ${i_try} of" \
            "${n_try})." >&2
        printf '%q ' "${cmd_wget[@]}" >&2
        echo >&2

        if ! \
            "${cmd_wget[@]}" > "${log_out}" 2> "${log_err}"
        then
            rm -f -- "${fil_prt}"
            echo_err "failed to download '${url}'; see '${log_err}'."
            return 1
        fi

        if \
            reason="$(verify_fastq "${fil_prt}" "${md5}" "${siz}")"
        then
            mv -f -- "${fil_prt}" "${fil_out}" || {
                rm -f -- "${fil_prt}"
                echo_err "failed to move '${fil_prt}' to '${fil_out}'."
                return 1
            }
            return 0
        fi

        if (( i_try < n_try )); then
            echo_warn \
                "download of '${url}' failed verification: ${reason};" \
                "downloading it again."
        fi
    done

    rm -f -- "${fil_prt}"
    echo_err \
        "download of '${url}' failed verification ${n_try} times: ${reason}."
    return 1
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
            "${md5_1}" \
            "${siz_1}" \
            "${dir_eo}/${nam_job}.${srr}_R1.stdout.txt" \
            "${dir_eo}/${nam_job}.${srr}_R1.stderr.txt" || return 1

        download_fastq \
            "${url_2}" \
            "${dir_out}/${srr}_R2.fastq.gz" \
            "${md5_2}" \
            "${siz_2}" \
            "${dir_eo}/${nam_job}.${srr}_R2.stdout.txt" \
            "${dir_eo}/${nam_job}.${srr}_R2.stderr.txt" || return 1
    else
        download_fastq \
            "${url_1}" \
            "${dir_out}/${srr}.fastq.gz" \
            "${md5_1}" \
            "${siz_1}" \
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

    validate_args || return 1
    run_downloads || return 1
}


main "$@"
