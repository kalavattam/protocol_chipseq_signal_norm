#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: check_unity.sh
#
# Copyright 2025-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
# GPT-5.6) were used in design, development, and documentation, with all output
# reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


# check_unity


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
{
    _dir_src_unity="$(
        cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
    )"

    # shellcheck source=lib/bash/core/source_helpers.sh
    source "${_dir_src_unity}/source_helpers.sh" || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source '${_dir_src_unity}/source_helpers.sh'." >&2

        if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
            return 1
        else
            exit 1
        fi
    }

    source_helpers "${_dir_src_unity}" \
        check_args check_numbers check_source format_outputs || {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to source required helper dependencies." >&2

            if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
                return 1
            else
                exit 1
            fi
        }

    unset _dir_src_unity
}


function check_unity() {
    local fil_in="${1}"            # bedGraph fil_in ('.gz' OK)
    local bnd_gt="${2:-0.999998}"  # Lower bound for unity check
    local bnd_lt="${3:-1.000002}"  # Upper bound for unity check
    local quiet="${4:-false}"      # Boolean-like: suppress success message
    local awk_unity                # AWK program for unity check
    local status=0                 # Exit status of AWK or gzip|AWK pipeline
    local pipefail_set=0           # Flag: caller had 'pipefail' enabled
    local show_help                # Help message/documentation

    show_help=$(cat << EOM
Usage
-----
  check_unity
    [--help] fil_in [bnd_gt] [bnd_lt] [quiet]

  Check whether the sum of column 4 in a bedGraph file is approximately 1 (unity), which is expected for "normalized coverage" (PMID 37160995). The input file may be plain text or gzipped.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  fil_in : file
    Input file path. Path to input bedGraph file ('.gz' OK).

  2  bnd_gt : float
    Lower bound for unity check (default: ${bnd_gt}).

  3  bnd_lt : float
    Upper bound for unity check (default: ${bnd_lt}).

  4  quiet : bool
    Suppress success message if true-like; accepted values are case-insensitive Boolean-like strings (default: ${quiet}).

Returns
-------
  0 if the sum is within the specified bounds; otherwise, 1 and an error message.

Notes
-----
  Runtime requirements:
    - awk
    - bash >= 4.4
    - gunzip (when processing a gzipped bedGraph)

  Accepted true-like values: true, t, yes, y, 1. Accepted false-like values: false, f, no, n, 0. Matching is case-insensitive; surrounding whitespace and empty values are invalid. Successful normalization emits 'true' or 'false'.

Examples
--------
  1. Check a normalized bedGraph with the default unity bounds.
    '''bash
    check_unity "coverage_track.bdg"
    '''

  2. Check a gzipped track with wider bounds and quiet success output.
    '''bash
    check_unity "coverage_track.bdg.gz" 0.99 1.01 "true"
    '''
EOM
    )

    #  Display help message if help flag or no arguments is given (latter with
    #+ error)
    if [[ "${fil_in}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${fil_in}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'fil_in', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ "${fil_in}" == -* ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "unknown option '${fil_in}'. This function accepts only '-h'," \
            "'--hlp', or '--help' as options; otherwise supply 'fil_in' as" \
            "positional argument 1."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  Validate input file
    if [[ ! -f "${fil_in}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "input file not found or not accessible: '${fil_in}'."
        return 1
    fi

    #  Validate bounds
    check_flt_pos "${bnd_gt}" "bnd_gt" || return 1
    check_flt_pos "${bnd_lt}" "bnd_lt" || return 1

    if ! \
        awk \
            -v bnd_gt="${bnd_gt}" \
            -v bnd_lt="${bnd_lt}" \
            'BEGIN { exit (bnd_gt < bnd_lt) ? 0 : 1 }'
    then
        echo_err_func "${FUNCNAME[0]}" \
            "'bnd_gt' must be less than 'bnd_lt', but got '${bnd_gt}' and" \
            "'${bnd_lt}'."
        return 1
    fi

    quiet="$(normalize_bool "${quiet}" "quiet")" || return 1

    #  Record whether the caller already had 'pipefail' enabled, as gzip-stream
    #+ failures should propagate when reading compressed input
    case ":${SHELLOPTS:-}:" in
        *:pipefail:*) pipefail_set=1 ;;
    esac

    set -o pipefail

    #  Define AWK program locally (as it is only used by this function)
    # shellcheck disable=SC2016
    awk_unity='
        /^[[:space:]]*$/        { next }
        /^[[:space:]]*#/        { next }
        /^[[:space:]]*track/    { next }
        /^[[:space:]]*browser/  { next }

        { sum += $4 }

        END {
            if (sum >= bnd_gt && sum <= bnd_lt) {
                if (quiet != "true") {
                    print "File sums to (approximately) unity:", sum
                }
                exit 0
            } else {
                print "File does not sum to unity:", sum \
                    > "/dev/stderr"
                exit 1
            }
        }
    '

    if [[ "${fil_in}" == *.gz ]]; then
        gunzip -cd -- "${fil_in}" \
            | awk \
                -v bnd_gt="${bnd_gt}" \
                -v bnd_lt="${bnd_lt}" \
                -v quiet="${quiet}" \
                "${awk_unity}"
        status=$?
    else
        awk \
            -v bnd_gt="${bnd_gt}" \
            -v bnd_lt="${bnd_lt}" \
            -v quiet="${quiet}" \
            "${awk_unity}" \
            "${fil_in}"
        status=$?
    fi

    if (( ! pipefail_set )); then
        set +o pipefail
    fi

    return "${status}"
}


#  Print an error message when function script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    err_source_only "${BASH_SOURCE[0]}"
fi
