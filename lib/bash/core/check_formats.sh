#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: check_formats.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


# canonicalize_fmt
# fmt_of_path
# check_fmt_path
# strip_fmt_suffix


# Require Bash >= 4.4 before defining functions.
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

# Source required helper functions if needed.
{
    _dir_src_fmt="$(
        cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
    )"

    # shellcheck source=lib/bash/core/source_helpers.sh
    source "${_dir_src_fmt}/source_helpers.sh" || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source '${_dir_src_fmt}/source_helpers.sh'." >&2

        if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
            return 1
        else
            exit 1
        fi
    }

    source_helpers "${_dir_src_fmt}" \
        check_source format_outputs || {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to source required helper dependencies." >&2

            if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
                return 1
            else
                exit 1
            fi
        }

    unset _dir_src_fmt
}


# Resolve a format name, in any letter case, to its canonical spelling.
function canonicalize_fmt() {
    local val="${1:-}"
    local key gz=""
    local show_help

    show_help=$(cat << EOM
Usage
-----
  canonicalize_fmt
    [--help] val

  Resolve a format name, in any letter case, to its canonical spelling and print it: 'bam', 'cram', 'bed', or 'bedGraph', with '.gz' kept for 'bed' and 'bedGraph'.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  val : str
    Format name as supplied, such as 'bedgraph.gz' or 'CRAM'.

Returns
-------
  0 and the canonical name on stdout if 'val' names a format; 1 with an error if 'val' is a short bedGraph spelling ('bdg', 'bg', optionally with '.gz') or is missing; 2, silently, if 'val' names no format, so the caller decides what to do.

Notes
-----
  Runtime requirements:
    bash >= 4.4

  Only 'bedGraph' names the bedGraph format; the short spellings are refused, not coerced.

Examples
--------
  1. Normalize a format name supplied in another letter case.
    '''bash
    typ_out="\$(canonicalize_fmt "BEDGRAPH.gz")"  # bedGraph.gz
    '''

  2. Refuse a short bedGraph spelling.
    '''bash
    canonicalize_fmt "bdg.gz" || echo "refused" >&2
    '''
EOM
    )

    if [[ "${val}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${val}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'val', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    key="${val,,}"
    if [[ "${key}" == *.gz ]]; then
        gz=".gz"
        key="${key%.gz}"
    fi

    case "${key}" in
        bdg|bg)
            echo_err_func "${FUNCNAME[0]}" \
                "'${val}' is not accepted; use 'bedGraph${gz}'."
            return 1
            ;;

        bedgraph) printf 'bedGraph%s\n' "${gz}" ;;
        bed)      printf 'bed%s\n' "${gz}" ;;

        bam|cram)
            if [[ -n "${gz}" ]]; then return 2; fi
            printf '%s\n' "${key}"
            ;;

        *) return 2 ;;
    esac
}


# Print a file's canonical format, recognized from its suffix in any case.
function fmt_of_path() {
    local pth="${1:-}"
    local base gz="" ext
    local show_help

    show_help=$(cat << EOM
Usage
-----
  fmt_of_path
    [--help] pth

  Recognize a file's format from its suffix, in any letter case, and print the canonical format: 'bam', 'cram', 'bed', or 'bedGraph', with '.gz' appended when the file name ends in it ('bed' and 'bedGraph' only).

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  pth : str
    File path, such as 'x.CRAM' or 'x.BEDGRAPH.gz'.

Returns
-------
  0 and the canonical format on stdout if the suffix names one; 1 with an error if the suffix is a short bedGraph spelling ('.bdg', '.bg', optionally with '.gz') or 'pth' is missing; 2, silently, if the suffix names no format.

Notes
-----
  Runtime requirements:
    bash >= 4.4

  The file is never renamed; only its format is recognized.

Examples
--------
  1. Recognize a compressed bedGraph track in any letter case.
    '''bash
    fmt="\$(fmt_of_path "track.BEDGRAPH.gz")"  # bedGraph.gz
    '''

  2. Refuse a short bedGraph suffix.
    '''bash
    fmt_of_path "track.bdg.gz" || echo "rename it" >&2
    '''
EOM
    )

    if [[ "${pth}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${pth}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'pth', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    check_fmt_path "${pth}" || return 1

    base="${pth##*/}"
    if [[ "${base,,}" == *.gz ]]; then
        gz=".gz"
        base="${base:0:${#base}-3}"
    fi

    if [[ "${base}" != *.* ]]; then return 2; fi
    ext="${base##*.}"

    case "${ext,,}" in
        bedgraph) printf 'bedGraph%s\n' "${gz}" ;;
        bed)      printf 'bed%s\n' "${gz}" ;;

        bam|cram)
            if [[ -n "${gz}" ]]; then return 2; fi
            printf '%s\n' "${ext,,}"
            ;;

        *) return 2 ;;
    esac
}


# Refuse a file path whose name uses a short bedGraph suffix.
function check_fmt_path() {
    local pth="${1:-}"
    local base gz="" ext
    local show_help

    show_help=$(cat << EOM
Usage
-----
  check_fmt_path
    [--help] pth

  Refuse a file path whose name ends in a short bedGraph suffix ('.bdg' or '.bg', optionally with '.gz', in any letter case), naming the accepted one.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  pth : str
    File path to check.

Returns
-------
  0 if the path does not use a short bedGraph suffix; 1 with an error if it does or if 'pth' is missing.

Notes
-----
  Runtime requirements:
    bash >= 4.4

  Any other name passes, so piped inputs and temporary files keep working.

Examples
--------
  1. Accept a bedGraph file.
    '''bash
    check_fmt_path "sample.bedGraph.gz"
    '''

  2. Refuse a short bedGraph suffix.
    '''bash
    check_fmt_path "sample.bdg.gz" || echo "rename it" >&2
    '''
EOM
    )

    if [[ "${pth}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${pth}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'pth', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    base="${pth##*/}"
    if [[ "${base,,}" == *.gz ]]; then
        gz="${base: -3}"
        base="${base:0:${#base}-3}"
    fi

    ext="${base##*.}"
    if [[ "${base}" == *.* && "${ext,,}" =~ ^(bdg|bg)$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'${pth}': '.${ext}${gz}' is not accepted; the file name must" \
            "end in '.bedGraph${gz:+.gz}'."
        return 1
    fi
}


# Print a file name without its format suffix, matched in any letter case.
function strip_fmt_suffix() {
    local nam="${1:-}"
    local ext
    local show_help

    show_help=$(cat << EOM
Usage
-----
  strip_fmt_suffix
    [--help] nam

  Print 'nam' without a trailing '.gz' and then one format suffix ('.bedGraph', '.bed', '.bam', '.cram', or '.sam'), each matched in any letter case.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  nam : str
    File name or path.

Returns
-------
  0 and the stripped name on stdout; 1 if 'nam' is missing.

Notes
-----
  Runtime requirements:
    bash >= 4.4

  '.sam' is stripped so that log names derived from SAM input match those of the other formats; SAM is not otherwise a recognized format.

Examples
--------
  1. Strip a bedGraph suffix supplied in another letter case.
    '''bash
    strip_fmt_suffix "IP_sample.BEDGRAPH.gz"  # IP_sample
    '''

  2. Derive a sample name from an output path.
    '''bash
    samp="\$(strip_fmt_suffix "\$(basename "\${fil_out}")")"
    '''
EOM
    )

    if [[ "${nam}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${nam}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'nam', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    if [[ "${nam,,}" == *.gz ]]; then nam="${nam:0:${#nam}-3}"; fi

    for ext in bedgraph bed bam cram sam; do
        if [[ "${nam,,}" == *."${ext}" ]]; then
            nam="${nam:0:${#nam}-${#ext}-1}"
            break
        fi
    done

    printf '%s\n' "${nam}"
}


# Print an error message when function script is executed directly.
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    err_source_only "${BASH_SOURCE[0]}"
fi
