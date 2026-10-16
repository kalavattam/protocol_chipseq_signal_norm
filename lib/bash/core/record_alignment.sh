#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: record_alignment.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


# print_record_commit
# add_record_run_lines


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
    _dir_src_rec="$(
        cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
    )"

    # shellcheck source=lib/bash/core/source_helpers.sh
    source "${_dir_src_rec}/source_helpers.sh" || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source '${_dir_src_rec}/source_helpers.sh'." >&2

        if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
            return 1
        else
            exit 1
        fi
    }

    source_helpers "${_dir_src_rec}" \
        check_source format_outputs || {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to source required helper dependencies." >&2

            if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
                return 1
            else
                exit 1
            fi
        }

    unset _dir_src_rec
}


# Print the commit of the checkout holding these libraries, for the pipeline
# record of an alignment file ('OUTPUT.ALIGNMENT.RECORD').
function print_record_commit() {
    local dir_rep commit
    local show_help

    show_help=$(cat << EOM
Usage
-----
  print_record_commit
    [--help]

  Print the commit of the Git checkout that holds these libraries, for the 'pipeline: commit=' line of an alignment file's pipeline record.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

Returns
-------
  Prints the full commit hash, with '-dirty' appended when tracked files have uncommitted changes, or 'unknown' when Git is missing or the libraries are not in a Git checkout; returns 0.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - dirname
    - git (optional; without it, the commit is 'unknown')

  - Untracked files do not mark the commit dirty.

Examples
--------
  1. Print the commit for a record line.
    '''bash
    printf '@CO\\tpipeline: commit=%s\\n' "\$(print_record_commit)"
    '''

  2. Display the function's argument contract.
    '''bash
    print_record_commit --help
    '''
EOM
    )

    if [[ "${1:-}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    dir_rep="$(dirname "${BASH_SOURCE[0]}")"

    if ! \
        commit="$(git -C "${dir_rep}" rev-parse HEAD 2> /dev/null)"
    then
        printf '%s\n' "unknown"
        return 0
    fi

    if [[ -n "$(
        git -C "${dir_rep}" status --porcelain --untracked-files=no \
            2> /dev/null
    )" ]]; then
        commit+="-dirty"
    fi

    printf '%s\n' "${commit}"
}


# Append the run lines of an alignment file's pipeline record
# ('OUTPUT.ALIGNMENT.RECORD') to a named array.
function add_record_run_lines() {
    local rec_nam="${1:-}"
    local fnc="${2:-}"
    local threads="${3:-}"
    local calls idx job key pth
    local -a arr_pth
    local show_help

    show_help=$(cat << EOM
Usage
-----
  add_record_run_lines
    [--help] rec_nam fnc threads [key path ...]

  Append the 'pipeline-run:' lines of an alignment file's pipeline record to the named array as '@CO' header lines: the thread count, the start time, the host, the call path, the environment, any Slurm job and task, and then the full path of each named file.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  rec_nam : str
    Name of the indexed array to append to, one header line per element.

  2  fnc : str
    The library function writing the record, as 'file.sh::function'; it ends the call path.

  3  threads : int
    Number of threads the chain uses.

  4+  key path : str
    Pairs of an option name without its leading dashes and the path given for it; a pair with an empty path is skipped, and a relative path is made absolute against the working directory.

Returns
-------
  Returns 0 after appending; prints an error and returns 1 if 'rec_nam' is missing or the pairs are incomplete.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - date
    - hostname

  - The start time is the time of the call, so call it when the chain starts.
  - The call path is 'PIPELINE_CALLS', to which each 'execute_*' and 'submit_*' script adds its name, then 'fnc', joined with ' > '; a direct library call gives 'fnc' alone.
  - In a Slurm array task, the job is the array's job ID ('SLURM_ARRAY_JOB_ID'), as in the log names and 'sacct', with the task ID beside it; otherwise it is 'SLURM_JOB_ID'.

Examples
--------
  1. Append the run lines for a filter run and print them.
    '''bash
    arr_rec_run=()
    add_record_run_lines arr_rec_run filter_alignment.sh::filter_alignment_sc 2 \\
        fil_in work/sample.bam \\
        fil_out work/sample.sc.bam
    printf '%s\\n' "\${arr_rec_run[@]}"
    '''

  2. Display the function's argument contract.
    '''bash
    add_record_run_lines --help
    '''
EOM
    )

    if [[ "${rec_nam}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${rec_nam}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'rec_nam', is missing."
        return 1
    elif (( ($# - 3) % 2 != 0 )); then
        echo_err_func "${FUNCNAME[0]}" \
            "each 'key' needs a 'path'."
        return 1
    fi

    local -n run_ref="${rec_nam}"
    shift 3

    calls="${PIPELINE_CALLS:+${PIPELINE_CALLS} > }${fnc}"

    run_ref+=(
        $'@CO\tpipeline-run: --threads='"${threads}"
        $'@CO\tpipeline-run: start='"$(date '+%Y-%m-%dT%H:%M:%S%z')"
        $'@CO\tpipeline-run: host='"$(hostname)"
        $'@CO\tpipeline-run: calls='"${calls}"
        $'@CO\tpipeline-run: env='"${CONDA_DEFAULT_ENV:-none}"
    )

    # In an array, give the array's job ID, which names the logs and 'sacct'
    # rows; each task's own 'SLURM_JOB_ID' differs from it.
    job="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-}}"
    if [[ -n "${job}" ]]; then
        run_ref+=( $'@CO\tpipeline-run: slurm_job='"${job}" )
    fi

    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        run_ref+=(
            $'@CO\tpipeline-run: slurm_task='"${SLURM_ARRAY_TASK_ID}"
        )
    fi

    arr_pth=( "$@" )
    for (( idx = 0; idx < ${#arr_pth[@]}; idx += 2 )); do
        key="${arr_pth[idx]}"
        pth="${arr_pth[idx+1]}"

        if [[ -z "${pth}" ]]; then continue; fi
        if [[ "${pth}" != /* ]]; then pth="${PWD}/${pth}"; fi

        run_ref+=( $'@CO\tpipeline-run: --'"${key}=${pth}" )
    done
}


# Print an error message when function script is executed directly.
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    err_source_only "${BASH_SOURCE[0]}"
fi
