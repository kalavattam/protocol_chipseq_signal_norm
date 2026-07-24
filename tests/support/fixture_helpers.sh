#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: fixture_helpers.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


#  Print a fixture-generator error and stop
function die() {
    local src="${BASH_SOURCE[1]:-${0}}"

    echo "error($(basename "${src}")):" "$@" >&2
    exit 1
}


#  Print a fixture-generator note
function note() {
    local src="${BASH_SOURCE[1]:-${0}}"

    echo "note($(basename "${src}")):" "$@" >&2
}


#  Print a fixture-generator success message
function succeed() {
    local src="${BASH_SOURCE[1]:-${0}}"

    echo "success($(basename "${src}")):" "$@" >&2
}


#  Require one command to be available
function require_cmd() {
    local cmd="${1:-}"
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  require_cmd
    [--help] cmd [context...]

  Require one executable command to be available in PATH.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  cmd : str
    Command name to resolve with 'command -v'.

  2+ context : str
    Optional diagnostic context appended to failure output.

Returns
-------
  Returns 0 when the command is available. Exits through 'die' when the command name is missing or unavailable.

  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Require gzip before generating compressed FASTQ fixtures.
    '''bash
    # shellcheck source=tests/support/fixture_helpers.sh
    source tests/support/fixture_helpers.sh
    require_cmd gzip
    '''

  2. Capture the expected fatal diagnostic for a missing fixture tool.
    '''bash
    # shellcheck source=tests/support/fixture_helpers.sh
    source tests/support/fixture_helpers.sh
    if ! (
        require_cmd definitely_missing_fixture_tool 'while building FASTQs'
    ); then
        printf '%s\n' 'missing tool rejected as expected'
    fi
    '''
EOM
    )

    if [[ "${cmd}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    shift || true

    if [[ -z "${cmd}" ]]; then
        die "'require_cmd' requires a command name."
    fi

    if ! \
        command -v "${cmd}" > /dev/null 2>&1
    then
        die "'${cmd}' must be available" "$@"
    fi
}


#  Require one or more commands to be available
function require_cmds() {
    local context="${1:-}"
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  require_cmds
    [--help] context cmd [cmd...]

  Require one or more executable commands to be available in PATH.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  context : str
    Diagnostic context passed to 'require_cmd'.

  2+ cmd : str
    Command names to resolve with 'command -v'.

Returns
-------
  Returns 0 when all commands are available. Exits through 'die' if no command names are supplied or any command is unavailable.

  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Require the Samtools commands used by an alignment-fixture recipe.
    '''bash
    # shellcheck source=tests/support/fixture_helpers.sh
    source tests/support/fixture_helpers.sh
    require_cmds 'while building BAM and CRAM fixtures' samtools gzip
    '''

  2. Capture the expected fatal diagnostic when no command names are supplied.
    '''bash
    # shellcheck source=tests/support/fixture_helpers.sh
    source tests/support/fixture_helpers.sh
    if ! (
        require_cmds 'empty dependency list'
    ); then
        printf '%s\n' 'empty command list rejected as expected'
    fi
    '''
EOM
    )

    if [[ "${context}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    shift || true

    local cmd=""

    if (( $# == 0 )); then
        die "'require_cmds' requires at least one command name."
    fi

    for cmd in "$@"; do
        require_cmd "${cmd}" "${context}"
    done
}


#  Require a named Conda environment to be active
function require_env() {
    local env_nam="${1:-}"
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  require_env
    [--help] env_nam [context...]

  Require a named Conda environment to be active before fixture
  generation continues.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  env_nam : str
    Conda environment to activate. Required environment name.

  2+ context : str
    Optional diagnostic context appended to failure output.

Returns
-------
  Returns 0 when the requested environment is active. Exits through 'die' when the environment name is missing or the active environment differs.

  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Continue when the requested fixture environment is active.
    '''bash
    # shellcheck source=tests/support/fixture_helpers.sh
    source tests/support/fixture_helpers.sh
    CONDA_DEFAULT_ENV=env_protocol
    require_env env_protocol
    '''

  2. Capture the expected fatal diagnostic for an environment mismatch.
    '''bash
    # shellcheck source=tests/support/fixture_helpers.sh
    source tests/support/fixture_helpers.sh
    CONDA_DEFAULT_ENV=base
    if ! (
        require_env env_protocol 'for alignment fixtures'
    ); then
        printf '%s\n' 'inactive fixture environment rejected as expected'
    fi
    '''
EOM
    )

    if [[ "${env_nam}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    shift || true

    if [[ -z "${env_nam}" ]]; then
        die "'require_env' requires an environment name."
    fi

    if [[ "${CONDA_DEFAULT_ENV:-}" != "${env_nam}" ]]; then
        die \
            "activate '${env_nam}' before generating fixtures;" \
            "current environment: '${CONDA_DEFAULT_ENV:-none}'" \
            "$@"
    fi
}


#  Create one or more directories
function mkdirs() {
    if (( $# == 0 )); then
        die "'mkdirs' requires at least one directory path."
    fi

    mkdir -p "$@"
}


#  Remove a generated fixture file inside a fixture root
function rm_file() {
    local dir_fix="${1:-}"
    local file="${2:-}"
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  rm_file
    [--help] dir_fix file

  Remove one generated fixture file after validating that it is inside the
  fixture root.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  dir_fix : dir
    Fixture root that bounds allowed removals.

  2  file : file
    Generated fixture file to remove.

Returns
-------
  Returns 0 after removal succeeds. Exits through 'die' if the path is empty, outside the fixture root, or a directory.

  Runtime requirements:
    - bash >= 4.4
    - rm

Examples
--------
  1. Remove one generated file contained by its fixture root.
    '''bash
    # shellcheck source=tests/support/fixture_helpers.sh
    source tests/support/fixture_helpers.sh
    dir_fix="$(mktemp -d)"
    : > "${dir_fix}/generated.bam"
    rm_file "${dir_fix}" "${dir_fix}/generated.bam"
    rmdir "${dir_fix}"
    '''

  2. Refuse an attempted removal outside the fixture root.
    '''bash
    # shellcheck source=tests/support/fixture_helpers.sh
    source tests/support/fixture_helpers.sh
    dir_fix="$(mktemp -d)"
    if ! (
        rm_file "${dir_fix}" /private/tmp/outside-fixture.bam
    ); then
        printf '%s\n' 'outside path rejected as expected'
    fi
    rmdir "${dir_fix}"
    '''
EOM
    )

    if [[ "${dir_fix}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    if [[ -z "${dir_fix}" ]]; then
        die "refusing to remove a file without a fixture root."
    elif [[ -z "${file}" ]]; then
        die "refusing to remove an empty file path."
    elif [[ "${file}" != "${dir_fix}/"* ]]; then
        die "refusing to remove path outside fixture directory: '${file}'."
    elif [[ -d "${file}" ]]; then
        die "refusing to remove directory: '${file}'."
    fi

    rm -f -- "${file}"
}


#  Remove generated fixture files inside a fixture root
function rm_files() {
    local dir_fix="${1:-}"
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  rm_files
    [--help] dir_fix file [file...]

  Remove generated fixture files after validating that each path is inside the
  fixture root.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  dir_fix : dir
    Fixture root that bounds allowed removals.

  2+ file : file
    Generated fixture files to remove.

Returns
-------
  Returns 0 after all removals succeed. Exits through 'die' if no files are supplied or any path is invalid.

  Runtime requirements:
    - bash >= 4.4
    - rm

Examples
--------
  1. Remove two generated fixture files in one call.
    '''bash
    # shellcheck source=tests/support/fixture_helpers.sh
    source tests/support/fixture_helpers.sh
    dir_fix="$(mktemp -d)"
    : > "${dir_fix}/a.bam"
    : > "${dir_fix}/b.bam"
    rm_files "${dir_fix}" "${dir_fix}/a.bam" "${dir_fix}/b.bam"
    rmdir "${dir_fix}"
    '''

  2. Capture the expected fatal diagnostic when no files are supplied.
    '''bash
    # shellcheck source=tests/support/fixture_helpers.sh
    source tests/support/fixture_helpers.sh
    dir_fix="$(mktemp -d)"
    if ! (
        rm_files "${dir_fix}"
    ); then
        printf '%s\n' 'empty removal list rejected as expected'
    fi
    rmdir "${dir_fix}"
    '''
EOM
    )

    if [[ "${dir_fix}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    shift || true

    local file=""

    if (( $# == 0 )); then
        die "rm_files requires at least one file path."
    fi

    for file in "$@"; do
        rm_file "${dir_fix}" "${file}"
    done
}


#  Register a cleanup function for shell exit
function register_cleanup() {
    local fnc="${1:-}"
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  register_cleanup
    [--help] fnc

  Register a named cleanup function for the fixture script EXIT trap.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  fnc : str
    Name of a defined shell function to run on exit.

Returns
-------
  Returns 0 after registering the trap. Exits through 'die' when the function name is missing or undefined.

  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Register a defined cleanup function in a bounded subshell.
    '''bash
    # shellcheck source=tests/support/fixture_helpers.sh
    source tests/support/fixture_helpers.sh
    (cleanup_example() { printf '%s\n' cleaned; }; register_cleanup cleanup_example)
    '''

  2. Capture the expected fatal diagnostic for an undefined cleanup function.
    '''bash
    # shellcheck source=tests/support/fixture_helpers.sh
    source tests/support/fixture_helpers.sh
    if ! (
        register_cleanup missing_cleanup
    ); then
        printf '%s\n' 'undefined cleanup rejected as expected'
    fi
    '''
EOM
    )

    if [[ "${fnc}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    if [[ -z "${fnc}" ]]; then
        die "register_cleanup requires a function name."
    elif ! \
        declare -F "${fnc}" > /dev/null
    then
        die "cleanup function is not defined: '${fnc}'."
    fi

    # shellcheck disable=SC2064
    trap "${fnc}" EXIT
}


#  Write one tab-delimited line/row
function write_line() {
    local field
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  write_line
    [--help] field [field...]

  Print one tab-delimited row to stdout.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1+ field : str
    Field values to print in order.

Returns
-------
  Writes one tab-delimited line to stdout. Exits through 'die' if no fields are supplied.

  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Write a two-field TSV row.
    '''bash
    # shellcheck source=tests/support/fixture_helpers.sh
    source tests/support/fixture_helpers.sh
    write_line sample_A 1.25
    '''

  2. Write a four-field SAM header row.
    '''bash
    # shellcheck source=tests/support/fixture_helpers.sh
    source tests/support/fixture_helpers.sh
    write_line @SQ SN:I LN:80 M5:ad01ff07883bde3965f91f71e009c5b0
    '''
EOM
    )

    if [[ "${1:-}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    if (( $# == 0 )); then
        die \
            "'write_line', 'write_tsv_row', and 'write_sam_line' require at" \
            "least one field."
    fi

    printf '%s' "${1}"
    shift

    for field in "${@}"; do
        printf '\t%s' "${field}"
    done

    printf '\n'
}


#  Write one tab-delimited TSV row (alias for 'write_line')
function write_tsv_row() { write_line "${@}"; }


#  Write one tab-delimited SAM line (alias for 'write_line')
function write_sam_line() { write_line "${@}"; }


#  Write gzip output with deterministic headers
function gzip_n() {
    local fil_in="${1:-}"
    local fil_out="${2:-}"
    local show_help

    show_help=$(cat << 'EOM'
Usage
-----
  gzip_n
    [--help] fil_in fil_out

  Write gzip-compressed output with deterministic gzip headers.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  fil_in : file
    Input file path. Input file to compress.

  2  fil_out : file
    Output file path. Output gzip file.

Returns
-------
  Returns the exit status of 'gzip -n -c'. Exits through 'die' if either path is missing.

  Runtime requirements:
    - bash >= 4.4
    - gzip

Examples
--------
  1. Compress a FASTQ with deterministic gzip headers and read it back.
    '''bash
    # shellcheck source=tests/support/fixture_helpers.sh
    source tests/support/fixture_helpers.sh
    tmp="$(mktemp -d)"
    printf '%s\n' '@read' ACGT + IIII > "${tmp}/read.fastq"
    gzip_n "${tmp}/read.fastq" "${tmp}/read.fastq.gz"
    gzip -cd "${tmp}/read.fastq.gz"
    rm -r -- "${tmp}"
    '''

  2. Propagate gzip failure for a missing input file.
    '''bash
    # shellcheck source=tests/support/fixture_helpers.sh
    source tests/support/fixture_helpers.sh
    tmp="$(mktemp -d)"
    if ! \
        gzip_n "${tmp}/missing.fastq" "${tmp}/read.fastq.gz"
    then
        printf '%s\n' 'missing FASTQ rejected as expected'
    fi
    rm -r -- "${tmp}"
    '''
EOM
    )

    if [[ "${fil_in}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    if [[ -z "${fil_in}" ]]; then
        die "'gzip_n' requires an input path."
    elif [[ -z "${fil_out}" ]]; then
        die "'gzip_n' requires an output path."
    fi

    gzip -n -c "${fil_in}" > "${fil_out}"
}
