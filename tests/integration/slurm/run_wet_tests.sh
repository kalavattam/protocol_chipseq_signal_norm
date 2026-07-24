#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: run_wet_tests.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


#  Require Bash >= 4.4 before parsing or running the wet-test interface
if [[ -z "${BASH_VERSION:-}" ]]; then
    echo "error(shell): this script requires Bash >= 4.4." >&2
    exit 1
elif ((
    BASH_VERSINFO[0] < 4 || ( BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] < 4 )
)); then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "this script requires Bash >= 4.4; current version is" \
        "'${BASH_VERSION}'." >&2
    exit 1
fi

set -euo pipefail


function show_help() {
    cat << 'EOM'
Usage
-----
  run_slurm_wet_tests.sh
    [--help] --config <file>

  Run the dedicated, bounded Slurm wet-integration suite and emit its structured result bundle.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  -cf, --config : file
    Runtime configuration produced by 'tests/integration/slurm/coordinator.py remote-launch'.

Returns
-------
  Returns 0 only when remote preflight succeeds and every required job reaches a successful terminal state with passing output assertions; otherwise, returns nonzero.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - python3 >= 3.11
    - sacct
    - sbatch
    - scontrol
    - squeue

  - Every scheduler action requires all three exact gates: 'RUN_SLURM=1', 'WAIT_SLURM=1', and 'CONFIRM_SLURM_WET=1'.
  - This runner is intentionally separate from 'tests/run_tests.sh'.
  - It uses two tiny committed-fixture jobs and writes only beneath the configured run directory.

Examples
--------
  1. Run a prepared wet-validation session after logging in to the Slurm host.
    '''bash
    RUN_SLURM=1 WAIT_SLURM=1 CONFIRM_SLURM_WET=1 \
        bash tests/integration/slurm/run_wet_tests.sh \
            --config ../runtime_config.json
    '''

  2. Display the wet-runner interface without checking gates or the scheduler.
    '''bash
    bash tests/integration/slurm/run_wet_tests.sh --help
    '''
EOM
}


function require_wet_gate() {
    local name="${1:-}"
    local value="${!name:-}"

    if [[ "${value}" != "1" ]]; then
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "${name}=1 is required before any wet Slurm action." >&2
        return 1
    fi
}


function parse_args() {
    config=""

    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
            -h|--help)
                show_help
                return 2
                ;;

            -cf|--config)
                if [[ -z "${2:-}" || "${2}" == -* ]]; then
                    echo "error($(basename "${BASH_SOURCE[0]}")):" \
                        "'${1}' requires a following file argument." >&2
                    return 1
                fi
                config="${2}"
                shift 2
                ;;

            *)
                echo "error($(basename "${BASH_SOURCE[0]}")):" \
                    "unknown option/parameter passed: '${1}'." >&2
                return 1
                ;;
        esac
    done

    if [[ -z "${config}" || ! -f "${config}" ]]; then
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "'--config' must name a readable runtime configuration file." >&2
        return 1
    fi
}


function main() {
    local config=""
    local status=0
    local dir_repo=""

    parse_args "$@" || {
        status=$?
        [[ "${status}" -eq 2 ]] && return 0
        return "${status}"
    }

    require_wet_gate RUN_SLURM
    require_wet_gate WAIT_SLURM
    require_wet_gate CONFIRM_SLURM_WET

    if [[ -z "${SLURM_REMOTE_PYTHON:-}" || \
        "${SLURM_REMOTE_PYTHON}" != /* || \
        ! -x "${SLURM_REMOTE_PYTHON}" ]]
    then
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "SLURM_REMOTE_PYTHON must name an executable absolute path." >&2
        return 1
    fi

    if ! "${SLURM_REMOTE_PYTHON}" -c \
        'import sys; raise SystemExit(sys.version_info < (3, 11))'
    then
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "Python >= 3.11 required." >&2
        return 1
    fi

    dir_repo="$(
        cd "$(dirname "${BASH_SOURCE[0]}")/../../.." > /dev/null 2>&1 && pwd
    )"
    PYTHONDONTWRITEBYTECODE=1 "${SLURM_REMOTE_PYTHON}" \
        "${dir_repo}/tests/integration/slurm/coordinator.py" \
        wet-run --config "${config}"
}


main "$@"
