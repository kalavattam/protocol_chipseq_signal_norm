#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_merge_bins_bdg.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in design, development,
# and documentation, with all output reviewed, edited, and approved by the
# author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="merge_bins_bdg.py smoke"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


function write_exact_fixture() {
    local out="${1:-}"

    cat > "${out}" << 'EOF'
chrI 0 10 1
chrI 10 20 1
track type=bedGraph
chrI 20 30 1.000
chrI 30 40 1.000
chrI 40 50 NaN
chrI 50 60 nan
EOF
}


function write_exact_expected() {
    local out="${1:-}"

    cat > "${out}" << 'EOF'
chrI	0	20	1
track type=bedGraph
chrI	20	40	1.000
chrI	40	50	nan
chrI	50	60	nan
EOF
}


function write_round_fixture() {
    local out="${1:-}"

    cat > "${out}" << 'EOF'
chrI 0 10 1.001
chrI 10 20 1.002
chrI 20 30 1.014
EOF
}


function write_round_expected() {
    local out="${1:-}"

    cat > "${out}" << 'EOF'
chrI	0	20	1.00
chrI	20	30	1.01
EOF
}


dir_tmp="${TEST_DIR_TMP}/merge_bins_bdg"
dir_log="${TEST_DIR_LOG}/merge_bins_bdg"

scr="${ROOT_REPO}/src/protocol_chipseq_signal_norm/cli/merge_bins_bdg.py"
in_exact="${dir_tmp}/exact.input.bdg"
exp_exact="${dir_tmp}/exact.expected.bdg"
obs_exact="${dir_tmp}/exact.observed.bdg"
in_round="${dir_tmp}/round.input.bdg"
exp_round="${dir_tmp}/round.expected.bdg"
obs_round="${dir_tmp}/round.observed.bdg"
bad_input="${dir_tmp}/malformed.input.bdg"
log_env_rslv="${dir_log}/resolve_env_python.log"
log_help="${dir_log}/help.log"
log_exact="${dir_log}/exact_stdin_stdout.log"
log_round="${dir_log}/round_file_output.log"
log_bad="${dir_log}/malformed_input.log"


print_section "${TEST_NAME}"

mkdir -p "${dir_tmp}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

if [[ -n "${CONDA_DEFAULT_ENV:-}" && "${CONDA_DEFAULT_ENV}" != "base" ]]; then
    if ! \
        py="$(find_python)"
    then
        record_fail "active project environment has no python/python3 on PATH"
        finish
        exit $?
    fi

    py_cmd=( "${py}" )
else
    if \
        run_capture \
            "resolve env python" \
            "${log_env_rslv}" \
            conda run -n "${env_nam}" python -c \
                'import sys; print(sys.executable)'
    then
        IFS= read -r py < "${log_env_rslv}"
    else
        record_fail \
            "failed to resolve python from '${env_nam}'; see" \
            "$(print_relpath "${log_env_rslv}")"
        finish
        exit $?
    fi

    if [[ -z "${py}" || ! -x "${py}" ]]; then
        record_fail \
            "resolved python is not executable; see" \
            "$(print_relpath "${log_env_rslv}")"
        finish
        exit $?
    fi

    py_cmd=( "${py}" )
fi

write_exact_fixture "${in_exact}"
write_exact_expected "${exp_exact}"
write_round_fixture "${in_round}"
write_round_expected "${exp_round}"
printf 'chrI 0 10\n' > "${bad_input}"

if \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONPYCACHEPREFIX="${TEST_DIR_OUT}/pycache" \
    PYTHONPATH="${ROOT_REPO}/src" \
    run_capture \
        "merge_bins_bdg help" \
        "${log_help}" \
        "${py_cmd[@]}" "${scr}" --help
then
    record_pass "merge_bins_bdg.py --help exits successfully"
    assert_pattern_found \
        "${log_help}" \
        'Usage:' \
        "merge_bins_bdg.py help includes Usage"
else
    record_fail "merge_bins_bdg.py --help failed; see $(print_relpath "${log_help}")"
fi

if \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONPYCACHEPREFIX="${TEST_DIR_OUT}/pycache" \
    PYTHONPATH="${ROOT_REPO}/src" \
    run_capture \
        "merge_bins_bdg exact stdin stdout" \
        "${log_exact}" \
        "${py_cmd[@]}" "${scr}" -fi - -fo - < "${in_exact}"
then
    cp "${log_exact}" "${obs_exact}"
    assert_files_equal \
        "${obs_exact}" \
        "${exp_exact}" \
        "merge_bins_bdg.py exact stdin/stdout output"
else
    record_fail \
        "merge_bins_bdg.py exact stdin/stdout failed; see" \
        "$(print_relpath "${log_exact}")"
fi

if \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONPYCACHEPREFIX="${TEST_DIR_OUT}/pycache" \
    PYTHONPATH="${ROOT_REPO}/src" \
    run_capture \
        "merge_bins_bdg rounded file output" \
        "${log_round}" \
        "${py_cmd[@]}" "${scr}" \
            -fi "${in_round}" \
            -fo "${obs_round}" \
            --dp 2
then
    assert_files_equal \
        "${obs_round}" \
        "${exp_round}" \
        "merge_bins_bdg.py rounded file output"
else
    record_fail \
        "merge_bins_bdg.py rounded file output failed; see" \
        "$(print_relpath "${log_round}")"
fi

if \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONPYCACHEPREFIX="${TEST_DIR_OUT}/pycache" \
    PYTHONPATH="${ROOT_REPO}/src" \
    run_capture \
        "merge_bins_bdg malformed input" \
        "${log_bad}" \
        "${py_cmd[@]}" "${scr}" -fi "${bad_input}" -fo -
then
    record_fail "merge_bins_bdg.py malformed input unexpectedly succeeded"
else
    record_pass "merge_bins_bdg.py malformed input exits nonzero"
    assert_pattern_found \
        "${log_bad}" \
        'Malformed bedGraph line' \
        "merge_bins_bdg.py malformed input reports parser error"
fi

finish
