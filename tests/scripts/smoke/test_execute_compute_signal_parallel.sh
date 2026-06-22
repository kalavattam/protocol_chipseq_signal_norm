#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_compute_signal_parallel.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="execute compute-signal GNU Parallel"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

print_section "${TEST_NAME}"

if ! \
    is_parallel_enabled
then
    record_skip \
        "GNU Parallel compute-signal check disabled; set RUN_PARALLEL=1 to" \
        "enable"
    finish
    exit $?
fi

#  Define fixture and output paths for lightweight execute-layer config checks
dir_fx="${ROOT_REPO}/tests/compute_signal/fixtures/bedgraph"
fil_A="${dir_fx}/ratio_A.bdg"
fil_B="${dir_fx}/ratio_B.bdg"

tmp="${TEST_DIR_TMP}/execute_compute_signal_parallel"
dir_in="${tmp}/in"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/compute_signal"
fil_A_1="${dir_in}/ratio_A_1.bdg"
fil_A_2="${dir_in}/ratio_A_2.bdg"

cfg_parallel="${dir_err}/test_execute_compute_parallel.config_parallel.txt"
fil_out_wet_1="${dir_out}/exec_parallel_wet_ratio_A_1.bdg"
fil_out_wet_2="${dir_out}/exec_parallel_wet_ratio_A_2.bdg"

log_env_parallel="${dir_log}/execute_compute_signal_parallel_env.log"
log_dry_run="${dir_log}/execute_compute_signal_parallel_dry_run.log"
log_wet_run="${dir_log}/execute_compute_signal_parallel_wet_run.log"

rm -rf "${tmp}"
mkdir -p "${dir_in}" "${dir_out}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty \
    "${fil_A}" \
    "${fil_B}" || {
    finish
    exit $?
}

cp "${fil_A}" "${fil_A_1}"
cp "${fil_A}" "${fil_A_2}"

require_files_nonempty \
    "${fil_A_1}" \
    "${fil_A_2}" || {
    finish
    exit $?
}

#  GNU Parallel dry-run config should invoke non-executable submit scripts
#+ through Bash
# shellcheck disable=SC2154
if ! \
    require_env_parallel \
        "${env_nam}" \
        "${log_env_parallel}"
then
    finish
    exit $?
fi

if \
    run_capture \
        "execute compute-signal GNU Parallel dry-run" \
        "${log_dry_run}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/execute_compute_signal.sh" \
            --dry_run \
            --threads 2 \
            --mode ratio \
            --method unadj \
            --csv_fil_A "${fil_A}" \
            --csv_fil_B "${fil_B}" \
            --dir_out "${dir_out}" \
            --typ_out bdg \
            --prefix exec_parallel \
            --eps 0 \
            --dp 3 \
            --err_out "${dir_err}" \
            --nam_job "test_execute_compute_parallel" \
            --max_job 2
then
    record_pass "execute_compute_signal.sh GNU Parallel dry-run exits 0"
else
    record_fail \
        "execute_compute_signal.sh GNU Parallel dry-run failed; see" \
        "$(print_relpath "${log_dry_run}")"
fi

assert_file_nonempty \
    "${cfg_parallel}" \
    "execute GNU Parallel config"

if [[ -s "${cfg_parallel}" ]]; then
    assert_pattern_found \
        "${cfg_parallel}" \
        "${TEST_BASH} ${ROOT_REPO}/scripts/submit_compute_signal.sh" \
        "execute GNU Parallel config uses Bash-prefixed submit command"
fi


#  Real GNU Parallel execution should produce one output per ratio job
if \
    run_capture \
        "execute compute-signal GNU Parallel wet run" \
        "${log_wet_run}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/execute_compute_signal.sh" \
            --threads 2 \
            --mode ratio \
            --method unadj \
            --csv_fil_A "${fil_A_1},${fil_A_2}" \
            --csv_fil_B "${fil_B},${fil_B}" \
            --dir_out "${dir_out}" \
            --typ_out bdg \
            --prefix exec_parallel_wet \
            --eps 0 \
            --dp 3 \
            --err_out "${dir_err}" \
            --nam_job "test_execute_compute_parallel_wet" \
            --max_job 2
then
    record_pass "execute_compute_signal.sh GNU Parallel wet run exits 0"
else
    record_fail \
        "execute_compute_signal.sh GNU Parallel wet run failed; see" \
        "$(print_relpath "${log_wet_run}")"
fi

assert_file_nonempty \
    "${fil_out_wet_1}" \
    "execute GNU Parallel wet output 1"

assert_file_nonempty \
    "${fil_out_wet_2}" \
    "execute GNU Parallel wet output 2"

for out in "${fil_out_wet_1}" "${fil_out_wet_2}"; do
    if [[ -s "${out}" ]]; then
        assert_pattern_found \
            "${out}" \
            $'^I\t0\t10\t2$' \
            "$(basename "${out}") has I:0-10 = 2"

        assert_pattern_found \
            "${out}" \
            $'^I\t40\t50\t4$' \
            "$(basename "${out}") has I:40-50 = 4"

        assert_pattern_found \
            "${out}" \
            $'^I\t60\t70\t0.333$' \
            "$(basename "${out}") has I:60-70 = 0.333"
    fi
done

finish
