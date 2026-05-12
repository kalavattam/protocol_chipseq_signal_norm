#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_compute_signal_parallel.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="execute compute-signal GNU Parallel"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

rec_section "${TEST_NAME}"

if ! is_integration; then
    rec_skip \
        "GNU Parallel dry-run/config integration check disabled;" \
        "set RUN_INTEGRATION=1 or RUN_CONDA_INTEGRATION=1 to enable"
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

rm -rf "${tmp}"
mkdir -p "${dir_in}" "${dir_out}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_exist "${fil_A}" "${fil_B}" || {
    finish
    exit $?
}

cp "${fil_A}" "${fil_A_1}"
cp "${fil_A}" "${fil_A_2}"

require_files_exist "${fil_A_1}" "${fil_A_2}" || {
    finish
    exit $?
}

#  GNU Parallel dry-run config should invoke non-executable submit scripts
#+ through Bash
# shellcheck disable=SC2154
if ! \
    require_parallel_env \
        "${env_nam}" \
        "${dir_log}/execute_compute_signal_parallel_env.log"
then
    finish
    exit $?
fi

config="${dir_err}/test_execute_compute_parallel.config_parallel.txt"
log="${dir_log}/execute_compute_signal_parallel_dry_run.log"

if \
    run_capture \
        "execute compute-signal GNU Parallel dry-run" "${log}" \
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
    rec_pass "execute_compute_signal.sh GNU Parallel dry-run exits 0"
else
    rec_fail \
        "execute_compute_signal.sh GNU Parallel dry-run failed; see" \
        "$(rec_relpath "${log}")"
fi

assert_file_nonempty "${config}" "execute GNU Parallel config"

if [[ -s "${config}" ]]; then
    assert_grep_pattern \
        "${config}" \
        "${TEST_BASH} ${ROOT_REPO}/scripts/submit_compute_signal.sh" \
        "execute GNU Parallel config uses Bash-prefixed submit command"
fi


#  Real GNU Parallel execution should produce one output per ratio job
outfile_1="${dir_out}/exec_parallel_wet_ratio_A_1.bdg"
outfile_2="${dir_out}/exec_parallel_wet_ratio_A_2.bdg"
log="${dir_log}/execute_compute_signal_parallel_wet_run.log"

if \
    run_capture \
        "execute compute-signal GNU Parallel wet run" "${log}" \
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
    rec_pass "execute_compute_signal.sh GNU Parallel wet run exits 0"
else
    rec_fail \
        "execute_compute_signal.sh GNU Parallel wet run failed; see" \
        "$(rec_relpath "${log}")"
fi

assert_file_nonempty "${outfile_1}" "execute GNU Parallel wet output 1"
assert_file_nonempty "${outfile_2}" "execute GNU Parallel wet output 2"

for outfile in "${outfile_1}" "${outfile_2}"; do
    if [[ -s "${outfile}" ]]; then
        assert_grep_pattern "${outfile}" $'^I\t0\t10\t2$' \
            "$(basename "${outfile}") has I:0-10 = 2"
        assert_grep_pattern "${outfile}" $'^I\t40\t50\t4$' \
            "$(basename "${outfile}") has I:40-50 = 4"
        assert_grep_pattern "${outfile}" $'^I\t60\t70\t0.333$' \
            "$(basename "${outfile}") has I:60-70 = 0.333"
    fi
done

finish
