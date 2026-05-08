#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_submit_compute_signal_ratio.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="submit compute-signal ratio"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

rec_section "${TEST_NAME}"

#  Define fixture and output paths for the simplest wet ratio-mode path
dir_fx="${ROOT_REPO}/tests/compute_signal/fixtures/bedgraph"
fil_A="${dir_fx}/ratio_A.bdg"
fil_B="${dir_fx}/ratio_B.bdg"

tmp="${TEST_DIR_TMP}/compute_signal_ratio"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
outfile="${dir_out}/ratio_unadj.dp3.bdg"
log="${TEST_DIR_LOG}/compute_signal/submit_compute_signal_ratio.log"

rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_err}" "$(dirname "${log}")"

#  Use the active project environment when available; otherwise require the
#+ conventional env_protocol environment so the wrapper can activate it
if [[ -n "${CONDA_DEFAULT_ENV:-}" && "${CONDA_DEFAULT_ENV}" != "base" ]]; then
    env_nam="${CONDA_DEFAULT_ENV}"
else
    env_nam="env_protocol"

    if ! \
        check_cmd_exists conda
    then
        rec_fail \
            "setup requires active project env or conda to activate" \
            "'env_protocol'"
        finish
        exit $?
    fi

    if ! \
        conda env list 2> /dev/null \
            | awk -v env="${env_nam}" '
                $1 == env { found = 1 } END { exit !found }
            '
    then
        rec_fail "setup requires Conda/Mamba environment '${env_nam}'"
        finish
        exit $?
    fi
fi

if [[ ! -s "${fil_A}" ]]; then
    rec_fail "missing fixture $(rec_relpath "${fil_A}")"
fi

if [[ ! -s "${fil_B}" ]]; then
    rec_fail "missing fixture $(rec_relpath "${fil_B}")"
fi

if (( TEST_FAIL > 0 )); then
    finish
    exit $?
fi

#  Run the local serial submit wrapper path into compute_signal_ratio.py
if \
    run_capture \
        "submit compute-signal ratio" "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/submit_compute_signal.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/scripts" \
            --threads 1 \
            --mode ratio \
            --method unadj \
            --csv_fil_A "${fil_A}" \
            --csv_fil_B "${fil_B}" \
            --csv_outfile "${outfile}" \
            --csv_scl_fct NA \
            --csv_dep_min NA \
            --csv_pseudo NA \
            --eps 0 \
            --skip_00 NA \
            --dp 3 \
            --err_out "${dir_err}" \
            --nam_job "test_compute_ratio_unadj"
then
    rec_pass "submit_compute_signal.sh ratio unadj exits 0"
else
    rec_fail \
        "submit_compute_signal.sh ratio unadj failed; see" \
        "$(rec_relpath "${log}")"
fi

if [[ -s "${outfile}" ]]; then
    rec_pass "ratio output exists and is non-empty"
else
    rec_fail "ratio output missing or empty: $(rec_relpath "${outfile}")"
fi

if [[ -s "${outfile}" ]]; then
    assert_grep_pattern "${outfile}" $'^I\t0\t10\t2$' \
        "ratio output has I:0-10 = 2"
    assert_grep_pattern "${outfile}" $'^I\t10\t20\t0$' \
        "ratio output has I:10-20 = 0"
    assert_grep_pattern "${outfile}" $'^I\t40\t50\t4$' \
        "ratio output has I:40-50 = 4"
    assert_grep_pattern "${outfile}" $'^I\t60\t70\t0.333$' \
        "ratio output has I:60-70 = 0.333"
    assert_grep_pattern "${outfile}" $'^I\t70\t80\t1$' \
        "ratio output has I:70-80 = 1"
fi

finish
