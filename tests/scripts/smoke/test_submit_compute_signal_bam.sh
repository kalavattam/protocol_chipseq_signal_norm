#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_submit_compute_signal_bam.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="submit compute-signal BAM"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

rec_section "${TEST_NAME}"

#  Define fixture and output paths for local serial BAM-backed smoke tests
dir_fx="${ROOT_REPO}/tests/compute_signal/fixtures/bam"
infile_se="${dir_fx}/tiny_se.bam"
infile_pe="${dir_fx}/tiny_pe.bam"

tmp="${TEST_DIR_TMP}/compute_signal_bam"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/compute_signal"

rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_err}" "${dir_log}"

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

if [[ ! -s "${infile_se}" ]]; then
    rec_fail "missing fixture $(rec_relpath "${infile_se}")"
fi

if [[ ! -s "${infile_pe}" ]]; then
    rec_fail "missing fixture $(rec_relpath "${infile_pe}")"
fi

if (( TEST_FAIL > 0 )); then
    finish
    exit $?
fi


#  Assert that a generated output file exists and is non-empty
function assert_file_nonempty() {
    local file="${1:-}"
    local label="${2:-output}"

    if [[ -s "${file}" ]]; then
        rec_pass "${label} exists and is non-empty"
    else
        rec_fail "${label} missing or empty: $(rec_relpath "${file}")"
    fi
}


#  Run a local serial submit-wrapper BAM case into compute_signal.py
function run_case_bam() {
    local nam_case="${1:-}"
    local mode="${2:-}"
    local infile_lcl="${3:-}"
    local outfile_lcl="${4:-}"
    local log_lcl="${5:-}"

    shift 5

    if \
        run_capture \
            "submit compute-signal BAM ${nam_case}" "${log_lcl}" \
            "${TEST_BASH}" "${ROOT_REPO}/scripts/submit_compute_signal.sh" \
                --env_nam "${env_nam}" \
                --dir_scr "${ROOT_REPO}/scripts" \
                --threads 1 \
                --mode "${mode}" \
                --csv_infile "${infile_lcl}" \
                --csv_outfile "${outfile_lcl}" \
                --csv_usr_frg NA \
                --err_out "${dir_err}" \
                --nam_job "test_compute_bam_${nam_case}" \
                "$@"
    then
        rec_pass "submit_compute_signal.sh ${mode} ${nam_case} exits 0"
    else
        rec_fail \
            "submit_compute_signal.sh ${mode} ${nam_case} failed; see" \
            "$(rec_relpath "${log_lcl}")"
    fi
}


#  Signal mode: two 10-bp SE alignments produce chromosome-I bedGraph bins
outfile="${dir_out}/tiny_se_signal_unadj.bdg"
log="${dir_log}/submit_compute_signal_bam_se_signal.log"

run_case_bam \
    "se_signal" "signal" "${infile_se}" "${outfile}" "${log}" \
    --method unadj \
    --siz_bin 10 \
    --csv_scl_fct NA \
    --dp 3

assert_file_nonempty "${outfile}" "signal bedGraph output"

if [[ -s "${outfile}" ]]; then
    assert_grep_pattern "${outfile}" $'^I\t0\t10\t10$' \
        "SE signal output has chromosome-I bin I:0-10"
    assert_grep_pattern "${outfile}" $'^I\t20\t30\t10$' \
        "SE signal output has chromosome-I bin I:20-30"
fi


#  Coord mode: the same SE alignments emit BED-like processed fragments
outfile="${dir_out}/tiny_se_coord.bed"
log="${dir_log}/submit_compute_signal_bam_se_coord.log"

run_case_bam \
    "se_coord" "coord" "${infile_se}" "${outfile}" "${log}" \
    --dp 3

assert_file_nonempty "${outfile}" "coord BED output"

if [[ -s "${outfile}" ]]; then
    assert_grep_pattern "${outfile}" $'^I\t0\t10\t10$' \
        "SE coord output has chromosome-I fragment I:0-10"
    assert_grep_pattern "${outfile}" $'^I\t20\t30\t10$' \
        "SE coord output has chromosome-I fragment I:20-30"
fi


#  Signal mode: two PE fragments cover bins from I:10-60
outfile="${dir_out}/tiny_pe_signal_unadj.bdg"
log="${dir_log}/submit_compute_signal_bam_pe_signal.log"

run_case_bam \
    "pe_signal" "signal" "${infile_pe}" "${outfile}" "${log}" \
    --method unadj \
    --siz_bin 10 \
    --csv_scl_fct NA \
    --dp 3

assert_file_nonempty "${outfile}" "PE signal bedGraph output"

if [[ -s "${outfile}" ]]; then
    assert_grep_pattern "${outfile}" $'^I\t10\t20\t10$' \
        "PE signal output has chromosome-I bin I:10-20"
    assert_grep_pattern "${outfile}" $'^I\t40\t50\t10$' \
        "PE signal output has chromosome-I bin I:40-50"
fi


#  Coord mode: PE output emits one BED-like row per leftmost proper pair
outfile="${dir_out}/tiny_pe_coord.bed"
log="${dir_log}/submit_compute_signal_bam_pe_coord.log"

run_case_bam \
    "pe_coord" "coord" "${infile_pe}" "${outfile}" "${log}" \
    --dp 3

assert_file_nonempty "${outfile}" "PE coord BED output"

if [[ -s "${outfile}" ]]; then
    assert_grep_pattern "${outfile}" $'^I\t10\t40\t30$' \
        "PE coord output has chromosome-I fragment I:10-40"
    assert_grep_pattern "${outfile}" $'^I\t40\t60\t20$' \
        "PE coord output has chromosome-I fragment I:40-60"
fi

finish
