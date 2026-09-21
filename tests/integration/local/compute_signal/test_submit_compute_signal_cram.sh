#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_submit_compute_signal_cram.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6);
# - Anthropic Claude Code (Opus 5).
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="submit compute-signal CRAM"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Define fixture and output paths for local serial CRAM-backed tests.
dir_fx="${ROOT_REPO}/tests/fixtures/compute_signal"
in_se="${dir_fx}/cram/se/tiny_se.cram"
in_pe="${dir_fx}/cram/pe/tiny_pe.cram"
ref_fa="${dir_fx}/reference/tiny.fa"

tmp="${TEST_DIR_TMP}/submit_compute_signal_cram"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/compute_signal"

fil_out_se_signal="${dir_out}/tiny_se_signal_unadj.bdg"
fil_out_se_coord="${dir_out}/tiny_se_coord.bed"
fil_out_pe_signal="${dir_out}/tiny_pe_signal_unadj.bdg"
fil_out_pe_coord="${dir_out}/tiny_pe_coord.bed"
fil_out_missing_ref="${dir_out}/tiny_se_missing_ref.bdg"

log_se_signal="${dir_log}/submit_compute_signal_cram_se_signal.log"
log_se_coord="${dir_log}/submit_compute_signal_cram_se_coord.log"
log_pe_signal="${dir_log}/submit_compute_signal_cram_pe_signal.log"
log_pe_coord="${dir_log}/submit_compute_signal_cram_pe_coord.log"
log_missing_ref="${dir_log}/submit_compute_signal_cram_missing_ref.log"


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty \
    "${in_se}" \
    "${in_pe}" \
    "${ref_fa}" || {
    finish
    exit $?
}


# Signal mode: two 10-bp SE alignments produce chromosome-I bedGraph bins.
run_case_compute_signal \
    submit \
    cram \
    "se_signal" \
    "signal" \
    "${in_se}" \
    "${fil_out_se_signal}" \
    "${log_se_signal}" \
    "${dir_out}" \
    "${dir_err}" \
    "${ref_fa}" \
    --method unadj \
    --siz_bin 10 \
    --csv_scl_fct NA \
    --dp 3

assert_file_nonempty \
    "${fil_out_se_signal}" \
    "CRAM SE signal bedGraph output"

assert_pattern_found \
    "${log_se_signal}" \
    "samp=tiny_se" \
    "submit CRAM SE signal strips .cram from sample name"

assert_pattern_absent \
    "${log_se_signal}" \
    "samp=tiny_se.cram" \
    "submit CRAM SE signal sample name omits .cram suffix"

if [[ -s "${fil_out_se_signal}" ]]; then
    assert_pattern_found \
        "${fil_out_se_signal}" \
        $'^I\t0\t10\t10$' \
        "CRAM SE signal output has chromosome-I bin I:0-10"

    assert_pattern_found \
        "${fil_out_se_signal}" \
        $'^I\t20\t30\t10$' \
        "CRAM SE signal output has chromosome-I bin I:20-30"
fi


# Coord mode: the same SE alignments emit BED-like processed fragments.
run_case_compute_signal \
    submit \
    cram \
    "se_coord" \
    "coord" \
    "${in_se}" \
    "${fil_out_se_coord}" \
    "${log_se_coord}" \
    "${dir_out}" \
    "${dir_err}" \
    "${ref_fa}" \
    --dp 3

assert_file_nonempty \
    "${fil_out_se_coord}" \
    "CRAM SE coord BED output"

if [[ -s "${fil_out_se_coord}" ]]; then
    assert_pattern_found \
        "${fil_out_se_coord}" \
        $'^I\t0\t10\t10$' \
        "CRAM SE coord output has chromosome-I fragment I:0-10"

    assert_pattern_found \
        "${fil_out_se_coord}" \
        $'^I\t20\t30\t10$' \
        "CRAM SE coord output has chromosome-I fragment I:20-30"
fi


# Signal mode: two PE fragments cover bins from I:10-60.
run_case_compute_signal \
    submit \
    cram \
    "pe_signal" \
    "signal" \
    "${in_pe}" \
    "${fil_out_pe_signal}" \
    "${log_pe_signal}" \
    "${dir_out}" \
    "${dir_err}" \
    "${ref_fa}" \
    --method unadj \
    --siz_bin 10 \
    --csv_scl_fct NA \
    --dp 3

assert_file_nonempty \
    "${fil_out_pe_signal}" \
    "CRAM PE signal bedGraph output"

if [[ -s "${fil_out_pe_signal}" ]]; then
    assert_pattern_found \
        "${fil_out_pe_signal}" \
        $'^I\t10\t20\t10$' \
        "CRAM PE signal output has chromosome-I bin I:10-20"

    assert_pattern_found \
        "${fil_out_pe_signal}" \
        $'^I\t40\t50\t10$' \
        "CRAM PE signal output has chromosome-I bin I:40-50"
fi


# Coord mode: PE output emits one BED-like row per leftmost proper pair.
run_case_compute_signal \
    submit \
    cram \
    "pe_coord" \
    "coord" \
    "${in_pe}" \
    "${fil_out_pe_coord}" \
    "${log_pe_coord}" \
    "${dir_out}" \
    "${dir_err}" \
    "${ref_fa}" \
    --dp 3

assert_file_nonempty \
    "${fil_out_pe_coord}" \
    "CRAM PE coord BED output"

if [[ -s "${fil_out_pe_coord}" ]]; then
    assert_pattern_found \
        "${fil_out_pe_coord}" \
        $'^I\t10\t40\t30$' \
        "CRAM PE coord output has chromosome-I fragment I:10-40"

    assert_pattern_found \
        "${fil_out_pe_coord}" \
        $'^I\t40\t60\t20$' \
        "CRAM PE coord output has chromosome-I fragment I:40-60"
fi


# CRAM input without a reference FASTA should fail clearly in the wrapper.
# shellcheck disable=SC2154
if \
    run_capture \
        "submit compute-signal CRAM missing_ref" \
        "${log_missing_ref}" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/submit_compute_signal.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/bin" \
            --threads 1 \
            --mode signal \
            --csv_fil_in "${in_se}" \
            --csv_fil_out "${fil_out_missing_ref}" \
            --csv_usr_frg NA \
            --dir_eo "${dir_err}" \
            --nam_job "test_compute_cram_missing_ref" \
            --method unadj \
            --siz_bin 10 \
            --csv_scl_fct NA \
            --dp 3
then
    record_fail \
        "submit_compute_signal.sh CRAM without --ref_fa unexpectedly passed"
else
    record_pass "submit_compute_signal.sh CRAM without --ref_fa fails"
fi

assert_pattern_found \
    "${log_missing_ref}" \
    "'--ref_fa' is required when '--csv_fil_in' contains CRAM" \
    "submit CRAM missing-ref error mentions --ref_fa"

if [[ ! -s "${fil_out_missing_ref}" ]]; then
    record_pass "submit CRAM missing-ref output is absent or empty"
else
    record_fail "submit CRAM missing-ref output was unexpectedly written"
fi


# Window engine over CRAM, multi-threaded, as previously CRAM suites had no
# '--siz_win' case, and nothing paired '--siz_win' with '--threads' above one,
# in which case window tasks are distributed across cores.
fil_out_pe_window="${dir_out}/tiny_pe_window_t2.bdg"
log_pe_window="${dir_log}/submit_compute_signal_cram_pe_window.log"

run_case_compute_signal \
    submit \
    cram \
    "pe_window" \
    "signal" \
    "${in_pe}" \
    "${fil_out_pe_window}" \
    "${log_pe_window}" \
    "${dir_out}" \
    "${dir_err}" \
    "${ref_fa}" \
    --method unadj \
    --siz_bin 10 \
    --engine window \
    --siz_win 20 \
    --threads 2 \
    --csv_scl_fct NA

assert_file_nonempty \
    "${fil_out_pe_window}" \
    "submit CRAM PE window-engine multi-threaded output"

log_window="${tmp}/logs/test_compute_cram_pe_window.tiny_pe_window_t2.stderr.txt"

assert_file_nonempty \
    "${log_window}" \
    "submit CRAM window stderr log"

if [[ -s "${log_window}" ]]; then
    assert_pattern_found \
        "${log_window}" \
        "--siz_win *20" \
        "submit CRAM window run forwards '--siz_win'"

    assert_pattern_found \
        "${log_window}" \
        "--engine *window" \
        "submit CRAM window run forwards '--engine window'"
fi


# Reports over CRAM: the counter shares the signal path's reference-backed
# iterator, so '--ref_fa' must not disturb the counts.
dir_rep="${tmp}/reports"
mkdir -p "${dir_rep}"

bash "${ROOT_REPO}/bin/submit_compute_signal.sh" \
    --mode signal \
    --csv_fil_in "${in_pe}" \
    --csv_fil_out "${dir_rep}/pe.bdg" \
    --csv_report_n_frg "${dir_rep}/pe.n_frg.txt" \
    --csv_report_n_ovlp "${dir_rep}/pe.n_ovlp.txt" \
    --ref_fa "${ref_fa}" \
    --dir_eo "${dir_err}" \
    --siz_bin 10 \
    --method unadj \
    > /dev/null 2>&1 || true

assert_file_nonempty "${dir_rep}/pe.n_frg.txt" \
    "submit CRAM report run writes the fragment count"

if \
    PYTHONDONTWRITEBYTECODE=1 \
    "${TEST_MANAGED_PYTHON}" \
        -m protocol_chipseq_signal_norm.cli.compute_signal \
        --fil_in "${dir_fx}/bam/pe/tiny_pe.bam" \
        --siz_bin 10 \
        --report_n_frg "${dir_rep}/bam.n_frg.txt" \
        --report_n_ovlp "${dir_rep}/bam.n_ovlp.txt" \
        > /dev/null 2>&1
then
    assert_files_equal \
        "${dir_rep}/pe.n_frg.txt" \
        "${dir_rep}/bam.n_frg.txt" \
        "CRAM fragment count equals the BAM count for the same alignments"

    assert_files_equal \
        "${dir_rep}/pe.n_ovlp.txt" \
        "${dir_rep}/bam.n_ovlp.txt" \
        "CRAM overlap count equals the BAM value for the same alignments"
else
    record_fail "BAM reference counts for the CRAM comparison failed"
fi


finish
