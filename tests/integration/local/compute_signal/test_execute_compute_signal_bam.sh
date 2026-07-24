#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_compute_signal_bam.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="execute compute-signal BAM"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


#  Define fixture and output paths for the execute-to-submit BAM path
dir_fx="${ROOT_REPO}/tests/fixtures/compute_signal"
chr_sizes="${dir_fx}/reference/tiny.fa.fai"
in_se="${dir_fx}/bam/se/tiny_se.bam"
in_pe="${dir_fx}/bam/pe/tiny_pe.bam"

tmp="${TEST_DIR_TMP}/execute_compute_signal_bam"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/compute_signal"

fil_out_se_signal="${dir_out}/tiny_se.bdg"
log_se_signal="${dir_log}/execute_compute_signal_bam_se_signal.log"

fil_out_se_coord="${dir_out}/tiny_se.bed"
log_se_coord="${dir_log}/execute_compute_signal_bam_se_coord.log"

fil_out_pe_signal="${dir_out}/tiny_pe.bdg"
log_pe_signal="${dir_log}/execute_compute_signal_bam_pe_signal.log"

fil_out_pe_coord="${dir_out}/tiny_pe.bed"
log_pe_coord="${dir_log}/execute_compute_signal_bam_pe_coord.log"

fil_out_se_signal_scaled="${dir_out}/scaled.tiny_se.bdg"
log_se_signal_scaled="${dir_log}/execute_compute_signal_bam_se_signal_scaled.log"

fil_out_se_signal_usr_frg="${dir_out}/usr_frg_signal.tiny_se.bdg"
log_se_signal_usr_frg="${dir_log}/execute_compute_signal_bam_se_signal_usr_frg.log"

fil_out_se_coord_usr_frg="${dir_out}/usr_frg_coord.tiny_se.bed"
log_se_coord_usr_frg="${dir_log}/execute_compute_signal_bam_se_coord_usr_frg.log"


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty \
    "${chr_sizes}" \
    "${in_se}" \
    "${in_pe}" || {
    finish
    exit $?
}


#  Signal mode: two 10-bp SE alignments produce chromosome-I bedGraph bins
run_case_compute_signal \
    execute \
    bam \
    "se_signal" \
    "signal" \
    "${in_se}" \
    "bdg" \
    "${log_se_signal}" \
    "${dir_out}" \
    "${dir_err}" \
    "" \
    --chr_sizes "${chr_sizes}" \
    --method unadj \
    --siz_bin 10 \
    --engine window \
    --csv_scl_fct NA

assert_file_nonempty \
    "${fil_out_se_signal}" \
    "execute SE signal bedGraph output"

if [[ -s "${fil_out_se_signal}" ]]; then
    assert_pattern_found \
        "${fil_out_se_signal}" \
        $'^I\t0\t10\t10$' \
        "execute SE signal output has chromosome-I bin I:0-10"

    assert_pattern_found \
        "${fil_out_se_signal}" \
        $'^I\t20\t30\t10$' \
        "execute SE signal output has chromosome-I bin I:20-30"
fi


#  Coord mode: the same SE alignments emit BED-like processed fragments
run_case_compute_signal \
    execute \
    bam \
    "se_coord" \
    "coord" \
    "${in_se}" \
    "bed" \
    "${log_se_coord}" \
    "${dir_out}" \
    "${dir_err}" \
    ""

assert_file_nonempty \
    "${fil_out_se_coord}" \
    "execute SE coord BED output"

if [[ -s "${fil_out_se_coord}" ]]; then
    assert_pattern_found \
        "${fil_out_se_coord}" \
        $'^I\t0\t10\t10$' \
        "execute SE coord output has chromosome-I fragment I:0-10"

    assert_pattern_found \
        "${fil_out_se_coord}" \
        $'^I\t20\t30\t10$' \
        "execute SE coord output has chromosome-I fragment I:20-30"
fi


#  Signal mode: two PE fragments cover bins from I:10-60
run_case_compute_signal \
    execute \
    bam \
    "pe_signal" \
    "signal" \
    "${in_pe}" \
    "bdg" \
    "${log_pe_signal}" \
    "${dir_out}" \
    "${dir_err}" \
    "" \
    --method unadj \
    --siz_bin 10 \
    --csv_scl_fct NA

assert_file_nonempty \
    "${fil_out_pe_signal}" \
    "execute PE signal bedGraph output"

if [[ -s "${fil_out_pe_signal}" ]]; then
    assert_pattern_found \
        "${fil_out_pe_signal}" \
        $'^I\t10\t20\t10$' \
        "execute PE signal output has chromosome-I bin I:10-20"

    assert_pattern_found \
        "${fil_out_pe_signal}" \
        $'^I\t40\t50\t10$' \
        "execute PE signal output has chromosome-I bin I:40-50"
fi


#  Coord mode: PE output emits one BED-like row per leftmost proper pair
run_case_compute_signal \
    execute \
    bam \
    "pe_coord" \
    "coord" \
    "${in_pe}" \
    "bed" \
    "${log_pe_coord}" \
    "${dir_out}" \
    "${dir_err}" \
    ""

assert_file_nonempty \
    "${fil_out_pe_coord}" \
    "execute PE coord BED output"

if [[ -s "${fil_out_pe_coord}" ]]; then
    assert_pattern_found \
        "${fil_out_pe_coord}" \
        $'^I\t10\t40\t30$' \
        "execute PE coord output has chromosome-I fragment I:10-40"

    assert_pattern_found \
        "${fil_out_pe_coord}" \
        $'^I\t40\t60\t20$' \
        "execute PE coord output has chromosome-I fragment I:40-60"
fi


#  Scaling factor and prefix propagation: raw 10-bp SE bins scaled by 2
run_case_compute_signal \
    execute \
    bam \
    "se_signal_scaled" \
    "signal" \
    "${in_se}" \
    "bdg" \
    "${log_se_signal_scaled}" \
    "${dir_out}" \
    "${dir_err}" \
    "" \
    --method unadj \
    --siz_bin 10 \
    --csv_scl_fct 2 \
    --prefix scaled

assert_file_nonempty \
    "${fil_out_se_signal_scaled}" \
    "execute scaled SE signal output"

if [[ -s "${fil_out_se_signal_scaled}" ]]; then
    assert_pattern_found \
        "${fil_out_se_signal_scaled}" \
        $'^I\t0\t10\t20$' \
        "execute scaled SE signal output has I:0-10 = 20"

    assert_pattern_found \
        "${fil_out_se_signal_scaled}" \
        $'^I\t20\t30\t20$' \
        "execute scaled SE signal output has I:20-30 = 20"
fi


#  Fixed SE fragment length extends reads to 20 bp in signal mode
run_case_compute_signal \
    execute \
    bam \
    "se_signal_usr_frg" \
    "signal" \
    "${in_se}" \
    "bdg" \
    "${log_se_signal_usr_frg}" \
    "${dir_out}" \
    "${dir_err}" \
    "" \
    --method unadj \
    --siz_bin 10 \
    --csv_scl_fct NA \
    --csv_usr_frg 20 \
    --prefix usr_frg_signal

assert_file_nonempty \
    "${fil_out_se_signal_usr_frg}" \
    "execute usr_frg SE signal output"

if [[ -s "${fil_out_se_signal_usr_frg}" ]]; then
    assert_pattern_found \
        "${fil_out_se_signal_usr_frg}" \
        $'^I\t0\t10\t10$' \
        "execute usr_frg SE signal output has I:0-10 = 10"

    assert_pattern_found \
        "${fil_out_se_signal_usr_frg}" \
        $'^I\t10\t20\t20$' \
        "execute usr_frg SE signal output has I:10-20 = 20"

    assert_pattern_found \
        "${fil_out_se_signal_usr_frg}" \
        $'^I\t20\t30\t10$' \
        "execute usr_frg SE signal output has I:20-30 = 10"
fi


#  Fixed SE fragment length is reflected in coord-mode BED intervals
run_case_compute_signal \
    execute \
    bam \
    "se_coord_usr_frg" \
    "coord" \
    "${in_se}" \
    "bed" \
    "${log_se_coord_usr_frg}" \
    "${dir_out}" \
    "${dir_err}" \
    "" \
    --csv_usr_frg 20 \
    --prefix usr_frg_coord

assert_file_nonempty \
    "${fil_out_se_coord_usr_frg}" \
    "execute usr_frg SE coord output"

if [[ -s "${fil_out_se_coord_usr_frg}" ]]; then
    assert_pattern_found \
        "${fil_out_se_coord_usr_frg}" \
        $'^I\t0\t20\t20$' \
        "execute usr_frg SE coord output has chromosome-I fragment I:0-20"

    assert_pattern_found \
        "${fil_out_se_coord_usr_frg}" \
        $'^I\t10\t30\t20$' \
        "execute usr_frg SE coord output has chromosome-I fragment I:10-30"
fi

finish
