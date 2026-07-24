#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_submit_compute_signal_bam.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="submit compute-signal BAM"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


#  Define fixture and output paths for local serial BAM-backed tests
dir_fx="${ROOT_REPO}/tests/fixtures/compute_signal"
chr_sizes="${dir_fx}/reference/tiny.fa.fai"
in_se="${dir_fx}/bam/se/tiny_se.bam"
in_pe="${dir_fx}/bam/pe/tiny_pe.bam"

tmp="${TEST_DIR_TMP}/submit_compute_signal_bam"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/compute_signal"

fil_out_se_signal="${dir_out}/tiny_se_signal_unadj.bdg"
fil_out_se_coord="${dir_out}/tiny_se_coord.bed"
fil_out_pe_signal="${dir_out}/tiny_pe_signal_unadj.bdg"
fil_out_pe_coord="${dir_out}/tiny_pe_coord.bed"
fil_out_se_signal_scaled="${dir_out}/tiny_se_signal_scaled.bdg"
fil_out_se_signal_frag="${dir_out}/tiny_se_signal_frag.bdg"
fil_out_se_signal_norm="${dir_out}/tiny_se_signal_norm.bdg"
fil_out_se_signal_usr_frg="${dir_out}/tiny_se_signal_usr_frg.bdg"
fil_out_se_coord_usr_frg="${dir_out}/tiny_se_coord_usr_frg.bed"

log_se_signal="${dir_log}/submit_compute_signal_bam_se_signal.log"
log_se_coord="${dir_log}/submit_compute_signal_bam_se_coord.log"
log_pe_signal="${dir_log}/submit_compute_signal_bam_pe_signal.log"
log_pe_coord="${dir_log}/submit_compute_signal_bam_pe_coord.log"
log_se_signal_scaled="${dir_log}/submit_compute_signal_bam_se_signal_scaled.log"
log_se_signal_frag="${dir_log}/submit_compute_signal_bam_se_signal_frag.log"
log_se_signal_norm="${dir_log}/submit_compute_signal_bam_se_signal_norm.log"
log_se_signal_usr_frg="${dir_log}/submit_compute_signal_bam_se_signal_usr_frg.log"
log_se_coord_usr_frg="${dir_log}/submit_compute_signal_bam_se_coord_usr_frg.log"


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
    submit \
    bam \
    "se_signal" \
    "signal" \
    "${in_se}" \
    "${fil_out_se_signal}" \
    "${log_se_signal}" \
    "${dir_out}" \
    "${dir_err}" \
    "" \
    --chr_sizes "${chr_sizes}" \
    --method unadj \
    --siz_bin 10 \
    --engine window \
    --csv_scl_fct NA \
    --dp 3

assert_file_nonempty \
    "${fil_out_se_signal}" \
    "signal bedGraph output"

if [[ -s "${fil_out_se_signal}" ]]; then
    assert_pattern_found \
        "${fil_out_se_signal}" \
        $'^I\t0\t10\t10$' \
        "SE signal output has chromosome-I bin I:0-10"

    assert_pattern_found \
        "${fil_out_se_signal}" \
        $'^I\t20\t30\t10$' \
        "SE signal output has chromosome-I bin I:20-30"
fi


#  Coord mode: the same SE alignments emit BED-like processed fragments
run_case_compute_signal \
    submit \
    bam \
    "se_coord" \
    "coord" \
    "${in_se}" \
    "${fil_out_se_coord}" \
    "${log_se_coord}" \
    "${dir_out}" \
    "${dir_err}" \
    "" \
    --dp 3

assert_file_nonempty \
    "${fil_out_se_coord}" \
    "coord BED output"

if [[ -s "${fil_out_se_coord}" ]]; then
    assert_pattern_found \
        "${fil_out_se_coord}" \
        $'^I\t0\t10\t10$' \
        "SE coord output has chromosome-I fragment I:0-10"

    assert_pattern_found \
        "${fil_out_se_coord}" \
        $'^I\t20\t30\t10$' \
        "SE coord output has chromosome-I fragment I:20-30"
fi


#  Signal mode: two PE fragments cover bins from I:10-60
run_case_compute_signal \
    submit \
    bam \
    "pe_signal" \
    "signal" \
    "${in_pe}" \
    "${fil_out_pe_signal}" \
    "${log_pe_signal}" \
    "${dir_out}" \
    "${dir_err}" \
    "" \
    --method unadj \
    --siz_bin 10 \
    --csv_scl_fct NA \
    --dp 3

assert_file_nonempty \
    "${fil_out_pe_signal}" \
    "PE signal bedGraph output"

if [[ -s "${fil_out_pe_signal}" ]]; then
    assert_pattern_found \
        "${fil_out_pe_signal}" \
        $'^I\t10\t20\t10$' \
        "PE signal output has chromosome-I bin I:10-20"

    assert_pattern_found \
        "${fil_out_pe_signal}" \
        $'^I\t40\t50\t10$' \
        "PE signal output has chromosome-I bin I:40-50"
fi


#  Coord mode: PE output emits one BED-like row per leftmost proper pair
run_case_compute_signal \
    submit \
    bam \
    "pe_coord" \
    "coord" \
    "${in_pe}" \
    "${fil_out_pe_coord}" \
    "${log_pe_coord}" \
    "${dir_out}" \
    "${dir_err}" \
    "" \
    --dp 3

assert_file_nonempty \
    "${fil_out_pe_coord}" \
    "PE coord BED output"

if [[ -s "${fil_out_pe_coord}" ]]; then
    assert_pattern_found \
        "${fil_out_pe_coord}" \
        $'^I\t10\t40\t30$' \
        "PE coord output has chromosome-I fragment I:10-40"

    assert_pattern_found \
        "${fil_out_pe_coord}" \
        $'^I\t40\t60\t20$' \
        "PE coord output has chromosome-I fragment I:40-60"
fi


#  Scaling factor and prefix propagation: raw 10-bp SE bins scaled by 2
run_case_compute_signal \
    submit \
    bam \
    "se_signal_scaled" \
    "signal" \
    "${in_se}" \
    "${fil_out_se_signal_scaled}" \
    "${log_se_signal_scaled}" \
    "${dir_out}" \
    "${dir_err}" \
    "" \
    --method unadj \
    --siz_bin 10 \
    --csv_scl_fct 2 \
    --dp 3

assert_file_nonempty \
    "${fil_out_se_signal_scaled}" \
    "scaled SE signal output"

if [[ -s "${fil_out_se_signal_scaled}" ]]; then
    assert_pattern_found \
        "${fil_out_se_signal_scaled}" \
        $'^I\t0\t10\t20$' \
        "scaled SE signal output has I:0-10 = 20"

    assert_pattern_found \
        "${fil_out_se_signal_scaled}" \
        $'^I\t20\t30\t20$' \
        "scaled SE signal output has I:20-30 = 20"
fi


#  Fragment-length normalization: each 10-bp SE fragment contributes 1
run_case_compute_signal \
    submit \
    bam \
    "se_signal_frag" \
    "signal" \
    "${in_se}" \
    "${fil_out_se_signal_frag}" \
    "${log_se_signal_frag}" \
    "${dir_out}" \
    "${dir_err}" \
    "" \
    --method frag \
    --siz_bin 10 \
    --csv_scl_fct NA \
    --dp 3

assert_file_nonempty \
    "${fil_out_se_signal_frag}" \
    "frag-normalized SE signal output"

if [[ -s "${fil_out_se_signal_frag}" ]]; then
    assert_pattern_found \
        "${fil_out_se_signal_frag}" \
        $'^I\t0\t10\t1$' \
        "frag-normalized SE signal output has I:0-10 = 1"

    assert_pattern_found \
        "${fil_out_se_signal_frag}" \
        $'^I\t20\t30\t1$' \
        "frag-normalized SE signal output has I:20-30 = 1"
fi


#  Normalized coverage divides fragment-normalized signal by total fragments
run_case_compute_signal \
    submit \
    bam \
    "se_signal_norm" \
    "signal" \
    "${in_se}" \
    "${fil_out_se_signal_norm}" \
    "${log_se_signal_norm}" \
    "${dir_out}" \
    "${dir_err}" \
    "" \
    --method norm \
    --siz_bin 10 \
    --csv_scl_fct NA \
    --dp 3

assert_file_nonempty \
    "${fil_out_se_signal_norm}" \
    "norm SE signal output"

if [[ -s "${fil_out_se_signal_norm}" ]]; then
    assert_pattern_found \
        "${fil_out_se_signal_norm}" \
        $'^I\t0\t10\t0.5$' \
        "norm SE signal output has I:0-10 = 0.5"

    assert_pattern_found \
        "${fil_out_se_signal_norm}" \
        $'^I\t20\t30\t0.5$' \
        "norm SE signal output has I:20-30 = 0.5"
fi


#  Fixed SE fragment length extends reads to 20 bp in signal mode
run_case_compute_signal \
    submit \
    bam \
    "se_signal_usr_frg" \
    "signal" \
    "${in_se}" \
    "${fil_out_se_signal_usr_frg}" \
    "${log_se_signal_usr_frg}" \
    "${dir_out}" \
    "${dir_err}" \
    "" \
    --method unadj \
    --siz_bin 10 \
    --csv_scl_fct NA \
    --csv_usr_frg 20 \
    --dp 3

assert_file_nonempty \
    "${fil_out_se_signal_usr_frg}" \
    "usr_frg SE signal output"

if [[ -s "${fil_out_se_signal_usr_frg}" ]]; then
    assert_pattern_found \
        "${fil_out_se_signal_usr_frg}" \
        $'^I\t0\t10\t10$' \
        "usr_frg SE signal output has I:0-10 = 10"

    assert_pattern_found \
        "${fil_out_se_signal_usr_frg}" \
        $'^I\t10\t20\t20$' \
        "usr_frg SE signal output has I:10-20 = 20"

    assert_pattern_found \
        "${fil_out_se_signal_usr_frg}" \
        $'^I\t20\t30\t10$' \
        "usr_frg SE signal output has I:20-30 = 10"
fi


#  Fixed SE fragment length is reflected in coord-mode BED intervals
run_case_compute_signal \
    submit \
    bam \
    "se_coord_usr_frg" \
    "coord" \
    "${in_se}" \
    "${fil_out_se_coord_usr_frg}" \
    "${log_se_coord_usr_frg}" \
    "${dir_out}" \
    "${dir_err}" \
    "" \
    --csv_usr_frg 20 \
    --dp 3

assert_file_nonempty \
    "${fil_out_se_coord_usr_frg}" \
    "usr_frg SE coord output"

if [[ -s "${fil_out_se_coord_usr_frg}" ]]; then
    assert_pattern_found \
        "${fil_out_se_coord_usr_frg}" \
        $'^I\t0\t20\t20$' \
        "usr_frg SE coord output has chromosome-I fragment I:0-20"

    assert_pattern_found \
        "${fil_out_se_coord_usr_frg}" \
        $'^I\t10\t30\t20$' \
        "usr_frg SE coord output has chromosome-I fragment I:10-30"
fi

finish
