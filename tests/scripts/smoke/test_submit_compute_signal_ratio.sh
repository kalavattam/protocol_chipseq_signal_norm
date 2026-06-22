#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_submit_compute_signal_ratio.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="submit compute-signal ratio"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

print_section "${TEST_NAME}"

dir_fx="${ROOT_REPO}/tests/compute_signal/fixtures/bedgraph"
fil_A="${dir_fx}/ratio_A.bdg"
fil_B="${dir_fx}/ratio_B.bdg"
fil_A_hdr="${dir_fx}/ratio_headers_A.bdg"
fil_B_hdr="${dir_fx}/ratio_headers_B.bdg"

tmp="${TEST_DIR_TMP}/submit_compute_signal_ratio"
dir_in="${tmp}/in"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/compute_signal"
fil_A_gz="${dir_in}/ratio_A.bdg.gz"
fil_B_gz="${dir_in}/ratio_B.bdg.gz"

fil_out_unadj="${dir_out}/ratio_unadj.dp3.bdg"
fil_out_scl_fct="${dir_out}/ratio_scl_fct_2_1.dp3.bdg"
fil_out_log2="${dir_out}/ratio_log2.dp3.bdg"
fil_out_unadj_r="${dir_out}/ratio_unadj_r.dp3.bdg"
fil_out_log2_r="${dir_out}/ratio_log2_r.dp3.bdg"
fil_out_dep_min="${dir_out}/ratio_dep_min_0p1.dp3.bdg"
fil_out_eps="${dir_out}/ratio_eps_0p05.dp3.bdg"
fil_out_pseudo="${dir_out}/ratio_pseudo_1_1.dp3.bdg"
fil_out_drp_nan="${dir_out}/ratio_drp_nan.dp3.bdg"
fil_out_skip_00="${dir_out}/ratio_skip_00_pre_scale.dp3.bdg"
fil_out_skip_00_post="${dir_out}/ratio_skip_00_post_scale.dp3.bdg"
fil_out_track="${dir_out}/ratio_track.dp3.bdg"
fil_out_rnd_alias="${dir_out}/ratio_rnd_alias.dp2.bdg"
fil_out_gzip_io="${dir_out}/ratio_gzip_io.dp3.bdg.gz"
fil_out_skp_pfx="${dir_out}/ratio_skp_pfx.dp3.bdg"

trackfile_track="${dir_out}/ratio_track.dp3.track.bdg"
outfile_txt_gzip_io="${dir_out}/ratio_gzip_io.dp3.bdg"

log_unadj="${dir_log}/submit_compute_signal_ratio_unadj.log"
log_scl_fct="${dir_log}/submit_compute_signal_ratio_scl_fct.log"
log_log2="${dir_log}/submit_compute_signal_ratio_log2.log"
log_unadj_r="${dir_log}/submit_compute_signal_ratio_unadj_r.log"
log_log2_r="${dir_log}/submit_compute_signal_ratio_log2_r.log"
log_dep_min="${dir_log}/submit_compute_signal_ratio_dep_min.log"
log_eps="${dir_log}/submit_compute_signal_ratio_eps.log"
log_pseudo="${dir_log}/submit_compute_signal_ratio_pseudo.log"
log_drp_nan="${dir_log}/submit_compute_signal_ratio_drp_nan.log"
log_skip_00="${dir_log}/submit_compute_signal_ratio_skip_00.log"
log_skip_00_post="${dir_log}/submit_compute_signal_ratio_skip_00_post_scale.log"
log_track="${dir_log}/submit_compute_signal_ratio_track.log"
log_rnd_alias="${dir_log}/submit_compute_signal_ratio_rnd_alias.log"
log_gzip_io="${dir_log}/submit_compute_signal_ratio_gzip_io.log"
log_skp_pfx="${dir_log}/submit_compute_signal_ratio_skp_pfx.log"


rm -rf "${tmp}"
mkdir -p "${dir_in}" "${dir_out}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty \
    "${fil_A}" \
    "${fil_B}" \
    "${fil_A_hdr}" \
    "${fil_B_hdr}" || {
    finish
    exit $?
}

if ! \
    check_cmd_exists gzip
then
    record_fail "gzip is required for ratio gzip I/O coverage"
    finish
    exit $?
fi

gzip -c "${fil_A}" > "${fil_A_gz}"
gzip -c "${fil_B}" > "${fil_B_gz}"

require_files_nonempty \
    "${fil_A_gz}" \
    "${fil_B_gz}" || {
    finish
    exit $?
}


#  Baseline unadjusted ratio with three-decimal rounding
run_case_compute_signal_ratio \
    submit \
    "unadj" \
    "${fil_out_unadj}" \
    "${log_unadj}" \
    "unadj" \
    "${fil_A}" \
    "${fil_B}" \
    "${dir_out}" \
    "${dir_err}" \
    --csv_scl_fct \
    NA \
    --csv_dep_min \
    NA \
    --csv_pseudo \
    NA \
    --eps \
    0 \
    --skip_00 \
    NA \
    --dp \
    3

assert_file_nonempty \
    "${fil_out_unadj}" \
    "baseline ratio output"

if [[ -s "${fil_out_unadj}" ]]; then
    assert_pattern_found \
        "${fil_out_unadj}" \
        $'^I\t0\t10\t2$' \
        "baseline ratio output has I:0-10 = 2"

    assert_pattern_found \
        "${fil_out_unadj}" \
        $'^I\t10\t20\t0$' \
        "baseline ratio output has I:10-20 = 0"

    assert_pattern_found \
        "${fil_out_unadj}" \
        $'^I\t40\t50\t4$' \
        "baseline ratio output has I:40-50 = 4"

    assert_pattern_found \
        "${fil_out_unadj}" \
        $'^I\t60\t70\t0.333$' \
        "baseline ratio output has I:60-70 = 0.333"

    assert_pattern_found \
        "${fil_out_unadj}" \
        $'^I\t70\t80\t1$' \
        "baseline ratio output has I:70-80 = 1"
fi


#  Scaling factors are applied before ratio calculation: (2 * A) / (1 * B)
run_case_compute_signal_ratio \
    submit \
    "scl_fct" \
    "${fil_out_scl_fct}" \
    "${log_scl_fct}" \
    "unadj" \
    "${fil_A}" \
    "${fil_B}" \
    "${dir_out}" \
    "${dir_err}" \
    --csv_scl_fct \
    2:1 \
    --csv_dep_min \
    NA \
    --csv_pseudo \
    NA \
    --eps \
    0 \
    --skip_00 \
    NA \
    --dp \
    3

assert_file_nonempty \
    "${fil_out_scl_fct}" \
    "scaling-factor ratio output"

if [[ -s "${fil_out_scl_fct}" ]]; then
    assert_pattern_found \
        "${fil_out_scl_fct}" \
        $'^I\t0\t10\t4$' \
        "scaling-factor ratio output has I:0-10 = 4"

    assert_pattern_found \
        "${fil_out_scl_fct}" \
        $'^I\t40\t50\t8$' \
        "scaling-factor ratio output has I:40-50 = 8"

    assert_pattern_found \
        "${fil_out_scl_fct}" \
        $'^I\t60\t70\t0.667$' \
        "scaling-factor ratio output has I:60-70 = 0.667"
fi


#  Log2 ratio: log2(4 / 2) = 1 and log2(2 / 0.5) = 2
run_case_compute_signal_ratio \
    submit \
    "log2" \
    "${fil_out_log2}" \
    "${log_log2}" \
    "log2" \
    "${fil_A}" \
    "${fil_B}" \
    "${dir_out}" \
    "${dir_err}" \
    --csv_scl_fct \
    NA \
    --csv_dep_min \
    NA \
    --csv_pseudo \
    NA \
    --eps \
    0 \
    --skip_00 \
    NA \
    --dp \
    3

assert_file_nonempty \
    "${fil_out_log2}" \
    "log2 ratio output"

if [[ -s "${fil_out_log2}" ]]; then
    assert_pattern_found \
        "${fil_out_log2}" \
        $'^I\t0\t10\t1$' \
        "log2 ratio output has I:0-10 = 1"

    assert_pattern_found \
        "${fil_out_log2}" \
        $'^I\t40\t50\t2$' \
        "log2 ratio output has I:40-50 = 2"

    assert_pattern_found \
        "${fil_out_log2}" \
        $'^I\t60\t70\t-1.585$' \
        "log2 ratio output has I:60-70 = -1.585"
fi


#  Reciprocal ratio: B / A gives 0.5, 0.25, and 3 for selected rows
run_case_compute_signal_ratio \
    submit \
    "unadj_r" \
    "${fil_out_unadj_r}" \
    "${log_unadj_r}" \
    "unadj_r" \
    "${fil_A}" \
    "${fil_B}" \
    "${dir_out}" \
    "${dir_err}" \
    --csv_scl_fct \
    NA \
    --csv_dep_min \
    NA \
    --csv_pseudo \
    NA \
    --eps \
    0 \
    --skip_00 \
    NA \
    --dp \
    3

assert_file_nonempty \
    "${fil_out_unadj_r}" \
    "reciprocal ratio output"

if [[ -s "${fil_out_unadj_r}" ]]; then
    assert_pattern_found \
        "${fil_out_unadj_r}" \
        $'^I\t0\t10\t0.5$' \
        "reciprocal ratio output has I:0-10 = 0.5"

    assert_pattern_found \
        "${fil_out_unadj_r}" \
        $'^I\t40\t50\t0.25$' \
        "reciprocal ratio output has I:40-50 = 0.25"

    assert_pattern_found \
        "${fil_out_unadj_r}" \
        $'^I\t60\t70\t3$' \
        "reciprocal ratio output has I:60-70 = 3"
fi


#  Reciprocal log2 ratio: log2(B / A)
run_case_compute_signal_ratio \
    submit \
    "log2_r" \
    "${fil_out_log2_r}" \
    "${log_log2_r}" \
    "log2_r" \
    "${fil_A}" \
    "${fil_B}" \
    "${dir_out}" \
    "${dir_err}" \
    --csv_scl_fct \
    NA \
    --csv_dep_min \
    NA \
    --csv_pseudo \
    NA \
    --eps \
    0 \
    --skip_00 \
    NA \
    --dp \
    3

assert_file_nonempty \
    "${fil_out_log2_r}" \
    "reciprocal log2 ratio output"

if [[ -s "${fil_out_log2_r}" ]]; then
    assert_pattern_found \
        "${fil_out_log2_r}" \
        $'^I\t0\t10\t-1$' \
        "reciprocal log2 ratio output has I:0-10 = -1"

    assert_pattern_found \
        "${fil_out_log2_r}" \
        $'^I\t40\t50\t-2$' \
        "reciprocal log2 ratio output has I:40-50 = -2"

    assert_pattern_found \
        "${fil_out_log2_r}" \
        $'^I\t60\t70\t1.585$' \
        "reciprocal log2 ratio output has I:60-70 = 1.585"
fi


#  Denominator floor: B=0.04 is floored to 0.1, so 1 / 0.1 = 10
run_case_compute_signal_ratio \
    submit \
    "dep_min" \
    "${fil_out_dep_min}" \
    "${log_dep_min}" \
    "unadj" \
    "${fil_A}" \
    "${fil_B}" \
    "${dir_out}" \
    "${dir_err}" \
    --csv_scl_fct \
    NA \
    --csv_dep_min \
    0.1 \
    --csv_pseudo \
    NA \
    --eps \
    0 \
    --skip_00 \
    NA \
    --dp \
    3

assert_file_nonempty \
    "${fil_out_dep_min}" \
    "dep_min ratio output"

if [[ -s "${fil_out_dep_min}" ]]; then
    assert_pattern_found \
        "${fil_out_dep_min}" \
        $'^I\t50\t60\t10$' \
        "dep_min ratio output has I:50-60 = 10"
fi


#  Epsilon guards denominator values at or below eps: B=0.04 <= 0.05 -> nan
run_case_compute_signal_ratio \
    submit \
    "eps" \
    "${fil_out_eps}" \
    "${log_eps}" \
    "unadj" \
    "${fil_A}" \
    "${fil_B}" \
    "${dir_out}" \
    "${dir_err}" \
    --csv_scl_fct \
    NA \
    --csv_dep_min \
    NA \
    --csv_pseudo \
    NA \
    --eps \
    0.05 \
    --skip_00 \
    NA \
    --dp \
    3

assert_file_nonempty \
    "${fil_out_eps}" \
    "eps ratio output"

if [[ -s "${fil_out_eps}" ]]; then
    assert_pattern_found \
        "${fil_out_eps}" \
        $'^I\t0\t10\t2$' \
        "eps ratio output retains I:0-10 = 2"

    assert_pattern_found \
        "${fil_out_eps}" \
        $'^I\t50\t60\tnan$' \
        "eps ratio output has I:50-60 = nan"

    assert_pattern_found \
        "${fil_out_eps}" \
        $'^I\t60\t70\t0.333$' \
        "eps ratio output retains I:60-70 = 0.333"
fi


#  Pseudocounts: (0 + 1) / (2 + 1) = 0.333 at three decimals
run_case_compute_signal_ratio \
    submit \
    "pseudo" \
    "${fil_out_pseudo}" \
    "${log_pseudo}" \
    "unadj" \
    "${fil_A}" \
    "${fil_B}" \
    "${dir_out}" \
    "${dir_err}" \
    --csv_scl_fct \
    NA \
    --csv_dep_min \
    NA \
    --csv_pseudo \
    1:1 \
    --eps \
    0 \
    --skip_00 \
    NA \
    --dp \
    3

assert_file_nonempty \
    "${fil_out_pseudo}" \
    "pseudo ratio output"

if [[ -s "${fil_out_pseudo}" ]]; then
    assert_pattern_found \
        "${fil_out_pseudo}" \
        $'^I\t10\t20\t0.333$' \
        "pseudo ratio output has I:10-20 = 0.333"
fi


#  Drop non-finite rows while preserving finite ratio rows
run_case_compute_signal_ratio \
    submit \
    "drp_nan" \
    "${fil_out_drp_nan}" \
    "${log_drp_nan}" \
    "unadj" \
    "${fil_A}" \
    "${fil_B}" \
    "${dir_out}" \
    "${dir_err}" \
    --csv_scl_fct \
    NA \
    --csv_dep_min \
    NA \
    --csv_pseudo \
    NA \
    --eps \
    0 \
    --skip_00 \
    NA \
    --drp_nan \
    --dp \
    3

assert_file_nonempty \
    "${fil_out_drp_nan}" \
    "drp_nan ratio output"

if [[ -s "${fil_out_drp_nan}" ]]; then
    assert_pattern_found \
        "${fil_out_drp_nan}" \
        $'^I\t0\t10\t2$' \
        "drp_nan ratio output retains I:0-10 = 2"

    assert_pattern_found \
        "${fil_out_drp_nan}" \
        $'^I\t60\t70\t0.333$' \
        "drp_nan ratio output retains I:60-70 = 0.333"

    assert_pattern_absent \
        "${fil_out_drp_nan}" \
        $'^I\t20\t30\t' \
        "drp_nan ratio output omits I:20-30"

    assert_pattern_absent \
        "${fil_out_drp_nan}" \
        $'^I\t30\t40\t' \
        "drp_nan ratio output omits I:30-40"
fi


#  Zero-zero skipping before scaling removes the A=0, B=0 bin
run_case_compute_signal_ratio \
    submit \
    "skip_00" \
    "${fil_out_skip_00}" \
    "${log_skip_00}" \
    "unadj" \
    "${fil_A}" \
    "${fil_B}" \
    "${dir_out}" \
    "${dir_err}" \
    --csv_scl_fct \
    NA \
    --csv_dep_min \
    NA \
    --csv_pseudo \
    NA \
    --eps \
    0 \
    --skip_00 \
    pre_scale \
    --dp \
    3

assert_file_nonempty \
    "${fil_out_skip_00}" \
    "skip_00 ratio output"

if [[ -s "${fil_out_skip_00}" ]]; then
    assert_pattern_absent \
        "${fil_out_skip_00}" \
        $'^I\t30\t40\t' \
        "skip_00 ratio output omits I:30-40"
fi


#  Post-scale zero-zero skipping can remove bins that are non-zero before
#+ scaling: I:50-60 has A=1 and B=0.04, then scales to 0.001 and 0.004
run_case_compute_signal_ratio \
    submit \
    "skip_00_post_scale" \
    "${fil_out_skip_00_post}" \
    "${log_skip_00_post}" \
    "unadj" \
    "${fil_A}" \
    "${fil_B}" \
    "${dir_out}" \
    "${dir_err}" \
    --csv_scl_fct \
    0.001:0.1 \
    --csv_dep_min \
    NA \
    --csv_pseudo \
    NA \
    --eps \
    0.005 \
    --skip_00 \
    post_scale \
    --dp \
    3

assert_file_nonempty \
    "${fil_out_skip_00_post}" \
    "skip_00 post_scale ratio output"

if [[ -s "${fil_out_skip_00_post}" ]]; then
    assert_pattern_found \
        "${fil_out_skip_00_post}" \
        $'^I\t0\t10\t0.02$' \
        "skip_00 post_scale ratio output has I:0-10 = 0.02"

    assert_pattern_found \
        "${fil_out_skip_00_post}" \
        $'^I\t40\t50\t0.04$' \
        "skip_00 post_scale ratio output has I:40-50 = 0.04"

    assert_pattern_absent \
        "${fil_out_skip_00_post}" \
        $'^I\t50\t60\t' \
        "skip_00 post_scale ratio output omits scaled zero-zero I:50-60"
fi


#  Track sidecar should be generated and should omit non-finite rows
run_case_compute_signal_ratio \
    submit \
    "track" \
    "${fil_out_track}" \
    "${log_track}" \
    "unadj" \
    "${fil_A}" \
    "${fil_B}" \
    "${dir_out}" \
    "${dir_err}" \
    --csv_scl_fct \
    NA \
    --csv_dep_min \
    NA \
    --csv_pseudo \
    NA \
    --eps \
    0 \
    --skip_00 \
    NA \
    --track \
    --dp \
    3

assert_file_nonempty \
    "${fil_out_track}" \
    "track main ratio output"

assert_file_nonempty \
    "${trackfile_track}" \
    "track sidecar output"

if [[ -s "${trackfile_track}" ]]; then
    assert_pattern_found \
        "${trackfile_track}" \
        $'^I\t0\t10\t2$' \
        "track sidecar retains I:0-10 = 2"

    assert_pattern_found \
        "${trackfile_track}" \
        $'^I\t60\t70\t0.333$' \
        "track sidecar retains I:60-70 = 0.333"

    assert_pattern_absent \
        "${trackfile_track}" \
        $'^I\t20\t30\t' \
        "track sidecar omits I:20-30"

    assert_pattern_absent \
        "${trackfile_track}" \
        $'^I\t30\t40\t' \
        "track sidecar omits I:30-40"
fi


#  Legacy rounding alias: 1 / 3 rounds to 0.33 at two decimals
run_case_compute_signal_ratio \
    submit \
    "rnd_alias" \
    "${fil_out_rnd_alias}" \
    "${log_rnd_alias}" \
    "unadj" \
    "${fil_A}" \
    "${fil_B}" \
    "${dir_out}" \
    "${dir_err}" \
    --csv_scl_fct \
    NA \
    --csv_dep_min \
    NA \
    --csv_pseudo \
    NA \
    --eps \
    0 \
    --skip_00 \
    NA \
    --rnd \
    2

assert_file_nonempty \
    "${fil_out_rnd_alias}" \
    "rnd alias ratio output"

if [[ -s "${fil_out_rnd_alias}" ]]; then
    assert_pattern_found \
        "${fil_out_rnd_alias}" \
        $'^I\t60\t70\t0.33$' \
        "rnd alias ratio output has I:60-70 = 0.33"
fi


#  Gzipped bedGraph input and output should round-trip through ratio mode
# shellcheck disable=SC2154
if \
    run_capture \
        "submit compute-signal ratio gzip_io" \
        "${log_gzip_io}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/submit_compute_signal.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/scripts" \
            --threads 1 \
            --mode ratio \
            --method unadj \
            --csv_fil_A "${fil_A_gz}" \
            --csv_fil_B "${fil_B_gz}" \
            --csv_outfile "${fil_out_gzip_io}" \
            --err_out "${dir_err}" \
            --nam_job "test_compute_ratio_gzip_io" \
            --csv_scl_fct NA \
            --csv_dep_min NA \
            --csv_pseudo NA \
            --eps 0 \
            --skip_00 NA \
            --dp 3
then
    record_pass "submit_compute_signal.sh ratio gzip_io exits 0"
else
    record_fail \
        "submit_compute_signal.sh ratio gzip_io failed; see" \
        "$(print_relpath "${log_gzip_io}")"
fi

assert_file_nonempty \
    "${fil_out_gzip_io}" \
    "gzip ratio output"

if [[ -s "${fil_out_gzip_io}" ]]; then
    if \
        gzip -t "${fil_out_gzip_io}"
    then
        record_pass "gzip ratio output passes gzip integrity check"
        gzip -cd "${fil_out_gzip_io}" > "${outfile_txt_gzip_io}"
    else
        record_fail "gzip ratio output fails gzip integrity check"
    fi
fi

if [[ -s "${outfile_txt_gzip_io}" ]]; then
    assert_pattern_found \
        "${outfile_txt_gzip_io}" \
        $'^I\t0\t10\t2$' \
        "gzip ratio output has I:0-10 = 2"

    assert_pattern_found \
        "${outfile_txt_gzip_io}" \
        $'^I\t40\t50\t4$' \
        "gzip ratio output has I:40-50 = 4"

    assert_pattern_found \
        "${outfile_txt_gzip_io}" \
        $'^I\t60\t70\t0.333$' \
        "gzip ratio output has I:60-70 = 0.333"
fi


#  Header/prefix skipping should ignore default and custom metadata lines

# shellcheck disable=SC2154
if \
    run_capture \
        "submit compute-signal ratio skp_pfx" \
        "${log_skp_pfx}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/submit_compute_signal.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/scripts" \
            --threads 1 \
            --mode ratio \
            --method unadj \
            --csv_fil_A "${fil_A_hdr}" \
            --csv_fil_B "${fil_B_hdr}" \
            --csv_outfile "${fil_out_skp_pfx}" \
            --err_out "${dir_err}" \
            --nam_job "test_compute_ratio_skp_pfx" \
            --csv_scl_fct NA \
            --csv_dep_min NA \
            --csv_pseudo NA \
            --eps 0 \
            --skip_00 NA \
            --skp_pfx "#,track,browser,customHeader" \
            --dp 3
then
    record_pass "submit_compute_signal.sh ratio skp_pfx exits 0"
else
    record_fail \
        "submit_compute_signal.sh ratio skp_pfx failed; see" \
        "$(print_relpath "${log_skp_pfx}")"
fi

assert_file_nonempty \
    "${fil_out_skp_pfx}" \
    "skp_pfx ratio output"

if [[ -s "${fil_out_skp_pfx}" ]]; then
    assert_pattern_found \
        "${fil_out_skp_pfx}" \
        $'^I\t0\t10\t2$' \
        "skp_pfx ratio output has I:0-10 = 2"

    assert_pattern_found \
        "${fil_out_skp_pfx}" \
        $'^I\t40\t50\t4$' \
        "skp_pfx ratio output has I:40-50 = 4"

    assert_pattern_found \
        "${fil_out_skp_pfx}" \
        $'^I\t60\t70\t0.333$' \
        "skp_pfx ratio output has I:60-70 = 0.333"

    assert_pattern_absent \
        "${fil_out_skp_pfx}" \
        '^track' \
        "skp_pfx ratio output omits track headers"

    assert_pattern_absent \
        "${fil_out_skp_pfx}" \
        '^browser' \
        "skp_pfx ratio output omits browser headers"

    assert_pattern_absent \
        "${fil_out_skp_pfx}" \
        '^customHeader' \
        "skp_pfx ratio output omits custom headers"

    assert_pattern_absent \
        "${fil_out_skp_pfx}" \
        '^#' \
        "skp_pfx ratio output omits comment headers"
fi

finish
