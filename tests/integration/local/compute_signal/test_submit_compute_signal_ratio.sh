#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_submit_compute_signal_ratio.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6);
# - Anthropic Claude Code (Opus 5, Opus 5.5).
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="submit compute-signal ratio"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


dir_fx="${ROOT_REPO}/tests/fixtures/compute_signal/bedgraph"
chr_siz="${ROOT_REPO}/tests/fixtures/compute_signal/reference/tiny.fa.fai"
fil_A="${dir_fx}/ratio_A.bdg"
fil_B="${dir_fx}/ratio_B.bdg"
fil_A_hdr="${dir_fx}/ratio_headers_A.bdg"
fil_B_hdr="${dir_fx}/ratio_headers_B.bdg"
fil_A_gz="${dir_fx}/ratio_A.bdg.gz"
fil_B_gz="${dir_fx}/ratio_B.bdg.gz"

tmp="${TEST_DIR_TMP}/submit_compute_signal_ratio"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/compute_signal"

fil_out_linear="${dir_out}/ratio_linear.dp3.bdg"
fil_out_scl_fct="${dir_out}/ratio_scl_fct_2_1.dp3.bdg"
fil_out_log2="${dir_out}/ratio_log2.dp3.bdg"
fil_out_linear_r="${dir_out}/ratio_linear_r.dp3.bdg"
fil_out_log2_r="${dir_out}/ratio_log2_r.dp3.bdg"
fil_out_dep_min="${dir_out}/ratio_dep_min_0p1.dp3.bdg"
fil_out_eps="${dir_out}/ratio_eps_0p05.dp3.bdg"
fil_out_pseudo="${dir_out}/ratio_pseudo_1_1.dp3.bdg"
fil_out_drp_nan="${dir_out}/ratio_drp_nan.dp3.bdg"
fil_out_skip_00="${dir_out}/ratio_skip_00_pre_scale.dp3.bdg"
fil_out_skip_00_post="${dir_out}/ratio_skip_00_post_scale.dp3.bdg"
fil_out_track="${dir_out}/ratio_track.dp3.bdg"
fil_out_dp_alias="${dir_out}/ratio_dp_alias.dp2.bdg"
fil_out_gzip_io="${dir_out}/ratio_gzip_io.dp3.bdg.gz"
fil_out_skp_pfx="${dir_out}/ratio_skp_pfx.dp3.bdg"

trackfile_track="${dir_out}/ratio_track.dp3.track.bdg"
outfile_txt_gzip_io="${dir_out}/ratio_gzip_io.dp3.bdg"

log_linear="${dir_log}/submit_compute_signal_ratio_linear.log"
log_scl_fct="${dir_log}/submit_compute_signal_ratio_scl_fct.log"
log_log2="${dir_log}/submit_compute_signal_ratio_log2.log"
log_linear_r="${dir_log}/submit_compute_signal_ratio_linear_r.log"
log_log2_r="${dir_log}/submit_compute_signal_ratio_log2_r.log"
log_dep_min="${dir_log}/submit_compute_signal_ratio_dep_min.log"
log_eps="${dir_log}/submit_compute_signal_ratio_eps.log"
log_pseudo="${dir_log}/submit_compute_signal_ratio_pseudo.log"
log_drp_nan="${dir_log}/submit_compute_signal_ratio_drp_nan.log"
log_skip_00="${dir_log}/submit_compute_signal_ratio_skip_00.log"
log_skip_00_post="${dir_log}/submit_compute_signal_ratio_skip_00_post_scale.log"
log_track="${dir_log}/submit_compute_signal_ratio_track.log"
log_dp_alias="${dir_log}/submit_compute_signal_ratio_dp_alias.log"
log_gzip_io="${dir_log}/submit_compute_signal_ratio_gzip_io.log"
log_skp_pfx="${dir_log}/submit_compute_signal_ratio_skp_pfx.log"


print_section "${TEST_NAME}"


rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty \
    "${chr_siz}" \
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

require_files_nonempty \
    "${fil_A_gz}" \
    "${fil_B_gz}" || {
    finish
    exit $?
}


# Baseline linear ratio with three-decimal rounding.
run_case_compute_signal_ratio \
    submit \
    "linear" \
    "${fil_out_linear}" \
    "${log_linear}" \
    "linear" \
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
    --chr_siz \
    "${chr_siz}" \
    --strict_bins \
    --dp \
    3

assert_file_nonempty \
    "${fil_out_linear}" \
    "baseline ratio output"

if [[ -s "${fil_out_linear}" ]]; then
    assert_pattern_found \
        "${fil_out_linear}" \
        $'^I\t0\t10\t2$' \
        "baseline ratio output has I:0-10 = 2"

    assert_pattern_found \
        "${fil_out_linear}" \
        $'^I\t10\t20\t0$' \
        "baseline ratio output has I:10-20 = 0"

    assert_pattern_found \
        "${fil_out_linear}" \
        $'^I\t40\t50\t4$' \
        "baseline ratio output has I:40-50 = 4"

    assert_pattern_found \
        "${fil_out_linear}" \
        $'^I\t60\t70\t0.333$' \
        "baseline ratio output has I:60-70 = 0.333"

    assert_pattern_found \
        "${fil_out_linear}" \
        $'^I\t70\t80\t1$' \
        "baseline ratio output has I:70-80 = 1"
fi


# Scaling factors are applied before ratio calculation: (2 * A) / (1 * B).
run_case_compute_signal_ratio \
    submit \
    "scl_fct" \
    "${fil_out_scl_fct}" \
    "${log_scl_fct}" \
    "linear" \
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


# Log2 ratio: log2(4 / 2) = 1 and log2(2 / 0.5) = 2.
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


# Reciprocal ratio: B / A gives 0.5, 0.25, and 3 for selected rows.
run_case_compute_signal_ratio \
    submit \
    "linear_r" \
    "${fil_out_linear_r}" \
    "${log_linear_r}" \
    "linear_r" \
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
    "${fil_out_linear_r}" \
    "reciprocal ratio output"

if [[ -s "${fil_out_linear_r}" ]]; then
    assert_pattern_found \
        "${fil_out_linear_r}" \
        $'^I\t0\t10\t0.5$' \
        "reciprocal ratio output has I:0-10 = 0.5"

    assert_pattern_found \
        "${fil_out_linear_r}" \
        $'^I\t40\t50\t0.25$' \
        "reciprocal ratio output has I:40-50 = 0.25"

    assert_pattern_found \
        "${fil_out_linear_r}" \
        $'^I\t60\t70\t3$' \
        "reciprocal ratio output has I:60-70 = 3"
fi


# Reciprocal log2 ratio: log2(B / A).
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


# Denominator floor: B=0.04 is floored to 0.1, so 1 / 0.1 = 10.
run_case_compute_signal_ratio \
    submit \
    "dep_min" \
    "${fil_out_dep_min}" \
    "${log_dep_min}" \
    "linear" \
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


# Epsilon guards denominator values at or below eps: B=0.04 <= 0.05 -> nan.
run_case_compute_signal_ratio \
    submit \
    "eps" \
    "${fil_out_eps}" \
    "${log_eps}" \
    "linear" \
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


# Pseudocounts: (0 + 1) / (2 + 1) = 0.333 at three decimals.
run_case_compute_signal_ratio \
    submit \
    "pseudo" \
    "${fil_out_pseudo}" \
    "${log_pseudo}" \
    "linear" \
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


# Drop non-finite rows while preserving finite ratio rows.
run_case_compute_signal_ratio \
    submit \
    "drp_nan" \
    "${fil_out_drp_nan}" \
    "${log_drp_nan}" \
    "linear" \
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


# Zero-zero skipping before scaling removes the A=0, B=0 bin.
run_case_compute_signal_ratio \
    submit \
    "skip_00" \
    "${fil_out_skip_00}" \
    "${log_skip_00}" \
    "linear" \
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


# Post-scale zero-zero skipping can remove bins that are non-zero before
# scaling: I:50-60 has A=1 and B=0.04, then scales to 0.001 and 0.004.
run_case_compute_signal_ratio \
    submit \
    "skip_00_post_scale" \
    "${fil_out_skip_00_post}" \
    "${log_skip_00_post}" \
    "linear" \
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


# Track sidecar should be generated and should omit non-finite rows.
run_case_compute_signal_ratio \
    submit \
    "track" \
    "${fil_out_track}" \
    "${log_track}" \
    "linear" \
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


# Legacy rounding alias: 1 / 3 rounds to 0.33 at two decimals.
run_case_compute_signal_ratio \
    submit \
    "dp_precision" \
    "${fil_out_dp_alias}" \
    "${log_dp_alias}" \
    "linear" \
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
    2

assert_file_nonempty \
    "${fil_out_dp_alias}" \
    "dp alias ratio output"

if [[ -s "${fil_out_dp_alias}" ]]; then
    assert_pattern_found \
        "${fil_out_dp_alias}" \
        $'^I\t60\t70\t0.33$' \
        "dp alias ratio output has I:60-70 = 0.33"
fi


# Gzipped bedGraph input and output should round-trip through ratio mode.
# shellcheck disable=SC2154
if \
    run_capture \
        "submit compute-signal ratio gzip_io" \
        "${log_gzip_io}" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/submit_compute_signal.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/bin" \
            --threads 1 \
            --mode ratio \
            --method linear \
            --csv_fil_A "${fil_A_gz}" \
            --csv_fil_B "${fil_B_gz}" \
            --csv_fil_out "${fil_out_gzip_io}" \
            --dir_eo "${dir_err}" \
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

assert_file_nonempty \
    "${outfile_txt_gzip_io}" \
    "submit decompressed gzip ratio output"

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


# Header/prefix skipping should ignore default and custom metadata lines.
# shellcheck disable=SC2154
if \
    run_capture \
        "submit compute-signal ratio skp_pfx" \
        "${log_skp_pfx}" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/submit_compute_signal.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/bin" \
            --threads 1 \
            --mode ratio \
            --method linear \
            --csv_fil_A "${fil_A_hdr}" \
            --csv_fil_B "${fil_B_hdr}" \
            --csv_fil_out "${fil_out_skp_pfx}" \
            --dir_eo "${dir_err}" \
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


# Ratio runs must never receive the signal-only window options; forwarding is
# covered by the execute suite.
log_rat_mode="${tmp}/logs/test_compute_ratio_skp_pfx.ratio_skp_pfx.dp3.stderr.txt"

assert_file_nonempty \
    "${log_rat_mode}" \
    "submit ratio stderr log for mode separation"

if [[ -s "${log_rat_mode}" ]]; then
    assert_pattern_absent \
        "${log_rat_mode}" \
        "--siz_win" \
        "submit ratio omits '--siz_win'"

    assert_pattern_absent \
        "${log_rat_mode}" \
        "--engine" \
        "submit ratio omits '--engine'"
fi


# '--prior_count' is refused wherever it could change no output. A direct
# submit run must refuse it too, not only the execute driver.
rc_pc=0
out_pc="$(
    "${TEST_BASH}" "${ROOT_REPO}/bin/submit_compute_signal.sh" \
        --env_nam "${env_nam}" \
        --dir_scr "${ROOT_REPO}/bin" \
        --threads 1 \
        --mode ratio \
        --method linear \
        --prior_count 1 \
        --csv_fil_A "${fil_A}" \
        --csv_fil_B "${fil_B}" \
        --csv_fil_out "${tmp}/ratio_pc_lit.bedGraph" \
        --dir_eo "${dir_err}" \
        --nam_job "test_compute_ratio_pc_lit" \
        --csv_scl_fct NA \
        --csv_dep_min NA \
        --csv_pseudo 1:1 2>&1
)" || rc_pc=$?

if [[
    "${rc_pc}" -ne 0
    && "${out_pc}" == *"'--prior_count' applies only to '--csv_pseudo edger'"*
    && ! -e "${tmp}/ratio_pc_lit.bedGraph"
]]; then
    record_pass "submit rejects '--prior_count' when no element is 'edger'"
else
    record_fail \
        "submit did not reject '--prior_count' with only literal" \
        "pseudocounts (exit ${rc_pc})"
fi

# With '--csv_pseudo' omitted, every element is the 'NA' sentinel.
rc_pc=0
out_pc="$(
    "${TEST_BASH}" "${ROOT_REPO}/bin/submit_compute_signal.sh" \
        --env_nam "${env_nam}" \
        --dir_scr "${ROOT_REPO}/bin" \
        --threads 1 \
        --mode ratio \
        --method linear \
        --prior_count 1 \
        --csv_fil_A "${fil_A}" \
        --csv_fil_B "${fil_B}" \
        --csv_fil_out "${tmp}/ratio_pc_none.bedGraph" \
        --dir_eo "${dir_err}" \
        --nam_job "test_compute_ratio_pc_none" 2>&1
)" || rc_pc=$?

if [[
    "${rc_pc}" -ne 0
    && "${out_pc}" == *"'--prior_count' applies only to '--csv_pseudo edger'"*
    && ! -e "${tmp}/ratio_pc_none.bedGraph"
]]; then
    record_pass "submit rejects '--prior_count' when '--csv_pseudo' is omitted"
else
    record_fail \
        "submit did not reject '--prior_count' without '--csv_pseudo'" \
        "(exit ${rc_pc})"
fi

bam_se="${ROOT_REPO}/tests/fixtures/compute_signal/bam/se/tiny_se.bam"

rc_pc=0
out_pc="$(
    "${TEST_BASH}" "${ROOT_REPO}/bin/submit_compute_signal.sh" \
        --env_nam "${env_nam}" \
        --dir_scr "${ROOT_REPO}/bin" \
        --threads 1 \
        --mode signal \
        --method norm \
        --prior_count 1 \
        --csv_fil_in "${bam_se}" \
        --csv_fil_out "${tmp}/signal_pc.bedGraph" \
        --dir_eo "${dir_err}" \
        --nam_job "test_compute_signal_pc" 2>&1
)" || rc_pc=$?

if [[
    "${rc_pc}" -ne 0
    && "${out_pc}" == *"'--prior_count' is for '--mode ratio'"*
]]; then
    record_pass "submit rejects '--prior_count' in signal mode"
else
    record_fail \
        "submit did not reject '--prior_count' in signal mode" \
        "(exit ${rc_pc})"
fi
unset bam_se rc_pc out_pc

# A bare literal would regularize file A alone and leave file B's zero bins
# undefined, so a direct submit run refuses it too.
rc_bare=0
out_bare="$(
    "${TEST_BASH}" "${ROOT_REPO}/bin/submit_compute_signal.sh" \
        --env_nam "${env_nam}" \
        --dir_scr "${ROOT_REPO}/bin" \
        --threads 1 \
        --mode ratio \
        --method linear \
        --csv_fil_A "${fil_A}" \
        --csv_fil_B "${fil_B}" \
        --csv_fil_out "${tmp}/ratio_bare.bedGraph" \
        --dir_eo "${dir_err}" \
        --nam_job "test_compute_ratio_bare" \
        --csv_pseudo 1 2>&1
)" || rc_bare=$?

if [[
    "${rc_bare}" -ne 0
    && "${out_bare}" == *"needs 'A:B'"*
    && ! -e "${tmp}/ratio_bare.bedGraph"
]]; then
    record_pass "submit refuses a bare '--csv_pseudo' literal, naming 'A:B'"
else
    record_fail \
        "submit did not refuse a bare '--csv_pseudo' literal" \
        "(exit ${rc_bare})"
fi
unset rc_bare out_bare

# '--prior_count' belongs to ratio mode, so a direct submit run in coord mode
# refuses it too.
bam_crd="${ROOT_REPO}/tests/fixtures/compute_signal/bam/se/tiny_se.bam"

rc_crd=0
out_crd="$(
    "${TEST_BASH}" "${ROOT_REPO}/bin/submit_compute_signal.sh" \
        --env_nam "${env_nam}" \
        --dir_scr "${ROOT_REPO}/bin" \
        --threads 1 \
        --mode coord \
        --prior_count 1 \
        --csv_fil_in "${bam_crd}" \
        --csv_fil_out "${tmp}/coord_pc.bed" \
        --dir_eo "${dir_err}" \
        --nam_job "test_compute_coord_pc" 2>&1
)" || rc_crd=$?

if [[
    "${rc_crd}" -ne 0
    && "${out_crd}" == *"'--prior_count' is for '--mode ratio'"*
    && ! -e "${tmp}/coord_pc.bed"
]]; then
    record_pass "submit rejects '--prior_count' in coord mode"
else
    record_fail \
        "submit did not reject '--prior_count' in coord mode (exit ${rc_crd})"
fi
unset bam_crd rc_crd out_crd

# A non-numeric '--prior_count' is refused before any work.
rc_nan=0
out_nan="$(
    "${TEST_BASH}" "${ROOT_REPO}/bin/submit_compute_signal.sh" \
        --env_nam "${env_nam}" \
        --dir_scr "${ROOT_REPO}/bin" \
        --threads 1 \
        --mode ratio \
        --method log2 \
        --typ_sig count \
        --prior_count two \
        --csv_fil_A "${fil_A}" \
        --csv_fil_B "${fil_B}" \
        --csv_fil_out "${tmp}/ratio_pc_nan.bedGraph" \
        --dir_eo "${dir_err}" \
        --nam_job "test_compute_ratio_pc_nan" \
        --csv_pseudo edger 2>&1
)" || rc_nan=$?

if [[
    "${rc_nan}" -ne 0
    && "${out_nan}" == *"'--prior_count' was assigned 'two'"*
    && ! -e "${tmp}/ratio_pc_nan.bedGraph"
]]; then
    record_pass "submit refuses a non-numeric '--prior_count' before any work"
else
    record_fail \
        "submit did not refuse a non-numeric '--prior_count'" \
        "(exit ${rc_nan})"
fi
unset rc_nan out_nan


# A mixed list given to submit whole, as under Slurm: each task must resolve
# only its own element.
# shellcheck source=lib/bash/core/format_outputs.sh
source "${ROOT_REPO}/lib/bash/core/format_outputs.sh"

dir_cnt="${ROOT_REPO}/tests/fixtures/compute_signal/bedgraph/count"
trk_se="${dir_cnt}/tiny_se.bedGraph"
trk_pe="${dir_cnt}/tiny_pe.bedGraph"
dir_mix="${tmp}/mixed"
out_mix_1="${dir_mix}/edg_se_pe.bedGraph"
out_mix_2="${dir_mix}/lit_pe_se.bedGraph"
mkdir -p "${dir_mix}"

exp_mix_1="$(
    "${TEST_MANAGED_PYTHON}" \
        -m protocol_chipseq_signal_norm.cli.compute_pseudo \
        --method edger \
        --typ_sig count \
        --prior_count 1 \
        --fil_A "${trk_se}" \
        --fil_B "${trk_pe}" \
        --n_frg_A "$(< "$(derive_report_path "${trk_se}" n_frg)")" \
        --n_frg_B "$(< "$(derive_report_path "${trk_pe}" n_frg)")" \
        --n_ovlp_A "$(< "$(derive_report_path "${trk_se}" n_ovlp)")" \
        --n_ovlp_B "$(< "$(derive_report_path "${trk_pe}" n_ovlp)")" \
        2>/dev/null
)"

rc_mix=0
"${TEST_BASH}" "${ROOT_REPO}/bin/submit_compute_signal.sh" \
    --env_nam "${env_nam}" \
    --dir_scr "${ROOT_REPO}/bin" \
    --threads 1 \
    --mode ratio \
    --method log2 \
    --typ_sig count \
    --prior_count 1 \
    --csv_fil_A "${trk_se},${trk_pe}" \
    --csv_fil_B "${trk_pe},${trk_se}" \
    --csv_fil_out "${out_mix_1},${out_mix_2}" \
    --dir_eo "${dir_mix}" \
    --nam_job "test_compute_ratio_mixed" \
    --csv_scl_fct NA,NA \
    --csv_dep_min NA,NA \
    --csv_pseudo "edger,1:1" > /dev/null 2>&1 || rc_mix=$?

# Compare each pair's logged '--pseudo' (printed as floats) to its spec
# numerically, so a value applied to the wrong pair cannot pass.
function call_has_pseudo() {
    local fil_out="${1}" spec="${2}"
    local log="${dir_mix}/test_compute_ratio_mixed.$(
        basename "${fil_out}" .bedGraph
    ).stderr.txt"

    [[ -s "${log}" ]] || return 1

    "${TEST_MANAGED_PYTHON}" - "${log}" "${spec}" << 'PY'
import math
import re
import sys

text = open(sys.argv[1], encoding="utf-8").read()
found = re.search(r"^--pseudo\s+(\S+):(\S+)$", text, re.MULTILINE)
want = [float(v) for v in sys.argv[2].split(":")]
got = [float(v) for v in found.groups()] if found else []
same = len(got) == 2 and all(
    math.isclose(g, w, rel_tol=1e-12) for g, w in zip(got, want)
)
raise SystemExit(0 if same else 1)
PY
}

if [[
    "${rc_mix}" -eq 0
    && -n "${exp_mix_1}"
    && -s "${out_mix_1}"
    && -s "${out_mix_2}"
]] \
    && call_has_pseudo "${out_mix_1}" "${exp_mix_1}" \
    && call_has_pseudo "${out_mix_2}" "1:1"
then
    record_pass \
        "submit given 'edger,1:1' whole applies '--prior_count' to the" \
        "'edger' pair only"
else
    record_fail \
        "submit given 'edger,1:1' whole did not give each pair its own" \
        "element (exit ${rc_mix}); see $(print_relpath "${dir_mix}")"
fi
unset dir_cnt trk_se trk_pe dir_mix out_mix_1 out_mix_2 exp_mix_1 rc_mix

finish
