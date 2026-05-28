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

print_section "${TEST_NAME}"

#  Define fixture and output paths for the simplest wet ratio-mode path
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
fil_out="${dir_out}/ratio_unadj.dp3.bdg"
log="${dir_log}/submit_compute_signal_ratio_unadj.log"

run_case_compute_signal_ratio \
    submit \
    "unadj" \
    "${fil_out}" \
    "${log}" \
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
    "${fil_out}" \
    "baseline ratio output"

if [[ -s "${fil_out}" ]]; then
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t0\t10\t2$' \
        "baseline ratio output has I:0-10 = 2"
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t10\t20\t0$' \
        "baseline ratio output has I:10-20 = 0"
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t40\t50\t4$' \
        "baseline ratio output has I:40-50 = 4"
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t60\t70\t0.333$' \
        "baseline ratio output has I:60-70 = 0.333"
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t70\t80\t1$' \
        "baseline ratio output has I:70-80 = 1"
fi


#  Scaling factors are applied before ratio calculation: (2 * A) / (1 * B)
fil_out="${dir_out}/ratio_scl_fct_2_1.dp3.bdg"
log="${dir_log}/submit_compute_signal_ratio_scl_fct.log"

run_case_compute_signal_ratio \
    submit \
    "scl_fct" \
    "${fil_out}" \
    "${log}" \
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
    "${fil_out}" \
    "scaling-factor ratio output"

if [[ -s "${fil_out}" ]]; then
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t0\t10\t4$' \
        "scaling-factor ratio output has I:0-10 = 4"
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t40\t50\t8$' \
        "scaling-factor ratio output has I:40-50 = 8"
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t60\t70\t0.667$' \
        "scaling-factor ratio output has I:60-70 = 0.667"
fi


#  Log2 ratio: log2(4 / 2) = 1 and log2(2 / 0.5) = 2
fil_out="${dir_out}/ratio_log2.dp3.bdg"
log="${dir_log}/submit_compute_signal_ratio_log2.log"

run_case_compute_signal_ratio \
    submit \
    "log2" \
    "${fil_out}" \
    "${log}" \
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
    "${fil_out}" \
    "log2 ratio output"

if [[ -s "${fil_out}" ]]; then
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t0\t10\t1$' \
        "log2 ratio output has I:0-10 = 1"
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t40\t50\t2$' \
        "log2 ratio output has I:40-50 = 2"
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t60\t70\t-1.585$' \
        "log2 ratio output has I:60-70 = -1.585"
fi


#  Reciprocal ratio: B / A gives 0.5, 0.25, and 3 for selected rows
fil_out="${dir_out}/ratio_unadj_r.dp3.bdg"
log="${dir_log}/submit_compute_signal_ratio_unadj_r.log"

run_case_compute_signal_ratio \
    submit \
    "unadj_r" \
    "${fil_out}" \
    "${log}" \
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
    "${fil_out}" \
    "reciprocal ratio output"

if [[ -s "${fil_out}" ]]; then
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t0\t10\t0.5$' \
        "reciprocal ratio output has I:0-10 = 0.5"
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t40\t50\t0.25$' \
        "reciprocal ratio output has I:40-50 = 0.25"
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t60\t70\t3$' \
        "reciprocal ratio output has I:60-70 = 3"
fi


#  Reciprocal log2 ratio: log2(B / A)
fil_out="${dir_out}/ratio_log2_r.dp3.bdg"
log="${dir_log}/submit_compute_signal_ratio_log2_r.log"

run_case_compute_signal_ratio \
    submit \
    "log2_r" \
    "${fil_out}" \
    "${log}" \
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
    "${fil_out}" \
    "reciprocal log2 ratio output"

if [[ -s "${fil_out}" ]]; then
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t0\t10\t-1$' \
        "reciprocal log2 ratio output has I:0-10 = -1"
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t40\t50\t-2$' \
        "reciprocal log2 ratio output has I:40-50 = -2"
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t60\t70\t1.585$' \
        "reciprocal log2 ratio output has I:60-70 = 1.585"
fi


#  Denominator floor: B=0.04 is floored to 0.1, so 1 / 0.1 = 10
fil_out="${dir_out}/ratio_dep_min_0p1.dp3.bdg"
log="${dir_log}/submit_compute_signal_ratio_dep_min.log"

run_case_compute_signal_ratio \
    submit \
    "dep_min" \
    "${fil_out}" \
    "${log}" \
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
    "${fil_out}" \
    "dep_min ratio output"

if [[ -s "${fil_out}" ]]; then
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t50\t60\t10$' \
        "dep_min ratio output has I:50-60 = 10"
fi


#  Epsilon guards denominator values at or below eps: B=0.04 <= 0.05 -> nan
fil_out="${dir_out}/ratio_eps_0p05.dp3.bdg"
log="${dir_log}/submit_compute_signal_ratio_eps.log"

run_case_compute_signal_ratio \
    submit \
    "eps" \
    "${fil_out}" \
    "${log}" \
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
    "${fil_out}" \
    "eps ratio output"

if [[ -s "${fil_out}" ]]; then
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t0\t10\t2$' \
        "eps ratio output retains I:0-10 = 2"
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t50\t60\tnan$' \
        "eps ratio output has I:50-60 = nan"
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t60\t70\t0.333$' \
        "eps ratio output retains I:60-70 = 0.333"
fi


#  Pseudocounts: (0 + 1) / (2 + 1) = 0.333 at three decimals
fil_out="${dir_out}/ratio_pseudo_1_1.dp3.bdg"
log="${dir_log}/submit_compute_signal_ratio_pseudo.log"

run_case_compute_signal_ratio \
    submit \
    "pseudo" \
    "${fil_out}" \
    "${log}" \
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
    "${fil_out}" \
    "pseudo ratio output"

if [[ -s "${fil_out}" ]]; then
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t10\t20\t0.333$' \
        "pseudo ratio output has I:10-20 = 0.333"
fi


#  Drop non-finite rows while preserving finite ratio rows
fil_out="${dir_out}/ratio_drp_nan.dp3.bdg"
log="${dir_log}/submit_compute_signal_ratio_drp_nan.log"

run_case_compute_signal_ratio \
    submit \
    "drp_nan" \
    "${fil_out}" \
    "${log}" \
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
    "${fil_out}" \
    "drp_nan ratio output"

if [[ -s "${fil_out}" ]]; then
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t0\t10\t2$' \
        "drp_nan ratio output retains I:0-10 = 2"
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t60\t70\t0.333$' \
        "drp_nan ratio output retains I:60-70 = 0.333"
    assert_pattern_absent \
        "${fil_out}" \
        $'^I\t20\t30\t' \
        "drp_nan ratio output omits I:20-30"
    assert_pattern_absent \
        "${fil_out}" \
        $'^I\t30\t40\t' \
        "drp_nan ratio output omits I:30-40"
fi


#  Zero-zero skipping before scaling removes the A=0, B=0 bin
fil_out="${dir_out}/ratio_skip_00_pre_scale.dp3.bdg"
log="${dir_log}/submit_compute_signal_ratio_skip_00.log"

run_case_compute_signal_ratio \
    submit \
    "skip_00" \
    "${fil_out}" \
    "${log}" \
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
    "${fil_out}" \
    "skip_00 ratio output"

if [[ -s "${fil_out}" ]]; then
    assert_pattern_absent \
        "${fil_out}" \
        $'^I\t30\t40\t' \
        "skip_00 ratio output omits I:30-40"
fi


#  Post-scale zero-zero skipping can remove bins that are non-zero before
#+ scaling: I:50-60 has A=1 and B=0.04, then scales to 0.001 and 0.004
fil_out="${dir_out}/ratio_skip_00_post_scale.dp3.bdg"
log="${dir_log}/submit_compute_signal_ratio_skip_00_post_scale.log"

run_case_compute_signal_ratio \
    submit \
    "skip_00_post_scale" \
    "${fil_out}" \
    "${log}" \
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
    "${fil_out}" \
    "skip_00 post_scale ratio output"

if [[ -s "${fil_out}" ]]; then
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t0\t10\t0.02$' \
        "skip_00 post_scale ratio output has I:0-10 = 0.02"
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t40\t50\t0.04$' \
        "skip_00 post_scale ratio output has I:40-50 = 0.04"
    assert_pattern_absent \
        "${fil_out}" \
        $'^I\t50\t60\t' \
        "skip_00 post_scale ratio output omits scaled zero-zero I:50-60"
fi


#  Track sidecar should be generated and should omit non-finite rows
fil_out="${dir_out}/ratio_track.dp3.bdg"
trackfile="${dir_out}/ratio_track.dp3.track.bdg"
log="${dir_log}/submit_compute_signal_ratio_track.log"

run_case_compute_signal_ratio \
    submit \
    "track" \
    "${fil_out}" \
    "${log}" \
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
    "${fil_out}" \
    "track main ratio output"
assert_file_nonempty \
    "${trackfile}" \
    "track sidecar output"

if [[ -s "${trackfile}" ]]; then
    assert_pattern_found \
        "${trackfile}" \
        $'^I\t0\t10\t2$' \
        "track sidecar retains I:0-10 = 2"
    assert_pattern_found \
        "${trackfile}" \
        $'^I\t60\t70\t0.333$' \
        "track sidecar retains I:60-70 = 0.333"
    assert_pattern_absent \
        "${trackfile}" \
        $'^I\t20\t30\t' \
        "track sidecar omits I:20-30"
    assert_pattern_absent \
        "${trackfile}" \
        $'^I\t30\t40\t' \
        "track sidecar omits I:30-40"
fi


#  Legacy rounding alias: 1 / 3 rounds to 0.33 at two decimals
fil_out="${dir_out}/ratio_rnd_alias.dp2.bdg"
log="${dir_log}/submit_compute_signal_ratio_rnd_alias.log"

run_case_compute_signal_ratio \
    submit \
    "rnd_alias" \
    "${fil_out}" \
    "${log}" \
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
    "${fil_out}" \
    "rnd alias ratio output"

if [[ -s "${fil_out}" ]]; then
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t60\t70\t0.33$' \
        "rnd alias ratio output has I:60-70 = 0.33"
fi


#  Gzipped bedGraph input and output should round-trip through ratio mode
fil_out="${dir_out}/ratio_gzip_io.dp3.bdg.gz"
outfile_txt="${dir_out}/ratio_gzip_io.dp3.bdg"
log="${dir_log}/submit_compute_signal_ratio_gzip_io.log"

# shellcheck disable=SC2154
if \
    run_capture \
        "submit compute-signal ratio gzip_io" \
        "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/submit_compute_signal.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/scripts" \
            --threads 1 \
            --mode ratio \
            --method unadj \
            --csv_fil_A "${fil_A_gz}" \
            --csv_fil_B "${fil_B_gz}" \
            --csv_outfile "${fil_out}" \
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
        "$(print_relpath "${log}")"
fi

assert_file_nonempty \
    "${fil_out}" \
    "gzip ratio output"

if [[ -s "${fil_out}" ]]; then
    if \
        gzip -t "${fil_out}"
    then
        record_pass "gzip ratio output passes gzip integrity check"
        gzip -cd "${fil_out}" > "${outfile_txt}"
    else
        record_fail "gzip ratio output fails gzip integrity check"
    fi
fi

if [[ -s "${outfile_txt}" ]]; then
    assert_pattern_found \
        "${outfile_txt}" \
        $'^I\t0\t10\t2$' \
        "gzip ratio output has I:0-10 = 2"
    assert_pattern_found \
        "${outfile_txt}" \
        $'^I\t40\t50\t4$' \
        "gzip ratio output has I:40-50 = 4"
    assert_pattern_found \
        "${outfile_txt}" \
        $'^I\t60\t70\t0.333$' \
        "gzip ratio output has I:60-70 = 0.333"
fi


#  Header/prefix skipping should ignore default and custom metadata lines
fil_out="${dir_out}/ratio_skp_pfx.dp3.bdg"
log="${dir_log}/submit_compute_signal_ratio_skp_pfx.log"

# shellcheck disable=SC2154
if \
    run_capture \
        "submit compute-signal ratio skp_pfx" \
        "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/submit_compute_signal.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/scripts" \
            --threads 1 \
            --mode ratio \
            --method unadj \
            --csv_fil_A "${fil_A_hdr}" \
            --csv_fil_B "${fil_B_hdr}" \
            --csv_outfile "${fil_out}" \
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
        "$(print_relpath "${log}")"
fi

assert_file_nonempty \
    "${fil_out}" \
    "skp_pfx ratio output"

if [[ -s "${fil_out}" ]]; then
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t0\t10\t2$' \
        "skp_pfx ratio output has I:0-10 = 2"
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t40\t50\t4$' \
        "skp_pfx ratio output has I:40-50 = 4"
    assert_pattern_found \
        "${fil_out}" \
        $'^I\t60\t70\t0.333$' \
        "skp_pfx ratio output has I:60-70 = 0.333"
    assert_pattern_absent \
        "${fil_out}" \
        '^track' \
        "skp_pfx ratio output omits track headers"
    assert_pattern_absent \
        "${fil_out}" \
        '^browser' \
        "skp_pfx ratio output omits browser headers"
    assert_pattern_absent \
        "${fil_out}" \
        '^customHeader' \
        "skp_pfx ratio output omits custom headers"
    assert_pattern_absent \
        "${fil_out}" \
        '^#' \
        "skp_pfx ratio output omits comment headers"
fi

finish
