#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_compute_signal_ratio.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="execute compute-signal ratio"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

rec_section "${TEST_NAME}"

#  Define fixture and output paths for the execute-to-submit ratio path
dir_fx="${ROOT_REPO}/tests/compute_signal/fixtures/bedgraph"
fil_A="${dir_fx}/ratio_A.bdg"
fil_B="${dir_fx}/ratio_B.bdg"
fil_hdr_A="${dir_fx}/ratio_headers_A.bdg"
fil_hdr_B="${dir_fx}/ratio_headers_B.bdg"

tmp="${TEST_DIR_TMP}/execute_compute_signal_ratio"
dir_in="${tmp}/in"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/compute_signal"
fil_gz_A="${dir_in}/ratio_A.bdg.gz"
fil_gz_B="${dir_in}/ratio_B.bdg.gz"

rm -rf "${tmp}"
mkdir -p "${dir_in}" "${dir_out}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_exist "${fil_A}" "${fil_B}" "${fil_hdr_A}" "${fil_hdr_B}" || {
    finish
    exit $?
}

if ! check_cmd_exists gzip; then
    rec_fail "gzip is required for ratio gzip I/O coverage"
    finish
    exit $?
fi

gzip -c "${fil_A}" > "${fil_gz_A}"
gzip -c "${fil_B}" > "${fil_gz_B}"

require_files_exist "${fil_gz_A}" "${fil_gz_B}" || {
    finish
    exit $?
}


#  Run a local serial execute-wrapper ratio case through submit and Python
function run_case_ratio() {
    local cas_nam="${1:-}"
    local log_lcl="${2:-}"
    local pfx_lcl="${3:-exec}"
    local method="${4:-unadj}"

    shift 4

    # shellcheck disable=SC2154
    if \
        run_capture \
            "execute compute-signal ratio ${cas_nam}" "${log_lcl}" \
            "${TEST_BASH}" "${ROOT_REPO}/scripts/execute_compute_signal.sh" \
                --threads 1 \
                --mode ratio \
                --method "${method}" \
                --csv_fil_A "${fil_A}" \
                --csv_fil_B "${fil_B}" \
                --dir_out "${dir_out}" \
                --typ_out bdg \
                --prefix "${pfx_lcl}" \
                --eps 0 \
                --dp 3 \
                --err_out "${dir_err}" \
                --nam_job "test_execute_compute_ratio_${cas_nam}" \
                --max_job 1 \
                "$@"
    then
        rec_pass "execute_compute_signal.sh ratio ${cas_nam} exits 0"
    else
        rec_fail \
            "execute_compute_signal.sh ratio ${cas_nam} failed; see" \
            "$(rec_relpath "${log_lcl}")"
    fi
}


#  Baseline unadjusted ratio with three-decimal rounding
outfile="${dir_out}/exec_ratio_A.bdg"
log="${dir_log}/execute_compute_signal_ratio_unadj.log"

run_case_ratio "unadj" "${log}" "exec" "unadj"

assert_file_nonempty "${outfile}" "execute ratio output"

if [[ -s "${outfile}" ]]; then
    assert_grep_pattern "${outfile}" $'^I\t0\t10\t2$' \
        "execute ratio output has I:0-10 = 2"
    assert_grep_pattern "${outfile}" $'^I\t10\t20\t0$' \
        "execute ratio output has I:10-20 = 0"
    assert_grep_pattern "${outfile}" $'^I\t40\t50\t4$' \
        "execute ratio output has I:40-50 = 4"
    assert_grep_pattern "${outfile}" $'^I\t60\t70\t0.333$' \
        "execute ratio output has I:60-70 = 0.333"
    assert_grep_pattern "${outfile}" $'^I\t70\t80\t1$' \
        "execute ratio output has I:70-80 = 1"
fi


#  Scaling factors are applied before ratio calculation: (2 * A) / (1 * B)
outfile="${dir_out}/exec_scl_fct_ratio_A.bdg"
log="${dir_log}/execute_compute_signal_ratio_scl_fct.log"

run_case_ratio \
    "scl_fct" "${log}" "exec_scl_fct" "unadj" \
    --csv_scl_fct 2:1

assert_file_nonempty "${outfile}" "execute scaling-factor ratio output"

if [[ -s "${outfile}" ]]; then
    assert_grep_pattern "${outfile}" $'^I\t0\t10\t4$' \
        "execute scaling-factor ratio output has I:0-10 = 4"
    assert_grep_pattern "${outfile}" $'^I\t40\t50\t8$' \
        "execute scaling-factor ratio output has I:40-50 = 8"
    assert_grep_pattern "${outfile}" $'^I\t60\t70\t0.667$' \
        "execute scaling-factor ratio output has I:60-70 = 0.667"
fi


#  Log2 ratio: log2(4 / 2) = 1 and log2(2 / 0.5) = 2
outfile="${dir_out}/exec_log2_ratio_A.bdg"
log="${dir_log}/execute_compute_signal_ratio_log2.log"

run_case_ratio "log2" "${log}" "exec_log2" "log2"

assert_file_nonempty "${outfile}" "execute log2 ratio output"

if [[ -s "${outfile}" ]]; then
    assert_grep_pattern "${outfile}" $'^I\t0\t10\t1$' \
        "execute log2 ratio output has I:0-10 = 1"
    assert_grep_pattern "${outfile}" $'^I\t40\t50\t2$' \
        "execute log2 ratio output has I:40-50 = 2"
    assert_grep_pattern "${outfile}" $'^I\t60\t70\t-1.585$' \
        "execute log2 ratio output has I:60-70 = -1.585"
fi


#  Reciprocal ratio: B / A gives 0.5, 0.25, and 3 for selected rows
outfile="${dir_out}/exec_unadj_r_ratio_A.bdg"
log="${dir_log}/execute_compute_signal_ratio_unadj_r.log"

run_case_ratio "unadj_r" "${log}" "exec_unadj_r" "unadj_r"

assert_file_nonempty "${outfile}" "execute reciprocal ratio output"

if [[ -s "${outfile}" ]]; then
    assert_grep_pattern "${outfile}" $'^I\t0\t10\t0.5$' \
        "execute reciprocal ratio output has I:0-10 = 0.5"
    assert_grep_pattern "${outfile}" $'^I\t40\t50\t0.25$' \
        "execute reciprocal ratio output has I:40-50 = 0.25"
    assert_grep_pattern "${outfile}" $'^I\t60\t70\t3$' \
        "execute reciprocal ratio output has I:60-70 = 3"
fi


#  Reciprocal log2 ratio: log2(B / A)
outfile="${dir_out}/exec_log2_r_ratio_A.bdg"
log="${dir_log}/execute_compute_signal_ratio_log2_r.log"

run_case_ratio "log2_r" "${log}" "exec_log2_r" "log2_r"

assert_file_nonempty "${outfile}" "execute reciprocal log2 ratio output"

if [[ -s "${outfile}" ]]; then
    assert_grep_pattern "${outfile}" $'^I\t0\t10\t-1$' \
        "execute reciprocal log2 ratio output has I:0-10 = -1"
    assert_grep_pattern "${outfile}" $'^I\t40\t50\t-2$' \
        "execute reciprocal log2 ratio output has I:40-50 = -2"
    assert_grep_pattern "${outfile}" $'^I\t60\t70\t1.585$' \
        "execute reciprocal log2 ratio output has I:60-70 = 1.585"
fi


#  Denominator floor: B=0.04 is floored to 0.1, so 1 / 0.1 = 10
outfile="${dir_out}/exec_dep_min_ratio_A.bdg"
log="${dir_log}/execute_compute_signal_ratio_dep_min.log"

run_case_ratio "dep_min" "${log}" "exec_dep_min" "unadj" --csv_dep_min 0.1

assert_file_nonempty "${outfile}" "execute dep_min ratio output"

if [[ -s "${outfile}" ]]; then
    assert_grep_pattern "${outfile}" $'^I\t50\t60\t10$' \
        "execute dep_min ratio output has I:50-60 = 10"
fi


#  Epsilon guards denominator values at or below eps: B=0.04 <= 0.05 -> nan
outfile="${dir_out}/exec_eps_ratio_A.bdg"
log="${dir_log}/execute_compute_signal_ratio_eps.log"

run_case_ratio "eps" "${log}" "exec_eps" "unadj" --eps 0.05

assert_file_nonempty "${outfile}" "execute eps ratio output"

if [[ -s "${outfile}" ]]; then
    assert_grep_pattern "${outfile}" $'^I\t0\t10\t2$' \
        "execute eps ratio output retains I:0-10 = 2"
    assert_grep_pattern "${outfile}" $'^I\t50\t60\tnan$' \
        "execute eps ratio output has I:50-60 = nan"
    assert_grep_pattern "${outfile}" $'^I\t60\t70\t0.333$' \
        "execute eps ratio output retains I:60-70 = 0.333"
fi


#  Pseudocounts: (0 + 1) / (2 + 1) = 0.333 at three decimals
outfile="${dir_out}/exec_pseudo_ratio_A.bdg"
log="${dir_log}/execute_compute_signal_ratio_pseudo.log"

run_case_ratio "pseudo" "${log}" "exec_pseudo" "unadj" --csv_pseudo 1:1

assert_file_nonempty "${outfile}" "execute pseudo ratio output"

if [[ -s "${outfile}" ]]; then
    assert_grep_pattern "${outfile}" $'^I\t10\t20\t0.333$' \
        "execute pseudo ratio output has I:10-20 = 0.333"
fi


#  Drop non-finite rows while preserving finite ratio rows
outfile="${dir_out}/exec_drp_nan_ratio_A.bdg"
log="${dir_log}/execute_compute_signal_ratio_drp_nan.log"

run_case_ratio "drp_nan" "${log}" "exec_drp_nan" "unadj" --drp_nan

assert_file_nonempty "${outfile}" "execute drp_nan ratio output"

if [[ -s "${outfile}" ]]; then
    assert_grep_pattern "${outfile}" $'^I\t0\t10\t2$' \
        "execute drp_nan ratio output retains I:0-10 = 2"
    assert_grep_pattern "${outfile}" $'^I\t60\t70\t0.333$' \
        "execute drp_nan ratio output retains I:60-70 = 0.333"
    assert_no_grep_pattern "${outfile}" $'^I\t20\t30\t' \
        "execute drp_nan ratio output omits I:20-30"
    assert_no_grep_pattern "${outfile}" $'^I\t30\t40\t' \
        "execute drp_nan ratio output omits I:30-40"
fi


#  Zero-zero skipping before scaling removes the A=0, B=0 bin
outfile="${dir_out}/exec_skip_00_ratio_A.bdg"
log="${dir_log}/execute_compute_signal_ratio_skip_00.log"

run_case_ratio "skip_00" "${log}" "exec_skip_00" "unadj" --skip_00 pre_scale

assert_file_nonempty "${outfile}" "execute skip_00 ratio output"

if [[ -s "${outfile}" ]]; then
    assert_no_grep_pattern "${outfile}" $'^I\t30\t40\t' \
        "execute skip_00 ratio output omits I:30-40"
fi


#  Post-scale zero-zero skipping can remove bins that are non-zero before
#+ scaling: I:50-60 has A=1 and B=0.04, then scales to 0.001 and 0.004
outfile="${dir_out}/exec_skip_00_post_scale_ratio_A.bdg"
log="${dir_log}/execute_compute_signal_ratio_skip_00_post_scale.log"

run_case_ratio \
    "skip_00_post_scale" "${log}" "exec_skip_00_post_scale" "unadj" \
    --csv_scl_fct 0.001:0.1 \
    --eps 0.005 \
    --skip_00 post_scale

assert_file_nonempty "${outfile}" "execute skip_00 post_scale ratio output"

if [[ -s "${outfile}" ]]; then
    assert_grep_pattern "${outfile}" $'^I\t0\t10\t0.02$' \
        "execute skip_00 post_scale ratio output has I:0-10 = 0.02"
    assert_grep_pattern "${outfile}" $'^I\t40\t50\t0.04$' \
        "execute skip_00 post_scale ratio output has I:40-50 = 0.04"
    assert_no_grep_pattern "${outfile}" $'^I\t50\t60\t' \
        "execute skip_00 post_scale ratio output omits scaled zero-zero I:50-60"
fi


#  Track sidecar should be generated and should omit non-finite rows
outfile="${dir_out}/exec_track_ratio_A.bdg"
trackfile="${dir_out}/exec_track_ratio_A.track.bdg"
log="${dir_log}/execute_compute_signal_ratio_track.log"

run_case_ratio "track" "${log}" "exec_track" "unadj" --track

assert_file_nonempty "${outfile}" "execute track main ratio output"
assert_file_nonempty "${trackfile}" "execute track sidecar output"

if [[ -s "${trackfile}" ]]; then
    assert_grep_pattern "${trackfile}" $'^I\t0\t10\t2$' \
        "execute track sidecar retains I:0-10 = 2"
    assert_grep_pattern "${trackfile}" $'^I\t60\t70\t0.333$' \
        "execute track sidecar retains I:60-70 = 0.333"
    assert_no_grep_pattern "${trackfile}" $'^I\t20\t30\t' \
        "execute track sidecar omits I:20-30"
    assert_no_grep_pattern "${trackfile}" $'^I\t30\t40\t' \
        "execute track sidecar omits I:30-40"
fi


#  Gzipped bedGraph input and output should round-trip through execute mode
outfile="${dir_out}/exec_gzip_io_ratio_A.bdg.gz"
outfile_txt="${dir_out}/exec_gzip_io_ratio_A.bdg"
log="${dir_log}/execute_compute_signal_ratio_gzip_io.log"

if \
    run_capture \
        "execute compute-signal ratio gzip_io" "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/execute_compute_signal.sh" \
            --threads 1 \
            --mode ratio \
            --method unadj \
            --csv_fil_A "${fil_gz_A}" \
            --csv_fil_B "${fil_gz_B}" \
            --dir_out "${dir_out}" \
            --typ_out bdg.gz \
            --prefix "exec_gzip_io" \
            --eps 0 \
            --dp 3 \
            --err_out "${dir_err}" \
            --nam_job "test_execute_compute_ratio_gzip_io" \
            --max_job 1 \
            --csv_scl_fct NA \
            --csv_dep_min NA \
            --csv_pseudo NA
then
    rec_pass "execute_compute_signal.sh ratio gzip_io exits 0"
else
    rec_fail \
        "execute_compute_signal.sh ratio gzip_io failed; see" \
        "$(rec_relpath "${log}")"
fi

assert_file_nonempty "${outfile}" "execute gzip ratio output"

if [[ -s "${outfile}" ]]; then
    if gzip -t "${outfile}"; then
        rec_pass "execute gzip ratio output passes gzip integrity check"
        gzip -cd "${outfile}" > "${outfile_txt}"
    else
        rec_fail "execute gzip ratio output fails gzip integrity check"
    fi
fi

if [[ -s "${outfile_txt}" ]]; then
    assert_grep_pattern "${outfile_txt}" $'^I\t0\t10\t2$' \
        "execute gzip ratio output has I:0-10 = 2"
    assert_grep_pattern "${outfile_txt}" $'^I\t40\t50\t4$' \
        "execute gzip ratio output has I:40-50 = 4"
    assert_grep_pattern "${outfile_txt}" $'^I\t60\t70\t0.333$' \
        "execute gzip ratio output has I:60-70 = 0.333"
fi


#  Header/prefix skipping should propagate through execute to submit/Python
outfile="${dir_out}/exec_skp_pfx_ratio_headers_A.bdg"
log="${dir_log}/execute_compute_signal_ratio_skp_pfx.log"

if \
    run_capture \
        "execute compute-signal ratio skp_pfx" "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/execute_compute_signal.sh" \
            --threads 1 \
            --mode ratio \
            --method unadj \
            --csv_fil_A "${fil_hdr_A}" \
            --csv_fil_B "${fil_hdr_B}" \
            --dir_out "${dir_out}" \
            --typ_out bdg \
            --prefix "exec_skp_pfx" \
            --eps 0 \
            --dp 3 \
            --err_out "${dir_err}" \
            --nam_job "test_execute_compute_ratio_skp_pfx" \
            --max_job 1 \
            --csv_scl_fct NA \
            --csv_dep_min NA \
            --csv_pseudo NA \
            --skp_pfx "#,track,browser,customHeader"
then
    rec_pass "execute_compute_signal.sh ratio skp_pfx exits 0"
else
    rec_fail \
        "execute_compute_signal.sh ratio skp_pfx failed; see" \
        "$(rec_relpath "${log}")"
fi

assert_file_nonempty "${outfile}" "execute skp_pfx ratio output"

if [[ -s "${outfile}" ]]; then
    assert_grep_pattern "${outfile}" $'^I\t0\t10\t2$' \
        "execute skp_pfx ratio output has I:0-10 = 2"
    assert_grep_pattern "${outfile}" $'^I\t40\t50\t4$' \
        "execute skp_pfx ratio output has I:40-50 = 4"
    assert_grep_pattern "${outfile}" $'^I\t60\t70\t0.333$' \
        "execute skp_pfx ratio output has I:60-70 = 0.333"
    assert_no_grep_pattern "${outfile}" '^track' \
        "execute skp_pfx ratio output omits track headers"
    assert_no_grep_pattern "${outfile}" '^browser' \
        "execute skp_pfx ratio output omits browser headers"
    assert_no_grep_pattern "${outfile}" '^customHeader' \
        "execute skp_pfx ratio output omits custom headers"
    assert_no_grep_pattern "${outfile}" '^#' \
        "execute skp_pfx ratio output omits comment headers"
fi


finish
