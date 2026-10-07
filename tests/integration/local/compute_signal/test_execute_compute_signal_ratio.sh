#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_compute_signal_ratio.sh
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

TEST_NAME="execute compute-signal ratio"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Define fixture and output paths for the execute-to-submit ratio path.
dir_fx="${ROOT_REPO}/tests/fixtures/compute_signal/bedgraph"
chr_siz="${ROOT_REPO}/tests/fixtures/compute_signal/reference/tiny.fa.fai"
fil_A="${dir_fx}/ratio_A.bdg"
fil_B="${dir_fx}/ratio_B.bdg"
fil_A_hdr="${dir_fx}/ratio_headers_A.bdg"
fil_B_hdr="${dir_fx}/ratio_headers_B.bdg"
fil_A_gz="${dir_fx}/ratio_A.bdg.gz"
fil_B_gz="${dir_fx}/ratio_B.bdg.gz"

tmp="${TEST_DIR_TMP}/execute_compute_signal_ratio"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/compute_signal"

fil_out_linear="${dir_out}/exec_ratio_A.bdg"
fil_out_scl_fct="${dir_out}/exec_scl_fct_ratio_A.bdg"
fil_out_log2="${dir_out}/exec_log2_ratio_A.bdg"
fil_out_linear_r="${dir_out}/exec_linear_r_ratio_A.bdg"
fil_out_log2_r="${dir_out}/exec_log2_r_ratio_A.bdg"
fil_out_dep_min="${dir_out}/exec_dep_min_ratio_A.bdg"
fil_out_eps="${dir_out}/exec_eps_ratio_A.bdg"
fil_out_pseudo="${dir_out}/exec_pseudo_ratio_A.bdg"
fil_out_drp_nan="${dir_out}/exec_drp_nan_ratio_A.bdg"
fil_out_skip_00="${dir_out}/exec_skip_00_ratio_A.bdg"
fil_out_skip_00_post_scale="${dir_out}/exec_skip_00_post_scale_ratio_A.bdg"
fil_out_track="${dir_out}/exec_track_ratio_A.bdg"
fil_out_gzip_io="${dir_out}/exec_gzip_io_ratio_A.bdg.gz"
fil_out_skp_pfx="${dir_out}/exec_skp_pfx_ratio_headers_A.bdg"

trackfile_track="${dir_out}/exec_track_ratio_A.track.bdg"
outfile_txt_gzip_io="${dir_out}/exec_gzip_io_ratio_A.bdg"

log_linear="${dir_log}/execute_compute_signal_ratio_linear.log"
log_scl_fct="${dir_log}/execute_compute_signal_ratio_scl_fct.log"
log_log2="${dir_log}/execute_compute_signal_ratio_log2.log"
log_linear_r="${dir_log}/execute_compute_signal_ratio_linear_r.log"
log_log2_r="${dir_log}/execute_compute_signal_ratio_log2_r.log"
log_dep_min="${dir_log}/execute_compute_signal_ratio_dep_min.log"
log_eps="${dir_log}/execute_compute_signal_ratio_eps.log"
log_pseudo="${dir_log}/execute_compute_signal_ratio_pseudo.log"
log_drp_nan="${dir_log}/execute_compute_signal_ratio_drp_nan.log"
log_skip_00="${dir_log}/execute_compute_signal_ratio_skip_00.log"
log_skip_00_post_scale="${dir_log}/execute_compute_signal_ratio_skip_00_post_scale.log"
log_track="${dir_log}/execute_compute_signal_ratio_track.log"
log_gzip_io="${dir_log}/execute_compute_signal_ratio_gzip_io.log"
log_skp_pfx="${dir_log}/execute_compute_signal_ratio_skp_pfx.log"


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_err}" "${dir_log}"

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

require_files_nonempty \
    "${chr_siz}" \
    "${fil_A_gz}" \
    "${fil_B_gz}" || {
    finish
    exit $?
}


# Baseline linear ratio with three-decimal rounding.
run_case_compute_signal_ratio \
    execute \
    "linear" \
    "exec" \
    "${log_linear}" \
    "linear" \
    "${fil_A}" \
    "${fil_B}" \
    "${dir_out}" \
    "${dir_err}" \
    --chr_siz "${chr_siz}" \
    --strict_bins

assert_file_nonempty \
    "${fil_out_linear}" \
    "execute ratio output"

if [[ -s "${fil_out_linear}" ]]; then
    assert_pattern_found \
        "${fil_out_linear}" \
        $'^I\t0\t10\t2$' \
        "execute ratio output has I:0-10 = 2"

    assert_pattern_found \
        "${fil_out_linear}" \
        $'^I\t10\t20\t0$' \
        "execute ratio output has I:10-20 = 0"

    assert_pattern_found \
        "${fil_out_linear}" \
        $'^I\t40\t50\t4$' \
        "execute ratio output has I:40-50 = 4"

    assert_pattern_found \
        "${fil_out_linear}" \
        $'^I\t60\t70\t0.333$' \
        "execute ratio output has I:60-70 = 0.333"

    assert_pattern_found \
        "${fil_out_linear}" \
        $'^I\t70\t80\t1$' \
        "execute ratio output has I:70-80 = 1"
fi


# Scaling factors are applied before ratio calculation: (2 * A) / (1 * B).
run_case_compute_signal_ratio \
    execute \
    "scl_fct" \
    "exec_scl_fct" \
    "${log_scl_fct}" \
    "linear" \
    "${fil_A}" \
    "${fil_B}" \
    "${dir_out}" \
    "${dir_err}" \
    --csv_scl_fct \
    2:1

assert_file_nonempty \
    "${fil_out_scl_fct}" \
    "execute scaling-factor ratio output"

if [[ -s "${fil_out_scl_fct}" ]]; then
    assert_pattern_found \
        "${fil_out_scl_fct}" \
        $'^I\t0\t10\t4$' \
        "execute scaling-factor ratio output has I:0-10 = 4"

    assert_pattern_found \
        "${fil_out_scl_fct}" \
        $'^I\t40\t50\t8$' \
        "execute scaling-factor ratio output has I:40-50 = 8"

    assert_pattern_found \
        "${fil_out_scl_fct}" \
        $'^I\t60\t70\t0.667$' \
        "execute scaling-factor ratio output has I:60-70 = 0.667"
fi


# Log2 ratio: log2(4 / 2) = 1 and log2(2 / 0.5) = 2.
run_case_compute_signal_ratio \
    execute \
    "log2" \
    "exec_log2" \
    "${log_log2}" \
    "log2" \
    "${fil_A}" \
    "${fil_B}" \
    "${dir_out}" \
    "${dir_err}"

assert_file_nonempty \
    "${fil_out_log2}" \
    "execute log2 ratio output"

if [[ -s "${fil_out_log2}" ]]; then
    assert_pattern_found \
        "${fil_out_log2}" \
        $'^I\t0\t10\t1$' \
        "execute log2 ratio output has I:0-10 = 1"

    assert_pattern_found \
        "${fil_out_log2}" \
        $'^I\t40\t50\t2$' \
        "execute log2 ratio output has I:40-50 = 2"

    assert_pattern_found \
        "${fil_out_log2}" \
        $'^I\t60\t70\t-1.585$' \
        "execute log2 ratio output has I:60-70 = -1.585"
fi


# Reciprocal ratio: B / A gives 0.5, 0.25, and 3 for selected rows.
run_case_compute_signal_ratio \
    execute \
    "linear_r" \
    "exec_linear_r" \
    "${log_linear_r}" \
    "linear_r" \
    "${fil_A}" \
    "${fil_B}" \
    "${dir_out}" \
    "${dir_err}"

assert_file_nonempty \
    "${fil_out_linear_r}" \
    "execute reciprocal ratio output"

if [[ -s "${fil_out_linear_r}" ]]; then
    assert_pattern_found \
        "${fil_out_linear_r}" \
        $'^I\t0\t10\t0.5$' \
        "execute reciprocal ratio output has I:0-10 = 0.5"

    assert_pattern_found \
        "${fil_out_linear_r}" \
        $'^I\t40\t50\t0.25$' \
        "execute reciprocal ratio output has I:40-50 = 0.25"

    assert_pattern_found \
        "${fil_out_linear_r}" \
        $'^I\t60\t70\t3$' \
        "execute reciprocal ratio output has I:60-70 = 3"
fi


# Reciprocal log2 ratio: log2(B / A).
run_case_compute_signal_ratio \
    execute \
    "log2_r" \
    "exec_log2_r" \
    "${log_log2_r}" \
    "log2_r" \
    "${fil_A}" \
    "${fil_B}" \
    "${dir_out}" \
    "${dir_err}"

assert_file_nonempty \
    "${fil_out_log2_r}" \
    "execute reciprocal log2 ratio output"

if [[ -s "${fil_out_log2_r}" ]]; then
    assert_pattern_found \
        "${fil_out_log2_r}" \
        $'^I\t0\t10\t-1$' \
        "execute reciprocal log2 ratio output has I:0-10 = -1"

    assert_pattern_found \
        "${fil_out_log2_r}" \
        $'^I\t40\t50\t-2$' \
        "execute reciprocal log2 ratio output has I:40-50 = -2"

    assert_pattern_found \
        "${fil_out_log2_r}" \
        $'^I\t60\t70\t1.585$' \
        "execute reciprocal log2 ratio output has I:60-70 = 1.585"
fi


# Denominator floor: B=0.04 is floored to 0.1, so 1 / 0.1 = 10.
run_case_compute_signal_ratio \
    execute \
    "dep_min" \
    "exec_dep_min" \
    "${log_dep_min}" \
    "linear" \
    "${fil_A}" \
    "${fil_B}" \
    "${dir_out}" \
    "${dir_err}" \
    --csv_dep_min \
    0.1

assert_file_nonempty \
    "${fil_out_dep_min}" \
    "execute dep_min ratio output"

if [[ -s "${fil_out_dep_min}" ]]; then
    assert_pattern_found \
        "${fil_out_dep_min}" \
        $'^I\t50\t60\t10$' \
        "execute dep_min ratio output has I:50-60 = 10"
fi


# Epsilon guards denominator values at or below eps: B=0.04 <= 0.05 -> nan.
run_case_compute_signal_ratio \
    execute \
    "eps" \
    "exec_eps" \
    "${log_eps}" \
    "linear" \
    "${fil_A}" \
    "${fil_B}" \
    "${dir_out}" \
    "${dir_err}" \
    --eps \
    0.05

assert_file_nonempty \
    "${fil_out_eps}" \
    "execute eps ratio output"

if [[ -s "${fil_out_eps}" ]]; then
    assert_pattern_found \
        "${fil_out_eps}" \
        $'^I\t0\t10\t2$' \
        "execute eps ratio output retains I:0-10 = 2"

    assert_pattern_found \
        "${fil_out_eps}" \
        $'^I\t50\t60\tnan$' \
        "execute eps ratio output has I:50-60 = nan"

    assert_pattern_found \
        "${fil_out_eps}" \
        $'^I\t60\t70\t0.333$' \
        "execute eps ratio output retains I:60-70 = 0.333"
fi


# Pseudocounts: (0 + 1) / (2 + 1) = 0.333 at three decimals.
run_case_compute_signal_ratio \
    execute \
    "pseudo" \
    "exec_pseudo" \
    "${log_pseudo}" \
    "linear" \
    "${fil_A}" \
    "${fil_B}" \
    "${dir_out}" \
    "${dir_err}" \
    --csv_pseudo \
    1:1

assert_file_nonempty \
    "${fil_out_pseudo}" \
    "execute pseudo ratio output"

if [[ -s "${fil_out_pseudo}" ]]; then
    assert_pattern_found \
        "${fil_out_pseudo}" \
        $'^I\t10\t20\t0.333$' \
        "execute pseudo ratio output has I:10-20 = 0.333"
fi


# Drop non-finite rows while preserving finite ratio rows.
run_case_compute_signal_ratio \
    execute \
    "drp_nan" \
    "exec_drp_nan" \
    "${log_drp_nan}" \
    "linear" \
    "${fil_A}" \
    "${fil_B}" \
    "${dir_out}" \
    "${dir_err}" \
    --drp_nan

assert_file_nonempty \
    "${fil_out_drp_nan}" \
    "execute drp_nan ratio output"

if [[ -s "${fil_out_drp_nan}" ]]; then
    assert_pattern_found \
        "${fil_out_drp_nan}" \
        $'^I\t0\t10\t2$' \
        "execute drp_nan ratio output retains I:0-10 = 2"

    assert_pattern_found \
        "${fil_out_drp_nan}" \
        $'^I\t60\t70\t0.333$' \
        "execute drp_nan ratio output retains I:60-70 = 0.333"

    assert_pattern_absent \
        "${fil_out_drp_nan}" \
        $'^I\t20\t30\t' \
        "execute drp_nan ratio output omits I:20-30"

    assert_pattern_absent \
        "${fil_out_drp_nan}" \
        $'^I\t30\t40\t' \
        "execute drp_nan ratio output omits I:30-40"
fi


# Zero-zero skipping before scaling removes the A=0, B=0 bin.
run_case_compute_signal_ratio \
    execute \
    "skip_00" \
    "exec_skip_00" \
    "${log_skip_00}" \
    "linear" \
    "${fil_A}" \
    "${fil_B}" \
    "${dir_out}" \
    "${dir_err}" \
    --skip_00 \
    pre_scale

assert_file_nonempty \
    "${fil_out_skip_00}" \
    "execute skip_00 ratio output"

if [[ -s "${fil_out_skip_00}" ]]; then
    assert_pattern_absent \
        "${fil_out_skip_00}" \
        $'^I\t30\t40\t' \
        "execute skip_00 ratio output omits I:30-40"
fi


# Post-scale zero-zero skipping can remove bins that are non-zero before
# scaling: I:50-60 has A=1 and B=0.04, then scales to 0.001 and 0.004.
run_case_compute_signal_ratio \
    execute \
    "skip_00_post_scale" \
    "exec_skip_00_post_scale" \
    "${log_skip_00_post_scale}" \
    "linear" \
    "${fil_A}" \
    "${fil_B}" \
    "${dir_out}" \
    "${dir_err}" \
    --csv_scl_fct \
    0.001:0.1 \
    --eps \
    0.005 \
    --skip_00 \
    post_scale

assert_file_nonempty \
    "${fil_out_skip_00_post_scale}" \
    "execute skip_00 post_scale ratio output"

if [[ -s "${fil_out_skip_00_post_scale}" ]]; then
    assert_pattern_found \
        "${fil_out_skip_00_post_scale}" \
        $'^I\t0\t10\t0.02$' \
        "execute skip_00 post_scale ratio output has I:0-10 = 0.02"

    assert_pattern_found \
        "${fil_out_skip_00_post_scale}" \
        $'^I\t40\t50\t0.04$' \
        "execute skip_00 post_scale ratio output has I:40-50 = 0.04"

    assert_pattern_absent \
        "${fil_out_skip_00_post_scale}" \
        $'^I\t50\t60\t' \
        "execute skip_00 post_scale ratio output omits scaled zero-zero I:50-60"
fi


# Track sidecar should be generated and should omit non-finite rows.
run_case_compute_signal_ratio \
    execute \
    "track" \
    "exec_track" \
    "${log_track}" \
    "linear" \
    "${fil_A}" \
    "${fil_B}" \
    "${dir_out}" \
    "${dir_err}" \
    --track

assert_file_nonempty \
    "${fil_out_track}" \
    "execute track main ratio output"

assert_file_nonempty \
    "${trackfile_track}" \
    "execute track sidecar output"

if [[ -s "${trackfile_track}" ]]; then
    assert_pattern_found \
        "${trackfile_track}" \
        $'^I\t0\t10\t2$' \
        "execute track sidecar retains I:0-10 = 2"

    assert_pattern_found \
        "${trackfile_track}" \
        $'^I\t60\t70\t0.333$' \
        "execute track sidecar retains I:60-70 = 0.333"

    assert_pattern_absent \
        "${trackfile_track}" \
        $'^I\t20\t30\t' \
        "execute track sidecar omits I:20-30"

    assert_pattern_absent \
        "${trackfile_track}" \
        $'^I\t30\t40\t' \
        "execute track sidecar omits I:30-40"
fi


# shellcheck disable=SC2154
# Gzipped bedGraph input and output should round-trip through execute mode.
if \
    run_capture \
        "execute compute-signal ratio gzip_io" \
        "${log_gzip_io}" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/execute_compute_signal.sh" \
            --env_nam "${env_nam}" \
            --threads 1 \
            --mode ratio \
            --method linear \
            --csv_fil_A "${fil_A_gz}" \
            --csv_fil_B "${fil_B_gz}" \
            --dir_out "${dir_out}" \
            --typ_out bdg.gz \
            --prefix "exec_gzip_io" \
            --eps 0 \
            --dp 3 \
            --dir_eo "${dir_err}" \
            --nam_job "test_execute_compute_ratio_gzip_io" \
            --max_job 1 \
            --csv_scl_fct NA \
            --csv_dep_min NA \
            --csv_pseudo NA
then
    record_pass "execute_compute_signal.sh ratio gzip_io exits 0"
else
    record_fail \
        "execute_compute_signal.sh ratio gzip_io failed; see" \
        "$(print_relpath "${log_gzip_io}")"
fi

assert_file_nonempty \
    "${fil_out_gzip_io}" \
    "execute gzip ratio output"

if [[ -s "${fil_out_gzip_io}" ]]; then
    if \
        gzip -t "${fil_out_gzip_io}"
    then
        record_pass "execute gzip ratio output passes gzip integrity check"
        gzip -cd "${fil_out_gzip_io}" > "${outfile_txt_gzip_io}"
    else
        record_fail "execute gzip ratio output fails gzip integrity check"
    fi
fi

assert_file_nonempty \
    "${outfile_txt_gzip_io}" \
    "execute decompressed gzip ratio output"

if [[ -s "${outfile_txt_gzip_io}" ]]; then
    assert_pattern_found \
        "${outfile_txt_gzip_io}" \
        $'^I\t0\t10\t2$' \
        "execute gzip ratio output has I:0-10 = 2"

    assert_pattern_found \
        "${outfile_txt_gzip_io}" \
        $'^I\t40\t50\t4$' \
        "execute gzip ratio output has I:40-50 = 4"

    assert_pattern_found \
        "${outfile_txt_gzip_io}" \
        $'^I\t60\t70\t0.333$' \
        "execute gzip ratio output has I:60-70 = 0.333"
fi


# Header/prefix skipping should propagate through execute to submit/Python.
if \
    run_capture \
        "execute compute-signal ratio skp_pfx" \
        "${log_skp_pfx}" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/execute_compute_signal.sh" \
            --env_nam "${env_nam}" \
            --threads 1 \
            --mode ratio \
            --method linear \
            --csv_fil_A "${fil_A_hdr}" \
            --csv_fil_B "${fil_B_hdr}" \
            --dir_out "${dir_out}" \
            --typ_out bdg \
            --prefix "exec_skp_pfx" \
            --eps 0 \
            --dp 3 \
            --dir_eo "${dir_err}" \
            --nam_job "test_execute_compute_ratio_skp_pfx" \
            --max_job 1 \
            --csv_scl_fct NA \
            --csv_dep_min NA \
            --csv_pseudo NA \
            --skp_pfx "#,track,browser,customHeader"
then
    record_pass "execute_compute_signal.sh ratio skp_pfx exits 0"
else
    record_fail \
        "execute_compute_signal.sh ratio skp_pfx failed; see" \
        "$(print_relpath "${log_skp_pfx}")"
fi

assert_file_nonempty \
    "${fil_out_skp_pfx}" \
    "execute skp_pfx ratio output"

if [[ -s "${fil_out_skp_pfx}" ]]; then
    assert_pattern_found \
        "${fil_out_skp_pfx}" \
        $'^I\t0\t10\t2$' \
        "execute skp_pfx ratio output has I:0-10 = 2"

    assert_pattern_found \
        "${fil_out_skp_pfx}" \
        $'^I\t40\t50\t4$' \
        "execute skp_pfx ratio output has I:40-50 = 4"

    assert_pattern_found \
        "${fil_out_skp_pfx}" \
        $'^I\t60\t70\t0.333$' \
        "execute skp_pfx ratio output has I:60-70 = 0.333"

    assert_pattern_absent \
        "${fil_out_skp_pfx}" \
        '^track' \
        "execute skp_pfx ratio output omits track headers"

    assert_pattern_absent \
        "${fil_out_skp_pfx}" \
        '^browser' \
        "execute skp_pfx ratio output omits browser headers"

    assert_pattern_absent \
        "${fil_out_skp_pfx}" \
        '^customHeader' \
        "execute skp_pfx ratio output omits custom headers"

    assert_pattern_absent \
        "${fil_out_skp_pfx}" \
        '^#' \
        "execute skp_pfx ratio output omits comment headers"
fi


# Mode separation: ratio runs take '--chr_siz', which they use to validate
# bedGraph bounds, and must never receive the signal-only window options.
log_rat_mode="${tmp}/logs/test_execute_compute_ratio_dep_min.exec_dep_min_ratio_A.stderr.txt"

assert_file_nonempty \
    "${log_rat_mode}" \
    "execute ratio stderr log for mode separation"

if [[ -s "${log_rat_mode}" ]]; then
    assert_pattern_found \
        "${log_rat_mode}" \
        "--chr_siz" \
        "execute ratio still forwards '--chr_siz'"

    assert_pattern_absent \
        "${log_rat_mode}" \
        "--siz_win" \
        "execute ratio omits '--siz_win'"

    assert_pattern_absent \
        "${log_rat_mode}" \
        "--engine" \
        "execute ratio omits '--engine'"
fi


# Match the method diagnostic, not exit status: an unrelated argument fault
# also exits non-zero, so a status check would pass for a valid method.
arr_mth_gone=(
    unadj unadj_r unadjusted r raw u s smp simple
    2 lg2 rr ur sr 2r l2r lg2_r bogus
)

for retired in "${arr_mth_gone[@]}"; do
    out_rej="$(
        bash "${ROOT_REPO}/bin/execute_compute_signal.sh" \
            --mode ratio \
            --method "${retired}" \
            --csv_fil_A "${fil_A}" \
            --csv_fil_B "${fil_B}" \
            --dir_out "${dir_out}" 2>&1 || true
    )"

    if [[ "${out_rej}" == *"invalid value for '--method'"* ]]; then
        record_pass "execute rejects retired ratio method '${retired}'"
    else
        record_fail "execute did not reject ratio method '${retired}'"
    fi
done


# Read the accepted vocabulary out of the tool itself, so a method added to
# 'METHOD_CANON' without a matching wrapper arm fails here.
mapfile -t arr_mth_keep < <(
    "${TEST_MANAGED_PYTHON}" -c \
        "from protocol_chipseq_signal_norm.cli.compute_signal_ratio import \
METHOD_CANON; print(chr(10).join(METHOD_CANON))"
)

if [[ "${#arr_mth_keep[@]}" -eq 0 ]]; then
    record_fail "could not read the ratio method vocabulary"
fi


# The same probe stays silent for every kept spelling, which is what makes the
# loop above discriminating rather than inert.
for kept in "${arr_mth_keep[@]}"; do
    out_keep="$(
        bash "${ROOT_REPO}/bin/execute_compute_signal.sh" \
            --mode ratio \
            --method "${kept}" \
            --csv_fil_A "${fil_A}" \
            --csv_fil_B "${fil_B}" \
            --dir_out "${dir_out}" 2>&1 || true
    )"

    if [[ "${out_keep}" == *"invalid value for '--method'"* ]]; then
        record_fail "execute wrongly rejected ratio method '${kept}'"
    else
        record_pass "execute accepts ratio method '${kept}'"
    fi
done


# 'edger' reads the counts beside each track; 'count' tracks, unlike 'norm',
# make an A/B swap change the value, so these checks can see one.
# shellcheck source=lib/bash/core/format_outputs.sh
source "${ROOT_REPO}/lib/bash/core/format_outputs.sh"

dir_edg="${tmp}/edger"
mkdir -p "${dir_edg}"

# The fixture 'count' tracks, with the counts the ratio stage must find for
# itself beside each.
dir_cnt="${ROOT_REPO}/tests/fixtures/compute_signal/bedgraph/count"
trk_A="${dir_cnt}/tiny_se.bedGraph"
trk_B="${dir_cnt}/tiny_pe.bedGraph"

require_files_nonempty \
    "${trk_A}" \
    "${trk_B}" \
    "$(derive_report_path "${trk_A}" n_frg)" \
    "$(derive_report_path "${trk_B}" n_frg)" \
    "$(derive_report_path "${trk_A}" n_ovlp)" \
    "$(derive_report_path "${trk_B}" n_ovlp)" || {
    finish
    exit $?
}

# Both arms run 'compute_pseudo', so this checks the plumbing (report paths,
# A/B pairing, precision, '--typ_sig' forwarding), not the arithmetic.
exp_pseudo="$(
    "${TEST_MANAGED_PYTHON}" \
        -m protocol_chipseq_signal_norm.cli.compute_pseudo \
        --method edger \
        --typ_sig count \
        --fil_A "${trk_A}" \
        --fil_B "${trk_B}" \
        --n_frg_A "$(< "$(derive_report_path "${trk_A}" n_frg)")" \
        --n_frg_B "$(< "$(derive_report_path "${trk_B}" n_frg)")" \
        --n_ovlp_A "$(< "$(derive_report_path "${trk_A}" n_ovlp)")" \
        --n_ovlp_B "$(< "$(derive_report_path "${trk_B}" n_ovlp)")" \
        2>/dev/null
)"

"${TEST_BASH}" "${ROOT_REPO}/bin/execute_compute_signal.sh" \
    --env_nam "${env_nam}" \
    --threads 1 \
    --mode ratio \
    --method log2 \
    --typ_sig count \
    --csv_fil_A "${trk_A}" \
    --csv_fil_B "${trk_B}" \
    --dir_out "${dir_edg}" \
    --typ_out bedGraph \
    --prefix derived \
    --dir_eo "${dir_err}" \
    --nam_job "edger_ratio" \
    --max_job 1 \
    --csv_scl_fct NA \
    --csv_dep_min NA \
    --csv_pseudo edger \
    --verbose > "${dir_edg}/ratio.log" 2>&1 || true

log_edg="${dir_err}/edger_ratio.derived_tiny_se.stderr.txt"

swap_pseudo="$(
    "${TEST_MANAGED_PYTHON}" \
        -m protocol_chipseq_signal_norm.cli.compute_pseudo \
        --method edger \
        --typ_sig count \
        --fil_A "${trk_A}" \
        --fil_B "${trk_B}" \
        --n_frg_A "$(< "$(derive_report_path "${trk_B}" n_frg)")" \
        --n_frg_B "$(< "$(derive_report_path "${trk_A}" n_frg)")" \
        --n_ovlp_A "$(< "$(derive_report_path "${trk_B}" n_ovlp)")" \
        --n_ovlp_B "$(< "$(derive_report_path "${trk_A}" n_ovlp)")" \
        2>/dev/null
)"

if [[ "${exp_pseudo}" != "${swap_pseudo}" ]]; then
    record_pass "the equivalence check can see an A/B pairing error"
else
    record_fail \
        "the equivalence check has no pairing power: swapping A and B gives" \
        "the same value, so a mis-wired pair would pass"
fi

# 'norm' uses pair means, so an A/B swap is harmless, but swapping 'N' and 'L'
# inverts 'k = L / N'; pin both.
norm_ok="$(
    "${TEST_MANAGED_PYTHON}" \
        -m protocol_chipseq_signal_norm.cli.compute_pseudo \
        --method edger \
        --typ_sig norm \
        --fil_A "${trk_A}" \
        --fil_B "${trk_B}" \
        --n_frg_A 3 \
        --n_frg_B 7 \
        --n_ovlp_A 11 \
        --n_ovlp_B 29 \
        2>/dev/null
)"
norm_swap="$(
    "${TEST_MANAGED_PYTHON}" \
        -m protocol_chipseq_signal_norm.cli.compute_pseudo \
        --method edger \
        --typ_sig norm \
        --fil_A "${trk_A}" \
        --fil_B "${trk_B}" \
        --n_frg_A 7 \
        --n_frg_B 3 \
        --n_ovlp_A 29 \
        --n_ovlp_B 11 \
        2>/dev/null
)"
norm_trans="$(
    "${TEST_MANAGED_PYTHON}" \
        -m protocol_chipseq_signal_norm.cli.compute_pseudo \
        --method edger \
        --typ_sig norm \
        --fil_A "${trk_A}" \
        --fil_B "${trk_B}" \
        --n_frg_A 11 \
        --n_frg_B 29 \
        --n_ovlp_A 3 \
        --n_ovlp_B 7 \
        2>/dev/null
)"

if [[ "${norm_ok}" == "${norm_ok#*:}:${norm_ok#*:}" ]]; then
    record_pass "'norm' gives both tracks the same pseudocount, as designed"
else
    record_fail "'norm' is no longer symmetric in A and B: '${norm_ok}'"
fi

if [[ "${norm_ok}" == "${norm_swap}" ]]; then
    record_pass "'norm' is unchanged by an A/B swap, so a swap cannot corrupt it"
else
    record_fail \
        "'norm' changed under an A/B swap: '${norm_ok}' became '${norm_swap}'"
fi

if [[ "${norm_ok}" != "${norm_trans}" ]]; then
    record_pass "'norm' sees an 'N'/'L' transposition, the error that would corrupt it"
else
    record_fail "'norm' cannot see an 'N'/'L' transposition; 'k = L / N' is unpinned"
fi

if grep -qF -- "--pseudo ${exp_pseudo}" "${log_edg}"; then
    record_pass "plumbing equivalence: the derived pseudocount matches a hand run"
else
    record_fail \
        "plumbing equivalence: expected '--pseudo ${exp_pseudo}' in the" \
        "ratio call; see $(print_relpath "${log_edg}")"
fi

# '--prior_count' reaches the derivation; the default run above already pins
# that omitting it leaves 'compute_pseudo' its own default.
exp_pc1="$(
    "${TEST_MANAGED_PYTHON}" \
        -m protocol_chipseq_signal_norm.cli.compute_pseudo \
        --method edger \
        --typ_sig count \
        --prior_count 1 \
        --fil_A "${trk_A}" \
        --fil_B "${trk_B}" \
        --n_frg_A "$(< "$(derive_report_path "${trk_A}" n_frg)")" \
        --n_frg_B "$(< "$(derive_report_path "${trk_B}" n_frg)")" \
        --n_ovlp_A "$(< "$(derive_report_path "${trk_A}" n_ovlp)")" \
        --n_ovlp_B "$(< "$(derive_report_path "${trk_B}" n_ovlp)")" \
        2>/dev/null
)"

if [[ -n "${exp_pc1}" && "${exp_pc1}" != "${exp_pseudo}" ]]; then
    record_pass "'--prior_count 1' and the default derive different values"
else
    record_fail \
        "'--prior_count 1' derives '${exp_pc1}', the same as the default" \
        "'${exp_pseudo}', so the forwarding check below has no power"
fi

"${TEST_BASH}" "${ROOT_REPO}/bin/execute_compute_signal.sh" \
    --env_nam "${env_nam}" \
    --threads 1 \
    --mode ratio \
    --method log2 \
    --typ_sig count \
    --prior_count 1 \
    --csv_fil_A "${trk_A}" \
    --csv_fil_B "${trk_B}" \
    --dir_out "${dir_edg}" \
    --typ_out bedGraph \
    --prefix pc1 \
    --dir_eo "${dir_err}" \
    --nam_job "edger_pc1" \
    --max_job 1 \
    --csv_scl_fct NA \
    --csv_dep_min NA \
    --csv_pseudo edger > /dev/null 2>&1 || true

log_pc1="${dir_err}/edger_pc1.pc1_tiny_se.stderr.txt"

if grep -qF -- "--pseudo ${exp_pc1}" "${log_pc1}" 2>/dev/null; then
    record_pass \
        "'--prior_count 1' reaches 'compute_pseudo' and its value is used"
else
    record_fail \
        "'--prior_count 1': expected '--pseudo ${exp_pc1}' in the ratio" \
        "call; see $(print_relpath "${log_pc1}")"
fi

# Print the 'A:B' spec 'compute_pseudo' derives for one pair from the counts
# beside its tracks, with an optional prior count.
function exp_edger() {
    local fil_A="${1}" fil_B="${2}" pc="${3:-}"
    local -a opt=()

    if [[ -n "${pc}" ]]; then opt=( --prior_count "${pc}" ); fi

    "${TEST_MANAGED_PYTHON}" \
        -m protocol_chipseq_signal_norm.cli.compute_pseudo \
        --method edger \
        --typ_sig count \
        --fil_A "${fil_A}" \
        --fil_B "${fil_B}" \
        --n_frg_A "$(< "$(derive_report_path "${fil_A}" n_frg)")" \
        --n_frg_B "$(< "$(derive_report_path "${fil_B}" n_frg)")" \
        --n_ovlp_A "$(< "$(derive_report_path "${fil_A}" n_ovlp)")" \
        --n_ovlp_B "$(< "$(derive_report_path "${fil_B}" n_ovlp)")" \
        "${opt[@]}" \
        2>/dev/null
}

# Mixed lists: each pair must get its own element, checked in both orders so an
# element applied to the wrong pair cannot pass.
arr_pair_A=( "${trk_A}" "${trk_B}" )
arr_pair_B=( "${trk_B}" "${trk_A}" )
arr_case_mix=(
    "lit_edg|1:1,edger|1"
    "edg_lit|edger,1:1|1"
    "lit_edg_def|1:1,edger|"
)

for case_mix in "${arr_case_mix[@]}"; do
    IFS='|' read -r tag_mix csv_mix pc_mix <<< "${case_mix}"
    IFS=',' read -r -a arr_mix <<< "${csv_mix}"
    arr_pc_mix=()

    if [[ -n "${pc_mix}" ]]; then arr_pc_mix=( --prior_count "${pc_mix}" ); fi

    rc_mix=0
    "${TEST_BASH}" "${ROOT_REPO}/bin/execute_compute_signal.sh" \
        --env_nam "${env_nam}" \
        --threads 1 \
        --mode ratio \
        --method log2 \
        --typ_sig count \
        "${arr_pc_mix[@]}" \
        --csv_fil_A "${trk_A},${trk_B}" \
        --csv_fil_B "${trk_B},${trk_A}" \
        --dir_out "${dir_edg}" \
        --typ_out bedGraph \
        --prefix "mix_${tag_mix}" \
        --dir_eo "${dir_err}" \
        --nam_job "edger_mix_${tag_mix}" \
        --max_job 1 \
        --csv_scl_fct NA,NA \
        --csv_dep_min NA,NA \
        --csv_pseudo "${csv_mix}" > /dev/null 2>&1 || rc_mix=$?

    # The driver's exit status does not see a failed pair, so each pair's track
    # and ratio call are checked directly.
    n_bad=0
    for idx in 0 1; do
        nam="mix_${tag_mix}_$(basename "${arr_pair_A[idx]}" .bedGraph)"
        log_mix="${dir_err}/edger_mix_${tag_mix}.${nam}.stderr.txt"
        exp_mix="${arr_mix[idx]}"

        if [[ "${exp_mix}" == "edger" ]]; then
            exp_mix="$(
                exp_edger "${arr_pair_A[idx]}" "${arr_pair_B[idx]}" "${pc_mix}"
            )"
        fi

        if [[
            ! -s "${dir_edg}/${nam}.bedGraph"
            || -z "${exp_mix}"
        ]] || ! grep -qF -- "--pseudo ${exp_mix} " "${log_mix}" 2>/dev/null
        then
            n_bad=$(( n_bad + 1 ))
        fi
    done

    if [[ "${rc_mix}" -eq 0 && "${n_bad}" -eq 0 ]]; then
        record_pass \
            "'${csv_mix}' with '--prior_count ${pc_mix:-<default>}' gives" \
            "each pair its own element"
    else
        record_fail \
            "'${csv_mix}' with '--prior_count ${pc_mix:-<default>}':" \
            "${n_bad} of 2 pairs missing a track or the expected '--pseudo'" \
            "(exit ${rc_mix}); see $(print_relpath "${dir_err}")"
    fi
done
unset arr_case_mix case_mix tag_mix csv_mix pc_mix arr_mix arr_pc_mix rc_mix
unset n_bad idx nam log_mix exp_mix

# A single value applies to every pair: one 'edger' still derives each pair's
# own pseudocount, and one scale factor and floor reach both.
for case_bc in \
    "edg|--csv_pseudo edger --prior_count 1" \
    "scl|--csv_scl_fct 2:1 --csv_dep_min 0.5 --csv_pseudo NA"
do
    IFS='|' read -r tag_bc opt_bc <<< "${case_bc}"
    read -r -a arr_opt_bc <<< "${opt_bc}"

    "${TEST_BASH}" "${ROOT_REPO}/bin/execute_compute_signal.sh" \
        --env_nam "${env_nam}" \
        --threads 1 \
        --mode ratio \
        --method log2 \
        --typ_sig count \
        --csv_fil_A "${trk_A},${trk_B}" \
        --csv_fil_B "${trk_B},${trk_A}" \
        --dir_out "${dir_edg}" \
        --typ_out bedGraph \
        --prefix "bc_${tag_bc}" \
        --dir_eo "${dir_err}" \
        --nam_job "edger_bc_${tag_bc}" \
        --max_job 1 \
        "${arr_opt_bc[@]}" > /dev/null 2>&1 || true

    n_bad=0
    for idx in 0 1; do
        nam="bc_${tag_bc}_$(basename "${arr_pair_A[idx]}" .bedGraph)"
        log_bc="${dir_err}/edger_bc_${tag_bc}.${nam}.stderr.txt"

        if [[ "${tag_bc}" == "edg" ]]; then
            want="--pseudo $(
                exp_edger "${arr_pair_A[idx]}" "${arr_pair_B[idx]}" 1
            ) "
        else
            want="--scl_fct 2:1 --dep_min 0.5 "
        fi

        if [[ ! -s "${dir_edg}/${nam}.bedGraph" ]] \
            || ! grep -qF -- "${want}" "${log_bc}" 2>/dev/null
        then
            n_bad=$(( n_bad + 1 ))
        fi
    done

    if [[ "${n_bad}" -eq 0 ]]; then
        record_pass "a single '${opt_bc}' applies to both pairs"
    else
        record_fail \
            "a single '${opt_bc}' did not reach ${n_bad} of 2 pairs; see" \
            "$(print_relpath "${dir_err}")"
    fi
done
unset case_bc tag_bc opt_bc arr_opt_bc n_bad idx nam log_bc want

# A single literal pseudocount is not applied to every pair.
rc_bl=0
out_bl="$(
    "${TEST_BASH}" "${ROOT_REPO}/bin/execute_compute_signal.sh" \
        --env_nam "${env_nam}" \
        --threads 1 \
        --mode ratio \
        --method log2 \
        --csv_fil_A "${trk_A},${trk_B}" \
        --csv_fil_B "${trk_B},${trk_A}" \
        --dir_out "${dir_edg}" \
        --typ_out bedGraph \
        --prefix bc_lit \
        --dir_eo "${dir_err}" \
        --nam_job "edger_bc_lit" \
        --max_job 1 \
        --csv_pseudo 1:1 2>&1
)" || rc_bl=$?

if [[
    "${rc_bl}" -ne 0
    && "${out_bl}" == *"single literal '--csv_pseudo'"*
    && ! -e "${dir_edg}/bc_lit_tiny_se.bedGraph"
]]; then
    record_pass "a single literal '--csv_pseudo' is not applied to all pairs"
else
    record_fail \
        "a single literal '--csv_pseudo' was applied to all pairs" \
        "(exit ${rc_bl})"
fi
unset rc_bl out_bl

# Separability: the pseudocount is added before scaling, so a scale factor is a
# pure additive offset in log2 space and never interacts with the prior.
for arm in plain scaled; do
    csv_scl="NA"

    if [[ "${arm}" == "scaled" ]]; then
        csv_scl="2:1"
    fi

    "${TEST_BASH}" "${ROOT_REPO}/bin/execute_compute_signal.sh" \
        --env_nam "${env_nam}" \
        --threads 1 \
        --mode ratio \
        --method log2 \
        --typ_sig count \
        --csv_fil_A "${trk_A}" \
        --csv_fil_B "${trk_B}" \
        --dir_out "${dir_edg}" \
        --typ_out bedGraph \
        --prefix "sep_${arm}" \
        --dir_eo "${dir_err}" \
        --nam_job "edger_sep_${arm}" \
        --max_job 1 \
        --dp 15 \
        --csv_scl_fct "${csv_scl}" \
        --csv_dep_min NA \
        --csv_pseudo edger > /dev/null 2>&1 || true
done
unset arm csv_scl

sep_plain="${dir_edg}/sep_plain_tiny_se.bedGraph"
sep_scaled="${dir_edg}/sep_scaled_tiny_se.bedGraph"

# The combination itself: a derived pseudocount alongside a scale factor.
assert_file_nonempty \
    "${sep_scaled}" \
    "'--csv_pseudo edger' with '--csv_scl_fct' produces a track"

if [[ -s "${sep_plain}" && -s "${sep_scaled}" ]]; then
    read -r n_chk n_bad < <(
        paste "${sep_plain}" "${sep_scaled}" \
            | awk -F'\t' '
                function abs(x) { return x < 0 ? -x : x }

                $4 !~ /^-?[0-9]/ || $8 !~ /^-?[0-9]/ { next }

                {
                    n_chk += 1

                    if (abs($8 - ($4 + 1.0)) > 1e-9) {
                        n_bad += 1
                    }
                }

                END { print (n_chk + 0) "\t" (n_bad + 0) }
            '
    )

    if [[ "${n_chk}" -gt 0 && "${n_bad}" -eq 0 ]]; then
        record_pass \
            "separability: 'scl_fct 2:1' shifts all ${n_chk} finite bins by" \
            "exactly 1 in log2 space"
    else
        record_fail \
            "separability: ${n_bad} of ${n_chk} finite bins do not satisfy" \
            "'log2 scaled == log2 plain + 1';" \
            "see $(print_relpath "${sep_scaled}")"
    fi

    unset n_chk n_bad
fi

# A missing report must name what is missing, not fail as a bare open error.
dir_miss="${tmp}/edger_missing"
mkdir -p "${dir_miss}"
cp "${trk_A}" "${dir_miss}/a.bedGraph"
cp "${trk_B}" "${dir_miss}/b.bedGraph"

out_miss="$(
    "${TEST_BASH}" "${ROOT_REPO}/bin/execute_compute_signal.sh" \
        --env_nam "${env_nam}" \
        --threads 1 \
        --mode ratio \
        --method log2 \
        --typ_sig count \
        --csv_fil_A "${dir_miss}/a.bedGraph" \
        --csv_fil_B "${dir_miss}/b.bedGraph" \
        --dir_out "${dir_miss}" \
        --typ_out bedGraph \
        --dir_eo "${dir_err}" \
        --nam_job "edger_missing" \
        --max_job 1 \
        --csv_scl_fct NA \
        --csv_dep_min NA \
        --csv_pseudo edger 2>&1
)" || true

# The descriptor is derived from the method and sample, so glob rather than
# guess it.
txt_miss="${out_miss}$(
    cat "${dir_err}"/edger_missing.*.stderr.txt 2>/dev/null || true
)"

if [[
    "${txt_miss}" == *"n_frg"* && "${txt_miss}" == *"--csv_pseudo edger"*
]]; then
    record_pass "a missing count names the report and the flag that needed it"
else
    record_fail \
        "a missing count did not name the report and the flag; see" \
        "$(print_relpath "${dir_err}")/edger_missing.*.stderr.txt"
fi

# A failed first pair must fail the serial run, and the second pair must still
# run; before, only the last task's status survived.
rc_ser=0
out_ser="$(
    "${TEST_BASH}" "${ROOT_REPO}/bin/execute_compute_signal.sh" \
        --env_nam "${env_nam}" \
        --threads 1 \
        --mode ratio \
        --method log2 \
        --typ_sig count \
        --csv_fil_A "${dir_miss}/a.bedGraph,${trk_A}" \
        --csv_fil_B "${dir_miss}/b.bedGraph,${trk_B}" \
        --dir_out "${dir_miss}" \
        --typ_out bedGraph \
        --prefix serial \
        --dir_eo "${dir_err}" \
        --nam_job "serial_status" \
        --max_job 1 \
        --csv_scl_fct NA,NA \
        --csv_dep_min NA,NA \
        --csv_pseudo "edger,1:1" 2>&1
)" || rc_ser=$?

if [[
    "${rc_ser}" -ne 0
    && "${out_ser}" == *"1 of 2 task(s) failed"*
    && ! -e "${dir_miss}/serial_a.bedGraph"
    && -s "${dir_miss}/serial_tiny_se.bedGraph"
]]; then
    record_pass "a failed first pair fails the serial run; the second still runs"
else
    record_fail \
        "a failed first pair was not reported by the serial run, or the" \
        "second pair did not run (exit ${rc_ser})"
fi
unset rc_ser out_ser

# An invalid argument stops at its first check rather than running on to fail
# somewhere unrelated.
rc_thr=0
out_thr="$(
    "${TEST_BASH}" "${ROOT_REPO}/bin/execute_compute_signal.sh" \
        --env_nam "${env_nam}" \
        --threads 0 \
        --mode ratio \
        --method linear \
        --csv_fil_A "${fil_A}" \
        --csv_fil_B "${fil_B}" \
        --dir_out "${dir_out}" \
        --typ_out bdg \
        --dir_eo "${dir_err}" \
        --max_job 1 2>&1
)" || rc_thr=$?

if [[
    "${rc_thr}" -ne 0
    && "${out_thr}" == *"'--threads' was assigned '0'"*
    && "${out_thr}" != *"print_parallel_info"*
]]; then
    record_pass "an invalid '--threads' stops at its own check"
else
    record_fail \
        "an invalid '--threads' did not stop at its own check" \
        "(exit ${rc_thr})"
fi
unset rc_thr out_thr

# '--typ_sig' belongs to ratio mode; '--method' owns that meaning in signal.
out_sig="$(
    "${TEST_BASH}" "${ROOT_REPO}/bin/execute_compute_signal.sh" \
        --env_nam "${env_nam}" \
        --threads 1 \
        --mode signal \
        --method norm \
        --typ_sig norm \
        --csv_fil_in "${ROOT_REPO}/tests/fixtures/compute_signal/bam/se/tiny_se.bam" \
        --dir_out "${dir_edg}" \
        --dir_eo "${dir_err}" \
        --siz_bin 10 2>&1
)" || true

if [[
    "${out_sig}" == *"'--typ_sig'"* && "${out_sig}" == *"'--method'"*
]]; then
    record_pass "'--typ_sig' is rejected in signal mode, naming '--method'"
else
    record_fail "'--typ_sig' was not rejected in signal mode by name"
fi

# '--prior_count' is refused wherever it could change no output, so a sweep
# over it cannot silently return identical tracks.
bam_se="${ROOT_REPO}/tests/fixtures/compute_signal/bam/se/tiny_se.bam"

rc_pc=0
out_pc="$(
    "${TEST_BASH}" "${ROOT_REPO}/bin/execute_compute_signal.sh" \
        --env_nam "${env_nam}" \
        --threads 1 \
        --mode signal \
        --method norm \
        --prior_count 1 \
        --csv_fil_in "${bam_se}" \
        --dir_out "${dir_edg}" \
        --dir_eo "${dir_err}" \
        --siz_bin 10 2>&1
)" || rc_pc=$?

if [[
    "${rc_pc}" -ne 0
    && "${out_pc}" == *"'--prior_count' is for '--mode ratio'"*
]]; then
    record_pass "'--prior_count' is rejected in signal mode"
else
    record_fail \
        "'--prior_count' was not rejected in signal mode (exit ${rc_pc})"
fi

rc_pc=0
out_pc="$(
    "${TEST_BASH}" "${ROOT_REPO}/bin/execute_compute_signal.sh" \
        --env_nam "${env_nam}" \
        --threads 1 \
        --mode ratio \
        --method log2 \
        --prior_count 1 \
        --csv_fil_A "${trk_A}" \
        --csv_fil_B "${trk_B}" \
        --dir_out "${dir_edg}" \
        --typ_out bedGraph \
        --prefix pc_lit \
        --dir_eo "${dir_err}" \
        --nam_job "edger_pc_lit" \
        --max_job 1 \
        --csv_scl_fct NA \
        --csv_dep_min NA \
        --csv_pseudo 1:1 2>&1
)" || rc_pc=$?

if [[
    "${rc_pc}" -ne 0
    && "${out_pc}" == *"'--prior_count' applies only to '--csv_pseudo edger'"*
    && ! -e "${dir_edg}/pc_lit_tiny_se.bedGraph"
]]; then
    record_pass "'--prior_count' is rejected when no element is 'edger'"
else
    record_fail \
        "'--prior_count' with only literal pseudocounts was not rejected" \
        "before writing a track (exit ${rc_pc})"
fi

# With '--csv_pseudo' omitted, every element is the 'NA' sentinel.
rc_pc=0
out_pc="$(
    "${TEST_BASH}" "${ROOT_REPO}/bin/execute_compute_signal.sh" \
        --env_nam "${env_nam}" \
        --threads 1 \
        --mode ratio \
        --method log2 \
        --prior_count 1 \
        --csv_fil_A "${trk_A}" \
        --csv_fil_B "${trk_B}" \
        --dir_out "${dir_edg}" \
        --typ_out bedGraph \
        --prefix pc_none \
        --dir_eo "${dir_err}" \
        --nam_job "edger_pc_none" \
        --max_job 1 2>&1
)" || rc_pc=$?

if [[
    "${rc_pc}" -ne 0
    && "${out_pc}" == *"'--prior_count' applies only to '--csv_pseudo edger'"*
    && ! -e "${dir_edg}/pc_none_tiny_se.bedGraph"
]]; then
    record_pass "'--prior_count' is rejected when '--csv_pseudo' is omitted"
else
    record_fail \
        "'--prior_count' without '--csv_pseudo' was not rejected before" \
        "writing a track (exit ${rc_pc})"
fi

rc_pc=0
out_pc="$(
    "${TEST_BASH}" "${ROOT_REPO}/bin/execute_compute_signal.sh" \
        --env_nam "${env_nam}" \
        --threads 1 \
        --mode ratio \
        --method log2 \
        --typ_sig count \
        --prior_count two \
        --csv_fil_A "${trk_A}" \
        --csv_fil_B "${trk_B}" \
        --dir_out "${dir_edg}" \
        --typ_out bedGraph \
        --prefix pc_bad \
        --dir_eo "${dir_err}" \
        --nam_job "edger_pc_bad" \
        --max_job 1 \
        --csv_scl_fct NA \
        --csv_dep_min NA \
        --csv_pseudo edger 2>&1
)" || rc_pc=$?

if [[
    "${rc_pc}" -ne 0
    && "${out_pc}" == *"'--prior_count' was assigned 'two'"*
    && ! -e "${dir_edg}/pc_bad_tiny_se.bedGraph"
]]; then
    record_pass "a non-numeric '--prior_count' is rejected before any work"
else
    record_fail \
        "a non-numeric '--prior_count' was not rejected (exit ${rc_pc})"
fi
unset bam_se rc_pc out_pc

# A bare literal would regularize file A alone and leave file B's zero bins
# undefined, so it is refused before any work.
rc_bare=0
out_bare="$(
    "${TEST_BASH}" "${ROOT_REPO}/bin/execute_compute_signal.sh" \
        --env_nam "${env_nam}" \
        --threads 1 \
        --mode ratio \
        --method log2 \
        --csv_fil_A "${fil_A}" \
        --csv_fil_B "${fil_B}" \
        --dir_out "${dir_out}" \
        --typ_out bdg \
        --prefix exec_bare \
        --dir_eo "${dir_err}" \
        --nam_job "test_execute_compute_ratio_bare" \
        --max_job 1 \
        --csv_pseudo 1 2>&1
)" || rc_bare=$?

if [[
    "${rc_bare}" -ne 0
    && "${out_bare}" == *"needs 'A:B'"*
    && ! -e "${dir_out}/exec_bare_ratio_A.bdg"
]]; then
    record_pass "a bare '--csv_pseudo' literal is refused, naming 'A:B'"
else
    record_fail \
        "a bare '--csv_pseudo' literal was not refused before any work" \
        "(exit ${rc_bare})"
fi
unset rc_bare out_bare

# '--prior_count' belongs to ratio mode, so coord mode refuses it too.
bam_crd="${ROOT_REPO}/tests/fixtures/compute_signal/bam/se/tiny_se.bam"

rc_crd=0
out_crd="$(
    "${TEST_BASH}" "${ROOT_REPO}/bin/execute_compute_signal.sh" \
        --env_nam "${env_nam}" \
        --threads 1 \
        --mode coord \
        --prior_count 1 \
        --csv_fil_in "${bam_crd}" \
        --dir_out "${dir_out}" \
        --typ_out bed \
        --prefix exec_crd_pc \
        --dir_eo "${dir_err}" 2>&1
)" || rc_crd=$?

if [[
    "${rc_crd}" -ne 0
    && "${out_crd}" == *"'--prior_count' is for '--mode ratio'"*
]]; then
    record_pass "'--prior_count' is rejected in coord mode"
else
    record_fail \
        "'--prior_count' was not rejected in coord mode (exit ${rc_crd})"
fi
unset bam_crd rc_crd out_crd

# A non-numeric side of 'A:B' must stop the run before dispatch; the exit
# status alone cannot show that, so also require that no task log exists.
rc_nan=0
out_nan="$(
    "${TEST_BASH}" "${ROOT_REPO}/bin/execute_compute_signal.sh" \
        --env_nam "${env_nam}" \
        --threads 1 \
        --mode ratio \
        --method log2 \
        --csv_fil_A "${fil_A}" \
        --csv_fil_B "${fil_B}" \
        --dir_out "${dir_out}" \
        --typ_out bdg \
        --prefix exec_nan_psd \
        --dir_eo "${dir_err}" \
        --nam_job "test_execute_compute_ratio_nan_psd" \
        --max_job 1 \
        --csv_pseudo x:1 2>&1
)" || rc_nan=$?

if [[
    "${rc_nan}" -ne 0
    && "${out_nan}" == *"'--csv_pseudo' was assigned 'x'"*
    && ! -e "${dir_out}/exec_nan_psd_ratio_A.bdg"
    && -z "$(
        find "${dir_err}" -name 'test_execute_compute_ratio_nan_psd.*' \
            2>/dev/null
    )"
]]; then
    record_pass "a non-numeric side of an 'A:B' pseudocount stops the run"
else
    record_fail \
        "a non-numeric side of an 'A:B' pseudocount did not stop the run" \
        "(exit ${rc_nan})"
fi
unset rc_nan out_nan


finish
