#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_submit_qsort_bam.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in design, development,
# and documentation, with all output reviewed, edited, and approved by the
# author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="submit qsort BAM/CRAM"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


#  Define fixture and output paths for local serial qsort tests
dir_fx="${ROOT_REPO}/tests/fixtures/calculate_scaling_factor"
in_bam="${dir_fx}/bam/pe/IP_WT_G1_Hho1_6336.sc.bam"
in_cram="${dir_fx}/cram/pe/IP_WT_G1_Hho1_6336.sc.cram"
ref_fa="${dir_fx}/reference/tiny.fa"

tmp="${TEST_DIR_TMP}/submit_qsort_bam"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/qsort_bam"

out_bam="${dir_out}/IP_WT_G1_Hho1_6336.sc.qnam.bam"
out_cram="${dir_out}/IP_WT_G1_Hho1_6336.sc.qnam.cram"

log_help="${dir_log}/submit_qsort_bam_help.log"
log_positional="${dir_log}/submit_qsort_bam_positional_fails.log"
log_unknown="${dir_log}/submit_qsort_bam_unknown_option.log"
log_missing_ref="${dir_log}/submit_qsort_bam_cram_missing_ref.log"
log_serial="${dir_log}/submit_qsort_bam_serial.log"
log_qc_bam="${dir_log}/submit_qsort_bam_bam_quickcheck.log"
log_qc_cram="${dir_log}/submit_qsort_bam_cram_quickcheck.log"


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty \
    "${in_bam}" \
    "${in_cram}" \
    "${ref_fa}" || {
    finish
    exit $?
}


#  Help should exit successfully
if \
    run_capture \
        "submit_qsort_bam.sh --help" \
        "${log_help}" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/submit_qsort_bam.sh" \
            --help
then
    record_pass "submit_qsort_bam.sh --help exits 0"
else
    record_fail \
        "submit_qsort_bam.sh --help failed;" \
        "see $(print_relpath "${log_help}")"
fi


#  Old positional invocation should fail after the keyarg migration
# shellcheck disable=SC2154
if \
    run_capture \
        "submit_qsort_bam.sh positional invocation fails" \
        "${log_positional}" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/submit_qsort_bam.sh" \
            "${env_nam}" \
            1 \
            "${in_bam}" \
            "${dir_out}" \
            "${dir_err}" \
            test_submit_qsort_positional
then
    record_fail "submit_qsort_bam.sh positional invocation passed"
else
    record_pass "submit_qsort_bam.sh positional invocation fails"
fi


#  Unknown option should fail clearly
if \
    run_capture \
        "submit_qsort_bam.sh unknown option fails" \
        "${log_unknown}" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/submit_qsort_bam.sh" \
            --not_an_option
then
    record_fail "submit_qsort_bam.sh unknown option passed"
else
    record_pass "submit_qsort_bam.sh unknown option fails"
fi


#  CRAM input without a reference should fail before execution
if \
    run_capture \
        "submit_qsort_bam.sh CRAM missing ref fails" \
        "${log_missing_ref}" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/submit_qsort_bam.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/bin" \
            --threads 1 \
            --csv_fil_in "${in_cram}" \
            --dir_out "${dir_out}" \
            --dir_eo "${dir_err}" \
            --nam_job test_submit_qsort_cram_missing_ref
then
    record_fail "submit_qsort_bam.sh CRAM without --ref_fa passed"
else
    record_pass "submit_qsort_bam.sh CRAM without --ref_fa fails"
fi


#  Serial BAM and CRAM qsort should write one output per input
if \
    run_capture \
        "submit_qsort_bam.sh serial BAM/CRAM" \
        "${log_serial}" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/submit_qsort_bam.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/bin" \
            --threads 1 \
            --csv_fil_in "${in_bam},${in_cram}" \
            --ref_fa "${ref_fa}" \
            --dir_out "${dir_out}" \
            --dir_eo "${dir_err}" \
            --nam_job test_submit_qsort_serial
then
    record_pass "submit_qsort_bam.sh serial BAM/CRAM exits 0"
else
    record_fail \
        "submit_qsort_bam.sh serial BAM/CRAM failed;" \
        "see $(print_relpath "${log_serial}")"
fi

assert_file_nonempty \
    "${out_bam}" \
    "qsort BAM output"

assert_file_nonempty \
    "${out_cram}" \
    "qsort CRAM output"

if [[ -s "${out_bam}" ]]; then
    if \
        run_capture \
            "quickcheck qsort BAM output" \
            "${log_qc_bam}" \
            run_samtools quickcheck "${out_bam}"
    then
        record_pass "qsort BAM output passes samtools quickcheck"
    else
        record_fail "qsort BAM output fails samtools quickcheck"
    fi
fi

if [[ -s "${out_cram}" ]]; then
    if \
        run_capture \
            "quickcheck qsort CRAM output" \
            "${log_qc_cram}" \
            run_samtools quickcheck "${out_cram}"
    then
        record_pass "qsort CRAM output passes samtools quickcheck"
    else
        record_fail "qsort CRAM output fails samtools quickcheck"
    fi
fi

finish
