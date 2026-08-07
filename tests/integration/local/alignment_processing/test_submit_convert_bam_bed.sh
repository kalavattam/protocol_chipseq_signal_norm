#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_submit_convert_bam_bed.sh
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

TEST_NAME="submit convert BAM/CRAM to BED"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


#  Define fixture and output paths for local serial conversion tests
dir_fx="${ROOT_REPO}/tests/fixtures/calculate_scaling_factor"
in_bam="${dir_fx}/bam/pe/IP_WT_G1_Hho1_6336.sc.bam"
in_cram="${dir_fx}/cram/pe/IP_WT_G1_Hho1_6337.sc.cram"
ref_fa="${dir_fx}/reference/tiny.fa"

tmp="${TEST_DIR_TMP}/submit_convert_bam_bed"
dir_qsort="${tmp}/qsort"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/convert_bam_bed"

q_bam="${dir_qsort}/IP_WT_G1_Hho1_6336.sc.qnam.bam"
q_cram="${dir_qsort}/IP_WT_G1_Hho1_6337.sc.qnam.cram"

out_bam_bed="${dir_out}/IP_WT_G1_Hho1_6336.sc.bed.gz"
out_cram_bed="${dir_out}/IP_WT_G1_Hho1_6337.sc.bed.gz"

vw_bam_bed="${dir_out}/IP_WT_G1_Hho1_6336.sc.bed.txt"
vw_cram_bed="${dir_out}/IP_WT_G1_Hho1_6337.sc.bed.txt"
vw_awk_bam_bed="${dir_out}/IP_WT_G1_Hho1_6336.sc.awk.bed.txt"

log_help="${dir_log}/submit_convert_bam_bed_help.log"
log_positional="${dir_log}/submit_convert_bam_bed_positional_fails.log"
log_unknown="${dir_log}/submit_convert_bam_bed_unknown_option.log"
log_missing_ref="${dir_log}/submit_convert_bam_bed_cram_missing_ref.log"
log_py_awk="${dir_log}/submit_convert_bam_bed_py_and_awk.log"
log_prepare="${dir_log}/submit_convert_bam_bed_prepare_qsort.log"
log_python="${dir_log}/submit_convert_bam_bed_python_serial.log"
log_awk="${dir_log}/submit_convert_bam_bed_awk_bam.log"


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${dir_qsort}" "${dir_out}" "${dir_err}" "${dir_log}"

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
        "submit_convert_bam_bed.sh --help" \
        "${log_help}" \
        "${TEST_BASH}" \
        "${ROOT_REPO}/bin/submit_convert_bam_bed.sh" \
            --help
then
    record_pass "submit_convert_bam_bed.sh --help exits 0"
else
    record_fail \
        "submit_convert_bam_bed.sh --help failed;" \
        "see $(print_relpath "${log_help}")"
fi


#  Old positional invocation should fail after the keyarg migration
# shellcheck disable=SC2154
if \
    run_capture \
        "submit_convert_bam_bed.sh positional invocation fails" \
        "${log_positional}" \
        "${TEST_BASH}" \
        "${ROOT_REPO}/bin/submit_convert_bam_bed.sh" \
            "${env_nam}" \
            1 \
            "${in_bam}" \
            "${ROOT_REPO}/src/protocol_chipseq_signal_norm/cli/compute_signal.py" \
            "${dir_out}" \
            "${dir_err}" \
            test_submit_convert_positional \
            false
then
    record_fail \
        "submit_convert_bam_bed.sh positional invocation passed"
else
    record_pass \
        "submit_convert_bam_bed.sh positional invocation fails"
fi


#  Unknown option should fail clearly
if \
    run_capture \
        "submit_convert_bam_bed.sh unknown option fails" \
        "${log_unknown}" \
        "${TEST_BASH}" \
        "${ROOT_REPO}/bin/submit_convert_bam_bed.sh" \
            --not_an_option
then
    record_fail "submit_convert_bam_bed.sh unknown option passed"
else
    record_pass "submit_convert_bam_bed.sh unknown option fails"
fi


#  CRAM input without a reference should fail before execution
if \
    run_capture \
        "submit_convert_bam_bed.sh CRAM missing ref fails" \
        "${log_missing_ref}" \
        "${TEST_BASH}" \
        "${ROOT_REPO}/bin/submit_convert_bam_bed.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/bin" \
            --threads 1 \
            --csv_fil_in "${in_cram}" \
            --dir_out "${dir_out}" \
            --dir_eo "${dir_err}" \
            --nam_job test_submit_convert_cram_missing_ref
then
    record_fail \
        "submit_convert_bam_bed.sh CRAM without --ref_fa passed"
else
    record_pass \
        "submit_convert_bam_bed.sh CRAM without --ref_fa fails"
fi


#  Explicit Python script and AWK branch are mutually exclusive
if \
    run_capture \
        "submit_convert_bam_bed.sh Python and AWK conflict" \
        "${log_py_awk}" \
        "${TEST_BASH}" \
        "${ROOT_REPO}/bin/submit_convert_bam_bed.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/bin" \
            --threads 1 \
            --csv_fil_in "${in_bam}" \
            --pth_scr_py "${ROOT_REPO}/src/protocol_chipseq_signal_norm/cli/compute_signal.py" \
            --dir_out "${dir_out}" \
            --dir_eo "${dir_err}" \
            --nam_job test_submit_convert_py_and_awk \
            --use_awk
then
    record_fail \
        "submit_convert_bam_bed.sh --pth_scr_py with --use_awk passed"
else
    record_pass \
        "submit_convert_bam_bed.sh --pth_scr_py with --use_awk fails"
fi


#  Prepare qname-sorted temporary inputs for conversion tests
if \
    run_capture \
        "prepare qname-sorted BAM/CRAM conversion inputs" \
        "${log_prepare}" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/submit_qsort_bam.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/bin" \
            --threads 1 \
            --csv_fil_in "${in_bam},${in_cram}" \
            --ref_fa "${ref_fa}" \
            --dir_out "${dir_qsort}" \
            --dir_eo "${dir_err}" \
            --nam_job test_submit_convert_prepare_qsort
then
    record_pass "prepared qname-sorted BAM/CRAM conversion inputs"
else
    record_fail \
        "failed to prepare qname-sorted conversion inputs;" \
        "see $(print_relpath "${log_prepare}")"
fi

require_files_nonempty \
    "${q_bam}" \
    "${q_cram}" || {
    finish
    exit $?
}


#  Serial Python conversion should write BED.GZ output for BAM and CRAM inputs
if \
    run_capture \
        "submit_convert_bam_bed.sh Python BAM/CRAM" \
        "${log_python}" \
        "${TEST_BASH}" \
        "${ROOT_REPO}/bin/submit_convert_bam_bed.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/bin" \
            --threads 1 \
            --csv_fil_in "${q_bam},${q_cram}" \
            --ref_fa "${ref_fa}" \
            --pth_scr_py "${ROOT_REPO}/src/protocol_chipseq_signal_norm/cli/compute_signal.py" \
            --dir_out "${dir_out}" \
            --dir_eo "${dir_err}" \
            --nam_job test_submit_convert_python
then
    record_pass "submit_convert_bam_bed.sh Python BAM/CRAM exits 0"
else
    record_fail \
        "submit_convert_bam_bed.sh Python BAM/CRAM failed;" \
        "see $(print_relpath "${log_python}")"
fi

assert_file_nonempty \
    "${out_bam_bed}" \
    "convert BAM BED.GZ output"

assert_file_nonempty \
    "${out_cram_bed}" \
    "convert CRAM BED.GZ output"

if [[ -s "${out_bam_bed}" ]]; then
    if \
        run_capture \
            "gzip view convert BAM BED.GZ" \
            "${vw_bam_bed}" \
            gzip -cd "${out_bam_bed}"
    then
        record_pass "convert BAM BED.GZ decompresses"
    else
        record_fail "convert BAM BED.GZ does not decompress"
    fi

    assert_pattern_found \
        "${vw_bam_bed}" \
        $'^[^\t][^\t]*\t[0-9][0-9]*\t[0-9][0-9]*\t[0-9][0-9]*$' \
        "convert BAM BED.GZ has BED-like rows"
fi

if [[ -s "${out_cram_bed}" ]]; then
    if \
        run_capture \
            "gzip view convert CRAM BED.GZ" \
            "${vw_cram_bed}" \
            gzip -cd "${out_cram_bed}"
    then
        record_pass "convert CRAM BED.GZ decompresses"
    else
        record_fail "convert CRAM BED.GZ does not decompress"
    fi

    assert_pattern_found \
        "${vw_cram_bed}" \
        $'^[^\t][^\t]*\t[0-9][0-9]*\t[0-9][0-9]*\t[0-9][0-9]*$' \
        "convert CRAM BED.GZ has BED-like rows"
fi


#  AWK branch should also handle qname-sorted PE BAM input
rm -f "${out_bam_bed}"
if \
    run_capture \
        "submit_convert_bam_bed.sh AWK BAM" \
        "${log_awk}" \
        "${TEST_BASH}" \
        "${ROOT_REPO}/bin/submit_convert_bam_bed.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/bin" \
            --threads 1 \
            --csv_fil_in "${q_bam}" \
            --dir_out "${dir_out}" \
            --dir_eo "${dir_err}" \
            --nam_job test_submit_convert_awk \
            --use_awk
then
    record_pass "submit_convert_bam_bed.sh AWK BAM exits 0"
else
    record_fail \
        "submit_convert_bam_bed.sh AWK BAM failed;" \
        "see $(print_relpath "${log_awk}")"
fi

assert_file_nonempty \
    "${out_bam_bed}" \
    "convert AWK BAM BED.GZ output"

if [[ -s "${out_bam_bed}" ]]; then
    if \
        run_capture \
            "gzip view convert AWK BAM BED.GZ" \
            "${vw_awk_bam_bed}" \
            gzip -cd "${out_bam_bed}"
    then
        record_pass "convert AWK BAM BED.GZ decompresses"
    else
        record_fail "convert AWK BAM BED.GZ does not decompress"
    fi

    assert_pattern_found \
        "${vw_awk_bam_bed}" \
        $'^[^\t][^\t]*\t[0-9][0-9]*\t[0-9][0-9]*\t[0-9][0-9]*$' \
        "convert AWK BAM BED.GZ has BED-like rows"
fi

finish
