#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_submit_convert_bam_bed.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5) were used in development.
#
# Distributed under the MIT license.


TEST_NAME="submit convert BAM/CRAM to BED"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

print_section "${TEST_NAME}"

#  Define fixture and output paths for local serial conversion smoke tests
dir_fx="${ROOT_REPO}/tests/calculate_scaling_factor/fixtures"
in_bam="${dir_fx}/bam/pe/IP_WT_G1_Hho1_6336.sc.bam"
in_cram="${dir_fx}/cram/pe/IP_WT_G1_Hho1_6337.sc.cram"
ref_fa="${dir_fx}/reference/tiny.fa"

tmp="${TEST_DIR_TMP}/submit_convert_bam_bed"
dir_qsort="${tmp}/qsort"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/convert_bam_bed"

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
log="${dir_log}/submit_convert_bam_bed_help.log"
if \
    run_capture \
        "submit_convert_bam_bed.sh --help" \
        "${log}" \
        "${TEST_BASH}" \
        "${ROOT_REPO}/scripts/submit_convert_bam_bed.sh" \
            --help
then
    record_pass "submit_convert_bam_bed.sh --help exits 0"
else
    record_fail \
        "submit_convert_bam_bed.sh --help failed;" \
        "see $(print_relpath "${log}")"
fi


#  Old positional invocation should fail after the keyarg migration
log="${dir_log}/submit_convert_bam_bed_positional_fails.log"
if \
    run_capture \
        "submit_convert_bam_bed.sh positional invocation fails" \
        "${log}" \
        "${TEST_BASH}" \
        "${ROOT_REPO}/scripts/submit_convert_bam_bed.sh" \
            "${env_nam}" \
            1 \
            "${in_bam}" \
            "${ROOT_REPO}/scripts/compute_signal.py" \
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
log="${dir_log}/submit_convert_bam_bed_unknown_option.log"
if \
    run_capture \
        "submit_convert_bam_bed.sh unknown option fails" \
        "${log}" \
        "${TEST_BASH}" \
        "${ROOT_REPO}/scripts/submit_convert_bam_bed.sh" \
            --not_an_option
then
    record_fail "submit_convert_bam_bed.sh unknown option passed"
else
    record_pass "submit_convert_bam_bed.sh unknown option fails"
fi


#  CRAM input without a reference should fail before execution
log="${dir_log}/submit_convert_bam_bed_cram_missing_ref.log"
if \
    run_capture \
        "submit_convert_bam_bed.sh CRAM missing ref fails" \
        "${log}" \
        "${TEST_BASH}" \
        "${ROOT_REPO}/scripts/submit_convert_bam_bed.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/scripts" \
            --threads 1 \
            --csv_infile "${in_cram}" \
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
log="${dir_log}/submit_convert_bam_bed_py_and_awk.log"
if \
    run_capture \
        "submit_convert_bam_bed.sh Python and AWK conflict" \
        "${log}" \
        "${TEST_BASH}" \
        "${ROOT_REPO}/scripts/submit_convert_bam_bed.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/scripts" \
            --threads 1 \
            --csv_infile "${in_bam}" \
            --pth_scr_py "${ROOT_REPO}/scripts/compute_signal.py" \
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
log="${dir_log}/submit_convert_bam_bed_prepare_qsort.log"
if \
    run_capture \
        "prepare qname-sorted BAM/CRAM conversion inputs" \
        "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/submit_qsort_bam.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/scripts" \
            --threads 1 \
            --csv_infile "${in_bam},${in_cram}" \
            --ref_fa "${ref_fa}" \
            --dir_out "${dir_qsort}" \
            --dir_eo "${dir_err}" \
            --nam_job test_submit_convert_prepare_qsort
then
    record_pass "prepared qname-sorted BAM/CRAM conversion inputs"
else
    record_fail \
        "failed to prepare qname-sorted conversion inputs;" \
        "see $(print_relpath "${log}")"
fi

q_bam="${dir_qsort}/IP_WT_G1_Hho1_6336.sc.qnam.bam"
q_cram="${dir_qsort}/IP_WT_G1_Hho1_6337.sc.qnam.cram"

require_files_nonempty \
    "${q_bam}" \
    "${q_cram}" || {
    finish
    exit $?
}


#  Serial Python conversion should write BED.GZ output for BAM and CRAM inputs
log="${dir_log}/submit_convert_bam_bed_python_serial.log"
if \
    run_capture \
        "submit_convert_bam_bed.sh Python BAM/CRAM" \
        "${log}" \
        "${TEST_BASH}" \
        "${ROOT_REPO}/scripts/submit_convert_bam_bed.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/scripts" \
            --threads 1 \
            --csv_infile "${q_bam},${q_cram}" \
            --ref_fa "${ref_fa}" \
            --pth_scr_py "${ROOT_REPO}/scripts/compute_signal.py" \
            --dir_out "${dir_out}" \
            --dir_eo "${dir_err}" \
            --nam_job test_submit_convert_python
then
    record_pass "submit_convert_bam_bed.sh Python BAM/CRAM exits 0"
else
    record_fail \
        "submit_convert_bam_bed.sh Python BAM/CRAM failed;" \
        "see $(print_relpath "${log}")"
fi

out_bam_bed="${dir_out}/IP_WT_G1_Hho1_6336.sc.bed.gz"
out_cram_bed="${dir_out}/IP_WT_G1_Hho1_6337.sc.bed.gz"

assert_file_nonempty \
    "${out_bam_bed}" \
    "convert BAM BED.GZ output"
assert_file_nonempty \
    "${out_cram_bed}" \
    "convert CRAM BED.GZ output"

if [[ -s "${out_bam_bed}" ]]; then
    vw="${dir_out}/IP_WT_G1_Hho1_6336.sc.bed.txt"

    if \
        run_capture \
            "gzip view convert BAM BED.GZ" \
            "${vw}" \
            gzip -cd "${out_bam_bed}"
    then
        record_pass "convert BAM BED.GZ decompresses"
    else
        record_fail "convert BAM BED.GZ does not decompress"
    fi

    assert_pattern_found \
        "${vw}" \
        $'^[^\t][^\t]*\t[0-9][0-9]*\t[0-9][0-9]*\t[0-9][0-9]*$' \
        "convert BAM BED.GZ has BED-like rows"
fi

if [[ -s "${out_cram_bed}" ]]; then
    vw="${dir_out}/IP_WT_G1_Hho1_6337.sc.bed.txt"

    if \
        run_capture \
            "gzip view convert CRAM BED.GZ" \
            "${vw}" \
            gzip -cd "${out_cram_bed}"
    then
        record_pass "convert CRAM BED.GZ decompresses"
    else
        record_fail "convert CRAM BED.GZ does not decompress"
    fi

    assert_pattern_found \
        "${vw}" \
        $'^[^\t][^\t]*\t[0-9][0-9]*\t[0-9][0-9]*\t[0-9][0-9]*$' \
        "convert CRAM BED.GZ has BED-like rows"
fi


#  AWK branch should also handle qname-sorted PE BAM input
rm -f "${out_bam_bed}"
log="${dir_log}/submit_convert_bam_bed_awk_bam.log"

if \
    run_capture \
        "submit_convert_bam_bed.sh AWK BAM" \
        "${log}" \
        "${TEST_BASH}" \
        "${ROOT_REPO}/scripts/submit_convert_bam_bed.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/scripts" \
            --threads 1 \
            --csv_infile "${q_bam}" \
            --dir_out "${dir_out}" \
            --dir_eo "${dir_err}" \
            --nam_job test_submit_convert_awk \
            --use_awk
then
    record_pass "submit_convert_bam_bed.sh AWK BAM exits 0"
else
    record_fail \
        "submit_convert_bam_bed.sh AWK BAM failed;" \
        "see $(print_relpath "${log}")"
fi

assert_file_nonempty \
    "${out_bam_bed}" \
    "convert AWK BAM BED.GZ output"

if [[ -s "${out_bam_bed}" ]]; then
    vw="${dir_out}/IP_WT_G1_Hho1_6336.sc.awk.bed.txt"

    if \
        run_capture \
            "gzip view convert AWK BAM BED.GZ" \
            "${vw}" \
            gzip -cd "${out_bam_bed}"
    then
        record_pass "convert AWK BAM BED.GZ decompresses"
    else
        record_fail "convert AWK BAM BED.GZ does not decompress"
    fi

    assert_pattern_found \
        "${vw}" \
        $'^[^\t][^\t]*\t[0-9][0-9]*\t[0-9][0-9]*\t[0-9][0-9]*$' \
        "convert AWK BAM BED.GZ has BED-like rows"
fi

finish
