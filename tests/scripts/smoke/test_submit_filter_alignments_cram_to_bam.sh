#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_submit_filter_alignments_cram_to_bam.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5) were used in development.
#
# Distributed under the MIT license.


TEST_NAME="submit filter-alignments CRAM to BAM"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

print_section "${TEST_NAME}"

dir_fx="${ROOT_REPO}/tests/filter_alignments/fixtures"
in_sam="${dir_fx}/sam/filter_sc_sp.sam"
ref_fa="${dir_fx}/reference/filter_sc_sp.fa"
ref_fai="${ref_fa}.fai"

tmp="${TEST_DIR_TMP}/submit_filter_alignments_cram_to_bam"
dir_in="${tmp}/input"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_missing="${tmp}/missing_ref"
dir_log="${TEST_DIR_LOG}/filter_alignments"
in_cram="${dir_in}/filter_sc_sp.cram"

rm -rf "${tmp}"
mkdir -p "${dir_in}" "${dir_out}" "${dir_err}" "${dir_missing}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

# shellcheck disable=SC2034
{
    arr_cmd_filter=(
        "${TEST_BASH}" "${ROOT_REPO}/scripts/submit_filter_alignments.sh"
            --env_nam "${env_nam}"
            --dir_scr "${ROOT_REPO}/scripts"
            --threads 1
            --csv_infile "${in_cram}"
            --dir_out "${dir_out}"
            --ref_fa "${ref_fa}"
            --err_out "${dir_err}"
    )
    nam_job_pfx="test_submit_filter_alignments_cram_to_bam"
    arr_arg_nil=()
    arr_arg_sp=( --tg --mtr --mito )
}

require_files_nonempty \
    "${in_sam}" \
    "${ref_fa}" \
    "${ref_fai}" || {
    finish
    exit $?
}


#  Build deterministic CRAM input from committed SAM and reference fixtures
log="${dir_log}/submit_filter_alignments_prepare_cram_to_bam.log"
if ! \
    build_filter_alignments_fixture_cram \
        "${in_sam}" \
        "${ref_fa}" \
        "${in_cram}" \
        "${log}" \
        "submit filter-alignments CRAM-to-BAM fixture"
then
    finish
    exit $?
fi


#  CRAM input without --ref_fa should fail clearly
log="${dir_log}/submit_filter_alignments_cram_to_bam_missing_ref.log"
if \
    run_capture \
        "submit filter-alignments CRAM to BAM missing ref" \
        "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/submit_filter_alignments.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/scripts" \
            --threads 1 \
            --csv_infile "${in_cram}" \
            --dir_out "${dir_missing}" \
            --retain sc \
            --err_out "${dir_err}" \
            --nam_job "test_submit_filter_alignments_cram_to_bam_missing_ref"
then
    record_fail "submit_filter_alignments.sh CRAM to BAM without --ref_fa unexpectedly passed"
else
    record_pass "submit_filter_alignments.sh CRAM to BAM without --ref_fa fails"
fi

assert_pattern_found \
    "${log}" \
    "'--ref_fa' is required" \
    "submit CRAM to BAM missing-ref error mentions --ref_fa"


#  S. cerevisiae filtering should retain only canonical SC chromosomes
fil_out="${dir_out}/filter_sc_sp.sc.bam"
stat="${dir_out}/filter_sc_sp.sc.idxstats.txt"
log="${dir_log}/submit_filter_alignments_cram_to_bam_sc.log"

run_case_filter \
    sc \
    sc \
    "${log}" \
    arr_cmd_filter \
    "${nam_job_pfx}" \
    arr_arg_nil \
    "submit_filter_alignments.sh CRAM to BAM retain=sc sc"

assert_file_nonempty \
    "${fil_out}" \
    "submit CRAM to BAM retain=sc BAM output"
assert_file_nonempty \
    "${fil_out}.bai" \
    "submit CRAM to BAM retain=sc BAM index"
assert_filter_alignments_pg_header \
    "${fil_out}" \
    "" \
    filter_alignment_sc \
    sc \
    bam \
    "${dir_out}/filter_sc_sp.sc.header.txt" \
    "submit CRAM to BAM retain=sc BAM output"
assert_file_exists \
    "${dir_err}/test_submit_filter_alignments_cram_to_bam_sc.filter_sc_sp.stdout.txt" \
    "submit CRAM to BAM retain=sc stdout log"
assert_file_exists \
    "${dir_err}/test_submit_filter_alignments_cram_to_bam_sc.filter_sc_sp.stderr.txt" \
    "submit CRAM to BAM retain=sc stderr log"

if [[ -s "${fil_out}" ]]; then
    if \
        run_capture \
            "quickcheck submit filter-alignments CRAM to BAM sc" \
            "${dir_log}/submit_filter_alignments_cram_to_bam_sc_quickcheck.log" \
            run_samtools quickcheck "${fil_out}"
    then
        record_pass "submit CRAM to BAM retain=sc BAM passes samtools quickcheck"
    else
        record_fail "submit CRAM to BAM retain=sc BAM fails samtools quickcheck"
    fi

    run_capture \
        "idxstats submit filter-alignments CRAM to BAM sc" \
        "${stat}" \
        run_samtools idxstats "${fil_out}"

    assert_pattern_found \
        "${stat}" \
        $'^I\t100\t1\t0$' \
        "submit CRAM to BAM retain=sc output keeps chromosome I"
    assert_pattern_absent \
        "${stat}" \
        $'^Mito\t' \
        "submit CRAM to BAM retain=sc output omits Mito without --mito"
    assert_pattern_absent \
        "${stat}" \
        $'^SP_I\t' \
        "submit CRAM to BAM retain=sc output omits S. pombe chromosomes"
fi


#  S. pombe filtering should honor optional TG, MTR, and mito contigs
fil_out="${dir_out}/filter_sc_sp.sp.bam"
stat="${dir_out}/filter_sc_sp.sp.idxstats.txt"
log="${dir_log}/submit_filter_alignments_cram_to_bam_sp.log"

run_case_filter \
    sp \
    sp \
    "${log}" \
    arr_cmd_filter \
    "${nam_job_pfx}" \
    arr_arg_sp \
    "submit_filter_alignments.sh CRAM to BAM retain=sp sp"

assert_file_nonempty \
    "${fil_out}" \
    "submit CRAM to BAM retain=sp BAM output"
assert_file_nonempty \
    "${fil_out}.bai" \
    "submit CRAM to BAM retain=sp BAM index"
assert_filter_alignments_pg_header \
    "${fil_out}" \
    "" \
    filter_alignment_sp \
    sp \
    bam \
    "${dir_out}/filter_sc_sp.sp.header.txt" \
    "submit CRAM to BAM retain=sp BAM output"
assert_file_exists \
    "${dir_err}/test_submit_filter_alignments_cram_to_bam_sp.filter_sc_sp.stdout.txt" \
    "submit CRAM to BAM retain=sp stdout log"
assert_file_exists \
    "${dir_err}/test_submit_filter_alignments_cram_to_bam_sp.filter_sc_sp.stderr.txt" \
    "submit CRAM to BAM retain=sp stderr log"

if [[ -s "${fil_out}" ]]; then
    if \
        run_capture \
            "quickcheck submit filter-alignments CRAM to BAM sp" \
            "${dir_log}/submit_filter_alignments_cram_to_bam_sp_quickcheck.log" \
            run_samtools quickcheck "${fil_out}"
    then
        record_pass "submit CRAM to BAM retain=sp BAM passes samtools quickcheck"
    else
        record_fail "submit CRAM to BAM retain=sp BAM fails samtools quickcheck"
    fi

    run_capture \
        "idxstats submit filter-alignments CRAM to BAM sp" \
        "${stat}" \
        run_samtools idxstats "${fil_out}"

    assert_pattern_found \
        "${stat}" \
        $'^SP_I\t100\t1\t0$' \
        "submit CRAM to BAM retain=sp output keeps SP_I"
    assert_pattern_found \
        "${stat}" \
        $'^SP_II_TG\t100\t1\t0$' \
        "submit CRAM to BAM retain=sp output keeps SP_II_TG"
    assert_pattern_found \
        "${stat}" \
        $'^SP_MTR\t100\t1\t0$' \
        "submit CRAM to BAM retain=sp output keeps SP_MTR"
    assert_pattern_found \
        "${stat}" \
        $'^SP_Mito\t100\t1\t0$' \
        "submit CRAM to BAM retain=sp output keeps SP_Mito"
    assert_pattern_absent \
        "${stat}" \
        $'^I\t' \
        "submit CRAM to BAM retain=sp output omits S. cerevisiae chromosomes"
    assert_pattern_absent \
        "${stat}" \
        $'^chrUn\t' \
        "submit CRAM to BAM retain=sp output omits unrelated contigs"
fi

finish
