#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_submit_filter_bams_cram_input_bam.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="submit filter-bams CRAM input BAM output"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

rec_section "${TEST_NAME}"

dir_fx="${ROOT_REPO}/tests/filter_bams/fixtures"
in_sam="${dir_fx}/sam/filter_sc_sp.sam"
ref_fa="${dir_fx}/reference/filter_sc_sp.fa"
ref_fai="${ref_fa}.fai"

tmp="${TEST_DIR_TMP}/submit_filter_bams_cram_input_bam"
dir_in="${tmp}/input"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_missing="${tmp}/missing_ref"
dir_log="${TEST_DIR_LOG}/filter_bams"
in_cram="${dir_in}/filter_sc_sp.cram"

rm -rf "${tmp}"
mkdir -p "${dir_in}" "${dir_out}" "${dir_err}" "${dir_missing}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_exist "${in_sam}" "${ref_fa}" "${ref_fai}" || {
    finish
    exit $?
}


#  Build deterministic CRAM input from committed SAM and reference fixtures
log="${dir_log}/submit_filter_bams_prepare_cram.log"
if \
    run_capture \
        "prepare submit filter-bams CRAM fixture" "${log}" \
        run_samtools view -C -T "${ref_fa}" -o "${in_cram}" "${in_sam}" \
    && run_capture \
        "index submit filter-bams CRAM fixture" "${log}.index" \
        run_samtools index "${in_cram}"
then
    rec_pass "filter CRAM input fixture is prepared"
else
    rec_fail "failed to prepare filter CRAM fixture; see $(rec_relpath "${log}")"
    finish
    exit $?
fi


function run_case_filter() {
    local retain="${1:-}"
    local nam_case="${2:-}"
    local log_lcl="${3:-}"

    shift 3

    # shellcheck disable=SC2154
    if \
        run_capture \
            "submit filter-bams CRAM ${nam_case}" "${log_lcl}" \
            "${TEST_BASH}" "${ROOT_REPO}/scripts/submit_filter_bams.sh" \
                --env_nam "${env_nam}" \
                --dir_scr "${ROOT_REPO}/scripts" \
                --threads 1 \
                --csv_infile "${in_cram}" \
                --dir_out "${dir_out}" \
                --retain "${retain}" \
                --ref_fa "${ref_fa}" \
                --err_out "${dir_err}" \
                --nam_job "test_submit_filter_bams_cram_${nam_case}" \
                "$@"
    then
        rec_pass "submit_filter_bams.sh CRAM retain=${retain} ${nam_case} exits 0"
    else
        rec_fail \
            "submit_filter_bams.sh CRAM retain=${retain} ${nam_case} failed; see" \
            "$(rec_relpath "${log_lcl}")"
    fi
}


#  CRAM input without --ref_fa should fail clearly
log="${dir_log}/submit_filter_bams_cram_missing_ref.log"
if \
    run_capture \
        "submit filter-bams CRAM missing ref" "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/submit_filter_bams.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/scripts" \
            --threads 1 \
            --csv_infile "${in_cram}" \
            --dir_out "${dir_missing}" \
            --retain sc \
            --err_out "${dir_err}" \
            --nam_job "test_submit_filter_bams_cram_missing_ref"
then
    rec_fail "submit_filter_bams.sh CRAM without --ref_fa unexpectedly passed"
else
    rec_pass "submit_filter_bams.sh CRAM without --ref_fa fails"
fi

assert_grep_pattern "${log}" "'--ref_fa' is required" \
    "submit CRAM missing-ref error mentions --ref_fa"


#  S. cerevisiae filtering should retain only canonical SC chromosomes
outfile="${dir_out}/filter_sc_sp.sc.bam"
idxstats="${dir_out}/filter_sc_sp.sc.idxstats.txt"
log="${dir_log}/submit_filter_bams_cram_sc.log"

run_case_filter "sc" "sc" "${log}"

assert_file_nonempty "${outfile}" "submit CRAM retain=sc BAM output"
assert_file_nonempty "${outfile}.bai" "submit CRAM retain=sc BAM index"
assert_file_exists \
    "${dir_err}/test_submit_filter_bams_cram_sc.filter_sc_sp.stdout.txt" \
    "submit CRAM retain=sc stdout log"
assert_file_exists \
    "${dir_err}/test_submit_filter_bams_cram_sc.filter_sc_sp.stderr.txt" \
    "submit CRAM retain=sc stderr log"

if [[ -s "${outfile}" ]]; then
    if \
        run_capture \
            "quickcheck submit filter-bams CRAM sc" \
            "${dir_log}/submit_filter_bams_cram_sc_quickcheck.log" \
            run_samtools quickcheck "${outfile}"
    then
        rec_pass "submit CRAM retain=sc BAM passes samtools quickcheck"
    else
        rec_fail "submit CRAM retain=sc BAM fails samtools quickcheck"
    fi

    run_capture \
        "idxstats submit filter-bams CRAM sc" "${idxstats}" \
        run_samtools idxstats "${outfile}"

    assert_grep_pattern "${idxstats}" $'^I\t100\t1\t0$' \
        "submit CRAM retain=sc output keeps chromosome I"
    assert_no_grep_pattern "${idxstats}" $'^Mito\t' \
        "submit CRAM retain=sc output omits Mito without --mito"
    assert_no_grep_pattern "${idxstats}" $'^SP_I\t' \
        "submit CRAM retain=sc output omits S. pombe chromosomes"
fi


#  S. pombe filtering should honor optional TG, MTR, and mito contigs
outfile="${dir_out}/filter_sc_sp.sp.bam"
idxstats="${dir_out}/filter_sc_sp.sp.idxstats.txt"
log="${dir_log}/submit_filter_bams_cram_sp.log"

run_case_filter "sp" "sp" "${log}" --tg --mtr --mito

assert_file_nonempty "${outfile}" "submit CRAM retain=sp BAM output"
assert_file_nonempty "${outfile}.bai" "submit CRAM retain=sp BAM index"
assert_file_exists \
    "${dir_err}/test_submit_filter_bams_cram_sp.filter_sc_sp.stdout.txt" \
    "submit CRAM retain=sp stdout log"
assert_file_exists \
    "${dir_err}/test_submit_filter_bams_cram_sp.filter_sc_sp.stderr.txt" \
    "submit CRAM retain=sp stderr log"

if [[ -s "${outfile}" ]]; then
    if \
        run_capture \
            "quickcheck submit filter-bams CRAM sp" \
            "${dir_log}/submit_filter_bams_cram_sp_quickcheck.log" \
            run_samtools quickcheck "${outfile}"
    then
        rec_pass "submit CRAM retain=sp BAM passes samtools quickcheck"
    else
        rec_fail "submit CRAM retain=sp BAM fails samtools quickcheck"
    fi

    run_capture \
        "idxstats submit filter-bams CRAM sp" "${idxstats}" \
        run_samtools idxstats "${outfile}"

    assert_grep_pattern "${idxstats}" $'^SP_I\t100\t1\t0$' \
        "submit CRAM retain=sp output keeps SP_I"
    assert_grep_pattern "${idxstats}" $'^SP_II_TG\t100\t1\t0$' \
        "submit CRAM retain=sp output keeps SP_II_TG"
    assert_grep_pattern "${idxstats}" $'^SP_MTR\t100\t1\t0$' \
        "submit CRAM retain=sp output keeps SP_MTR"
    assert_grep_pattern "${idxstats}" $'^SP_Mito\t100\t1\t0$' \
        "submit CRAM retain=sp output keeps SP_Mito"
    assert_no_grep_pattern "${idxstats}" $'^I\t' \
        "submit CRAM retain=sp output omits S. cerevisiae chromosomes"
    assert_no_grep_pattern "${idxstats}" $'^chrUn\t' \
        "submit CRAM retain=sp output omits unrelated contigs"
fi

finish
