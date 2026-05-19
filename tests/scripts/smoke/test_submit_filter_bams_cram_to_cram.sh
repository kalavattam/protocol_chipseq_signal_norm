#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_submit_filter_bams_cram_to_cram.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="submit filter-bams CRAM to CRAM"

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

tmp="${TEST_DIR_TMP}/submit_filter_bams_cram_to_cram"
dir_in="${tmp}/input"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/filter_bams"
in_cram="${dir_in}/filter_sc_sp.cram"

rm -rf "${tmp}"
mkdir -p "${dir_in}" "${dir_out}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_exist "${in_sam}" "${ref_fa}" "${ref_fai}" || {
    finish
    exit $?
}


#  Build deterministic CRAM input from committed SAM and reference fixtures
log="${dir_log}/submit_filter_bams_prepare_cram_to_cram.log"
if ! \
    prepare_filter_bams_cram_fixture \
        "${in_sam}" "${ref_fa}" "${in_cram}" "${log}" \
        "submit filter-bams CRAM fixture for CRAM output"
then
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
            "submit filter-bams CRAM to CRAM ${nam_case}" "${log_lcl}" \
            "${TEST_BASH}" "${ROOT_REPO}/scripts/submit_filter_bams.sh" \
                --env_nam "${env_nam}" \
                --dir_scr "${ROOT_REPO}/scripts" \
                --threads 1 \
                --csv_infile "${in_cram}" \
                --dir_out "${dir_out}" \
                --out_ext cram \
                --retain "${retain}" \
                --ref_fa "${ref_fa}" \
                --err_out "${dir_err}" \
                --nam_job "test_submit_filter_bams_cram_to_cram_${nam_case}" \
                "$@"
    then
        rec_pass "submit_filter_bams.sh CRAM to CRAM retain=${retain} ${nam_case} exits 0"
    else
        rec_fail \
            "submit_filter_bams.sh CRAM to CRAM retain=${retain} ${nam_case} failed; see" \
            "$(rec_relpath "${log_lcl}")"
    fi
}


#  S. cerevisiae filtering should retain only canonical SC chromosomes
outfile="${dir_out}/filter_sc_sp.sc.cram"
log="${dir_log}/submit_filter_bams_cram_to_cram_sc.log"

run_case_filter "sc" "sc" "${log}"

assert_file_nonempty "${outfile}" "submit CRAM-to-CRAM retain=sc output"
assert_cram_index "${outfile}" "submit CRAM-to-CRAM retain=sc CRAI index"
assert_filter_bams_pg_header \
    "${outfile}" "${ref_fa}" filter_bam_sc sc cram \
    "${dir_out}/filter_sc_sp.sc.header.txt" \
    "submit CRAM-to-CRAM retain=sc output"
assert_file_exists \
    "${dir_err}/test_submit_filter_bams_cram_to_cram_sc.filter_sc_sp.stdout.txt" \
    "submit CRAM-to-CRAM retain=sc stdout log"
assert_file_exists \
    "${dir_err}/test_submit_filter_bams_cram_to_cram_sc.filter_sc_sp.stderr.txt" \
    "submit CRAM-to-CRAM retain=sc stderr log"

if [[ -s "${outfile}" ]]; then
    if \
        run_capture \
            "quickcheck submit filter-bams CRAM to CRAM sc" \
            "${dir_log}/submit_filter_bams_cram_to_cram_sc_quickcheck.log" \
            run_samtools quickcheck "${outfile}"
    then
        rec_pass "submit CRAM-to-CRAM retain=sc passes samtools quickcheck"
    else
        rec_fail "submit CRAM-to-CRAM retain=sc fails samtools quickcheck"
    fi

    assert_cram_count "${outfile}" "${ref_fa}" I 1 \
        "submit CRAM-to-CRAM retain=sc has one read on chromosome I" \
        "${dir_out}/filter_sc_sp.sc.I.count.txt"
    assert_cram_count "${outfile}" "${ref_fa}" Mito 0 \
        "submit CRAM-to-CRAM retain=sc omits Mito without --mito" \
        "${dir_out}/filter_sc_sp.sc.Mito.count.txt"
    assert_cram_count "${outfile}" "${ref_fa}" SP_I 0 \
        "submit CRAM-to-CRAM retain=sc omits S. pombe chromosomes" \
        "${dir_out}/filter_sc_sp.sc.SP_I.count.txt"
fi


#  S. pombe filtering should honor optional TG, MTR, and mito contigs
outfile="${dir_out}/filter_sc_sp.sp.cram"
log="${dir_log}/submit_filter_bams_cram_to_cram_sp.log"

run_case_filter "sp" "sp" "${log}" --tg --mtr --mito

assert_file_nonempty "${outfile}" "submit CRAM-to-CRAM retain=sp output"
assert_cram_index "${outfile}" "submit CRAM-to-CRAM retain=sp CRAI index"
assert_filter_bams_pg_header \
    "${outfile}" "${ref_fa}" filter_bam_sp sp cram \
    "${dir_out}/filter_sc_sp.sp.header.txt" \
    "submit CRAM-to-CRAM retain=sp output"
assert_file_exists \
    "${dir_err}/test_submit_filter_bams_cram_to_cram_sp.filter_sc_sp.stdout.txt" \
    "submit CRAM-to-CRAM retain=sp stdout log"
assert_file_exists \
    "${dir_err}/test_submit_filter_bams_cram_to_cram_sp.filter_sc_sp.stderr.txt" \
    "submit CRAM-to-CRAM retain=sp stderr log"

if [[ -s "${outfile}" ]]; then
    if \
        run_capture \
            "quickcheck submit filter-bams CRAM to CRAM sp" \
            "${dir_log}/submit_filter_bams_cram_to_cram_sp_quickcheck.log" \
            run_samtools quickcheck "${outfile}"
    then
        rec_pass "submit CRAM-to-CRAM retain=sp passes samtools quickcheck"
    else
        rec_fail "submit CRAM-to-CRAM retain=sp fails samtools quickcheck"
    fi

    assert_cram_count "${outfile}" "${ref_fa}" SP_I 1 \
        "submit CRAM-to-CRAM retain=sp keeps SP_I" \
        "${dir_out}/filter_sc_sp.sp.SP_I.count.txt"
    assert_cram_count "${outfile}" "${ref_fa}" SP_II_TG 1 \
        "submit CRAM-to-CRAM retain=sp keeps SP_II_TG" \
        "${dir_out}/filter_sc_sp.sp.SP_II_TG.count.txt"
    assert_cram_count "${outfile}" "${ref_fa}" SP_MTR 1 \
        "submit CRAM-to-CRAM retain=sp keeps SP_MTR" \
        "${dir_out}/filter_sc_sp.sp.SP_MTR.count.txt"
    assert_cram_count "${outfile}" "${ref_fa}" SP_Mito 1 \
        "submit CRAM-to-CRAM retain=sp keeps SP_Mito" \
        "${dir_out}/filter_sc_sp.sp.SP_Mito.count.txt"
    assert_cram_count "${outfile}" "${ref_fa}" I 0 \
        "submit CRAM-to-CRAM retain=sp omits S. cerevisiae chromosomes" \
        "${dir_out}/filter_sc_sp.sp.I.count.txt"
    assert_cram_count "${outfile}" "${ref_fa}" chrUn 0 \
        "submit CRAM-to-CRAM retain=sp omits unrelated contigs" \
        "${dir_out}/filter_sc_sp.sp.chrUn.count.txt"
fi

finish
