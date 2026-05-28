#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_submit_filter_alignments_bam_to_cram.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="submit filter-alignments BAM to CRAM"

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

# shellcheck disable=SC2034
{
    arr_rgn_i=( I )
    arr_rgn_mito=( Mito )
    arr_rgn_sp_i=( SP_I )
    arr_rgn_sp_ii_tg=( SP_II_TG )
    arr_rgn_sp_mtr=( SP_MTR )
    arr_rgn_sp_mito=( SP_Mito )
    arr_rgn_chrun=( chrUn )
}

tmp="${TEST_DIR_TMP}/submit_filter_alignments_bam_to_cram"
dir_in="${tmp}/input"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_missing="${tmp}/missing_ref"
dir_log="${TEST_DIR_LOG}/filter_alignments"
in_bam="${dir_in}/filter_sc_sp.bam"

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
            --csv_infile "${in_bam}"
            --dir_out "${dir_out}"
            --out_ext cram
            --ref_fa "${ref_fa}"
            --err_out "${dir_err}"
    )
    nam_job_pfx="test_submit_filter_alignments_bam_to_cram"
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


#  Build deterministic BAM input from committed SAM fixture
log="${dir_log}/submit_filter_alignments_prepare_bam_to_cram.log"
if ! \
    build_filter_alignments_bam_fixture \
        "${in_sam}" \
        "${in_bam}" \
        "${log}" \
        "submit filter-alignments BAM fixture for CRAM output"
then
    finish
    exit $?
fi


#  CRAM output without --ref_fa should fail clearly
log="${dir_log}/submit_filter_alignments_bam_to_cram_missing_ref.log"
if \
    run_capture \
        "submit filter-alignments BAM to CRAM missing ref" \
        "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/submit_filter_alignments.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/scripts" \
            --threads 1 \
            --csv_infile "${in_bam}" \
            --dir_out "${dir_missing}" \
            --out_ext cram \
            --retain sc \
            --err_out "${dir_err}" \
            --nam_job "test_submit_filter_alignments_bam_to_cram_missing_ref"
then
    record_fail \
        "submit_filter_alignments.sh --out_ext cram without --ref_fa" \
        "unexpectedly passed"
else
    record_pass \
        "submit_filter_alignments.sh --out_ext cram without --ref_fa fails"
fi

assert_pattern_found \
    "${log}" \
    "'--ref_fa' is required" \
    "submit BAM-to-CRAM missing-ref error mentions --ref_fa"


#  S. cerevisiae filtering should retain only canonical SC chromosomes
fil_out="${dir_out}/filter_sc_sp.sc.cram"
log="${dir_log}/submit_filter_alignments_bam_to_cram_sc.log"

run_case_filter \
    sc \
    sc \
    "${log}" \
    arr_cmd_filter \
    "${nam_job_pfx}" \
    arr_arg_nil \
    "submit_filter_alignments.sh BAM to CRAM retain=sc sc"

assert_file_nonempty \
    "${fil_out}" \
    "submit retain=sc CRAM output"
assert_cram_index \
    "${fil_out}" \
    "submit retain=sc CRAI index"
assert_filter_alignments_pg_header \
    "${fil_out}" \
    "${ref_fa}" \
    filter_alignment_sc \
    sc \
    cram \
    "${dir_out}/filter_sc_sp.sc.header.txt" \
    "submit retain=sc CRAM output"
assert_file_exists \
    "${dir_err}/test_submit_filter_alignments_bam_to_cram_sc.filter_sc_sp.stdout.txt" \
    "submit retain=sc CRAM stdout log"
assert_file_exists \
    "${dir_err}/test_submit_filter_alignments_bam_to_cram_sc.filter_sc_sp.stderr.txt" \
    "submit retain=sc CRAM stderr log"

if [[ -s "${fil_out}" ]]; then
    if \
        run_capture \
            "quickcheck submit filter-alignments BAM to CRAM sc" \
            "${dir_log}/submit_filter_alignments_bam_to_cram_sc_quickcheck.log" \
            run_samtools quickcheck "${fil_out}"
    then
        record_pass "submit retain=sc CRAM passes samtools quickcheck"
    else
        record_fail "submit retain=sc CRAM fails samtools quickcheck"
    fi

    assert_cram_count \
        "${fil_out}" \
        "${ref_fa}" \
        1 \
        "${dir_out}/filter_sc_sp.sc.I.count.txt" \
        arr_rgn_i \
        "submit retain=sc CRAM has one read on chromosome I"
    assert_cram_count \
        "${fil_out}" \
        "${ref_fa}" \
        0 \
        "${dir_out}/filter_sc_sp.sc.Mito.count.txt" \
        arr_rgn_mito \
        "submit retain=sc CRAM omits Mito without --mito"
    assert_cram_count \
        "${fil_out}" \
        "${ref_fa}" \
        0 \
        "${dir_out}/filter_sc_sp.sc.SP_I.count.txt" \
        arr_rgn_sp_i \
        "submit retain=sc CRAM omits S. pombe chromosomes"
fi


#  S. pombe filtering should honor optional TG, MTR, and mito contigs
fil_out="${dir_out}/filter_sc_sp.sp.cram"
log="${dir_log}/submit_filter_alignments_bam_to_cram_sp.log"

run_case_filter \
    sp \
    sp \
    "${log}" \
    arr_cmd_filter \
    "${nam_job_pfx}" \
    arr_arg_sp \
    "submit_filter_alignments.sh BAM to CRAM retain=sp sp"

assert_file_nonempty \
    "${fil_out}" \
    "submit retain=sp CRAM output"
assert_cram_index \
    "${fil_out}" \
    "submit retain=sp CRAI index"
assert_filter_alignments_pg_header \
    "${fil_out}" \
    "${ref_fa}" \
    filter_alignment_sp \
    sp \
    cram \
    "${dir_out}/filter_sc_sp.sp.header.txt" \
    "submit retain=sp CRAM output"
assert_file_exists \
    "${dir_err}/test_submit_filter_alignments_bam_to_cram_sp.filter_sc_sp.stdout.txt" \
    "submit retain=sp CRAM stdout log"
assert_file_exists \
    "${dir_err}/test_submit_filter_alignments_bam_to_cram_sp.filter_sc_sp.stderr.txt" \
    "submit retain=sp CRAM stderr log"

if [[ -s "${fil_out}" ]]; then
    if \
        run_capture \
            "quickcheck submit filter-alignments BAM to CRAM sp" \
            "${dir_log}/submit_filter_alignments_bam_to_cram_sp_quickcheck.log" \
            run_samtools quickcheck "${fil_out}"
    then
        record_pass "submit retain=sp CRAM passes samtools quickcheck"
    else
        record_fail "submit retain=sp CRAM fails samtools quickcheck"
    fi

    assert_cram_count \
        "${fil_out}" \
        "${ref_fa}" \
        1 \
        "${dir_out}/filter_sc_sp.sp.SP_I.count.txt" \
        arr_rgn_sp_i \
        "submit retain=sp CRAM keeps SP_I"
    assert_cram_count \
        "${fil_out}" \
        "${ref_fa}" \
        1 \
        "${dir_out}/filter_sc_sp.sp.SP_II_TG.count.txt" \
        arr_rgn_sp_ii_tg \
        "submit retain=sp CRAM keeps SP_II_TG"
    assert_cram_count \
        "${fil_out}" \
        "${ref_fa}" \
        1 \
        "${dir_out}/filter_sc_sp.sp.SP_MTR.count.txt" \
        arr_rgn_sp_mtr \
        "submit retain=sp CRAM keeps SP_MTR"
    assert_cram_count \
        "${fil_out}" \
        "${ref_fa}" \
        1 \
        "${dir_out}/filter_sc_sp.sp.SP_Mito.count.txt" \
        arr_rgn_sp_mito \
        "submit retain=sp CRAM keeps SP_Mito"
    assert_cram_count \
        "${fil_out}" \
        "${ref_fa}" \
        0 \
        "${dir_out}/filter_sc_sp.sp.I.count.txt" \
        arr_rgn_i \
        "submit retain=sp CRAM omits S. cerevisiae chromosomes"
    assert_cram_count \
        "${fil_out}" \
        "${ref_fa}" \
        0 \
        "${dir_out}/filter_sc_sp.sp.chrUn.count.txt" \
        arr_rgn_chrun \
        "submit retain=sp CRAM omits unrelated contigs"
fi

finish
