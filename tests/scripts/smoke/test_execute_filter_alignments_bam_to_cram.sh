#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_filter_alignments_bam_to_cram.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5) were used in development.
#
# Distributed under the MIT license.


TEST_NAME="execute filter-alignments BAM to CRAM"

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

tmp="${TEST_DIR_TMP}/execute_filter_alignments_bam_to_cram"
dir_in="${tmp}/input"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/filter_alignments"
in_bam="${dir_in}/filter_sc_sp.bam"

# shellcheck disable=SC2034
{
    arr_cmd_filter=(
        "${TEST_BASH}" "${ROOT_REPO}/scripts/execute_filter_alignments.sh"
            --threads 1
            --csv_infile "${in_bam}"
            --dir_out "${dir_out}"
            --out_ext cram
            --ref_fa "${ref_fa}"
            --err_out "${dir_err}"
            --max_job 1
    )
    nam_job_pfx="test_execute_filter_alignments_bam_to_cram"
    arr_arg_nil=()
    arr_arg_sp=( --tg --mtr --mito )
}

rm -rf "${tmp}"
mkdir -p "${dir_in}" "${dir_out}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty \
    "${in_sam}" \
    "${ref_fa}" \
    "${ref_fai}" || {
    finish
    exit $?
}


#  Build deterministic BAM input from committed SAM fixture
log="${dir_log}/execute_filter_alignments_prepare_bam_to_cram.log"
if ! \
    build_filter_alignments_fixture_bam \
        "${in_sam}" \
        "${in_bam}" \
        "${log}" \
        "execute filter-alignments BAM fixture for CRAM output"
then
    finish
    exit $?
fi


#  S. cerevisiae filtering should retain only canonical SC chromosomes
fil_out="${dir_out}/filter_sc_sp.sc.cram"
log="${dir_log}/execute_filter_alignments_bam_to_cram_sc.log"

run_case_filter \
    sc \
    sc \
    "${log}" \
    arr_cmd_filter \
    "${nam_job_pfx}" \
    arr_arg_nil \
    "execute_filter_alignments.sh BAM to CRAM retain=sc sc"

assert_file_nonempty \
    "${fil_out}" \
    "execute retain=sc CRAM output"
assert_cram_index \
    "${fil_out}" \
    "execute retain=sc CRAI index"
assert_file_exists \
    "${dir_err}/test_execute_filter_alignments_bam_to_cram_sc_ser.stdout.txt" \
    "execute retain=sc serial stdout log"
assert_file_exists \
    "${dir_err}/test_execute_filter_alignments_bam_to_cram_sc_ser.stderr.txt" \
    "execute retain=sc serial stderr log"
assert_file_exists \
    "${dir_err}/test_execute_filter_alignments_bam_to_cram_sc.filter_sc_sp.stdout.txt" \
    "execute retain=sc submit stdout log"
assert_file_exists \
    "${dir_err}/test_execute_filter_alignments_bam_to_cram_sc.filter_sc_sp.stderr.txt" \
    "execute retain=sc submit stderr log"

if [[ -s "${fil_out}" ]]; then
    if \
        run_capture \
            "quickcheck execute filter-alignments BAM to CRAM sc" \
            "${dir_log}/execute_filter_alignments_bam_to_cram_sc_quickcheck.log" \
            run_samtools quickcheck "${fil_out}"
    then
        record_pass "execute retain=sc CRAM passes samtools quickcheck"
    else
        record_fail "execute retain=sc CRAM fails samtools quickcheck"
    fi

    assert_cram_count \
        "${fil_out}" \
        "${ref_fa}" \
        1 \
        "${dir_out}/filter_sc_sp.sc.I.count.txt" \
        arr_rgn_i \
        "execute retain=sc CRAM has one read on chromosome I"
    assert_cram_count \
        "${fil_out}" \
        "${ref_fa}" \
        0 \
        "${dir_out}/filter_sc_sp.sc.Mito.count.txt" \
        arr_rgn_mito \
        "execute retain=sc CRAM omits Mito without --mito"
    assert_cram_count \
        "${fil_out}" \
        "${ref_fa}" \
        0 \
        "${dir_out}/filter_sc_sp.sc.SP_I.count.txt" \
        arr_rgn_sp_i \
        "execute retain=sc CRAM omits S. pombe chromosomes"
fi


#  S. pombe filtering should honor optional TG, MTR, and mito contigs
fil_out="${dir_out}/filter_sc_sp.sp.cram"
log="${dir_log}/execute_filter_alignments_bam_to_cram_sp.log"

run_case_filter \
    sp \
    sp \
    "${log}" \
    arr_cmd_filter \
    "${nam_job_pfx}" \
    arr_arg_sp \
    "execute_filter_alignments.sh BAM to CRAM retain=sp sp"

assert_file_nonempty \
    "${fil_out}" \
    "execute retain=sp CRAM output"
assert_cram_index \
    "${fil_out}" \
    "execute retain=sp CRAI index"
assert_file_exists \
    "${dir_err}/test_execute_filter_alignments_bam_to_cram_sp_ser.stdout.txt" \
    "execute retain=sp serial stdout log"
assert_file_exists \
    "${dir_err}/test_execute_filter_alignments_bam_to_cram_sp_ser.stderr.txt" \
    "execute retain=sp serial stderr log"
assert_file_exists \
    "${dir_err}/test_execute_filter_alignments_bam_to_cram_sp.filter_sc_sp.stdout.txt" \
    "execute retain=sp submit stdout log"
assert_file_exists \
    "${dir_err}/test_execute_filter_alignments_bam_to_cram_sp.filter_sc_sp.stderr.txt" \
    "execute retain=sp submit stderr log"

if [[ -s "${fil_out}" ]]; then
    if \
        run_capture \
            "quickcheck execute filter-alignments BAM to CRAM sp" \
            "${dir_log}/execute_filter_alignments_bam_to_cram_sp_quickcheck.log" \
            run_samtools quickcheck "${fil_out}"
    then
        record_pass "execute retain=sp CRAM passes samtools quickcheck"
    else
        record_fail "execute retain=sp CRAM fails samtools quickcheck"
    fi

    assert_cram_count \
        "${fil_out}" \
        "${ref_fa}" \
        1 \
        "${dir_out}/filter_sc_sp.sp.SP_I.count.txt" \
        arr_rgn_sp_i \
        "execute retain=sp CRAM keeps SP_I"
    assert_cram_count \
        "${fil_out}" \
        "${ref_fa}" \
        1 \
        "${dir_out}/filter_sc_sp.sp.SP_II_TG.count.txt" \
        arr_rgn_sp_ii_tg \
        "execute retain=sp CRAM keeps SP_II_TG"
    assert_cram_count \
        "${fil_out}" \
        "${ref_fa}" \
        1 \
        "${dir_out}/filter_sc_sp.sp.SP_MTR.count.txt" \
        arr_rgn_sp_mtr \
        "execute retain=sp CRAM keeps SP_MTR"
    assert_cram_count \
        "${fil_out}" \
        "${ref_fa}" \
        1 \
        "${dir_out}/filter_sc_sp.sp.SP_Mito.count.txt" \
        arr_rgn_sp_mito \
        "execute retain=sp CRAM keeps SP_Mito"
    assert_cram_count \
        "${fil_out}" \
        "${ref_fa}" \
        0 \
        "${dir_out}/filter_sc_sp.sp.I.count.txt" \
        arr_rgn_i \
        "execute retain=sp CRAM omits S. cerevisiae chromosomes"
    assert_cram_count \
        "${fil_out}" \
        "${ref_fa}" \
        0 \
        "${dir_out}/filter_sc_sp.sp.chrUn.count.txt" \
        arr_rgn_chrun \
        "execute retain=sp CRAM omits unrelated contigs"
fi

finish
