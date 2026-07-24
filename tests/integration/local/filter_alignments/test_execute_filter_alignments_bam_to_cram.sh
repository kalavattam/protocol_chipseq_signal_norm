#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_filter_alignments_bam_to_cram.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="execute filter-alignments BAM to CRAM"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


dir_fx="${ROOT_REPO}/tests/fixtures/filter_alignments"
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

log_prepare="${dir_log}/execute_filter_alignments_prepare_bam_to_cram.log"

fil_out_sc="${dir_out}/filter_sc_sp.sc.cram"
log_sc="${dir_log}/execute_filter_alignments_bam_to_cram_sc.log"
log_quick_sc="${dir_log}/execute_filter_alignments_bam_to_cram_sc_quickcheck.log"
log_pfx_sc="${dir_err}/test_execute_filter_alignments_bam_to_cram_sc"
log_ser_out_sc="${log_pfx_sc}_ser.stdout.txt"
log_ser_err_sc="${log_pfx_sc}_ser.stderr.txt"
log_sub_out_sc="${log_pfx_sc}.filter_sc_sp.stdout.txt"
log_sub_err_sc="${log_pfx_sc}.filter_sc_sp.stderr.txt"

cnt_sc_i="${dir_out}/filter_sc_sp.sc.I.count.txt"
cnt_sc_mito="${dir_out}/filter_sc_sp.sc.Mito.count.txt"
cnt_sc_sp_i="${dir_out}/filter_sc_sp.sc.SP_I.count.txt"

fil_out_sp="${dir_out}/filter_sc_sp.sp.cram"
log_sp="${dir_log}/execute_filter_alignments_bam_to_cram_sp.log"
log_quick_sp="${dir_log}/execute_filter_alignments_bam_to_cram_sp_quickcheck.log"
log_pfx_sp="${dir_err}/test_execute_filter_alignments_bam_to_cram_sp"
log_ser_out_sp="${log_pfx_sp}_ser.stdout.txt"
log_ser_err_sp="${log_pfx_sp}_ser.stderr.txt"
log_sub_out_sp="${log_pfx_sp}.filter_sc_sp.stdout.txt"
log_sub_err_sp="${log_pfx_sp}.filter_sc_sp.stderr.txt"

cnt_sp_i="${dir_out}/filter_sc_sp.sp.I.count.txt"
cnt_sp_sp_i="${dir_out}/filter_sc_sp.sp.SP_I.count.txt"
cnt_sp_sp_ii_tg="${dir_out}/filter_sc_sp.sp.SP_II_TG.count.txt"
cnt_sp_sp_mtr="${dir_out}/filter_sc_sp.sp.SP_MTR.count.txt"
cnt_sp_sp_mito="${dir_out}/filter_sc_sp.sp.SP_Mito.count.txt"
cnt_sp_chrun="${dir_out}/filter_sc_sp.sp.chrUn.count.txt"

# shellcheck disable=SC2034
{
    arr_cmd_filter=(
        "${TEST_BASH}" "${ROOT_REPO}/bin/execute_filter_alignments.sh"
            --threads 1
            --csv_fil_in "${in_bam}"
            --dir_out "${dir_out}"
            --out_ext cram
            --ref_fa "${ref_fa}"
            --dir_eo "${dir_err}"
            --max_job 1
    )
    nam_job_pfx="test_execute_filter_alignments_bam_to_cram"
    arr_arg_nil=()
    arr_arg_sp=( --tg --mtr --mito )
}

print_section "${TEST_NAME}"

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
if ! \
    build_filter_alignments_fixture_bam \
        "${in_sam}" \
        "${in_bam}" \
        "${log_prepare}" \
        "execute filter-alignments BAM fixture for CRAM output"
then
    finish
    exit $?
fi


#  S. cerevisiae filtering should retain only canonical SC chromosomes
run_case_filter \
    sc \
    sc \
    "${log_sc}" \
    arr_cmd_filter \
    "${nam_job_pfx}" \
    arr_arg_nil \
    "execute_filter_alignments.sh BAM to CRAM retain=sc sc"

assert_file_nonempty \
    "${fil_out_sc}" \
    "execute retain=sc CRAM output"

assert_cram_index \
    "${fil_out_sc}" \
    "execute retain=sc CRAI index"

assert_file_exists \
    "${log_ser_out_sc}" \
    "execute retain=sc serial stdout log"

assert_file_exists \
    "${log_ser_err_sc}" \
    "execute retain=sc serial stderr log"

assert_file_exists \
    "${log_sub_out_sc}" \
    "execute retain=sc submit stdout log"

assert_file_exists \
    "${log_sub_err_sc}" \
    "execute retain=sc submit stderr log"

if [[ -s "${fil_out_sc}" ]]; then
    if \
        run_capture \
            "quickcheck execute filter-alignments BAM to CRAM sc" \
            "${log_quick_sc}" \
            run_samtools quickcheck "${fil_out_sc}"
    then
        record_pass "execute retain=sc CRAM passes samtools quickcheck"
    else
        record_fail "execute retain=sc CRAM fails samtools quickcheck"
    fi

    assert_cram_count \
        "${fil_out_sc}" \
        "${ref_fa}" \
        1 \
        "${cnt_sc_i}" \
        arr_rgn_i \
        "execute retain=sc CRAM has one read on chromosome I"

    assert_cram_count \
        "${fil_out_sc}" \
        "${ref_fa}" \
        0 \
        "${cnt_sc_mito}" \
        arr_rgn_mito \
        "execute retain=sc CRAM omits Mito without --mito"

    assert_cram_count \
        "${fil_out_sc}" \
        "${ref_fa}" \
        0 \
        "${cnt_sc_sp_i}" \
        arr_rgn_sp_i \
        "execute retain=sc CRAM omits S. pombe chromosomes"
fi


#  S. pombe filtering should honor optional TG, MTR, and mito contigs
run_case_filter \
    sp \
    sp \
    "${log_sp}" \
    arr_cmd_filter \
    "${nam_job_pfx}" \
    arr_arg_sp \
    "execute_filter_alignments.sh BAM to CRAM retain=sp sp"

assert_file_nonempty \
    "${fil_out_sp}" \
    "execute retain=sp CRAM output"

assert_cram_index \
    "${fil_out_sp}" \
    "execute retain=sp CRAI index"

assert_file_exists \
    "${log_ser_out_sp}" \
    "execute retain=sp serial stdout log"

assert_file_exists \
    "${log_ser_err_sp}" \
    "execute retain=sp serial stderr log"

assert_file_exists \
    "${log_sub_out_sp}" \
    "execute retain=sp submit stdout log"

assert_file_exists \
    "${log_sub_err_sp}" \
    "execute retain=sp submit stderr log"

if [[ -s "${fil_out_sp}" ]]; then
    if \
        run_capture \
            "quickcheck execute filter-alignments BAM to CRAM sp" \
            "${log_quick_sp}" \
            run_samtools quickcheck "${fil_out_sp}"
    then
        record_pass "execute retain=sp CRAM passes samtools quickcheck"
    else
        record_fail "execute retain=sp CRAM fails samtools quickcheck"
    fi

    assert_cram_count \
        "${fil_out_sp}" \
        "${ref_fa}" \
        1 \
        "${cnt_sp_sp_i}" \
        arr_rgn_sp_i \
        "execute retain=sp CRAM keeps SP_I"

    assert_cram_count \
        "${fil_out_sp}" \
        "${ref_fa}" \
        1 \
        "${cnt_sp_sp_ii_tg}" \
        arr_rgn_sp_ii_tg \
        "execute retain=sp CRAM keeps SP_II_TG"

    assert_cram_count \
        "${fil_out_sp}" \
        "${ref_fa}" \
        1 \
        "${cnt_sp_sp_mtr}" \
        arr_rgn_sp_mtr \
        "execute retain=sp CRAM keeps SP_MTR"

    assert_cram_count \
        "${fil_out_sp}" \
        "${ref_fa}" \
        1 \
        "${cnt_sp_sp_mito}" \
        arr_rgn_sp_mito \
        "execute retain=sp CRAM keeps SP_Mito"

    assert_cram_count \
        "${fil_out_sp}" \
        "${ref_fa}" \
        0 \
        "${cnt_sp_i}" \
        arr_rgn_i \
        "execute retain=sp CRAM omits S. cerevisiae chromosomes"

    assert_cram_count \
        "${fil_out_sp}" \
        "${ref_fa}" \
        0 \
        "${cnt_sp_chrun}" \
        arr_rgn_chrun \
        "execute retain=sp CRAM omits unrelated contigs"
fi

finish
