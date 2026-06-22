#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_filter_alignments_cram_to_bam.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="execute filter-alignments CRAM to BAM"

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

tmp="${TEST_DIR_TMP}/execute_filter_alignments_cram_to_bam"
dir_in="${tmp}/input"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/filter_alignments"
in_cram="${dir_in}/filter_sc_sp.cram"

log_prepare="${dir_log}/execute_filter_alignments_prepare_cram_to_bam.log"

fil_out_sc="${dir_out}/filter_sc_sp.sc.bam"
bai_sc="${fil_out_sc}.bai"
stat_sc="${dir_out}/filter_sc_sp.sc.idxstats.txt"

log_sc="${dir_log}/execute_filter_alignments_cram_to_bam_sc.log"
log_quick_sc="${dir_log}/execute_filter_alignments_cram_to_bam_sc_quickcheck.log"
log_pfx_sc="${dir_err}/test_execute_filter_alignments_cram_to_bam_sc"
log_ser_out_sc="${log_pfx_sc}_ser.stdout.txt"
log_ser_err_sc="${log_pfx_sc}_ser.stderr.txt"
log_sub_out_sc="${log_pfx_sc}.filter_sc_sp.stdout.txt"
log_sub_err_sc="${log_pfx_sc}.filter_sc_sp.stderr.txt"

fil_out_sp="${dir_out}/filter_sc_sp.sp.bam"
bai_sp="${fil_out_sp}.bai"
stat_sp="${dir_out}/filter_sc_sp.sp.idxstats.txt"

log_sp="${dir_log}/execute_filter_alignments_cram_to_bam_sp.log"
log_quick_sp="${dir_log}/execute_filter_alignments_cram_to_bam_sp_quickcheck.log"
log_pfx_sp="${dir_err}/test_execute_filter_alignments_cram_to_bam_sp"
log_ser_out_sp="${log_pfx_sp}_ser.stdout.txt"
log_ser_err_sp="${log_pfx_sp}_ser.stderr.txt"
log_sub_out_sp="${log_pfx_sp}.filter_sc_sp.stdout.txt"
log_sub_err_sp="${log_pfx_sp}.filter_sc_sp.stderr.txt"

# shellcheck disable=SC2034
{
    arr_cmd_filter=(
        "${TEST_BASH}" "${ROOT_REPO}/scripts/execute_filter_alignments.sh"
            --threads 1
            --csv_infile "${in_cram}"
            --dir_out "${dir_out}"
            --ref_fa "${ref_fa}"
            --err_out "${dir_err}"
            --max_job 1
    )
    nam_job_pfx="test_execute_filter_alignments_cram_to_bam"
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


#  Build deterministic CRAM input from committed SAM and reference fixtures
if ! \
    build_filter_alignments_fixture_cram \
        "${in_sam}" \
        "${ref_fa}" \
        "${in_cram}" \
        "${log_prepare}" \
        "execute filter-alignments CRAM-to-BAM fixture"
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
    "execute_filter_alignments.sh CRAM to BAM retain=sc sc"

assert_file_nonempty \
    "${fil_out_sc}" \
    "execute CRAM to BAM retain=sc BAM output"

assert_file_nonempty \
    "${bai_sc}" \
    "execute CRAM to BAM retain=sc BAM index"

assert_file_exists \
    "${log_ser_out_sc}" \
    "execute CRAM to BAM retain=sc serial stdout log"

assert_file_exists \
    "${log_ser_err_sc}" \
    "execute CRAM to BAM retain=sc serial stderr log"

assert_file_exists \
    "${log_sub_out_sc}" \
    "execute CRAM to BAM retain=sc submit stdout log"

assert_file_exists \
    "${log_sub_err_sc}" \
    "execute CRAM to BAM retain=sc submit stderr log"

if [[ -s "${fil_out_sc}" ]]; then
    if \
        run_capture \
            "quickcheck execute filter-alignments CRAM to BAM sc" \
            "${log_quick_sc}" \
            run_samtools quickcheck "${fil_out_sc}"
    then
        record_pass "execute CRAM to BAM retain=sc BAM passes samtools quickcheck"
    else
        record_fail "execute CRAM to BAM retain=sc BAM fails samtools quickcheck"
    fi

    run_capture \
        "idxstats execute filter-alignments CRAM to BAM sc" \
        "${stat_sc}" \
        run_samtools idxstats "${fil_out_sc}"

    assert_pattern_found \
        "${stat_sc}" \
        $'^I\t100\t1\t0$' \
        "execute CRAM to BAM retain=sc output keeps chromosome I"

    assert_pattern_absent \
        "${stat_sc}" \
        $'^Mito\t' \
        "execute CRAM to BAM retain=sc output omits Mito without --mito"

    assert_pattern_absent \
        "${stat_sc}" \
        $'^SP_I\t' \
        "execute CRAM to BAM retain=sc output omits S. pombe chromosomes"
fi


#  S. pombe filtering should honor optional TG, MTR, and mito contigs
run_case_filter \
    sp \
    sp \
    "${log_sp}" \
    arr_cmd_filter \
    "${nam_job_pfx}" \
    arr_arg_sp \
    "execute_filter_alignments.sh CRAM to BAM retain=sp sp"

assert_file_nonempty \
    "${fil_out_sp}" \
    "execute CRAM to BAM retain=sp BAM output"

assert_file_nonempty \
    "${bai_sp}" \
    "execute CRAM to BAM retain=sp BAM index"

assert_file_exists \
    "${log_ser_out_sp}" \
    "execute CRAM to BAM retain=sp serial stdout log"

assert_file_exists \
    "${log_ser_err_sp}" \
    "execute CRAM to BAM retain=sp serial stderr log"

assert_file_exists \
    "${log_sub_out_sp}" \
    "execute CRAM to BAM retain=sp submit stdout log"

assert_file_exists \
    "${log_sub_err_sp}" \
    "execute CRAM to BAM retain=sp submit stderr log"

if [[ -s "${fil_out_sp}" ]]; then
    if \
        run_capture \
            "quickcheck execute filter-alignments CRAM to BAM sp" \
            "${log_quick_sp}" \
            run_samtools quickcheck "${fil_out_sp}"
    then
        record_pass "execute CRAM to BAM retain=sp BAM passes samtools quickcheck"
    else
        record_fail "execute CRAM to BAM retain=sp BAM fails samtools quickcheck"
    fi

    run_capture \
        "idxstats execute filter-alignments CRAM to BAM sp" \
        "${stat_sp}" \
        run_samtools idxstats "${fil_out_sp}"

    assert_pattern_found \
        "${stat_sp}" \
        $'^SP_I\t100\t1\t0$' \
        "execute CRAM to BAM retain=sp output keeps SP_I"

    assert_pattern_found \
        "${stat_sp}" \
        $'^SP_II_TG\t100\t1\t0$' \
        "execute CRAM to BAM retain=sp output keeps SP_II_TG"

    assert_pattern_found \
        "${stat_sp}" \
        $'^SP_MTR\t100\t1\t0$' \
        "execute CRAM to BAM retain=sp output keeps SP_MTR"

    assert_pattern_found \
        "${stat_sp}" \
        $'^SP_Mito\t100\t1\t0$' \
        "execute CRAM to BAM retain=sp output keeps SP_Mito"

    assert_pattern_absent \
        "${stat_sp}" \
        $'^I\t' \
        "execute CRAM to BAM retain=sp output omits S. cerevisiae chromosomes"

    assert_pattern_absent \
        "${stat_sp}" \
        $'^chrUn\t' \
        "execute CRAM to BAM retain=sp output omits unrelated contigs"
fi

finish
