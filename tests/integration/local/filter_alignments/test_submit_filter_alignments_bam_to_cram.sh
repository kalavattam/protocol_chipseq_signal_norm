#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_submit_filter_alignments_bam_to_cram.sh
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

TEST_NAME="submit filter-alignments BAM to CRAM"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


dir_fx="${ROOT_REPO}/tests/fixtures/filter_alignments"
in_sam="${dir_fx}/sam/filter_sc_sp.sam"
ref_fa="${dir_fx}/reference/filter_sc_sp.fa"
ref_fai="${ref_fa}.fai"

tmp="${TEST_DIR_TMP}/submit_filter_alignments_bam_to_cram"
dir_in="${tmp}/input"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_missing="${tmp}/missing_ref"
dir_log="${TEST_DIR_LOG}/filter_alignments"
in_bam="${dir_in}/filter_sc_sp.bam"

log_prepare="${dir_log}/submit_filter_alignments_prepare_bam_to_cram.log"
log_missing_ref="${dir_log}/submit_filter_alignments_bam_to_cram_missing_ref.log"

fil_out_sc="${dir_out}/filter_sc_sp.sc.cram"
hdr_sc="${dir_out}/filter_sc_sp.sc.header.txt"

log_sc="${dir_log}/submit_filter_alignments_bam_to_cram_sc.log"
log_quick_sc="${dir_log}/submit_filter_alignments_bam_to_cram_sc_quickcheck.log"
log_pfx_sc="${dir_err}/test_submit_filter_alignments_bam_to_cram_sc"
log_sub_out_sc="${log_pfx_sc}.filter_sc_sp.stdout.txt"
log_sub_err_sc="${log_pfx_sc}.filter_sc_sp.stderr.txt"

cnt_sc_i="${dir_out}/filter_sc_sp.sc.I.count.txt"
cnt_sc_mito="${dir_out}/filter_sc_sp.sc.Mito.count.txt"
cnt_sc_sp_i="${dir_out}/filter_sc_sp.sc.SP_I.count.txt"

fil_out_sp="${dir_out}/filter_sc_sp.sp.cram"
hdr_sp="${dir_out}/filter_sc_sp.sp.header.txt"

log_sp="${dir_log}/submit_filter_alignments_bam_to_cram_sp.log"
log_quick_sp="${dir_log}/submit_filter_alignments_bam_to_cram_sp_quickcheck.log"
log_pfx_sp="${dir_err}/test_submit_filter_alignments_bam_to_cram_sp"
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
    arr_rgn_i=( I )
    arr_rgn_mito=( Mito )
    arr_rgn_sp_i=( SP_I )
    arr_rgn_sp_ii_tg=( SP_II_TG )
    arr_rgn_sp_mtr=( SP_MTR )
    arr_rgn_sp_mito=( SP_Mito )
    arr_rgn_chrun=( chrUn )
}

print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${dir_in}" "${dir_out}" "${dir_err}" "${dir_missing}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

# shellcheck disable=SC2034
{
    arr_cmd_filter=(
        "${TEST_BASH}" "${ROOT_REPO}/bin/submit_filter_alignments.sh"
            --env_nam "${env_nam}"
            --dir_scr "${ROOT_REPO}/bin"
            --threads 1
            --csv_fil_in "${in_bam}"
            --dir_out "${dir_out}"
            --out_ext cram
            --ref_fa "${ref_fa}"
            --dir_eo "${dir_err}"
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
if ! \
    build_filter_alignments_fixture_bam \
        "${in_sam}" \
        "${in_bam}" \
        "${log_prepare}" \
        "submit filter-alignments BAM fixture for CRAM output"
then
    finish
    exit $?
fi


#  CRAM output without --ref_fa should fail clearly
if \
    run_capture \
        "submit filter-alignments BAM to CRAM missing ref" \
        "${log_missing_ref}" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/submit_filter_alignments.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/bin" \
            --threads 1 \
            --csv_fil_in "${in_bam}" \
            --dir_out "${dir_missing}" \
            --out_ext cram \
            --retain sc \
            --dir_eo "${dir_err}" \
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
    "${log_missing_ref}" \
    "'--ref_fa' is required" \
    "submit BAM-to-CRAM missing-ref error mentions --ref_fa"


#  S. cerevisiae filtering should retain only canonical SC chromosomes
run_case_filter \
    sc \
    sc \
    "${log_sc}" \
    arr_cmd_filter \
    "${nam_job_pfx}" \
    arr_arg_nil \
    "submit_filter_alignments.sh BAM to CRAM retain=sc sc"

assert_file_nonempty \
    "${fil_out_sc}" \
    "submit retain=sc CRAM output"

assert_cram_index \
    "${fil_out_sc}" \
    "submit retain=sc CRAI index"

assert_filter_alignments_pg_header \
    "${fil_out_sc}" \
    "${ref_fa}" \
    filter_alignment_sc \
    sc \
    cram \
    "${hdr_sc}" \
    "submit retain=sc CRAM output"

assert_file_exists \
    "${log_sub_out_sc}" \
    "submit retain=sc CRAM stdout log"

assert_file_exists \
    "${log_sub_err_sc}" \
    "submit retain=sc CRAM stderr log"

if [[ -s "${fil_out_sc}" ]]; then
    if \
        run_capture \
            "quickcheck submit filter-alignments BAM to CRAM sc" \
            "${log_quick_sc}" \
            run_samtools quickcheck "${fil_out_sc}"
    then
        record_pass "submit retain=sc CRAM passes samtools quickcheck"
    else
        record_fail "submit retain=sc CRAM fails samtools quickcheck"
    fi

    assert_cram_count \
        "${fil_out_sc}" \
        "${ref_fa}" \
        1 \
        "${cnt_sc_i}" \
        arr_rgn_i \
        "submit retain=sc CRAM has one read on chromosome I"

    assert_cram_count \
        "${fil_out_sc}" \
        "${ref_fa}" \
        0 \
        "${cnt_sc_mito}" \
        arr_rgn_mito \
        "submit retain=sc CRAM omits Mito without --mito"

    assert_cram_count \
        "${fil_out_sc}" \
        "${ref_fa}" \
        0 \
        "${cnt_sc_sp_i}" \
        arr_rgn_sp_i \
        "submit retain=sc CRAM omits S. pombe chromosomes"
fi


#  S. pombe filtering should honor optional TG, MTR, and mito contigs
run_case_filter \
    sp \
    sp \
    "${log_sp}" \
    arr_cmd_filter \
    "${nam_job_pfx}" \
    arr_arg_sp \
    "submit_filter_alignments.sh BAM to CRAM retain=sp sp"

assert_file_nonempty \
    "${fil_out_sp}" \
    "submit retain=sp CRAM output"

assert_cram_index \
    "${fil_out_sp}" \
    "submit retain=sp CRAI index"

assert_filter_alignments_pg_header \
    "${fil_out_sp}" \
    "${ref_fa}" \
    filter_alignment_sp \
    sp \
    cram \
    "${hdr_sp}" \
    "submit retain=sp CRAM output"

assert_file_exists \
    "${log_sub_out_sp}" \
    "submit retain=sp CRAM stdout log"

assert_file_exists \
    "${log_sub_err_sp}" \
    "submit retain=sp CRAM stderr log"

if [[ -s "${fil_out_sp}" ]]; then
    if \
        run_capture \
            "quickcheck submit filter-alignments BAM to CRAM sp" \
            "${log_quick_sp}" \
            run_samtools quickcheck "${fil_out_sp}"
    then
        record_pass "submit retain=sp CRAM passes samtools quickcheck"
    else
        record_fail "submit retain=sp CRAM fails samtools quickcheck"
    fi

    assert_cram_count \
        "${fil_out_sp}" \
        "${ref_fa}" \
        1 \
        "${cnt_sp_sp_i}" \
        arr_rgn_sp_i \
        "submit retain=sp CRAM keeps SP_I"

    assert_cram_count \
        "${fil_out_sp}" \
        "${ref_fa}" \
        1 \
        "${cnt_sp_sp_ii_tg}" \
        arr_rgn_sp_ii_tg \
        "submit retain=sp CRAM keeps SP_II_TG"

    assert_cram_count \
        "${fil_out_sp}" \
        "${ref_fa}" \
        1 \
        "${cnt_sp_sp_mtr}" \
        arr_rgn_sp_mtr \
        "submit retain=sp CRAM keeps SP_MTR"

    assert_cram_count \
        "${fil_out_sp}" \
        "${ref_fa}" \
        1 \
        "${cnt_sp_sp_mito}" \
        arr_rgn_sp_mito \
        "submit retain=sp CRAM keeps SP_Mito"

    assert_cram_count \
        "${fil_out_sp}" \
        "${ref_fa}" \
        0 \
        "${cnt_sp_i}" \
        arr_rgn_i \
        "submit retain=sp CRAM omits S. cerevisiae chromosomes"

    assert_cram_count \
        "${fil_out_sp}" \
        "${ref_fa}" \
        0 \
        "${cnt_sp_chrun}" \
        arr_rgn_chrun \
        "submit retain=sp CRAM omits unrelated contigs"
fi

finish
