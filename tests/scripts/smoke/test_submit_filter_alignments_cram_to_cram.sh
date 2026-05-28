#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_submit_filter_alignments_cram_to_cram.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="submit filter-alignments CRAM to CRAM"

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

tmp="${TEST_DIR_TMP}/submit_filter_alignments_cram_to_cram"
dir_in="${tmp}/input"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/filter_alignments"
in_cram="${dir_in}/filter_sc_sp.cram"

rm -rf "${tmp}"
mkdir -p "${dir_in}" "${dir_out}" "${dir_err}" "${dir_log}"

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
            --out_ext cram
            --ref_fa "${ref_fa}"
            --err_out "${dir_err}"
    )
    nam_job_pfx="test_submit_filter_alignments_cram_to_cram"
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
log="${dir_log}/submit_filter_alignments_prepare_cram_to_cram.log"
if ! \
    build_filter_alignments_cram_fixture \
        "${in_sam}" \
        "${ref_fa}" \
        "${in_cram}" \
        "${log}" \
        "submit filter-alignments CRAM fixture for CRAM output"
then
    finish
    exit $?
fi


#  S. cerevisiae filtering should retain only canonical SC chromosomes
fil_out="${dir_out}/filter_sc_sp.sc.cram"
log="${dir_log}/submit_filter_alignments_cram_to_cram_sc.log"

run_case_filter \
    sc \
    sc \
    "${log}" \
    arr_cmd_filter \
    "${nam_job_pfx}" \
    arr_arg_nil \
    "submit_filter_alignments.sh CRAM to CRAM retain=sc sc"

assert_file_nonempty \
    "${fil_out}" \
    "submit CRAM-to-CRAM retain=sc output"
assert_cram_index \
    "${fil_out}" \
    "submit CRAM-to-CRAM retain=sc CRAI index"
assert_filter_alignments_pg_header \
    "${fil_out}" \
    "${ref_fa}" \
    filter_alignment_sc \
    sc \
    cram \
    "${dir_out}/filter_sc_sp.sc.header.txt" \
    "submit CRAM-to-CRAM retain=sc output"
assert_file_exists \
    "${dir_err}/test_submit_filter_alignments_cram_to_cram_sc.filter_sc_sp.stdout.txt" \
    "submit CRAM-to-CRAM retain=sc stdout log"
assert_file_exists \
    "${dir_err}/test_submit_filter_alignments_cram_to_cram_sc.filter_sc_sp.stderr.txt" \
    "submit CRAM-to-CRAM retain=sc stderr log"

if [[ -s "${fil_out}" ]]; then
    if \
        run_capture \
            "quickcheck submit filter-alignments CRAM to CRAM sc" \
            "${dir_log}/submit_filter_alignments_cram_to_cram_sc_quickcheck.log" \
            run_samtools quickcheck "${fil_out}"
    then
        record_pass "submit CRAM-to-CRAM retain=sc passes samtools quickcheck"
    else
        record_fail "submit CRAM-to-CRAM retain=sc fails samtools quickcheck"
    fi

    assert_cram_count \
        "${fil_out}" \
        "${ref_fa}" \
        1 \
        "${dir_out}/filter_sc_sp.sc.I.count.txt" \
        arr_rgn_i \
        "submit CRAM-to-CRAM retain=sc has one read on chromosome I"
    assert_cram_count \
        "${fil_out}" \
        "${ref_fa}" \
        0 \
        "${dir_out}/filter_sc_sp.sc.Mito.count.txt" \
        arr_rgn_mito \
        "submit CRAM-to-CRAM retain=sc omits Mito without --mito"
    assert_cram_count \
        "${fil_out}" \
        "${ref_fa}" \
        0 \
        "${dir_out}/filter_sc_sp.sc.SP_I.count.txt" \
        arr_rgn_sp_i \
        "submit CRAM-to-CRAM retain=sc omits S. pombe chromosomes"
fi


#  S. pombe filtering should honor optional TG, MTR, and mito contigs
fil_out="${dir_out}/filter_sc_sp.sp.cram"
log="${dir_log}/submit_filter_alignments_cram_to_cram_sp.log"

run_case_filter \
    sp \
    sp \
    "${log}" \
    arr_cmd_filter \
    "${nam_job_pfx}" \
    arr_arg_sp \
    "submit_filter_alignments.sh CRAM to CRAM retain=sp sp"

assert_file_nonempty \
    "${fil_out}" \
    "submit CRAM-to-CRAM retain=sp output"
assert_cram_index \
    "${fil_out}" \
    "submit CRAM-to-CRAM retain=sp CRAI index"
assert_filter_alignments_pg_header \
    "${fil_out}" \
    "${ref_fa}" \
    filter_alignment_sp \
    sp \
    cram \
    "${dir_out}/filter_sc_sp.sp.header.txt" \
    "submit CRAM-to-CRAM retain=sp output"
assert_file_exists \
    "${dir_err}/test_submit_filter_alignments_cram_to_cram_sp.filter_sc_sp.stdout.txt" \
    "submit CRAM-to-CRAM retain=sp stdout log"
assert_file_exists \
    "${dir_err}/test_submit_filter_alignments_cram_to_cram_sp.filter_sc_sp.stderr.txt" \
    "submit CRAM-to-CRAM retain=sp stderr log"

if [[ -s "${fil_out}" ]]; then
    if \
        run_capture \
            "quickcheck submit filter-alignments CRAM to CRAM sp" \
            "${dir_log}/submit_filter_alignments_cram_to_cram_sp_quickcheck.log" \
            run_samtools quickcheck "${fil_out}"
    then
        record_pass "submit CRAM-to-CRAM retain=sp passes samtools quickcheck"
    else
        record_fail "submit CRAM-to-CRAM retain=sp fails samtools quickcheck"
    fi

    assert_cram_count \
        "${fil_out}" \
        "${ref_fa}" \
        1 \
        "${dir_out}/filter_sc_sp.sp.SP_I.count.txt" \
        arr_rgn_sp_i \
        "submit CRAM-to-CRAM retain=sp keeps SP_I"
    assert_cram_count \
        "${fil_out}" \
        "${ref_fa}" \
        1 \
        "${dir_out}/filter_sc_sp.sp.SP_II_TG.count.txt" \
        arr_rgn_sp_ii_tg \
        "submit CRAM-to-CRAM retain=sp keeps SP_II_TG"
    assert_cram_count \
        "${fil_out}" \
        "${ref_fa}" \
        1 \
        "${dir_out}/filter_sc_sp.sp.SP_MTR.count.txt" \
        arr_rgn_sp_mtr \
        "submit CRAM-to-CRAM retain=sp keeps SP_MTR"
    assert_cram_count \
        "${fil_out}" \
        "${ref_fa}" \
        1 \
        "${dir_out}/filter_sc_sp.sp.SP_Mito.count.txt" \
        arr_rgn_sp_mito \
        "submit CRAM-to-CRAM retain=sp keeps SP_Mito"
    assert_cram_count \
        "${fil_out}" \
        "${ref_fa}" \
        0 \
        "${dir_out}/filter_sc_sp.sp.I.count.txt" \
        arr_rgn_i \
        "submit CRAM-to-CRAM retain=sp omits S. cerevisiae chromosomes"
    assert_cram_count \
        "${fil_out}" \
        "${ref_fa}" \
        0 \
        "${dir_out}/filter_sc_sp.sp.chrUn.count.txt" \
        arr_rgn_chrun \
        "submit CRAM-to-CRAM retain=sp omits unrelated contigs"
fi

finish
