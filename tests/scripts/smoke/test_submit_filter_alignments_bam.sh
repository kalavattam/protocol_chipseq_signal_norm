#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_submit_filter_alignments_bam.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="submit filter-alignments BAM"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

print_section "${TEST_NAME}"

#  Define fixture and output paths for BAM input filtering
dir_fx="${ROOT_REPO}/tests/filter_alignments/fixtures/sam"
in_sam="${dir_fx}/filter_sc_sp.sam"

tmp="${TEST_DIR_TMP}/submit_filter_alignments_bam"
dir_in="${tmp}/input"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/filter_alignments"
in_bam="${dir_in}/filter_sc_sp.bam"

rm -rf "${tmp}"
mkdir -p "${dir_in}" "${dir_out}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

#  Configure submit_filter_alignments.sh for retain-mode cases
# shellcheck disable=SC2034
{
    arr_cmd_filter=(
        "${TEST_BASH}" "${ROOT_REPO}/scripts/submit_filter_alignments.sh"
            --env_nam "${env_nam}"
            --dir_scr "${ROOT_REPO}/scripts"
            --threads 1
            --csv_infile "${in_bam}"
            --dir_out "${dir_out}"
            --err_out "${dir_err}"
    )
    nam_job_pfx="test_submit_filter_alignments"
    arr_arg_nil=()
    arr_arg_sp=( --tg --mtr --mito )
}

require_files_nonempty \
    "${in_sam}" || {
    finish
    exit $?
}

#  Build a deterministic BAM input from the committed SAM fixture
log="${dir_log}/submit_filter_alignments_prepare_bam.log"
if ! \
    build_filter_alignments_bam_fixture \
        "${in_sam}" \
        "${in_bam}" \
        "${log}" \
        "submit filter-alignments BAM fixture"
then
    finish
    exit $?
fi


#  S. cerevisiae filtering should retain only canonical SC chromosomes
fil_out="${dir_out}/filter_sc_sp.sc.bam"
stat="${dir_out}/filter_sc_sp.sc.idxstats.txt"
log="${dir_log}/submit_filter_alignments_sc.log"

run_case_filter \
    sc \
    sc \
    "${log}" \
    arr_cmd_filter \
    "${nam_job_pfx}" \
    arr_arg_nil \
    "submit_filter_alignments.sh retain=sc sc"

assert_file_nonempty \
    "${fil_out}" \
    "submit retain=sc BAM output"
assert_file_nonempty \
    "${fil_out}.bai" \
    "submit retain=sc BAM index"
assert_filter_alignments_pg_header \
    "${fil_out}" \
    "" \
    filter_alignment_sc \
    sc \
    bam \
    "${dir_out}/filter_sc_sp.sc.header.txt" \
    "submit retain=sc BAM output"

if [[ -s "${fil_out}" ]]; then
    run_capture \
        "idxstats submit filter-alignments sc" \
        "${stat}" \
        run_samtools idxstats "${fil_out}"

    assert_pattern_found \
        "${stat}" \
        $'^I\t100\t1\t0$' \
        "submit retain=sc output keeps chromosome I"
    assert_pattern_absent \
        "${stat}" \
        $'^Mito\t' \
        "submit retain=sc output omits Mito without --mito"
    assert_pattern_absent \
        "${stat}" \
        $'^SP_I\t' \
        "submit retain=sc output omits S. pombe chromosomes"
fi


#  S. pombe filtering should honor optional TG, MTR, and mito contigs
fil_out="${dir_out}/filter_sc_sp.sp.bam"
stat="${dir_out}/filter_sc_sp.sp.idxstats.txt"
log="${dir_log}/submit_filter_alignments_sp.log"

run_case_filter \
    sp \
    sp \
    "${log}" \
    arr_cmd_filter \
    "${nam_job_pfx}" \
    arr_arg_sp \
    "submit_filter_alignments.sh retain=sp sp"

assert_file_nonempty \
    "${fil_out}" \
    "submit retain=sp BAM output"
assert_file_nonempty \
    "${fil_out}.bai" \
    "submit retain=sp BAM index"
assert_filter_alignments_pg_header \
    "${fil_out}" \
    "" \
    filter_alignment_sp \
    sp \
    bam \
    "${dir_out}/filter_sc_sp.sp.header.txt" \
    "submit retain=sp BAM output"

if [[ -s "${fil_out}" ]]; then
    run_capture \
        "idxstats submit filter-alignments sp" \
        "${stat}" \
        run_samtools idxstats "${fil_out}"

    assert_pattern_found \
        "${stat}" \
        $'^SP_I\t100\t1\t0$' \
        "submit retain=sp output keeps SP_I"
    assert_pattern_found \
        "${stat}" \
        $'^SP_II_TG\t100\t1\t0$' \
        "submit retain=sp output keeps SP_II_TG"
    assert_pattern_found \
        "${stat}" \
        $'^SP_MTR\t100\t1\t0$' \
        "submit retain=sp output keeps SP_MTR"
    assert_pattern_found \
        "${stat}" \
        $'^SP_Mito\t100\t1\t0$' \
        "submit retain=sp output keeps SP_Mito"
    assert_pattern_absent \
        "${stat}" \
        $'^I\t' \
        "submit retain=sp output omits S. cerevisiae chromosomes"
    assert_pattern_absent \
        "${stat}" \
        $'^chrUn\t' \
        "submit retain=sp output omits unrelated contigs"
fi

finish
