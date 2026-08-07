#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_submit_filter_alignments_bam_to_bam.sh
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

TEST_NAME="submit filter-alignments BAM to BAM"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


#  Define fixture and output paths for BAM input filtering
dir_fx="${ROOT_REPO}/tests/fixtures/filter_alignments/sam"
in_sam="${dir_fx}/filter_sc_sp.sam"

tmp="${TEST_DIR_TMP}/submit_filter_alignments_bam_to_bam"
dir_in="${tmp}/input"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/filter_alignments"
in_bam="${dir_in}/filter_sc_sp.bam"

log_prep="${dir_log}/submit_filter_alignments_prepare_bam_to_bam.log"

fil_out_sc="${dir_out}/filter_sc_sp.sc.bam"
bai_sc="${fil_out_sc}.bai"
stat_sc="${dir_out}/filter_sc_sp.sc.idxstats.txt"
hdr_sc="${dir_out}/filter_sc_sp.sc.header.txt"
log_sc="${dir_log}/submit_filter_alignments_bam_to_bam_sc.log"

fil_out_sp="${dir_out}/filter_sc_sp.sp.bam"
bai_sp="${fil_out_sp}.bai"
stat_sp="${dir_out}/filter_sc_sp.sp.idxstats.txt"
hdr_sp="${dir_out}/filter_sc_sp.sp.header.txt"
log_sp="${dir_log}/submit_filter_alignments_bam_to_bam_sp.log"


print_section "${TEST_NAME}"

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
        "${TEST_BASH}" "${ROOT_REPO}/bin/submit_filter_alignments.sh"
            --env_nam "${env_nam}"
            --dir_scr "${ROOT_REPO}/bin"
            --threads 1
            --csv_fil_in "${in_bam}"
            --dir_out "${dir_out}"
            --dir_eo "${dir_err}"
    )
    nam_job_pfx="test_submit_filter_alignments_bam_to_bam"
    arr_arg_nil=()
    arr_arg_sp=( --tg --mtr --mito )
}

require_files_nonempty \
    "${in_sam}" || {
    finish
    exit $?
}

#  Build a deterministic BAM input from the committed SAM fixture
if ! \
    build_filter_alignments_fixture_bam \
        "${in_sam}" \
        "${in_bam}" \
        "${log_prep}" \
        "submit filter-alignments BAM-to-BAM fixture"
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
    "submit_filter_alignments.sh BAM to BAM retain=sc sc"

assert_file_nonempty \
    "${fil_out_sc}" \
    "submit BAM to BAM retain=sc BAM output"

assert_file_nonempty \
    "${bai_sc}" \
    "submit BAM to BAM retain=sc BAM index"

assert_filter_alignments_pg_header \
    "${fil_out_sc}" \
    "" \
    filter_alignment_sc \
    sc \
    bam \
    "${hdr_sc}" \
    "submit BAM to BAM retain=sc BAM output"

if [[ -s "${fil_out_sc}" ]]; then
    run_capture \
        "idxstats submit filter-alignments BAM to BAM sc" \
        "${stat_sc}" \
        run_samtools idxstats "${fil_out_sc}"

    assert_pattern_found \
        "${stat_sc}" \
        $'^I\t100\t1\t0$' \
        "submit BAM to BAM retain=sc output keeps chromosome I"

    assert_pattern_absent \
        "${stat_sc}" \
        $'^Mito\t' \
        "submit BAM to BAM retain=sc output omits Mito without --mito"

    assert_pattern_absent \
        "${stat_sc}" \
        $'^SP_I\t' \
        "submit BAM to BAM retain=sc output omits S. pombe chromosomes"
fi


#  S. pombe filtering should honor optional TG, MTR, and mito contigs
run_case_filter \
    sp \
    sp \
    "${log_sp}" \
    arr_cmd_filter \
    "${nam_job_pfx}" \
    arr_arg_sp \
    "submit_filter_alignments.sh BAM to BAM retain=sp sp"

assert_file_nonempty \
    "${fil_out_sp}" \
    "submit BAM to BAM retain=sp BAM output"

assert_file_nonempty \
    "${bai_sp}" \
    "submit BAM to BAM retain=sp BAM index"

assert_filter_alignments_pg_header \
    "${fil_out_sp}" \
    "" \
    filter_alignment_sp \
    sp \
    bam \
    "${hdr_sp}" \
    "submit BAM to BAM retain=sp BAM output"

if [[ -s "${fil_out_sp}" ]]; then
    run_capture \
        "idxstats submit filter-alignments BAM to BAM sp" \
        "${stat_sp}" \
        run_samtools idxstats "${fil_out_sp}"

    assert_pattern_found \
        "${stat_sp}" \
        $'^SP_I\t100\t1\t0$' \
        "submit BAM to BAM retain=sp output keeps SP_I"

    assert_pattern_found \
        "${stat_sp}" \
        $'^SP_II_TG\t100\t1\t0$' \
        "submit BAM to BAM retain=sp output keeps SP_II_TG"

    assert_pattern_found \
        "${stat_sp}" \
        $'^SP_MTR\t100\t1\t0$' \
        "submit BAM to BAM retain=sp output keeps SP_MTR"

    assert_pattern_found \
        "${stat_sp}" \
        $'^SP_Mito\t100\t1\t0$' \
        "submit BAM to BAM retain=sp output keeps SP_Mito"

    assert_pattern_absent \
        "${stat_sp}" \
        $'^I\t' \
        "submit BAM to BAM retain=sp output omits S. cerevisiae chromosomes"

    assert_pattern_absent \
        "${stat_sp}" \
        $'^chrUn\t' \
        "submit BAM to BAM retain=sp output omits unrelated contigs"
fi

finish
