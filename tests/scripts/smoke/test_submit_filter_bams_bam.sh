#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_submit_filter_bams_bam.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="submit filter-bams BAM"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

rec_section "${TEST_NAME}"

#  Define fixture and output paths for BAM filtering
dir_fx="${ROOT_REPO}/tests/filter_bams/fixtures/sam"
in_sam="${dir_fx}/filter_sc_sp.sam"

tmp="${TEST_DIR_TMP}/submit_filter_bams_bam"
dir_in="${tmp}/input"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/filter_bams"
in_bam="${dir_in}/filter_sc_sp.bam"

rm -rf "${tmp}"
mkdir -p "${dir_in}" "${dir_out}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_exist "${in_sam}" || {
    finish
    exit $?
}


#  Build a deterministic BAM input from the committed SAM fixture
log="${dir_log}/submit_filter_bams_prepare_bam.log"
if ! \
    prepare_filter_bams_bam_fixture \
        "${in_sam}" "${in_bam}" "${log}" "submit filter-bams BAM fixture"
then
    finish
    exit $?
fi


#  Run submit_filter_bams.sh for one retain mode
function run_case_filter() {
    local retain="${1:-}"
    local nam_case="${2:-}"
    local log_lcl="${3:-}"

    shift 3

    # shellcheck disable=SC2154
    if \
        run_capture \
            "submit filter-bams ${nam_case}" "${log_lcl}" \
            "${TEST_BASH}" "${ROOT_REPO}/scripts/submit_filter_bams.sh" \
                --env_nam "${env_nam}" \
                --dir_scr "${ROOT_REPO}/scripts" \
                --threads 1 \
                --csv_infile "${in_bam}" \
                --dir_out "${dir_out}" \
                --retain "${retain}" \
                --err_out "${dir_err}" \
                --nam_job "test_submit_filter_bams_${nam_case}" \
                "$@"
    then
        rec_pass "submit_filter_bams.sh retain=${retain} ${nam_case} exits 0"
    else
        rec_fail \
            "submit_filter_bams.sh retain=${retain} ${nam_case} failed; see" \
            "$(rec_relpath "${log_lcl}")"
    fi
}


#  S. cerevisiae filtering should retain only canonical SC chromosomes
outfile="${dir_out}/filter_sc_sp.sc.bam"
idxstats="${dir_out}/filter_sc_sp.sc.idxstats.txt"
log="${dir_log}/submit_filter_bams_sc.log"

run_case_filter "sc" "sc" "${log}"

assert_file_nonempty "${outfile}" "submit retain=sc BAM output"
assert_file_nonempty "${outfile}.bai" "submit retain=sc BAM index"
assert_filter_bams_pg_header \
    "${outfile}" "" filter_bam_sc sc bam \
    "${dir_out}/filter_sc_sp.sc.header.txt" \
    "submit retain=sc BAM output"

if [[ -s "${outfile}" ]]; then
    run_capture \
        "idxstats submit filter-bams sc" "${idxstats}" \
        run_samtools idxstats "${outfile}"

    assert_grep_pattern "${idxstats}" $'^I\t100\t1\t0$' \
        "submit retain=sc output keeps chromosome I"
    assert_no_grep_pattern "${idxstats}" $'^Mito\t' \
        "submit retain=sc output omits Mito without --mito"
    assert_no_grep_pattern "${idxstats}" $'^SP_I\t' \
        "submit retain=sc output omits S. pombe chromosomes"
fi


#  S. pombe filtering should honor optional TG, MTR, and mito contigs
outfile="${dir_out}/filter_sc_sp.sp.bam"
idxstats="${dir_out}/filter_sc_sp.sp.idxstats.txt"
log="${dir_log}/submit_filter_bams_sp.log"

run_case_filter "sp" "sp" "${log}" --tg --mtr --mito

assert_file_nonempty "${outfile}" "submit retain=sp BAM output"
assert_file_nonempty "${outfile}.bai" "submit retain=sp BAM index"
assert_filter_bams_pg_header \
    "${outfile}" "" filter_bam_sp sp bam \
    "${dir_out}/filter_sc_sp.sp.header.txt" \
    "submit retain=sp BAM output"

if [[ -s "${outfile}" ]]; then
    run_capture \
        "idxstats submit filter-bams sp" "${idxstats}" \
        run_samtools idxstats "${outfile}"

    assert_grep_pattern "${idxstats}" $'^SP_I\t100\t1\t0$' \
        "submit retain=sp output keeps SP_I"
    assert_grep_pattern "${idxstats}" $'^SP_II_TG\t100\t1\t0$' \
        "submit retain=sp output keeps SP_II_TG"
    assert_grep_pattern "${idxstats}" $'^SP_MTR\t100\t1\t0$' \
        "submit retain=sp output keeps SP_MTR"
    assert_grep_pattern "${idxstats}" $'^SP_Mito\t100\t1\t0$' \
        "submit retain=sp output keeps SP_Mito"
    assert_no_grep_pattern "${idxstats}" $'^I\t' \
        "submit retain=sp output omits S. cerevisiae chromosomes"
    assert_no_grep_pattern "${idxstats}" $'^chrUn\t' \
        "submit retain=sp output omits unrelated contigs"
fi

finish
