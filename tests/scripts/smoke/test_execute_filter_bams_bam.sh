#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_filter_bams_bam.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="execute filter-bams BAM"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

rec_section "${TEST_NAME}"

#  Define fixture and output paths for the execute-to-submit BAM filter path
dir_fx="${ROOT_REPO}/tests/filter_bams/fixtures/sam"
in_sam="${dir_fx}/filter_sc_sp.sam"

tmp="${TEST_DIR_TMP}/execute_filter_bams_bam"
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


#  Run samtools from the active shell or the resolved project environment
function run_samtools() {
    if check_cmd_exists samtools; then
        samtools "$@"
    else
        conda run -n "${env_nam}" samtools "$@"
    fi
}


#  Build a deterministic BAM input from the committed SAM fixture
log="${dir_log}/execute_filter_bams_prepare_bam.log"
if \
    run_capture \
        "prepare execute filter-bams BAM fixture" "${log}" \
        run_samtools view -bS -o "${in_bam}" "${in_sam}" \
    && run_capture \
        "index execute filter-bams BAM fixture" "${log}.index" \
        run_samtools index "${in_bam}"
then
    rec_pass "filter BAM input fixture is prepared"
else
    rec_fail "failed to prepare filter BAM fixture; see $(rec_relpath "${log}")"
    finish
    exit $?
fi


#  Run execute_filter_bams.sh for one retain mode through submit_filter_bams.sh
function run_case_filter() {
    local retain="${1:-}"
    local nam_case="${2:-}"
    local log_lcl="${3:-}"

    shift 3

    if \
        run_capture \
            "execute filter-bams ${nam_case}" "${log_lcl}" \
            "${TEST_BASH}" "${ROOT_REPO}/scripts/execute_filter_bams.sh" \
                --threads 1 \
                --csv_infile "${in_bam}" \
                --dir_out "${dir_out}" \
                --retain "${retain}" \
                --err_out "${dir_err}" \
                --nam_job "test_execute_filter_bams_${nam_case}" \
                --max_job 1 \
                "$@"
    then
        rec_pass "execute_filter_bams.sh retain=${retain} ${nam_case} exits 0"
    else
        rec_fail \
            "execute_filter_bams.sh retain=${retain} ${nam_case} failed; see" \
            "$(rec_relpath "${log_lcl}")"
    fi
}


#  S. cerevisiae filtering should retain only canonical SC chromosomes
outfile="${dir_out}/filter_sc_sp.sc.bam"
idxstats="${dir_out}/filter_sc_sp.sc.idxstats.txt"
log="${dir_log}/execute_filter_bams_sc.log"

run_case_filter "sc" "sc" "${log}"

assert_file_nonempty "${outfile}" "execute retain=sc BAM output"
assert_file_nonempty "${outfile}.bai" "execute retain=sc BAM index"

if [[ -s "${outfile}" ]]; then
    run_capture \
        "idxstats execute filter-bams sc" "${idxstats}" \
        run_samtools idxstats "${outfile}"

    assert_grep_pattern "${idxstats}" $'^I\t100\t1\t0$' \
        "execute retain=sc output keeps chromosome I"
    assert_no_grep_pattern "${idxstats}" $'^Mito\t' \
        "execute retain=sc output omits Mito without --mito"
    assert_no_grep_pattern "${idxstats}" $'^SP_I\t' \
        "execute retain=sc output omits S. pombe chromosomes"
fi


#  S. pombe filtering should honor optional TG, MTR, and mito contigs
outfile="${dir_out}/filter_sc_sp.sp.bam"
idxstats="${dir_out}/filter_sc_sp.sp.idxstats.txt"
log="${dir_log}/execute_filter_bams_sp.log"

run_case_filter "sp" "sp" "${log}" --tg --mtr --mito

assert_file_nonempty "${outfile}" "execute retain=sp BAM output"
assert_file_nonempty "${outfile}.bai" "execute retain=sp BAM index"

if [[ -s "${outfile}" ]]; then
    run_capture \
        "idxstats execute filter-bams sp" "${idxstats}" \
        run_samtools idxstats "${outfile}"

    assert_grep_pattern "${idxstats}" $'^SP_I\t100\t1\t0$' \
        "execute retain=sp output keeps SP_I"
    assert_grep_pattern "${idxstats}" $'^SP_II_TG\t100\t1\t0$' \
        "execute retain=sp output keeps SP_II_TG"
    assert_grep_pattern "${idxstats}" $'^SP_MTR\t100\t1\t0$' \
        "execute retain=sp output keeps SP_MTR"
    assert_grep_pattern "${idxstats}" $'^SP_Mito\t100\t1\t0$' \
        "execute retain=sp output keeps SP_Mito"
    assert_no_grep_pattern "${idxstats}" $'^I\t' \
        "execute retain=sp output omits S. cerevisiae chromosomes"
    assert_no_grep_pattern "${idxstats}" $'^chrUn\t' \
        "execute retain=sp output omits unrelated contigs"
fi

finish
