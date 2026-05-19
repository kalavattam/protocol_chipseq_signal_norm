#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_filter_bams_cram_input_bam.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="execute filter-bams CRAM input BAM output"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

rec_section "${TEST_NAME}"

dir_fx="${ROOT_REPO}/tests/filter_bams/fixtures"
in_sam="${dir_fx}/sam/filter_sc_sp.sam"
ref_fa="${dir_fx}/reference/filter_sc_sp.fa"
ref_fai="${ref_fa}.fai"

tmp="${TEST_DIR_TMP}/execute_filter_bams_cram_input_bam"
dir_in="${tmp}/input"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/filter_bams"
in_cram="${dir_in}/filter_sc_sp.cram"

rm -rf "${tmp}"
mkdir -p "${dir_in}" "${dir_out}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_exist "${in_sam}" "${ref_fa}" "${ref_fai}" || {
    finish
    exit $?
}


#  Build deterministic CRAM input from committed SAM and reference fixtures
log="${dir_log}/execute_filter_bams_prepare_cram.log"
if \
    run_capture \
        "prepare execute filter-bams CRAM fixture" "${log}" \
        run_samtools view -C -T "${ref_fa}" -o "${in_cram}" "${in_sam}" \
    && run_capture \
        "index execute filter-bams CRAM fixture" "${log}.index" \
        run_samtools index "${in_cram}"
then
    rec_pass "filter CRAM input fixture is prepared"
else
    rec_fail "failed to prepare filter CRAM fixture; see $(rec_relpath "${log}")"
    finish
    exit $?
fi


function run_case_filter() {
    local retain="${1:-}"
    local nam_case="${2:-}"
    local log_lcl="${3:-}"

    shift 3

    if \
        run_capture \
            "execute filter-bams CRAM ${nam_case}" "${log_lcl}" \
            "${TEST_BASH}" "${ROOT_REPO}/scripts/execute_filter_bams.sh" \
                --threads 1 \
                --csv_infile "${in_cram}" \
                --dir_out "${dir_out}" \
                --retain "${retain}" \
                --ref_fa "${ref_fa}" \
                --err_out "${dir_err}" \
                --nam_job "test_execute_filter_bams_cram_${nam_case}" \
                --max_job 1 \
                "$@"
    then
        rec_pass "execute_filter_bams.sh CRAM retain=${retain} ${nam_case} exits 0"
    else
        rec_fail \
            "execute_filter_bams.sh CRAM retain=${retain} ${nam_case} failed; see" \
            "$(rec_relpath "${log_lcl}")"
    fi
}


#  S. cerevisiae filtering should retain only canonical SC chromosomes
outfile="${dir_out}/filter_sc_sp.sc.bam"
idxstats="${dir_out}/filter_sc_sp.sc.idxstats.txt"
log="${dir_log}/execute_filter_bams_cram_sc.log"

run_case_filter "sc" "sc" "${log}"

assert_file_nonempty "${outfile}" "execute CRAM retain=sc BAM output"
assert_file_nonempty "${outfile}.bai" "execute CRAM retain=sc BAM index"
assert_file_exists \
    "${dir_err}/test_execute_filter_bams_cram_sc_ser.stdout.txt" \
    "execute CRAM retain=sc serial stdout log"
assert_file_exists \
    "${dir_err}/test_execute_filter_bams_cram_sc_ser.stderr.txt" \
    "execute CRAM retain=sc serial stderr log"
assert_file_exists \
    "${dir_err}/test_execute_filter_bams_cram_sc.filter_sc_sp.stdout.txt" \
    "execute CRAM retain=sc submit stdout log"
assert_file_exists \
    "${dir_err}/test_execute_filter_bams_cram_sc.filter_sc_sp.stderr.txt" \
    "execute CRAM retain=sc submit stderr log"

if [[ -s "${outfile}" ]]; then
    if \
        run_capture \
            "quickcheck execute filter-bams CRAM sc" \
            "${dir_log}/execute_filter_bams_cram_sc_quickcheck.log" \
            run_samtools quickcheck "${outfile}"
    then
        rec_pass "execute CRAM retain=sc BAM passes samtools quickcheck"
    else
        rec_fail "execute CRAM retain=sc BAM fails samtools quickcheck"
    fi

    run_capture \
        "idxstats execute filter-bams CRAM sc" "${idxstats}" \
        run_samtools idxstats "${outfile}"

    assert_grep_pattern "${idxstats}" $'^I\t100\t1\t0$' \
        "execute CRAM retain=sc output keeps chromosome I"
    assert_no_grep_pattern "${idxstats}" $'^Mito\t' \
        "execute CRAM retain=sc output omits Mito without --mito"
    assert_no_grep_pattern "${idxstats}" $'^SP_I\t' \
        "execute CRAM retain=sc output omits S. pombe chromosomes"
fi


#  S. pombe filtering should honor optional TG, MTR, and mito contigs
outfile="${dir_out}/filter_sc_sp.sp.bam"
idxstats="${dir_out}/filter_sc_sp.sp.idxstats.txt"
log="${dir_log}/execute_filter_bams_cram_sp.log"

run_case_filter "sp" "sp" "${log}" --tg --mtr --mito

assert_file_nonempty "${outfile}" "execute CRAM retain=sp BAM output"
assert_file_nonempty "${outfile}.bai" "execute CRAM retain=sp BAM index"
assert_file_exists \
    "${dir_err}/test_execute_filter_bams_cram_sp_ser.stdout.txt" \
    "execute CRAM retain=sp serial stdout log"
assert_file_exists \
    "${dir_err}/test_execute_filter_bams_cram_sp_ser.stderr.txt" \
    "execute CRAM retain=sp serial stderr log"
assert_file_exists \
    "${dir_err}/test_execute_filter_bams_cram_sp.filter_sc_sp.stdout.txt" \
    "execute CRAM retain=sp submit stdout log"
assert_file_exists \
    "${dir_err}/test_execute_filter_bams_cram_sp.filter_sc_sp.stderr.txt" \
    "execute CRAM retain=sp submit stderr log"

if [[ -s "${outfile}" ]]; then
    if \
        run_capture \
            "quickcheck execute filter-bams CRAM sp" \
            "${dir_log}/execute_filter_bams_cram_sp_quickcheck.log" \
            run_samtools quickcheck "${outfile}"
    then
        rec_pass "execute CRAM retain=sp BAM passes samtools quickcheck"
    else
        rec_fail "execute CRAM retain=sp BAM fails samtools quickcheck"
    fi

    run_capture \
        "idxstats execute filter-bams CRAM sp" "${idxstats}" \
        run_samtools idxstats "${outfile}"

    assert_grep_pattern "${idxstats}" $'^SP_I\t100\t1\t0$' \
        "execute CRAM retain=sp output keeps SP_I"
    assert_grep_pattern "${idxstats}" $'^SP_II_TG\t100\t1\t0$' \
        "execute CRAM retain=sp output keeps SP_II_TG"
    assert_grep_pattern "${idxstats}" $'^SP_MTR\t100\t1\t0$' \
        "execute CRAM retain=sp output keeps SP_MTR"
    assert_grep_pattern "${idxstats}" $'^SP_Mito\t100\t1\t0$' \
        "execute CRAM retain=sp output keeps SP_Mito"
    assert_no_grep_pattern "${idxstats}" $'^I\t' \
        "execute CRAM retain=sp output omits S. cerevisiae chromosomes"
    assert_no_grep_pattern "${idxstats}" $'^chrUn\t' \
        "execute CRAM retain=sp output omits unrelated contigs"
fi

finish
