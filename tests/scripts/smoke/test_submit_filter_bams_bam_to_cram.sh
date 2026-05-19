#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_submit_filter_bams_bam_to_cram.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="submit filter-bams BAM to CRAM"

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

tmp="${TEST_DIR_TMP}/submit_filter_bams_bam_to_cram"
dir_in="${tmp}/input"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_missing="${tmp}/missing_ref"
dir_log="${TEST_DIR_LOG}/filter_bams"
in_bam="${dir_in}/filter_sc_sp.bam"

rm -rf "${tmp}"
mkdir -p "${dir_in}" "${dir_out}" "${dir_err}" "${dir_missing}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_exist "${in_sam}" "${ref_fa}" "${ref_fai}" || {
    finish
    exit $?
}


#  Build deterministic BAM input from committed SAM fixture
log="${dir_log}/submit_filter_bams_prepare_bam_to_cram.log"
if \
    run_capture \
        "prepare submit filter-bams BAM fixture for CRAM output" "${log}" \
        run_samtools view -bS -o "${in_bam}" "${in_sam}" \
    && run_capture \
        "index submit filter-bams BAM fixture for CRAM output" "${log}.index" \
        run_samtools index "${in_bam}"
then
    rec_pass "filter BAM input fixture is prepared"
else
    rec_fail "failed to prepare filter BAM fixture; see $(rec_relpath "${log}")"
    finish
    exit $?
fi


function run_case_filter() {
    local retain="${1:-}"
    local nam_case="${2:-}"
    local log_lcl="${3:-}"

    shift 3

    # shellcheck disable=SC2154
    if \
        run_capture \
            "submit filter-bams BAM to CRAM ${nam_case}" "${log_lcl}" \
            "${TEST_BASH}" "${ROOT_REPO}/scripts/submit_filter_bams.sh" \
                --env_nam "${env_nam}" \
                --dir_scr "${ROOT_REPO}/scripts" \
                --threads 1 \
                --csv_infile "${in_bam}" \
                --dir_out "${dir_out}" \
                --out_ext cram \
                --retain "${retain}" \
                --ref_fa "${ref_fa}" \
                --err_out "${dir_err}" \
                --nam_job "test_submit_filter_bams_bam_to_cram_${nam_case}" \
                "$@"
    then
        rec_pass "submit_filter_bams.sh BAM to CRAM retain=${retain} ${nam_case} exits 0"
    else
        rec_fail \
            "submit_filter_bams.sh BAM to CRAM retain=${retain} ${nam_case} failed; see" \
            "$(rec_relpath "${log_lcl}")"
    fi
}


#  CRAM output without --ref_fa should fail clearly
log="${dir_log}/submit_filter_bams_bam_to_cram_missing_ref.log"
if \
    run_capture \
        "submit filter-bams BAM to CRAM missing ref" "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/submit_filter_bams.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/scripts" \
            --threads 1 \
            --csv_infile "${in_bam}" \
            --dir_out "${dir_missing}" \
            --out_ext cram \
            --retain sc \
            --err_out "${dir_err}" \
            --nam_job "test_submit_filter_bams_bam_to_cram_missing_ref"
then
    rec_fail "submit_filter_bams.sh --out_ext cram without --ref_fa unexpectedly passed"
else
    rec_pass "submit_filter_bams.sh --out_ext cram without --ref_fa fails"
fi

assert_grep_pattern "${log}" "'--ref_fa' is required" \
    "submit BAM-to-CRAM missing-ref error mentions --ref_fa"


#  S. cerevisiae filtering should retain only canonical SC chromosomes
outfile="${dir_out}/filter_sc_sp.sc.cram"
log="${dir_log}/submit_filter_bams_bam_to_cram_sc.log"

run_case_filter "sc" "sc" "${log}"

assert_file_nonempty "${outfile}" "submit retain=sc CRAM output"
assert_cram_index "${outfile}" "submit retain=sc CRAI index"
assert_file_exists \
    "${dir_err}/test_submit_filter_bams_bam_to_cram_sc.filter_sc_sp.stdout.txt" \
    "submit retain=sc CRAM stdout log"
assert_file_exists \
    "${dir_err}/test_submit_filter_bams_bam_to_cram_sc.filter_sc_sp.stderr.txt" \
    "submit retain=sc CRAM stderr log"

if [[ -s "${outfile}" ]]; then
    if \
        run_capture \
            "quickcheck submit filter-bams BAM to CRAM sc" \
            "${dir_log}/submit_filter_bams_bam_to_cram_sc_quickcheck.log" \
            run_samtools quickcheck "${outfile}"
    then
        rec_pass "submit retain=sc CRAM passes samtools quickcheck"
    else
        rec_fail "submit retain=sc CRAM fails samtools quickcheck"
    fi

    assert_cram_count "${outfile}" "${ref_fa}" I 1 \
        "submit retain=sc CRAM has one read on chromosome I" \
        "${dir_out}/filter_sc_sp.sc.I.count.txt"
    assert_cram_count "${outfile}" "${ref_fa}" Mito 0 \
        "submit retain=sc CRAM omits Mito without --mito" \
        "${dir_out}/filter_sc_sp.sc.Mito.count.txt"
    assert_cram_count "${outfile}" "${ref_fa}" SP_I 0 \
        "submit retain=sc CRAM omits S. pombe chromosomes" \
        "${dir_out}/filter_sc_sp.sc.SP_I.count.txt"
fi


#  S. pombe filtering should honor optional TG, MTR, and mito contigs
outfile="${dir_out}/filter_sc_sp.sp.cram"
log="${dir_log}/submit_filter_bams_bam_to_cram_sp.log"

run_case_filter "sp" "sp" "${log}" --tg --mtr --mito

assert_file_nonempty "${outfile}" "submit retain=sp CRAM output"
assert_cram_index "${outfile}" "submit retain=sp CRAI index"
assert_file_exists \
    "${dir_err}/test_submit_filter_bams_bam_to_cram_sp.filter_sc_sp.stdout.txt" \
    "submit retain=sp CRAM stdout log"
assert_file_exists \
    "${dir_err}/test_submit_filter_bams_bam_to_cram_sp.filter_sc_sp.stderr.txt" \
    "submit retain=sp CRAM stderr log"

if [[ -s "${outfile}" ]]; then
    if \
        run_capture \
            "quickcheck submit filter-bams BAM to CRAM sp" \
            "${dir_log}/submit_filter_bams_bam_to_cram_sp_quickcheck.log" \
            run_samtools quickcheck "${outfile}"
    then
        rec_pass "submit retain=sp CRAM passes samtools quickcheck"
    else
        rec_fail "submit retain=sp CRAM fails samtools quickcheck"
    fi

    assert_cram_count "${outfile}" "${ref_fa}" SP_I 1 \
        "submit retain=sp CRAM keeps SP_I" \
        "${dir_out}/filter_sc_sp.sp.SP_I.count.txt"
    assert_cram_count "${outfile}" "${ref_fa}" SP_II_TG 1 \
        "submit retain=sp CRAM keeps SP_II_TG" \
        "${dir_out}/filter_sc_sp.sp.SP_II_TG.count.txt"
    assert_cram_count "${outfile}" "${ref_fa}" SP_MTR 1 \
        "submit retain=sp CRAM keeps SP_MTR" \
        "${dir_out}/filter_sc_sp.sp.SP_MTR.count.txt"
    assert_cram_count "${outfile}" "${ref_fa}" SP_Mito 1 \
        "submit retain=sp CRAM keeps SP_Mito" \
        "${dir_out}/filter_sc_sp.sp.SP_Mito.count.txt"
    assert_cram_count "${outfile}" "${ref_fa}" I 0 \
        "submit retain=sp CRAM omits S. cerevisiae chromosomes" \
        "${dir_out}/filter_sc_sp.sp.I.count.txt"
    assert_cram_count "${outfile}" "${ref_fa}" chrUn 0 \
        "submit retain=sp CRAM omits unrelated contigs" \
        "${dir_out}/filter_sc_sp.sp.chrUn.count.txt"
fi

finish
