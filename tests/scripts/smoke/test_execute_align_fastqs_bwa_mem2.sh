#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_align_fastqs_bwa_mem2.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="execute align-fastqs bwa-mem2 MEM"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

rec_section "${TEST_NAME}"

#  Define fixture and output paths for execute-to-submit bwa-mem2 MEM alignment
dir_fx="${ROOT_REPO}/tests/align_fastqs/fixtures"
infile_se="${dir_fx}/fastq/tiny_se.atria.fastq.gz"
infile_pe_1="${dir_fx}/fastq/tiny_pe_R1.atria.fastq.gz"
infile_pe_2="${dir_fx}/fastq/tiny_pe_R2.atria.fastq.gz"
infile_pe="${infile_pe_1},${infile_pe_2}"
index="${dir_fx}/bwa-mem2/tiny.fa"

tmp="${TEST_DIR_TMP}/execute_align_fastqs_bwa_mem2"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/align_fastqs"

outfile_se="${dir_out}/tiny_se.bam"
idxstats_se="${dir_out}/tiny_se.idxstats.txt"
view_se="${dir_out}/tiny_se.view.txt"
outfile_pe="${dir_out}/tiny_pe.bam"
idxstats_pe="${dir_out}/tiny_pe.idxstats.txt"
view_pe="${dir_out}/tiny_pe.view.txt"

rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_exist \
    "${infile_se}" \
    "${infile_pe_1}" \
    "${infile_pe_2}" \
    "${index}" \
    "${index}.0123" \
    "${index}.amb" \
    "${index}.ann" \
    "${index}.bwt.2bit.64" \
    "${index}.pac" \
    || {
        finish
        exit $?
    }


#  Run samtools from the active shell or the resolved project environment
# shellcheck disable=SC2154
function run_samtools() {
    if check_cmd_exists samtools; then
        samtools "$@"
    else
        conda run -n "${env_nam}" samtools "$@"
    fi
}


#  Run execute_align_fastqs.sh through submit_align_fastqs.sh for SE input
log="${dir_log}/execute_align_fastqs_bwa_mem2_se.log"

if \
    run_capture \
        "execute align-fastqs bwa-mem2 mem se" "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/execute_align_fastqs.sh" \
            --threads 1 \
            --aligner bwa-mem2 \
            --mapq 0 \
            --index "${index}" \
            --csv_infile "${infile_se}" \
            --dir_out "${dir_out}" \
            --out_ext bam \
            --sfx_se ".atria.fastq.gz" \
            --sfx_pe "_R1.atria.fastq.gz" \
            --err_out "${dir_err}" \
            --nam_job "test_execute_align_bwa_mem2_se" \
            --max_job 1
then
    rec_pass "execute_align_fastqs.sh bwa-mem2 MEM SE BAM exits 0"
else
    rec_fail \
        "execute_align_fastqs.sh bwa-mem2 MEM SE BAM failed; see" \
        "$(rec_relpath "${log}")"
fi

assert_file_nonempty "${outfile_se}" "execute bwa-mem2 MEM SE BAM output"
assert_file_nonempty "${outfile_se}.bai" "execute bwa-mem2 MEM SE BAM index"

if [[ -s "${outfile_se}" ]]; then
    if run_capture \
        "quickcheck execute align-fastqs bwa-mem2 MEM SE BAM" \
        "${dir_log}/execute_align_fastqs_bwa_mem2_se_quickcheck.log" \
        run_samtools quickcheck "${outfile_se}"
    then
        rec_pass "execute bwa-mem2 MEM SE BAM passes samtools quickcheck"
    else
        rec_fail "execute bwa-mem2 MEM SE BAM fails samtools quickcheck"
    fi

    run_capture \
        "idxstats execute align-fastqs bwa-mem2 MEM SE BAM" "${idxstats_se}" \
        run_samtools idxstats "${outfile_se}"
    run_capture \
        "view execute align-fastqs bwa-mem2 MEM SE BAM" "${view_se}" \
        run_samtools view "${outfile_se}"

    assert_grep_pattern "${idxstats_se}" $'^I\t108\t1\t0$' \
        "execute bwa-mem2 MEM SE BAM has one mapped read on chromosome I"
    assert_grep_pattern "${view_se}" $'^tiny_se_read_1\t' \
        "execute bwa-mem2 MEM SE BAM contains expected read name"
fi


#  Run execute_align_fastqs.sh through submit_align_fastqs.sh for PE input
log="${dir_log}/execute_align_fastqs_bwa_mem2_pe.log"

if \
    run_capture \
        "execute align-fastqs bwa-mem2 mem pe" "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/execute_align_fastqs.sh" \
            --threads 1 \
            --aligner bwa-mem2 \
            --mapq 0 \
            --index "${index}" \
            --csv_infile "${infile_pe}" \
            --dir_out "${dir_out}" \
            --out_ext bam \
            --sfx_se ".atria.fastq.gz" \
            --sfx_pe "_R1.atria.fastq.gz" \
            --err_out "${dir_err}" \
            --nam_job "test_execute_align_bwa_mem2_pe" \
            --max_job 1
then
    rec_pass "execute_align_fastqs.sh bwa-mem2 MEM PE BAM exits 0"
else
    rec_fail \
        "execute_align_fastqs.sh bwa-mem2 MEM PE BAM failed; see" \
        "$(rec_relpath "${log}")"
fi

assert_file_nonempty "${outfile_pe}" "execute bwa-mem2 MEM PE BAM output"
assert_file_nonempty "${outfile_pe}.bai" "execute bwa-mem2 MEM PE BAM index"

if [[ -s "${outfile_pe}" ]]; then
    if run_capture \
        "quickcheck execute align-fastqs bwa-mem2 MEM PE BAM" \
        "${dir_log}/execute_align_fastqs_bwa_mem2_pe_quickcheck.log" \
        run_samtools quickcheck "${outfile_pe}"
    then
        rec_pass "execute bwa-mem2 MEM PE BAM passes samtools quickcheck"
    else
        rec_fail "execute bwa-mem2 MEM PE BAM fails samtools quickcheck"
    fi

    run_capture \
        "idxstats execute align-fastqs bwa-mem2 MEM PE BAM" "${idxstats_pe}" \
        run_samtools idxstats "${outfile_pe}"
    run_capture \
        "view execute align-fastqs bwa-mem2 MEM PE BAM" "${view_pe}" \
        run_samtools view "${outfile_pe}"

    assert_grep_pattern "${idxstats_pe}" $'^I\t108\t2\t0$' \
        "execute bwa-mem2 MEM PE BAM has two mapped reads on chromosome I"
    assert_grep_pattern "${view_pe}" $'^tiny_pe_pair_1\t' \
        "execute bwa-mem2 MEM PE BAM contains expected read name"
fi

finish
