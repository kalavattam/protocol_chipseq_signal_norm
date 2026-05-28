#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_submit_align_fastqs_bwa_mem2.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="submit align-fastqs bwa-mem2 MEM"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

print_section "${TEST_NAME}"

#  Define fixture and output paths for bwa-mem2 MEM alignment
dir_fx="${ROOT_REPO}/tests/align_fastqs/fixtures"
in_se="${dir_fx}/fastq/tiny_se.atria.fastq.gz"
in_pe_1="${dir_fx}/fastq/tiny_pe_R1.atria.fastq.gz"
in_pe_2="${dir_fx}/fastq/tiny_pe_R2.atria.fastq.gz"
in_pe="${in_pe_1},${in_pe_2}"
idx="${dir_fx}/bwa-mem2/tiny.fa"

tmp="${TEST_DIR_TMP}/submit_align_fastqs_bwa_mem2"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/align_fastqs"

out_se="${dir_out}/tiny_se.bam"
stat_se="${dir_out}/tiny_se.idxstats.txt"
vw_se="${dir_out}/tiny_se.view.txt"

out_pe="${dir_out}/tiny_pe.bam"
stat_pe="${dir_out}/tiny_pe.idxstats.txt"
vw_pe="${dir_out}/tiny_pe.view.txt"

rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty \
    "${in_se}" \
    "${in_pe_1}" \
    "${in_pe_2}" \
    "${idx}" \
    "${idx}.0123" \
    "${idx}.amb" \
    "${idx}.ann" \
    "${idx}.bwt.2bit.64" \
    "${idx}.pac" \
    || {
        finish
        exit $?
    }


#  Align one SE FASTQ fixture with bwa-mem2 MEM and emit BAM
log="${dir_log}/submit_align_fastqs_bwa_mem2_se.log"

# shellcheck disable=SC2154
if \
    run_capture \
        "submit align-fastqs bwa-mem2 mem se" \
        "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/submit_align_fastqs.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/scripts" \
            --threads 1 \
            --aligner bwa-mem2 \
            --mapq 0 \
            --index "${idx}" \
            --csv_infile "${in_se}" \
            --dir_out "${dir_out}" \
            --out_ext bam \
            --sfx_se ".atria.fastq.gz" \
            --sfx_pe "_R1.atria.fastq.gz" \
            --err_out "${dir_err}" \
            --nam_job "test_submit_align_bwa_mem2_se"
then
    record_pass "submit_align_fastqs.sh bwa-mem2 MEM SE BAM exits 0"
else
    record_fail \
        "submit_align_fastqs.sh bwa-mem2 MEM SE BAM failed; see" \
        "$(print_relpath "${log}")"
fi

assert_file_nonempty \
    "${out_se}" \
    "submit bwa-mem2 MEM SE BAM output"
assert_file_nonempty \
    "${out_se}.bai" \
    "submit bwa-mem2 MEM SE BAM index"

if [[ -s "${out_se}" ]]; then
    if \
        run_capture \
            "quickcheck submit align-fastqs bwa-mem2 MEM SE BAM" \
            "${dir_log}/submit_align_fastqs_bwa_mem2_se_quickcheck.log" \
            run_samtools quickcheck "${out_se}"
    then
        record_pass "submit bwa-mem2 MEM SE BAM passes samtools quickcheck"
    else
        record_fail "submit bwa-mem2 MEM SE BAM fails samtools quickcheck"
    fi

    run_capture \
        "idxstats submit align-fastqs bwa-mem2 MEM SE BAM" \
        "${stat_se}" \
        run_samtools idxstats "${out_se}"
    run_capture \
        "view submit align-fastqs bwa-mem2 MEM SE BAM" \
        "${vw_se}" \
        run_samtools view "${out_se}"

    assert_pattern_found \
        "${stat_se}" \
        $'^I\t108\t1\t0$' \
        "submit bwa-mem2 MEM SE BAM has one mapped read on chromosome I"
    assert_pattern_found \
        "${vw_se}" \
        $'^tiny_se_read_1\t' \
        "submit bwa-mem2 MEM SE BAM contains expected read name"
fi


#  Align one PE FASTQ fixture with bwa-mem2 MEM and emit BAM
log="${dir_log}/submit_align_fastqs_bwa_mem2_pe.log"

# shellcheck disable=SC2154
if \
    run_capture \
        "submit align-fastqs bwa-mem2 mem pe" \
        "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/submit_align_fastqs.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/scripts" \
            --threads 1 \
            --aligner bwa-mem2 \
            --mapq 0 \
            --index "${idx}" \
            --csv_infile "${in_pe}" \
            --dir_out "${dir_out}" \
            --out_ext bam \
            --sfx_se ".atria.fastq.gz" \
            --sfx_pe "_R1.atria.fastq.gz" \
            --err_out "${dir_err}" \
            --nam_job "test_submit_align_bwa_mem2_pe"
then
    record_pass "submit_align_fastqs.sh bwa-mem2 MEM PE BAM exits 0"
else
    record_fail \
        "submit_align_fastqs.sh bwa-mem2 MEM PE BAM failed; see" \
        "$(print_relpath "${log}")"
fi

assert_file_nonempty \
    "${out_pe}" \
    "submit bwa-mem2 MEM PE BAM output"
assert_file_nonempty \
    "${out_pe}.bai" \
    "submit bwa-mem2 MEM PE BAM index"

if [[ -s "${out_pe}" ]]; then
    if \
        run_capture \
            "quickcheck submit align-fastqs bwa-mem2 MEM PE BAM" \
            "${dir_log}/submit_align_fastqs_bwa_mem2_pe_quickcheck.log" \
            run_samtools quickcheck "${out_pe}"
    then
        record_pass "submit bwa-mem2 MEM PE BAM passes samtools quickcheck"
    else
        record_fail "submit bwa-mem2 MEM PE BAM fails samtools quickcheck"
    fi

    run_capture \
        "idxstats submit align-fastqs bwa-mem2 MEM PE BAM" \
        "${stat_pe}" \
        run_samtools idxstats "${out_pe}"
    run_capture \
        "view submit align-fastqs bwa-mem2 MEM PE BAM" \
        "${vw_pe}" \
        run_samtools view "${out_pe}"

    assert_pattern_found \
        "${stat_pe}" \
        $'^I\t108\t2\t0$' \
        "submit bwa-mem2 MEM PE BAM has two mapped reads on chromosome I"
    assert_pattern_found \
        "${vw_pe}" \
        $'^tiny_pe_pair_1\t' \
        "submit bwa-mem2 MEM PE BAM contains expected read name"
fi

finish
