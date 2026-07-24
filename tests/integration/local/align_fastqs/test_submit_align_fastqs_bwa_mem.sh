#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_submit_align_fastqs_bwa_mem.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="submit align-fastqs BWA MEM"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


#  Define fixture and output paths for BWA MEM alignment
dir_fx="${ROOT_REPO}/tests/fixtures/align_fastqs"
in_se="${dir_fx}/fastq/se/tiny_se.atria.fastq.gz"
in_pe_1="${dir_fx}/fastq/pe/tiny_pe_R1.atria.fastq.gz"
in_pe_2="${dir_fx}/fastq/pe/tiny_pe_R2.atria.fastq.gz"
in_pe="${in_pe_1},${in_pe_2}"
idx_bwa="${dir_fx}/bwa/tiny.fa"

tmp="${TEST_DIR_TMP}/submit_align_fastqs_bwa_mem"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/align_fastqs"

out_se="${dir_out}/tiny_se.bam"
bai_se="${out_se}.bai"
stat_se="${dir_out}/tiny_se.idxstats.txt"
vw_se="${dir_out}/tiny_se.view.txt"

out_pe="${dir_out}/tiny_pe.bam"
bai_pe="${out_pe}.bai"
stat_pe="${dir_out}/tiny_pe.idxstats.txt"
vw_pe="${dir_out}/tiny_pe.view.txt"

log_se="${dir_log}/submit_align_fastqs_bwa_mem_se.log"
log_pe="${dir_log}/submit_align_fastqs_bwa_mem_pe.log"
log_se_qc="${dir_log}/submit_align_fastqs_bwa_mem_se_quickcheck.log"
log_pe_qc="${dir_log}/submit_align_fastqs_bwa_mem_pe_quickcheck.log"


print_section "${TEST_NAME}"

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
    "${idx_bwa}" \
    "${idx_bwa}.amb" \
    "${idx_bwa}.ann" \
    "${idx_bwa}.bwt" \
    "${idx_bwa}.pac" \
    "${idx_bwa}.sa" \
    || {
        finish
        exit $?
    }


#  Align one SE FASTQ fixture with BWA MEM and emit BAM
# shellcheck disable=SC2154
if \
    run_capture \
        "submit align-fastqs bwa mem se" \
        "${log_se}" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/submit_align_fastqs.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/bin" \
            --threads 1 \
            --aligner bwa \
            --bwa_alg mem \
            --mapq 0 \
            --index "${idx_bwa}" \
            --csv_fil_in "${in_se}" \
            --dir_out "${dir_out}" \
            --out_ext bam \
            --sfx_se ".atria.fastq.gz" \
            --sfx_pe "_R1.atria.fastq.gz" \
            --dir_eo "${dir_err}" \
            --nam_job "test_submit_align_bwa_mem_se"
then
    record_pass "submit_align_fastqs.sh BWA MEM SE BAM exits 0"
else
    record_fail \
        "submit_align_fastqs.sh BWA MEM SE BAM failed; see" \
        "$(print_relpath "${log_se}")"
fi

assert_file_nonempty \
    "${out_se}" \
    "submit BWA MEM SE BAM output"

assert_file_nonempty \
    "${bai_se}" \
    "submit BWA MEM SE BAM index"

if [[ -s "${out_se}" ]]; then
    if \
        run_capture \
            "quickcheck submit align-fastqs BWA MEM SE BAM" \
            "${log_se_qc}" \
            run_samtools quickcheck "${out_se}"
    then
        record_pass "submit BWA MEM SE BAM passes samtools quickcheck"
    else
        record_fail "submit BWA MEM SE BAM fails samtools quickcheck"
    fi

    run_capture \
        "idxstats submit align-fastqs BWA MEM SE BAM" \
        "${stat_se}" \
        run_samtools idxstats "${out_se}"

    run_capture \
        "view submit align-fastqs BWA MEM SE BAM" \
        "${vw_se}" \
        run_samtools view "${out_se}"

    assert_pattern_found \
        "${stat_se}" \
        $'^I\t108\t1\t0$' \
        "submit BWA MEM SE BAM has one mapped read on chromosome I"

    assert_pattern_found \
        "${vw_se}" \
        $'^tiny_se_read_1\t' \
        "submit BWA MEM SE BAM contains expected read name"
fi


#  Align one PE FASTQ fixture with BWA MEM and emit BAM
# shellcheck disable=SC2154
if \
    run_capture \
        "submit align-fastqs bwa mem pe" \
        "${log_pe}" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/submit_align_fastqs.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/bin" \
            --threads 1 \
            --aligner bwa \
            --bwa_alg mem \
            --mapq 0 \
            --index "${idx_bwa}" \
            --csv_fil_in "${in_pe}" \
            --dir_out "${dir_out}" \
            --out_ext bam \
            --sfx_se ".atria.fastq.gz" \
            --sfx_pe "_R1.atria.fastq.gz" \
            --dir_eo "${dir_err}" \
            --nam_job "test_submit_align_bwa_mem_pe"
then
    record_pass "submit_align_fastqs.sh BWA MEM PE BAM exits 0"
else
    record_fail \
        "submit_align_fastqs.sh BWA MEM PE BAM failed; see" \
        "$(print_relpath "${log_pe}")"
fi

assert_file_nonempty \
    "${out_pe}" \
    "submit BWA MEM PE BAM output"

assert_file_nonempty \
    "${bai_pe}" \
    "submit BWA MEM PE BAM index"

if [[ -s "${out_pe}" ]]; then
    if \
        run_capture \
            "quickcheck submit align-fastqs BWA MEM PE BAM" \
            "${log_pe_qc}" \
            run_samtools quickcheck "${out_pe}"
    then
        record_pass "submit BWA MEM PE BAM passes samtools quickcheck"
    else
        record_fail "submit BWA MEM PE BAM fails samtools quickcheck"
    fi

    run_capture \
        "idxstats submit align-fastqs BWA MEM PE BAM" \
        "${stat_pe}" \
        run_samtools idxstats "${out_pe}"

    run_capture \
        "view submit align-fastqs BWA MEM PE BAM" \
        "${vw_pe}" \
        run_samtools view "${out_pe}"

    assert_pattern_found \
        "${stat_pe}" \
        $'^I\t108\t2\t0$' \
        "submit BWA MEM PE BAM has two mapped reads on chromosome I"

    assert_pattern_found \
        "${vw_pe}" \
        $'^tiny_pe_pair_1\t' \
        "submit BWA MEM PE BAM contains expected read name"
fi

finish
