#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_align_fastqs_bwa_aln.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6);
# - Anthropic Claude Code (Opus 5.5).
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="execute align-fastqs BWA-backtrack"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Define fixture and output paths for execute-to-submit BWA-backtrack
# alignment.
dir_fix="${ROOT_REPO}/tests/fixtures/align_fastqs"
in_se="${dir_fix}/fastq/se/tiny_se.atria.fastq.gz"
in_pe_1="${dir_fix}/fastq/pe/tiny_pe_R1.atria.fastq.gz"
in_pe_2="${dir_fix}/fastq/pe/tiny_pe_R2.atria.fastq.gz"
in_pe="${in_pe_1},${in_pe_2}"
idx_bwa="${dir_fix}/bwa/tiny.fa"

tmp="${TEST_DIR_TMP}/execute_align_fastqs_bwa_aln"
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

log_se="${dir_log}/execute_align_fastqs_bwa_aln_se.log"
log_pe="${dir_log}/execute_align_fastqs_bwa_aln_pe.log"
log_se_qc="${dir_log}/execute_align_fastqs_bwa_aln_se_quickcheck.log"
log_pe_qc="${dir_log}/execute_align_fastqs_bwa_aln_pe_quickcheck.log"


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


# Run 'execute_align_fastqs.sh' through 'submit_align_fastqs.sh' for SE input.
if \
    run_capture \
        "execute align-fastqs bwa aln se" \
        "${log_se}" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/execute_align_fastqs.sh" \
            --env_nam "${env_nam}" \
            --threads 1 \
            --aligner bwa \
            --bwa_alg aln \
            --mapq 0 \
            --index "${idx_bwa}" \
            --csv_fil_in "${in_se}" \
            --dir_out "${dir_out}" \
            --out_ext bam \
            --sfx_se ".atria.fastq.gz" \
            --sfx_pe "_R1.atria.fastq.gz" \
            --dir_eo "${dir_err}" \
            --nam_job "test_execute_align_bwa_aln_se" \
            --max_job 1
then
    record_pass "execute_align_fastqs.sh BWA-backtrack SE BAM exits 0"
else
    record_fail \
        "execute_align_fastqs.sh BWA-backtrack SE BAM failed; see" \
        "$(print_relpath "${log_se}")"
fi

assert_file_nonempty \
    "${out_se}" \
    "execute BWA-backtrack SE BAM output"

assert_file_nonempty \
    "${bai_se}" \
    "execute BWA-backtrack SE BAM index"

if [[ -s "${out_se}" ]]; then
    if \
        run_capture \
            "quickcheck execute align-fastqs BWA-backtrack SE BAM" \
            "${log_se_qc}" \
            run_samtools quickcheck "${out_se}"
    then
        record_pass "execute BWA-backtrack SE BAM passes samtools quickcheck"
    else
        record_fail "execute BWA-backtrack SE BAM fails samtools quickcheck"
    fi

    run_capture \
        "idxstats execute align-fastqs BWA-backtrack SE BAM" \
        "${stat_se}" \
        run_samtools idxstats "${out_se}"

    run_capture \
        "view execute align-fastqs BWA-backtrack SE BAM" \
        "${vw_se}" \
        run_samtools view "${out_se}"

    assert_pattern_found \
        "${stat_se}" \
        $'^I\t108\t1\t0$' \
        "execute BWA-backtrack SE BAM has one mapped read on chromosome I"

    assert_pattern_found \
        "${vw_se}" \
        $'^tiny_se_read_1\t' \
        "execute BWA-backtrack SE BAM contains expected read name"
fi


# Run 'execute_align_fastqs.sh' through 'submit_align_fastqs.sh' for PE input.
if \
    run_capture \
        "execute align-fastqs bwa aln pe" \
        "${log_pe}" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/execute_align_fastqs.sh" \
            --env_nam "${env_nam}" \
            --threads 1 \
            --aligner bwa \
            --bwa_alg aln \
            --mapq 0 \
            --index "${idx_bwa}" \
            --csv_fil_in "${in_pe}" \
            --dir_out "${dir_out}" \
            --out_ext bam \
            --sfx_se ".atria.fastq.gz" \
            --sfx_pe "_R1.atria.fastq.gz" \
            --dir_eo "${dir_err}" \
            --nam_job "test_execute_align_bwa_aln_pe" \
            --max_job 1
then
    record_pass "execute_align_fastqs.sh BWA-backtrack PE BAM exits 0"
else
    record_fail \
        "execute_align_fastqs.sh BWA-backtrack PE BAM failed; see" \
        "$(print_relpath "${log_pe}")"
fi

assert_file_nonempty \
    "${out_pe}" \
    "execute BWA-backtrack PE BAM output"

assert_file_nonempty \
    "${bai_pe}" \
    "execute BWA-backtrack PE BAM index"

if [[ -s "${out_pe}" ]]; then
    if \
        run_capture \
            "quickcheck execute align-fastqs BWA-backtrack PE BAM" \
            "${log_pe_qc}" \
            run_samtools quickcheck "${out_pe}"
    then
        record_pass "execute BWA-backtrack PE BAM passes samtools quickcheck"
    else
        record_fail "execute BWA-backtrack PE BAM fails samtools quickcheck"
    fi

    run_capture \
        "idxstats execute align-fastqs BWA-backtrack PE BAM" \
        "${stat_pe}" \
        run_samtools idxstats "${out_pe}"

    run_capture \
        "view execute align-fastqs BWA-backtrack PE BAM" \
        "${vw_pe}" \
        run_samtools view "${out_pe}"

    assert_pattern_found \
        "${stat_pe}" \
        $'^I\t108\t2\t0$' \
        "execute BWA-backtrack PE BAM has two mapped reads on chromosome I"

    assert_pattern_found \
        "${vw_pe}" \
        $'^tiny_pe_pair_1\t' \
        "execute BWA-backtrack PE BAM contains expected read name"

    assert_pattern_found \
        "${vw_pe}" \
        $'^tiny_pe_pair_1\t99\tI\t17\t' \
        "execute BWA-backtrack PE BAM has proper-pair R1 flag and start"

    assert_pattern_found \
        "${vw_pe}" \
        $'^tiny_pe_pair_1\t147\tI\t70\t' \
        "execute BWA-backtrack PE BAM has proper-pair R2 flag and start"
fi


finish
