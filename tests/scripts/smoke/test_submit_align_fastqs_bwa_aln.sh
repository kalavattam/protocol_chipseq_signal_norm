#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_submit_align_fastqs_bwa_aln.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="submit align-fastqs BWA ALN"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

rec_section "${TEST_NAME}"

#  Define fixture and output paths for BWA ALN alignment
dir_fx="${ROOT_REPO}/tests/align_fastqs/fixtures"
infile_se="${dir_fx}/fastq/tiny_se.atria.fastq.gz"
infile_pe_1="${dir_fx}/fastq/tiny_pe_R1.atria.fastq.gz"
infile_pe_2="${dir_fx}/fastq/tiny_pe_R2.atria.fastq.gz"
infile_pe="${infile_pe_1},${infile_pe_2}"
index="${dir_fx}/bwa/tiny.fa"

tmp="${TEST_DIR_TMP}/submit_align_fastqs_bwa_aln"
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
    "${index}.amb" \
    "${index}.ann" \
    "${index}.bwt" \
    "${index}.pac" \
    "${index}.sa" \
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


#  Align one SE FASTQ fixture with BWA ALN and emit BAM
log="${dir_log}/submit_align_fastqs_bwa_aln_se.log"

# shellcheck disable=SC2154
if \
    run_capture \
        "submit align-fastqs bwa aln se" "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/submit_align_fastqs.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/scripts" \
            --threads 1 \
            --aligner bwa \
            --bwa_alg aln \
            --mapq 0 \
            --index "${index}" \
            --csv_infile "${infile_se}" \
            --dir_out "${dir_out}" \
            --out_ext bam \
            --sfx_se ".atria.fastq.gz" \
            --sfx_pe "_R1.atria.fastq.gz" \
            --err_out "${dir_err}" \
            --nam_job "test_submit_align_bwa_aln_se"
then
    rec_pass "submit_align_fastqs.sh BWA ALN SE BAM exits 0"
else
    rec_fail \
        "submit_align_fastqs.sh BWA ALN SE BAM failed; see" \
        "$(rec_relpath "${log}")"
fi

assert_file_nonempty "${outfile_se}" "submit BWA ALN SE BAM output"
assert_file_nonempty "${outfile_se}.bai" "submit BWA ALN SE BAM index"

if [[ -s "${outfile_se}" ]]; then
    if run_capture \
        "quickcheck submit align-fastqs BWA ALN SE BAM" \
        "${dir_log}/submit_align_fastqs_bwa_aln_se_quickcheck.log" \
        run_samtools quickcheck "${outfile_se}"
    then
        rec_pass "submit BWA ALN SE BAM passes samtools quickcheck"
    else
        rec_fail "submit BWA ALN SE BAM fails samtools quickcheck"
    fi

    run_capture \
        "idxstats submit align-fastqs BWA ALN SE BAM" "${idxstats_se}" \
        run_samtools idxstats "${outfile_se}"
    run_capture \
        "view submit align-fastqs BWA ALN SE BAM" "${view_se}" \
        run_samtools view "${outfile_se}"

    assert_grep_pattern "${idxstats_se}" $'^I\t108\t1\t0$' \
        "submit BWA ALN SE BAM has one mapped read on chromosome I"
    assert_grep_pattern "${view_se}" $'^tiny_se_read_1\t' \
        "submit BWA ALN SE BAM contains expected read name"
fi


#  Align one PE FASTQ fixture with BWA ALN and emit BAM
log="${dir_log}/submit_align_fastqs_bwa_aln_pe.log"

# shellcheck disable=SC2154
if \
    run_capture \
        "submit align-fastqs bwa aln pe" "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/submit_align_fastqs.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/scripts" \
            --threads 1 \
            --aligner bwa \
            --bwa_alg aln \
            --mapq 0 \
            --index "${index}" \
            --csv_infile "${infile_pe}" \
            --dir_out "${dir_out}" \
            --out_ext bam \
            --sfx_se ".atria.fastq.gz" \
            --sfx_pe "_R1.atria.fastq.gz" \
            --err_out "${dir_err}" \
            --nam_job "test_submit_align_bwa_aln_pe"
then
    rec_pass "submit_align_fastqs.sh BWA ALN PE BAM exits 0"
else
    rec_fail \
        "submit_align_fastqs.sh BWA ALN PE BAM failed; see" \
        "$(rec_relpath "${log}")"
fi

assert_file_nonempty "${outfile_pe}" "submit BWA ALN PE BAM output"
assert_file_nonempty "${outfile_pe}.bai" "submit BWA ALN PE BAM index"

if [[ -s "${outfile_pe}" ]]; then
    if run_capture \
        "quickcheck submit align-fastqs BWA ALN PE BAM" \
        "${dir_log}/submit_align_fastqs_bwa_aln_pe_quickcheck.log" \
        run_samtools quickcheck "${outfile_pe}"
    then
        rec_pass "submit BWA ALN PE BAM passes samtools quickcheck"
    else
        rec_fail "submit BWA ALN PE BAM fails samtools quickcheck"
    fi

    run_capture \
        "idxstats submit align-fastqs BWA ALN PE BAM" "${idxstats_pe}" \
        run_samtools idxstats "${outfile_pe}"
    run_capture \
        "view submit align-fastqs BWA ALN PE BAM" "${view_pe}" \
        run_samtools view "${outfile_pe}"

    assert_grep_pattern "${idxstats_pe}" $'^I\t108\t2\t0$' \
        "submit BWA ALN PE BAM has two mapped reads on chromosome I"
    assert_grep_pattern "${view_pe}" $'^tiny_pe_pair_1\t' \
        "submit BWA ALN PE BAM contains expected read name"
    assert_grep_pattern "${view_pe}" $'^tiny_pe_pair_1\t99\tI\t17\t' \
        "submit BWA ALN PE BAM has proper-pair R1 flag and start"
    assert_grep_pattern "${view_pe}" $'^tiny_pe_pair_1\t147\tI\t70\t' \
        "submit BWA ALN PE BAM has proper-pair R2 flag and start"
fi

finish
