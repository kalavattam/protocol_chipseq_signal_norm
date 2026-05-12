#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_submit_align_fastqs_bowtie2.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="submit align-fastqs Bowtie2"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

rec_section "${TEST_NAME}"

#  Define fixture and output paths for single-end Bowtie2 alignment
dir_fx="${ROOT_REPO}/tests/align_fastqs/fixtures"
infile="${dir_fx}/fastq/tiny_se.atria.fastq.gz"
index="${dir_fx}/bowtie2/tiny"

tmp="${TEST_DIR_TMP}/submit_align_fastqs_bowtie2"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/align_fastqs"

outfile="${dir_out}/tiny_se.bam"
idxstats="${dir_out}/tiny_se.idxstats.txt"
view_out="${dir_out}/tiny_se.view.txt"

rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_exist \
    "${infile}" \
    "${index}.1.bt2" \
    "${index}.2.bt2" \
    "${index}.3.bt2" \
    "${index}.4.bt2" \
    "${index}.rev.1.bt2" \
    "${index}.rev.2.bt2" \
    || {
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


#  Align one SE FASTQ fixture with Bowtie2 and emit BAM
log="${dir_log}/submit_align_fastqs_bowtie2.log"

# shellcheck disable=SC2154
if \
    run_capture \
        "submit align-fastqs bowtie2" "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/submit_align_fastqs.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/scripts" \
            --threads 1 \
            --aligner bowtie2 \
            --bt2_aln global \
            --bwa_alg mem \
            --mapq 0 \
            --index "${index}" \
            --csv_infile "${infile}" \
            --dir_out "${dir_out}" \
            --out_ext bam \
            --sfx_se ".atria.fastq.gz" \
            --sfx_pe "_R1.atria.fastq.gz" \
            --err_out "${dir_err}" \
            --nam_job "test_submit_align_bowtie2"
then
    rec_pass "submit_align_fastqs.sh Bowtie2 BAM exits 0"
else
    rec_fail \
        "submit_align_fastqs.sh Bowtie2 BAM failed; see" \
        "$(rec_relpath "${log}")"
fi

assert_file_nonempty "${outfile}" "submit Bowtie2 BAM output"
assert_file_nonempty "${outfile}.bai" "submit Bowtie2 BAM index"

if [[ -s "${outfile}" ]]; then
    if run_capture \
        "quickcheck submit align-fastqs Bowtie2 BAM" \
        "${dir_log}/submit_align_fastqs_bowtie2_quickcheck.log" \
        run_samtools quickcheck "${outfile}"
    then
        rec_pass "submit Bowtie2 BAM passes samtools quickcheck"
    else
        rec_fail "submit Bowtie2 BAM fails samtools quickcheck"
    fi

    run_capture \
        "idxstats submit align-fastqs Bowtie2 BAM" "${idxstats}" \
        run_samtools idxstats "${outfile}"
    run_capture \
        "view submit align-fastqs Bowtie2 BAM" "${view_out}" \
        run_samtools view "${outfile}"

    assert_grep_pattern "${idxstats}" $'^I\t108\t1\t0$' \
        "submit Bowtie2 BAM has one mapped read on chromosome I"
    assert_grep_pattern "${view_out}" $'^tiny_se_read_1\t' \
        "submit Bowtie2 BAM contains expected read name"
fi

finish
