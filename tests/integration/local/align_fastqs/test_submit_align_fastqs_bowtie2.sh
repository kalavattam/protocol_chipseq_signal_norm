#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_submit_align_fastqs_bowtie2.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in design, development,
# and documentation, with all output reviewed, edited, and approved by the
# author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="submit align-fastqs Bowtie2"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


#  Define fixture and output paths for Bowtie2 alignment
dir_fx="${ROOT_REPO}/tests/fixtures/align_fastqs"
ref_src="${dir_fx}/reference/tiny.fa"
in_se="${dir_fx}/fastq/se/tiny_se.atria.fastq.gz"
in_pe_1="${dir_fx}/fastq/pe/tiny_pe_R1.atria.fastq.gz"
in_pe_2="${dir_fx}/fastq/pe/tiny_pe_R2.atria.fastq.gz"
in_pe="${in_pe_1},${in_pe_2}"
idx_bt2="${dir_fx}/bowtie2/tiny"

tmp="${TEST_DIR_TMP}/submit_align_fastqs_bowtie2"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/align_fastqs"
ref_fa="${tmp}/tiny.fa"

out_se="${dir_out}/tiny_se.bam"
bai_se="${out_se}.bai"
stat_se="${dir_out}/tiny_se.idxstats.txt"
vw_se="${dir_out}/tiny_se.view.txt"

out_pe="${dir_out}/tiny_pe.bam"
bai_pe="${out_pe}.bai"
stat_pe="${dir_out}/tiny_pe.idxstats.txt"
vw_pe="${dir_out}/tiny_pe.view.txt"

out_se_cram="${dir_out}/tiny_se.cram"
crai_se_cram="${out_se_cram}.crai"
count_se_cram="${dir_out}/tiny_se.cram.count.txt"
vw_se_cram="${dir_out}/tiny_se.cram.view.txt"

out_pe_cram="${dir_out}/tiny_pe.cram"
crai_pe_cram="${out_pe_cram}.crai"
count_pe_cram="${dir_out}/tiny_pe.cram.count.txt"
vw_pe_cram="${dir_out}/tiny_pe.cram.view.txt"

log_se="${dir_log}/submit_align_fastqs_bowtie2_se.log"
log_se_qc="${dir_log}/submit_align_fastqs_bowtie2_quickcheck.log"
log_se_cram="${dir_log}/submit_align_fastqs_bowtie2_se_cram.log"
log_se_qc_cram="${dir_log}/submit_align_fastqs_bowtie2_se_cram_quickcheck.log"

log_pe="${dir_log}/submit_align_fastqs_bowtie2_pe.log"
log_pe_qc="${dir_log}/submit_align_fastqs_bowtie2_pe_quickcheck.log"
log_pe_cram="${dir_log}/submit_align_fastqs_bowtie2_pe_cram.log"
log_pe_qc_cram="${dir_log}/submit_align_fastqs_bowtie2_pe_cram_quickcheck.log"


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty \
    "${ref_src}" \
    "${in_se}" \
    "${in_pe_1}" \
    "${in_pe_2}" \
    "${idx_bt2}.1.bt2" \
    "${idx_bt2}.2.bt2" \
    "${idx_bt2}.3.bt2" \
    "${idx_bt2}.4.bt2" \
    "${idx_bt2}.rev.1.bt2" \
    "${idx_bt2}.rev.2.bt2" \
    || {
        finish
        exit $?
    }

cp "${ref_src}" "${ref_fa}"


#  Align one SE FASTQ fixture with Bowtie2 and emit BAM
# shellcheck disable=SC2154
if \
    run_capture \
        "submit align-fastqs bowtie2" \
        "${log_se}" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/submit_align_fastqs.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/bin" \
            --threads 1 \
            --aligner bowtie2 \
            --bt2_mode global \
            --mapq 0 \
            --index "${idx_bt2}" \
            --csv_fil_in "${in_se}" \
            --dir_out "${dir_out}" \
            --out_ext bam \
            --sfx_se ".atria.fastq.gz" \
            --sfx_pe "_R1.atria.fastq.gz" \
            --dir_eo "${dir_err}" \
            --nam_job "test_submit_align_bowtie2"
then
    record_pass "submit_align_fastqs.sh Bowtie2 BAM exits 0"
else
    record_fail \
        "submit_align_fastqs.sh Bowtie2 BAM failed; see" \
        "$(print_relpath "${log_se}")"
fi

assert_file_nonempty \
    "${out_se}" \
    "submit Bowtie2 SE BAM output"

assert_file_nonempty \
    "${bai_se}" \
    "submit Bowtie2 SE BAM index"

if [[ -s "${out_se}" ]]; then
    if \
        run_capture \
            "quickcheck submit align-fastqs Bowtie2 BAM" \
            "${log_se_qc}" \
            run_samtools quickcheck "${out_se}"
    then
        record_pass "submit Bowtie2 SE BAM passes samtools quickcheck"
    else
        record_fail "submit Bowtie2 SE BAM fails samtools quickcheck"
    fi

    run_capture \
        "idxstats submit align-fastqs Bowtie2 SE BAM" \
        "${stat_se}" \
        run_samtools idxstats "${out_se}"

    run_capture \
        "view submit align-fastqs Bowtie2 SE BAM" \
        "${vw_se}" \
        run_samtools view "${out_se}"

    assert_pattern_found \
        "${stat_se}" \
        $'^I\t108\t1\t0$' \
        "submit Bowtie2 SE BAM has one mapped read on chromosome I"

    assert_pattern_found \
        "${vw_se}" \
        $'^tiny_se_read_1\t' \
        "submit Bowtie2 SE BAM contains expected read name"
fi


#  Align one PE FASTQ fixture with Bowtie2, requiring proper pairs
# shellcheck disable=SC2154
if \
    run_capture \
        "submit align-fastqs bowtie2 pe" \
        "${log_pe}" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/submit_align_fastqs.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/bin" \
            --threads 1 \
            --aligner bowtie2 \
            --bt2_mode global \
            --mapq 0 \
            --req_flg \
            --index "${idx_bt2}" \
            --csv_fil_in "${in_pe}" \
            --dir_out "${dir_out}" \
            --out_ext bam \
            --sfx_se ".atria.fastq.gz" \
            --sfx_pe "_R1.atria.fastq.gz" \
            --dir_eo "${dir_err}" \
            --nam_job "test_submit_align_bowtie2_pe"
then
    record_pass "submit_align_fastqs.sh Bowtie2 PE BAM exits 0"
else
    record_fail \
        "submit_align_fastqs.sh Bowtie2 PE BAM failed; see" \
        "$(print_relpath "${log_pe}")"
fi

assert_file_nonempty \
    "${out_pe}" \
    "submit Bowtie2 PE BAM output"

assert_file_nonempty \
    "${bai_pe}" \
    "submit Bowtie2 PE BAM index"

if [[ -s "${out_pe}" ]]; then
    if \
        run_capture \
            "quickcheck submit align-fastqs Bowtie2 PE BAM" \
            "${log_pe_qc}" \
            run_samtools quickcheck "${out_pe}"
    then
        record_pass "submit Bowtie2 PE BAM passes samtools quickcheck"
    else
        record_fail "submit Bowtie2 PE BAM fails samtools quickcheck"
    fi

    run_capture \
        "idxstats submit align-fastqs Bowtie2 PE BAM" \
        "${stat_pe}" \
        run_samtools idxstats "${out_pe}"

    run_capture \
        "view submit align-fastqs Bowtie2 PE BAM" \
        "${vw_pe}" \
        run_samtools view "${out_pe}"

    assert_pattern_found \
        "${stat_pe}" \
        $'^I\t108\t2\t0$' \
        "submit Bowtie2 PE BAM has two mapped reads on chromosome I"

    assert_pattern_found \
        "${vw_pe}" \
        $'^tiny_pe_pair_1\t' \
        "submit Bowtie2 PE BAM contains expected read name"

    assert_pattern_found \
        "${vw_pe}" \
        $'^tiny_pe_pair_1\t99\tI\t17\t' \
        "submit Bowtie2 PE BAM has proper-pair R1 flag and start"

    assert_pattern_found \
        "${vw_pe}" \
        $'^tiny_pe_pair_1\t147\tI\t70\t' \
        "submit Bowtie2 PE BAM has proper-pair R2 flag and start"
fi


#  Align one SE FASTQ fixture with Bowtie2 and emit CRAM
# shellcheck disable=SC2154
if \
    run_capture \
        "submit align-fastqs bowtie2 se cram" \
        "${log_se_cram}" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/submit_align_fastqs.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/bin" \
            --threads 1 \
            --aligner bowtie2 \
            --bt2_mode global \
            --mapq 0 \
            --index "${idx_bt2}" \
            --ref_fa "${ref_fa}" \
            --csv_fil_in "${in_se}" \
            --dir_out "${dir_out}" \
            --out_ext cram \
            --sfx_se ".atria.fastq.gz" \
            --sfx_pe "_R1.atria.fastq.gz" \
            --dir_eo "${dir_err}" \
            --nam_job "test_submit_align_bowtie2_se_cram"
then
    record_pass "submit_align_fastqs.sh Bowtie2 SE CRAM exits 0"
else
    record_fail \
        "submit_align_fastqs.sh Bowtie2 SE CRAM failed; see" \
        "$(print_relpath "${log_se_cram}")"
fi

assert_file_nonempty \
    "${out_se_cram}" \
    "submit Bowtie2 SE CRAM output"

assert_file_nonempty \
    "${crai_se_cram}" \
    "submit Bowtie2 SE CRAI index"

if [[ -s "${out_se_cram}" ]]; then
    if \
        run_capture \
            "quickcheck submit align-fastqs Bowtie2 SE CRAM" \
            "${log_se_qc_cram}" \
            run_samtools quickcheck "${out_se_cram}"
    then
        record_pass "submit Bowtie2 SE CRAM passes samtools quickcheck"
    else
        record_fail "submit Bowtie2 SE CRAM fails samtools quickcheck"
    fi

    run_capture \
        "count submit align-fastqs Bowtie2 SE CRAM" \
        "${count_se_cram}" \
        run_samtools view -T "${ref_fa}" -c "${out_se_cram}" I

    run_capture \
        "view submit align-fastqs Bowtie2 SE CRAM" \
        "${vw_se_cram}" \
        run_samtools view -T "${ref_fa}" "${out_se_cram}"

    assert_pattern_found \
        "${count_se_cram}" \
        $'^1$' \
        "submit Bowtie2 SE CRAM has one mapped read on chromosome I"

    assert_pattern_found \
        "${vw_se_cram}" \
        $'^tiny_se_read_1\t' \
        "submit Bowtie2 SE CRAM contains expected read name"
fi


#  Align one PE FASTQ fixture with Bowtie2 and emit CRAM
# shellcheck disable=SC2154
if \
    run_capture \
        "submit align-fastqs bowtie2 pe cram" \
        "${log_pe_cram}" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/submit_align_fastqs.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/bin" \
            --threads 1 \
            --aligner bowtie2 \
            --bt2_mode global \
            --mapq 0 \
            --req_flg \
            --index "${idx_bt2}" \
            --ref_fa "${ref_fa}" \
            --csv_fil_in "${in_pe}" \
            --dir_out "${dir_out}" \
            --out_ext cram \
            --sfx_se ".atria.fastq.gz" \
            --sfx_pe "_R1.atria.fastq.gz" \
            --dir_eo "${dir_err}" \
            --nam_job "test_submit_align_bowtie2_pe_cram"
then
    record_pass "submit_align_fastqs.sh Bowtie2 PE CRAM exits 0"
else
    record_fail \
        "submit_align_fastqs.sh Bowtie2 PE CRAM failed; see" \
        "$(print_relpath "${log_pe_cram}")"
fi

assert_file_nonempty \
    "${out_pe_cram}" \
    "submit Bowtie2 PE CRAM output"

assert_file_nonempty \
    "${crai_pe_cram}" \
    "submit Bowtie2 PE CRAI index"

if [[ -s "${out_pe_cram}" ]]; then
    if \
        run_capture \
            "quickcheck submit align-fastqs Bowtie2 PE CRAM" \
            "${log_pe_qc_cram}" \
            run_samtools quickcheck "${out_pe_cram}"
    then
        record_pass "submit Bowtie2 PE CRAM passes samtools quickcheck"
    else
        record_fail "submit Bowtie2 PE CRAM fails samtools quickcheck"
    fi

    run_capture \
        "count submit align-fastqs Bowtie2 PE CRAM" \
        "${count_pe_cram}" \
        run_samtools view -T "${ref_fa}" -c "${out_pe_cram}" I

    run_capture \
        "view submit align-fastqs Bowtie2 PE CRAM" \
        "${vw_pe_cram}" \
        run_samtools view -T "${ref_fa}" "${out_pe_cram}"

    assert_pattern_found \
        "${count_pe_cram}" \
        $'^2$' \
        "submit Bowtie2 PE CRAM has two mapped reads on chromosome I"

    assert_pattern_found \
        "${vw_pe_cram}" \
        $'^tiny_pe_pair_1\t' \
        "submit Bowtie2 PE CRAM contains expected read name"

    assert_pattern_found \
        "${vw_pe_cram}" \
        $'^tiny_pe_pair_1\t99\tI\t17\t' \
        "submit Bowtie2 PE CRAM has proper-pair R1 flag and start"

    assert_pattern_found \
        "${vw_pe_cram}" \
        $'^tiny_pe_pair_1\t147\tI\t70\t' \
        "submit Bowtie2 PE CRAM has proper-pair R2 flag and start"
fi

finish
