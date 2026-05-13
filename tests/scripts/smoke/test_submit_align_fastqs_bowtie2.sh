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

#  Define fixture and output paths for Bowtie2 alignment
dir_fx="${ROOT_REPO}/tests/align_fastqs/fixtures"
ref_src="${dir_fx}/reference/tiny.fa"
infile_se="${dir_fx}/fastq/tiny_se.atria.fastq.gz"
infile_pe_1="${dir_fx}/fastq/tiny_pe_R1.atria.fastq.gz"
infile_pe_2="${dir_fx}/fastq/tiny_pe_R2.atria.fastq.gz"
infile_pe="${infile_pe_1},${infile_pe_2}"
index="${dir_fx}/bowtie2/tiny"

tmp="${TEST_DIR_TMP}/submit_align_fastqs_bowtie2"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/align_fastqs"
ref_fa="${tmp}/tiny.fa"

outfile_se="${dir_out}/tiny_se.bam"
idxstats_se="${dir_out}/tiny_se.idxstats.txt"
view_se="${dir_out}/tiny_se.view.txt"
outfile_pe="${dir_out}/tiny_pe.bam"
idxstats_pe="${dir_out}/tiny_pe.idxstats.txt"
view_pe="${dir_out}/tiny_pe.view.txt"
outfile_se_cram="${dir_out}/tiny_se.cram"
count_se_cram="${dir_out}/tiny_se.cram.count.txt"
view_se_cram="${dir_out}/tiny_se.cram.view.txt"
outfile_pe_cram="${dir_out}/tiny_pe.cram"
count_pe_cram="${dir_out}/tiny_pe.cram.count.txt"
view_pe_cram="${dir_out}/tiny_pe.cram.view.txt"

rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_exist \
    "${ref_src}" \
    "${infile_se}" \
    "${infile_pe_1}" \
    "${infile_pe_2}" \
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

cp "${ref_src}" "${ref_fa}"


#  Run samtools from the active shell or the resolved project environment
function run_samtools() {
    if check_cmd_exists samtools; then
        samtools "$@"
    else
        conda run -n "${env_nam}" samtools "$@"
    fi
}


#  Align one SE FASTQ fixture with Bowtie2 and emit BAM
log="${dir_log}/submit_align_fastqs_bowtie2_se.log"

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
            --mapq 0 \
            --index "${index}" \
            --csv_infile "${infile_se}" \
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

assert_file_nonempty "${outfile_se}" "submit Bowtie2 SE BAM output"
assert_file_nonempty "${outfile_se}.bai" "submit Bowtie2 SE BAM index"

if [[ -s "${outfile_se}" ]]; then
    if run_capture \
        "quickcheck submit align-fastqs Bowtie2 BAM" \
        "${dir_log}/submit_align_fastqs_bowtie2_quickcheck.log" \
        run_samtools quickcheck "${outfile_se}"
    then
        rec_pass "submit Bowtie2 SE BAM passes samtools quickcheck"
    else
        rec_fail "submit Bowtie2 SE BAM fails samtools quickcheck"
    fi

    run_capture \
        "idxstats submit align-fastqs Bowtie2 SE BAM" "${idxstats_se}" \
        run_samtools idxstats "${outfile_se}"
    run_capture \
        "view submit align-fastqs Bowtie2 SE BAM" "${view_se}" \
        run_samtools view "${outfile_se}"

    assert_grep_pattern "${idxstats_se}" $'^I\t108\t1\t0$' \
        "submit Bowtie2 SE BAM has one mapped read on chromosome I"
    assert_grep_pattern "${view_se}" $'^tiny_se_read_1\t' \
        "submit Bowtie2 SE BAM contains expected read name"
fi


#  Align one PE FASTQ fixture with Bowtie2, requiring proper pairs
log="${dir_log}/submit_align_fastqs_bowtie2_pe.log"

# shellcheck disable=SC2154
if \
    run_capture \
        "submit align-fastqs bowtie2 pe" "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/submit_align_fastqs.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/scripts" \
            --threads 1 \
            --aligner bowtie2 \
            --bt2_aln global \
            --mapq 0 \
            --req_flg \
            --index "${index}" \
            --csv_infile "${infile_pe}" \
            --dir_out "${dir_out}" \
            --out_ext bam \
            --sfx_se ".atria.fastq.gz" \
            --sfx_pe "_R1.atria.fastq.gz" \
            --err_out "${dir_err}" \
            --nam_job "test_submit_align_bowtie2_pe"
then
    rec_pass "submit_align_fastqs.sh Bowtie2 PE BAM exits 0"
else
    rec_fail \
        "submit_align_fastqs.sh Bowtie2 PE BAM failed; see" \
        "$(rec_relpath "${log}")"
fi

assert_file_nonempty "${outfile_pe}" "submit Bowtie2 PE BAM output"
assert_file_nonempty "${outfile_pe}.bai" "submit Bowtie2 PE BAM index"

if [[ -s "${outfile_pe}" ]]; then
    if run_capture \
        "quickcheck submit align-fastqs Bowtie2 PE BAM" \
        "${dir_log}/submit_align_fastqs_bowtie2_pe_quickcheck.log" \
        run_samtools quickcheck "${outfile_pe}"
    then
        rec_pass "submit Bowtie2 PE BAM passes samtools quickcheck"
    else
        rec_fail "submit Bowtie2 PE BAM fails samtools quickcheck"
    fi

    run_capture \
        "idxstats submit align-fastqs Bowtie2 PE BAM" "${idxstats_pe}" \
        run_samtools idxstats "${outfile_pe}"
    run_capture \
        "view submit align-fastqs Bowtie2 PE BAM" "${view_pe}" \
        run_samtools view "${outfile_pe}"

    assert_grep_pattern "${idxstats_pe}" $'^I\t108\t2\t0$' \
        "submit Bowtie2 PE BAM has two mapped reads on chromosome I"
    assert_grep_pattern "${view_pe}" $'^tiny_pe_pair_1\t' \
        "submit Bowtie2 PE BAM contains expected read name"
    assert_grep_pattern "${view_pe}" $'^tiny_pe_pair_1\t99\tI\t17\t' \
        "submit Bowtie2 PE BAM has proper-pair R1 flag and start"
    assert_grep_pattern "${view_pe}" $'^tiny_pe_pair_1\t147\tI\t70\t' \
        "submit Bowtie2 PE BAM has proper-pair R2 flag and start"
fi


#  Align one SE FASTQ fixture with Bowtie2 and emit CRAM
log="${dir_log}/submit_align_fastqs_bowtie2_se_cram.log"

# shellcheck disable=SC2154
if \
    run_capture \
        "submit align-fastqs bowtie2 se cram" "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/submit_align_fastqs.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/scripts" \
            --threads 1 \
            --aligner bowtie2 \
            --bt2_aln global \
            --mapq 0 \
            --index "${index}" \
            --ref "${ref_fa}" \
            --csv_infile "${infile_se}" \
            --dir_out "${dir_out}" \
            --out_ext cram \
            --sfx_se ".atria.fastq.gz" \
            --sfx_pe "_R1.atria.fastq.gz" \
            --err_out "${dir_err}" \
            --nam_job "test_submit_align_bowtie2_se_cram"
then
    rec_pass "submit_align_fastqs.sh Bowtie2 SE CRAM exits 0"
else
    rec_fail \
        "submit_align_fastqs.sh Bowtie2 SE CRAM failed; see" \
        "$(rec_relpath "${log}")"
fi

assert_file_nonempty "${outfile_se_cram}" "submit Bowtie2 SE CRAM output"
assert_file_nonempty "${outfile_se_cram}.crai" "submit Bowtie2 SE CRAI index"

if [[ -s "${outfile_se_cram}" ]]; then
    if run_capture \
        "quickcheck submit align-fastqs Bowtie2 SE CRAM" \
        "${dir_log}/submit_align_fastqs_bowtie2_se_cram_quickcheck.log" \
        run_samtools quickcheck "${outfile_se_cram}"
    then
        rec_pass "submit Bowtie2 SE CRAM passes samtools quickcheck"
    else
        rec_fail "submit Bowtie2 SE CRAM fails samtools quickcheck"
    fi

    run_capture \
        "count submit align-fastqs Bowtie2 SE CRAM" "${count_se_cram}" \
        run_samtools view -T "${ref_fa}" -c "${outfile_se_cram}" I
    run_capture \
        "view submit align-fastqs Bowtie2 SE CRAM" "${view_se_cram}" \
        run_samtools view -T "${ref_fa}" "${outfile_se_cram}"

    assert_grep_pattern "${count_se_cram}" $'^1$' \
        "submit Bowtie2 SE CRAM has one mapped read on chromosome I"
    assert_grep_pattern "${view_se_cram}" $'^tiny_se_read_1\t' \
        "submit Bowtie2 SE CRAM contains expected read name"
fi


#  Align one PE FASTQ fixture with Bowtie2 and emit CRAM
log="${dir_log}/submit_align_fastqs_bowtie2_pe_cram.log"

# shellcheck disable=SC2154
if \
    run_capture \
        "submit align-fastqs bowtie2 pe cram" "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/submit_align_fastqs.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/scripts" \
            --threads 1 \
            --aligner bowtie2 \
            --bt2_aln global \
            --mapq 0 \
            --req_flg \
            --index "${index}" \
            --ref "${ref_fa}" \
            --csv_infile "${infile_pe}" \
            --dir_out "${dir_out}" \
            --out_ext cram \
            --sfx_se ".atria.fastq.gz" \
            --sfx_pe "_R1.atria.fastq.gz" \
            --err_out "${dir_err}" \
            --nam_job "test_submit_align_bowtie2_pe_cram"
then
    rec_pass "submit_align_fastqs.sh Bowtie2 PE CRAM exits 0"
else
    rec_fail \
        "submit_align_fastqs.sh Bowtie2 PE CRAM failed; see" \
        "$(rec_relpath "${log}")"
fi

assert_file_nonempty "${outfile_pe_cram}" "submit Bowtie2 PE CRAM output"
assert_file_nonempty "${outfile_pe_cram}.crai" "submit Bowtie2 PE CRAI index"

if [[ -s "${outfile_pe_cram}" ]]; then
    if run_capture \
        "quickcheck submit align-fastqs Bowtie2 PE CRAM" \
        "${dir_log}/submit_align_fastqs_bowtie2_pe_cram_quickcheck.log" \
        run_samtools quickcheck "${outfile_pe_cram}"
    then
        rec_pass "submit Bowtie2 PE CRAM passes samtools quickcheck"
    else
        rec_fail "submit Bowtie2 PE CRAM fails samtools quickcheck"
    fi

    run_capture \
        "count submit align-fastqs Bowtie2 PE CRAM" "${count_pe_cram}" \
        run_samtools view -T "${ref_fa}" -c "${outfile_pe_cram}" I
    run_capture \
        "view submit align-fastqs Bowtie2 PE CRAM" "${view_pe_cram}" \
        run_samtools view -T "${ref_fa}" "${outfile_pe_cram}"

    assert_grep_pattern "${count_pe_cram}" $'^2$' \
        "submit Bowtie2 PE CRAM has two mapped reads on chromosome I"
    assert_grep_pattern "${view_pe_cram}" $'^tiny_pe_pair_1\t' \
        "submit Bowtie2 PE CRAM contains expected read name"
    assert_grep_pattern "${view_pe_cram}" $'^tiny_pe_pair_1\t99\tI\t17\t' \
        "submit Bowtie2 PE CRAM has proper-pair R1 flag and start"
    assert_grep_pattern "${view_pe_cram}" $'^tiny_pe_pair_1\t147\tI\t70\t' \
        "submit Bowtie2 PE CRAM has proper-pair R2 flag and start"
fi

finish
