#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_align_fastqs_parallel.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="execute align-fastqs GNU Parallel"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

print_section "${TEST_NAME}"

if ! \
    is_parallel_enabled
then
    record_skip \
        "GNU Parallel align-fastqs check disabled; set RUN_PARALLEL=1 to" \
        "enable"
    finish
    exit $?
fi

#  Define fixture and output paths for a local GNU Parallel Bowtie2 wet run
dir_fx="${ROOT_REPO}/tests/align_fastqs/fixtures"
in_se="${dir_fx}/fastq/se/tiny_se.atria.fastq.gz"
idx_bt2="${dir_fx}/bowtie2/tiny"

tmp="${TEST_DIR_TMP}/execute_align_fastqs_parallel"
dir_in="${tmp}/in"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/align_fastqs"

in_1="${dir_in}/tiny_se_1.atria.fastq.gz"
in_2="${dir_in}/tiny_se_2.atria.fastq.gz"
csv_in="${in_1};${in_2}"

out_1="${dir_out}/tiny_se_1.bam"
out_2="${dir_out}/tiny_se_2.bam"
bai_1="${out_1}.bai"
bai_2="${out_2}.bai"
stat_1="${dir_out}/tiny_se_1.idxstats.txt"
stat_2="${dir_out}/tiny_se_2.idxstats.txt"
vw_1="${dir_out}/tiny_se_1.view.txt"
vw_2="${dir_out}/tiny_se_2.view.txt"
cfg="${dir_err}/test_execute_align_parallel.config_parallel.txt"

log_env_parallel="${dir_log}/execute_align_fastqs_parallel_env.log"
log_run="${dir_log}/execute_align_fastqs_parallel_bowtie2.log"
log_qc_1="${dir_log}/execute_align_fastqs_parallel_bam_1_quickcheck.log"
log_qc_2="${dir_log}/execute_align_fastqs_parallel_bam_2_quickcheck.log"

rm -rf "${tmp}"
mkdir -p "${dir_in}" "${dir_out}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty \
    "${in_se}" \
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

cp "${in_se}" "${in_1}"
cp "${in_se}" "${in_2}"

require_files_nonempty \
    "${in_1}" \
    "${in_2}" || {
    finish
    exit $?
}

# shellcheck disable=SC2154
if ! \
    require_env_parallel \
        "${env_nam}" \
        "${log_env_parallel}"
then
    finish
    exit $?
fi


#  Run execute_align_fastqs.sh through local GNU Parallel for two SE inputs
if \
    run_capture \
        "execute align-fastqs GNU Parallel Bowtie2 wet run" \
        "${log_run}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/execute_align_fastqs.sh" \
            --threads 2 \
            --aligner bowtie2 \
            --bt2_aln global \
            --mapq 0 \
            --index "${idx_bt2}" \
            --csv_infile "${csv_in}" \
            --dir_out "${dir_out}" \
            --out_ext bam \
            --sfx_se ".atria.fastq.gz" \
            --sfx_pe "_R1.atria.fastq.gz" \
            --err_out "${dir_err}" \
            --nam_job "test_execute_align_parallel" \
            --max_job 2
then
    record_pass "execute_align_fastqs.sh GNU Parallel Bowtie2 wet run exits 0"
else
    record_fail \
        "execute_align_fastqs.sh GNU Parallel Bowtie2 wet run failed; see" \
        "$(print_relpath "${log_run}")"
fi

assert_file_nonempty \
    "${cfg}" \
    "execute align-fastqs GNU Parallel config"

assert_file_nonempty \
    "${out_1}" \
    "execute GNU Parallel Bowtie2 BAM output 1"

assert_file_nonempty \
    "${bai_1}" \
    "execute GNU Parallel Bowtie2 BAM index 1"

assert_file_nonempty \
    "${out_2}" \
    "execute GNU Parallel Bowtie2 BAM output 2"

assert_file_nonempty \
    "${bai_2}" \
    "execute GNU Parallel Bowtie2 BAM index 2"

if [[ -s "${cfg}" ]]; then
    assert_pattern_found \
        "${cfg}" \
        "${TEST_BASH} ${ROOT_REPO}/scripts/submit_align_fastqs.sh" \
        "execute align-fastqs GNU Parallel config uses Bash-prefixed submit" \
        "command"
fi

if [[ -s "${out_1}" ]]; then
    if \
        run_capture \
            "quickcheck execute align-fastqs GNU Parallel Bowtie2 BAM 1" \
            "${log_qc_1}" \
            run_samtools quickcheck "${out_1}"
    then
        record_pass \
            "execute GNU Parallel Bowtie2 BAM 1 passes samtools quickcheck"
    else
        record_fail \
            "execute GNU Parallel Bowtie2 BAM 1 fails samtools quickcheck"
    fi

    run_capture \
        "idxstats execute align-fastqs GNU Parallel Bowtie2 BAM 1" \
        "${stat_1}" \
        run_samtools idxstats "${out_1}"

    run_capture \
        "view execute align-fastqs GNU Parallel Bowtie2 BAM 1" \
        "${vw_1}" \
        run_samtools view "${out_1}"

    assert_pattern_found \
        "${stat_1}" \
        $'^I\t108\t1\t0$' \
        "execute GNU Parallel Bowtie2 BAM 1 has one mapped read on" \
        "chromosome I"

    assert_pattern_found \
        "${vw_1}" \
        $'^tiny_se_read_1\t' \
        "execute GNU Parallel Bowtie2 BAM 1 contains expected read name"
fi

if [[ -s "${out_2}" ]]; then
    if \
        run_capture \
            "quickcheck execute align-fastqs GNU Parallel Bowtie2 BAM 2" \
            "${log_qc_2}" \
            run_samtools quickcheck "${out_2}"
    then
        record_pass \
            "execute GNU Parallel Bowtie2 BAM 2 passes samtools quickcheck"
    else
        record_fail \
            "execute GNU Parallel Bowtie2 BAM 2 fails samtools quickcheck"
    fi

    run_capture \
        "idxstats execute align-fastqs GNU Parallel Bowtie2 BAM 2" \
        "${stat_2}" \
        run_samtools idxstats "${out_2}"

    run_capture \
        "view execute align-fastqs GNU Parallel Bowtie2 BAM 2" \
        "${vw_2}" \
        run_samtools view "${out_2}"

    assert_pattern_found \
        "${stat_2}" \
        $'^I\t108\t1\t0$' \
        "execute GNU Parallel Bowtie2 BAM 2 has one mapped read on" \
        "chromosome I"

    assert_pattern_found \
        "${vw_2}" \
        $'^tiny_se_read_1\t' \
        "execute GNU Parallel Bowtie2 BAM 2 contains expected read name"
fi

finish
