#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_align_fastqs_slurm.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="execute align-fastqs Slurm"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

print_section "${TEST_NAME}"

#  This smoke test must be run on a Slurm-capable system
#+
#+ Remote submission check:
#+
#+   RUN_SLURM=1 bash tests/scripts/smoke/test_execute_align_fastqs_slurm.sh
#+
#+ Optional output polling after submission:
#+
#+   RUN_SLURM=1 SLURM_WAIT=1 \
#+       bash tests/scripts/smoke/test_execute_align_fastqs_slurm.sh
#+
if [[ "${RUN_SLURM:-0}" != "1" ]]; then
    record_skip \
        "Slurm align-fastqs integration check disabled; set RUN_SLURM=1 on a" \
        "Slurm-capable system to enable"
    finish
    exit $?
fi

#  Define fixture and output paths for a Slurm Bowtie2 submission smoke test
dir_fx="${ROOT_REPO}/tests/align_fastqs/fixtures"
in_se="${dir_fx}/fastq/tiny_se.atria.fastq.gz"
idx="${dir_fx}/bowtie2/tiny"

tmp="${TEST_DIR_TMP}/execute_align_fastqs_slurm"
dir_in="${tmp}/in"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/align_fastqs"

in_1="${dir_in}/tiny_se_1.atria.fastq.gz"
in_2="${dir_in}/tiny_se_2.atria.fastq.gz"
csv_in="${in_1};${in_2}"

out_1="${dir_out}/tiny_se_1.bam"
out_2="${dir_out}/tiny_se_2.bam"
stat_1="${dir_out}/tiny_se_1.idxstats.txt"
stat_2="${dir_out}/tiny_se_2.idxstats.txt"
vw_1="${dir_out}/tiny_se_1.view.txt"
vw_2="${dir_out}/tiny_se_2.view.txt"

rm -rf "${tmp}"
mkdir -p "${dir_in}" "${dir_out}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty \
    "${in_se}" \
    "${idx}.1.bt2" \
    "${idx}.2.bt2" \
    "${idx}.3.bt2" \
    "${idx}.4.bt2" \
    "${idx}.rev.1.bt2" \
    "${idx}.rev.2.bt2" \
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

if \
    check_cmd_exists sbatch
then
    record_pass "sbatch is available for Slurm submission"
else
    record_fail "RUN_SLURM=1 requires sbatch in PATH on a Slurm-capable system"
    finish
    exit $?
fi


#  Submit execute_align_fastqs.sh as a two-task Slurm array
log="${dir_log}/execute_align_fastqs_slurm_bowtie2.log"

if \
    run_capture \
        "execute align-fastqs Slurm Bowtie2 submission" \
        "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/execute_align_fastqs.sh" \
            --threads 1 \
            --aligner bowtie2 \
            --bt2_aln global \
            --mapq 0 \
            --index "${idx}" \
            --csv_infile "${csv_in}" \
            --dir_out "${dir_out}" \
            --out_ext bam \
            --sfx_se ".atria.fastq.gz" \
            --sfx_pe "_R1.atria.fastq.gz" \
            --err_out "${dir_err}" \
            --nam_job "test_execute_align_slurm" \
            --max_job 2 \
            --time "05:00" \
            --slurm
then
    record_pass "execute_align_fastqs.sh Slurm Bowtie2 submission exits 0"
else
    record_fail \
        "execute_align_fastqs.sh Slurm Bowtie2 submission failed; see" \
        "$(print_relpath "${log}")"
fi

assert_file_nonempty \
    "${log}" \
    "execute align-fastqs Slurm submission log"

if [[ -s "${log}" ]]; then
    assert_pattern_found \
        "${log}" \
        '^Submitted batch job [0-9][0-9]*$' \
        "execute align-fastqs Slurm submission reports batch job ID"
fi

if [[ "${SLURM_WAIT:-0}" != "1" ]]; then
    record_skip \
        "Slurm output polling disabled; set SLURM_WAIT=1 to wait for BAM" \
        "products after submission"
    finish
    exit $?
fi


#  Optional, cluster-agnostic output polling
#+
#+ This does not inspect Slurm state; it waits only for expected products to
#+ appear under the test output directory
sec_wait="${SLURM_sec_wait:-120}"
sec_slp="${SLURM_POLL_SECONDS:-5}"
elapsed=0

while (( elapsed < sec_wait )); do
    if \
           [[ -s "${out_1}" ]] \
        && [[ -s "${out_1}.bai" ]] \
        && [[ -s "${out_2}" ]] \
        && [[ -s "${out_2}.bai" ]]
    then
        break
    fi

    sleep "${sec_slp}"
    elapsed=$(( elapsed + sec_slp ))
done

assert_file_nonempty \
    "${out_1}" \
    "execute Slurm Bowtie2 BAM output 1"
assert_file_nonempty \
    "${out_1}.bai" \
    "execute Slurm Bowtie2 BAM index 1"
assert_file_nonempty \
    "${out_2}" \
    "execute Slurm Bowtie2 BAM output 2"
assert_file_nonempty \
    "${out_2}.bai" \
    "execute Slurm Bowtie2 BAM index 2"

if [[ -s "${out_1}" ]]; then
    if \
        run_capture \
            "quickcheck execute align-fastqs Slurm Bowtie2 BAM 1" \
            "${dir_log}/execute_align_fastqs_slurm_bam_1_quickcheck.log" \
            run_samtools quickcheck "${out_1}"
    then
        record_pass "execute Slurm Bowtie2 BAM 1 passes samtools quickcheck"
    else
        record_fail "execute Slurm Bowtie2 BAM 1 fails samtools quickcheck"
    fi

    run_capture \
        "idxstats execute align-fastqs Slurm Bowtie2 BAM 1" \
        "${stat_1}" \
        run_samtools idxstats "${out_1}"
    run_capture \
        "view execute align-fastqs Slurm Bowtie2 BAM 1" \
        "${vw_1}" \
        run_samtools view "${out_1}"

    assert_pattern_found \
        "${stat_1}" \
        $'^I\t108\t1\t0$' \
        "execute Slurm Bowtie2 BAM 1 has one mapped read on chromosome I"
    assert_pattern_found \
        "${vw_1}" \
        $'^tiny_se_read_1\t' \
        "execute Slurm Bowtie2 BAM 1 contains expected read name"
fi

if [[ -s "${out_2}" ]]; then
    if \
        run_capture \
            "quickcheck execute align-fastqs Slurm Bowtie2 BAM 2" \
            "${dir_log}/execute_align_fastqs_slurm_bam_2_quickcheck.log" \
            run_samtools quickcheck "${out_2}"
    then
        record_pass "execute Slurm Bowtie2 BAM 2 passes samtools quickcheck"
    else
        record_fail "execute Slurm Bowtie2 BAM 2 fails samtools quickcheck"
    fi

    run_capture \
        "idxstats execute align-fastqs Slurm Bowtie2 BAM 2" \
        "${stat_2}" \
        run_samtools idxstats "${out_2}"
    run_capture \
        "view execute align-fastqs Slurm Bowtie2 BAM 2" \
        "${vw_2}" \
        run_samtools view "${out_2}"

    assert_pattern_found \
        "${stat_2}" \
        $'^I\t108\t1\t0$' \
        "execute Slurm Bowtie2 BAM 2 has one mapped read on chromosome I"
    assert_pattern_found \
        "${vw_2}" \
        $'^tiny_se_read_1\t' \
        "execute Slurm Bowtie2 BAM 2 contains expected read name"
fi

finish
