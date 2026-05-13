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

rec_section "${TEST_NAME}"

#  This smoke test must be run on a Slurm-capable system
#+
#+ Remote submission check:
#+
#+   RUN_SLURM=1 bash tests/scripts/smoke/test_execute_align_fastqs_slurm.sh
#+
#+ Optional output polling after submission:
#+
#+   RUN_SLURM=1 SLURM_WAIT=1 bash tests/scripts/smoke/test_execute_align_fastqs_slurm.sh
#+
if [[ "${RUN_SLURM:-0}" != "1" ]]; then
    rec_skip \
        "Slurm align-fastqs integration check disabled;" \
        "set RUN_SLURM=1 on a Slurm-capable system to enable"
    finish
    exit $?
fi

#  Define fixture and output paths for a Slurm Bowtie2 submission smoke test
dir_fx="${ROOT_REPO}/tests/align_fastqs/fixtures"
infile_se="${dir_fx}/fastq/tiny_se.atria.fastq.gz"
index="${dir_fx}/bowtie2/tiny"

tmp="${TEST_DIR_TMP}/execute_align_fastqs_slurm"
dir_in="${tmp}/in"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/align_fastqs"

infile_1="${dir_in}/tiny_se_1.atria.fastq.gz"
infile_2="${dir_in}/tiny_se_2.atria.fastq.gz"
csv_infile="${infile_1};${infile_2}"

outfile_1="${dir_out}/tiny_se_1.bam"
outfile_2="${dir_out}/tiny_se_2.bam"
idxstats_1="${dir_out}/tiny_se_1.idxstats.txt"
idxstats_2="${dir_out}/tiny_se_2.idxstats.txt"
view_1="${dir_out}/tiny_se_1.view.txt"
view_2="${dir_out}/tiny_se_2.view.txt"

rm -rf "${tmp}"
mkdir -p "${dir_in}" "${dir_out}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_exist \
    "${infile_se}" \
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

cp "${infile_se}" "${infile_1}"
cp "${infile_se}" "${infile_2}"

require_files_exist "${infile_1}" "${infile_2}" || {
    finish
    exit $?
}

if check_cmd_exists sbatch; then
    rec_pass "sbatch is available for Slurm submission"
else
    rec_fail "RUN_SLURM=1 requires sbatch in PATH on a Slurm-capable system"
    finish
    exit $?
fi


#  Run samtools from the active shell or the resolved project environment
# shellcheck disable=SC2154
function run_samtools() {
    if check_cmd_exists samtools; then
        samtools "$@"
    else
        conda run -n "${env_nam}" samtools "$@"
    fi
}


#  Submit execute_align_fastqs.sh as a two-task Slurm array
log="${dir_log}/execute_align_fastqs_slurm_bowtie2.log"

if \
    run_capture \
        "execute align-fastqs Slurm Bowtie2 submission" "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/execute_align_fastqs.sh" \
            --threads 1 \
            --aligner bowtie2 \
            --bt2_aln global \
            --mapq 0 \
            --index "${index}" \
            --csv_infile "${csv_infile}" \
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
    rec_pass "execute_align_fastqs.sh Slurm Bowtie2 submission exits 0"
else
    rec_fail \
        "execute_align_fastqs.sh Slurm Bowtie2 submission failed; see" \
        "$(rec_relpath "${log}")"
fi

assert_file_nonempty "${log}" "execute align-fastqs Slurm submission log"

if [[ -s "${log}" ]]; then
    assert_grep_pattern "${log}" '^Submitted batch job [0-9][0-9]*$' \
        "execute align-fastqs Slurm submission reports batch job ID"
fi

if [[ "${SLURM_WAIT:-0}" != "1" ]]; then
    rec_skip \
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
           [[ -s "${outfile_1}" ]] \
        && [[ -s "${outfile_1}.bai" ]] \
        && [[ -s "${outfile_2}" ]] \
        && [[ -s "${outfile_2}.bai" ]]
    then
        break
    fi

    sleep "${sec_slp}"
    elapsed=$(( elapsed + sec_slp ))
done

assert_file_nonempty "${outfile_1}" "execute Slurm Bowtie2 BAM output 1"
assert_file_nonempty "${outfile_1}.bai" "execute Slurm Bowtie2 BAM index 1"
assert_file_nonempty "${outfile_2}" "execute Slurm Bowtie2 BAM output 2"
assert_file_nonempty "${outfile_2}.bai" "execute Slurm Bowtie2 BAM index 2"

if [[ -s "${outfile_1}" ]]; then
    if run_capture \
        "quickcheck execute align-fastqs Slurm Bowtie2 BAM 1" \
        "${dir_log}/execute_align_fastqs_slurm_bam_1_quickcheck.log" \
        run_samtools quickcheck "${outfile_1}"
    then
        rec_pass "execute Slurm Bowtie2 BAM 1 passes samtools quickcheck"
    else
        rec_fail "execute Slurm Bowtie2 BAM 1 fails samtools quickcheck"
    fi

    run_capture \
        "idxstats execute align-fastqs Slurm Bowtie2 BAM 1" "${idxstats_1}" \
        run_samtools idxstats "${outfile_1}"
    run_capture \
        "view execute align-fastqs Slurm Bowtie2 BAM 1" "${view_1}" \
        run_samtools view "${outfile_1}"

    assert_grep_pattern "${idxstats_1}" $'^I\t108\t1\t0$' \
        "execute Slurm Bowtie2 BAM 1 has one mapped read on chromosome I"
    assert_grep_pattern "${view_1}" $'^tiny_se_read_1\t' \
        "execute Slurm Bowtie2 BAM 1 contains expected read name"
fi

if [[ -s "${outfile_2}" ]]; then
    if run_capture \
        "quickcheck execute align-fastqs Slurm Bowtie2 BAM 2" \
        "${dir_log}/execute_align_fastqs_slurm_bam_2_quickcheck.log" \
        run_samtools quickcheck "${outfile_2}"
    then
        rec_pass "execute Slurm Bowtie2 BAM 2 passes samtools quickcheck"
    else
        rec_fail "execute Slurm Bowtie2 BAM 2 fails samtools quickcheck"
    fi

    run_capture \
        "idxstats execute align-fastqs Slurm Bowtie2 BAM 2" "${idxstats_2}" \
        run_samtools idxstats "${outfile_2}"
    run_capture \
        "view execute align-fastqs Slurm Bowtie2 BAM 2" "${view_2}" \
        run_samtools view "${outfile_2}"

    assert_grep_pattern "${idxstats_2}" $'^I\t108\t1\t0$' \
        "execute Slurm Bowtie2 BAM 2 has one mapped read on chromosome I"
    assert_grep_pattern "${view_2}" $'^tiny_se_read_1\t' \
        "execute Slurm Bowtie2 BAM 2 contains expected read name"
fi

finish
