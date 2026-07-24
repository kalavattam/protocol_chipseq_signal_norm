#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_trim_fastqs_pe.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="execute trim-fastqs PE"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


if ! \
    is_atria_enabled
then
    record_skip \
        "Atria execute trim-fastqs PE check disabled; set RUN_ATRIA=1 to" \
        "enable"
    finish
    exit $?
fi


#  Define fixture and output paths for local serial execute-layer Atria run
dir_fx="${ROOT_REPO}/tests/fixtures/trim_fastqs"
in_r1="${dir_fx}/fastq/pe/tiny_pe_R1.fastq.gz"
in_r2="${dir_fx}/fastq/pe/tiny_pe_R2.fastq.gz"
csv_in="${in_r1},${in_r2}"

tmp="${TEST_DIR_TMP}/execute_trim_fastqs_pe"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/trim_fastqs"

vw_r1="${dir_out}/tiny_pe_R1.trimmed.fastq"
vw_r2="${dir_out}/tiny_pe_R2.trimmed.fastq"
cnt_r1="${dir_out}/tiny_pe_R1.read_count.txt"
cnt_r2="${dir_out}/tiny_pe_R2.read_count.txt"

log_env="${dir_log}/execute_trim_fastqs_pe_env.log"
log_run="${dir_log}/execute_trim_fastqs_pe.log"
log_out_ser="${dir_err}/test_execute_trim_pe_ser.stdout.txt"
log_err_ser="${dir_err}/test_execute_trim_pe_ser.stderr.txt"
log_out_pe="${dir_err}/test_execute_trim_pe.tiny_pe.stdout.txt"
log_err_pe="${dir_err}/test_execute_trim_pe.tiny_pe.stderr.txt"


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_err}" "${dir_log}"

require_files_nonempty \
    "${in_r1}" \
    "${in_r2}" || {
    finish
    exit $?
}

require_env_project env_nam || {
    finish
    exit $?
}

if ! \
    require_env_atria \
        "${env_nam}" \
        "${log_env}"
then
    finish
    exit $?
fi


#  Run execute_trim_fastqs.sh serially on one PE FASTQ pair
if \
    run_capture \
        "execute trim-fastqs PE Atria wet run" \
        "${log_run}" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/execute_trim_fastqs.sh" \
            --env_nam "${env_nam}" \
            --threads 1 \
            --csv_fil_in "${csv_in}" \
            --dir_out "${dir_out}" \
            --sfx_se ".fastq.gz" \
            --sfx_pe "_R1.fastq.gz" \
            --dir_eo "${dir_err}" \
            --nam_job "test_execute_trim_pe" \
            --max_job 1
then
    record_pass "execute_trim_fastqs.sh PE Atria wet run exits 0"
else
    record_fail \
        "execute_trim_fastqs.sh PE Atria wet run failed; see" \
        "$(print_relpath "${log_run}")"
fi

# shellcheck disable=SC2034
{
    mapfile -t trimmed_outputs_r1 < <(
        find "${dir_out}" -maxdepth 1 -type f \
            \( -name 'tiny_pe*R1*.fastq.gz' -o -name 'tiny_pe*R1*.fq.gz' \) \
            | sort
    )
    mapfile -t trimmed_outputs_r2 < <(
        find "${dir_out}" -maxdepth 1 -type f \
            \( -name 'tiny_pe*R2*.fastq.gz' -o -name 'tiny_pe*R2*.fq.gz' \) \
            | sort
    )
}

assert_path_found \
    trimmed_outputs_r1 \
    "compressed PE R1 trimmed FASTQ output" \
    "${dir_out}" \
    out_r1

assert_path_found \
    trimmed_outputs_r2 \
    "compressed PE R2 trimmed FASTQ output" \
    "${dir_out}" \
    out_r2

# shellcheck disable=SC2154
{
    assert_fastq_gzip \
        "${out_r1}" \
        '^@tiny_trim_pe_pair_1' \
        "${vw_r1}" \
        "${cnt_r1}" \
        "" \
        "execute trim-fastqs PE R1 FASTQ output"

    assert_fastq_gzip \
        "${out_r2}" \
        '^@tiny_trim_pe_pair_1' \
        "${vw_r2}" \
        "${cnt_r2}" \
        "" \
        "execute trim-fastqs PE R2 FASTQ output"
}

assert_file_exists \
    "${log_out_ser}" \
    "execute trim-fastqs PE serial stdout log exists"

assert_file_exists \
    "${log_err_ser}" \
    "execute trim-fastqs PE serial stderr log exists"

assert_file_exists \
    "${log_out_pe}" \
    "execute trim-fastqs PE submit stdout log exists"

assert_file_exists \
    "${log_err_pe}" \
    "execute trim-fastqs PE submit stderr log exists"

finish
