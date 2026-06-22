#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_trim_fastqs_se.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="execute trim-fastqs SE"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

print_section "${TEST_NAME}"

if ! \
    is_atria_enabled
then
    record_skip \
        "Atria execute trim-fastqs SE check disabled; set RUN_ATRIA=1 to" \
        "enable"
    finish
    exit $?
fi


#  Define fixture and output paths for local serial execute-layer Atria run
dir_fx="${ROOT_REPO}/tests/trim_fastqs/fixtures"
in_se="${dir_fx}/fastq/se/tiny_se.fastq.gz"

tmp="${TEST_DIR_TMP}/execute_trim_fastqs_se"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/trim_fastqs"

vw_fq="${dir_out}/tiny_se.trimmed.fastq"
cnt="${dir_out}/tiny_se.read_count.txt"

log_env="${dir_log}/execute_trim_fastqs_se_env.log"
log_run="${dir_log}/execute_trim_fastqs_se.log"
log_out_ser="${dir_err}/test_execute_trim_se_ser.stdout.txt"
log_err_ser="${dir_err}/test_execute_trim_se_ser.stderr.txt"
log_out_se="${dir_err}/test_execute_trim_se.tiny_se.stdout.txt"
log_err_se="${dir_err}/test_execute_trim_se.tiny_se.stderr.txt"

rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_err}" "${dir_log}"

require_files_nonempty \
    "${in_se}" || {
    finish
    exit $?
}

if ! \
    require_env_atria "env_protocol" "${log_env}"
then
    finish
    exit $?
fi


#  Run execute_trim_fastqs.sh serially on one SE FASTQ fixture
if \
    run_capture \
        "execute trim-fastqs SE Atria wet run" \
        "${log_run}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/execute_trim_fastqs.sh" \
            --threads 1 \
            --csv_infile "${in_se}" \
            --dir_out "${dir_out}" \
            --sfx_se ".fastq.gz" \
            --sfx_pe "_R1.fastq.gz" \
            --err_out "${dir_err}" \
            --nam_job "test_execute_trim_se" \
            --max_job 1
then
    record_pass "execute_trim_fastqs.sh SE Atria wet run exits 0"
else
    record_fail \
        "execute_trim_fastqs.sh SE Atria wet run failed; see" \
        "$(print_relpath "${log_run}")"
fi

# shellcheck disable=2034
mapfile -t trimmed_outputs < <(
    find "${dir_out}" -maxdepth 1 -type f \
        \( -name 'tiny_se*.fastq.gz' -o -name 'tiny_se*.fq.gz' \) \
        | sort
)

assert_path_found \
    trimmed_outputs \
    "compressed SE trimmed FASTQ output" \
    "${dir_out}" \
    fil_out

if [[ -n "${fil_out}" ]]; then
    assert_file_nonempty \
        "${fil_out}" \
        "execute trim-fastqs SE FASTQ output"

    if \
        gzip -t "${fil_out}"
    then
        record_pass "execute trim-fastqs SE FASTQ output passes gzip integrity"
    else
        record_fail "execute trim-fastqs SE FASTQ output fails gzip integrity"
    fi

    if \
        gzip -cd "${fil_out}" > "${vw_fq}"
    then
        record_pass "execute trim-fastqs SE FASTQ output can be decompressed"
    else
        record_fail "execute trim-fastqs SE FASTQ output cannot be decompressed"
    fi

    if [[ -s "${vw_fq}" ]]; then
        assert_pattern_found \
            "${vw_fq}" \
            '^@tiny_trim_se_read_1$' \
            "execute trim-fastqs SE FASTQ contains expected read name"

        awk 'NR % 4 == 1 { n++ } END { print n + 0 }' \
            "${vw_fq}" > "${cnt}"
        assert_pattern_found \
            "${cnt}" \
            '^1$' \
            "execute trim-fastqs SE FASTQ contains one read"
    fi
fi

assert_file_exists \
    "${log_out_ser}" \
    "execute trim-fastqs SE serial stdout log exists"

assert_file_exists \
    "${log_err_ser}" \
    "execute trim-fastqs SE serial stderr log exists"

assert_file_exists \
    "${log_out_se}" \
    "execute trim-fastqs SE submit stdout log exists"

assert_file_exists \
    "${log_err_se}" \
    "execute trim-fastqs SE submit stderr log exists"

finish
