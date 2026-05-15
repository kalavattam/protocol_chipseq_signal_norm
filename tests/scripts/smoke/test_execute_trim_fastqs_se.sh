#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_trim_fastqs_se.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="execute trim-fastqs SE"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

rec_section "${TEST_NAME}"

if ! is_atria_enabled; then
    rec_skip \
        "Atria execute trim-fastqs SE check disabled;" \
        "set RUN_ATRIA=1 to enable"
    finish
    exit $?
fi


#  Define fixture and output paths for local serial execute-layer Atria run
dir_fx="${ROOT_REPO}/tests/trim_fastqs/fixtures"
infile_se="${dir_fx}/fastq/tiny_se.fastq.gz"

tmp="${TEST_DIR_TMP}/execute_trim_fastqs_se"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/trim_fastqs"
view_fastq="${dir_out}/tiny_se.trimmed.fastq"
count_reads="${dir_out}/tiny_se.read_count.txt"

rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_err}" "${dir_log}"

require_files_exist "${infile_se}" || {
    finish
    exit $?
}

if ! \
    require_atria_env \
        "env_protocol" \
        "${dir_log}/execute_trim_fastqs_se_env.log"
then
    finish
    exit $?
fi


#  Run execute_trim_fastqs.sh serially on one SE FASTQ fixture
log="${dir_log}/execute_trim_fastqs_se.log"

if \
    run_capture \
        "execute trim-fastqs SE Atria wet run" "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/execute_trim_fastqs.sh" \
            --threads 1 \
            --csv_infile "${infile_se}" \
            --dir_out "${dir_out}" \
            --sfx_se ".fastq.gz" \
            --sfx_pe "_R1.fastq.gz" \
            --err_out "${dir_err}" \
            --nam_job "test_execute_trim_se" \
            --max_job 1
then
    rec_pass "execute_trim_fastqs.sh SE Atria wet run exits 0"
else
    rec_fail \
        "execute_trim_fastqs.sh SE Atria wet run failed; see" \
        "$(rec_relpath "${log}")"
fi

mapfile -t trimmed_outputs < <(
    find "${dir_out}" -maxdepth 1 -type f \
        \( -name 'tiny_se*.fastq.gz' -o -name 'tiny_se*.fq.gz' \) \
        | sort
)

assert_one_path_found \
    trimmed_outputs "compressed SE trimmed FASTQ output" "${dir_out}" outfile

if [[ -n "${outfile}" ]]; then
    assert_file_nonempty "${outfile}" "execute trim-fastqs SE FASTQ output"

    if gzip -t "${outfile}"; then
        rec_pass "execute trim-fastqs SE FASTQ output passes gzip integrity"
    else
        rec_fail "execute trim-fastqs SE FASTQ output fails gzip integrity"
    fi

    if gzip -cd "${outfile}" > "${view_fastq}"; then
        rec_pass "execute trim-fastqs SE FASTQ output can be decompressed"
    else
        rec_fail "execute trim-fastqs SE FASTQ output cannot be decompressed"
    fi

    if [[ -s "${view_fastq}" ]]; then
        assert_grep_pattern "${view_fastq}" '^@tiny_trim_se_read_1$' \
            "execute trim-fastqs SE FASTQ contains expected read name"

        awk 'NR % 4 == 1 { n++ } END { print n + 0 }' \
            "${view_fastq}" > "${count_reads}"
        assert_grep_pattern "${count_reads}" '^1$' \
            "execute trim-fastqs SE FASTQ contains one read"
    fi
fi

assert_file_exists \
    "${dir_err}/test_execute_trim_se_ser.stdout.txt" \
    "execute trim-fastqs SE serial stdout log exists"
assert_file_exists \
    "${dir_err}/test_execute_trim_se_ser.stderr.txt" \
    "execute trim-fastqs SE serial stderr log exists"
assert_file_exists \
    "${dir_err}/test_execute_trim_se.tiny_se.stdout.txt" \
    "execute trim-fastqs SE submit stdout log exists"
assert_file_exists \
    "${dir_err}/test_execute_trim_se.tiny_se.stderr.txt" \
    "execute trim-fastqs SE submit stderr log exists"

finish
