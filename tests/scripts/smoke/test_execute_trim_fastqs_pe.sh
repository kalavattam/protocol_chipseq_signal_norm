#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_trim_fastqs_pe.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="execute trim-fastqs PE"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

rec_section "${TEST_NAME}"

if ! is_atria_enabled; then
    rec_skip \
        "Atria execute trim-fastqs PE check disabled;" \
        "set RUN_ATRIA=1 to enable"
    finish
    exit $?
fi


function assert_trimmed_fastq() {
    local outfile="${1:-}"
    local mate="${2:-}"
    local view_fastq="${3:-}"
    local count_reads="${4:-}"

    if [[ -z "${outfile}" ]]; then
        return 0
    fi

    assert_file_nonempty "${outfile}" "execute trim-fastqs PE ${mate} FASTQ output"

    if gzip -t "${outfile}"; then
        rec_pass "execute trim-fastqs PE ${mate} FASTQ output passes gzip integrity"
    else
        rec_fail "execute trim-fastqs PE ${mate} FASTQ output fails gzip integrity"
    fi

    if gzip -cd "${outfile}" > "${view_fastq}"; then
        rec_pass "execute trim-fastqs PE ${mate} FASTQ output can be decompressed"
    else
        rec_fail "execute trim-fastqs PE ${mate} FASTQ output cannot be decompressed"
    fi

    if [[ -s "${view_fastq}" ]]; then
        assert_grep_pattern "${view_fastq}" '^@tiny_trim_pe_pair_1' \
            "execute trim-fastqs PE ${mate} FASTQ contains expected read name"

        awk 'NR % 4 == 1 { n++ } END { print n + 0 }' \
            "${view_fastq}" > "${count_reads}"
        assert_grep_pattern "${count_reads}" '^1$' \
            "execute trim-fastqs PE ${mate} FASTQ contains one read"
    fi
}


#  Define fixture and output paths for local serial execute-layer Atria run
dir_fx="${ROOT_REPO}/tests/trim_fastqs/fixtures"
infile_r1="${dir_fx}/fastq/tiny_pe_R1.fastq.gz"
infile_r2="${dir_fx}/fastq/tiny_pe_R2.fastq.gz"
csv_infile="${infile_r1},${infile_r2}"

tmp="${TEST_DIR_TMP}/execute_trim_fastqs_pe"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/trim_fastqs"
view_fastq_r1="${dir_out}/tiny_pe_R1.trimmed.fastq"
view_fastq_r2="${dir_out}/tiny_pe_R2.trimmed.fastq"
count_reads_r1="${dir_out}/tiny_pe_R1.read_count.txt"
count_reads_r2="${dir_out}/tiny_pe_R2.read_count.txt"

rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_err}" "${dir_log}"

require_files_exist "${infile_r1}" "${infile_r2}" || {
    finish
    exit $?
}

if ! \
    require_atria_env \
        "env_protocol" \
        "${dir_log}/execute_trim_fastqs_pe_env.log"
then
    finish
    exit $?
fi


#  Run execute_trim_fastqs.sh serially on one PE FASTQ pair
log="${dir_log}/execute_trim_fastqs_pe.log"

if \
    run_capture \
        "execute trim-fastqs PE Atria wet run" "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/execute_trim_fastqs.sh" \
            --threads 1 \
            --csv_infile "${csv_infile}" \
            --dir_out "${dir_out}" \
            --sfx_se ".fastq.gz" \
            --sfx_pe "_R1.fastq.gz" \
            --err_out "${dir_err}" \
            --nam_job "test_execute_trim_pe" \
            --max_job 1
then
    rec_pass "execute_trim_fastqs.sh PE Atria wet run exits 0"
else
    rec_fail \
        "execute_trim_fastqs.sh PE Atria wet run failed; see" \
        "$(rec_relpath "${log}")"
fi

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

assert_one_path_found \
    trimmed_outputs_r1 "compressed PE R1 trimmed FASTQ output" \
    "${dir_out}" outfile_r1
assert_one_path_found \
    trimmed_outputs_r2 "compressed PE R2 trimmed FASTQ output" \
    "${dir_out}" outfile_r2

assert_trimmed_fastq "${outfile_r1}" "R1" "${view_fastq_r1}" "${count_reads_r1}"
assert_trimmed_fastq "${outfile_r2}" "R2" "${view_fastq_r2}" "${count_reads_r2}"

assert_file_exists \
    "${dir_err}/test_execute_trim_pe_ser.stdout.txt" \
    "execute trim-fastqs PE serial stdout log exists"
assert_file_exists \
    "${dir_err}/test_execute_trim_pe_ser.stderr.txt" \
    "execute trim-fastqs PE serial stderr log exists"
assert_file_exists \
    "${dir_err}/test_execute_trim_pe.tiny_pe.stdout.txt" \
    "execute trim-fastqs PE submit stdout log exists"
assert_file_exists \
    "${dir_err}/test_execute_trim_pe.tiny_pe.stderr.txt" \
    "execute trim-fastqs PE submit stderr log exists"

finish
