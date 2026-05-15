#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_trim_fastqs_parallel.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="execute trim-fastqs GNU Parallel"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

rec_section "${TEST_NAME}"

if ! is_atria_enabled || ! is_parallel_enabled; then
    rec_skip \
        "Atria GNU Parallel execute trim-fastqs check disabled;" \
        "set RUN_ATRIA=1 RUN_PARALLEL=1 to enable"
    finish
    exit $?
fi


function assert_trimmed_fastq() {
    local outfile="${1:-}"
    local label="${2:-trimmed FASTQ}"
    local read_pattern="${3:-}"
    local view_fastq="${4:-}"
    local count_reads="${5:-}"

    if [[ -z "${outfile}" ]]; then
        return 0
    fi

    assert_file_nonempty "${outfile}" "${label}"

    if gzip -t "${outfile}"; then
        rec_pass "${label} passes gzip integrity"
    else
        rec_fail "${label} fails gzip integrity"
    fi

    if gzip -cd "${outfile}" > "${view_fastq}"; then
        rec_pass "${label} can be decompressed"
    else
        rec_fail "${label} cannot be decompressed"
    fi

    if [[ -s "${view_fastq}" ]]; then
        assert_grep_pattern "${view_fastq}" "${read_pattern}" \
            "${label} contains expected read name"

        awk 'NR % 4 == 1 { n++ } END { print n + 0 }' \
            "${view_fastq}" > "${count_reads}"
        assert_grep_pattern "${count_reads}" '^1$' \
            "${label} contains one read"
    fi
}


#  Define fixture and output paths for a local GNU Parallel Atria wet run
dir_fx="${ROOT_REPO}/tests/trim_fastqs/fixtures"
infile_se="${dir_fx}/fastq/tiny_se.fastq.gz"
infile_r1="${dir_fx}/fastq/tiny_pe_R1.fastq.gz"
infile_r2="${dir_fx}/fastq/tiny_pe_R2.fastq.gz"
csv_infile="${infile_se};${infile_r1},${infile_r2}"

tmp="${TEST_DIR_TMP}/execute_trim_fastqs_parallel"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/trim_fastqs"
config="${dir_err}/test_execute_trim_parallel.config_parallel.txt"
view_fastq_se="${dir_out}/tiny_se.trimmed.fastq"
view_fastq_r1="${dir_out}/tiny_pe_R1.trimmed.fastq"
view_fastq_r2="${dir_out}/tiny_pe_R2.trimmed.fastq"
count_reads_se="${dir_out}/tiny_se.read_count.txt"
count_reads_r1="${dir_out}/tiny_pe_R1.read_count.txt"
count_reads_r2="${dir_out}/tiny_pe_R2.read_count.txt"

rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_err}" "${dir_log}"

require_files_exist "${infile_se}" "${infile_r1}" "${infile_r2}" || {
    finish
    exit $?
}

if ! \
    require_atria_env \
        "env_protocol" \
        "${dir_log}/execute_trim_fastqs_parallel_atria_env.log"
then
    finish
    exit $?
fi

if ! \
    require_parallel_env \
        "env_protocol" \
        "${dir_log}/execute_trim_fastqs_parallel_parallel_env.log"
then
    finish
    exit $?
fi


#  Run execute_trim_fastqs.sh through local GNU Parallel for one SE and one PE
#+ input entry. With two entries and '--max_job 2', the execute wrapper writes
#+ a GNU Parallel config and dispatches the per-entry submit commands through
#+ 'parallel' instead of the serial branch.
log="${dir_log}/execute_trim_fastqs_parallel.log"

if \
    run_capture \
        "execute trim-fastqs GNU Parallel Atria wet run" "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/execute_trim_fastqs.sh" \
            --threads 2 \
            --csv_infile "${csv_infile}" \
            --dir_out "${dir_out}" \
            --sfx_se ".fastq.gz" \
            --sfx_pe "_R1.fastq.gz" \
            --err_out "${dir_err}" \
            --nam_job "test_execute_trim_parallel" \
            --max_job 2
then
    rec_pass "execute_trim_fastqs.sh GNU Parallel Atria wet run exits 0"
else
    rec_fail \
        "execute_trim_fastqs.sh GNU Parallel Atria wet run failed; see" \
        "$(rec_relpath "${log}")"
fi

assert_file_nonempty "${config}" "execute trim-fastqs GNU Parallel config"

if [[ -s "${config}" ]]; then
    assert_grep_pattern \
        "${config}" \
        "${TEST_BASH} ${ROOT_REPO}/scripts/submit_trim_fastqs.sh" \
        "execute trim-fastqs GNU Parallel config uses Bash-prefixed submit command"
fi

mapfile -t trimmed_outputs_se < <(
    find "${dir_out}" -maxdepth 1 -type f \
        \( -name 'tiny_se*.fastq.gz' -o -name 'tiny_se*.fq.gz' \) \
        | sort
)
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
    trimmed_outputs_se "compressed SE trimmed FASTQ output" \
    "${dir_out}" outfile_se
assert_one_path_found \
    trimmed_outputs_r1 "compressed PE R1 trimmed FASTQ output" \
    "${dir_out}" outfile_r1
assert_one_path_found \
    trimmed_outputs_r2 "compressed PE R2 trimmed FASTQ output" \
    "${dir_out}" outfile_r2

assert_trimmed_fastq \
    "${outfile_se}" "execute trim-fastqs GNU Parallel SE FASTQ output" \
    '^@tiny_trim_se_read_1$' "${view_fastq_se}" "${count_reads_se}"
assert_trimmed_fastq \
    "${outfile_r1}" "execute trim-fastqs GNU Parallel PE R1 FASTQ output" \
    '^@tiny_trim_pe_pair_1' "${view_fastq_r1}" "${count_reads_r1}"
assert_trimmed_fastq \
    "${outfile_r2}" "execute trim-fastqs GNU Parallel PE R2 FASTQ output" \
    '^@tiny_trim_pe_pair_1' "${view_fastq_r2}" "${count_reads_r2}"

assert_file_exists \
    "${dir_err}/test_execute_trim_parallel.tiny_se.stdout.txt" \
    "execute trim-fastqs GNU Parallel SE submit stdout log exists"
assert_file_exists \
    "${dir_err}/test_execute_trim_parallel.tiny_se.stderr.txt" \
    "execute trim-fastqs GNU Parallel SE submit stderr log exists"
assert_file_exists \
    "${dir_err}/test_execute_trim_parallel.tiny_pe.stdout.txt" \
    "execute trim-fastqs GNU Parallel PE submit stdout log exists"
assert_file_exists \
    "${dir_err}/test_execute_trim_parallel.tiny_pe.stderr.txt" \
    "execute trim-fastqs GNU Parallel PE submit stderr log exists"

finish
