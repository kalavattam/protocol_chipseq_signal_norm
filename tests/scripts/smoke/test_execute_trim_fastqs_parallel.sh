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

print_section "${TEST_NAME}"

if ! {
    is_atria_enabled && is_parallel_enabled
}; then
    record_skip \
        "Atria GNU Parallel execute trim-fastqs check disabled;" \
        "set RUN_ATRIA=1 RUN_PARALLEL=1 to enable"
    finish
    exit $?
fi


#  Define fixture and output paths for a local GNU Parallel Atria wet run
dir_fx="${ROOT_REPO}/tests/trim_fastqs/fixtures"
in_se="${dir_fx}/fastq/tiny_se.fastq.gz"
in_r1="${dir_fx}/fastq/tiny_pe_R1.fastq.gz"
in_r2="${dir_fx}/fastq/tiny_pe_R2.fastq.gz"
csv_in="${in_se};${in_r1},${in_r2}"

tmp="${TEST_DIR_TMP}/execute_trim_fastqs_parallel"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/trim_fastqs"

cfg="${dir_err}/test_execute_trim_parallel.config_parallel.txt"

vw_se="${dir_out}/tiny_se.trimmed.fastq"
vw_r1="${dir_out}/tiny_pe_R1.trimmed.fastq"
vw_r2="${dir_out}/tiny_pe_R2.trimmed.fastq"

cnt_se="${dir_out}/tiny_se.read_count.txt"
cnt_r1="${dir_out}/tiny_pe_R1.read_count.txt"
cnt_r2="${dir_out}/tiny_pe_R2.read_count.txt"

rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_err}" "${dir_log}"

require_files_nonempty \
    "${in_se}" \
    "${in_r1}" \
    "${in_r2}" || {
    finish
    exit $?
}

if ! \
    require_env_atria \
        "env_protocol" \
        "${dir_log}/execute_trim_fastqs_parallel_atria_env.log"
then
    finish
    exit $?
fi

if ! \
    require_env_parallel \
        "env_protocol" \
        "${dir_log}/execute_trim_fastqs_parallel_parallel_env.log"
then
    finish
    exit $?
fi


#  Run execute_trim_fastqs.sh through local GNU Parallel for one SE and one PE
#+ input entry; with two entries and '--max_job 2', the execute wrapper writes
#+ a GNU Parallel config and dispatches the per-entry submit commands through
#+ 'parallel' instead of the serial branch
log="${dir_log}/execute_trim_fastqs_parallel.log"

if \
    run_capture \
        "execute trim-fastqs GNU Parallel Atria wet run" \
        "${log}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/execute_trim_fastqs.sh" \
            --threads 2 \
            --csv_infile "${csv_in}" \
            --dir_out "${dir_out}" \
            --sfx_se ".fastq.gz" \
            --sfx_pe "_R1.fastq.gz" \
            --err_out "${dir_err}" \
            --nam_job "test_execute_trim_parallel" \
            --max_job 2
then
    record_pass "execute_trim_fastqs.sh GNU Parallel Atria wet run exits 0"
else
    record_fail \
        "execute_trim_fastqs.sh GNU Parallel Atria wet run failed; see" \
        "$(print_relpath "${log}")"
fi

assert_file_nonempty \
    "${cfg}" \
    "execute trim-fastqs GNU Parallel config"

if [[ -s "${cfg}" ]]; then
    assert_pattern_found \
        "${cfg}" \
        "${TEST_BASH} ${ROOT_REPO}/scripts/submit_trim_fastqs.sh" \
        "execute trim-fastqs GNU Parallel config uses Bash-prefixed submit command"
fi

# shellcheck disable=2034
{
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
}

assert_path_found \
    trimmed_outputs_se \
    "compressed SE trimmed FASTQ output" \
    "${dir_out}" \
    out_se
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

# shellcheck disable=2154
{
    assert_fastq_gzip \
        "${out_se}" \
        '^@tiny_trim_se_read_1$' \
        "${vw_se}" \
        "${cnt_se}" \
        "" \
        "execute trim-fastqs GNU Parallel SE FASTQ output"
    assert_fastq_gzip \
        "${out_r1}" \
        '^@tiny_trim_pe_pair_1' \
        "${vw_r1}" \
        "${cnt_r1}" \
        "" \
        "execute trim-fastqs GNU Parallel PE R1 FASTQ output"
    assert_fastq_gzip \
        "${out_r2}" \
        '^@tiny_trim_pe_pair_1' \
        "${vw_r2}" \
        "${cnt_r2}" \
        "" \
        "execute trim-fastqs GNU Parallel PE R2 FASTQ output"
}

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
