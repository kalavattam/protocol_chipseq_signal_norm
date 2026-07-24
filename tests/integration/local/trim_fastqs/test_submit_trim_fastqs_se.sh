#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_submit_trim_fastqs_se.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="submit trim-fastqs SE"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


if ! \
    is_atria_enabled
then
    record_skip \
        "Atria trim-fastqs SE check disabled; set RUN_ATRIA=1 to enable"
    finish
    exit $?
fi


#  Define fixture and output paths for a minimal submit-layer Atria wet run
dir_fx="${ROOT_REPO}/tests/fixtures/trim_fastqs"
in_se="${dir_fx}/fastq/se/tiny_se.fastq.gz"

tmp="${TEST_DIR_TMP}/submit_trim_fastqs_se"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/trim_fastqs"

vw_fq="${dir_out}/tiny_se.trimmed.fastq"
cnt="${dir_out}/tiny_se.read_count.txt"

log_env="${dir_log}/submit_trim_fastqs_se_env.log"
log_run="${dir_log}/submit_trim_fastqs_se.log"
log_out_se="${dir_err}/test_submit_trim_se.tiny_se.stdout.txt"
log_err_se="${dir_err}/test_submit_trim_se.tiny_se.stderr.txt"


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty \
    "${in_se}" || {
    finish
    exit $?
}

# shellcheck disable=SC2154
if ! \
    require_env_atria \
        "${env_nam}" \
        "${log_env}"
then
    finish
    exit $?
fi


#  Run submit_trim_fastqs.sh on one SE FASTQ fixture through real Atria
if \
    run_capture \
        "submit trim-fastqs SE Atria wet run" \
        "${log_run}" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/submit_trim_fastqs.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/bin" \
            --threads 1 \
            --csv_fil_in "${in_se}" \
            --dir_out "${dir_out}" \
            --sfx_se ".fastq.gz" \
            --sfx_pe "_R1.fastq.gz" \
            --dir_eo "${dir_err}" \
            --nam_job "test_submit_trim_se"
then
    record_pass "submit_trim_fastqs.sh SE Atria wet run exits 0"
else
    record_fail \
        "submit_trim_fastqs.sh SE Atria wet run failed; see" \
        "$(print_relpath "${log_run}")"
fi

# shellcheck disable=SC2034
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
        "submit trim-fastqs SE FASTQ output"

    if \
        gzip -t "${fil_out}"
    then
        record_pass "submit trim-fastqs SE FASTQ output passes gzip integrity"
    else
        record_fail "submit trim-fastqs SE FASTQ output fails gzip integrity"
    fi

    if \
        gzip -cd "${fil_out}" > "${vw_fq}"
    then
        record_pass "submit trim-fastqs SE FASTQ output can be decompressed"
    else
        record_fail "submit trim-fastqs SE FASTQ output cannot be decompressed"
    fi

    if [[ -s "${vw_fq}" ]]; then
        assert_pattern_found \
            "${vw_fq}" \
            '^@tiny_trim_se_read_1$' \
            "submit trim-fastqs SE FASTQ contains expected read name"

        awk 'NR % 4 == 1 { n++ } END { print n + 0 }' "${vw_fq}" > "${cnt}"
        assert_pattern_found \
            "${cnt}" \
            '^1$' \
            "submit trim-fastqs SE FASTQ contains one read"
    fi
fi

assert_file_exists \
    "${log_out_se}" \
    "submit trim-fastqs SE stdout log exists"

assert_file_exists \
    "${log_err_se}" \
    "submit trim-fastqs SE stderr log exists"

finish
