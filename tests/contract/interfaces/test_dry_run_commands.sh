#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_dry_run_commands.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="dry-run and expected failure"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


function assert_fixed_before() {
    local file="${1}"
    local first="${2}"
    local second="${3}"
    local label="${4}"
    local first_line
    local second_line

    first_line="$(
        grep -n -m 1 -F -- "${first}" "${file}" | cut -d: -f1
    )"
    second_line="$(
        grep -n -m 1 -F -- "${second}" "${file}" | cut -d: -f1
    )"
    if \
        [[ -n "${first_line}" ]] \
        && [[ -n "${second_line}" ]] \
        && (( first_line < second_line ))
    then
        record_pass "${label}"
    else
        record_fail "${label}"
    fi
}


#  Create small temporary inputs for dry-run smoke checks
tmp="${TEST_DIR_TMP}/dry_run"
fil_sample="${tmp}/in/sample.txt"
fil_header="${tmp}/out/header.tsv"
dir_compress="${tmp}/compress"
fil_compress="${dir_compress}/large.log"
fil_empty="${dir_compress}/empty.log"

log_symlink="${TEST_DIR_LOG}/dry_run/symlink_files.log"
log_header="${TEST_DIR_LOG}/dry_run/write_header.log"
log_compress="${TEST_DIR_LOG}/dry_run/compress_remove_files.log"
log_compress_short="${TEST_DIR_LOG}/dry_run/compress_remove_files_short.log"
log_compress_help="${TEST_DIR_LOG}/dry_run/compress_remove_files_help.log"
log_compress_mutual="${TEST_DIR_LOG}/dry_run/compress_mutual.log"
log_symlink_missing="${TEST_DIR_LOG}/expected_fail/symlink_missing_args.log"
log_find_missing="${TEST_DIR_LOG}/expected_fail/find_missing_args.log"


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${tmp}/in" "${tmp}/out" "${tmp}/logs" "${dir_compress}"
printf 'smoke\n' > "${fil_sample}"
printf '%02048d' 0 > "${fil_compress}"
: > "${fil_empty}"

if \
    run_capture \
        "compress help" \
        "${log_compress_help}" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/compress_remove_files.sh" --help
then
    assert_pattern_found \
        "${log_compress_help}" \
        '^  -dr, --dry_run : flag$' \
        "compress help exposes the public dry-run pair"
    assert_pattern_found \
        "${log_compress_help}" \
        '^    \[--help\] \[--dry_run\]$' \
        "compress Usage places dry-run after help"
    assert_fixed_before \
        "${log_compress_help}" \
        "  -dr, --dry_run : flag" \
        "  -t, --thr, --threads : int" \
        "compress Parameters places dry-run before ordinary options"
else
    record_fail "compress_remove_files.sh --help failed"
fi
assert_fixed_before \
    "${ROOT_REPO}/bin/compress_remove_files.sh" \
    "            -dr|--dry_run)" \
    "            -t|--thr|--threads)" \
    "compress parser places dry-run before ordinary options"

if \
    run_capture \
        "symlink dry-run" \
        "${log_symlink}" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/symlink_files.sh" \
            --dry_run \
            --csv_fil_in "${fil_sample}" \
            --dir_out "${tmp}/out"
then
    record_pass "symlink_files.sh minimal --dry_run"
else
    record_fail \
        "symlink_files.sh minimal --dry_run; see" \
        "$(print_relpath "${log_symlink}")"
fi

#  Check that write_header.sh dry-run reports work without writing output
if \
    run_capture \
        "write_header dry-run" \
        "${log_header}" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/write_header.sh" \
            --dry_run \
            --mode siq \
            --fil_out "${fil_header}"
then
    if [[ -e "${fil_header}" ]]; then
        record_fail \
            "write_header.sh --dry_run created an output file; see" \
            "$(print_relpath "${log_header}")"
    else
        assert_pattern_found \
            "${log_header}" \
            "Dry run: would create" \
            "write_header.sh minimal --dry_run"
    fi
else
    record_fail \
        "write_header.sh minimal --dry_run; see" \
        "$(print_relpath "${log_header}")"
fi

#  Check canonical compress dry run without compressing or deleting files
compress_before="$(< "${fil_compress}")"
if \
    run_capture \
        "compress dry-run" \
        "${log_compress}" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/compress_remove_files.sh" \
            --dir_fnd "${dir_compress}" \
            --pattern '*.log' \
            --size 1 \
            --dry_run
then
    assert_pattern_found \
        "${log_compress}" \
        '^## Call to find for files larger than 1k ##$' \
        "compress dry run prints the plan"
    assert_pattern_found \
        "${log_compress}" \
        '^## Results of find command for files larger than 1k ##$' \
        "compress dry run prints results"
    if [[ ! -f "${fil_compress}" || -e "${fil_compress}.gz" ]]; then
        record_fail "compress dry run compressed the matching file"
    elif [[ "$(< "${fil_compress}")" != "${compress_before}" ]]; then
        record_fail "compress dry run changed the matching file"
    else
        record_pass "compress dry run leaves matching content unchanged"
    fi
    if [[ -e "${fil_empty}" ]]; then
        record_pass "compress dry run does not delete empty files"
    else
        record_fail "compress dry run deleted an empty file"
    fi
else
    record_fail \
        "compress_remove_files.sh --dry_run failed; see" \
        "$(print_relpath "${log_compress}")"
fi

# Prove the public short form selects the same non-mutating behavior.
if \
    run_capture \
        "compress short dry-run" \
        "${log_compress_short}" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/compress_remove_files.sh" \
            -dr \
            --dir_fnd "${dir_compress}" \
            --pattern '*.log' \
            --size 1
then
    assert_pattern_found \
        "${log_compress_short}" \
        '^## Results of find command for files larger than 1k ##$' \
        "compress -dr selects dry-run behavior"
else
    record_fail "compress_remove_files.sh -dr failed"
fi
if [[ -f "${fil_compress}" && ! -e "${fil_compress}.gz" && -e "${fil_empty}" ]]
then
    record_pass "compress -dr is non-mutating"
else
    record_fail "compress -dr changed target files"
fi

# Prove every retired spelling and bounded near miss is rejected.
for alias in \
    -ce \
    -cu \
    --chk_exc \
    --chk-exc \
    --chk_exu \
    --chk-exu \
    --dry \
    --dry-run \
    --dry__run \
    --dry_-run \
    --dry-_run \
    --dry_runs \
    --dry_run_ \
    --dryrun \
    --Dry_run \
    -dry_run \
    ---dry_run
do
    safe_alias="${alias//[^A-Za-z0-9]/_}"
    alias_log="${TEST_DIR_LOG}/dry_run/compress_alias_${safe_alias}.log"
    if \
        run_capture \
            "compress alias ${alias}" \
            "${alias_log}" \
            "${TEST_BASH}" \
            "${ROOT_REPO}/bin/compress_remove_files.sh" \
                "${alias}" \
                --not_a_real_option
    then
        record_fail "compress retired alias ${alias} unexpectedly succeeded"
    else
        assert_pattern_found \
            "${alias_log}" \
            "unknown option/parameter passed: '${alias}'" \
            "compress parser rejects ${alias}"
    fi
done

if \
    run_capture \
        "compress mutually exclusive modes" \
        "${log_compress_mutual}" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/compress_remove_files.sh" \
            --dir_fnd "${dir_compress}" \
            --pattern '*.log' \
            --size 1 \
            --chk_con \
            --dry_run
then
    record_fail "compress --chk_con --dry_run unexpectedly succeeded"
else
    assert_pattern_found \
        "${log_compress_mutual}" \
        "only one of '--chk_con' or '--dry_run'" \
        "compress check and dry-run modes remain mutually exclusive"
fi
if [[ -f "${fil_compress}" && ! -e "${fil_compress}.gz" && -e "${fil_empty}" ]]
then
    record_pass "compress mutual-exclusion failure is non-mutating"
else
    record_fail "compress mutual-exclusion failure changed target files"
fi

#  Check that selected wrappers fail clearly when required arguments are absent
if \
    run_capture \
        "symlink missing args" \
        "${log_symlink_missing}" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/symlink_files.sh" --dry_run
then
    record_fail "symlink_files.sh missing required args unexpectedly succeeded"
else
    assert_pattern_found \
        "${log_symlink_missing}" \
        'error' \
        "symlink_files.sh missing args emits useful error"
fi

if \
    run_capture \
        "find missing args" \
        "${log_find_missing}" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/find_files.sh" --pattern '*.sh'
then
    record_fail "find_files.sh missing required args unexpectedly succeeded"
else
    assert_pattern_found \
        "${log_find_missing}" \
        'error' \
        "find_files.sh missing args emits useful error"
fi

record_skip \
    "execute_* dry-run wrappers require realistic files and/or environment" \
    "checks; covered later by integration fixtures"

record_skip \
    "submit_* dry-run wrappers may activate Conda/Mamba or require Slurm" \
    "context; covered later by integration fixtures"

finish
