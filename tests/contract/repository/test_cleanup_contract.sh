#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_cleanup_contract.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="scoped cleanup contract"

# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


function initialize_fixture_repo() {
    local root="${1}"
    local workflow=""

    mkdir -p \
        "${root}/tests/support" \
        "${root}/artifacts/tests/logs" \
        "${root}/artifacts/tests/package" \
        "${root}/artifacts/tests/pycache" \
        "${root}/artifacts/tests/pytest_cache" \
        "${root}/artifacts/tests/tmp"
    cp "${ROOT_REPO}/tests/support/clean_artifacts.sh" \
        "${root}/tests/support/clean_artifacts.sh"
    for workflow in \
        align_fastqs \
        calculate_scaling_factor \
        compute_signal \
        download_fastqs \
        filter_alignments \
        trim_fastqs
    do
        mkdir -p "${root}/tests/fixtures/${workflow}"
        printf '# fixture\n' > "${root}/tests/fixtures/${workflow}/README.md"
        printf '#!/usr/bin/env bash\n' > "${root}/tests/fixtures/${workflow}/make.sh"
    done
    printf '%s\n' \
        '/artifacts/tests/' \
        'tests/fixtures/*/*' \
        '!tests/fixtures/*/README.md' \
        '!tests/fixtures/*/make.sh' \
        > "${root}/.gitignore"
    git -C "${root}" init -q
    git -C "${root}" add .gitignore tests
}


function cleanup_command() {
    local root="${1}"

    shift 1
    "${TEST_BASH}" "${root}/tests/support/clean_artifacts.sh" "$@"
}


print_section "${TEST_NAME}"

fixture_repo="$(mktemp -d "${TEST_DIR_TMP}/cleanup_repo.XXXXXX")"
initialize_fixture_repo "${fixture_repo}"
fixture_generated="${fixture_repo}/tests/fixtures/align_fastqs/generated.dat"
output_generated="${fixture_repo}/artifacts/tests/logs/generated.log"
retained_checkpoint="${fixture_repo}/artifacts/tests/checkpoint/retained.json"
retained_transfer="${fixture_repo}/artifacts/tests/transfer/retained.json"
retained_legacy="${fixture_repo}/artifacts/tests/legacy_migration/retained.json"
retained_migration="${fixture_repo}/artifacts/tests/remaining_standards_migration/retained.json"
printf 'fixture\n' > "${fixture_generated}"
printf 'output\n' > "${output_generated}"
mkdir -p \
    "$(dirname "${retained_checkpoint}")" \
    "$(dirname "${retained_transfer}")" \
    "$(dirname "${retained_legacy}")" \
    "$(dirname "${retained_migration}")"
printf 'retained\n' > "${retained_checkpoint}"
printf 'retained\n' > "${retained_transfer}"
printf 'retained\n' > "${retained_legacy}"
printf 'retained\n' > "${retained_migration}"

if cleanup_command "${fixture_repo}" --dry_run > /dev/null 2>&1; then
    record_fail "dry-run without an explicit selector succeeded"
else
    record_pass "dry-run requires an explicit selector"
fi

if cleanup_command "${fixture_repo}" --outputs > /dev/null 2>&1; then
    record_fail "destructive cleanup without confirmation succeeded"
elif [[ -f "${output_generated}" ]]; then
    record_pass "missing confirmation preserves ignored output"
else
    record_fail "missing confirmation removed ignored output"
fi

if CONFIRM_TEST_CLEANUP=true \
    cleanup_command "${fixture_repo}" --outputs > /dev/null 2>&1
then
    record_fail "true-like cleanup confirmation was accepted"
elif [[ -f "${output_generated}" ]]; then
    record_pass "true-like cleanup confirmation is rejected"
else
    record_fail "invalid confirmation removed ignored output"
fi

if CONFIRM_TEST_CLEANUP=1 \
    cleanup_command "${fixture_repo}" --outputs > /dev/null 2>&1
then
    if [[ \
        ! -e "${output_generated}" && \
        -e "${fixture_generated}" && \
        -e "${retained_checkpoint}" && \
        -e "${retained_transfer}" && \
        -e "${retained_legacy}" && \
        -e "${retained_migration}"
    ]]; then
        record_pass \
            "output cleanup deletes disposable output and preserves retained evidence"
    else
        record_fail "output cleanup crossed its selected boundary"
    fi
else
    record_fail "exactly confirmed output cleanup failed"
fi

if CONFIRM_TEST_CLEANUP=1 \
    cleanup_command "${fixture_repo}" --fixtures > /dev/null 2>&1
then
    if [[ \
        ! -e "${fixture_generated}" && \
        -e "${fixture_repo}/tests/fixtures/align_fastqs/README.md" && \
        -e "${fixture_repo}/tests/fixtures/align_fastqs/make.sh"
    ]]; then
        record_pass "fixture cleanup preserves tracked recipe sources"
    else
        record_fail "fixture cleanup removed or retained the wrong paths"
    fi
else
    record_fail "exactly confirmed fixture cleanup failed"
fi

unexpected="${fixture_repo}/tests/fixtures/align_fastqs/unexpected.txt"
printf 'tracked\n' > "${unexpected}"
git -C "${fixture_repo}" add -f \
    tests/fixtures/align_fastqs/unexpected.txt
if cleanup_command "${fixture_repo}" --dry_run --fixtures > /dev/null 2>&1; then
    record_fail "fixture cleanup accepted an unexpected tracked path"
else
    record_pass "fixture cleanup refuses unexpected tracked paths"
fi

tracked_output="${fixture_repo}/artifacts/tests/logs/tracked.txt"
mkdir -p "$(dirname "${tracked_output}")"
printf 'tracked\n' > "${tracked_output}"
git -C "${fixture_repo}" add -f artifacts/tests/logs/tracked.txt
if cleanup_command "${fixture_repo}" --dry_run --outputs > /dev/null 2>&1; then
    record_fail "output cleanup accepted a tracked output path"
else
    record_pass "output cleanup refuses tracked output paths"
fi

mismatch_root="$(mktemp -d "${TEST_DIR_TMP}/mismatch_root.XXXXXX")"
mismatch_parent="${mismatch_root}/mismatch_parent"
git -C "${mismatch_root}" init -q
mkdir -p "${mismatch_parent}/tests/support"
cp "${ROOT_REPO}/tests/support/clean_artifacts.sh" \
    "${mismatch_parent}/tests/support/clean_artifacts.sh"
if "${TEST_BASH}" "${mismatch_parent}/tests/support/clean_artifacts.sh" \
    --dry_run --outputs > /dev/null 2>&1
then
    record_fail "cleanup accepted a mismatched repository root"
else
    record_pass "cleanup rejects a mismatched repository root"
fi

finish
