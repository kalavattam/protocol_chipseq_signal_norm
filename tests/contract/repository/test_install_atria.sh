#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_install_atria.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="Atria installer"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"

script="${ROOT_REPO}/install/scripts/install_atria.sh"
tool="${ROOT_REPO}/tests/fixtures/install_atria/tool"
log="${TEST_DIR_LOG}/install_atria.log"

print_section "${TEST_NAME}"

for mode in fail reuse update; do
    mode_log="${TEST_DIR_LOG}/install_atria_${mode}.log"
    if \
        run_capture \
            "Atria ${mode} dry run" \
            "${mode_log}" \
            "${TEST_BASH}" "${script}" \
                --dry_run \
                --if_exists "${mode}" \
                --dir_install "${TEST_DIR_TMP}/${mode}"
    then
        record_pass "Atria ${mode} dry run exits 0"
    else
        record_fail "Atria ${mode} dry run failed"
    fi
done

latest_log="${TEST_DIR_LOG}/install_atria_latest.log"
if run_capture "Atria latest dry run" "${latest_log}" \
    "${TEST_BASH}" "${script}" --dry_run --if_exists update \
    --v_atria latest --dir_install "${TEST_DIR_TMP}/latest"
then
    assert_pattern_found \
        "${latest_log}" \
        "would resolve the latest stable" \
        "latest dry run does not query upstream"
else
    record_fail "Atria latest dry run failed"
fi

install="${TEST_DIR_TMP}/install"
fake_log="${TEST_DIR_LOG}/install_atria_fake.log"
rm -rf "${install}"
: > "${fake_log}"
if PATH="${tool}:${PATH}" ATRIA_FAKE_LOG="${fake_log}" \
    CONDA_DEFAULT_ENV=env_protocol \
    "${TEST_BASH}" "${script}" --if_exists update --v_atria latest \
        --dir_install "${install}" > "${log}" 2>&1
then
    assert_pattern_found \
        "${fake_log}" \
        'git ls-remote --tags --refs' \
        "latest resolves tags through Git"
    assert_pattern_found \
        "${fake_log}" \
        'git checkout --detach v4.1.5' \
        "latest checks out the highest stable tag"
    assert_pattern_found \
        "${log}" \
        'Atria tag v4.1.5 (reported v4.1.5)' \
        "summary separates Atria tag and reported version"
else
    record_fail \
        "fake latest installation failed; see $(print_relpath "${log}")"
fi

: > "${fake_log}"
if PATH="${tool}:${PATH}" ATRIA_FAKE_LOG="${fake_log}" \
    CONDA_DEFAULT_ENV=env_protocol \
    "${TEST_BASH}" "${script}" --if_exists update --v_atria latest \
        --dir_install "${install}" > "${log}" 2>&1
then
    assert_pattern_absent \
        "${fake_log}" \
        '^tar ' \
        "matching update does not recreate Julia"
    assert_pattern_absent \
        "${fake_log}" \
        '^git checkout' \
        "matching update does not check out or rebuild Atria"
else
    record_fail "matching fake update failed; see $(print_relpath "${log}")"
fi

finish "$@"
