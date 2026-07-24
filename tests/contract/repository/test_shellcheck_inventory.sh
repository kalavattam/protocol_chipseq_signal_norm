#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_shellcheck_inventory.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="ShellCheck inventory provenance"

# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


function json_value() {
    local source="${1}"
    local expression="${2}"

    "${TEST_PYTHON}" -c \
        'import json,sys; value=json.load(open(sys.argv[1])); print(eval(sys.argv[2], {"value": value}))' \
        "${source}" "${expression}"
}


function make_fake_prefix() {
    local prefix="${1}"

    mkdir -p "${prefix}/bin"
    cp "${ROOT_REPO}/tests/fixtures/shellcheck/fake_shellcheck.fixture" \
        "${prefix}/bin/shellcheck"
    chmod +x "${prefix}/bin/shellcheck"
    ln -sf "${TEST_PYTHON}" "${prefix}/bin/python"
}


function run_fake() {
    local prefix="${1}"
    local status="${2}"
    local output="${3}"
    local invocation_log="${4}"
    shift 4

    CONDA_PREFIX="${prefix}" \
    CONDA_DEFAULT_ENV=env_protocol \
    FAKE_SHELLCHECK_STATUS="${status}" \
    FAKE_SHELLCHECK_LOG="${invocation_log}" \
        "${TEST_BASH}" "${runner}" \
            --output-dir "${output}" -- "$@"
}


function make_discovery_repo() {
    local root="${1}"

    mkdir -p \
        "${root}/bin" \
        "${root}/dev/audit" \
        "${root}/install/scripts" \
        "${root}/tests/contract"
    cp "${runner}" "${root}/dev/audit/run_shellcheck.sh"
    cp "${ROOT_REPO}/tests/fixtures/shellcheck/bash.fixture" \
        "${root}/bin/discovered.sh"
    cp "${ROOT_REPO}/tests/fixtures/shellcheck/posix.fixture" \
        "${root}/install/scripts/install_envs_entrypoint.sh"
    cp "${ROOT_REPO}/tests/fixtures/shellcheck/bash.fixture" \
        "${root}/tests/contract/discovered.sh"
    cp "${ROOT_REPO}/tests/fixtures/shellcheck/bash.fixture" \
        "${root}/tests/contract/deleted.sh"
    git -C "${root}" init -q
    git -C "${root}" add bin dev install tests
    rm "${root}/tests/contract/deleted.sh"
}


print_section "${TEST_NAME}"

runner="${ROOT_REPO}/dev/audit/run_shellcheck.sh"
environment="${ROOT_REPO}/install/envs/env_protocol.yml"
actual_log="${TEST_DIR_LOG}/shellcheck/provenance.log"
nonactivated_log="${TEST_DIR_LOG}/shellcheck/nonactivated.log"
mkdir -p "$(dirname "${actual_log}")"

assert_pattern_found \
    "${environment}" \
    '^  - shellcheck=0\.10\.0' \
    "env_protocol pins ShellCheck 0.10.0"

if [[ ! -x "${TEST_ENV_PREFIX}/bin/shellcheck" ]]; then
    record_fail "env_protocol ShellCheck is not executable"
else
    record_pass "env_protocol ShellCheck is an explicit executable"
fi

if run_capture \
    "activated-prefix ShellCheck provenance" \
    "${actual_log}" \
    env CONDA_PREFIX="${TEST_ENV_PREFIX}" CONDA_DEFAULT_ENV=env_protocol \
        PATH="/opt/homebrew/bin:${PATH}" \
        "${TEST_BASH}" "${runner}" --provenance-only
then
    assert_pattern_found \
        "${actual_log}" \
        "^ShellCheck executable: ${TEST_ENV_PREFIX}/bin/shellcheck$" \
        "activated caller uses the explicit prefix executable"
else
    record_fail "activated-prefix ShellCheck provenance failed"
fi

if run_capture \
    "nonactivated ShellCheck provenance" \
    "${nonactivated_log}" \
    conda run --no-capture-output -n env_protocol /bin/bash \
        "${runner}" --provenance-only
then
    assert_pattern_found \
        "${nonactivated_log}" \
        "^ShellCheck executable: ${TEST_ENV_PREFIX}/bin/shellcheck$" \
        "nonactivated caller uses the explicit prefix executable"
    assert_pattern_absent \
        "${nonactivated_log}" \
        '/opt/homebrew/bin/shellcheck' \
        "nonactivated caller does not select Homebrew ShellCheck"
else
    record_fail "nonactivated ShellCheck provenance failed"
fi

fake_prefix="${TEST_DIR_TMP}/shellcheck_prefix"
make_fake_prefix "${fake_prefix}"
bash_fixture="${ROOT_REPO}/tests/fixtures/shellcheck/bash.fixture"
bootstrap="${ROOT_REPO}/install/scripts/install_envs_entrypoint.sh"

discovery_root="${TEST_DIR_TMP}/shellcheck_discovery"
make_discovery_repo "${discovery_root}"
discovery_log="${TEST_DIR_TMP}/shellcheck_discovery.log"
if CONDA_PREFIX="${fake_prefix}" CONDA_DEFAULT_ENV=env_protocol \
    "${TEST_BASH}" \
        "${discovery_root}/dev/audit/run_shellcheck.sh" --list \
        > "${discovery_log}" 2>&1
then
    if grep -Fq $'bash\t'"${discovery_root}/bin/discovered.sh" \
        "${discovery_log}" && \
        grep -Fq $'bash\t'"${discovery_root}/tests/contract/discovered.sh" \
            "${discovery_log}" && \
        grep -Fq $'sh\t'"${discovery_root}/install/scripts/install_envs_entrypoint.sh" \
            "${discovery_log}" && \
        ! grep -Fq 'tests/contract/deleted.sh' "${discovery_log}"
    then
        record_pass "default discovery assigns maintained Bash and POSIX fixtures"
    else
        record_fail "default discovery omitted or misclassified a fixture"
    fi
else
    record_fail "fixture-backed default discovery failed"
fi

for expected_status in 0 1 3; do
    output="${TEST_DIR_TMP}/shellcheck_status_${expected_status}"
    invocation_log="${output}/invocations.log"
    mkdir -p "${output}"
    set +e
    run_fake \
        "${fake_prefix}" "${expected_status}" "${output}" \
        "${invocation_log}" "${bash_fixture}" \
        > "${output}/runner.log" 2>&1
    observed_status="$?"
    set -e
    summary="${output}/summary.json"
    if [[ ! -s "${summary}" ]]; then
        record_fail "status ${expected_status} did not retain a summary"
        continue
    fi
    if (( expected_status == 0 && observed_status == 0 )); then
        record_pass "status 0 is clean"
    elif (( expected_status == 1 && observed_status == 1 )); then
        record_pass "status 1 completes with findings"
    elif (( expected_status > 1 && observed_status > 1 )); then
        record_pass "status ${expected_status} is an infrastructure failure"
    else
        record_fail \
            "status ${expected_status} returned ${observed_status}"
    fi
done

status_one_summary="${TEST_DIR_TMP}/shellcheck_status_1/summary.json"
if [[ "$(json_value "${status_one_summary}" 'value["finding_count"]')" == "1" ]]; then
    record_pass "status-1 structured evidence retains its finding count"
else
    record_fail "status-1 structured evidence has the wrong finding count"
fi

split_output="${TEST_DIR_TMP}/shellcheck_split"
split_log="${split_output}/invocations.log"
mkdir -p "${split_output}"
if run_fake \
    "${fake_prefix}" 0 "${split_output}" "${split_log}" \
    "${bash_fixture}" "${bootstrap}" \
    > "${split_output}/runner.log" 2>&1
then
    if grep -Fq -- '--shell=bash' "${split_log}" && \
        grep -Fq -- '--shell=sh' "${split_log}" && \
        [[ "$(json_value "${split_output}/summary.json" \
            'value["languages"]["bash"]["file_count"]')" == "1" ]] && \
        [[ "$(json_value "${split_output}/summary.json" \
            'value["languages"]["sh"]["file_count"]')" == "1" ]]
    then
        record_pass "Bash and the POSIX bootstrap use separate language modes"
    else
        record_fail "ShellCheck language-mode discovery was not split"
    fi
else
    record_fail "fixture-backed language split failed"
fi

if grep -Fq -- "-x -P ${ROOT_REPO}" "${split_log}"; then
    record_pass "ShellCheck retains external sources and repository source path"
else
    record_fail "ShellCheck invocation lost -x or the repository source path"
fi

finish
