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


function run_fake_with_runner() {
    local shell_path="${1}"
    local runner_path="${2}"
    local prefix="${3}"
    local status="${4}"
    local output="${5}"
    local invocation_log="${6}"
    shift 6

    CONDA_PREFIX="${prefix}" \
    CONDA_DEFAULT_ENV=env_protocol \
    FAKE_SHELLCHECK_STATUS="${status}" \
    FAKE_SHELLCHECK_LOG="${invocation_log}" \
        "${shell_path}" "${runner_path}" \
            --output-dir "${output}" -- "$@"
}


function run_fake() {
    run_fake_with_runner "${TEST_BASH}" "${runner}" "$@"
}


function run_fake_system_bash() {
    run_fake_with_runner /bin/bash "${runner}" "$@"
}


function invocation_count() {
    local log="${1}"
    local pattern="${2}"

    grep -Fc -- "${pattern}" "${log}" 2> /dev/null || true
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


function make_legacy_runner() {
    local destination="${1}"

    "${TEST_PYTHON}" - "${runner}" "${destination}" << 'PY'
from pathlib import Path
import sys


source = Path(sys.argv[1]).read_text(encoding="utf-8")
destination = Path(sys.argv[2])
replacements = {
    '''    if (( ${#arr_bash_paths[@]} > 0 )); then
        for resolved in "${arr_bash_paths[@]}"; do
            printf 'bash\\t%s\\n' "${resolved}"
        done
    fi
''': '''    for resolved in "${arr_bash_paths[@]}"; do
        printf 'bash\\t%s\\n' "${resolved}"
    done
''',
    '''    if (( ${#arr_sh_paths[@]} > 0 )); then
        for resolved in "${arr_sh_paths[@]}"; do
            printf 'sh\\t%s\\n' "${resolved}"
        done
    fi
''': '''    for resolved in "${arr_sh_paths[@]}"; do
        printf 'sh\\t%s\\n' "${resolved}"
    done
''',
    '''if (( ${#arr_bash_paths[@]} > 0 )); then
    run_language bash "${bash_raw}" "${arr_bash_paths[@]}"
else
    run_language bash "${bash_raw}"
fi
''': '''run_language bash "${bash_raw}" "${arr_bash_paths[@]}"
''',
    '''if (( ${#arr_sh_paths[@]} > 0 )); then
    run_language sh "${sh_raw}" "${arr_sh_paths[@]}"
else
    run_language sh "${sh_raw}"
fi
''': '''run_language sh "${sh_raw}" "${arr_sh_paths[@]}"
''',
}
for expected, replacement in replacements.items():
    if expected not in source:
        raise SystemExit(f"missing guarded expansion: {expected!r}")
    source = source.replace(expected, replacement, 1)
destination.write_text(source, encoding="utf-8")
PY
    chmod +x "${destination}"
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

legacy_runner="${TEST_DIR_TMP}/run_shellcheck_legacy.sh"
legacy_output="${TEST_DIR_TMP}/shellcheck_legacy"
legacy_log="${legacy_output}/invocations.log"
mkdir -p "${legacy_output}"
make_legacy_runner "${legacy_runner}"
set +e
run_fake_with_runner /bin/bash "${legacy_runner}" \
    "${fake_prefix}" 0 "${legacy_output}" "${legacy_log}" "${bash_fixture}" \
    > "${legacy_output}/runner.log" 2>&1
legacy_status="$?"
set -e
if (( legacy_status > 0 )) && \
    grep -Fq 'unbound variable' "${legacy_output}/runner.log"
then
    record_pass "legacy Bash-only runner exposes the nounset array fault"
else
    record_fail "legacy Bash-only runner did not expose the nounset array fault"
fi

bash_only_output="${TEST_DIR_TMP}/shellcheck_bash_only"
bash_only_log="${bash_only_output}/invocations.log"
mkdir -p "${bash_only_output}"
if run_fake_system_bash \
    "${fake_prefix}" 0 "${bash_only_output}" "${bash_only_log}" \
    "${bash_fixture}" > "${bash_only_output}/runner.log" 2>&1
then
    if [[ "$(json_value "${bash_only_output}/summary.json" \
        'value["languages"]["bash"]["file_count"]')" == "1" ]] && \
        [[ "$(json_value "${bash_only_output}/summary.json" \
        'value["languages"]["sh"]["file_count"]')" == "0" ]] && \
        [[ "$(invocation_count "${bash_only_log}" '--shell=bash')" == "1" ]] && \
        [[ "$(invocation_count "${bash_only_log}" '--shell=sh')" == "0" ]] && \
        [[ -s "${bash_only_output}/bash_findings.json" ]] && \
        [[ -s "${bash_only_output}/sh_findings.json" ]]
    then
        record_pass "Bash-only scan preserves empty POSIX evidence under nounset"
    else
        record_fail "Bash-only scan lost language counts, calls, or raw evidence"
    fi
else
    record_fail "Bash-only fixture-backed scan failed"
fi

posix_only_output="${TEST_DIR_TMP}/shellcheck_posix_only"
posix_only_log="${posix_only_output}/invocations.log"
mkdir -p "${posix_only_output}"
if run_fake_system_bash \
    "${fake_prefix}" 0 "${posix_only_output}" "${posix_only_log}" \
    "${bootstrap}" > "${posix_only_output}/runner.log" 2>&1
then
    if [[ "$(json_value "${posix_only_output}/summary.json" \
        'value["languages"]["bash"]["file_count"]')" == "0" ]] && \
        [[ "$(json_value "${posix_only_output}/summary.json" \
        'value["languages"]["sh"]["file_count"]')" == "1" ]] && \
        [[ "$(invocation_count "${posix_only_log}" '--shell=bash')" == "0" ]] && \
        [[ "$(invocation_count "${posix_only_log}" '--shell=sh')" == "1" ]] && \
        [[ -s "${posix_only_output}/bash_findings.json" ]] && \
        [[ -s "${posix_only_output}/sh_findings.json" ]]
    then
        record_pass "POSIX-only scan preserves empty Bash evidence under nounset"
    else
        record_fail "POSIX-only scan lost language counts, calls, or raw evidence"
    fi
else
    record_fail "POSIX-only fixture-backed scan failed"
fi

bash_list="${TEST_DIR_TMP}/shellcheck_bash_only.list"
posix_list="${TEST_DIR_TMP}/shellcheck_posix_only.list"
if CONDA_PREFIX="${fake_prefix}" CONDA_DEFAULT_ENV=env_protocol \
    /bin/bash "${runner}" --list -- "${bash_fixture}" > "${bash_list}" 2>&1 && \
    CONDA_PREFIX="${fake_prefix}" CONDA_DEFAULT_ENV=env_protocol \
    /bin/bash "${runner}" --list -- "${bootstrap}" > "${posix_list}" 2>&1 && \
    grep -Eq '^bash\t' "${bash_list}" && \
    ! grep -Eq '^sh\t' "${bash_list}" && \
    grep -Eq '^sh\t' "${posix_list}" && \
    ! grep -Eq '^bash\t' "${posix_list}"
then
    record_pass "single-dialect list loops are nounset-safe"
else
    record_fail "single-dialect list loops failed under nounset"
fi

default_output="${TEST_DIR_TMP}/shellcheck_default"
default_log="${default_output}/invocations.log"
mkdir -p "${default_output}"
if run_fake_with_runner /bin/bash \
    "${discovery_root}/dev/audit/run_shellcheck.sh" "${fake_prefix}" 0 \
    "${default_output}" "${default_log}" \
    > "${default_output}/runner.log" 2>&1
then
    if [[ "$(json_value "${default_output}/summary.json" \
        'value["languages"]["bash"]["file_count"]')" == "2" ]] && \
        [[ "$(json_value "${default_output}/summary.json" \
        'value["languages"]["sh"]["file_count"]')" == "1" ]] && \
        [[ "$(invocation_count "${default_log}" '--shell=bash')" == "1" ]] && \
        [[ "$(invocation_count "${default_log}" '--shell=sh')" == "1" ]] && \
        ! grep -Fq 'tests/contract/deleted.sh' "${default_log}"
    then
        record_pass "default scan executes both discovered language passes"
    else
        record_fail "default scan lost discovery or language-pass behavior"
    fi
else
    record_fail "fixture-backed default scan failed"
fi

missing_output="${TEST_DIR_TMP}/shellcheck_missing"
missing_log="${missing_output}/invocations.log"
mkdir -p "${missing_output}"
set +e
run_fake_system_bash \
    "${fake_prefix}" 0 "${missing_output}" "${missing_log}" \
    "${TEST_DIR_TMP}/not_a_shell_file.sh" > "${missing_output}/runner.log" 2>&1
missing_status="$?"
set -e
if (( missing_status == 2 )) && \
    grep -Fq 'ShellCheck path is not a file' "${missing_output}/runner.log" && \
    [[ ! -e "${missing_log}" ]]
then
    record_pass "missing explicit path remains an infrastructure failure"
else
    record_fail "missing explicit path changed ShellCheck failure behavior"
fi

duplicate_output="${TEST_DIR_TMP}/shellcheck_duplicates"
duplicate_log="${duplicate_output}/invocations.log"
mkdir -p "${duplicate_output}"
if run_fake_system_bash \
    "${fake_prefix}" 0 "${duplicate_output}" "${duplicate_log}" \
    "${bash_fixture}" "${bash_fixture}" \
    > "${duplicate_output}/runner.log" 2>&1
then
    duplicate_paths="$(grep -oF -- "${bash_fixture}" "${duplicate_log}" | wc -l)"
    if [[ "$(json_value "${duplicate_output}/summary.json" \
        'value["languages"]["bash"]["file_count"]')" == "2" ]] && \
        [[ "${duplicate_paths//[[:space:]]/}" == "2" ]] && \
        [[ "$(invocation_count "${duplicate_log}" '--shell=bash')" == "1" ]]
    then
        record_pass "duplicate explicit paths retain count and order"
    else
        record_fail "duplicate explicit paths were deduplicated or reordered"
    fi
else
    record_fail "duplicate explicit path scan failed"
fi

failure_output="${TEST_DIR_TMP}/shellcheck_failure_propagation"
failure_log="${failure_output}/invocations.log"
mkdir -p "${failure_output}"
set +e
run_fake_system_bash \
    "${fake_prefix}" 3 "${failure_output}" "${failure_log}" \
    "${bash_fixture}" "${bootstrap}" \
    > "${failure_output}/runner.log" 2>&1
failure_status="$?"
set -e
if (( failure_status > 1 )) && \
    [[ -s "${failure_output}/summary.json" ]] && \
    [[ "$(json_value "${failure_output}/summary.json" 'value["status"]')" == "2" ]] && \
    [[ "$(invocation_count "${failure_log}" '--shell=bash')" == "1" ]] && \
    [[ "$(invocation_count "${failure_log}" '--shell=sh')" == "1" ]]
then
    record_pass "mixed infrastructure failures propagate after both passes"
else
    record_fail "mixed infrastructure failure propagation changed"
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
