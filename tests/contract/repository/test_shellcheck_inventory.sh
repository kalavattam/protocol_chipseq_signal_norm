#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_shellcheck_inventory.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.6);
# - Anthropic Claude Code (Opus 5, Opus 5.5).
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="ShellCheck inventory provenance"

# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Print the value of a Python expression over a JSON file's 'value'.
function print_json_value() {
    local source="${1}"
    local expression="${2}"

    "${TEST_PYTHON}" \
        "${ROOT_REPO}/tests/support/json_value.py" \
        "${source}" \
        "${expression}"
}


# Print how many files of one language a scan's summary counts.
function print_file_count() {
    print_json_value \
        "${1}/summary.json" \
        "value[\"languages\"][\"${2}\"][\"file_count\"]"
}


# Print how many lines of a file contain a fixed string; a missing file prints
# nothing.
function print_match_count() {
    local fil="${1}"
    local pattern="${2}"

    grep -Fc -- "${pattern}" "${fil}" 2> /dev/null || true
}


# Print how many lines of a file are exactly a given line.
function print_line_count() {
    local fil="${1}"
    local line="${2}"

    grep -Fxc -- "${line}" "${fil}" 2> /dev/null || true
}


# Print the languages that a '--list' output names, sorted and joined by
# spaces. Its provenance lines have no tab and are skipped.
function print_languages() {
    awk -F '\t' 'NF == 2 { print $1 }' "${1}" | sort -u | paste -sd ' ' -
}


# Make a prefix whose 'shellcheck' is the fixture stub and whose 'python' is
# the test interpreter.
function make_fake_prefix() {
    local prefix="${1}"

    mkdir -p "${prefix}/bin"
    cp \
        "${ROOT_REPO}/tests/fixtures/shellcheck/tool/fake_shellcheck.sh" \
        "${prefix}/bin/shellcheck"
    chmod +x "${prefix}/bin/shellcheck"
    ln -sf "${TEST_PYTHON}" "${prefix}/bin/python"
}


# Run a runner under a shell against the fake prefix, with the stub's exit
# status and call log set.
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


# Run the maintained runner under the test Bash.
function run_fake() {
    run_fake_with_runner "${TEST_BASH}" "${runner}" "$@"
}


# Run the maintained runner under the system Bash.
function run_fake_system_bash() {
    run_fake_with_runner /bin/bash "${runner}" "$@"
}


# Print the maintained runner's '--list' output for the given paths, under the
# system Bash and the fake prefix.
function run_list() {
    CONDA_PREFIX="${fake_prefix}" \
    CONDA_DEFAULT_ENV=env_protocol \
        /bin/bash "${runner}" --list -- "$@"
}


# Make a small Git repository that holds the runner, its report module, one
# Bash script in each discovered directory, the POSIX bootstrap, and a tracked
# script deleted from the working tree.
function make_discovery_repo() {
    local root="${1}"

    mkdir -p \
        "${root}/bin" \
        "${root}/dev/audit" \
        "${root}/install/scripts" \
        "${root}/tests/contract"

    cp "${runner}" "${root}/dev/audit/run_shellcheck.sh"
    cp \
        "${ROOT_REPO}/dev/audit/shellcheck_report.py" \
        "${root}/dev/audit/shellcheck_report.py"
    cp \
        "${ROOT_REPO}/tests/fixtures/shellcheck/script/bash.sh" \
        "${root}/bin/discovered.sh"
    cp \
        "${ROOT_REPO}/tests/fixtures/shellcheck/script/posix.sh" \
        "${root}/install/scripts/install_envs_entrypoint.sh"
    cp \
        "${ROOT_REPO}/tests/fixtures/shellcheck/script/bash.sh" \
        "${root}/tests/contract/discovered.sh"
    cp \
        "${ROOT_REPO}/tests/fixtures/shellcheck/script/bash.sh" \
        "${root}/tests/contract/deleted.sh"

    git -C "${root}" init -q
    git -C "${root}" add bin dev install tests
    rm "${root}/tests/contract/deleted.sh"
}


# Write the Bash-only runner that predates the nounset fix.
function make_legacy_runner() {
    local destination="${1}"

    "${TEST_PYTHON}" \
        "${ROOT_REPO}/tests/support/make_legacy_shellcheck_runner.py" \
        "${runner}" \
        "${destination}"
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

if \
    run_capture \
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

if \
    run_capture \
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
bash_fixture="${ROOT_REPO}/tests/fixtures/shellcheck/script/bash.sh"
bootstrap="${ROOT_REPO}/install/scripts/install_envs_entrypoint.sh"
discovery_root="${TEST_DIR_TMP}/shellcheck_discovery"
discovery_log="${TEST_DIR_TMP}/shellcheck_discovery.log"

make_fake_prefix "${fake_prefix}"
make_discovery_repo "${discovery_root}"

# Name each expected discovery entry so the assertion stays scannable. Each
# entry must be a whole line, since 'sh<TAB><path>' is also part of
# 'bash<TAB><path>'.
entry_bin=$'bash\t'"${discovery_root}/bin/discovered.sh"
entry_test=$'bash\t'"${discovery_root}/tests/contract/discovered.sh"
entry_boot=$'sh\t'"${discovery_root}/install/scripts"
entry_boot+="/install_envs_entrypoint.sh"
deleted='tests/contract/deleted.sh'

if \
    CONDA_PREFIX="${fake_prefix}" CONDA_DEFAULT_ENV=env_protocol \
        "${TEST_BASH}" "${discovery_root}/dev/audit/run_shellcheck.sh" \
            --list > "${discovery_log}" 2>&1
then
    if [[
        "$(print_line_count "${discovery_log}" "${entry_bin}")" == "1"
            && "$(print_line_count "${discovery_log}" "${entry_test}")" == "1"
            && "$(print_line_count "${discovery_log}" "${entry_boot}")" == "1"
            && "$(print_match_count "${discovery_log}" "${deleted}")" == "0"
    ]]; then
        record_pass \
            "default discovery assigns maintained Bash and POSIX fixtures"
    else
        record_fail "default discovery omitted or misclassified a fixture"
    fi
else
    record_fail "fixture-backed default discovery failed"
fi

legacy_runner="${TEST_DIR_TMP}/run_shellcheck_legacy.sh"
legacy_output="${TEST_DIR_TMP}/shellcheck_legacy"
legacy_log="${legacy_output}/invocations.log"

rm -rf "${legacy_output}"
mkdir -p "${legacy_output}"

make_legacy_runner "${legacy_runner}"

legacy_status=0
run_fake_with_runner /bin/bash "${legacy_runner}" \
    "${fake_prefix}" 0 "${legacy_output}" "${legacy_log}" "${bash_fixture}" \
    > "${legacy_output}/runner.log" 2>&1 \
    || legacy_status=$?

if [[
    "${legacy_status}" -gt 0
        && "$(
            print_match_count "${legacy_output}/runner.log" 'unbound variable'
        )" -ge 1
]]; then
    record_pass "legacy Bash-only runner exposes the nounset array fault"
else
    record_fail \
        "legacy Bash-only runner did not expose the nounset array fault"
fi

bash_only_output="${TEST_DIR_TMP}/shellcheck_bash_only"
bash_only_log="${bash_only_output}/invocations.log"

rm -rf "${bash_only_output}"
mkdir -p "${bash_only_output}"

if \
    run_fake_system_bash \
        "${fake_prefix}" 0 "${bash_only_output}" "${bash_only_log}" \
        "${bash_fixture}" > "${bash_only_output}/runner.log" 2>&1
then
    if [[
        "$(print_file_count "${bash_only_output}" bash)" == "1"
            && "$(print_file_count "${bash_only_output}" sh)" == "0"
            && "$(print_match_count "${bash_only_log}" '--shell=bash')" == "1"
            && "$(print_match_count "${bash_only_log}" '--shell=sh')" == "0"
            && -s "${bash_only_output}/bash_findings.json"
            && -s "${bash_only_output}/sh_findings.json"
    ]]; then
        record_pass \
            "Bash-only scan preserves empty POSIX evidence under nounset"
    else
        record_fail \
            "Bash-only scan lost language counts, calls, or raw evidence"
    fi
else
    record_fail "Bash-only fixture-backed scan failed"
fi

posix_only_output="${TEST_DIR_TMP}/shellcheck_posix_only"
posix_only_log="${posix_only_output}/invocations.log"

rm -rf "${posix_only_output}"
mkdir -p "${posix_only_output}"

if \
    run_fake_system_bash \
        "${fake_prefix}" 0 "${posix_only_output}" "${posix_only_log}" \
        "${bootstrap}" > "${posix_only_output}/runner.log" 2>&1
then
    if [[
        "$(print_file_count "${posix_only_output}" bash)" == "0"
            && "$(print_file_count "${posix_only_output}" sh)" == "1"
            && "$(print_match_count "${posix_only_log}" '--shell=bash')" == "0"
            && "$(print_match_count "${posix_only_log}" '--shell=sh')" == "1"
            && -s "${posix_only_output}/bash_findings.json"
            && -s "${posix_only_output}/sh_findings.json"
    ]]; then
        record_pass \
            "POSIX-only scan preserves empty Bash evidence under nounset"
    else
        record_fail \
            "POSIX-only scan lost language counts, calls, or raw evidence"
    fi
else
    record_fail "POSIX-only fixture-backed scan failed"
fi

bash_list="${TEST_DIR_TMP}/shellcheck_bash_only.list"
posix_list="${TEST_DIR_TMP}/shellcheck_posix_only.list"

list_status=0
run_list "${bash_fixture}" > "${bash_list}" 2>&1 || list_status=$?
run_list "${bootstrap}" > "${posix_list}" 2>&1 || list_status=$?

if [[
    "${list_status}" -eq 0
        && "$(print_languages "${bash_list}")" == "bash"
        && "$(print_languages "${posix_list}")" == "sh"
]]; then
    record_pass "single-dialect list loops are nounset-safe"
else
    record_fail "single-dialect list loops failed under nounset"
fi

default_output="${TEST_DIR_TMP}/shellcheck_default"
default_log="${default_output}/invocations.log"

rm -rf "${default_output}"
mkdir -p "${default_output}"

if \
    run_fake_with_runner /bin/bash \
        "${discovery_root}/dev/audit/run_shellcheck.sh" "${fake_prefix}" 0 \
        "${default_output}" "${default_log}" \
        > "${default_output}/runner.log" 2>&1
then
    if [[
        "$(print_file_count "${default_output}" bash)" == "2"
            && "$(print_file_count "${default_output}" sh)" == "1"
            && "$(print_match_count "${default_log}" '--shell=bash')" == "1"
            && "$(print_match_count "${default_log}" '--shell=sh')" == "1"
            && "$(print_match_count "${default_log}" "${deleted}")" == "0"
    ]]; then
        record_pass "default scan executes both discovered language passes"
    else
        record_fail "default scan lost discovery or language-pass behavior"
    fi
else
    record_fail "fixture-backed default scan failed"
fi

missing_output="${TEST_DIR_TMP}/shellcheck_missing"
missing_log="${missing_output}/invocations.log"
missing_msg='ShellCheck path is not a file'

rm -rf "${missing_output}"
mkdir -p "${missing_output}"

missing_status=0
run_fake_system_bash \
    "${fake_prefix}" 0 "${missing_output}" "${missing_log}" \
    "${TEST_DIR_TMP}/not_a_shell_file.sh" \
    > "${missing_output}/runner.log" 2>&1 \
    || missing_status=$?

if [[
    "${missing_status}" -eq 2
        && "$(
            print_match_count "${missing_output}/runner.log" "${missing_msg}"
        )" -ge 1
        && ! -e "${missing_log}"
]]; then
    record_pass "missing explicit path remains an infrastructure failure"
else
    record_fail "missing explicit path changed ShellCheck failure behavior"
fi

duplicate_output="${TEST_DIR_TMP}/shellcheck_duplicates"
duplicate_log="${duplicate_output}/invocations.log"

rm -rf "${duplicate_output}"
mkdir -p "${duplicate_output}"

if \
    run_fake_system_bash \
        "${fake_prefix}" 0 "${duplicate_output}" "${duplicate_log}" \
        "${bash_fixture}" "${bash_fixture}" \
        > "${duplicate_output}/runner.log" 2>&1
then
    duplicate_paths="$(
        grep -oF -- "${bash_fixture}" "${duplicate_log}" | wc -l
    )"

    if [[
        "$(print_file_count "${duplicate_output}" bash)" == "2"
            && "${duplicate_paths//[[:space:]]/}" == "2"
            && "$(print_match_count "${duplicate_log}" '--shell=bash')" == "1"
    ]]; then
        record_pass "duplicate explicit paths retain count and order"
    else
        record_fail "duplicate explicit paths were deduplicated or reordered"
    fi
else
    record_fail "duplicate explicit path scan failed"
fi

failure_output="${TEST_DIR_TMP}/shellcheck_failure_propagation"
failure_log="${failure_output}/invocations.log"
failure_summary="${failure_output}/summary.json"

rm -rf "${failure_output}"
mkdir -p "${failure_output}"

failure_status=0
run_fake_system_bash \
    "${fake_prefix}" 3 "${failure_output}" "${failure_log}" \
    "${bash_fixture}" "${bootstrap}" \
    > "${failure_output}/runner.log" 2>&1 \
    || failure_status=$?

if [[
    "${failure_status}" -gt 1
        && -s "${failure_summary}"
        && "$(print_json_value "${failure_summary}" 'value["status"]')" == "2"
        && "$(print_match_count "${failure_log}" '--shell=bash')" == "1"
        && "$(print_match_count "${failure_log}" '--shell=sh')" == "1"
]]; then
    record_pass "mixed infrastructure failures propagate after both passes"
else
    record_fail "mixed infrastructure failure propagation changed"
fi

for expected_status in 0 1 3; do
    output="${TEST_DIR_TMP}/shellcheck_status_${expected_status}"
    invocation_log="${output}/invocations.log"

    mkdir -p "${output}"

    observed_status=0
    run_fake \
        "${fake_prefix}" "${expected_status}" "${output}" \
        "${invocation_log}" "${bash_fixture}" \
        > "${output}/runner.log" 2>&1 \
        || observed_status=$?

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

if [[
    "$(print_json_value "${status_one_summary}" 'value["finding_count"]')" \
        == "1"
]]; then
    record_pass "status-1 structured evidence retains its finding count"
else
    record_fail "status-1 structured evidence has the wrong finding count"
fi

split_output="${TEST_DIR_TMP}/shellcheck_split"
split_log="${split_output}/invocations.log"

rm -rf "${split_output}"
mkdir -p "${split_output}"

if \
    run_fake \
        "${fake_prefix}" 0 "${split_output}" "${split_log}" \
        "${bash_fixture}" "${bootstrap}" \
        > "${split_output}/runner.log" 2>&1
then
    if [[
        "$(print_match_count "${split_log}" '--shell=bash')" -ge 1
            && "$(print_match_count "${split_log}" '--shell=sh')" -ge 1
            && "$(print_file_count "${split_output}" bash)" == "1"
            && "$(print_file_count "${split_output}" sh)" == "1"
    ]]; then
        record_pass "Bash and the POSIX bootstrap use separate language modes"
    else
        record_fail "ShellCheck language-mode discovery was not split"
    fi
else
    record_fail "fixture-backed language split failed"
fi

if [[ "$(print_match_count "${split_log}" "-x -P ${ROOT_REPO}")" -ge 1 ]]; then
    record_pass \
        "ShellCheck retains external sources and repository source path"
else
    record_fail "ShellCheck invocation lost -x or the repository source path"
fi


finish
