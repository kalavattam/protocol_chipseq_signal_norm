#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_runner_selection.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="safe runner selection"
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"
cd "${ROOT_REPO}"
print_section "${TEST_NAME}"

canonical_artifact_root="${ROOT_REPO}/artifacts/tests"
case "${TEST_ARTIFACT_ROOT}" in
    "${canonical_artifact_root}")
        record_pass "runner uses the canonical repository artifact root"
        ;;

    "${ROOT_REPO}"|"${ROOT_REPO}"/*)
        record_fail "runner uses an arbitrary in-repository artifact root"
        ;;

    /*)
        record_pass "runner uses an approved absolute external artifact root"
        ;;

    *)
        record_fail "runner artifact root is not absolute"
        ;;
esac
if [[ "${PYTHONPYCACHEPREFIX}" == "${TEST_ARTIFACT_ROOT}/pycache" ]]; then
    record_pass "runner derives bytecode state from the exact artifact root"
else
    record_fail "runner bytecode state diverges from the artifact root"
fi
if [[ "${PYTEST_ADDOPTS}" == \
    *"cache_dir=${TEST_ARTIFACT_ROOT}/pytest_cache"* ]]
then
    record_pass "runner derives pytest state from the exact artifact root"
else
    record_fail "runner pytest state diverges from the artifact root"
fi

selection="$(bash tests/run_tests.sh --list all-safe)"
[[ -n "${selection}" ]] || record_fail "all-safe selection is empty"
default_selection="$(bash tests/run_tests.sh --list)"
if [[ "${default_selection}" == "${selection}" ]]; then
    record_pass "default selection equals explicit all-safe"
else
    record_fail "default selection differs from explicit all-safe"
fi
if [[ -n "$(printf '%s\n' "${selection}" | sort | uniq -d)" ]]; then
    record_fail "all-safe selection contains duplicate paths"
else
    record_pass "all-safe selection is unique"
fi

if grep -Fqx 'tests/integration/slurm/run_wet_tests.sh' <<< "${selection}"; then
    record_fail "safe selection includes the wet Slurm coordinator"
else
    record_pass "wet Slurm coordinator is excluded"
fi

overlap="$(bash tests/run_tests.sh --list all-safe contract all-safe)"
if [[ -n "$(printf '%s\n' "${overlap}" | sort | uniq -d)" ]]; then
    record_fail "overlapping groups contain duplicate paths"
else
    record_pass "overlapping groups are deduplicated"
fi

for group in \
    unit \
    contract \
    integration-local \
    integration-parallel \
    integration-slurm
do
    if bash tests/run_tests.sh --list "${group}" > /dev/null; then
        record_pass "group lists independently: ${group}"
    else
        record_fail "group failed independent listing: ${group}"
    fi
done

if RUN_PARALLEL=invalid RUN_SLURM=invalid \
    bash tests/run_tests.sh --list all-safe > /dev/null
then
    record_pass "list mode does not evaluate optional gates"
else
    record_fail "list mode evaluated an optional gate"
fi

runner_source="${ROOT_REPO}/tests/run_tests.sh"
download_gate_line="$(
    grep -n 'if gate_enabled RUN_DOWNLOAD; then' "${runner_source}" \
        | cut -d: -f1
)"
download_fixture_line="$(
    grep -n 'ensure_fixture download_fastqs' "${runner_source}" \
        | cut -d: -f1
)"
if [[ \
    "${download_gate_line}" =~ ^[0-9]+$ \
    && "${download_fixture_line}" =~ ^[0-9]+$ \
    && "${download_gate_line}" -lt "${download_fixture_line}" \
]]; then
    record_pass "download fixture preparation follows the explicit gate"
else
    record_fail "download fixture preparation is not explicitly gated"
fi

invalid_log="${TEST_DIR_TMP}/runner_invalid_group.log"
if bash tests/run_tests.sh --list contract invalid-group \
    > "${invalid_log}" 2>&1
then
    record_fail "runner accepted an invalid group"
elif grep -Eq '^(unit|tests/)' "${invalid_log}"; then
    record_fail "runner discovered tests before rejecting an invalid group"
else
    record_pass "runner rejects invalid groups before discovery"
fi

before="$(find artifacts/tests -type f -print 2> /dev/null | sort)"
unset TEST_ARTIFACT_ROOT
bash tests/run_tests.sh --list unit > /dev/null
after="$(find artifacts/tests -type f -print 2> /dev/null | sort)"
if [[ "${before}" == "${after}" ]]; then
    record_pass "list mode does not create repository evidence"
else
    record_fail "list mode changed repository evidence"
fi

if TEST_ARTIFACT_ROOT="${ROOT_REPO}/invalid-test-root" \
    bash tests/run_tests.sh unit > /dev/null 2>&1
then
    record_fail "runner accepted an in-repository artifact override"
else
    record_pass "runner rejects in-repository artifact overrides"
fi

if TEST_ARTIFACT_ROOT="" bash tests/run_tests.sh unit > /dev/null 2>&1; then
    record_fail "runner accepted an empty artifact override"
else
    record_pass "runner rejects an empty artifact override"
fi

artifact_link_parent="$(
    mktemp -d "${TEST_DIR_TMP}/artifact_root_link.XXXXXX"
)"
artifact_link="${artifact_link_parent}/root"
ln -s "${ROOT_REPO}" "${artifact_link}"
if TEST_ARTIFACT_ROOT="${artifact_link}" \
    bash tests/run_tests.sh unit > /dev/null 2>&1
then
    record_fail "runner accepted a symlink artifact override"
else
    record_pass "runner rejects a symlink artifact override"
fi

finish
