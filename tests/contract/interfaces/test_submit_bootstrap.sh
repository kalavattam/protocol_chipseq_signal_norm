#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_submit_bootstrap.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="submit shebang and bootstrap policy"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


#  Evaluate the canonical major/minor boundary without mutating Bash internals
function check_version_boundary() {
    local major="${1:-0}"
    local minor="${2:-0}"
    local expected="${3:-reject}"
    local actual="accept"

    if (( major < 4 || ( major == 4 && minor < 4 ) )); then
        actual="reject"
    fi

    if [[ "${actual}" == "${expected}" ]]; then
        record_pass "Bash ${major}.${minor} boundary=${expected}"
    else
        record_fail "Bash ${major}.${minor} boundary expected ${expected}, got ${actual}"
    fi
}


#  Require direct execution of a probe to resolve through one entry shebang
function check_shebang_resolution() {
    local script="${1:-}"
    local probe="${TEST_DIR_TMP}/submit_bootstrap/shebang_probe.sh"
    local actual=""

    {
        sed -n '1p' "${script}"
        # shellcheck disable=SC2016  # Emit a literal probe expression.
        printf '%s\n' 'printf "%s\n" "${BASH}"'
    } > "${probe}"
    chmod 700 "${probe}"

    actual="$(
        PATH="$(dirname "${TEST_BASH}"):/usr/bin:/bin" "${probe}"
    )"
    if [[ "${actual}" == "${TEST_BASH}" ]]; then
        record_pass "direct shebang resolution $(print_relpath "${script}")"
    else
        record_fail \
            "direct shebang resolution $(print_relpath "${script}") expected" \
            "'${TEST_BASH}', got '${actual}'"
    fi
}


dir_tmp="${TEST_DIR_TMP}/submit_bootstrap"
log_audit="${TEST_DIR_LOG}/submit_bootstrap/audit.log"
log_direct="${TEST_DIR_LOG}/submit_bootstrap/source_helper_direct.log"
log_resolve_help="${TEST_DIR_LOG}/submit_bootstrap/source_helper_resolve_help.log"
out_resolve_help="${dir_tmp}/source_helper_resolve_help.stdout"
log_sourced="${TEST_DIR_LOG}/submit_bootstrap/source_helper_sourced.log"
marker="${dir_tmp}/source_count"
helper="${dir_tmp}/counted_helper.sh"


print_section "${TEST_NAME}"

mkdir -p "${dir_tmp}" "$(dirname "${log_audit}")"

# shellcheck disable=SC2016  # Expand in the child authoritative Bash.
if \
    PYTHONDONTWRITEBYTECODE=1 \
    "${TEST_BASH}" -c 'python="$1"; shift; "${python}" "$@"' \
        _ \
        "$(dirname "${TEST_BASH}")/python" \
        -m dev.audit.source_policy \
        --root "${ROOT_REPO}" \
        --submit-bootstrap \
        > "${log_audit}"
then
    record_pass "repository submit-bootstrap audit"
else
    record_fail \
        "repository submit-bootstrap audit; see $(print_relpath "${log_audit}")"
fi

while IFS= read -r script; do
    check_shebang_resolution "${script}"
done < <(find "${ROOT_REPO}/bin" -maxdepth 1 -name 'submit_*.sh' -print | sort)

check_version_boundary 3 2 reject
check_version_boundary 4 3 reject
check_version_boundary 4 4 accept
check_version_boundary 4 9 accept
check_version_boundary 5 0 accept

#  The live macOS system Bash must reach the intended diagnostic before later
#+ Bash >= 4.4 syntax, not fail while parsing that syntax.
set +e
/bin/bash "${ROOT_REPO}/bin/submit_align_fastqs.sh" --help \
    > "${log_direct}" 2>&1
status=$?
set -e
if (( status != 0 )) \
    && grep -Fq "requires Bash >= 4.4" "${log_direct}" \
    && ! grep -Fq "syntax error" "${log_direct}"
then
    record_pass "macOS Bash 3.2 reaches the submit guard before unsupported syntax"
else
    record_fail "macOS Bash 3.2 submit guard ordering"
fi

#  A source-only helper guard returns to its caller under unsupported Bash.
set +e
# shellcheck disable=SC2016  # Expand in the child macOS Bash.
/bin/bash -c '
    source "$1"
    status=$?
    printf "survived:%s\n" "${status}"
    exit 0
' _ "${ROOT_REPO}/lib/bash/core/source_helpers.sh" \
    > "${log_sourced}" 2>&1
status_sourced=$?
/bin/bash "${ROOT_REPO}/lib/bash/core/source_helpers.sh" \
    > "${log_direct}" 2>&1
status_direct=$?
set -e
if [[ "${status_sourced}" -eq 0 ]] \
    && grep -Fq "survived:1" "${log_sourced}"
then
    record_pass "source-only guard returns to its caller"
else
    record_fail "source-only guard did not return to its caller"
fi
if (( status_direct != 0 )) && ! grep -Fq "survived" "${log_direct}"; then
    record_pass "direct source-library invocation exits nonzero"
else
    record_fail "direct source-library invocation did not exit nonzero"
fi

#  Resolver help uses stderr without producing a path or changing the registry.
# shellcheck disable=SC2016  # Expand in the child authoritative Bash.
if "${TEST_BASH}" -c '
    source "$1"
    registry_before="$(declare -p __SOURCED_HELPERS)"
    _source_helper_resolve --help > "$2" 2> "$3"
    status=$?
    registry_after="$(declare -p __SOURCED_HELPERS)"
    (( status == 0 )) && [[ "${registry_before}" == "${registry_after}" ]]
' _ \
    "${ROOT_REPO}/lib/bash/core/source_helpers.sh" \
    "${out_resolve_help}" \
    "${log_resolve_help}" \
    && [[ ! -s "${out_resolve_help}" ]] \
    && grep -Fxq "  _source_helper_resolve" "${log_resolve_help}"
then
    record_pass "source-helper resolver help channel and registry isolation"
else
    record_fail "source-helper resolver help channel or registry isolation"
fi

#  The canonical source registry suppresses a repeated helper bootstrap.
# shellcheck disable=SC2016  # Emit a literal helper expression.
printf '%s\n' 'COUNTED_HELPER_LOADS=$(( ${COUNTED_HELPER_LOADS:-0} + 1 ))' \
    > "${helper}"
# shellcheck disable=SC2016  # Expand in the child authoritative Bash.
"${TEST_BASH}" -c '
    source "$1"
    source_once "$2"
    source_once "$2"
    printf "%s\n" "${COUNTED_HELPER_LOADS}"
' _ "${ROOT_REPO}/lib/bash/core/source_helpers.sh" "${helper}" > "${marker}"
if [[ "$(< "${marker}")" == "1" ]]; then
    record_pass "helper bootstrap runs exactly once"
else
    record_fail "helper bootstrap did not run exactly once"
fi

#  Direct submit help and sourced help-owner behavior remain available.
if "${TEST_BASH}" "${ROOT_REPO}/bin/submit_compute_signal.sh" --help \
    > /dev/null 2>&1
then
    record_pass "direct submit help"
else
    record_fail "direct submit help"
fi
# shellcheck disable=SC2016  # Expand in the child authoritative Bash.
if "${TEST_BASH}" -c '
    source "$1"
    help_submit_compute_signal > /dev/null 2>&1
' _ "${ROOT_REPO}/lib/bash/help/help_submit_compute_signal.sh"
then
    record_pass "sourced submit-help owner"
else
    record_fail "sourced submit-help owner"
fi

#  Compatibility delegators must preserve target help and failure status.
for legacy in submit_filter_bams.sh submit_filter_crams.sh; do
    set +e
    out_legacy="$("${TEST_BASH}" "${ROOT_REPO}/bin/${legacy}" --invalid 2>&1)"
    status_legacy=$?
    out_target="$(
        "${TEST_BASH}" "${ROOT_REPO}/bin/submit_filter_alignments.sh" \
            --invalid 2>&1
    )"
    status_target=$?
    set -e
    if [[ "${status_legacy}" -eq "${status_target}" ]] \
        && [[ "${out_legacy}" == "${out_target}" ]]
    then
        record_pass "${legacy} preserves arguments, output, and status"
    else
        record_fail "${legacy} compatibility delegation drift"
    fi
done

# shellcheck disable=SC2119  # finish intentionally receives no script arguments.
finish
