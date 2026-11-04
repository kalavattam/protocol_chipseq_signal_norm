#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_embedded_programs.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="embedded programs and recorded inputs"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


log_ep="${TEST_DIR_LOG}/embedded_programs/findings.tsv"
mkdir -p "$(dirname "${log_ep}")"

print_section "${TEST_NAME}"

# Each blocking finding fails. Each advisory diagnostic with cases left warns
# once with its count, so the migration's progress shows without one warning
# per case.
rc=0
(
    cd "${ROOT_REPO}"
    PYTHONDONTWRITEBYTECODE=1 "${TEST_PYTHON}" -m dev.audit.embedded_programs \
        --root .
) > "${log_ep}" 2>&1 || rc=$?

if (( rc > 1 )); then
    record_fail \
        "embedded-program checker crashed; see $(print_relpath "${log_ep}")"
fi

while IFS=$'\t' read -r sev rule loc msg; do
    case "${sev}" in
        FAIL)
            record_fail "${rule} ${loc}: ${msg}"
            ;;

        COUNT)
            if [[ "${msg}" == "advisory" && "${loc}" != "0" ]]; then
                record_warn "${rule}: ${loc} case(s) left to convert"
            fi
            ;;
    esac
done < "${log_ep}"

if (( rc == 0 )); then
    record_pass "no blocking embedded-program or recorded-input finding"
fi

finish
