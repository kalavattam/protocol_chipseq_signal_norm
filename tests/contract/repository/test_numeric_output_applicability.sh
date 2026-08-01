#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_numeric_output_applicability.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Distributed under the MIT license.

set -euo pipefail

TEST_NAME="numeric output applicability"

# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


print_section "${TEST_NAME}"

inventory="$(mktemp "${TMPDIR:-/tmp}/numeric_inventory.XXXXXX.json")"
trap 'rm -f "${inventory}"' EXIT

env PYTHONDONTWRITEBYTECODE=1 conda run -n env_protocol python \
    -m dev.audit.numeric_emission_inventory --root "${ROOT_REPO}" \
    --json-out "${inventory}"

env PYTHONDONTWRITEBYTECODE=1 conda run -n env_protocol python \
    -m dev.audit.numeric_output_applicability --root "${ROOT_REPO}" \
    --inventory "${inventory}"
record_pass "numeric applicability evidence and readiness contract pass"

finish
