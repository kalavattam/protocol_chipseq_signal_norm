#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_ai_attribution_gate.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="AI attribution gate"

# Source shared test helpers. The repository root is derived from BASH_SOURCE
# rather than from Git, because Git is not a test-runtime dependency.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/../../.." > /dev/null 2>&1 && pwd
)/tests/support/test_helpers.sh"


# The header rules were applied across the maintained tree, but nothing ran
# the checker over that tree afterward, so drift went unreported until it was
# found by eye. This gate exists to make that impossible: it runs the
# registered command repository-wide and proves the scan was not vacuous.
inspected_floor=200

print_section "${TEST_NAME}"

inventory="${TEST_DIR_OUT}/tmp/ai_attribution_inventory.json"
mkdir -p "$(dirname "${inventory}")"

if run_capture \
    "maintained source attribution" \
    "${TEST_DIR_OUT}/logs/ai_attribution_gate.log" \
    python3 -m dev.audit.ai_attribution \
        --root "${ROOT_REPO}" \
        --current-year 2026 \
        --inventory-output "${inventory}"
then
    record_pass "maintained source headers carry declared attribution"
else
    record_fail "maintained source headers report attribution findings"
fi

# A gate that inspects nothing passes for the wrong reason. Assert the
# inspected count so a discovery regression fails here instead of going quiet.
if [[ -s "${inventory}" ]]; then
    inspected="$(
        python3 -c 'import json,sys; print(len(json.load(open(sys.argv[1]))))' \
            "${inventory}"
    )"
else
    inspected=0
fi

if (( inspected >= inspected_floor )); then
    record_pass "attribution scan inspected ${inspected} maintained source(s)"
else
    record_fail \
        "attribution scan inspected only ${inspected} source(s); a" \
        "changed-file scope cannot report drift from an earlier commit"
fi

finish
