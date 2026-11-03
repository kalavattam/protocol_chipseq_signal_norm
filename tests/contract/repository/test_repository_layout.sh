#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_repository_layout.sh
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

TEST_NAME="repository layout"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


cd "${ROOT_REPO}"

print_section "${TEST_NAME}"

for path in \
    bin \
    lib/bash/core \
    lib/bash/dispatch \
    lib/bash/workflows \
    lib/bash/help \
    src/protocol_chipseq_signal_norm/cli \
    src/protocol_chipseq_signal_norm/utilities \
    docs/standards
do
    test -d "${path}"
done

for retired in \
    scripts \
    analysis \
    scraps \
    protocol_chipseq_signal_norm \
    lib/bash/blog \
    docs/style
do
    test ! -e "${retired}"
done

python "${ROOT_REPO}/tests/support/check_python_package.py" \
    packages \
    "${ROOT_REPO}"

record_pass "required roots, retired roots, and package discovery conform"

finish
