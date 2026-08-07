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
# - Anthropic Claude Code (Opus 5).
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="repository layout"
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"
cd "${ROOT_REPO}"
print_section "${TEST_NAME}"

for path in bin lib/bash/core lib/bash/dispatch lib/bash/workflows \
    lib/bash/help src/protocol_chipseq_signal_norm/cli \
    src/protocol_chipseq_signal_norm/utilities docs/standards; do
    test -d "${path}"
done

for retired in scripts analysis scraps protocol_chipseq_signal_norm; do
    test ! -e "${retired}"
done
test ! -e lib/bash/blog
test ! -e docs/style

python - << PY
from setuptools.discovery import PackageFinder

found = set(PackageFinder.find("src"))
expected = {
    "protocol_chipseq_signal_norm",
    "protocol_chipseq_signal_norm.cli",
    "protocol_chipseq_signal_norm.utilities",
}
raise SystemExit(0 if found == expected else f"unexpected packages: {sorted(found)}")
PY
record_pass "required roots, retired roots, and package discovery conform"
finish
