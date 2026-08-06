#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_python_boundaries.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1

TEST_NAME="Python package boundaries"
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"
cd "${ROOT_REPO}"
print_section "${TEST_NAME}"

python - << PY
from __future__ import annotations

import ast
from pathlib import Path

root = Path("src/protocol_chipseq_signal_norm")
findings: list[str] = []
for path in sorted(root.rglob("*.py")):
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    role = path.relative_to(root).parts[0]
    for node in ast.walk(tree):
        names: list[str] = []
        if isinstance(node, ast.Import):
            names = [alias.name for alias in node.names]
        elif isinstance(node, ast.ImportFrom):
            names = [node.module or ""]
        for name in names:
            if name == "scripts" or name.startswith("scripts."):
                findings.append(f"{path}:{node.lineno}: retired scripts import")
            if role == "utilities" and ".cli" in name:
                findings.append(f"{path}:{node.lineno}: utility imports CLI layer")
            if role == "cli" and name.startswith(
                "protocol_chipseq_signal_norm.cli."
            ):
                findings.append(f"{path}:{node.lineno}: CLI imports peer CLI")

if findings:
    raise SystemExit("\n".join(findings))
PY
record_pass "package import directions avoid retired and inverted layers"
finish
