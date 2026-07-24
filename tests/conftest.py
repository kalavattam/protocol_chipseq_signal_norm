#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: conftest.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


from __future__ import annotations

import sys
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
for source_root in (
    REPOSITORY_ROOT / "src",
    REPOSITORY_ROOT,
    REPOSITORY_ROOT / "dev" / "audit",
):
    value = str(source_root)
    if value not in sys.path:
        sys.path.insert(0, value)


def write_text(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")
