#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_python_help_format.py
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


"""
Test bounded Python help-literal formatting.
"""

from __future__ import annotations

import ast
from pathlib import Path

from dev.tools.python_help_format import format_source

ROOT = Path(__file__).resolve().parents[3]
FIXTURES = ROOT / "tests/fixtures/python_source_policy"


def rendered(source: str) -> str:
    tree = ast.parse(source)
    keyword = next(
        keyword
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        for keyword in node.keywords
        if keyword.arg == "help"
    )
    return ast.literal_eval(keyword.value)


def test_formatter_is_greedy_value_preserving_and_idempotent() -> None:
    source = (FIXTURES / "format/help_input.py").read_text()
    expected = (FIXTURES / "format/help_expected.py").read_text()
    actual = format_source(source)

    assert actual == expected
    assert rendered(actual) == rendered(source)
    assert format_source(actual) == actual


def test_formatter_preserves_unicode_literals_and_packs_them_greedily() -> (
    None
):
    source = (FIXTURES / "format/help_unicode_input.py").read_text()
    expected = (
        FIXTURES / "format/help_unicode_expected.py"
    ).read_text()
    actual = format_source(source)

    assert actual == expected
    assert rendered(actual) == rendered(source)
    assert format_source(actual) == actual


def test_accepted_pilot_help_literals_are_already_canonical() -> None:
    for name in ("compute_input_floor.py", "compute_pseudo.py"):
        source = (
            ROOT / "src/protocol_chipseq_signal_norm/cli" / name
        ).read_text(encoding="utf-8")

        assert format_source(source) == source
