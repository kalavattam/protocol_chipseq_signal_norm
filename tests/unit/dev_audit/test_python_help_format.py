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


def _table_fixture(*rows: str) -> str:
    """
    Build one help group whose paragraph is a pipe table.
    """

    literals = "\n".join(f'            "{row}\\n"' for row in rows)

    return f'''"""
Provide one table fixture.
"""

import argparse


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse fixture arguments.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--value",
        help=(
            "Lead in.\\n"
            "\\n"
{literals}
        ),
    )
    return parser.parse_args(argv)
'''


def test_formatter_leaves_a_table_row_intact() -> None:
    """
    The approved formatter must not undo the form the standard mandates.

    'PY.CLI.HELP.LAYOUT' gives each table row its own literal. Reflowing a row
    across two literals preserves the rendered value but destroys the source
    table, so the formatter would silently fight the rule it realizes.
    """

    header = "| name" + (" " * 55) + "| value |"
    delimiter = "| :---" + (" " * 55) + "| :---  |"
    source = _table_fixture(header, delimiter)
    row_line = next(
        line for line in source.splitlines() if "| name" in line
    )

    assert len(row_line) > 79
    assert format_source(source) == source
    assert rendered(format_source(source)) == rendered(source)


def test_formatter_still_reflows_a_lone_pipe_line() -> None:
    """
    One row is not a table, so it keeps the ordinary greedy treatment.
    """

    lone = "| name" + (" " * 55) + "| value |"
    source = _table_fixture(lone)

    assert format_source(source) != source
    assert rendered(format_source(source)) == rendered(source)


def test_formatter_keeps_a_quoted_span_whole() -> None:
    """
    The formatter splits on the same units the checker measures.

    Splitting a quoted span to satisfy the width budget would emit source the
    checker accepts while the expression reads as broken, so the two must
    agree on what a unit is.
    """

    source = _table_fixture(
        "Filter values by threshold and report each retained bin",
        "'|x| >= eps' trailing.",
    )
    formatted = format_source(source)

    assert "'|x| >= eps'" in source
    assert "'|x| >= eps'" in formatted
    assert rendered(formatted) == rendered(source)
