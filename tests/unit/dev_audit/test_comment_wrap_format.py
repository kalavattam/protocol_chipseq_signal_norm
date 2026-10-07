#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_comment_wrap_format.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Test the Shell and Python comment-wrap formatter.
"""

from __future__ import annotations

import ast
import shutil
from pathlib import Path

from dev.audit.shell_source_form import check_text
from dev.tools.comment_wrap_format import format_text, main

ROOT = Path(__file__).resolve().parents[3]
SHELL = ROOT / "tests" / "fixtures" / "shell_source_form"
PYTHON = ROOT / "tests" / "fixtures" / "python_source_policy"


def wrap_lines(text: str) -> list[int]:
    """
    Return the lines that the Shell checker's wrap facet reports.
    """

    return [
        finding.line
        for finding in check_text(text, "tests/fixture.sh")
        if finding.message.startswith("comment breaks before")
    ]


def test_shell_pair_is_greedy_idempotent_and_leaves_review_lines() -> None:
    """
    Refill the Shell input to the expected output and leave only review lines.

    The lines the formatter leaves for review are exactly the lines the checker
    still reports in the expected output: the hyphenated paragraph.
    """

    source = (SHELL / "format/comment_wrap_input.sh").read_text()
    expected = (SHELL / "format/comment_wrap_expected.sh").read_text()
    actual, review = format_text(source, "shell")

    assert actual == expected
    assert format_text(actual, "shell")[0] == actual
    assert review == wrap_lines(expected)
    assert review


def test_python_pair_is_greedy_idempotent_and_code_preserving() -> None:
    """
    Refill the Python input to the expected output and change no code.

    Comments are not part of the syntax tree, so equal trees prove that only
    comments moved, and the comment-like lines in the string stayed as written.
    """

    source = (PYTHON / "format/comment_wrap_input.py").read_text()
    expected = (PYTHON / "format/comment_wrap_expected.py").read_text()
    actual, review = format_text(source, "python")

    assert actual == expected
    assert format_text(actual, "python")[0] == actual
    assert review == []
    assert ast.dump(ast.parse(actual)) == ast.dump(ast.parse(source))


def test_main_previews_by_default_and_writes_only_when_asked(
    tmp_path: Path,
) -> None:
    """
    Report a file that would change, and rewrite it only with '--write'.
    """

    target = tmp_path / "comment_wrap_input.sh"
    shutil.copyfile(SHELL / "format/comment_wrap_input.sh", target)
    original = target.read_text()
    expected = (SHELL / "format/comment_wrap_expected.sh").read_text()

    assert main([str(target)]) == 1
    assert target.read_text() == original

    assert main([str(target), "--write"]) == 0
    assert target.read_text() == expected
