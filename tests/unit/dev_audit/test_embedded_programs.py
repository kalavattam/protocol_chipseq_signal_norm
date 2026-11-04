#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_embedded_programs.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Test the embedded-program and record-input checker on its fixtures.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from dev.audit.embedded_programs import (
    BLOCKING,
    main,
    next_quote_state,
    scan_repository,
)

ROOT = Path(__file__).resolve().parents[3]
FIXTURES = ROOT / "tests/fixtures/embedded_programs"

# Each rejected case: its diagnostic, its file, and its line.
REJECTED = {
    ("SHELL.EMBED.PYTHON", "bin/long_programs.sh", 2),
    ("SHELL.EMBED.PYTHON", "bin/long_programs.sh", 5),
    ("SHELL.EMBED.BASH_C", "bin/long_programs.sh", 7),
    ("SHELL.EMBED.AWK", "bin/long_programs.sh", 9),
    ("SHELL.EMBED.SED", "bin/long_programs.sh", 15),
    ("SHELL.EMBED.SED_GNU", "bin/long_programs.sh", 19),
    ("SHELL.EMBED.PERL", "bin/long_programs.sh", 20),
    ("TEST.FIXTURES.RECORD_INPUT", "tests/contract/test_writes.sh", 2),
    ("TEST.FIXTURES.RECORD_INPUT", "tests/contract/test_writes.sh", 3),
    ("TEST.FIXTURES.PERTURBATION_REASON", "tests/contract/test_writes.sh", 4),
    ("TEST.FIXTURES.RECORD_INPUT", "tests/contract/test_writes.sh", 5),
    ("TEST.FIXTURES.RECORD_INPUT", "tests/contract/test_writes.sh", 6),
    ("PY.EMBED.STRING_INPUT", "tests/unit/test_literal.py", 5),
    ("PY.EMBED.STRING_INPUT", "tests/unit/test_literal.py", 15),
    ("TEST.FIXTURES.RECORD_INPUT", "tests/unit/test_literal.py", 10),
    ("TEST.FIXTURES.INLINE_STUB", "tests/fixtures/demo/make.sh", 7),
    ("TEST.FIXTURES.INLINE_STUB", "tests/fixtures/demo/make.sh", 9),
}


def found_at(verdict: str) -> set[tuple[str, str, int]]:
    """
    Return the findings in one verdict's fixture root.
    """

    return {
        (item.rule_id, item.path, item.line)
        for item in scan_repository(FIXTURES / verdict)
    }


def test_rejected_cases_are_each_found_once() -> None:
    found = scan_repository(FIXTURES / "rejected")

    assert {(item.rule_id, item.path, item.line) for item in found} == (
        REJECTED
    )
    assert len(found) == len(REJECTED)


@pytest.mark.parametrize("verdict", ("accepted", "boundary", "nonapplicable"))
def test_conforming_roots_have_no_findings(verdict: str) -> None:
    assert found_at(verdict) == set()


def test_main_exits_one_only_on_a_blocking_finding(
    capsys: pytest.CaptureFixture[str],
) -> None:
    assert main(["--root", str(FIXTURES / "rejected")]) == 1

    capsys.readouterr()

    assert main(["--root", str(FIXTURES / "accepted")]) == 0

    out = capsys.readouterr().out

    assert "COUNT\tSHELL.EMBED.PYTHON\t0\tblocking" in out
    assert "COUNT\tSHELL.EMBED.BASH_C\t0\tadvisory" in out


def test_json_output_counts_each_diagnostic(
    capsys: pytest.CaptureFixture[str],
) -> None:
    main(["--root", str(FIXTURES / "rejected"), "--json"])

    report = json.loads(capsys.readouterr().out)

    assert report["counts"]["TEST.FIXTURES.RECORD_INPUT"] == 5
    assert report["counts"]["SHELL.EMBED.PYTHON"] == 2
    assert set(report["blocking"]) == set(BLOCKING)
    assert len(report["findings"]) == len(REJECTED)


def test_quote_state_reopens_commands_inside_substitution() -> None:
    assert next_quote_state('out="$(', []) == ['"', "("]
    assert next_quote_state(')"', ['"', "("]) == []
    assert next_quote_state("exp='a", []) == ["'"]
    assert next_quote_state("echo $'it\\'s'", []) == []
    assert next_quote_state("echo x  # it's a comment", []) == []
