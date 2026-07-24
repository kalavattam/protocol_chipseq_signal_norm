#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_check_args_parser.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


"""Focused parser regression for 'check_arg_supplied' assignment aliases."""

from __future__ import annotations

import subprocess
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
CHECK_ARGS = REPO_ROOT / "lib/bash/core/check_args.sh"


class CheckArgSuppliedParserTest(unittest.TestCase):
    """Keep the documented assignment spelling synchronized with parsing."""

    def call(self, *arguments: str) -> subprocess.CompletedProcess[str]:
        command = 'source "$1"; shift; check_arg_supplied "$@"'
        return subprocess.run(
            ["bash", "-c", command, "check-arg-supplied", str(CHECK_ARGS), *arguments],
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
        )

    def test_asgmt_is_accepted(self) -> None:
        result = self.call("--asgmt", "value", "--name", "required_arg")
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_asmgt_transposition_is_rejected(self) -> None:
        result = self.call("--asmgt", "value", "--name", "required_arg")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(
            "unknown option/parameter passed: '--asmgt'.",
            result.stderr,
        )

    def test_help_uses_asgmt_and_omits_transposition(self) -> None:
        result = self.call("--help")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--asgmt", result.stderr)
        self.assertNotIn("--asmgt", result.stderr)


if __name__ == "__main__":
    unittest.main()
