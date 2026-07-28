#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_python_version_policy.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


"""
Focused regressions for the repository-wide Python version policy.
"""

from __future__ import annotations

import unittest
from pathlib import Path

from dev.audit.python_version_policy import (
    MINIMUM_PYTHON,
    RULE_ID,
    environment_python_version,
    is_supported_version,
)


class PythonVersionPolicyTest(unittest.TestCase):
    """
    Require one source-derived Python >= 3.11 policy architecture.
    """

    def test_authoritative_minimum_is_python_3_11(self) -> None:
        self.assertEqual(MINIMUM_PYTHON, (3, 11))
        self.assertEqual(RULE_ID, "PY.VERSION.FLOOR")

    def test_supported_version_boundary(self) -> None:
        self.assertFalse(is_supported_version((3, 10)))
        self.assertTrue(is_supported_version((3, 11)))
        self.assertTrue(is_supported_version((3, 99)))
        self.assertTrue(is_supported_version((4, 0)))

    def test_environment_requirement_accepts_supported_newer_versions(
        self,
    ) -> None:
        self.assertEqual(
            environment_python_version("dependencies:\n  - python=3.11\n"),
            (3, 11),
        )
        self.assertEqual(
            environment_python_version("dependencies:\n  - python>=3.12.2\n"),
            (3, 12),
        )
        self.assertEqual(
            environment_python_version("dependencies:\n  - python=3.10\n"),
            (3, 10),
        )
        self.assertIsNone(
            environment_python_version("dependencies:\n  - pip\n"),
        )

    def test_floor_parser_rejects_syntax_newer_than_python_3_11(self) -> None:
        from dev.audit.python_version_policy import parse_floor_syntax

        with self.assertRaises(SyntaxError):
            parse_floor_syntax("type Alias = int\n", "python312_only.py")

        parse_floor_syntax("Alias = int\n", "python311.py")

    def test_inventory_classifies_maintained_entry_points_and_modules(
        self,
    ) -> None:
        from dev.audit.python_version_policy import inventory_repository

        root = Path(__file__).resolve().parents[3]

        inventory = inventory_repository(root)

        self.assertTrue(inventory)
        self.assertTrue(any(row.role == "entry_point" for row in inventory))
        self.assertTrue(
            any(row.role == "imported_module" for row in inventory),
        )
        self.assertTrue(
            all(row.minimum == "Python >= 3.11" for row in inventory),
        )

    def test_guard_policy_rejects_python_3_10_and_accepts_python_3_11(
        self,
    ) -> None:
        from dev.audit.python_version_policy import guard_findings, local_guard

        stale = (
            'assert sys.version_info >= (3, 10), "Python >= 3.10 required."\n'
        )
        current = (
            'assert sys.version_info >= (3, 11), "Python >= 3.11 required."\n'
        )

        self.assertTrue(guard_findings("stale.py", stale))
        self.assertFalse(guard_findings("current.py", current))
        self.assertEqual(local_guard(stale), "Python >= 3.10")
        self.assertEqual(local_guard(current), "Python >= 3.11")
        self.assertEqual(
            local_guard(f"STALE_FIXTURE = {stale!r}\n"),
            "none; covered centrally",
        )

    def test_repository_parity_has_one_stable_rule_id(self) -> None:
        from dev.audit.python_version_policy import scan_repository

        root = Path(__file__).resolve().parents[3]

        findings, inventory = scan_repository(root)

        self.assertTrue(inventory)
        self.assertEqual(
            {item.rule_id for item in findings},
            set() if not findings else {"PY.VERSION.FLOOR"},
        )


if __name__ == "__main__":
    unittest.main()
