#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_boolean_contracts.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Focused regressions for Boolean-like contract enforcement.
"""

from __future__ import annotations

import unittest
from collections import Counter
from pathlib import Path

from dev.audit.boolean_contracts import (
    FALSE_TOKENS,
    TRUE_TOKENS,
    inventory_repository,
    scan_repository,
)

REPO_ROOT = Path(__file__).resolve().parents[3]


class BooleanContractsTest(unittest.TestCase):
    """
    Require one strict policy without changing presence-only flags.
    """

    def test_canonical_token_sets_are_exact_and_disjoint(self) -> None:
        self.assertEqual(TRUE_TOKENS, ("true", "t", "yes", "y", "1"))
        self.assertEqual(FALSE_TOKENS, ("false", "f", "no", "n", "0"))
        self.assertFalse(set(TRUE_TOKENS) & set(FALSE_TOKENS))

    def test_inventory_is_deterministic_and_fully_dispositioned(self) -> None:
        first = inventory_repository(REPO_ROOT)
        second = inventory_repository(REPO_ROOT)

        self.assertEqual(first, second)
        self.assertTrue(first)
        self.assertNotIn(
            "ambiguous",
            Counter(row.classification for row in first),
        )

        identities = [
            (row.owner_identity, row.line, row.name, row.classification)
            for row in first
        ]

        self.assertEqual(len(identities), len(set(identities)))

    def test_presence_flags_remain_value_free(self) -> None:
        flags = [
            row
            for row in inventory_repository(REPO_ROOT)
            if row.classification == "presence-only flag"
        ]

        self.assertTrue(flags)
        self.assertTrue(
            all(row.value_form == "no-value flag" for row in flags),
        )
        self.assertTrue(all(not row.accepted_true_like for row in flags))
        self.assertTrue(all(not row.accepted_false_like for row in flags))

    def test_non_boolean_modes_are_not_normalized(self) -> None:
        modes = [
            row
            for row in inventory_repository(REPO_ROOT)
            if row.classification == "non-Boolean enum or mode"
        ]

        self.assertTrue(modes)
        self.assertTrue(
            all(row.normalization_helper == "none" for row in modes),
        )

    def test_current_repository_has_no_boolean_contract_findings(self) -> None:
        findings, inventory = scan_repository(REPO_ROOT)

        self.assertTrue(inventory)
        self.assertEqual(
            [],
            findings,
            "\n".join(row.format() for row in findings),
        )


if __name__ == "__main__":
    unittest.main()
