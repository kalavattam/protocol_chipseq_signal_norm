#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_source_policy.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Focused regressions for bounded ShellCheck source policy.
"""

from __future__ import annotations

import unittest
from pathlib import Path

from dev.audit.source_policy import (
    CANONICAL_BASH_SHEBANG,
    bootstrap_help_spacing_warnings,
    classify_source,
    discover_submit_interfaces,
    submit_bootstrap_findings,
    suppression_mismatches,
)


class SourcePolicyTest(unittest.TestCase):
    """
    Keep source mapping and bootstrap spacing intentionally bounded.
    """

    def test_static_source_requires_exact_repository_relative_mapping(
        self,
    ) -> None:
        record = classify_source(
            "lib/bash/core/check_args.sh",
            'source "${_dir_src_args}/source_helpers.sh"',
        )

        self.assertEqual(record.classification, "static")
        self.assertEqual(record.target, "lib/bash/core/source_helpers.sh")
        self.assertEqual(record.required_suppression, None)

    def test_runtime_selected_source_retains_sc1090(self) -> None:
        record = classify_source(
            "bin/execute_align_fastqs.sh",
            'source "${fnc_src}"',
        )

        self.assertEqual(record.classification, "runtime_selected")
        self.assertEqual(record.required_suppression, "SC1090")

    def test_conda_activation_source_retains_sc1091(self) -> None:
        record = classify_source(
            "lib/bash/core/handle_env.sh",
            'source "$(conda info --base)/etc/profile.d/conda.sh"',
        )

        self.assertEqual(record.classification, "dynamic_activation")
        self.assertEqual(record.required_suppression, "SC1091")

    def test_semantic_test_path_maps_shared_test_helpers(self) -> None:
        record = classify_source(
            "tests/contract/repository/test_shell_syntax.sh",
            (
                'source "$(git rev-parse '
                '--show-toplevel)/tests/support/test_helpers.sh"'
            ),
        )

        self.assertEqual(record.classification, "static")
        self.assertEqual(record.target, "tests/support/test_helpers.sh")

    def test_fixture_recipe_maps_shared_fixture_helpers(self) -> None:
        record = classify_source(
            "tests/fixtures/compute_signal/make.sh",
            (
                'source "$(git rev-parse '
                '--show-toplevel)/tests/support/fixture_helpers.sh"'
            ),
        )

        self.assertEqual(record.classification, "static")
        self.assertEqual(record.target, "tests/support/fixture_helpers.sh")

    def test_stale_sc1090_for_static_source_is_reported(self) -> None:
        findings = suppression_mismatches(
            "bin/submit_align_fastqs.sh",
            "# shellcheck disable=SC1090\n"
            "source "
            '"${_dir_scr_help}/functions/help/help_submit_align_fastqs.sh'
            '"\n',
        )

        self.assertEqual(len(findings), 2)
        self.assertTrue(any("SC1090" in item.message for item in findings))
        self.assertTrue(any("source=" in item.message for item in findings))

    def test_bootstrap_help_source_requires_one_blank_before_directive(
        self,
    ) -> None:
        text = (
            '_dir_scr_help="$(pwd)"\n'
            "# shellcheck source=lib/bash/help/help_submit_align_fastqs.sh\n"
            "source "
            '"${_dir_scr_help}/functions/help/help_submit_align_fastqs.sh'
            '"\n'
        )

        findings = bootstrap_help_spacing_warnings(
            "bin/submit_align_fastqs.sh",
            text,
        )

        self.assertEqual(len(findings), 1)
        self.assertIn("exactly one blank", findings[0].message)

    def test_bootstrap_help_source_with_one_blank_passes(self) -> None:
        text = (
            '_dir_scr_help="$(pwd)"\n\n'
            "# shellcheck source=lib/bash/help/help_submit_align_fastqs.sh\n"
            "source "
            '"${_dir_scr_help}/functions/help/help_submit_align_fastqs.sh'
            '"\n'
        )

        self.assertEqual(
            bootstrap_help_spacing_warnings(
                "bin/submit_align_fastqs.sh",
                text,
            ),
            [],
        )

    def test_submit_inventory_is_source_derived_and_unambiguous(self) -> None:
        root = Path(__file__).resolve().parents[3]

        inventory = discover_submit_interfaces(root)

        self.assertTrue(inventory)
        self.assertEqual(
            {row.classification for row in inventory},
            {"compatibility delegator", "worker-only interface"},
        )
        self.assertTrue(all(not row.source_callers for row in inventory))
        self.assertTrue(all(not row.executable_bit for row in inventory))

    def test_submit_bootstrap_repository_policy_is_green(self) -> None:
        root = Path(__file__).resolve().parents[3]

        findings, inventory = submit_bootstrap_findings(root)

        self.assertEqual(findings, [])
        self.assertTrue(
            all(
                row.current_shebang == CANONICAL_BASH_SHEBANG
                for row in inventory
            ),
        )


if __name__ == "__main__":
    unittest.main()
