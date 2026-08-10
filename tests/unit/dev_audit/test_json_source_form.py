#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_json_source_form.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Focused regressions for the canonical JSON rendering.
"""

from __future__ import annotations

import json
import unittest
from pathlib import Path

from dev.audit.json_source_form import (
    BUDGET,
    canonical,
    check_text,
    maintained_paths,
    render,
)

ROOT = Path(__file__).resolve().parents[3]
FIXTURES = ROOT / "tests" / "fixtures" / "json_source_form" / "source"

# Every generated fixture, with the finding class it exists to exercise. A
# fixture is consumed by hard failure rather than by 'skipTest', so a broken
# generation step fails loudly instead of turning the suite green.
FIXTURE_CLASSES = (
    ("canonical.json", None),
    ("inline_overflow.json", "past the 79-column budget"),
    ("expanded_fits.json", "fits the 79-column budget inline"),
    ("hybrid_delimiter.json", "delimiter of an expanded structure"),
    ("wrong_indent.json", "indentation is"),
    ("tab_indent.json", "line contains a tab"),
    ("no_trailing_newline.json", "must end with exactly one newline"),
    ("duplicate_key.json", "unreadable JSON"),
    ("unreadable.json", "unreadable JSON"),
)


def fixture_text(name: str) -> str:
    """
    Read one generated fixture, failing loudly when it is absent.
    """

    path = FIXTURES / name

    if not path.is_file():
        raise AssertionError(
            f"missing generated fixture {path}; run "
            "'bash tests/fixtures/json_source_form/make.sh'",
        )

    return path.read_text(encoding="utf-8")


class FixtureClassTests(unittest.TestCase):
    """
    Prove each fixture reports its own class and no other.
    """

    def test_every_fixture_reports_its_class(self) -> None:
        """
        Report the expected class for each generated fixture.
        """

        for name, expected in FIXTURE_CLASSES:
            with self.subTest(fixture=name):
                findings = check_text(fixture_text(name), name)

                if expected is None:
                    self.assertEqual(
                        [item.format() for item in findings],
                        [],
                    )

                    continue

                self.assertTrue(findings, f"{name} reported nothing")
                self.assertTrue(
                    any(expected in item.message for item in findings),
                    f"{name} reported {[i.message for i in findings]}",
                )

    def test_fixture_classes_do_not_overlap(self) -> None:
        """
        Prove no negative fixture reports another fixture's class.
        """

        overflow = check_text(
            fixture_text("inline_overflow.json"),
            "inline_overflow.json",
        )
        fits = check_text(
            fixture_text("expanded_fits.json"),
            "expanded_fits.json",
        )

        self.assertNotIn(
            "fits the 79-column budget inline",
            " ".join(item.message for item in overflow),
        )
        self.assertNotIn(
            "past the 79-column budget",
            " ".join(item.message for item in fits),
        )

    def test_every_finding_carries_the_owner_id(self) -> None:
        """
        Prove no finding escapes with a different stable identifier.
        """

        for name, _ in FIXTURE_CLASSES:
            with self.subTest(fixture=name):
                for finding in check_text(fixture_text(name), name):
                    self.assertEqual(finding.rule_id, "JSON.SOURCE.FORM")


class CanonicalRenderingTests(unittest.TestCase):
    """
    Prove the rendering is total, stable, and value-preserving.
    """

    def test_canonical_fixture_is_a_fixpoint(self) -> None:
        """
        Prove re-rendering the canonical fixture changes nothing.
        """

        text = fixture_text("canonical.json")

        self.assertEqual(canonical(json.loads(text)), text)

    def test_rewriting_preserves_the_parsed_value(self) -> None:
        """
        Prove the rendering never alters a document's value.
        """

        for name, _ in FIXTURE_CLASSES:
            if name in {"unreadable.json", "duplicate_key.json"}:
                continue

            with self.subTest(fixture=name):
                value = json.loads(fixture_text(name))

                self.assertEqual(json.loads(canonical(value)), value)

    def test_rewriting_reaches_zero_findings(self) -> None:
        """
        Prove every repairable fixture is clean after one rewrite.
        """

        for name, _ in FIXTURE_CLASSES:
            if name in {"unreadable.json", "duplicate_key.json"}:
                continue

            with self.subTest(fixture=name):
                rewritten = canonical(json.loads(fixture_text(name)))

                self.assertEqual(check_text(rewritten, name), [])

    def test_budget_decides_structure_not_line_length(self) -> None:
        """
        Prove an indivisible scalar past the budget is left alone.
        """

        scalar = "x" * (BUDGET * 2)
        text = canonical({"note": scalar})

        self.assertEqual(check_text(text, "scalar.json"), [])
        self.assertTrue(
            any(len(line) > BUDGET for line in text.split("\n")),
        )

    def test_a_long_key_does_not_cascade_into_its_members(self) -> None:
        """
        Prove members indent from the line, not from the key's end column.
        """

        key = "a_deliberately_long_key_name_used_to_force_expansion"
        text = canonical({key: ["first value", "second value", "third value"]})
        member = next(
            line for line in text.split("\n") if "first value" in line
        )

        self.assertEqual(len(member) - len(member.lstrip(" ")), 4)

    def test_record_per_line_rows_survive_when_they_fit(self) -> None:
        """
        Prove a tabular array whose rows fit is preserved inline.
        """

        rows = [
            {"callable": "atria", "conceptual_names": ["Atria"]},
            {"callable": "bwa", "conceptual_names": ["BWA"]},
        ]
        text = canonical({"commands": rows})

        self.assertIn(
            '    {"callable": "atria", "conceptual_names": ["Atria"]},',
            text,
        )


class PerturbationTests(unittest.TestCase):
    """
    Prove each assertion still fails when the input stops conforming.
    """

    def test_widening_a_fitting_structure_forces_expansion(self) -> None:
        """
        Prove the budget boundary is live rather than incidental.
        """

        fits = ["x" * 30]
        wide = ["x" * 90]

        self.assertEqual(len(render(fits, 0, 0)), 1)
        self.assertGreater(len(render(wide, 0, 0)), 1)

    def test_one_added_space_is_reported(self) -> None:
        """
        Prove an indentation finding is not a rounding artifact.
        """

        text = canonical({"paths": ["first value", "second value"] * 6})
        lines = text.split("\n")
        lines[1] = " " + lines[1]
        perturbed = "\n".join(lines)

        self.assertTrue(
            any(
                "indentation is" in item.message
                for item in check_text(perturbed, "perturbed.json")
            ),
        )


class DiscoveryTests(unittest.TestCase):
    """
    Prove discovery reaches maintained JSON and nothing else.
    """

    def test_discovery_inspects_a_nonzero_corpus(self) -> None:
        """
        Prove a zero-finding report is not a report over zero paths.
        """

        self.assertGreater(len(maintained_paths(ROOT)), 0)

    def test_generated_fixtures_stay_invisible(self) -> None:
        """
        Prove an ignored negative fixture is never discovered as source.
        """

        discovered = maintained_paths(ROOT)

        self.assertEqual(
            [path for path in discovered if "json_source_form" in path],
            [],
        )

    def test_serialized_inventories_are_excluded(self) -> None:
        """
        Prove producer-owned inventories stay outside the author rule.
        """

        discovered = maintained_paths(ROOT)

        self.assertEqual(
            [path for path in discovered if path.startswith("dev/audit/")],
            [],
        )


if __name__ == "__main__":
    unittest.main()
