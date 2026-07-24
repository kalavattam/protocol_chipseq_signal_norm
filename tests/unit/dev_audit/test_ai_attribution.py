#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_ai_attribution.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


"""Focused regressions for bounded, model-generic source attribution."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

AUDIT_DIR = Path(__file__).resolve().parents[3] / "dev" / "audit"
sys.path.insert(0, str(AUDIT_DIR))

from ai_attribution import (
    attribution_block,
    check_attribution_source,
    load_applicability_manifest,
    model_tokens,
    normalize_attribution_source,
    source_header_inventory,
)


def source_file(
    attribution: str | None,
    *,
    basename: str = "tool.sh",
    shebang: str = "#!/usr/bin/env bash",
    copyright_line: str = "# Copyright 2024-2026 by Kris Alavattam",
) -> str:
    """Wrap an optional attribution block in a complete header fixture."""

    rows = [
        shebang,
        "# -*- coding: utf-8 -*-",
        "#",
        f"# Script: {basename}",
        "#",
        copyright_line,
        "# Email: kalavattam@gmail.com",
        "#",
    ]
    if attribution is not None:
        rows.extend((attribution, "#"))
    rows.extend(("# Distributed under the MIT license.", "", "", ""))
    return "\n".join(rows)


def flattened_attribution(text: str) -> str:
    """Return one source's attribution as a single logical line."""

    block = attribution_block(text)
    assert block is not None
    return " ".join(
        line.removeprefix("#").strip() for line in block[2].splitlines()
    )


class AiAttributionTest(unittest.TestCase):
    """Prove applicability, exact model tokens, tools, and domains."""

    def test_recognized_generic_and_explicit_model_forms_pass(self) -> None:
        headers = (
            "# OpenAI Codex (GPT-5-series models) was used in development.",
            "# OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:\n# GPT-5.6) were used in development and documentation.",
            "# OpenAI ChatGPT and Codex (GPT-5-series models; most recent: GPT-5.6) were\n# used in development and documentation.",
            "# OpenAI ChatGPT (GPT-5-series models; most recent: GPT-5.5, GPT-5.6)\n# was used in documentation.",
            "# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in development and\n# documentation.",
            "# OpenAI ChatGPT and Codex (GPT-7.2-preview) were used in development.",
        )
        for header in headers:
            with self.subTest(header=header):
                self.assertEqual(
                    check_attribution_source(source_file(header), "tool.sh"),
                    [],
                )

    def test_model_declaration_rejects_residue_and_bad_separators(self) -> None:
        invalid = (
            "banana, GPT-5.6, ???",
            "GPT-5.6!",
            "GPT-5.5, , GPT-5.6",
            ", GPT-5.6",
            "GPT-5.6,",
            "GPT-5.5,GPT-5.6",
            "GPT-5.5; GPT-5.6",
            "GPT-5-series models, GPT-5.6",
            "GPT-5-series models;most recent: GPT-5.6",
            "GPT-5-series models; most recent:",
            "GPT-5-series models; most recent: banana",
            "GPT-5-series models; most recent: GPT-5.6; GPT-5.7",
            "GPT-4- and GPT-5-series models; most recent: GPT-5.6,",
        )
        for models in invalid:
            header = (
                f"# OpenAI Codex ({models}) was used in development."
            )
            with self.subTest(models=models):
                findings = check_attribution_source(
                    source_file(header),
                    "tool.sh",
                )
                self.assertIn(
                    "SOURCE.HEADER.AI_ATTRIBUTION",
                    {item.rule_id for item in findings},
                )

    def test_standing_validation_accepts_coherent_no_ai_profile(self) -> None:
        self.assertEqual(
            check_attribution_source(source_file(None), "tool.sh"),
            [],
        )

    def test_explicit_applicability_requires_attribution(self) -> None:
        findings = check_attribution_source(
            source_file(None),
            "tool.sh",
            required_models=("GPT-5.6",),
            required_contribution_domain="development",
            required_tools="codex",
        )
        self.assertEqual(
            {item.rule_id for item in findings},
            {"SOURCE.HEADER.AI_ATTRIBUTION"},
        )

    def test_malformed_attribution_like_header_fails_standing_scan(self) -> None:
        text = source_file(None).replace(
            "# Distributed under the MIT license.",
            "# OpenAI tool GPT-series helped.\n#\n"
            "# Distributed under the MIT license.",
        )
        self.assertIn(
            "SOURCE.HEADER.AI_ATTRIBUTION",
            {item.rule_id for item in check_attribution_source(text, "tool.sh")},
        )

    def test_single_and_combined_tool_forms_are_truthful_and_accepted(self) -> None:
        headers = (
            "# OpenAI ChatGPT (GPT-5.6) was used in documentation.",
            "# OpenAI Codex (GPT-5.6) was used in development.",
            "# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and\n# documentation.",
        )
        for header in headers:
            with self.subTest(header=header):
                self.assertEqual(
                    check_attribution_source(source_file(header), "tool.sh"),
                    [],
                )

    def test_wrong_tool_verb_form_is_rejected(self) -> None:
        findings = check_attribution_source(
            source_file(
                "# OpenAI Codex (GPT-5.6) were used in development."
            ),
            "tool.sh",
        )
        self.assertIn(
            "SOURCE.HEADER.AI_ATTRIBUTION",
            {item.rule_id for item in findings},
        )

    def test_required_domain_and_tools_must_match_observed_wording(self) -> None:
        text = source_file(
            "# OpenAI Codex (GPT-5.6) was used in development."
        )
        self.assertEqual(
            check_attribution_source(
                text,
                "tool.sh",
                required_models=("GPT-5.6",),
                required_contribution_domain="development",
                required_tools="codex",
            ),
            [],
        )
        messages = {
            item.message
            for item in check_attribution_source(
                text,
                "tool.sh",
                required_models=("GPT-5.6",),
                required_contribution_domain="documentation",
                required_tools="chatgpt",
            )
        }
        self.assertIn(
            "observed contribution domain does not match explicit applicability",
            messages,
        )
        self.assertIn(
            "observed OpenAI tool set does not match explicit applicability",
            messages,
        )

    def test_exact_model_counters_do_not_confuse_gpt5_and_gpt56(self) -> None:
        text = source_file(
            "# OpenAI Codex (GPT-5.6) was used in development."
        )
        findings = check_attribution_source(
            text,
            "tool.sh",
            required_models=("GPT-5", "GPT-5.6"),
        )
        self.assertEqual(model_tokens("GPT-5, GPT-5.6"), ("GPT-5", "GPT-5.6"))
        self.assertTrue(
            any("GPT-5'" in item.message for item in findings)
        )
        self.assertFalse(
            any("GPT-5.6' must appear" in item.message for item in findings)
        )

    def test_duplicate_exact_model_token_fails_unique_owner(self) -> None:
        findings = check_attribution_source(
            source_file(
                "# OpenAI ChatGPT and Codex (GPT-5.6, GPT-5.6) were used in\n# development and documentation."
            ),
            "tool.sh",
        )
        self.assertIn(
            "SOURCE.HEADER.AI_ATTRIBUTION.UNIQUE",
            {item.rule_id for item in findings},
        )

    def test_generic_series_tokens_are_parsed_without_substrings(self) -> None:
        self.assertEqual(
            model_tokens(
                "GPT-4- and GPT-5-series models; most recent: GPT-5.6"
            ),
            ("GPT-4-series", "GPT-5-series", "GPT-5.6"),
        )
        findings = check_attribution_source(
            source_file(
                "# OpenAI ChatGPT (GPT-series) was used in documentation."
            ),
            "tool.sh",
        )
        self.assertTrue(
            any(item.rule_id == "SOURCE.HEADER.AI_ATTRIBUTION" for item in findings)
        )

    def test_repeatable_and_future_required_models_are_exact_data(self) -> None:
        text = source_file(
            "# OpenAI Codex (GPT-6.1, GPT-7.2-preview) was used in development."
        )
        self.assertEqual(
            check_attribution_source(
                text,
                "tool.sh",
                required_models=("GPT-6.1", "GPT-7.2-preview"),
            ),
            [],
        )

    def test_normalization_is_canonical_and_idempotent(self) -> None:
        cases = (
            (
                "# OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models) were used in\n# development and documentation.",
                "GPT-4- and GPT-5-series models; most recent: GPT-5.6",
            ),
            (
                "# OpenAI Codex (GPT-5-series models; most recent: GPT-5.5) was\n# used in development.",
                "GPT-5-series models; most recent: GPT-5.5, GPT-5.6",
            ),
            (
                "# OpenAI Codex (GPT-5.5) was used in development.",
                "GPT-5.5, GPT-5.6",
            ),
        )
        for header, expected_models in cases:
            with self.subTest(header=header):
                normalized = normalize_attribution_source(
                    source_file(header),
                    required_models=("GPT-5.6",),
                )
                self.assertIn(
                    f"({expected_models})",
                    flattened_attribution(normalized),
                )
                self.assertEqual(
                    normalize_attribution_source(
                        normalized,
                        required_models=("GPT-5.6",),
                    ),
                    normalized,
                )

    def test_normalization_does_not_invent_a_missing_series(self) -> None:
        text = source_file(
            "# OpenAI Codex (GPT-5.6) was used in development."
        )
        with self.assertRaisesRegex(ValueError, "generic-series requirement"):
            normalize_attribution_source(
                text,
                required_models=("GPT-6-series",),
            )

    def test_missing_attribution_requires_domain_models_and_tools(self) -> None:
        text = source_file(None)
        with self.assertRaisesRegex(ValueError, "explicit contribution domain"):
            normalize_attribution_source(text, required_models=("GPT-5.6",))
        with self.assertRaisesRegex(ValueError, "explicit OpenAI tool set"):
            normalize_attribution_source(
                text,
                contribution_domain="documentation",
                required_models=("GPT-5.6",),
            )
        normalized = normalize_attribution_source(
            text,
            contribution_domain="documentation",
            attribution_tools="chatgpt",
            required_models=("GPT-5.6",),
        )
        self.assertIn(
            "OpenAI ChatGPT (GPT-5.6) was used in documentation.",
            normalized,
        )
        self.assertNotIn("Codex", normalized)

    def test_body_attribution_does_not_change_no_ai_header_profile(self) -> None:
        text = source_file(None) + (
            "printf 'body\\n'\n"
            "# OpenAI Codex (GPT-5.6) was used in documentation.\n"
        )
        self.assertEqual(check_attribution_source(text, "tool.sh"), [])

    def test_header_structure_basename_width_year_and_profiles(self) -> None:
        valid = source_file(
            "# OpenAI Codex (GPT-5.6) was used in development."
        )
        cases = (
            valid.replace("# Script: tool.sh\n#\n", ""),
            valid.replace("# Script: tool.sh", "# Script: other.sh"),
            valid.replace(
                "# Copyright 2024-2026 by Kris Alavattam\n"
                "# Email: kalavattam@gmail.com\n",
                "# Email: kalavattam@gmail.com\n"
                "# Copyright 2024-2026 by Kris Alavattam\n",
            ),
            valid.replace("# Script: tool.sh\n#\n", "# Script: tool.sh\n# \n"),
            source_file(
                "# OpenAI ChatGPT and Codex (GPT-5.4, GPT-5.5, GPT-5.6) were used in development and documentation."
            ),
            source_file(
                "# OpenAI Codex (GPT-5.6) was used in development.",
                copyright_line="# Copyright 2024-2025 by Kris Alavattam",
            ),
        )
        for text in cases:
            with self.subTest(text=text):
                self.assertTrue(check_attribution_source(text, "tool.sh"))

    def test_boundary_profiles_and_body_spacing(self) -> None:
        python = source_file(
            "# OpenAI Codex (GPT-5.6) was used in development.",
            basename="tool.py",
            shebang="#!/usr/bin/env python3",
        )
        posix = source_file(
            None,
            basename="install_envs_entrypoint.sh",
            shebang="#!/bin/sh",
            copyright_line="# Copyright 2026 by Kris Alavattam",
        )
        self.assertEqual(check_attribution_source(python, "tool.py"), [])
        self.assertEqual(
            check_attribution_source(posix, "install_envs_entrypoint.sh"),
            [],
        )
        wrong = python.replace("#!/usr/bin/env python3", "#!/bin/python")
        self.assertIn(
            "SOURCE.HEADER.STRUCTURE",
            {item.rule_id for item in check_attribution_source(wrong, "tool.py")},
        )
        one_blank = python.rstrip() + "\n\nprint('body')\n"
        self.assertIn(
            "SOURCE.HEADER.STRUCTURE",
            {item.rule_id for item in check_attribution_source(one_blank, "tool.py")},
        )

    def test_duplicate_attribution_blocks_use_generic_unique_id(self) -> None:
        header = "# OpenAI Codex (GPT-5.6) was used in development."
        text = source_file(header).replace(
            f"{header}\n#\n",
            f"{header}\n{header}\n#\n",
        )
        self.assertIn(
            "SOURCE.HEADER.AI_ATTRIBUTION.UNIQUE",
            {item.rule_id for item in check_attribution_source(text, "tool.sh")},
        )

    def test_manifest_and_inventory_compare_observed_requirements(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = source_file(
                "# OpenAI Codex (GPT-5.6) was used in development."
            )
            (root / "tool.sh").write_text(source, encoding="utf-8")
            manifest = root / "cohort.json"
            manifest.write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "assessment_year": 2026,
                        "sources": [
                            {
                                "path": "tool.sh",
                                "required_models": ["GPT-5.6"],
                                "contribution_domain": "development",
                                "tools": "codex",
                            }
                        ],
                    }
                ),
                encoding="utf-8",
            )
            year, requirements = load_applicability_manifest(root, manifest)
            rows = source_header_inventory(
                root,
                ["tool.sh"],
                current_year=year,
                requirements=requirements,
            )
        self.assertEqual(requirements["tool.sh"].required_models, ("GPT-5.6",))
        self.assertEqual(rows[0]["observed_contribution_domain"], "development")
        self.assertEqual(rows[0]["required_contribution_domain"], "development")
        self.assertTrue(rows[0]["contribution_domain_agrees"])
        self.assertTrue(rows[0]["required_models_agree"])
        self.assertTrue(rows[0]["tools_agree"])

    def test_no_ai_inventory_reports_null_observed_attribution(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "tool.sh").write_text(
                source_file(None),
                encoding="utf-8",
            )
            rows = source_header_inventory(root, ["tool.sh"])
        self.assertEqual(len(rows), 1)
        self.assertIsNone(rows[0]["attribution_style"])
        self.assertEqual(rows[0]["observed_models"], [])
        self.assertIsNone(rows[0]["observed_contribution_domain"])
        self.assertIsNone(rows[0]["observed_tools"])


if __name__ == "__main__":
    unittest.main()
