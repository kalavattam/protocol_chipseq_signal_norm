#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_help_heredoc_reflow.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Focused regressions for the changed-content help-heredoc prose gate.
"""

from __future__ import annotations

import subprocess
import tempfile
import unittest
from pathlib import Path

from dev.audit.help_heredoc_reflow import (
    RULE_EXAMPLE_COLLAPSED,
    RULE_EXAMPLE_FENCE,
    RULE_EXAMPLE_FINAL_CONTINUATION,
    RULE_EXAMPLE_INDENT,
    RULE_EXAMPLE_RENDERED_STRUCTURE,
    RULE_EXAMPLE_SOURCE_BACKSLASH,
    RULE_EXAMPLE_TAB_INDENT,
    RULE_EXAMPLE_TRAILING_WHITESPACE,
    RULE_ID,
    audit_rendered_example_commands,
    audit_source_example_commands,
    extract_help_heredocs,
    extract_rendered_examples_text,
    fix_source_example_commands,
    scan_repository,
)

REPO_ROOT = Path(__file__).resolve().parents[3]

WRAPPED = """
#!/usr/bin/env bash
function show_command_help() {
    cat << 'EOM'
Notes
-----
Returns 0 when the command is available. Exits through 'die' when the
command name is missing or unavailable.
EOM
}
""".removeprefix("\n")


class HelpHeredocReflowTest(unittest.TestCase):
    """
    Prove strict failures and bounded structured-content exclusions.
    """

    def make_repo(self, root: Path) -> Path:
        subprocess.run(["git", "init", "-q", str(root)], check=True)
        subprocess.run(
            [
                "git",
                "-C",
                str(root),
                "config",
                "user.email",
                "style@example.test",
            ],
            check=True,
        )
        subprocess.run(
            ["git", "-C", str(root), "config", "user.name", "Style Test"],
            check=True,
        )
        (root / "seed.txt").write_text("seed\n", encoding="utf-8")
        subprocess.run(["git", "-C", str(root), "add", "seed.txt"], check=True)
        subprocess.run(
            ["git", "-C", str(root), "commit", "-qm", "seed"],
            check=True,
        )

        return root

    def commit_shell(self, root: Path, path: str, text: str) -> Path:
        target = root / path
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(text, encoding="utf-8")
        subprocess.run(["git", "-C", str(root), "add", path], check=True)
        subprocess.run(
            ["git", "-C", str(root), "commit", "-qm", path],
            check=True,
        )

        return target

    def scan_text(
        self,
        text: str,
        path: str = "new_help.sh",
    ) -> list[object]:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = self.make_repo(Path(temp_dir) / "repo")
            target = root / path
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(text, encoding="utf-8")

            return scan_repository(root)

    def notes_help(self, *lines: str) -> str:
        return "\n".join(
            [
                "#!/usr/bin/env bash",
                "function show_notes_help() {",
                "    cat << 'EOM'",
                "Notes",
                "-----",
                *lines,
                "EOM",
                "}",
                "",
            ],
        )

    def test_rule_identity_is_stable(self) -> None:
        self.assertEqual(RULE_ID, "HELP.HEREDOC.SOURCE_REFLOW")

    def test_rejects_exact_wrapped_command_contract(self) -> None:
        findings = self.scan_text(WRAPPED)

        self.assertEqual(len(findings), 1)
        self.assertEqual(findings[0].physical_lines, (6, 7))
        self.assertEqual(
            findings[0].rendered_prose,
            (
                "Returns 0 when the command is available. Exits through 'die' "
                "when the command name is missing or unavailable."
            ),
        )

    def test_rejects_existing_one_line_paragraph_reflowed_near_boundary(
        self,
    ) -> None:
        one_line = WRAPPED.replace(
            (
                "Returns 0 when the command is available. Exits through 'die' "
                "when the\ncommand name is missing or unavailable."
            ),
            (
                "This existing paragraph has enough ordinary words to sit "
                "near the preferred source width before its final clause "
                "remains attached."
            ),
        )
        reflowed = one_line.replace(
            (
                "This existing paragraph has enough ordinary words to sit "
                "near the preferred source width before its final clause "
                "remains attached."
            ),
            (
                "This existing paragraph has enough ordinary words to sit "
                "near the preferred\nsource width before its final clause "
                "remains attached."
            ),
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            root = self.make_repo(Path(temp_dir) / "repo")
            target = self.commit_shell(root, "tracked.sh", one_line)
            target.write_text(reflowed, encoding="utf-8")

            findings = scan_repository(root)

        self.assertEqual(len(findings), 1)
        self.assertEqual(findings[0].path, "tracked.sh")

    def test_rejects_new_wrapped_prose_in_untracked_shell_file(self) -> None:
        findings = self.scan_text(WRAPPED, "new/untracked_help.sh")

        self.assertEqual(
            [item.path for item in findings],
            ["new/untracked_help.sh"],
        )

    def test_rejects_wrapped_prose_in_changed_externally_owned_help_file(
        self,
    ) -> None:
        clean = WRAPPED.replace(
            (
                "Returns 0 when the command is available. Exits through 'die' "
                "when the\ncommand name is missing or unavailable."
            ),
            "Returns 0 when the command is available.",
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            root = self.make_repo(Path(temp_dir) / "repo")
            target = self.commit_shell(root, "external/owned/help.sh", clean)
            target.write_text(WRAPPED, encoding="utf-8")

            findings = scan_repository(root)

        self.assertEqual(
            [item.path for item in findings],
            ["external/owned/help.sh"],
        )

    def test_rejects_uppercase_proper_name_continuation(self) -> None:
        findings = self.scan_text(
            self.notes_help(
                (
                    "This paragraph explains the configured execution "
                    "environment before"
                ),
                "Samtools performs the requested alignment inspection.",
            ),
        )

        self.assertEqual(len(findings), 1)

    def test_rejects_wrapped_prose_with_short_first_line(self) -> None:
        findings = self.scan_text(
            self.notes_help(
                "This short prose line continues",
                "onto another ordinary physical source line.",
            ),
        )

        self.assertEqual(len(findings), 1)

    def test_rejects_wrapped_prose_with_long_first_line(self) -> None:
        findings = self.scan_text(
            self.notes_help(
                (
                    "This deliberately long ordinary prose line already "
                    "extends far beyond eighty-eight source characters before "
                    "it"
                ),
                "continues mechanically on the next physical line.",
            ),
        )

        self.assertEqual(len(findings), 1)

    def test_rejects_consecutive_sentences_after_terminal_punctuation(
        self,
    ) -> None:
        findings = self.scan_text(
            self.notes_help(
                "The first ordinary sentence ends on this physical line.",
                (
                    "The next ordinary sentence remains part of the same "
                    "prose paragraph."
                ),
            ),
        )

        self.assertEqual(len(findings), 1)

    def test_rejects_every_pair_in_three_line_ordinary_paragraph(self) -> None:
        findings = self.scan_text(
            self.notes_help(
                "This ordinary paragraph starts on one physical line",
                "and continues on a second physical line",
                "before ending on a third physical line.",
            ),
        )

        self.assertEqual(
            [item.physical_lines for item in findings],
            [(6, 7), (7, 8)],
        )

    def test_accepts_prose_paragraphs_separated_by_blank_line(self) -> None:
        text = self.notes_help(
            "This is the first complete ordinary prose paragraph.",
            "",
            "This is the second complete ordinary prose paragraph.",
        )

        self.assertEqual(self.scan_text(text), [])

    def test_accepts_one_line_prose_and_renders_it_on_one_output_line(
        self,
    ) -> None:
        one_line = WRAPPED.replace(
            (
                "Returns 0 when the command is available. Exits through 'die' "
                "when the\ncommand name is missing or unavailable."
            ),
            (
                "Returns 0 when the command is available. Exits through 'die' "
                "when the command name is missing or unavailable."
            ),
        )

        self.assertEqual(self.scan_text(one_line), [])

        with tempfile.TemporaryDirectory() as temp_dir:
            script = Path(temp_dir) / "render.sh"
            script.write_text(
                one_line + "\nshow_command_help\n",
                encoding="utf-8",
            )

            result = subprocess.run(
                ["bash", str(script)],
                text=True,
                capture_output=True,
                check=True,
            )

        self.assertIn(
            (
                "Returns 0 when the command is available. Exits through 'die' "
                "when the command name is missing or unavailable."
            ),
            result.stdout.splitlines(),
        )

    def test_accepts_recognized_structured_help_content(self) -> None:
        structured = """
#!/usr/bin/env bash
function help_structured() {
    cat << 'EOM'
Usage
-----
  tool_name
    [--input <file>] [--output <file>]

Parameters
----------
  -i, --input : file
    Input path with a structurally distinct indented description.

  -o, --output : file
    Output path with a structurally distinct indented description.

Notes
-----
  - First bullet with enough content to extend near the ordinary source-width boundary
    while retaining deliberate list indentation.
  1. First numbered item with enough content to extend near the ordinary source-width
     while retaining deliberate list indentation.

  | Name | Meaning      |
  | :--- | :---         |
  | foo  | first value  |
  | bar  | second value |

Examples
--------
  '''bash
  command_name --input input.txt
  printf '%s\\n' "intentionally structured literal output"
  '''
EOM
}
"""

        self.assertEqual(self.scan_text(structured), [])

    def test_accepts_long_one_line_prose_beyond_eighty_characters(
        self,
    ) -> None:
        long_line = WRAPPED.replace(
            (
                "Returns 0 when the command is available. Exits through 'die' "
                "when the\ncommand name is missing or unavailable."
            ),
            (
                "This ordinary help-heredoc paragraph deliberately remains on "
                "one physical source line even though its complete rendered "
                "wording is much longer than eighty characters."
            ),
        )

        self.assertEqual(self.scan_text(long_line), [])

    def test_shell_line_length_gate_excludes_long_help_prose(self) -> None:
        long_line = WRAPPED.replace(
            (
                "Returns 0 when the command is available. Exits through 'die' "
                "when the\ncommand name is missing or unavailable."
            ),
            (
                "This ordinary help-heredoc paragraph deliberately remains on "
                "one physical source line even though its complete rendered "
                "wording is much longer than eighty characters."
            ),
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            script = Path(temp_dir) / "long_help.sh"
            script.write_text(long_line, encoding="utf-8")

            result = subprocess.run(
                [
                    "bash",
                    str(
                        REPO_ROOT
                        / (
                            "tests/contract/repository/test_shell_line_length."
                            "sh"
                        ),
                    ),
                    str(script),
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
            )

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertNotIn("long_help.sh:", result.stdout + result.stderr)


class HelpExampleCommandTest(unittest.TestCase):
    """
    Prove source and rendered multiline-command invariants.
    """

    def source_help(
        self,
        *code: str,
        opener: str = "cat << EOM",
        owner: str = "help_demo",
        entries: int = 1,
    ) -> str:
        """
        Build one source help heredoc with literal established fences.
        """

        examples: list[str] = []

        for number in range(1, entries + 1):
            examples.extend(
                [
                    f"  {number}. Demonstrate invocation {number}.",
                    "    '''bash",
                    *(f"    {line}" for line in code),
                    "    '''",
                ],
            )

            if number != entries:
                examples.append("")

        return "\n".join(
            [
                "#!/usr/bin/env bash",
                f"function {owner}() {{",
                f"    {opener}",
                "Examples",
                "--------",
                *examples,
                "EOM",
                "}",
                "",
            ],
        )

    def rendered_help(self, *code: str, entries: int = 1) -> str:
        """
        Build one rendered Examples document.
        """

        examples: list[str] = []

        for number in range(1, entries + 1):
            examples.extend(
                [
                    f"  {number}. Demonstrate invocation {number}.",
                    "    '''bash",
                    *(f"    {line}" for line in code),
                    "    '''",
                ],
            )

            if number != entries:
                examples.append("")

        return "\n".join(["Examples", "--------", *examples, ""])

    def source_rules(self, text: str) -> set[str]:
        return {
            finding.rule_id
            for finding in audit_source_example_commands(text).findings
        }

    def rendered_rules(self, text: str) -> set[str]:
        return {
            finding.rule_id
            for finding in audit_rendered_example_commands(text).findings
        }

    def test_extracts_supported_quoted_unquoted_and_strip_tab_openers(
        self,
    ) -> None:
        forms = (
            ("cat <<EOM", False, False),
            ("cat << EOM", False, False),
            ("cat <<'EOM'", True, False),
            ("cat << 'EOM'", True, False),
            ("cat <<-EOM", False, True),
            ("cat <<-'EOM'", True, True),
        )

        for opener, quoted, strips_tabs in forms:
            with self.subTest(opener=opener):
                heredocs = extract_help_heredocs(
                    self.source_help("bash demo.sh", opener=opener),
                )

                self.assertEqual(len(heredocs), 1)
                self.assertEqual(heredocs[0].quoted, quoted)
                self.assertEqual(heredocs[0].strips_tabs, strips_tabs)

    def test_source_valid_unquoted_double_backslashes_pass(self) -> None:
        text = self.source_help(
            r"bash demo.sh \\",
            r"    --mode ${mode} \\",
            "    --output result.txt",
            entries=2,
        )
        audit = audit_source_example_commands(text)

        self.assertEqual(audit.findings, [])
        self.assertEqual(
            [row.classification for row in audit.inventory],
            ["checked_simple_command", "checked_simple_command"],
        )

    def test_source_valid_quoted_single_backslashes_pass(self) -> None:
        text = self.source_help(
            "bash demo.sh \\",
            "    --mode one \\",
            "    --output result.txt",
            opener="cat << 'EOM'",
        )

        self.assertEqual(audit_source_example_commands(text).findings, [])

    def test_source_deep_relative_indentation_and_single_line_pass(
        self,
    ) -> None:
        deep = self.source_help(
            r"    bash demo.sh \\",
            r"        --mode one \\",
            "        --output result.txt",
        )
        single = self.source_help("bash demo.sh --mode one")

        self.assertEqual(audit_source_example_commands(deep).findings, [])
        self.assertEqual(audit_source_example_commands(single).findings, [])

    def test_source_multiple_independent_single_line_commands_pass(
        self,
    ) -> None:
        text = self.source_help(
            "bash demo.sh --mode one",
            "printf '%s\\n' complete",
        )
        audit = audit_source_example_commands(text)

        self.assertEqual(audit.findings, [])
        self.assertEqual(
            [row.classification for row in audit.inventory],
            ["checked_simple_command", "checked_simple_command"],
        )

    def test_source_backslash_counts_fail_for_unquoted_nonfinal_lines(
        self,
    ) -> None:
        fixtures = {
            "one": "bash demo.sh \\",
            "zero": "bash demo.sh",
            "three": "bash demo.sh " + "\\" * 3,
        }

        for label, first in fixtures.items():
            with self.subTest(label=label):
                rules = self.source_rules(
                    self.source_help(first, "    --mode one"),
                )

                self.assertIn(RULE_EXAMPLE_SOURCE_BACKSLASH, rules)

    def test_source_relative_indentation_matrix_fails(self) -> None:
        for spaces in (0, 1, 2, 3, 5):
            with self.subTest(spaces=spaces):
                rules = self.source_rules(
                    self.source_help(
                        r"bash demo.sh \\",
                        f"{' ' * spaces}--mode one",
                    ),
                )

                self.assertIn(RULE_EXAMPLE_INDENT, rules)

    def test_source_inconsistent_indentation_fails(self) -> None:
        rules = self.source_rules(
            self.source_help(
                r"bash demo.sh \\",
                r"    --mode one \\",
                "     --output result.txt",
            ),
        )

        self.assertIn(RULE_EXAMPLE_INDENT, rules)

    def test_source_tab_indentation_fails(self) -> None:
        rules = self.source_rules(
            self.source_help(r"bash demo.sh \\", "\t--mode one"),
        )

        self.assertIn(RULE_EXAMPLE_TAB_INDENT, rules)

    def test_source_trailing_spaces_and_tabs_after_backslashes_fail(
        self,
    ) -> None:
        for suffix in ("  ", "\t"):
            with self.subTest(suffix=repr(suffix)):
                rules = self.source_rules(
                    self.source_help(
                        r"bash demo.sh \\" + suffix,
                        "    --mode one",
                    ),
                )

                self.assertIn(RULE_EXAMPLE_TRAILING_WHITESPACE, rules)

    def test_source_final_line_continuation_fails(self) -> None:
        rules = self.source_rules(
            self.source_help(
                r"bash demo.sh \\",
                r"    --mode one \\",
            ),
        )

        self.assertIn(RULE_EXAMPLE_FINAL_CONTINUATION, rules)

    def test_source_collapsed_layout_fails(self) -> None:
        rules = self.source_rules(
            self.source_help(r"bash demo.sh \\      --mode one"),
        )

        self.assertIn(RULE_EXAMPLE_COLLAPSED, rules)

    def test_source_malformed_or_unterminated_fences_fail(self) -> None:
        valid = self.source_help("bash demo.sh")
        fixtures = (
            valid.replace("    '''bash", "    ```bash"),
            valid.replace("    '''\nEOM", "EOM"),
            valid.replace("    '''bash", "    \\```bash"),
        )

        for text in fixtures:
            with self.subTest(text=text):
                self.assertIn(RULE_EXAMPLE_FENCE, self.source_rules(text))

    def test_explicit_fixer_is_conservative_and_idempotent(self) -> None:
        text = self.source_help(
            "bash demo.sh",
            "  --mode ${mode}",
            "  --output result.txt \\",
        )
        fixed, count, refused = fix_source_example_commands(text)

        self.assertEqual(count, 1)
        self.assertEqual(refused, ())
        self.assertIn("cat << EOM", fixed)
        self.assertIn("'''bash", fixed)
        self.assertIn("--mode ${mode}", fixed)
        self.assertEqual(audit_source_example_commands(fixed).findings, [])

        second, second_count, second_refused = fix_source_example_commands(
            fixed,
        )

        self.assertEqual(second, fixed)
        self.assertEqual(second_count, 0)
        self.assertEqual(second_refused, ())

    def test_fixer_preserves_complex_and_refuses_collapsed_or_malformed(
        self,
    ) -> None:
        complex_text = self.source_help(
            "if demo; then",
            "    printf '%s\\n' ok",
            "fi",
        )
        fixed, count, refused = fix_source_example_commands(complex_text)

        self.assertEqual(fixed, complex_text)
        self.assertEqual(count, 0)
        self.assertEqual(refused, ())

        for text in (
            self.source_help(r"bash demo.sh \\      --mode one"),
            self.source_help("bash demo.sh").replace(
                "    '''bash",
                "    ```bash",
            ),
        ):
            with self.subTest(text=text):
                fixed, count, refused = fix_source_example_commands(text)

                self.assertEqual(fixed, text)
                self.assertEqual(count, 0)
                self.assertTrue(refused)

    def test_complex_snippets_are_explicitly_excluded(self) -> None:
        snippets = (
            ("if demo; then", "    printf '%s\\n' ok", "fi"),
            ("for item in one two; do", '    demo "${item}"', "done"),
            ("values=(one two)", "printf '%s\\n' \"${values[@]}\""),
            ("demo one |", "    sed 's/one/two/'"),
            ("(demo one", "    demo two)"),
        )

        for code in snippets:
            with self.subTest(code=code):
                audit = audit_source_example_commands(self.source_help(*code))

                self.assertEqual(audit.findings, [])
                self.assertEqual(
                    [row.classification for row in audit.inventory],
                    ["deliberately_excluded_complex_snippet"],
                )
                self.assertTrue(audit.inventory[0].reason)

    def test_rendered_valid_multiline_structure_passes(self) -> None:
        text = self.rendered_help(
            "bash demo.sh \\",
            "    --mode one \\",
            "    --output result.txt",
            entries=2,
        )

        self.assertEqual(audit_rendered_example_commands(text).findings, [])

    def test_rendered_examples_extractor_preserves_exact_physical_text(
        self,
    ) -> None:
        examples = self.rendered_help(
            "bash demo.sh \\",
            "    --mode one \\",
            "    --output result.txt",
            entries=2,
        )
        text = "Usage\n-----\n  demo.sh\n\n" + examples

        self.assertEqual(extract_rendered_examples_text(text), examples)

        with self.assertRaises(ValueError):
            extract_rendered_examples_text("Usage\n-----\n  demo.sh\n")

    def test_rendered_backslash_and_physical_line_failures(self) -> None:
        fixtures = (
            self.rendered_help(
                "bash demo.sh \\      --mode one \\",
                "    --output result.txt",
            ),
            self.rendered_help("bash demo.sh", "    --mode one"),
            self.rendered_help(
                r"bash demo.sh \\",
                "    --mode one",
            ),
            self.rendered_help("bash demo.sh \\  ", "    --mode one"),
            self.rendered_help("bash demo.sh \\", "  --mode one"),
            self.rendered_help("bash demo.sh \\", "    --mode one \\"),
        )

        for text in fixtures:
            with self.subTest(text=text):
                self.assertTrue(self.rendered_rules(text))

    def test_rendered_collapsed_and_missing_fence_are_distinguished(
        self,
    ) -> None:
        collapsed = self.rendered_help(
            "bash demo.sh \\      --mode one --output result.txt",
        )
        missing_fence = self.rendered_help("bash demo.sh").replace(
            "    '''bash",
            "    ```bash",
        )

        self.assertIn(RULE_EXAMPLE_COLLAPSED, self.rendered_rules(collapsed))
        self.assertIn(RULE_EXAMPLE_FENCE, self.rendered_rules(missing_fence))
        self.assertIn(
            RULE_EXAMPLE_RENDERED_STRUCTURE,
            self.rendered_rules(collapsed),
        )

    def test_numbered_entries_require_blank_line_separation(self) -> None:
        source = self.source_help("bash demo.sh", entries=2).replace(
            "    '''\n\n  2.",
            "    '''\n  2.",
        )
        rendered = self.rendered_help("bash demo.sh", entries=2).replace(
            "    '''\n\n  2.",
            "    '''\n  2.",
        )

        self.assertIn(RULE_EXAMPLE_FENCE, self.source_rules(source))
        self.assertIn(
            RULE_EXAMPLE_RENDERED_STRUCTURE,
            self.rendered_rules(rendered),
        )


if __name__ == "__main__":
    unittest.main()
