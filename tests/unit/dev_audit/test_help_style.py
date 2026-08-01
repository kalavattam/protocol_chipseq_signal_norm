#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_help_style.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


"""
Focused regressions for structural shell-help conventions.
"""

from __future__ import annotations

import unittest
from pathlib import Path

from dev.audit import help_style
from dev.audit.help_style import (
    check_help_source,
    check_python_docstrings,
    check_rendered_help,
    scan_repository,
)

REPO_ROOT = Path(__file__).resolve().parents[3]
REGISTRY = REPO_ROOT / "dev/config/command_names.json"


def source(body: str, parser: str = "", owner: str = "show_help") -> str:
    """
    Wrap one rendered help body in a narrowly recognized heredoc.
    """

    return "\n".join(
        [
            "#!/usr/bin/env bash",
            parser,
            f"function {owner}() {{",
            "    cat << 'EOM'",
            body,
            "EOM",
            "}",
            "",
        ],
    )


def rule_ids(body: str, parser: str = "") -> set[str]:
    """
    Return rule identifiers emitted for one fixture.
    """

    return {
        finding.rule_id
        for finding in check_help_source(
            source(body, parser),
            registry_path=REGISTRY,
        )
    }


def findings(
    body: str,
    parser: str = "",
    owner: str = "show_help",
) -> list[help_style.Finding]:
    """
    Return full findings for one focused fixture.
    """

    return check_help_source(
        source(body, parser, owner),
        registry_path=REGISTRY,
    )


def function_body(*sections: str) -> str:
    """
    Build a minimally complete function-help fixture.
    """

    return "\n\n".join(sections)


class HelpStyleTest(unittest.TestCase):
    """
    Prove strict structure while retaining bounded exclusions.
    """

    def test_canonical_documented_help_spelling_passes(self) -> None:
        body = """
Usage
-----
  command_name
    [--help]

Parameters
----------
  -h, --help : flag
    Display this help message and exit.
""".removesuffix("\n")

        self.assertEqual(rule_ids(body), set())

    def test_parameter_row_with_public_short_and_long_aliases_passes(
        self,
    ) -> None:
        body = """
Parameters
----------
  -dr, --dry, --dry_run : flag
    Print the resolved command without running it.
""".removesuffix("\n")

        self.assertEqual(rule_ids(body), set())

    def test_duplicate_parameter_alias_fails(self) -> None:
        body = """
Parameters
----------
  -dr, --dry, -dr, --dry_run : flag
    Print the resolved command without running it.
""".removesuffix("\n")

        self.assertNotIn("HELP.PARAMETER.ALIAS_DUPLICATE", rule_ids(body))

    def test_positional_parameter_row_is_unaffected_by_alias_rules(
        self,
    ) -> None:
        body = """
Parameters
----------
  1  fil_in : file
    Input file.
""".removesuffix("\n")

        self.assertEqual(rule_ids(body), set())

    def test_hidden_help_aliases_in_public_usage_fail(self) -> None:
        body = """
Usage
-----
  command_name
    [-h|--hlp|--help]
""".removesuffix("\n")

        self.assertIn("HELP.DOC.CANONICAL_HELP", rule_ids(body))

    def test_hidden_parser_aliases_do_not_fail(self) -> None:
        body = """
Usage
-----
  command_name
    [--help]
""".removesuffix("\n")
        parser = (
            'function parse_args() { case "${1:-}" in -h|--hlp|--help) return '
            "0 ;; esac; }"
        )

        self.assertEqual(rule_ids(body, parser), set())

    def test_usage_name_alone_with_one_continuation_passes(self) -> None:
        body = """
Usage
-----
  require_optarg
    [--help] opt val [func]
""".removesuffix("\n")

        self.assertEqual(rule_ids(body), set())

    def test_usage_name_and_arguments_on_one_line_fail(self) -> None:
        body = """
Usage
-----
  require_optarg [--help] opt val [func]
""".removesuffix("\n")

        self.assertIn("HELP.USAGE.STRUCTURE", rule_ids(body))

    def test_usage_semantically_grouped_continuations_pass(self) -> None:
        body = """
Usage
-----
  filter_alignment_sc
    [--help] [--threads <int>]
    --fil_in <file> --fil_out <file> [--ref_fa <file>]
    [--mito] [--chk_chr]
""".removesuffix("\n")

        self.assertEqual(rule_ids(body), set())

    def test_usage_misindented_continuation_fails(self) -> None:
        body = """
Usage
-----
  require_optarg
   [--help] opt val [func]
""".removesuffix("\n")

        self.assertIn("HELP.USAGE.STRUCTURE", rule_ids(body))

    def test_parameter_entry_indentation_and_blank_line_pass(self) -> None:
        body = """
Parameters
----------
  1  var_nam : str
    Variable name to use in error messages.

  2  var_val : str
    Value to check.
""".removesuffix("\n")

        self.assertEqual(rule_ids(body), set())

    def test_adjacent_parameter_entries_without_blank_fail(self) -> None:
        body = """
Parameters
----------
  1  var_nam : str
    Variable name to use in error messages.
  2  var_val : str
    Value to check.
""".removesuffix("\n")

        self.assertIn("HELP.ENTRY.BLANK", rule_ids(body))

    def test_under_and_over_indented_parameter_content_fails(self) -> None:
        under = """
Parameters
----------
  1  var_nam : str
   Variable name to use in error messages.
""".removesuffix("\n")
        over = """
Parameters
----------
    1  var_nam : str
      Variable name to use in error messages.
""".removesuffix("\n")

        self.assertIn("HELP.ENTRY.INDENT", rule_ids(under))
        self.assertIn("HELP.ENTRY.INDENT", rule_ids(over))

    def test_inline_parameter_description_fails(self) -> None:
        body = """
Parameters
----------
  --help : flag    Display this help message and exit.
""".removesuffix("\n")

        self.assertIn("HELP.ENTRY.INDENT", rule_ids(body))

    def test_relative_runtime_direct_bullets_pass(self) -> None:
        for body in (
            (
                "  Runtime requirements:\n    - bash >= 4.4\n    - samtools "
                ">= 1.21"
            ),
        ):
            with self.subTest(body=body):
                self.assertEqual(rule_ids("Notes\n-----\n" + body), set())

    def test_runtime_requirements_outside_notes_fail(self) -> None:
        for body in (
            """
  Runtime requirements:
    bash >= 4.4

Usage
-----
  show_help
    [--help]
""".removesuffix("\n"),
            """
Returns
-------
  Returns 0.

  Runtime requirements:
    bash >= 4.4
""".removesuffix("\n"),
        ):
            with self.subTest(body=body):
                self.assertIn("HELP.RUNTIME.PARENT", rule_ids(body))

    def test_runtime_parent_tracks_later_section_boundary(self) -> None:
        body = """
Notes
-----
  Runtime requirements:
    bash >= 4.4

Returns
-------
  Returns 0.
""".removesuffix("\n")

        self.assertEqual(rule_ids(body), set())

    def test_runtime_parent_check_is_idempotent(self) -> None:
        body = """
Returns
-------
  Returns 0.

  Runtime requirements:
    bash >= 4.4
""".removesuffix("\n")

        self.assertEqual(findings(body), findings(body))

    def test_same_and_over_indented_runtime_bullets_fail(self) -> None:
        same = "Notes\n-----\n  Runtime requirements:\n  - bash >= 4.4"
        over = "Notes\n-----\n  Runtime requirements:\n      - bash >= 4.4"

        self.assertIn("HELP.RUNTIME.INDENT", rule_ids(same))
        self.assertIn("HELP.RUNTIME.INDENT", rule_ids(over))

    def test_runtime_block_ends_before_following_section(self) -> None:
        body = """
Notes
-----
  Runtime requirements:
    bash >= 4.4

Returns
-------
  - Prints help.
""".removesuffix("\n")

        self.assertEqual(rule_ids(body), set())

    def test_flat_program_only_runtime_list_passes(self) -> None:
        body = """
Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - grep
    - samtools >= 1.21
""".removesuffix("\n")

        self.assertEqual(rule_ids(body), set())

    def test_redundant_sole_programs_category_fails(self) -> None:
        body = """
Notes
-----
  Runtime requirements:
    - Programs
      + bash >= 4.4
      + samtools >= 1.21
""".removesuffix("\n")

        self.assertIn("HELP.RUNTIME.REDUNDANT_PROGRAMS", rule_ids(body))

    def test_multiple_genuine_runtime_categories_pass(self) -> None:
        body = """
Notes
-----
  Runtime requirements:
    - Programs
      + bash >= 4.4
      + samtools >= 1.21
    - Environment
      + Conda environment 'env_protocol'
    - Files
      + Reference FASTA and index
""".removesuffix("\n")

        self.assertEqual(rule_ids(body), set())

    def test_registered_executable_spellings_pass(self) -> None:
        body = """
Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - samtools
    - bamCompare
""".removesuffix("\n")

        self.assertEqual(rule_ids(body), set())

    def test_slurm_allocation_is_a_resource_not_a_callable(self) -> None:
        body = """
Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - Slurm allocation (when run as a Slurm array task)
""".removesuffix("\n")

        self.assertEqual(rule_ids(body), set())

    def test_project_capitalization_in_runtime_context_fails(self) -> None:
        body = """
Notes
-----
  Runtime requirements:
    - Bash >= 4.4
    - SAMtools
""".removesuffix("\n")

        self.assertIn("HELP.RUNTIME.EXACT_CALLABLE", rule_ids(body))

    def test_conceptual_project_prose_is_excluded(self) -> None:
        body = """
Notes
-----
  Bash is the implementation language used by this project.
""".removesuffix("\n")

        self.assertEqual(rule_ids(body), set())

    def test_typed_expected_and_generated_globals_pass(self) -> None:
        body = """
Expected globals
----------------
  dir_eo : dir
    Existing log directory.

Generated globals
-----------------
  cmd_bld : array of str
    Constructed command arguments.
""".removesuffix("\n")

        self.assertEqual(rule_ids(body), set())

    def test_grouped_global_row_passes(self) -> None:
        body = """
Expected globals
----------------
  scr_sub, ref_fa : file
    Submission-script path and optional reference-FASTA path, respectively.
""".removesuffix("\n")

        self.assertEqual(rule_ids(body), set())

    def test_grouped_global_bad_comma_spacing_fails(self) -> None:
        body = """
Expected globals
----------------
  scr_sub,ref_fa : file
    Submission-script path and optional reference-FASTA path, respectively.
""".removesuffix("\n")

        self.assertIn("HELP.GLOBAL.GROUP_SYNTAX", rule_ids(body))

    def test_grouped_global_invalid_identifier_fails(self) -> None:
        body = """
Expected globals
----------------
  scr_sub, 9ref_fa : file
    Submission-script path and optional reference-FASTA path, respectively.
""".removesuffix("\n")

        self.assertIn("HELP.GLOBAL.GROUP_SYNTAX", rule_ids(body))

    def test_grouped_global_without_one_shared_type_fails(self) -> None:
        body = """
Expected globals
----------------
  scr_sub : file, ref_fa : path
    Submission-script path and optional reference-FASTA path, respectively.
""".removesuffix("\n")

        self.assertIn("HELP.GLOBAL.GROUP_SYNTAX", rule_ids(body))

    def test_repeated_same_type_and_description_requires_grouping(
        self,
    ) -> None:
        body = """
Expected globals
----------------
  scr_sub : file
    Input path.

  ref_fa : file
    Input path.
""".removesuffix("\n")

        self.assertIn("HELP.GLOBAL.REPEATED_DESCRIPTION", rule_ids(body))

    def test_same_type_with_different_descriptions_passes_separately(
        self,
    ) -> None:
        body = """
Expected globals
----------------
  scr_sub : file
    Submission-script path.

  ref_fa : file
    Optional reference-FASTA path.
""".removesuffix("\n")

        self.assertEqual(rule_ids(body), set())

    def test_grouped_global_ambiguous_mapping_is_reported_for_review(
        self,
    ) -> None:
        body = """
Expected globals
----------------
  scr_sub, ref_fa : file
    Input paths used by the workflow.
""".removesuffix("\n")

        self.assertIn("HELP.GLOBAL.MAPPING_REVIEW", rule_ids(body))

    def test_globals_without_types_fail(self) -> None:
        body = """
Expected globals
----------------
  dir_eo
    Existing log directory.

Generated globals
-----------------
  cmd_bld
    Constructed command arguments.
""".removesuffix("\n")

        self.assertIn("HELP.GLOBAL.TYPE", rule_ids(body))

    def test_globals_require_blank_line_between_entries(self) -> None:
        body = """
Expected globals
----------------
  dir_eo : dir
    Existing log directory.
  nam_job : str
    Job name prefix.
""".removesuffix("\n")

        self.assertIn("HELP.ENTRY.BLANK", rule_ids(body))

    def test_notes_pseudo_header_prose_at_plus_two_passes(self) -> None:
        body = """
Notes
-----
  Runtime requirements:
    bc
""".removesuffix("\n")

        self.assertEqual(rule_ids(body), set())

    def test_notes_pseudo_header_same_indent_prose_fails(self) -> None:
        body = """
Notes
-----
  Runtime requirements:
  bc
""".removesuffix("\n")

        self.assertIn("HELP.NOTES.PSEUDO_HEADER_INDENT", rule_ids(body))

    def test_top_level_help_content_uses_two_space_indentation(self) -> None:
        valid = """
Returns
-------
  Returns 0.
""".removesuffix("\n")
        invalid = """
Returns
-------
Returns 0.
""".removesuffix("\n")

        self.assertEqual(rule_ids(valid), set())
        self.assertIn("HELP.SECTION.INDENT", rule_ids(invalid))

    def test_single_unbulleted_runtime_requirement_passes(self) -> None:
        body = """
Notes
-----
  Runtime requirements:
    bc
""".removesuffix("\n")

        self.assertEqual(rule_ids(body), set())

    def test_single_bulleted_runtime_requirement_fails(self) -> None:
        body = """
Notes
-----
  Runtime requirements:
    - bc
""".removesuffix("\n")

        self.assertIn("HELP.RUNTIME.CARDINALITY", rule_ids(body))

    def test_multiple_bulleted_runtime_requirements_pass(self) -> None:
        body = """
Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - grep
""".removesuffix("\n")

        self.assertEqual(rule_ids(body), set())

    def test_dependency_pseudo_sections_fail(self) -> None:
        for heading in ("Dependency:", "Dependencies:"):
            with self.subTest(heading=heading):
                body = f"""
Returns
-------
  Returns 0.

{heading}
  - bc
""".removesuffix("\n")

                self.assertIn(
                    "HELP.SECTION.UNSUPPORTED_DEPENDENCY",
                    rule_ids(body),
                )

    def test_duplicate_notes_headings_fail(self) -> None:
        body = """
Notes
-----
  First note.

Notes
-----
  Second note.
""".removesuffix("\n")

        self.assertIn("HELP.SECTION.DUPLICATE", rule_ids(body))

    def test_any_duplicate_top_level_heading_fails(self) -> None:
        for heading, underline in (
            ("Usage", "-----"),
            ("Parameters", "----------"),
            ("Returns", "-------"),
            ("Output", "------"),
            ("Examples", "--------"),
        ):
            with self.subTest(heading=heading):
                body = f"""
{heading}
{underline}
  First.

{heading}
{underline}
  Second.
""".removesuffix("\n")

                self.assertIn("HELP.SECTION.DUPLICATE", rule_ids(body))

    def test_canonical_function_heading_order_passes(self) -> None:
        body = function_body(
            """
Usage
-----
  sample_func
    [--help] value
""".removesuffix("\n"),
            """
Parameters
----------
  1  value : str
    Value to inspect.
""".removesuffix("\n"),
            """
Expected globals
----------------
  mode : str
    Active mode.
""".removesuffix("\n"),
            """
Generated globals
-----------------
  result : str
    Generated value.
""".removesuffix("\n"),
            """
Returns
-------
  Returns 0.
""".removesuffix("\n"),
            """
Output
------
  Prints the generated value.
""".removesuffix("\n"),
            """
Notes
-----
  The value is inspected once.
""".removesuffix("\n"),
            """
Examples
--------
  '''bash
  sample_func value
  '''
""".removesuffix("\n"),
        )

        self.assertEqual(
            findings(body, owner="sample_func"),
            [],
        )

    def test_each_out_of_order_function_heading_pair_fails_precisely(
        self,
    ) -> None:
        ordered = [
            ("Usage", "-----", "  sample_func\n    [--help] value"),
            ("Parameters", "----------", "  1  value : str\n    Value."),
            (
                "Expected globals",
                "----------------",
                "  mode : str\n    Mode.",
            ),
            (
                "Generated globals",
                "-----------------",
                "  result : str\n    Result.",
            ),
            ("Returns", "-------", "  Returns 0."),
            ("Output", "------", "  Prints output."),
            ("Notes", "-----", "  One note."),
            ("Examples", "--------", "  '''bash\n  sample_func value\n  '''"),
        ]

        for index in range(len(ordered) - 1):
            swapped = list(ordered)
            swapped[index], swapped[index + 1] = (
                swapped[index + 1],
                swapped[index],
            )
            body = "\n\n".join(
                f"{name}\n{underline}\n{content}"
                for name, underline, content in swapped
            )

            with self.subTest(pair=(ordered[index][0], ordered[index + 1][0])):
                order_findings = [
                    item
                    for item in findings(body, owner="sample_func")
                    if item.rule_id == "HELP.SECTION.ORDER"
                ]

                self.assertTrue(order_findings)
                self.assertIn("must precede", order_findings[0].message)

    def test_function_help_without_returns_fails(self) -> None:
        body = """
Usage
-----
  sample_func
    [--help]
""".removesuffix("\n")

        self.assertIn(
            "HELP.SECTION.REQUIRED",
            {item.rule_id for item in findings(body, owner="sample_func")},
        )

    def test_returns_plus_output_are_distinct_and_pass(self) -> None:
        body = function_body(
            """
Usage
-----
  sample_func
    [--help]
""".removesuffix("\n"),
            """
Returns
-------
  Returns 0.
""".removesuffix("\n"),
            """
Output
------
  Prints one record.
""".removesuffix("\n"),
        )

        self.assertEqual(findings(body, owner="sample_func"), [])

    def test_output_is_not_a_returns_duplicate(self) -> None:
        body = function_body(
            """
Usage
-----
  sample_func
    [--help]
""".removesuffix("\n"),
            """
Returns
-------
  Returns 0.
""".removesuffix("\n"),
            """
Output
------
  Prints one record.
""".removesuffix("\n"),
        )

        self.assertNotIn(
            "HELP.SECTION.DUPLICATE",
            {item.rule_id for item in findings(body, owner="sample_func")},
        )

    def test_examples_last_passes_and_later_section_fails(self) -> None:
        valid = function_body(
            """
Usage
-----
  sample_func
    [--help]
""".removesuffix("\n"),
            """
Returns
-------
  Returns 0.
""".removesuffix("\n"),
            """
Examples
--------
  '''bash
  sample_func
  '''
""".removesuffix("\n"),
        )
        invalid = valid + "\n\nNotes\n-----\n  Too late."

        self.assertEqual(findings(valid, owner="sample_func"), [])
        self.assertIn(
            "HELP.SECTION.ORDER",
            {item.rule_id for item in findings(invalid, owner="sample_func")},
        )

    def test_usage_with_only_long_options_passes(self) -> None:
        body = """
Usage
-----
  print_banner_pretty
    [--help] --text <str> [--wrap <str>] [--pad <int>] [--cols <int>]
""".removesuffix("\n")

        self.assertEqual(rule_ids(body), set())

    def test_usage_advertising_short_and_long_aliases_fails(self) -> None:
        body = """
Usage
-----
  print_banner_pretty
    [--help] -tx|--text <str> [-w|--wrap <str>]
""".removesuffix("\n")

        self.assertIn("HELP.DOC.LONG_OPTIONS", rule_ids(body))

    def test_hidden_parser_aliases_remain_allowed_with_long_help(self) -> None:
        body = """
Usage
-----
  print_banner_pretty
    [--help] --text <str>
""".removesuffix("\n")
        parser = (
            'function parse_args() { case "${1:-}" in -tx|--txt|--text) '
            "return 0 ;; esac; }"
        )

        self.assertEqual(rule_ids(body, parser), set())

    def test_flat_typed_globals_pass(self) -> None:
        body = """
Expected globals
----------------
  PY_INVOKE : {'module', 'file'}
    Python invocation mode.

  dir_rep : dir
    Optional repository root.
""".removesuffix("\n")

        self.assertEqual(rule_ids(body), set())

    def test_read_subgroup_beneath_expected_globals_fails(self) -> None:
        body = """
Expected globals
----------------
  Read:
    PY_INVOKE : {'module', 'file'}
      Python invocation mode.
""".removesuffix("\n")

        self.assertIn("HELP.GLOBAL.SUBGROUP", rule_ids(body))

    def test_embedded_tables_use_canonical_markdown_pipe_source(self) -> None:
        canonical = """
Notes
-----
  | Name | Meaning      |
  | :--- | :---         |
  | foo  | first value  |
  | bar  | second value |
""".removesuffix("\n")
        noncanonical = """
Notes
-----
  Name | Meaning
  -----|--------
  foo  | first value
  bar  | second value
""".removesuffix("\n")

        self.assertNotIn("MD.TABLE.CANONICAL", rule_ids(canonical))
        self.assertIn("MD.TABLE.CANONICAL", rule_ids(noncanonical))

    def test_python_docstring_tables_use_canonical_pipe_source(self) -> None:
        canonical = '''
def sample() -> None:
    """
    Describe one table.

    Notes
    -----
    | Name | Meaning      |
    | :--- | :---         |
    | foo  | first value  |
    | bar  | second value |
    """
'''
        noncanonical = '''
def sample() -> None:
    """
    Describe one table.

    Notes
    -----
    Name | Meaning
    -----|--------
    foo  | first value
    bar  | second value
    """
'''

        canonical_rules = {
            item.rule_id
            for item in check_python_docstrings(canonical, "sample.py", None)
        }
        noncanonical_rules = {
            item.rule_id
            for item in check_python_docstrings(
                noncanonical,
                "sample.py",
                None,
            )
        }

        self.assertNotIn("MD.TABLE.CANONICAL", canonical_rules)
        self.assertIn("MD.TABLE.CANONICAL", noncanonical_rules)

    def test_case_arm_spacing_advisory_is_available_and_bounded(self) -> None:
        self.assertTrue(hasattr(help_style, "case_arm_spacing_warnings"))

        if not hasattr(help_style, "case_arm_spacing_warnings"):
            return

        text = """
case \"${value}\" in
    loaded)
        return 0
        ;;
    loading)
        return 0
        ;;
esac
"""

        warnings = help_style.case_arm_spacing_warnings(text, "fixture.sh")

        self.assertEqual(len(warnings), 1)

    def test_clear_prose_backticks_fail(self) -> None:
        body = """
Notes
-----
  Use `module` mode for dotted imports.
""".removesuffix("\n")

        self.assertIn("HELP.TOKEN.QUOTING.SHELL", rule_ids(body))

    def test_literal_backtick_in_examples_is_excluded(self) -> None:
        body = """
Examples
--------
  '''bash
  printf '%s\\n' '`literal syntax`'
  '''
""".removesuffix("\n")

        self.assertNotIn("HELP.PROSE.STRAIGHT_QUOTES", rule_ids(body))

    def test_rendered_checker_reuses_multiline_command_auditor(self) -> None:
        body = """
Examples
--------
  1. Demonstrate collapsed output.
    '''bash
    bash demo.sh \\      --mode one
    '''
""".removesuffix("\n")

        observed = {
            item.rule_id
            for item in check_rendered_help(
                body,
                registry_path=REGISTRY,
            )
        }

        self.assertIn("HELP.EXAMPLES.COMMAND.COLLAPSED", observed)
        self.assertIn("HELP.EXAMPLES.COMMAND.RENDERED_STRUCTURE", observed)

    def test_python_docstring_prose_backticks_fail(self) -> None:
        for token in ("`module`", "``module``"):
            with self.subTest(token=token):
                text = f'''
def sample():
    """Use {token} mode for dotted imports."""
'''

                self.assertTrue(
                    check_python_docstrings(text, "sample.py", None),
                )

    def test_python_docstring_straight_single_quotes_pass(self) -> None:
        text = '''
def sample():
    """Use '--write' to update selected files."""
'''

        self.assertEqual(check_python_docstrings(text, "sample.py", None), [])

    def test_python_docstring_example_backticks_are_excluded(self) -> None:
        text = '''
def sample():
    """Show literal syntax.

    Examples
    --------
    >>> print("`literal syntax`")
    """
'''

        self.assertEqual(check_python_docstrings(text, "sample.py", None), [])

    def test_docstring_doctest_continuations_and_tilde_fences_are_excluded(
        self,
    ) -> None:
        text = '''
def sample():
    """Show literal syntax.

    ... print("`doctest continuation`")

    ~~~text
    `fenced literal`
    ~~~
    """
'''

        self.assertEqual(check_python_docstrings(text, "sample.py", None), [])

    def test_canonical_unknown_option_error_passes(self) -> None:
        parser = """
function parse_args() {
    case \"${1:-}\" in
        *) echo_err \"unknown option/parameter passed: '${1}'.\"; return 1 ;;
    esac
}
""".removesuffix("\n")

        self.assertNotIn(
            "SHELL.PARSER.UNKNOWN_ERROR",
            rule_ids("Notes\n-----\n  Parser fixture.", parser),
        )

    def test_obsolete_unknown_argument_banners_fail(self) -> None:
        for punctuation in ("", "."):
            with self.subTest(punctuation=punctuation):
                parser = (
                    "function parse_args() {\n"
                    '    case "${1:-}" in\n'
                    '        *) echo "## Unknown argument passed: '
                    f"'${{1}}'{punctuation} ##\" >&2; return 1 ;;\n"
                    "    esac\n"
                    "}"
                )

                self.assertIn(
                    "SHELL.PARSER.UNKNOWN_ERROR",
                    rule_ids("Notes\n-----\n  Parser fixture.", parser),
                )

    def test_manually_corrected_source_units_pass(self) -> None:
        for relative in (
            "lib/bash/workflows/align_fastqs.sh",
            "lib/bash/workflows/calculate_scaling_factor.sh",
            "lib/bash/core/check_args.sh",
            "lib/bash/core/format_outputs.sh",
            "lib/bash/core/handle_env.sh",
            "lib/bash/dispatch/run_python.sh",
            "lib/bash/core/wrap_cmd.sh",
        ):
            with self.subTest(path=relative):
                self.assertEqual(
                    scan_repository(
                        REPO_ROOT,
                        [relative],
                        registry_path=REGISTRY,
                    ),
                    [],
                )


if __name__ == "__main__":
    unittest.main()
