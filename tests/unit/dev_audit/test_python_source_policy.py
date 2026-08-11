#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_python_source_policy.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.6);
# - Anthropic Claude Code (Opus 5).
#
# Distributed under the MIT license.


"""
Test deterministic portions of the bounded Python source-policy checker.
"""

from __future__ import annotations

import ast
import hashlib
import tomllib
from pathlib import Path

import pytest
from dev.audit.python_naming_vocabulary import (
    CONTEXTS,
    MATCH_KINDS,
    STATUSES,
    evidence_candidate_segments,
    load_vocabulary,
    ordinary_short_words,
    prohibited_internal_segments,
)
from dev.audit.python_source_policy import (
    RULE_ANNOTATIONS,
    RULE_CLI_HELP_LAYOUT,
    RULE_COMMENTS,
    RULE_DOC_LAYOUT,
    RULE_DOC_NUMPY,
    RULE_HELP_SENTENCES,
    RULE_MULTILINE,
    RULE_NAMING,
    RULE_PROSE_WRAP,
    RULE_STRINGS,
    RULE_TOPOLOGY,
    analyze_text,
    parse_args,
)
from dev.tools.markdown_format import _rebase_details
from dev.tools.python_help_format import format_source

ROOT = Path(__file__).resolve().parents[3]
FIXTURES = ROOT / "tests" / "fixtures" / "python_source_policy"
CORRECTION_TOOLS = (
    ROOT
    / "artifacts"
    / "tests"
    / "python_source_migration_20260726_corrected"
    / "tools"
)


def test_governed_naming_vocabulary_preserves_policy_projections() -> None:
    """
    Preserve accepted, protected, prohibited, and review-only vocabulary.
    """

    vocabulary_path = ROOT / "dev" / "config" / "python_naming_vocabulary.toml"
    raw_vocabulary = tomllib.loads(
        vocabulary_path.read_text(encoding="utf-8"),
    )
    entries = load_vocabulary(vocabulary_path)
    prohibited = {
        entry.spelling: entry.replacement
        for entry in entries
        if entry.status == "prohibited_internal_shorthand"
    }
    protected = {
        entry.spelling
        for entry in entries
        if entry.status == "protected_external_or_interface_spelling"
    }

    def projection_hash(values: frozenset[str]) -> str:
        payload = "".join(f"{value}\n" for value in sorted(values))
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()

    assert raw_vocabulary["schema_version"] == 2
    assert len(entries) == 123
    assert {entry.status for entry in entries} == set(STATUSES)
    assert {entry.match_kind for entry in entries} == set(MATCH_KINDS)
    assert {context for entry in entries for context in entry.contexts} == set(
        CONTEXTS
    )
    assert all(entry.contexts for entry in entries)
    assert all(
        entry.evidence_candidate
        == (
            entry.match_kind == "casefold_segment"
            and entry.status
            in {
                "allowed_domain_term",
                "prohibited_internal_shorthand",
                "review_candidate",
            }
        )
        for entry in entries
    )
    assert all(
        set(raw_entry)
        == {
            "spelling",
            "status",
            "match_kind",
            "contexts",
            "evidence_candidate",
        }
        | (
            {"replacement"}
            if raw_entry["status"] == "prohibited_internal_shorthand"
            else set()
        )
        for raw_entry in raw_vocabulary["entry"]
    )

    ordinary = ordinary_short_words()
    evidence = evidence_candidate_segments()
    prohibited_projection = prohibited_internal_segments()

    assert len(ordinary) == 25
    assert len(evidence) == 77
    assert prohibited_projection == frozenset(prohibited)
    assert projection_hash(ordinary) == (
        "5d3389826ca59f5e638b14ae80314de7263c20ea2581cebea26c66b8c41d9629"
    )
    assert projection_hash(evidence) == (
        "2b9aabb7c42cf07089cf5bafc47b6644d86815f4a449466cd61816e1fb42bd0b"
    )
    assert projection_hash(prohibited_projection) == (
        "8bc06877dbeff1e556833a33cdb192356feab95bff987b4764f822e6f6380188"
    )
    assert prohibited == {
        "cfg": "configuration",
        "cmb": "combined",
        "col": "column",
        "cvg": "coverage",
        "dp": "decimal_places",
        "ext": "extension",
        "fh": "handle",
        "fmt": "format",
        "py": "python",
        "rc": "return_code",
        "sb": "start_bins",
        "src": "source",
        "str": "text",
        "xs": "values",
        "ys": "values",
    }
    assert protected == {
        "fil_in",
        "fil_out",
        "siz_bin",
        "skp_pfx",
        "scl_fct",
        "frg",
        "psc",
        "qntl",
    }


def test_comment_registry_preserves_shared_delegation() -> None:
    """
    Preserve the shared-to-Python comment checker delegation.
    """

    registry = tomllib.loads(
        (ROOT / "dev" / "config" / "rules.toml").read_text(
            encoding="utf-8",
        ),
    )
    entries = {entry["rule_id"]: entry for entry in registry["rule"]}

    assert "SOURCE.COMMENT.ATTACHMENT" not in entries

    python_entry = entries["PY.COMMENT.FORM"]
    assert python_entry["command"] == [
        "conda",
        "run",
        "-n",
        "env_protocol",
        "python",
        "-m",
        "dev.audit.python_source_policy",
        "--root",
        ".",
        "--all-maintained",
        "--rule",
        "PY.COMMENT.FORM",
    ]
    assert set(python_entry["applicable_paths"]) == {
        "docs/standards/source_layout.md",
        "docs/standards/python.md",
        "dev/audit/python_source_policy.py",
        "tests/fixtures/python_source_policy/**",
        "tests/unit/dev_audit/test_python_source_policy.py",
        "src/**/*.py",
        "tests/**/*.py",
        "dev/**/*.py",
        (
            "artifacts/tests/"
            "python_source_migration_20260727_semantic_review/tools/**/*.py"
        ),
    }
    assert "SOURCE.COMMENT.ATTACHMENT" in python_entry["covered_scope"]
    assert "SOURCE.COMMENT.ATTACHMENT" in python_entry["remaining_scope"]

    shared_owner = (
        ROOT / "docs" / "standards" / "source_layout.md"
    ).read_text(encoding="utf-8")
    router = (ROOT / "docs" / "standards" / "README.md").read_text(
        encoding="utf-8",
    )

    assert "No dedicated registry entry exists for this shared owner." in (
        shared_owner
    )
    assert "registered `PY.COMMENT.FORM` execution delegates" in shared_owner
    assert "source_layout.md#comment-attachment" in router
    assert "--rule PY.COMMENT.FORM <paths>" in router


def rule_messages(text: str, rule_id: str) -> list[str]:
    """
    Return messages for one rule and source string.
    """

    return [
        finding.message
        for finding in analyze_text(text, "<fixture>").findings
        if finding.rule_id == rule_id
    ]


def help_source(*parts: str) -> str:
    """
    Build one canonical parse_args source around literal parts.
    """

    literals = "\n".join(f'            "{part}"' for part in parts)
    return f'''
"""
Provide one help-layout fixture.
"""

import argparse


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse fixture arguments.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--value",
        help=(
{literals}
        ),
    )
    return parser.parse_args(argv)
'''


def help_value(text: str) -> str:
    """
    Return the constant value of one fixture help keyword.
    """

    tree = ast.parse(text)
    return next(
        node.value.value
        for node in ast.walk(tree)
        if isinstance(node, ast.keyword)
        and node.arg == "help"
        and isinstance(node.value, ast.Constant)
    )


def test_fixture_cohorts_cover_positive_negative_and_exceptions() -> None:
    """
    Keep canonical and literal-exception fixtures clean and negatives broad.
    """

    for name in (
        "accepted/canonical.py",
        "accepted/exceptions.py",
        "boundary/source_form.py",
    ):
        text = (FIXTURES / name).read_text(encoding="utf-8")

        assert analyze_text(text, name).findings == ()

    negative = "rejected/deterministic_owners.py"
    text = (FIXTURES / negative).read_text(encoding="utf-8")
    rule_ids = {
        finding.rule_id for finding in analyze_text(text, negative).findings
    }

    # RULE_NAMING and RULE_TOPOLOGY are asserted here rather than in the
    # fixture recipe. The recipe used to run the checker over its own output
    # and require eleven owners, nine of which this set already covered;
    # 'make.sh' generates, and deciding which owners a fixture must provoke is
    # a unit-level property of the checker.
    assert {
        RULE_HELP_SENTENCES,
        RULE_ANNOTATIONS,
        RULE_CLI_HELP_LAYOUT,
        RULE_COMMENTS,
        RULE_DOC_LAYOUT,
        RULE_DOC_NUMPY,
        RULE_MULTILINE,
        RULE_NAMING,
        RULE_PROSE_WRAP,
        RULE_STRINGS,
        RULE_TOPOLOGY,
    } <= rule_ids


def test_correction_tools_do_not_exempt_their_own_source_form() -> None:
    """
    Apply the maintained Python policy to every successor correction tool.
    """

    paths = sorted(CORRECTION_TOOLS.glob("*.py"))
    findings: dict[str, tuple[object, ...]] = {}

    for path in paths:
        analysis = analyze_text(
            path.read_text(encoding="utf-8"),
            path.relative_to(ROOT).as_posix(),
        )
        # These tools are immutable migration evidence under 'artifacts/', so
        # they answer to the facets that existed when they were accepted.
        # 'SOURCE.PROSE.WRAP' arrived afterward and is advisory; its residual
        # belongs to the recorded migration rather than to this contract.
        retained = tuple(
            finding
            for finding in analysis.findings
            if finding.rule_id != RULE_PROSE_WRAP
        )

        if retained:
            findings[path.name] = retained

    assert paths
    assert findings == {}


def test_details_rebasing_preserves_relative_heading_levels() -> None:
    """
    Keep heading rebasing identical while retaining successful matches.
    """

    source = [
        "## Parent",
        "<details>",
        "<summary>More</summary>",
        "## Child",
        "### Grandchild",
        "</details>",
    ]
    expected = [
        "## Parent",
        "<details>",
        "<summary>More</summary>",
        "### Child",
        "#### Grandchild",
        "</details>",
    ]

    assert _rebase_details(source) == expected


def test_docstring_layout_prefix_summary_and_spacing_boundaries() -> None:
    """
    Reject bad summaries and unnecessary raw prefixes.
    """

    no_full_stop = '''
"""
Describe this module
"""
'''
    unnecessary_raw = r'''
r"""
Describe this module.
"""
'''

    summary_messages = rule_messages(
        no_full_stop,
        RULE_DOC_LAYOUT,
    )
    raw_messages = rule_messages(unnecessary_raw, RULE_DOC_LAYOUT)

    assert "docstring summary must end with a full stop" in summary_messages
    assert "raw docstring prefix requires literal backslash content" in (
        raw_messages
    )


def test_numpy_docstring_sections_accept_canonical_rendering() -> None:
    """
    Accept exact headings, order, underlines, and body indentation.
    """

    source = '''
"""
Provide one module.
"""


def render(value: str) -> str:
    """
    Render one value.

    Parameters
    ----------
    value : str
        Value to render.

    Returns
    -------
    rendered : str
        Rendered value.

    Notes
    -----
    The fixture keeps each section at canonical indentation.
    """

    return value
'''

    assert rule_messages(source, RULE_DOC_NUMPY) == []


def test_numpy_multiline_type_closer_has_one_diagnostic_owner() -> None:
    """
    Align textual type closers and report only through the NumPy adapter.
    """

    accepted = '''
"""
Provide one module.
"""


def render() -> tuple[int, str]:
    """
    Render two values.

    Returns
    -------
    result : tuple[
        int, str
    ]
        Rendered values.
    """

    return 1, "value"
'''
    rejected = accepted.replace("\n    ]\n", "\n        ]\n")
    expected = (
        "multiline NumPy type closing delimiter must align with the entry line"
    )

    assert rule_messages(accepted, RULE_DOC_NUMPY) == []
    assert rule_messages(rejected, RULE_DOC_NUMPY) == [expected]
    assert expected not in rule_messages(rejected, RULE_MULTILINE)


def test_numpy_docstring_sections_require_typed_named_entries() -> None:
    """
    Reject untyped parameters and unnamed return or yield values.
    """

    source = '''
"""
Provide one module.
"""


def render(value: str):
    """
    Render one value.

    Parameters
    ----------
    value
        Value to render.

    Returns
    -------
    str
        Rendered value.

    Yields
    ------
    bytes
        Encoded values.
    """

    yield value.encode()
'''

    messages = rule_messages(source, RULE_DOC_NUMPY)

    assert "NumPy section 'Parameters' entry requires 'name : type'" in (
        messages
    )
    assert "NumPy section 'Returns' entry requires 'name : type'" in messages
    assert "NumPy section 'Yields' entry requires 'name : type'" in messages


def test_numpy_docstring_types_match_callable_annotations() -> None:
    """
    Reject documented types that disagree with callable annotations.
    """

    source = '''
"""
Provide one module.
"""


def render(value: str) -> int:
    """
    Render one value.

    Parameters
    ----------
    value : bytes
        Value to render.

    Returns
    -------
    rendered : str
        Rendered value.
    """

    return len(value)
'''

    messages = rule_messages(source, RULE_DOC_NUMPY)

    mismatch = "does not match annotation"
    parameter_message = (
        f"NumPy section 'Parameters' type 'bytes' {mismatch} 'str'"
    )
    return_message = f"NumPy section 'Returns' type 'str' {mismatch} 'int'"

    assert parameter_message in messages
    assert return_message in messages


def test_numpy_docstring_names_match_stable_direct_returns() -> None:
    """
    Require stable direct result names without guessing computed identities.
    """

    source = '''
"""
Provide one module.
"""


def prepare(header: list[str]) -> tuple[list[str], list[dict[str, str]]]:
    """
    Prepare one header and its rows.

    Parameters
    ----------
    header : list[str]
        Input header.

    Returns
    -------
    result : tuple[list[str], list[dict[str, str]]]
        Prepared header and rows.
    """

    rows: list[dict[str, str]] = []

    return header, rows
'''

    messages = rule_messages(source, RULE_DOC_NUMPY)

    assert (
        "NumPy section 'Returns' names 'result' do not match stable returned "
        "names 'header, rows'"
    ) in messages


def test_numpy_docstring_names_do_not_guess_branch_dependent_results() -> None:
    """
    Leave branch-dependent result identities to semantic review.
    """

    source = '''
"""
Provide one module.
"""


def choose(use_first: bool, first: str, second: str) -> str:
    """
    Choose one value.

    Parameters
    ----------
    use_first : bool
        Whether to select the first value.
    first : str
        First candidate.
    second : str
        Second candidate.

    Returns
    -------
    selected : str
        Selected candidate.
    """

    if use_first:
        return first

    return second
'''

    assert rule_messages(source, RULE_DOC_NUMPY) == []


def test_numpy_parameters_match_annotated_signature_membership() -> None:
    """
    Reject unknown entries and omitted annotated parameters.
    """

    source = '''
"""
Provide one module.
"""


def render(value: str, mode: str) -> str:
    """
    Render one value.

    Parameters
    ----------
    item : str
        Value to render.

    Returns
    -------
    rendered : str
        Rendered value.
    """

    return value if mode else ""
'''

    messages = rule_messages(source, RULE_DOC_NUMPY)

    assert (
        "NumPy section 'Parameters' documents unknown parameter(s): item"
    ) in messages
    assert (
        "NumPy section 'Parameters' omits annotated parameter(s): mode, value"
    ) in messages


def test_numpy_docstring_sections_reject_malformed_boundaries() -> None:
    """
    Reject pseudo-headings, bad underlines, order, and indentation.
    """

    source = '''
"""
Provide one module.
"""


def render(value: str) -> str:
    """
    Render one value.

    Returns
    ----
        str
            Rendered value.

    Parameters
    ----------
    value : str
        Value to render.

    Exits:
        Never.
    """

    return value
'''

    messages = rule_messages(source, RULE_DOC_NUMPY)

    assert "NumPy section 'Returns' requires an exact underline" in messages
    assert "NumPy section 'Returns' body is overindented" in messages
    assert "NumPy section 'Parameters' is out of canonical order" in messages
    assert "noncanonical NumPy pseudo-section 'Exits:'" in messages


def test_numpy_docstring_sections_require_blank_boundaries() -> None:
    """
    Reject a section heading attached directly to prior section prose.
    """

    source = '''
"""
Provide one module.
"""


def render(value: str) -> str:
    """
    Render one value.

    Parameters
    ----------
    value : str
        Value to render.
    Returns
    -------
    rendered : str
        Rendered value.
    """

    return value
'''

    messages = rule_messages(source, RULE_DOC_NUMPY)

    expected_message = (
        "NumPy section 'Returns' requires a blank line before its heading"
    )

    assert expected_message in messages


def test_numpy_docstring_sections_reject_empty_and_overlong_summary() -> None:
    """
    Reject an empty section and expose summary-line width candidates.
    """

    source = (
        '"""\nProvide one module.\n"""\n\n\n'
        "def render() -> None:\n"
        '    """\n'
        f"    {'Describe one deliberately overlong summary ' * 2}.\n"
        "\n"
        "    Notes\n"
        "    -----\n"
        '    """\n'
        "\n"
        "    return None\n"
    )

    messages = rule_messages(source, RULE_DOC_NUMPY)

    assert "docstring summary exceeds 79 columns" in messages
    assert "NumPy section 'Notes' requires substantive content" in messages


def test_adjacent_prose_literals_avoid_orphaned_words() -> None:
    """
    Reject a safely movable final prose word while preserving its value.
    """

    fragmented = """
def render(fact_count: int) -> str:
    return (
        f"No deterministic finding; {fact_count} normalized "
        f"evidence record(s) support that "
        f"result."
    )
"""
    canonical = """
def render(fact_count: int) -> str:
    return (
        f"No deterministic finding; {fact_count} normalized "
        f"evidence record(s) support that result."
    )
"""

    expected = (
        "adjacent prose literal breaks before a word or punctuation fragment "
        "that would still fit within 79 columns"
    )

    assert expected in rule_messages(fragmented, RULE_MULTILINE)
    assert rule_messages(canonical, RULE_MULTILINE) == []

    fragmented_namespace: dict[str, object] = {}
    canonical_namespace: dict[str, object] = {}
    exec(fragmented, fragmented_namespace)
    exec(canonical, canonical_namespace)
    fragmented_render = fragmented_namespace["render"]
    canonical_render = canonical_namespace["render"]

    assert callable(fragmented_render)
    assert callable(canonical_render)
    assert fragmented_render(3) == canonical_render(3)

    structured = """
def render() -> str:
    return (
        "Notes\\n"
        "-----\\n"
    )
"""

    assert rule_messages(structured, RULE_MULTILINE) == []

    orphaned_punctuation = """
def render() -> str:
    return (
        "Done"
        "."
    )
"""

    expected = (
        "adjacent prose literal breaks before a word or punctuation fragment "
        "that would still fit within 79 columns"
    )

    assert expected in rule_messages(orphaned_punctuation, RULE_MULTILINE)


def test_docstring_indentation_and_single_token_are_exact() -> None:
    """
    Reject overindented content and implicitly concatenated docstrings.
    """

    overindented = '''
"""
Provide one module.
"""


def run() -> None:
    """
        Run one value.
        """

    return None
'''

    messages = rule_messages(overindented, RULE_DOC_LAYOUT)

    assert "docstring summary must align with the opener" in messages
    assert "closing docstring delimiter must align with the opener" in messages

    concatenated = '''
(
"""
Provide one module.
"""
"""
Add implicit content.
"""
)
'''

    assert "docstring must be one triple-double-quoted token" in rule_messages(
        concatenated,
        RULE_DOC_LAYOUT,
    )


def test_annotations_exclude_only_actual_method_receivers() -> None:
    """
    Keep a real receiver exempt without exempting a free self parameter.
    """

    method = '''
"""
Provide one receiver fixture.
"""


class Value:
    """
    Hold one value.
    """

    def get(self) -> int:
        """
        Return one value.
        """

        return 1
'''
    free = '''
"""
Provide one free-function fixture.
"""


def get(self) -> int:
    """
    Return one value.
    """

    return 1
'''

    method_messages = rule_messages(method, RULE_ANNOTATIONS)
    annotation_messages = rule_messages(
        free,
        RULE_ANNOTATIONS,
    )

    assert method_messages == []
    assert "parameter 'self' requires an annotation" in annotation_messages[0]


def test_annotations_distinguish_class_static_and_free_parameters() -> None:
    """
    Exclude only the receiver implied by the callable's actual role.
    """

    source = '''
"""
Provide receiver-role fixtures.
"""


class Value:
    """
    Hold one value.
    """

    @classmethod
    def from_value(cls) -> "Value":
        """
        Return one value.
        """

        return cls()

    @staticmethod
    def static(self) -> int:
        """
        Return one value.
        """

        return 1
'''

    messages = rule_messages(source, RULE_ANNOTATIONS)

    assert messages == [
        "<module>.Value.static parameter 'self' requires an annotation",
    ]


def test_string_quote_and_literal_content_boundaries() -> None:
    """
    Distinguish canonical quotes from bounded literal-content exceptions.
    """

    good = (
        '"""\n'
        "Provide one string fixture.\n"
        '"""\n\n'
        "value = 'Contains \"double quotes\".'\n"
        "block = '''\n"
        'Contains a """ delimiter.\n'
        "'''\n"
    )
    bad = (
        '"""\n'
        "Provide one string fixture.\n"
        '"""\n\n'
        "value = 'plain'\n"
        "block = '''one line with "
        '""" inside'
        "'''\n"
    )

    assert rule_messages(good, RULE_STRINGS) == []
    assert len(rule_messages(bad, RULE_STRINGS)) == 2


def test_comments_cover_markers_separators_inline_spacing_and_width() -> None:
    """
    Exercise exact comment forms and 79/80-column boundaries.
    """

    accepted = "    # " + ("x" * 73)
    rejected = "    # " + ("x" * 74)

    assert len(accepted) == 79
    assert len(rejected) == 80

    good = f'''
"""
Provide one comment fixture.
"""


def run() -> None:
    """
    Run one fixture.
    """

{accepted}
    #
    # First paragraph.
    #
    # Second paragraph.
    # This wrapped sentence ends on its short
    # continuation.
    #
    # value = computed + 1
    value = 1  # Inline sentence.
    _ = value
'''
    bad = f'''
"""
Provide one comment fixture.
"""


def run() -> None:
    """
    Run one fixture.
    """

{rejected}
    #
    # lowercase prose has no full stop
    value = 1 # Bad inline spacing.
    #
    _ = value  # Inline prose has no full stop
'''

    good_messages = rule_messages(good, RULE_COMMENTS)
    messages = rule_messages(bad, RULE_COMMENTS)
    expected = "bare '#' is allowed only between full-line comment paragraphs"
    punctuation_count = messages.count(
        "ordinary comment prose must end with terminal punctuation",
    )

    assert good_messages == []
    assert "ordinary comment source line exceeds 79 columns" in messages
    assert expected in messages
    assert "inline comment requires exactly two preceding spaces" in messages
    assert (
        "ordinary comment prose must begin with a capital letter" in messages
    )
    assert punctuation_count == 2


def test_comment_prose_uses_one_space_between_sentences() -> None:
    """
    Reject repeated sentence spacing without duplicating marker findings.
    """

    accepted = '''
"""
Provide one comment fixture.
"""

# First sentence. Second sentence.
VALUE = 1
'''
    rejected = accepted.replace(
        "First sentence. Second",
        "First sentence.  Second",
    )

    assert rule_messages(accepted, RULE_COMMENTS) == []
    assert rule_messages(rejected, RULE_COMMENTS) == [
        "ordinary comment prose must use one space between sentences",
    ]


def docstring_source(*body: str) -> str:
    """
    Return one module whose callable docstring holds the supplied body.
    """

    rows = "\n".join(f"    {line}" if line else "" for line in body)
    return f'def subject() -> None:\n    """\n{rows}\n    """\n\n    return None\n'


def test_docstring_prose_wraps_greedily_through_column_79() -> None:
    """
    Require the next whole word to move whenever it still fits.
    """

    filler = "w" * 60
    accepted = docstring_source(
        "Summary.",
        "",
        f"{filler} abcdefghijklmn",
        "x.",
    )
    premature = docstring_source("Summary.", "", f"{filler} abcdefghijk", "x.")
    accepted_line = next(
        line for line in accepted.splitlines() if filler in line
    )
    premature_line = next(
        line for line in premature.splitlines() if filler in line
    )

    assert len(accepted_line) == 79
    assert len(premature_line) == 76
    assert rule_messages(accepted, RULE_PROSE_WRAP) == []
    assert rule_messages(premature, RULE_PROSE_WRAP) == [
        "docstring prose breaks before a word that would still fit within 79 "
        "columns",
    ]


def test_docstring_structural_boundaries_are_not_prose_breaks() -> None:
    """
    Exclude entry headers, dedents, textual types, and doctest rows.
    """

    structural = docstring_source(
        "Summary.",
        "",
        "Parameters",
        "----------",
        "first : int",
        "    A short description.",
        "second : int",
        "    A short description.",
        "",
        "Returns",
        "-------",
        "result : tuple[",
        "    int,",
        "]",
        "    A short description.",
        "",
        "Examples",
        "--------",
        ">>> subject()",
        "None",
    )

    assert rule_messages(structural, RULE_PROSE_WRAP) == []


def test_aligned_list_continuations_keep_the_fit_test() -> None:
    """
    Judge an aligned list item by its wrap point, not its indentation.

    A list item whose continuation is aligned under the item text rather than
    under the marker's own continuation column is still wrapped prose. The
    alignment is an indentation choice; only the break itself is checked.
    """

    premature = docstring_source(
        "Summary.",
        "",
        "- alpha: " + ("w" * 57),
        "               tail.",
    )
    accepted = docstring_source(
        "Summary.",
        "",
        "- beta: " + ("w" * 63),
        "              tail.",
    )
    premature_line = next(
        line for line in premature.splitlines() if "alpha" in line
    )
    accepted_line = next(
        line for line in accepted.splitlines() if "beta" in line
    )

    assert len(premature_line) == 70
    assert len(accepted_line) == 75
    assert rule_messages(accepted, RULE_PROSE_WRAP) == []
    assert rule_messages(premature, RULE_PROSE_WRAP) == [
        "docstring prose breaks before a word that would still fit within 79 "
        "columns",
    ]


def test_quoted_units_spanning_words_are_indivisible() -> None:
    """
    Keep a break that would split one quoted multi-word unit unreported.

    A following word that opens a quote without closing it begins a quoted
    formula or condition. Filling it forward would move only its first word and
    split the unit, so the break is not reported even though the word fits. A
    complete quoted token is a different case: moving it splits nothing, so it
    is reported like any other word.
    """

    spanning = docstring_source(
        "Summary.",
        "",
        "A paragraph line that ends before a quoted " + ("w" * 20),
        "'a / b' is the unit.",
    )
    complete = docstring_source(
        "Summary.",
        "",
        "A paragraph line that ends before a quoted " + ("w" * 20),
        "'ab' is the unit.",
    )

    assert rule_messages(spanning, RULE_PROSE_WRAP) == []
    assert rule_messages(complete, RULE_PROSE_WRAP) == [
        "docstring prose breaks before a word that would still fit within 79 "
        "columns",
    ]


def test_hyphen_broken_words_rejoin_without_a_separator() -> None:
    """
    Measure a mid-word hyphen break without charging for a space.

    A trailing hyphen touching the preceding character splits one word, so
    rejoining inserts no space and the width test must not reserve one. A
    hyphen preceded by a space is a minus sign or dash and keeps its separator.
    """

    # Rejoining the broken word gives exactly 79 columns without a separator
    # and 80 with one, so the two readings disagree on this line.
    mid_word = docstring_source(
        "Summary.",
        "",
        ("w" * 54) + " whitespace-",
        "separated columns.",
    )
    operator = docstring_source(
        "Summary.",
        "",
        ("w" * 64) + " -",
        "separated columns.",
    )
    mid_word_line = next(
        line for line in mid_word.splitlines() if "whitespace-" in line
    )
    operator_line = next(
        line for line in operator.splitlines() if line.rstrip().endswith(" -")
    )

    # A suspended hyphen joins two words sharing one hyphen. Rejoining it would
    # produce 'zero-or', so width cannot decide it.
    suspended = docstring_source(
        "Summary.",
        "",
        ("w" * 54) + " whitespace-",
        "or negative columns.",
    )

    assert len(mid_word_line) == len(operator_line) == 70
    assert len(mid_word_line) + len("separated") == 79
    assert len(operator_line) + 1 + len("separated") == 80
    assert rule_messages(suspended, RULE_PROSE_WRAP) == []
    assert rule_messages(mid_word, RULE_PROSE_WRAP) == [
        "docstring prose breaks before a word that would still fit within 79 "
        "columns",
    ]
    assert rule_messages(operator, RULE_PROSE_WRAP) == []


def test_fenced_blocks_are_verbatim_content() -> None:
    """
    Leave a fenced pseudocode block entirely outside the fit test.

    A line holding only ''', ```, or ~~~ opens and closes verbatim content
    whose line breaks carry meaning. Joining them would rewrite the pseudocode
    rather than restore greedy wrapping.
    """

    fenced = docstring_source(
        "Summary.",
        "",
        "All of the following are true:",
        "'''",
        "    read.is_paired",
        "and read.is_proper_pair",
        "and read.reference_id == read.next_reference_id",
        "'''",
        "A trailing paragraph.",
    )

    assert rule_messages(fenced, RULE_PROSE_WRAP) == []


def test_deeper_indentation_alone_is_not_wrapped_prose() -> None:
    """
    Keep a non-list line followed by a deeper line outside the fit test.

    Only a list marker licenses the wider continuation column. Without this
    boundary, every entry header followed by its indented description would be
    read as one wrapped paragraph.
    """

    nested = docstring_source(
        "Summary.",
        "",
        "A short paragraph line",
        "        a deeper continuation line.",
    )

    assert rule_messages(nested, RULE_PROSE_WRAP) == []


def test_help_literals_use_greedy_wrapping_and_preserve_value() -> None:
    """
    Require the next whole word to move whenever it still fits.
    """

    accepted = help_source(("x" * 64) + " ", "a continuation.")
    premature = help_source(("x" * 62) + " ", "a continuation.")
    paragraph = help_source("Short paragraph.\\n\\n", "Next paragraph.")
    accepted_line = next(
        line for line in accepted.splitlines() if ("x" * 20) in line
    )

    assert len(accepted_line) == 79
    assert rule_messages(accepted, RULE_CLI_HELP_LAYOUT) == []
    assert rule_messages(paragraph, RULE_CLI_HELP_LAYOUT) == []
    assert rule_messages(premature, RULE_CLI_HELP_LAYOUT) == [
        "parse_args help line breaks before a prose word that would still fit "
        "within 79 columns",
    ]

    before = help_source(
        "Output file path. Path to the output bedGraph file. ",
        "When a real path is provided.",
    )
    after = help_source(
        "Output file path. Path to the output bedGraph file. When a ",
        "real path is provided.",
    )

    assert help_value(before) == help_value(after)


def test_unicode_help_formatter_matches_cli_width_policy() -> None:
    """
    Require Unicode-preserving formatter widths to match emitted literals.
    """

    source = (FIXTURES / "format/help_unicode_input.py").read_text(
        encoding="utf-8",
    )
    expected = (
        FIXTURES / "format/help_unicode_expected.py"
    ).read_text(encoding="utf-8")
    joined = expected.replace('"\n            "available', "available")

    assert format_source(source) == expected
    assert rule_messages(source, RULE_CLI_HELP_LAYOUT) == [
        "parse_args help line breaks before a prose word that would still fit "
        "within 79 columns",
    ]
    assert rule_messages(expected, RULE_CLI_HELP_LAYOUT) == []
    assert rule_messages(joined, RULE_CLI_HELP_LAYOUT) == [
        "splittable static parse_args help source line exceeds 79 columns",
    ]
    assert help_value(source) == help_value(expected)
    assert format_source(expected) == expected


def test_single_static_help_boundaries_and_dynamic_exclusion() -> None:
    """
    Accept fitting static prose, reject overflow, and exclude dynamic help.
    """

    fitting = "x" * 65
    overlong = fitting + "x"
    accepted = help_source(fitting)
    rejected = help_source(overlong)
    accepted_line = next(
        line for line in accepted.splitlines() if fitting in line
    )
    rejected_line = next(
        line for line in rejected.splitlines() if overlong in line
    )

    assert len(accepted_line) < 80
    assert len(rejected_line) == len(accepted_line) + 1
    assert rule_messages(accepted, RULE_CLI_HELP_LAYOUT) == []
    assert rule_messages(rejected, RULE_CLI_HELP_LAYOUT) == [
        "splittable static parse_args help source line exceeds 79 columns",
    ]

    dynamic = help_source("value")
    dynamic = dynamic.replace(
        'help=(\n            "value"\n        )',
        'help=build_help("value", "other")',
    )

    assert rule_messages(dynamic, RULE_CLI_HELP_LAYOUT) == []


def test_static_help_prose_uses_complete_sentences() -> None:
    """
    Require sentence capitalization and terminal punctuation in static help.
    """

    accepted = help_source("Write one output file.")
    case_sensitive = help_source("siQ-ChIP metadata table.")
    rejected = help_source("write one output file")

    assert rule_messages(accepted, RULE_HELP_SENTENCES) == []
    assert rule_messages(case_sensitive, RULE_HELP_SENTENCES) == []
    assert rule_messages(rejected, RULE_HELP_SENTENCES) == [
        "help prose must begin with sentence capitalization",
        "help prose must end with terminal punctuation",
    ]


def test_static_help_sentences_exclude_structured_displays() -> None:
    """
    Exclude literal lists and formulas from prose-sentence enforcement.
    """

    choices = help_source("Choices:\\n  - first\\n  - second")
    formula = help_source("dep_min := max(dep_min, floor)")

    assert rule_messages(choices, RULE_HELP_SENTENCES) == []
    assert rule_messages(formula, RULE_HELP_SENTENCES) == []


def test_docstring_detection_excludes_strings_in_first_calls() -> None:
    """
    Do not mistake a string argument in the first statement for a docstring.
    """

    source = '''
"""
Provide one undocumented override fixture.
"""


class Parser:
    """
    Provide one parser fixture.
    """

    def configure(self) -> None:
        self.values.setdefault("formatter", object)
'''

    assert rule_messages(source, RULE_DOC_LAYOUT) == []


def test_naming_allows_framework_protocols_and_private_types() -> None:
    """
    Accept framework-owned spellings, unused parameters, and private types.
    """

    source = '''
"""
Provide one naming fixture.
"""


class _Reader:
    """
    Provide one private reader.
    """


class Case:
    """
    Provide one test case.
    """

    def setUp(self) -> None:
        self.ready = True

    def tearDown(self) -> None:
        self.ready = False

    def inspect(self, **_: object) -> None:
        return None
'''

    assert rule_messages(source, "PY.NAMING.IDENTIFIERS") == []


def test_naming_rejects_opaque_implementation_shorthand() -> None:
    """
    Reject private shorthand while preserving explicit boundary vocabulary.
    """

    source = '''
"""
Provide one naming fixture.
"""


def render(cfg: dict[str, object], fil_out: str) -> None:
    """
    Render one fixture.
    """

    cvg = cfg
    output_path = fil_out
    _ = cvg, output_path
'''

    messages = rule_messages(source, RULE_NAMING)

    assert messages == [
        "identifier 'cfg' contains opaque shorthand segment(s): cfg",
        "identifier 'cvg' contains opaque shorthand segment(s): cvg",
    ]


def test_topology_counts_blank_lines_before_attached_comments() -> None:
    """
    Keep an attached phase comment with the definition it governs.
    """

    accepted = '''
"""
Provide one topology fixture.
"""

VALUE = 1


# Explain the phase.
def run() -> None:
    """
    Run one fixture.
    """
'''

    assert rule_messages(accepted, RULE_TOPOLOGY) == []


def test_multiline_parameter_lists_and_compact_nested_tuples() -> None:
    """
    Check parameter-list boundaries without expanding compact tuple items.
    """

    accepted = '''
"""
Provide multiline parameter fixtures.
"""

PAIRS = (
    ("left", "right"),
    ("up", "down"),
)


async def run(
    first: int,
    /,
    second: int = 2,
    *values: int,
    third: int,
    **named: int,
) -> None:
    """
    Run one fixture.
    """

    return None


def bare(
    first: int,
    *,
    second: int,
) -> None:
    """
    Run one bare-marker fixture.
    """

    return None
'''

    assert rule_messages(accepted, RULE_MULTILINE) == []

    rejected = '''
"""
Provide one rejected parameter fixture.
"""


def run(
    first: int, second: int
  ) -> None:
    """
    Run one fixture.
    """

    return None
'''

    messages = rule_messages(rejected, RULE_MULTILINE)
    expected = (
        "parameter-list closing delimiter must align with the definition line"
    )

    assert "multiline parameter list requires a trailing comma" in messages
    assert "multiline parameter list requires one item per line" in messages
    assert expected in messages


def test_simple_return_annotations_remain_on_the_definition_line() -> None:
    """
    Reject formatter-unstable parentheses around one simple return type.
    """

    dense = '''
"""
Provide one dense return-annotation fixture.
"""


def result() -> (
    None
):
    return None
'''
    canonical = '''
"""
Provide one canonical return-annotation fixture.
"""


def result() -> None:
    return None
'''

    assert rule_messages(dense, RULE_MULTILINE) == [
        "simple return annotation must remain on the definition line",
    ]
    assert rule_messages(canonical, RULE_MULTILINE) == []


def test_multiline_forms_exclude_sole_generator_calls() -> None:
    """
    Enforce recognized calls without rewriting sole generator arguments.
    """

    source = '''
"""
Provide one multiline fixture.
"""


def values(items: list[int]) -> tuple[bool, list[int]]:
    """
    Return generator and list results.
    """

    matched = all(
        item > 0
        for item in items
    )
    selected = [
        1,
        2,
    ]
    return matched, selected
'''

    assert rule_messages(source, RULE_MULTILINE) == []


def test_multiline_forms_keep_comparisons_with_compact_call_rows() -> None:
    """
    Keep comparisons attached after a bounded compact call payload.
    """

    accepted = '''
"""
Provide one accepted comparison fixture.
"""

assert (
    calculate(
        4.0, 2.0, None, False
    ) == 4.5
)
'''
    rejected = '''
"""
Provide one rejected comparison fixture.
"""

assert (
    calculate(4.0, 2.0, None, False)
    == 4.5
)
'''

    assert rule_messages(accepted, RULE_MULTILINE) == []
    assert rule_messages(rejected, RULE_MULTILINE) == [
        "comparison operator must remain with its left operand",
    ]


def test_compound_conditions_extract_multiline_collection_predicates() -> None:
    """
    Keep a compound condition's boolean relationship visible.
    """

    accepted = '''
"""
Provide one accepted compound-condition fixture.
"""


def classify(value: str, enabled: bool) -> bool:
    """
    Return whether the value is enabled and recognized.
    """

    recognized = value in {"first", "second"}

    if enabled and recognized:
        return True

    return False
'''
    rejected = '''
"""
Provide one rejected compound-condition fixture.
"""


def classify(value: str, enabled: bool) -> bool:
    """
    Return whether the value is enabled and recognized.
    """

    if enabled and value in {
        "first",
        "second",
    }:
        return True

    return False
'''

    assert rule_messages(accepted, RULE_MULTILINE) == []
    assert rule_messages(rejected, RULE_MULTILINE) == [
        "compound condition with a multiline collection requires named "
        "predicates",
    ]


def test_multiline_forms_separate_all_block_literal_boundaries() -> None:
    """
    Separate every block literal from its opening and closing delimiters.
    """

    accepted = '''
source = """
function demo() {
    return 0
}
"""
'''
    rejected = '''
write(
    """function demo() {
        return 0
    }""",
)
'''

    assert rule_messages(accepted, RULE_MULTILINE) == []
    assert rule_messages(rejected, RULE_MULTILINE) == [
        "multiline string closer must be on its own line",
        "multiline string opener must be on its own line",
    ]


def test_multiline_forms_measure_dict_unpacking_from_the_stars() -> None:
    """
    Measure a dictionary-unpacking item from its leading stars.
    """

    source = '''
"""
Provide one dictionary-unpacking fixture.
"""

BASE = {
    "left": 1,
}
MERGED = {
    **BASE,
    "right": 2,
}
'''

    assert rule_messages(source, RULE_MULTILINE) == []


def test_multiline_forms_measure_parenthesized_items_from_the_group() -> None:
    """
    Measure a parenthesized item from its leading grouping delimiter.
    """

    source = '''
"""
Provide one parenthesized-item fixture.
"""

VALUES = [
    ("one"),
    (
        "two-"
        "continued"
    ),
]
'''

    assert rule_messages(source, RULE_MULTILINE) == []


def test_module_imports_precede_executable_calls() -> None:
    """
    Detect late imports even where Ruff deliberately exempts importorskip.
    """

    late = """
import pytest

pytest.importorskip("yaml")

from package import value
"""
    canonical = """
import pytest
from package import value
"""

    assert rule_messages(late, RULE_TOPOLOGY) == [
        "module import must precede executable module statements",
    ]
    assert rule_messages(canonical, RULE_TOPOLOGY) == []


def test_noncompact_control_flow_has_visible_sibling_boundaries() -> None:
    """
    Separate adjacent compound phases and completed compact guards.
    """

    dense = '''
"""
Provide one control-flow-boundary fixture.
"""


def status(remote: bool) -> int:
    value = 0
    if remote:
        value = 1
        print(value)
    if value:
        print(value)
        return value
    return 0
'''
    canonical = '''
"""
Provide one control-flow-boundary fixture.
"""


def status(remote: bool) -> int:
    value = 0

    if remote:
        value = 1
        print(value)

    if value:
        print(value)

        return value

    return 0
'''
    compact_guards = '''
"""
Provide one compact-guard fixture.
"""


def validate(first: bool, second: bool) -> None:
    if first:
        raise ValueError("first")

    if second:
        raise ValueError("second")
'''

    boundary_message = (
        "noncompact control-flow phase requires a blank-line boundary from "
        "its sibling statement"
    )

    assert rule_messages(dense, RULE_TOPOLOGY) == [
        boundary_message,
        boundary_message,
        "terminal transfer after a substantive phase requires a blank-line "
        "boundary",
        "terminal transfer after a substantive phase requires a blank-line "
        "boundary",
    ]
    assert rule_messages(canonical, RULE_TOPOLOGY) == []
    assert rule_messages(compact_guards, RULE_TOPOLOGY) == []


def test_compact_transfer_guards_end_with_visible_boundaries() -> None:
    """
    Separate each completed guard from its following sibling statement.
    """

    dense = '''
"""
Provide compact transfer-guard fixtures.
"""


def select(value: int | None) -> int:
    if value is None:
        return 0
    return value


def collect(values: list[int | None]) -> list[int]:
    selected = []

    for value in values:
        if value is None:
            continue
        selected.append(value)

    return selected
'''
    canonical = '''
"""
Provide compact transfer-guard fixtures.
"""


def select(value: int | None) -> int:
    if value is None:
        return 0

    return value


def collect(values: list[int | None]) -> list[int]:
    selected = []

    for value in values:
        if value is None:
            continue

        selected.append(value)

    return selected
'''
    expected = (
        "completed compact transfer guard requires a blank-line boundary from "
        "the following statement"
    )

    assert rule_messages(dense, RULE_TOPOLOGY) == [expected, expected]
    assert rule_messages(canonical, RULE_TOPOLOGY) == []


def test_one_action_branches_are_not_compact_transfer_guards() -> None:
    """
    Separate ordinary guarded work while retaining direct transfers.
    """

    dense = '''
"""
Provide one guarded-action fixture.
"""


def collect(enabled: bool, values: list[int]) -> int:
    current = 1
    if enabled:
        values.append(current)
    current += 1
    return current
'''
    canonical = '''
"""
Provide one canonical guarded-action fixture.
"""


def collect(enabled: bool, values: list[int]) -> int:
    current = 1

    if enabled:
        values.append(current)

    current += 1

    return current
'''
    boundary_message = (
        "noncompact control-flow phase requires a blank-line boundary from "
        "its sibling statement"
    )

    assert rule_messages(dense, RULE_TOPOLOGY) == [
        boundary_message,
        boundary_message,
        "terminal transfer after a substantive phase requires a blank-line "
        "boundary",
    ]
    assert rule_messages(canonical, RULE_TOPOLOGY) == []


def test_noncompact_control_flow_is_separated_from_plain_siblings() -> None:
    """
    Enforce both visible edges of one noncompact control-flow phase.
    """

    dense = '''
"""
Provide one mixed sibling-boundary fixture.
"""


def inspect(enabled: bool) -> int:
    value = 0
    if enabled:
        value = 1
        print(value)
    result = value + 1
    return result
'''
    canonical = '''
"""
Provide one canonical mixed sibling-boundary fixture.
"""


def inspect(enabled: bool) -> int:
    value = 0

    if enabled:
        value = 1
        print(value)

    result = value + 1

    return result
'''
    boundary_message = (
        "noncompact control-flow phase requires a blank-line boundary from "
        "its sibling statement"
    )

    assert rule_messages(dense, RULE_TOPOLOGY) == [
        boundary_message,
        boundary_message,
        "terminal transfer after a substantive phase requires a blank-line "
        "boundary",
    ]
    assert rule_messages(canonical, RULE_TOPOLOGY) == []


def test_substantive_phases_are_separated_from_every_transfer_form() -> None:
    """
    Enforce transfer boundaries while retaining one-step result preparation.
    """

    dense = '''
"""
Provide one terminal-transfer fixture.
"""


def return_value() -> int:
    first = 1
    second = 2
    return first + second


def raise_error() -> None:
    message = "bad"
    print(message)
    raise ValueError(message)


def yield_value() -> object:
    first = 1
    second = 2
    yield first + second


def stop(values: list[int]) -> None:
    for value in values:
        print(value)
        value += 1
        continue


def leave() -> None:
    first = 1
    second = 2
    sys.exit(first + second)
'''
    canonical = (
        dense.replace(
            "    return first + second",
            "\n    return first + second",
        )
        .replace(
            "    raise ValueError(message)",
            "\n    raise ValueError(message)",
        )
        .replace(
            "    yield first + second",
            "\n    yield first + second",
        )
        .replace(
            "        continue",
            "\n        continue",
        )
        .replace(
            "    sys.exit(first + second)",
            "\n    sys.exit(first + second)",
        )
    )
    compact = '''
"""
Provide one compact terminal-transfer fixture.
"""


def build() -> int:
    value = compute()
    return value


def validate(invalid: bool) -> None:
    if invalid:
        raise ValueError("invalid")
'''

    expected_message = (
        "terminal transfer after a substantive phase requires a blank-line "
        "boundary"
    )

    assert rule_messages(dense, RULE_TOPOLOGY) == [expected_message] * 5
    assert rule_messages(canonical, RULE_TOPOLOGY) == []
    assert rule_messages(compact, RULE_TOPOLOGY) == []


def test_setup_action_and_assertions_have_visible_phases() -> None:
    """
    Separate direct test assertions from setup and later mutations.
    """

    dense = '''
"""
Provide one dense test-phase fixture.
"""


def test_value() -> None:
    value = 1
    assert value == 1
    value += 1
    assert value == 2
'''
    canonical = '''
"""
Provide one canonical test-phase fixture.
"""


def test_value() -> None:
    value = 1

    assert value == 1

    value += 1

    assert value == 2
'''

    assert rule_messages(dense, RULE_TOPOLOGY) == [
        "semantic setup, action, validation, and assertion phases require a "
        "blank-line boundary",
        "semantic setup, action, validation, and assertion phases require a "
        "blank-line boundary",
        "semantic setup, action, validation, and assertion phases require a "
        "blank-line boundary",
    ]
    assert rule_messages(canonical, RULE_TOPOLOGY) == []


def test_actions_large_constructions_and_loops_have_visible_phases() -> None:
    """
    Separate fixture setup, baseline actions, case matrices, and loops.
    """

    dense = '''
"""
Provide one dense test-phase fixture.
"""


def test_cases() -> None:
    fixture = load_fixture()
    expected = expected_result()
    validate_fixture(fixture)
    cases = {
        "one": (
            1,
            2,
            3,
            4,
            5,
        ),
    }
    for label, values in cases.items():
        validate_case(label, values)
'''
    canonical = '''
"""
Provide one canonical test-phase fixture.
"""


def test_cases() -> None:
    fixture = load_fixture()
    expected = expected_result()

    validate_fixture(fixture)

    cases = {
        "one": (
            1,
            2,
            3,
            4,
            5,
        ),
    }

    for label, values in cases.items():
        validate_case(label, values)
'''

    expected = (
        "semantic setup, action, validation, and assertion phases require a "
        "blank-line boundary"
    )

    assert rule_messages(dense, RULE_TOPOLOGY) == [
        expected,
        expected,
        expected,
    ]
    assert rule_messages(canonical, RULE_TOPOLOGY) == []


def test_semantic_phase_boundaries_apply_inside_nested_suites() -> None:
    """
    Keep arrange, act, and assert phases visible inside test-owned suites.
    """

    dense = '''
"""
Provide one dense nested test fixture.
"""


def test_nested() -> None:
    with fixture():
        value = observe()
        assert value
'''
    canonical = '''
"""
Provide one canonical nested test fixture.
"""


def test_nested() -> None:
    with fixture():
        value = observe()

        assert value
'''

    expected = (
        "semantic setup, action, validation, and assertion phases require a "
        "blank-line boundary"
    )

    assert rule_messages(dense, RULE_TOPOLOGY) == [expected]
    assert rule_messages(canonical, RULE_TOPOLOGY) == []


def test_semantic_phase_boundaries_apply_outside_tests() -> None:
    """
    Keep setup, validation, and later work visible in production helpers.
    """

    dense = '''
"""
Provide one dense production fixture.
"""


def load_and_validate() -> None:
    first = load_first()
    second = load_second()
    validate_pair(first, second)
    persist_pair(first, second)
'''
    canonical = '''
"""
Provide one canonical production fixture.
"""


def load_and_validate() -> None:
    first = load_first()
    second = load_second()

    validate_pair(first, second)

    persist_pair(first, second)
'''

    expected = (
        "semantic setup, action, validation, and assertion phases require a "
        "blank-line boundary"
    )

    assert rule_messages(dense, RULE_TOPOLOGY) == [expected, expected]
    assert rule_messages(canonical, RULE_TOPOLOGY) == []


def test_long_repeated_validation_calls_remain_one_phase() -> None:
    """
    Keep consecutive validations together despite multiline call spans.
    """

    source = '''
"""
Provide one validation-run fixture.
"""


def validate_values(first: int, second: int) -> None:
    validate_comparison(
        first,
        "gt",
        0,
        "first",
        allow_none=False,
    )
    validate_comparison(
        second,
        "gt",
        0,
        "second",
        allow_none=False,
    )
'''

    assert rule_messages(source, RULE_TOPOLOGY) == []


def test_result_index_assignment_has_visible_return_boundary() -> None:
    """
    Separate a computed result index from its terminal return.
    """

    dense = '''
"""
Provide one dense result fixture.
"""


def selected_values(values: list[int]) -> list[int]:
    selected_indices = find_selected_indices(values)
    return [values[index] for index in selected_indices]
'''
    canonical = '''
"""
Provide one canonical result fixture.
"""


def selected_values(values: list[int]) -> list[int]:
    selected_indices = find_selected_indices(values)

    return [values[index] for index in selected_indices]
'''

    expected = (
        "semantic setup, action, validation, and assertion phases require a "
        "blank-line boundary"
    )

    assert rule_messages(dense, RULE_TOPOLOGY) == [expected]
    assert rule_messages(canonical, RULE_TOPOLOGY) == []


def test_operation_result_assignments_define_visible_phases() -> None:
    """
    Separate arrange, result capture, mutation, and verification topics.
    """

    dense = '''
"""
Provide one dense operation-result fixture.
"""


def test_report() -> None:
    root = fixture_root()
    target = root / "input.txt"
    target.write_text("input")
    report = scan_target(target)
    report["status"] = "reviewed"
    findings = validate_report(report)
    assert not findings


def production_report() -> dict[str, object]:
    root = fixture_root()
    target = root / "input.txt"
    target.write_text("input")
    report = scan_target(target)
    normalize_report(report)
    return report
'''
    canonical = '''
"""
Provide one canonical operation-result fixture.
"""


def test_report() -> None:
    root = fixture_root()
    target = root / "input.txt"
    target.write_text("input")

    report = scan_target(target)

    report["status"] = "reviewed"

    findings = validate_report(report)

    assert not findings


def production_report() -> dict[str, object]:
    root = fixture_root()
    target = root / "input.txt"
    target.write_text("input")

    report = scan_target(target)

    normalize_report(report)

    return report
'''

    expected = (
        "semantic setup, action, validation, and assertion phases require a "
        "blank-line boundary"
    )

    transfer = (
        "terminal transfer after a substantive phase requires a blank-line "
        "boundary"
    )

    assert rule_messages(dense, RULE_TOPOLOGY) == [
        expected,
        expected,
        expected,
        expected,
        expected,
        expected,
        transfer,
    ]
    assert rule_messages(canonical, RULE_TOPOLOGY) == []


def test_checker_is_idempotent_and_rule_selector_rejects_near_misses() -> None:
    """
    Preserve stable reports and reject unknown rule selectors.
    """

    text = (FIXTURES / "accepted/canonical.py").read_text(encoding="utf-8")

    assert analyze_text(text, "accepted/canonical.py") == analyze_text(
        text,
        "accepted/canonical.py",
    )

    with pytest.raises(SystemExit):
        parse_args(["--rule", "PY.DOCSTRING.LAYOU"])
