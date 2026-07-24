#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_markdown_policy.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


"""Unit tests for pure Markdown checker and formatter primitives."""

from __future__ import annotations

from pathlib import Path

import pytest
from dev.audit.markdown_policy import (
    check_explicit_links,
    check_standard_sections,
    check_text,
)
from dev.tools.markdown_format import (
    convert_delimited,
    format_deterministic,
    format_document,
)
from dev.tools.markdown_format import (
    main as formatter_main,
)

ROOT = Path(__file__).resolve().parents[3]
FIXTURES = ROOT / "tests/fixtures/markdown"


def read(relative: str) -> str:
    """Read one Markdown fixture."""

    return (FIXTURES / relative).read_text(encoding="utf-8")


@pytest.mark.parametrize(
    "fixture",
    [
        "accepted/basic.md",
        "accepted/anchors.md",
        "accepted/section_boundaries.md",
    ],
)
def test_accepted_fixture_has_no_deterministic_findings(
    fixture: str,
) -> None:
    findings = check_text(read(fixture))
    assert not [
        item for item in findings
        if item.classification == "deterministic"
    ]


def section_count(text: str) -> int:
    """Return deterministic section-boundary finding count."""

    return sum(
        item.rule_id == "MD.SECTION.BREAK"
        for item in check_text(text)
    )


def section_signature(text: str) -> tuple[tuple[str, str], ...]:
    """Return line-number-independent section decisions."""

    return tuple(
        (item.message, item.classification)
        for item in check_text(text)
        if item.rule_id == "MD.SECTION.BREAK"
    )


def anchored(heading: str, identifier: str, present: bool) -> str:
    """Return one heading with an optional canonical explicit anchor."""

    if not present:
        return heading
    return f'<a id="{identifier}"></a>\n{heading}'


@pytest.mark.parametrize("with_anchor", [False, True])
def test_ordinary_boundary_anchor_transparency(
    with_anchor: bool,
) -> None:
    valid = (
        "\n# Title\nContent.\n\n<br />\n\n"
        + anchored("## Later", "later", with_anchor)
        + "\nBody.\n"
    )
    invalid = (
        "\n# Title\nContent.\n"
        + anchored("## Later", "later", with_anchor)
        + "\nBody.\n"
    )
    assert section_count(valid) == 0
    assert section_count(invalid) == 1


@pytest.mark.parametrize("parent_anchor", [False, True])
@pytest.mark.parametrize("child_anchor", [False, True])
def test_contentless_parent_child_anchor_matrix(
    parent_anchor: bool,
    child_anchor: bool,
) -> None:
    valid = (
        "\n# Title\n"
        + anchored("## Parent", "parent", parent_anchor)
        + "\n"
        + anchored("### Child", "child", child_anchor)
        + "\nBody.\n"
    )
    invalid = valid.replace(
        anchored("### Child", "child", child_anchor),
        "\n<br />\n\n"
        + anchored("### Child", "child", child_anchor),
    )
    assert section_count(valid) == 0
    assert section_count(invalid) == 1


@pytest.mark.parametrize("with_anchor", [False, True])
def test_parent_content_requires_ordinary_boundary(
    with_anchor: bool,
) -> None:
    heading = anchored("### Child", "child", with_anchor)
    valid = (
        "\n# Title\n## Parent\nContent.\n\n<br />\n\n"
        f"{heading}\nBody.\n"
    )
    invalid = (
        "\n# Title\n## Parent\nContent.\n"
        f"{heading}\nBody.\n"
    )
    assert section_count(valid) == 0
    assert section_count(invalid) == 1


@pytest.mark.parametrize(
    ("later", "identifier"),
    [("## Sibling", "sibling"), ("# Shallower", "shallower")],
)
@pytest.mark.parametrize("with_anchor", [False, True])
def test_empty_sibling_and_shallower_use_ordinary_boundary(
    later: str,
    identifier: str,
    with_anchor: bool,
) -> None:
    heading = anchored(later, identifier, with_anchor)
    valid = (
        "\n# Title\n## Empty\n\n<br />\n\n"
        f"{heading}\nBody.\n"
    )
    invalid = f"\n# Title\n## Empty\n{heading}\nBody.\n"
    assert section_count(valid) == 0
    assert section_count(invalid) == 1


@pytest.mark.parametrize("with_anchor", [False, True])
def test_details_close_boundary_anchor_transparency(
    with_anchor: bool,
) -> None:
    heading = anchored("## Later", "later", with_anchor)
    prefix = (
        "\n# Title\n<details>\n<summary>More</summary>\n"
        "Content.\n</details>\n<br />\n\n"
    )
    valid = f"{prefix}{heading}\nBody.\n"
    invalid = f"{prefix}<br />\n\n{heading}\nBody.\n"
    assert section_count(valid) == 0
    assert section_count(invalid) == 1


@pytest.mark.parametrize("with_anchor", [False, True])
def test_doubled_ordinary_break_is_prohibited(
    with_anchor: bool,
) -> None:
    heading = anchored("## Later", "later", with_anchor)
    source = (
        "\n# Title\nContent.\n\n<br />\n\n<br />\n\n"
        f"{heading}\nBody.\n"
    )
    assert section_count(source) == 1


@pytest.mark.parametrize(
    "template",
    [
        "\n# Title\nContent.\n\n<br />\n\n{unit}\nBody.\n",
        "\n# Title\nContent.\n{unit}\nBody.\n",
        (
            "\n# Title\n## Parent\nContent.\n\n<br />\n\n"
            "{unit}\nBody.\n"
        ),
        (
            "\n# Title\n## Empty\n\n<br />\n\n"
            "{unit}\nBody.\n"
        ),
        (
            "\n# Title\n<details>\n<summary>More</summary>\n"
            "Content.\n</details>\n<br />\n\n{unit}\nBody.\n"
        ),
        (
            "\n# Title\nContent.\n\n<br />\n\n<br />\n\n"
            "{unit}\nBody.\n"
        ),
    ],
)
def test_anchor_presence_keeps_boundary_decisions_identical(
    template: str,
) -> None:
    without = template.format(unit="## Later")
    with_anchor = template.format(
        unit='<a id="later"></a>\n## Later'
    )
    assert section_signature(without) == section_signature(with_anchor)


def test_parent_and_child_anchor_matrix_has_identical_decisions() -> None:
    signatures = []
    for parent_anchor in (False, True):
        for child_anchor in (False, True):
            source = (
                "\n# Title\n"
                + anchored("## Parent", "parent", parent_anchor)
                + "\n"
                + anchored("### Child", "child", child_anchor)
                + "\nBody.\n"
            )
            signatures.append(section_signature(source))
    assert signatures == [signatures[0]] * 4


def test_noncanonical_adjacent_anchor_has_independent_diagnostic() -> None:
    without = "\n# Title\nContent.\n\n<br />\n\n## Later\nBody.\n"
    with_invalid = without.replace(
        "## Later",
        '<a name="later"></a>\n## Later',
    )
    assert section_count(without) == section_count(with_invalid) == 0
    assert "MD.ANCHOR.CANONICAL" in {
        item.rule_id for item in check_text(with_invalid)
    }


def test_informal_heading_recognition_does_not_require_conformance() -> None:
    source = (
        "\n# Title\n\n<br />\n\n###### Parent\n\n<br />\n\n"
        "**Child**\nBody.\n"
    )
    findings = check_text(source)
    assert section_count(source) == 1
    assert not [
        item for item in findings
        if item.rule_id == "MD.HEADING.INFORMAL"
    ]


def test_ambiguous_emphasis_remains_visible_advisory() -> None:
    source = "\n# Title\nContent.\n*Possible H8 or prose.*\n"
    findings = check_text(source)
    assert section_count(source) == 0
    assert [
        item for item in findings
        if item.rule_id == "MD.HEADING.INFORMAL"
        and item.classification == "advisory"
    ]


def test_section_fixture_covers_accepted_and_rejected_forms() -> None:
    assert section_count(read("accepted/section_boundaries.md")) == 0
    assert section_count(read("rejected/section_boundaries.md")) == 6


def test_anchor_fixture_reports_independent_defects() -> None:
    findings = [
        item
        for item in check_text(read("rejected/anchors.md"))
        if item.rule_id == "MD.ANCHOR.CANONICAL"
    ]
    messages = {item.message for item in findings}
    assert 'heading anchor must use exact <a id="ID"></a> source' in messages
    assert (
        "heading anchor must be immediately followed by one heading"
        in messages
    )
    assert "consecutive heading anchors are prohibited" in messages
    assert any("duplicated" in message for message in messages)


def test_explicit_links_are_narrow_and_exact(tmp_path: Path) -> None:
    target = tmp_path / "target.md"
    source = tmp_path / "source.md"
    texts = {
        target: (
            "\n# Target\n<a id=\"exact\"></a>\n"
            "## First\n\n<br />\n\n<a id=\"exact\"></a>\n## Second\n"
        ),
        source: (
            "\n# Source\n[Explicit](target.md#exact) and "
            "[generated](target.md#renderer-fragment).\n"
        ),
    }
    findings = check_explicit_links(texts)
    assert len(findings[source.resolve()]) == 1
    assert "explicit-ID link is ambiguous" in (
        findings[source.resolve()][0].message
    )


def test_fenced_anchor_and_boundary_examples_are_literal() -> None:
    source = (
        "\n# Title\n```markdown\n<a id=\"same\"></a>\n"
        "<a id=\"same\"></a>\n## Literal\n```\n"
    )
    owned = {
        "MD.ANCHOR.CANONICAL",
        "MD.SECTION.BREAK",
    }
    assert not [
        item for item in check_text(source)
        if item.rule_id in owned
    ]


@pytest.mark.parametrize(
    ("leading", "has_finding"),
    [("", True), ("\n", False), ("\n\n", True), ("\n\n\n", True)],
)
def test_file_boundary_complete_leading_line_feed_matrix(
    leading: str,
    has_finding: bool,
) -> None:
    """Cover zero, one, and two-or-more leading line feeds directly."""

    findings = [
        item
        for item in check_text(f"{leading}# Title\nBody.\n")
        if item.rule_id == "MD.FILE.BOUNDARY" and item.line == 1
    ]
    assert bool(findings) is has_finding


@pytest.mark.parametrize(
    ("terminal", "has_finding"),
    [("", True), ("\n", False), ("\n\n", True), ("\n\n\n", True)],
)
def test_file_boundary_complete_terminal_line_feed_matrix(
    terminal: str,
    has_finding: bool,
) -> None:
    """Cover zero, one, and two-or-more terminal line feeds directly."""

    findings = [
        item
        for item in check_text(f"\n# Title\nBody.{terminal}")
        if item.rule_id == "MD.FILE.BOUNDARY" and item.line != 1
    ]
    assert bool(findings) is has_finding


def test_colon_structure_finding_has_exact_owned_identity() -> None:
    """Keep constructor order and the deterministic diagnostic stable."""

    findings = [
        item
        for item in check_text("\n# Title\nIntro:\n\n- item\n")
        if item.rule_id == "MD.COLON.STRUCTURE"
    ]
    assert findings == [
        type(findings[0])(
            rule_id="MD.COLON.STRUCTURE",
            line=3,
            message=(
                "colon needs one blank before a table and none before a list or fence"
            ),
            classification="deterministic",
        )
    ]


@pytest.mark.parametrize(
    ("fixture", "rule_id"),
    [
        ("rejected/unclosed_fence.md", "MD.FENCE.CLOSED"),
        ("rejected/spacing.md", "MD.FILE.BOUNDARY"),
        ("rejected/spacing.md", "MD.HEADING.SPACING"),
        ("rejected/spacing.md", "MD.COLON.STRUCTURE"),
        ("rejected/spacing.md", "MD.LIST.SPACING"),
        ("rejected/anchors.md", "MD.ANCHOR.CANONICAL"),
        ("rejected/section_boundaries.md", "MD.SECTION.BREAK"),
    ],
)
def test_rejected_fixtures_report_owned_rules(
    fixture: str,
    rule_id: str,
) -> None:
    assert rule_id in {item.rule_id for item in check_text(read(fixture))}


def test_formatter_matches_expected_and_is_idempotent() -> None:
    expected = read("format/expected.md")
    actual = format_document(read("format/input.md"))
    assert actual == expected
    assert format_document(actual) == actual


def test_deterministic_formatter_preserves_deferred_behaviors() -> None:
    source = (
        "\n# Title\n\nHard-wrapped prose that must\nremain wrapped.\n\n"
        "<details>\n<summary>More</summary>\n### Child\n</details>\n"
        "\nValues:\n| A | B |\n| --- | ---: |\n| x | 界 |\n"
    )
    expected = (
        "\n# Title\nHard-wrapped prose that must\nremain wrapped.\n\n"
        "<details>\n<summary>More</summary>\n\n<br />\n\n"
        "### Child\n</details>\n"
        "\nValues:\n\n| A    | B    |\n| :--- | ---: |\n"
        "| x    | 界   |\n"
    )
    actual = format_deterministic(source)
    assert actual == expected
    assert format_deterministic(actual) == actual


def test_deterministic_formatter_preserves_fenced_literals() -> None:
    source = "\n# Title\n```markdown\n## Literal\nA | B\n--- | ---\n```\n"
    assert format_deterministic(source) == source


def test_deterministic_formatter_handles_anchor_boundary_cluster() -> None:
    source = (
        "\n# Title\nBody.\n\n<a id=\"later\"></a>\n\n<br />\n\n"
        "<a id=\"later\"></a>\n## Later\nBody.\n"
    )
    expected = (
        "\n# Title\nBody.\n\n<br />\n\n"
        "<a id=\"later\"></a>\n## Later\nBody.\n"
    )
    actual = format_deterministic(source)
    assert actual == expected
    assert format_deterministic(actual) == actual


def test_deterministic_formatter_does_not_guess_distinct_ids() -> None:
    source = (
        "\n# Title\nBody.\n\n<a id=\"first\"></a>\n\n<br />\n\n"
        "<a id=\"second\"></a>\n## Later\nBody.\n"
    )
    assert format_deterministic(source) == source


def test_proposed_mode_is_preview_only(tmp_path: Path) -> None:
    path = tmp_path / "sample.md"
    path.write_text("\n# Title\nHard\nwrapped.\n", encoding="utf-8")
    with pytest.raises(SystemExit, match="2"):
        formatter_main(["--mode", "proposed", "--write", str(path)])
    assert path.read_text(encoding="utf-8") == "\n# Title\nHard\nwrapped.\n"


def test_delimited_conversion_requires_explicit_equal_width_rows() -> None:
    output = convert_delimited("Name\tValue\nA\t界\n", "\t")
    assert "| Name | Value |" in output
    with pytest.raises(ValueError, match="equal column counts"):
        convert_delimited("a,b\n1\n", ",")


def test_standard_rule_sections_require_canonical_fields_and_order() -> None:
    """Require all canonical standards fields in their source order."""
    assert not check_standard_sections(
        read("accepted/standard_section.md")
    )
    findings = check_standard_sections(
        read("rejected/standard_section.md")
    )
    messages = {item.message for item in findings}
    assert "rule fields must appear in canonical order" in messages
    assert (
        "**Classification:** must immediately follow the heading" in messages
    )
    assert "rule section is missing **Scope:**" in messages
    assert {item.rule_id for item in findings} == {"MD.STANDARD.SECTION"}


def test_standard_rule_sections_reject_idless_h2() -> None:
    """Do not let an ID-less maintained H2 disappear from review."""
    findings = check_standard_sections(
        "\n# Standard\nIntro.\n\n<br />\n\n## Hidden policy\nText.\n"
    )
    assert [item.message for item in findings] == [
        "maintained standards H2 is missing a rule ID"
    ]
