#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: help_examples.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


"""
Enforce repository-wide strict Examples structure from current source.
"""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
import re
import shlex
import sys
from collections.abc import Iterable
from pathlib import Path

from dev.audit.help_aliases import (
    alias_chunks,
    file_alias_chunks,
    function_bodies,
    hidden_alias,
    parameter_rows,
)
from dev.audit.help_aliases import (
    scan_repository as scan_alias_repository,
)
from dev.audit.help_heredoc_reflow import Heredoc, extract_help_heredocs
from dev.audit.help_style import Section, sections

RULE_REQUIRED = "HELP.EXAMPLES.REQUIRED"
RULE_COUNT = "HELP.EXAMPLES.COUNT"
RULE_ENTRY = "HELP.EXAMPLES.ENTRY"
RULE_NUMBERING = "HELP.EXAMPLES.NUMBERING"
RULE_FINAL = "HELP.EXAMPLES.FINAL"
RULE_SECTION_COUNT = "HELP.EXAMPLES.SECTION_COUNT"
RULE_CODE_BLOCK = "HELP.EXAMPLES.CODE_BLOCK"
RULE_OWNER = "HELP.EXAMPLES.OWNER_INVOCATION"
RULE_ALIAS_VISIBILITY = "HELP.EXAMPLES.ALIAS_VISIBILITY"
RULE_ALIAS_ACCEPTANCE = "HELP.EXAMPLES.ALIAS_ACCEPTANCE"
RULE_COMPLETE = "HELP.EXAMPLES.STRUCTURAL_COMPLETE"
RULE_DUPLICATE = "HELP.EXAMPLES.DUPLICATE"
RULE_SIGNATURE_DUPLICATE = "HELP.EXAMPLES.SIGNATURE_DUPLICATE"
RULE_REVIEW = "HELP.EXAMPLES.REVIEW.MATERIAL_DISTINCTNESS"
RULE_REVIEW_UNSAFE = "HELP.EXAMPLES.REVIEW.UNSAFE"
RULE_REVIEW_INDIRECT = "HELP.EXAMPLES.REVIEW.INDIRECT_OWNER"
RULE_OWNERSHIP = "HELP.EXAMPLES.OWNERSHIP"
RULE_REVIEW_PENDING = "HELP.EXAMPLES.REVIEW.UNDISPOSITIONED"

ENTRY = re.compile(r"^  (?P<number>\d+)\. (?P<description>\S.*)$")
OPTION = re.compile(r"(?<![A-Za-z0-9_-])--?[A-Za-z][A-Za-z0-9_-]*")
DELEGATE = re.compile(
    r"exec\s+\"?\$\{?BASH\}?\"?\s+"
    r"\"?[^\n]*?/(?P<target>(?:execute|submit)_[A-Za-z0-9_]+\.sh)\"?\s+"
    r"\"?\$@\"?",
)
SELECTOR_OPTIONS = {
    "--aligner",
    "--aln_typ",
    "--bt2_mode",
    "--bwa_alg",
    "--engine",
    "--eqn",
    "--layout",
    "--method",
    "--mode",
    "--out_ext",
    "--retain",
    "--skip_00",
    "--typ_out",
}


@dataclasses.dataclass(frozen=True)
class Finding:
    """
    One strict repository-wide Examples finding.
    """

    rule_id: str
    path: str
    line: int
    owner: str
    message: str

    def format(self) -> str:
        """
        Render one stable line-oriented diagnostic.
        """

        return (
            f"{self.rule_id}: {self.path}:{self.line}: "
            f"owner={self.owner}; {self.message}"
        )


@dataclasses.dataclass(frozen=True)
class Review:
    """
    One non-mechanical semantic-review candidate.
    """

    rule_id: str
    path: str
    line: int
    owner: str
    examples: tuple[int, ...]
    message: str

    def as_dict(self) -> dict[str, object]:
        """
        Return a deterministic JSON-ready record.
        """

        value = dataclasses.asdict(self)
        payload = json.dumps(value, sort_keys=True, separators=(",", ":"))

        value["signature"] = (
            "sha256:" + hashlib.sha256(payload.encode()).hexdigest()
        )

        return value


@dataclasses.dataclass(frozen=True)
class Example:
    """
    One numbered example and its deterministic invocation facts.
    """

    number: int
    line: int
    description: str
    code: tuple[str, ...]
    normalized_invocation: str
    signature: tuple[object, ...]
    option_names: tuple[str, ...]
    positional_shape: tuple[str, ...]
    setup_global_state_shape: tuple[str, ...]
    temporary_resource_shape: tuple[str, ...]
    control_flow_shape: tuple[str, ...]
    expected_outcome_shape: tuple[str, ...]
    skip_failure_shape: tuple[str, ...]
    cleanup_shape: tuple[str, ...]

    def as_dict(self) -> dict[str, object]:
        """
        Return a deterministic JSON-ready record.
        """

        value = dataclasses.asdict(self)
        value["signature"] = list(self.signature)

        return value


@dataclasses.dataclass(frozen=True)
class Analysis:
    """
    Strict findings, semantic candidates, and parsed examples for one owner.
    """

    findings: list[Finding]
    reviews: list[Review]
    examples: list[Example]
    status: str


@dataclasses.dataclass(frozen=True)
class Ownership:
    """
    One execute/submit wrapper full-help classification.
    """

    path: str
    classification: str
    full_help_owner: str
    surface: str
    target: str = ""

    def as_dict(self) -> dict[str, str]:
        """
        Return a deterministic JSON-ready record.
        """

        return dataclasses.asdict(self)


@dataclasses.dataclass(frozen=True)
class RepositoryResult:
    """
    Complete wrapper/function Examples audit result.
    """

    findings: list[Finding]
    reviews: list[Review]
    inventory: list[dict[str, object]]
    ownership: list[Ownership]
    alias_findings: list[object]


InvocationFacts = tuple[
    str,
    tuple[object, ...],
    tuple[str, ...],
    tuple[str, ...],
    tuple[str, ...],
]


def rendered_heredoc(text: str, owner: str) -> Heredoc:
    """
    Represent one rendered sectioned document as a bounded heredoc.
    """

    lines = tuple(enumerate(text.splitlines(), 1))
    return Heredoc(
        owner=owner,
        delimiter="<rendered>",
        start_line=1,
        end_line=max(1, len(lines)),
        lines=lines,
    )


def normalized_code(lines: Iterable[str]) -> str:
    """
    Remove comments/help continuations and collapse cosmetic whitespace.
    """

    cleaned = [
        re.sub(
            r"\\+[ \t]*$",
            "",
            re.sub(r"\s+#.*$", "", line),
        ).strip()
        for line in lines
    ]

    return " ".join(" ".join(cleaned).split())


def positional_kind(value: str) -> str:
    """
    Classify one bounded positional token without preserving placeholders.
    """

    stripped = value.strip("\"'")
    if re.fullmatch(r"[-+]?\d+(?:\.\d+)?", stripped):
        return "number"

    if stripped.lower() in {"true", "false"}:
        return "bool"

    if stripped.startswith(("${", "$")):
        return "expansion"

    if "/" in stripped or re.search(
        r"\.(?:bam|cram|bed|bdg|gz|fq|fastq|txt)$",
        stripped,
    ):
        return "path"

    return "value"


def option_value_shape(value: str) -> str:
    """
    Classify notable list and file-format shape for one option value.
    """

    stripped = value.strip("\"'")
    delimiters = ""

    if ";" in stripped:
        delimiters += "semicolon_"

    if "," in stripped:
        delimiters += "comma_"

    extension = re.search(
        r"\.(bam|cram|fastq|fq|bdg|bed)(?:\.gz)?(?:[,;]|$)",
        stripped,
    )
    suffix = (
        f"{extension.group(1)}" if extension else positional_kind(stripped)
    )

    return f"{delimiters}{suffix}"


def invocation_facts(
    code: tuple[str, ...],
    owner: str,
) -> InvocationFacts:
    """
    Derive the deterministic bounded invocation signature.

    Parameters
    ----------
    code : tuple[str, ...]
        Shell source lines comprising one recognized invocation.
    owner : str
        Callable or command identity invoked by the source.

    Returns
    -------
    facts : InvocationFacts
        Normalized source, comparison signature, public options, positional
        shapes, and expected outcome shapes.
    """

    semantic_code = tuple(
        line for line in code if not line.lstrip().startswith("#")
    )
    normalized = normalized_code(semantic_code)

    try:
        tokens = shlex.split(normalized, comments=True, posix=True)
    except ValueError:
        tokens = normalized.split()

    owner_index = next(
        (
            index
            for index, token in enumerate(tokens)
            if Path(token.strip("\"'")).name == owner
            or token.strip("\"'") == owner
        ),
        -1,
    )
    tail = tokens[owner_index + 1 :] if owner_index >= 0 else []
    options: list[str] = []
    selectors: list[tuple[str, str]] = []
    value_shapes: list[tuple[str, str]] = []
    positionals: list[str] = []
    positional_selectors: list[str] = []
    index = 0

    while index < len(tail):
        token = tail[index]

        if token in {"|", "||", "&&", ";", "then", "fi"} or token.startswith(
            (">", "<"),
        ):
            index += 1

            continue

        if OPTION.fullmatch(token):
            options.append(token)

            if (
                token in SELECTOR_OPTIONS
                and index + 1 < len(tail)
                and not tail[index + 1].startswith("-")
            ):
                selectors.append((token, tail[index + 1].strip("\"'")))

            if index + 1 < len(tail) and not tail[index + 1].startswith("-"):
                value_shapes.append(
                    (token, option_value_shape(tail[index + 1])),
                )
                index += 2
            else:
                index += 1

            continue

        positionals.append(positional_kind(token))
        literal = token.strip("\"'")

        if literal.upper() in {"NA", "SE", "PE", "BAM", "CRAM", "F", "D"}:
            positional_selectors.append(literal.upper())

        index += 1

    operators = tuple(
        operator
        for operator in ("|", "||", "&&", ">", ">>", "<")
        if operator in normalized
    )

    owner_line = next(
        (index for index, line in enumerate(semantic_code) if owner in line),
        0,
    )
    setup_lines = semantic_code[:owner_line]
    setup_shape = tuple(
        shape
        for shape, pattern in (
            ("source", r"\bsource\b|^\s*\.\s"),
            ("assignment", r"^\s*[A-Za-z_][A-Za-z0-9_]*="),
            ("array_declaration", r"\b(?:declare|local)\s+-[A-Za-z]*a\b"),
            (
                "function_definition",
                r"^\s*(?:function\s+)?[A-Za-z_][A-Za-z0-9_]*\s*\(\)\s*\{",
            ),
            ("conditional", r"\b(?:if|unless)\b"),
            ("loop", r"\b(?:for|while|until)\b"),
        )
        if any(re.search(pattern, line) for line in setup_lines)
    )

    temporary_resource_shape = tuple(
        shape
        for shape, pattern in (
            ("temporary_path", r"\bmktemp\b"),
            ("temporary_directory", r"\bmkdir\s+-p\b"),
            ("loopback_resource", r"127\.0\.0\.1|\blocalhost\b"),
            ("background_process", r"(?:^|\s)&(?:\s|$)"),
        )
        if re.search(pattern, normalized)
    )
    control_shape = tuple(
        keyword
        for keyword in ("if", "for", "while", "until", "case")
        if re.search(rf"\b{keyword}\b", normalized)
    )

    owner_source = next((line for line in semantic_code if owner in line), "")
    owner_prefix, _, owner_suffix = owner_source.partition(owner)

    if re.search(r"\b(?:if|elif)\s+!\s+", owner_prefix):
        expected_outcome_shape = ("guarded_rejection",)
    elif re.search(r"(?:^|[;&|]\s*)!\s+", owner_prefix):
        expected_outcome_shape = ("negated_rejection",)
    elif re.search(r"\b(?:if|elif)\s+", owner_prefix):
        expected_outcome_shape = ("guarded_success",)
    elif "||" in owner_suffix:
        expected_outcome_shape = ("failure_fallback",)
    elif "&&" in owner_suffix:
        expected_outcome_shape = ("success_chain",)
    else:
        expected_outcome_shape = ("direct",)

    skip_failure_shape = tuple(
        shape
        for shape, pattern in (
            ("skip", r"\brecord_skip\b|\bSKIP:"),
            (
                "expected_failure",
                r"\b(?:if|elif)\s+!\s+|(?:^|[;&|]\s*)!\s+|\|\|",
            ),
            ("status_capture", r"\b(?:rc|status)=\$\?|\bset\s+\+e\b"),
        )
        if re.search(pattern, normalized)
    )
    cleanup_shape = tuple(
        shape
        for shape, pattern in (
            ("exit_trap", r"\btrap\b[^\n]*\bEXIT\b"),
            ("path_cleanup", r"\brm\s+-[A-Za-z]+\b"),
            ("process_cleanup", r"\b(?:kill|wait|cleanup_server_http)\b"),
        )
        if re.search(pattern, normalized)
    )

    signature: tuple[object, ...] = (
        owner,
        tuple(positionals),
        tuple(positional_selectors),
        tuple(sorted(set(options))),
        tuple(sorted(selectors)),
        tuple(sorted(value_shapes)),
        operators,
        setup_shape,
        control_shape,
        expected_outcome_shape,
        temporary_resource_shape,
        skip_failure_shape,
        cleanup_shape,
    )

    return (
        normalized,
        signature,
        tuple(sorted(set(options))),
        tuple(positionals),
        expected_outcome_shape,
    )


def _finding(
    rule_id: str,
    path: str,
    line: int,
    owner: str,
    message: str,
    offset: int,
) -> Finding:
    """
    Build one finding with a source-line offset.
    """

    return Finding(rule_id, path, line + offset, owner, message)


def _parse_example_entry(
    entry_lines: list[tuple[int, str]],
    match: re.Match[str],
    source_line: int,
    *,
    owner: str,
    accepted_aliases: set[str],
    public_aliases: set[str],
    hidden_aliases: set[str],
    path: str,
    line_offset: int,
    strict_findings: list[Finding],
    reviews: list[Review],
) -> Example | None:
    """
    Validate and parse one numbered help example.

    Parameters
    ----------
    entry_lines : list[tuple[int, str]]
        Physical rows governed by the numbered entry.
    match : re.Match[str]
        Parsed entry heading.
    source_line : int
        Entry source line.
    owner : str
        Documented command owner.
    accepted_aliases : set[str]
        Parser-accepted aliases.
    public_aliases : set[str]
        Publicly documented aliases.
    hidden_aliases : set[str]
        Accepted compatibility aliases.
    path : str
        Repository-relative source path.
    line_offset : int
        Rendered-to-source line offset.
    strict_findings : list[Finding]
        Mutable deterministic findings.
    reviews : list[Review]
        Mutable semantic-review candidates.

    Returns
    -------
    example : Example | None
        Parsed example, or None after an unusable code-block finding.
    """

    openings = [
        index
        for index, (_, line) in enumerate(entry_lines)
        if line == "    '''bash"
    ]
    closings = [
        index
        for index, (_, line) in enumerate(entry_lines)
        if line == "    '''"
    ]

    if len(openings) != 1 or len(closings) != 1:
        strict_findings.append(
            _finding(
                RULE_CODE_BLOCK,
                path,
                source_line,
                owner,
                "each numbered entry requires exactly one "
                "'''bash ... ''' block",
                line_offset,
            ),
        )

        return None

    opening, closing = openings[0], closings[0]

    if closing <= opening + 1:
        strict_findings.append(
            _finding(
                RULE_COMPLETE,
                path,
                source_line,
                owner,
                "example code block must be closed and nonempty",
                line_offset,
            ),
        )

        return None

    code_rows = entry_lines[opening + 1 : closing]

    if any(line and not line.startswith("    ") for _, line in code_rows):
        strict_findings.append(
            _finding(
                RULE_COMPLETE,
                path,
                source_line,
                owner,
                (
                    "example code must use the established four-space "
                    "indentation"
                ),
                line_offset,
            ),
        )

    code = tuple(
        line[4:] if line.startswith("    ") else line for _, line in code_rows
    )
    owner_pattern = re.compile(
        rf"(?<![A-Za-z0-9_]){re.escape(owner)}(?![A-Za-z0-9_])",
    )
    owner_hits = [line for line in code if owner_pattern.search(line)]

    if not owner_hits:
        strict_findings.append(
            _finding(
                RULE_OWNER,
                path,
                source_line,
                owner,
                "example code block must invoke the documented owner",
                line_offset,
            ),
        )

    observed_options = set(OPTION.findall("\n".join(code)))
    has_nonhelp_parser_options = bool(
        accepted_aliases
        - {
            "-h",
            "--hlp",
            "--help",
            "-d",
            "--details",
            "-ah",
            "--all_hlp",
            "--all_help",
            "--all-hlp",
            "--all-help",
        },
    )

    if not has_nonhelp_parser_options:
        observed_options &= accepted_aliases

    for alias in sorted(observed_options & hidden_aliases):
        strict_findings.append(
            _finding(
                RULE_ALIAS_VISIBILITY,
                path,
                source_line,
                owner,
                f"hidden compatibility alias '{alias}' may not appear "
                "in examples",
                line_offset,
            ),
        )

    for alias in sorted(observed_options - accepted_aliases):
        strict_findings.append(
            _finding(
                RULE_ALIAS_ACCEPTANCE,
                path,
                source_line,
                owner,
                f"option '{alias}' is not accepted by the relevant parser",
                line_offset,
            ),
        )

    for alias in sorted(
        (observed_options & accepted_aliases) - public_aliases,
    ):
        if alias not in hidden_aliases:
            strict_findings.append(
                _finding(
                    RULE_ALIAS_VISIBILITY,
                    path,
                    source_line,
                    owner,
                    (
                        f"accepted nonpublic alias '{alias}' may not appear "
                        "in examples"
                    ),
                    line_offset,
                ),
            )

    (
        normalized,
        signature,
        option_names,
        positional_shape,
        expected_outcome_shape,
    ) = invocation_facts(code, owner)
    example = Example(
        number=int(match.group("number")),
        line=source_line + line_offset,
        description=match.group("description"),
        code=code,
        normalized_invocation=normalized,
        signature=signature,
        option_names=option_names,
        positional_shape=positional_shape,
        setup_global_state_shape=signature[7],
        temporary_resource_shape=signature[10],
        control_flow_shape=signature[8],
        expected_outcome_shape=expected_outcome_shape,
        skip_failure_shape=signature[11],
        cleanup_shape=signature[12],
    )

    if re.search(r"\brm\s+-rf\b|\b(?:sudo|mkfs|shutdown)\b", normalized):
        reviews.append(
            Review(
                RULE_REVIEW_UNSAFE,
                path,
                source_line + line_offset,
                owner,
                (example.number,),
                "example contains a potentially unsafe command",
            ),
        )

    return example


def _analyze_examples_section(
    section: Section,
    *,
    owner: str,
    accepted_aliases: set[str],
    public_aliases: set[str],
    hidden_aliases: set[str],
    path: str,
    line_offset: int,
    strict_findings: list[Finding],
    reviews: list[Review],
) -> list[Example]:
    """
    Validate one Examples section and return its parsed entries.
    """

    body = list(section.lines)
    starts = [
        (index, number, ENTRY.fullmatch(line))
        for index, (number, line) in enumerate(body)
        if ENTRY.fullmatch(line)
    ]

    if not starts and any(line.strip() for _, line in body):
        strict_findings.append(
            _finding(
                RULE_ENTRY,
                path,
                section.heading_line,
                owner,
                "Examples content must use top-level numbered entries",
                line_offset,
            ),
        )
    elif starts and any(
        line.strip() and index < starts[0][0]
        for index, (_, line) in enumerate(body)
    ):
        strict_findings.append(
            _finding(
                RULE_ENTRY,
                path,
                section.heading_line,
                owner,
                "unnumbered content may not precede the first example entry",
                line_offset,
            ),
        )

    observed_numbers = [
        int(match.group("number")) for _, _, match in starts if match
    ]

    if observed_numbers != list(range(1, len(observed_numbers) + 1)):
        strict_findings.append(
            _finding(
                RULE_NUMBERING,
                path,
                starts[0][1] if starts else section.heading_line,
                owner,
                "numbering must be sequential and begin at 1",
                line_offset,
            ),
        )

    if len(starts) < 2:
        strict_findings.append(
            _finding(
                RULE_COUNT,
                path,
                section.heading_line,
                owner,
                "strict full help requires at least two numbered examples",
                line_offset,
            ),
        )

    parsed: list[Example] = []

    for position, (start, source_line, match) in enumerate(starts):
        assert match is not None

        end = (
            starts[position + 1][0]
            if position + 1 < len(starts)
            else len(body)
        )
        example = _parse_example_entry(
            body[start + 1 : end],
            match,
            source_line,
            owner=owner,
            accepted_aliases=accepted_aliases,
            public_aliases=public_aliases,
            hidden_aliases=hidden_aliases,
            path=path,
            line_offset=line_offset,
            strict_findings=strict_findings,
            reviews=reviews,
        )

        if example is not None:
            parsed.append(example)

    return parsed


def _compare_examples(
    examples: list[Example],
    *,
    owner: str,
    path: str,
    line_offset: int,
    strict_findings: list[Finding],
    reviews: list[Review],
) -> None:
    """
    Record duplicate and material-distinctness decisions for example pairs.
    """

    for index, first in enumerate(examples):
        for second in examples[index + 1 :]:
            if first.normalized_invocation == second.normalized_invocation:
                strict_findings.append(
                    _finding(
                        RULE_DUPLICATE,
                        path,
                        second.line - line_offset,
                        owner,
                        (
                            f"examples {first.number} and {second.number} "
                            "repeat one invocation"
                        ),
                        line_offset,
                    ),
                )

            if first.signature == second.signature:
                strict_findings.append(
                    _finding(
                        RULE_SIGNATURE_DUPLICATE,
                        path,
                        second.line - line_offset,
                        owner,
                        (
                            f"examples {first.number} and {second.number} "
                            "have identical normalized signatures"
                        ),
                        line_offset,
                    ),
                )

            same_branch_evidence = (
                first.signature[4] == second.signature[4]
                and first.signature[5] == second.signature[5]
                and first.signature[6] == second.signature[6]
            )
            outcomes_differ = (
                first.expected_outcome_shape != second.expected_outcome_shape
            )
            explicitly_distinct_outcomes = (
                same_branch_evidence and outcomes_differ
            )

            if (
                first.option_names == second.option_names
                and first.positional_shape == second.positional_shape
                and first.signature[2] == second.signature[2]
                and first.signature[5] == second.signature[5]
                and not explicitly_distinct_outcomes
            ):
                reviews.append(
                    Review(
                        RULE_REVIEW,
                        path,
                        second.line,
                        owner,
                        (first.number, second.number),
                        (
                            "examples share option and positional shape; "
                            "confirm that branch or domain meaning is "
                            "materially distinct"
                        ),
                    ),
                )


def analyze_help_document(
    text: str,
    *,
    owner: str,
    accepted_aliases: set[str],
    public_aliases: set[str],
    hidden_aliases: set[str],
    path: str = "<memory>",
    line_offset: int = 0,
) -> Analysis:
    """
    Analyze one full help document with no general shell parsing.

    Parameters
    ----------
    text : str
        Source text to inspect or normalize.
    owner : str
        Help owner whose source and rendered forms are analyzed.
    accepted_aliases : set[str]
        Option aliases accepted by the documented command.
    public_aliases : set[str]
        Option aliases documented as public command surfaces.
    hidden_aliases : set[str]
        Option aliases accepted but omitted from public help.
    path : str
        Repository-relative path associated with the source.
    line_offset : int
        Source-line offset applied to rendered findings.

    Returns
    -------
    analysis : Analysis
        Findings and semantic review facts for the help document.
    """

    heredoc = rendered_heredoc(text, owner)
    parsed = sections(heredoc)
    example_sections = [
        section for section in parsed if section.name == "Examples"
    ]
    strict_findings: list[Finding] = []
    reviews: list[Review] = []
    examples_parsed: list[Example] = []

    if len(example_sections) != 1:
        rule = RULE_REQUIRED if not example_sections else RULE_SECTION_COUNT
        strict_findings.append(
            _finding(
                rule,
                path,
                example_sections[0].heading_line if example_sections else 1,
                owner,
                "strict full help requires exactly one Examples section",
                line_offset,
            ),
        )

    if example_sections and parsed[-1].name != "Examples":
        strict_findings.append(
            _finding(
                RULE_FINAL,
                path,
                example_sections[-1].heading_line,
                owner,
                "Examples must be the final top-level section",
                line_offset,
            ),
        )

    if example_sections:
        examples_parsed.extend(
            _analyze_examples_section(
                example_sections[0],
                owner=owner,
                accepted_aliases=accepted_aliases,
                public_aliases=public_aliases,
                hidden_aliases=hidden_aliases,
                path=path,
                line_offset=line_offset,
                strict_findings=strict_findings,
                reviews=reviews,
            ),
        )

    _compare_examples(
        examples_parsed,
        owner=owner,
        path=path,
        line_offset=line_offset,
        strict_findings=strict_findings,
        reviews=reviews,
    )

    status = "strict_green" if not strict_findings else "strict_violation"
    return Analysis(strict_findings, reviews, examples_parsed, status)


def classify_wrapper_source(path: str, text: str) -> Ownership:
    """
    Classify one execute/submit wrapper's bounded full-help surface.
    """

    delegate = DELEGATE.search(text)
    if delegate:
        return Ownership(
            path,
            "compatibility_delegate",
            "",
            "delegated canonical full help",
            delegate.group("target"),
        )

    stem = Path(path).stem

    if re.search(r"--details", text) and re.search(
        rf"\bdetail_{re.escape(stem)}\b",
        text,
    ):
        surface = "--details"

        if "--all_help" in text or "--all[_-]h[e]?lp" in text:
            surface += "; --all_help renders short plus full help"

        return Ownership(
            path,
            "details_full_document",
            f"detail_{stem}",
            surface,
        )

    help_owner = f"help_{stem}"
    if re.search(rf"\b{re.escape(help_owner)}\b", text):
        return Ownership(
            path,
            "default_help_full_document",
            help_owner,
            "--help",
        )

    return Ownership(path, "no_valid_full_help_owner", "", "")


def short_help_advertises_details(text: str, owner: str) -> bool:
    """
    Return whether concise help publicly names its full-help surface.
    """

    return any(
        "--details" in row.aliases
        for row in parameter_rows(text)
        if row.owner == owner
    )


def undispositioned_review_findings(
    identity: str,
    reviews: list[Review],
) -> list[Finding]:
    """
    Block strict-green status while required semantic candidates remain.
    """

    return [
        Finding(
            RULE_REVIEW_PENDING,
            review.path,
            review.line,
            identity,
            (
                f"semantic candidate {review.as_dict()['signature']} requires "
                f"an explicit disposition"
            ),
        )
        for review in reviews
    ]


def shell_help_paths(root: Path) -> list[Path]:
    """
    Return repository shell sources that can own function help.
    """

    paths = list((root / "bin").glob("*.sh"))
    paths.extend((root / "lib/bash").rglob("*.sh"))
    paths.extend((root / "install" / "scripts").glob("*.sh"))
    paths.extend((root / "tests").rglob("*.sh"))

    return sorted(
        path
        for path in paths
        if path.is_file() and "outputs" not in path.relative_to(root).parts
    )


def top_level_shell_paths(root: Path) -> list[Path]:
    """
    Return maintained script-interface candidates, excluding helper modules.
    """

    paths = [
        *(root / "bin").glob("*.sh"),
        *(root / "install/scripts").glob("*.sh"),
    ]
    paths.extend((root / "tests").glob("*.sh"))
    paths.extend(
        (
            root / "tests/support/clean_artifacts.sh",
            root / "tests/integration/slurm/run_wet_tests.sh",
        ),
    )

    return sorted(path for path in paths if path.is_file())


def centralized_script_help_units(
    root: Path,
) -> dict[str, tuple[Path, Path, Heredoc]]:
    """
    Map top-level script identities to externally owned full-help heredocs.
    """

    result: dict[str, tuple[Path, Path, Heredoc]] = {}

    for source in top_level_shell_paths(root):
        relative = str(source.relative_to(root))
        classified = (
            classify_wrapper_source(
                relative,
                source.read_text(encoding="utf-8"),
            )
            if source.name.startswith(("execute_", "submit_"))
            else None
        )
        if (
            classified
            and classified.classification == "compatibility_delegate"
        ):
            continue

        help_path = root / "lib/bash/help" / f"help_{source.name}"
        if not help_path.is_file():
            continue

        docs = extract_help_heredocs(help_path.read_text(encoding="utf-8"))
        preferred = (
            classified.full_help_owner
            if classified and classified.full_help_owner
            else f"help_{source.stem}"
        )
        heredoc = next(
            (item for item in docs if item.owner == preferred),
            None,
        )

        if heredoc is None and docs:
            heredoc = docs[-1]

        if heredoc is not None:
            result[relative] = (source, help_path, heredoc)

    return result


def inline_script_help_units(
    root: Path,
) -> dict[str, tuple[Path, Path, Heredoc]]:
    """
    Map top-level scripts with an inline full-help heredoc.
    """

    centralized = centralized_script_help_units(root)
    result: dict[str, tuple[Path, Path, Heredoc]] = {}

    for source in top_level_shell_paths(root):
        relative = str(source.relative_to(root))
        if relative in centralized:
            continue

        text = source.read_text(encoding="utf-8")

        if source.name.startswith(("execute_", "submit_")):
            classified = classify_wrapper_source(relative, text)
            if classified.classification == "compatibility_delegate":
                continue

        docs = extract_help_heredocs(text)
        candidates = [
            item
            for item in docs
            if item.owner == "<file>"
            or item.owner.startswith(("help_", "show_help", "detail_"))
            or item.owner in {"main", "check_args_light"}
        ]

        if candidates:
            result[relative] = (source, source, candidates[-1])

    return result


def script_help_units(root: Path) -> dict[str, tuple[Path, Path, Heredoc]]:
    """
    Return every maintained top-level shell-script full-help owner.
    """

    return {
        **centralized_script_help_units(root),
        **inline_script_help_units(root),
    }


def function_help_units(root: Path) -> dict[str, tuple[Path, Heredoc]]:
    """
    Return stable identities for all recognized non-main function docs.
    """

    result: dict[str, tuple[Path, Heredoc]] = {}
    centralized = {
        (str(help_path.relative_to(root)), heredoc.owner, heredoc.start_line)
        for _, help_path, heredoc in script_help_units(root).values()
    }

    for path in shell_help_paths(root):
        relative = str(path.relative_to(root))

        for heredoc in extract_help_heredocs(path.read_text(encoding="utf-8")):
            if (
                heredoc.owner == "<file>"
                or (relative, heredoc.owner, heredoc.start_line) in centralized
            ):
                continue

            result[f"{relative}::{heredoc.owner}"] = (path, heredoc)

    return result


def repository_crosswalk(root: Path) -> dict[str, object]:
    """
    Return one deterministic help-owner universe and relationship crosswalk.
    """

    root = root.resolve()
    parsed: list[str] = []
    shared: list[str] = []

    for path in shell_help_paths(root):
        relative = str(path.relative_to(root))

        for heredoc in extract_help_heredocs(path.read_text(encoding="utf-8")):
            identity = f"{relative}::{heredoc.owner}@{heredoc.start_line}"
            parsed.append(identity)

            if heredoc.owner == "<file>":
                shared.append(identity)

    functions = function_help_units(root)
    scripts = script_help_units(root)
    centralized = centralized_script_help_units(root)
    wrappers = [
        item
        for item in (
            classify_wrapper_source(
                str(path.relative_to(root)),
                path.read_text(encoding="utf-8"),
            )
            for path in wrapper_paths(root)
        )
        if item.classification != "compatibility_delegate"
    ]
    delegators = [
        item
        for item in (
            classify_wrapper_source(
                str(path.relative_to(root)),
                path.read_text(encoding="utf-8"),
            )
            for path in wrapper_paths(root)
        )
        if item.classification == "compatibility_delegate"
    ]

    def category(
        name: str,
        definition: str,
        identities: list[str],
        relationship: str,
        governed: bool,
    ) -> dict[str, object]:
        return {
            "category": name,
            "definition": definition,
            "count": len(identities),
            "identities": sorted(identities),
            "overlap_or_delegation": relationship,
            "governed_by_two_example_rule": governed,
        }

    categories = [
        category(
            "all_parsed_help_bearing_source_units",
            (
                "Every bounded sectioned help heredoc parsed from maintained "
                "shell sources."
            ),
            parsed,
            (
                "Universe row; owner and non-owner categories below partition "
                "or overlap it explicitly."
            ),
            False,
        ),
        category(
            "shell_function_heredoc_owners",
            (
                "Parsed heredocs lexically owned by shell functions after "
                "centralized script help and shared fragments are removed."
            ),
            list(functions),
            "Disjoint from top-level script owners and shared heredocs.",
            True,
        ),
        category(
            "top_level_shell_script_help_owners",
            (
                "Maintained top-level script interfaces with one inline or "
                "centralized full-help document."
            ),
            list(scripts),
            (
                "Disjoint owner class; centralized functions delegate "
                "ownership to these script identities."
            ),
            True,
        ),
        category(
            "centralized_help_functions",
            (
                "Functions in lib/bash/help that render a top-level script's "
                "help."
            ),
            [
                f"{help.relative_to(root)}::{doc.owner}"
                for _, help, doc in centralized.values()
            ],
            (
                "Overlaps parsed units but delegates ownership to "
                "top_level_shell_script_help_owners."
            ),
            False,
        ),
        category(
            "execute_submit_rendered_wrapper_surfaces",
            (
                "Canonical execute/submit wrapper full-help surfaces "
                "classified from current parser dispatch."
            ),
            [f"{item.path}::{item.full_help_owner}" for item in wrappers],
            (
                "Subset of top_level_shell_script_help_owners; compatibility "
                "delegators are separate."
            ),
            True,
        ),
        category(
            "compatibility_delegators",
            (
                "Execute/submit shims that exec one canonical wrapper and own "
                "no separate full-help document."
            ),
            [f"{item.path}->{item.target}" for item in delegators],
            "Delegate to the target wrapper and do not add a second owner.",
            False,
        ),
        category(
            "shared_heredoc_nonowners",
            (
                "Parsed section fragments outside a function that are "
                "interpolated into one or more owned help documents."
            ),
            shared,
            (
                "Overlap parsed universe only; not independently callable "
                "help owners."
            ),
            False,
        ),
        category(
            "excluded_test_or_synthetic_fixtures",
            (
                "Generated outputs and Python string fixtures excluded by "
                "bounded repository discovery."
            ),
            [],
            (
                "Outside the parsed source universe; real tests shell "
                "functions remain governed."
            ),
            False,
        ),
    ]

    return {
        "schema_version": 1,
        "categories": categories,
        "totals": {
            "parsed_help_units": len(parsed),
            "script_help_owners": len(scripts),
            "function_help_owners": len(functions),
        },
    }


def accepted_function_aliases(text: str, owner: str) -> set[str]:
    """
    Return bounded parser-accepted spellings from one complete function body.
    """

    accepted = {
        alias
        for chunk in alias_chunks(text, {owner})
        for alias in chunk.aliases
    }
    body = function_bodies(text).get(owner, "")

    if re.search(r"-h\|--h\[e\]\?lp|-h\|--hlp\|--help", body):
        accepted.update({"-h", "--hlp", "--help"})

    return accepted


def documented_function_aliases(text: str, owner: str) -> set[str]:
    """
    Return public aliases from the function's Parameters rows.
    """

    return {
        alias
        for row in parameter_rows(text)
        if row.owner == owner
        for alias in row.aliases
    }


def wrapper_paths(root: Path) -> list[Path]:
    """
    Return all execute and submit wrappers, including compatibility shims.
    """

    return sorted(
        [
            *(root / "bin").glob("execute_*.sh"),
            *(root / "bin").glob("submit_*.sh"),
        ],
    )


def wrapper_aliases(
    root: Path,
    wrapper: Path,
) -> tuple[set[str], set[str], set[str]]:
    """
    Return accepted, public, and hidden aliases for one wrapper parser.
    """

    text = wrapper.read_text(encoding="utf-8")
    accepted = {
        alias
        for chunk in alias_chunks(
            text,
            {"parse_args", "resolve_dir_scr", "check_args_light", "main"},
        )
        for alias in chunk.aliases
    }

    if "--details" in text:
        accepted.update({"-d", "--details"})

    if "--all[_-]h[e]?lp" in text or "--all_help" in text:
        accepted.update({"-ah", "--all_help", "--all-help"})

    if re.search(r"-h\|--h\[e\]\?lp|-h\|--hlp\|--help", text):
        accepted.update({"-h", "--hlp", "--help"})

    hidden = {
        alias
        for alias in accepted
        if hidden_alias(alias, tuple(sorted(accepted)))
    }
    help_path = root / "lib/bash/help" / f"help_{wrapper.name}"
    public = set()

    if help_path.is_file():
        public = {
            alias
            for row in parameter_rows(help_path.read_text())
            for alias in row.aliases
        }

    return accepted, public, hidden


def script_aliases(
    source: Path,
    documentation: Path,
) -> tuple[set[str], set[str], set[str]]:
    """
    Return bounded parser and documentation aliases for one script owner.
    """

    text = source.read_text(encoding="utf-8")
    chunks = alias_chunks(
        text,
        {"parse_args", "resolve_dir_scr", "check_args_light", "main"},
    )

    if not chunks:
        chunks = file_alias_chunks(text)

    accepted = {alias for chunk in chunks for alias in chunk.aliases}

    if re.search(r"-h\|--h\[e\]\?lp|-h\|--hlp\|--help", text):
        accepted.update({"-h", "--hlp", "--help"})

    public = {
        alias
        for row in parameter_rows(documentation.read_text(encoding="utf-8"))
        for alias in row.aliases
    }
    hidden = {
        alias
        for alias in accepted
        if hidden_alias(alias, tuple(sorted(accepted)))
    }

    return accepted, public, hidden


def scan_repository(root: Path) -> RepositoryResult:
    """
    Audit every source-derived script and function help owner directly.

    Parameters
    ----------
    root : Path
        Repository root containing maintained help owners and registrations.

    Returns
    -------
    result : RepositoryResult
        Complete owner, finding, alias, and registration reconciliation.
    """

    root = root.resolve()
    findings: list[Finding] = []
    reviews: list[Review] = []
    inventory: list[dict[str, object]] = []
    ownership = [
        classify_wrapper_source(
            str(path.relative_to(root)),
            path.read_text(encoding="utf-8"),
        )
        for path in wrapper_paths(root)
    ]
    ownership_by_name = {Path(item.path).name: item for item in ownership}
    analyzed_canonical: dict[str, Analysis] = {}

    for item in ownership:
        canonical = item

        if item.classification == "compatibility_delegate":
            canonical = ownership_by_name.get(
                item.target,
                Ownership(item.path, "no_valid_full_help_owner", "", ""),
            )

        if canonical.classification == "no_valid_full_help_owner":
            findings.append(
                Finding(
                    RULE_OWNERSHIP,
                    item.path,
                    1,
                    Path(item.path).name,
                    "no valid full-help owner exists",
                ),
            )
            inventory.append(
                {**item.as_dict(), "status": "ownership_violation"},
            )

            continue

        canonical_wrapper = root / canonical.path
        canonical_name = canonical_wrapper.name

        if canonical_name not in analyzed_canonical:
            help_path = root / "lib/bash/help" / f"help_{canonical_name}"
            docs = extract_help_heredocs(help_path.read_text(encoding="utf-8"))
            heredoc = next(
                (
                    doc
                    for doc in docs
                    if doc.owner == canonical.full_help_owner
                ),
                None,
            )

            if heredoc is None:
                analysis = Analysis(
                    [
                        Finding(
                            RULE_OWNERSHIP,
                            str(help_path.relative_to(root)),
                            1,
                            canonical.full_help_owner,
                            "classified full-help heredoc was not found",
                        ),
                    ],
                    [],
                    [],
                    "strict_violation",
                )
            else:
                accepted, public, hidden = wrapper_aliases(
                    root,
                    canonical_wrapper,
                )
                analysis = analyze_help_document(
                    "\n".join(line for _, line in heredoc.lines),
                    owner=canonical_name,
                    accepted_aliases=accepted,
                    public_aliases=public,
                    hidden_aliases=hidden,
                    path=str(help_path.relative_to(root)),
                    line_offset=heredoc.lines[0][0] - 1,
                )

                if canonical.classification == "details_full_document":
                    short_owner = f"help_{canonical_wrapper.stem}"

                    if not short_help_advertises_details(
                        help_path.read_text(encoding="utf-8"),
                        short_owner,
                    ):
                        analysis = Analysis(
                            [
                                *analysis.findings,
                                Finding(
                                    RULE_OWNERSHIP,
                                    str(help_path.relative_to(root)),
                                    1,
                                    short_owner,
                                    (
                                        "concise --help must advertise the "
                                        "valid --details full-help surface"
                                    ),
                                ),
                            ],
                            analysis.reviews,
                            analysis.examples,
                            "strict_violation",
                        )

            pending = undispositioned_review_findings(
                canonical.path,
                analysis.reviews,
            )

            if pending:
                analysis = Analysis(
                    [*analysis.findings, *pending],
                    analysis.reviews,
                    analysis.examples,
                    "semantic_review_required",
                )

            analyzed_canonical[canonical_name] = analysis
            findings.extend(analysis.findings)
            reviews.extend(analysis.reviews)

        analysis = analyzed_canonical[canonical_name]
        inventory.append(
            {
                **item.as_dict(),
                "canonical_wrapper": canonical_name,
                "status": analysis.status,
                "examples": [
                    example.as_dict() for example in analysis.examples
                ],
            },
        )

    script_units = script_help_units(root)

    for identity, (source, help_path, heredoc) in sorted(script_units.items()):
        canonical = analyzed_canonical.get(source.name)

        if canonical is None:
            accepted, public, hidden = script_aliases(source, help_path)
            canonical = analyze_help_document(
                "\n".join(line for _, line in heredoc.lines),
                owner=source.name,
                accepted_aliases=accepted,
                public_aliases=public,
                hidden_aliases=hidden,
                path=str(help_path.relative_to(root)),
                line_offset=heredoc.lines[0][0] - 1,
            )
            pending = undispositioned_review_findings(
                identity,
                canonical.reviews,
            )

            if pending:
                canonical = Analysis(
                    [*canonical.findings, *pending],
                    canonical.reviews,
                    canonical.examples,
                    "semantic_review_required",
                )

            findings.extend(canonical.findings)
            reviews.extend(canonical.reviews)

        inventory.append(
            {
                "path": str(source.relative_to(root)),
                "documentation_path": str(help_path.relative_to(root)),
                "owner": heredoc.owner,
                "identity": identity,
                "kind": "script",
                "status": canonical.status,
                "examples": [
                    example.as_dict() for example in canonical.examples
                ],
            },
        )

    function_units = function_help_units(root)

    for identity, (path, heredoc) in sorted(function_units.items()):
        text = path.read_text(encoding="utf-8")
        accepted = accepted_function_aliases(text, heredoc.owner)
        public = documented_function_aliases(text, heredoc.owner)
        hidden = {
            alias
            for alias in accepted
            if hidden_alias(alias, tuple(sorted(accepted)))
        }

        analysis = analyze_help_document(
            "\n".join(line for _, line in heredoc.lines),
            owner=heredoc.owner,
            accepted_aliases=accepted,
            public_aliases=public,
            hidden_aliases=hidden,
            path=str(path.relative_to(root)),
            line_offset=heredoc.lines[0][0] - 1,
        )
        findings.extend(analysis.findings)
        status = analysis.status
        pending = undispositioned_review_findings(identity, analysis.reviews)

        if pending:
            status = "semantic_review_required"
            findings.extend(pending)

        reviews.extend(analysis.reviews)
        inventory.append(
            {
                "path": str(path.relative_to(root)),
                "owner": heredoc.owner,
                "identity": identity,
                "kind": "function",
                "status": status,
                "examples": [
                    example.as_dict() for example in analysis.examples
                ],
            },
        )

    expected = set(script_units) | set(function_units)
    observed = {
        str(row["identity"])
        for row in inventory
        if row.get("kind") in {"script", "function"}
    }

    for identity in sorted(expected - observed):
        findings.append(
            Finding(
                RULE_OWNERSHIP,
                identity.split("::", 1)[0],
                1,
                identity,
                (
                    "source-derived help owner is absent from the strict "
                    "inventory"
                ),
            ),
        )

    alias_findings, _ = scan_alias_repository(root)

    return RepositoryResult(
        findings,
        reviews,
        inventory,
        ownership,
        list(alias_findings),
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse command-line arguments.

    Parameters
    ----------
    argv : list[str] | None
        Explicit arguments, or None to read the process arguments.

    Returns
    -------
    arguments : argparse.Namespace
        Parsed Examples-audit paths and output options.
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--inventory-output", type=Path)
    parser.add_argument("--semantic-output", type=Path)
    parser.add_argument("--ownership-output", type=Path)
    parser.add_argument("--crosswalk-output", type=Path)

    return parser.parse_args(argv)


def write_json(path: Path | None, value: object) -> None:
    """
    Write one deterministic diagnostic artifact when requested.
    """

    if path is not None:
        path.write_text(
            json.dumps(value, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )


def compliance_summary(result: RepositoryResult) -> dict[str, object]:
    """
    Compute final repository-wide compliance from source-derived results.
    """

    scripts = [row for row in result.inventory if row.get("kind") == "script"]
    functions = [
        row for row in result.inventory if row.get("kind") == "function"
    ]
    script_green = sum(row["status"] == "strict_green" for row in scripts)
    function_green = sum(row["status"] == "strict_green" for row in functions)
    remaining = len(scripts) + len(functions) - script_green - function_green
    structural = [
        finding
        for finding in result.findings
        if finding.rule_id != RULE_REVIEW_PENDING
    ]
    compliant = (
        script_green == len(scripts)
        and function_green == len(functions)
        and remaining == 0
        and not structural
        and not result.alias_findings
        and not result.reviews
    )

    return {
        "script_green": script_green,
        "script_total": len(scripts),
        "function_green": function_green,
        "function_total": len(functions),
        "remaining": remaining,
        "structural_findings": len(structural),
        "alias_findings": len(result.alias_findings),
        "required_semantic_candidates": len(result.reviews),
        "global_compliance": compliant,
    }


def main(argv: list[str] | None = None) -> int:
    """
    Run the authoritative repository Examples audit.

    Parameters
    ----------
    argv : list[str] | None
        Explicit arguments, or None to read the process arguments.

    Returns
    -------
    status : int
        Zero when Examples contracts pass and one when findings remain.
    """

    args = parse_args(argv)
    result = scan_repository(args.root)

    write_json(args.inventory_output, result.inventory)
    write_json(
        args.semantic_output,
        [review.as_dict() for review in result.reviews],
    )
    write_json(
        args.ownership_output,
        [owner.as_dict() for owner in result.ownership],
    )
    write_json(args.crosswalk_output, repository_crosswalk(args.root))

    for finding in result.findings:
        print(finding.format())

    for finding in result.alias_findings:
        print(finding.format())

    summary = compliance_summary(result)
    script_score = f"{summary['script_green']}/{summary['script_total']}"
    function_score = f"{summary['function_green']}/{summary['function_total']}"
    print(f"script help owners: {script_score} strict-green")
    print(f"function help owners: {function_score} strict-green")
    print(f"remaining owners: {summary['remaining']}")
    answer = "yes" if summary["global_compliance"] else "no"
    print(f"global Examples compliance: {answer}")

    if result.findings or result.alias_findings:
        count = len(result.findings) + len(result.alias_findings)
        print(f"HELP.EXAMPLES: {count} strict violation(s)")

        return 1

    print("HELP.EXAMPLES: pass (repository-wide strict owner set green)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
