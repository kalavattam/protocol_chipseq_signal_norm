#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: ai_attribution.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


"""Validate and normalize bounded language-neutral source headers."""

from __future__ import annotations

import argparse
import dataclasses
import json
import re
import sys
import textwrap
from collections import Counter
from pathlib import Path

from help_heredoc_reflow import run_git, shell_paths

RULE_ATTRIBUTION = "SOURCE.HEADER.AI_ATTRIBUTION"
RULE_ATTRIBUTION_UNIQUE = "SOURCE.HEADER.AI_ATTRIBUTION.UNIQUE"
RULE_HEADER = "SOURCE.HEADER.STRUCTURE"
RULE_BASENAME = "SOURCE.HEADER.BASENAME"
RULE_WIDTH = "SOURCE.HEADER.WIDTH"
RULE_YEAR = "SOURCE.HEADER.YEAR"
CURRENT_YEAR = 2026
ATTRIBUTION_START = re.compile(
    r"^# OpenAI (?:ChatGPT and Codex|ChatGPT|Codex)\b"
)
ATTRIBUTION_LIKE = re.compile(
    r"^# .*\b(?:AI|OpenAI|ChatGPT|Codex|GPT-)\b",
    re.IGNORECASE,
)
ATTRIBUTION_FORM = re.compile(
    r"^OpenAI (?P<tools>ChatGPT and Codex|ChatGPT|Codex) "
    r"\((?P<models>[^()]*)\) (?P<verb>was|were) used in "
    r"(?P<domain>development and documentation|development|documentation)\.$"
)
MODEL_IDENTIFIER_TEXT = (
    r"GPT-\d+(?:\.\d+)?(?:-(?!series\b)[A-Za-z0-9]+)?"
)
SERIES_IDENTIFIER_TEXT = r"GPT-\d+(?:\.\d+)?-series"
MODEL_TOKEN = re.compile(
    rf"(?<![A-Za-z0-9.-]){MODEL_IDENTIFIER_TEXT}(?![A-Za-z0-9.-])"
)
EXPLICIT_MODEL = re.compile(MODEL_IDENTIFIER_TEXT)
EXPLICIT_MODEL_LIST = re.compile(
    rf"{MODEL_IDENTIFIER_TEXT}(?:, {MODEL_IDENTIFIER_TEXT})*"
)
SINGLE_SERIES = re.compile(
    r"GPT-(?P<series>\d+(?:\.\d+)?)-series models"
)
COMBINED_SERIES = re.compile(
    r"GPT-(?P<first>\d+(?:\.\d+)?)- and "
    r"GPT-(?P<second>\d+(?:\.\d+)?)-series models"
)
COPYRIGHT = re.compile(
    r"^# Copyright (?P<start>\d{4})(?:-(?P<end>\d{4}))? by Kris Alavattam$"
)
EMAIL = "# Email: kalavattam@gmail.com"
LICENSE = "# Distributed under the MIT license."
ENCODING = "# -*- coding: utf-8 -*-"


@dataclasses.dataclass(frozen=True)
class Finding:
    """One bounded source-header finding."""

    rule_id: str
    path: str
    line: int
    message: str

    def format(self) -> str:
        """Render one stable diagnostic."""

        return f"{self.rule_id}: {self.path}:{self.line}: {self.message}"


@dataclasses.dataclass(frozen=True)
class AttributionRequirement:
    """Path-scoped attribution data supplied by an applicability manifest."""

    path: str
    required_models: tuple[str, ...]
    contribution_domain: str
    tools: str


@dataclasses.dataclass(frozen=True)
class AttributionObservation:
    """Parsed facts from one recognized OpenAI attribution block."""

    tools: str
    models: tuple[str, ...]
    contribution_domain: str
    attribution_style: str


@dataclasses.dataclass(frozen=True)
class ModelDeclaration:
    """One parsed explicit-list or generic-series model declaration."""

    style: str
    tokens: tuple[str, ...]


def parse_model_declaration(text: str) -> ModelDeclaration | None:
    """Parse one bounded model declaration without accepting residue."""

    if EXPLICIT_MODEL_LIST.fullmatch(text):
        return ModelDeclaration("explicit_model_list", tuple(text.split(", ")))

    separator = "; most recent: "
    if separator in text:
        if text.count(separator) != 1:
            return None
        series_text, recent_text = text.split(separator, 1)
        if not EXPLICIT_MODEL_LIST.fullmatch(recent_text):
            return None
        recent = tuple(recent_text.split(", "))
    else:
        if ";" in text:
            return None
        series_text = text
        recent = ()

    combined = COMBINED_SERIES.fullmatch(series_text)
    if combined is not None:
        series = (
            f"GPT-{combined.group('first')}-series",
            f"GPT-{combined.group('second')}-series",
        )
    else:
        single = SINGLE_SERIES.fullmatch(series_text)
        if single is None:
            return None
        series = (f"GPT-{single.group('series')}-series",)
    return ModelDeclaration("generic_series", series + recent)


def model_tokens(text: str) -> tuple[str, ...]:
    """Return exact explicit and normalized generic-series model tokens."""

    tokens: list[tuple[int, str]] = []
    excluded: list[tuple[int, int]] = []
    for match in COMBINED_SERIES.finditer(text):
        tokens.extend(
            (
                (match.start(), f"GPT-{match.group('first')}-series"),
                (match.start() + 1, f"GPT-{match.group('second')}-series"),
            )
        )
        excluded.append(match.span())
    for match in SINGLE_SERIES.finditer(text):
        if any(start <= match.start() < end for start, end in excluded):
            continue
        tokens.append((match.start(), f"GPT-{match.group('series')}-series"))
        excluded.append(match.span())
    for match in MODEL_TOKEN.finditer(text):
        if any(start <= match.start() < end for start, end in excluded):
            continue
        tokens.append((match.start(), match.group(0)))
    return tuple(token for _, token in sorted(tokens))


def valid_required_model(model: str) -> bool:
    """Return whether *model* is exactly one supported model token."""

    return bool(
        EXPLICIT_MODEL.fullmatch(model)
        or re.fullmatch(SERIES_IDENTIFIER_TEXT, model)
    )


def parse_attribution(rendered: str) -> AttributionObservation | None:
    """Parse one approved bounded OpenAI attribution form."""

    normalized = " ".join(
        line.removeprefix("#").strip() for line in rendered.splitlines()
    )
    match = ATTRIBUTION_FORM.fullmatch(normalized)
    if match is None:
        return None
    declaration = parse_model_declaration(match.group("models"))
    if declaration is None:
        return None
    tools = match.group("tools")
    expected_verb = "were" if tools == "ChatGPT and Codex" else "was"
    if match.group("verb") != expected_verb:
        return None
    parsed_tools = {
        "ChatGPT": "chatgpt",
        "Codex": "codex",
        "ChatGPT and Codex": "both",
    }[tools]
    return AttributionObservation(
        tools=parsed_tools,
        models=declaration.tokens,
        contribution_domain=(
            "both"
            if match.group("domain") == "development and documentation"
            else match.group("domain")
        ),
        attribution_style=declaration.style,
    )


def load_applicability_manifest(
    root: Path,
    path: Path,
) -> tuple[int, dict[str, AttributionRequirement]]:
    """Validate and return one explicit attribution applicability manifest."""

    root = root.resolve()
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("schema_version") != 1:
        raise ValueError("applicability manifest requires schema_version 1")
    assessment_year = data.get("assessment_year")
    if not isinstance(assessment_year, int):
        raise ValueError("applicability manifest requires integer assessment_year")
    requirements: dict[str, AttributionRequirement] = {}
    for row in data.get("sources", []):
        relative = row.get("path")
        if not isinstance(relative, str) or not relative:
            raise ValueError("manifest source path must be a nonempty string")
        candidate = Path(relative)
        if candidate.is_absolute() or ".." in candidate.parts:
            raise ValueError(f"unsafe manifest source path: {relative!r}")
        normalized = candidate.as_posix()
        if normalized != relative or normalized in requirements:
            raise ValueError(f"duplicate or noncanonical manifest path: {relative!r}")
        resolved = (root / candidate).resolve()
        if root != resolved and root not in resolved.parents:
            raise ValueError(f"manifest path escapes root: {relative!r}")
        models = row.get("required_models")
        if (
            not isinstance(models, list)
            or not models
            or not all(isinstance(model, str) and model for model in models)
            or len(set(models)) != len(models)
            or not all(valid_required_model(model) for model in models)
        ):
            raise ValueError(f"invalid required_models for {relative!r}")
        contribution_domain = row.get("contribution_domain")
        if contribution_domain not in {"development", "documentation", "both"}:
            raise ValueError(f"invalid contribution_domain for {relative!r}")
        tools = row.get("tools")
        if tools not in {"chatgpt", "codex", "both"}:
            raise ValueError(f"invalid tools for {relative!r}")
        requirements[normalized] = AttributionRequirement(
            path=normalized,
            required_models=tuple(models),
            contribution_domain=contribution_domain,
            tools=tools,
        )
    return assessment_year, requirements


def attribution_block(text: str) -> tuple[int, int, str] | None:
    """Return the recognized attribution comment block."""

    lines = text.splitlines()
    upper = min(len(lines), 30)
    license_index = next(
        (index for index, line in enumerate(lines[:upper]) if line == LICENSE),
        None,
    )
    body_index = next(
        (
            index
            for index, line in enumerate(lines[1:upper], 1)
            if line and not line.startswith("#")
        ),
        None,
    )
    limits = [value for value in (license_index, body_index) if value is not None]
    upper = min(limits) if limits else upper
    for index, line in enumerate(lines[:upper]):
        if not ATTRIBUTION_START.match(line):
            continue
        end = index
        while end + 1 < len(lines) and lines[end + 1].startswith("#"):
            candidate = lines[end + 1]
            if candidate == "#":
                break
            end += 1
        return index, end, "\n".join(lines[index : end + 1])
    return None


def expected_basename(path: str, lines: list[str]) -> str | None:
    """Return the basename required by one explicit or in-memory source."""

    if path not in {"", "<memory>"}:
        return Path(path).name
    script = next((line for line in lines if line.startswith("# Script: ")), None)
    return script.removeprefix("# Script: ") if script else None


def check_attribution_source(
    text: str,
    path: str = "<memory>",
    *,
    current_year: int = CURRENT_YEAR,
    required_models: tuple[str, ...] = (),
    required_contribution_domain: str | None = None,
    required_tools: str | None = None,
) -> list[Finding]:
    """Return source-header and declared-attribution findings for one source."""

    lines = text.splitlines()
    if not lines or not lines[0].startswith("#!"):
        return []
    basename = expected_basename(path, lines)
    findings: list[Finding] = []
    license_index = next(
        (
            index
            for index, line in enumerate(lines[:30])
            if line == LICENSE
        ),
        None,
    )
    header_end = license_index if license_index is not None else min(len(lines) - 1, 20)
    header = lines[: header_end + 1]

    if basename == "install_envs_entrypoint.sh":
        if lines[0] != "#!/bin/sh":
            findings.append(Finding(RULE_HEADER, path, 1, "install_envs_entrypoint.sh must retain '#!/bin/sh'"))
    elif basename and basename.endswith(".py") and lines[0] != "#!/usr/bin/env python3":
        findings.append(Finding(RULE_HEADER, path, 1, "Python source requires '#!/usr/bin/env python3'"))
    elif basename and basename.endswith(".sh") and lines[0] != "#!/usr/bin/env bash":
        findings.append(Finding(RULE_HEADER, path, 1, "Bash source requires '#!/usr/bin/env bash'"))

    script_row = f"# Script: {basename}" if basename else None
    expected_rows = [ENCODING, script_row, EMAIL, LICENSE]
    for expected in expected_rows:
        if expected is None:
            continue
        matches = [index for index, line in enumerate(header) if line == expected]
        if len(matches) != 1:
            rule = RULE_BASENAME if expected.startswith("# Script:") else RULE_HEADER
            findings.append(Finding(rule, path, 1, f"source header requires exactly one exact row: {expected}"))
    copyright_rows = [
        (index, COPYRIGHT.fullmatch(line))
        for index, line in enumerate(header)
        if line.startswith("# Copyright ")
    ]
    if len(copyright_rows) != 1 or copyright_rows[0][1] is None:
        findings.append(Finding(RULE_HEADER, path, 1, "source header requires one canonical copyright row"))
        copyright_index = None
    else:
        copyright_index, match = copyright_rows[0]
        assert match is not None
        ending_year = int(match.group("end") or match.group("start"))
        if ending_year != current_year:
            findings.append(Finding(RULE_YEAR, path, copyright_index + 1, f"changed source copyright must end in {current_year}"))

    block = attribution_block(text)
    attribution_start = block[0] if block else None
    attribution_end = block[1] if block else None
    exact_positions = {
        1: ENCODING,
        2: "#",
        3: script_row,
        4: "#",
        6: EMAIL,
        7: "#",
    }
    for index, expected in exact_positions.items():
        if expected is not None and (
            index >= len(lines) or lines[index] != expected
        ):
            rule = RULE_BASENAME if index == 3 else RULE_HEADER
            findings.append(
                Finding(
                    rule,
                    path,
                    index + 1,
                    f"source-header row must be exactly: {expected}",
                )
            )
    if copyright_index != 5:
        findings.append(
            Finding(
                RULE_HEADER,
                path,
                6,
                "copyright must be the sixth header row",
            )
        )
    if attribution_start is not None and attribution_start != 8:
        findings.append(
            Finding(
                RULE_HEADER,
                path,
                attribution_start + 1,
                "attribution block must follow the email separator",
            )
        )
    if attribution_end is not None:
        separator = attribution_end + 1
        expected_license = attribution_end + 2
        if separator >= len(lines) or lines[separator] != "#":
            findings.append(Finding(RULE_HEADER, path, separator + 1, "attribution must be followed by exact '#' separator"))
        if expected_license >= len(lines) or lines[expected_license] != LICENSE:
            findings.append(Finding(RULE_HEADER, path, expected_license + 1, "license must immediately follow the attribution separator"))
        elif license_index != expected_license:
            findings.append(Finding(RULE_HEADER, path, expected_license + 1, "license row is out of order"))
    elif license_index != 8:
        findings.append(
            Finding(
                RULE_HEADER,
                path,
                9,
                "no-AI source profile requires the license after the email separator",
            )
        )

    if license_index is not None:
        body_index = license_index + 1
        blanks = 0
        while body_index < len(lines) and not lines[body_index]:
            blanks += 1
            body_index += 1
        if body_index < len(lines) and blanks != 2:
            findings.append(
                Finding(
                    RULE_HEADER,
                    path,
                    license_index + 2,
                    "source body must follow exactly two ordinary blank lines",
                )
            )

    for index, line in enumerate(header):
        if len(line) > 79:
            findings.append(Finding(RULE_WIDTH, path, index + 1, "source-header line exceeds 79 characters"))

    attribution_starts = [
        index for index, line in enumerate(header) if ATTRIBUTION_START.match(line)
    ]
    if len(attribution_starts) > 1:
        findings.append(
            Finding(
                RULE_ATTRIBUTION_UNIQUE,
                path,
                attribution_starts[1] + 1,
                "source header contains more than one attribution block",
            )
        )

    block = attribution_block(text)
    attribution_required = bool(
        required_models or required_contribution_domain or required_tools
    )
    if block is None:
        malformed = next(
            (
                index
                for index, line in enumerate(header)
                if ATTRIBUTION_LIKE.match(line)
            ),
            None,
        )
        if malformed is not None:
            findings.append(
                Finding(
                    RULE_ATTRIBUTION,
                    path,
                    malformed + 1,
                    "opening header contains an unrecognized attribution-like block",
                )
            )
        elif attribution_required:
            findings.append(
                Finding(
                    RULE_ATTRIBUTION,
                    path,
                    1,
                    "explicit applicability requires one recognized AI attribution block",
                )
            )
        return findings
    start, _, rendered = block
    observation = parse_attribution(rendered)
    if observation is None:
        findings.append(
            Finding(
                RULE_ATTRIBUTION,
                path,
                start + 1,
                "attribution must use an approved OpenAI tool, model, verb, and contribution-domain form",
            )
        )
        return findings
    declared_counts = Counter(observation.models)
    if not declared_counts:
        findings.append(
            Finding(
                RULE_ATTRIBUTION,
                path,
                start + 1,
                "recognized attribution must declare at least one model token",
            )
        )
    for model in sorted(set(declared_counts) | set(required_models)):
        count = declared_counts[model]
        if count > 1:
            findings.append(
                Finding(
                    RULE_ATTRIBUTION_UNIQUE,
                    path,
                    start + 1,
                    f"declared model {model!r} must appear exactly once",
                )
            )
        elif model in required_models and count == 0:
            findings.append(
                Finding(
                    RULE_ATTRIBUTION,
                    path,
                    start + 1,
                    f"applicability requires declared model {model!r}",
                )
            )
    if (
        required_contribution_domain is not None
        and observation.contribution_domain != required_contribution_domain
    ):
        findings.append(
            Finding(
                RULE_ATTRIBUTION,
                path,
                start + 1,
                "observed contribution domain does not match explicit applicability",
            )
        )
    if required_tools is not None and observation.tools != required_tools:
        findings.append(
            Finding(
                RULE_ATTRIBUTION,
                path,
                start + 1,
                "observed OpenAI tool set does not match explicit applicability",
            )
        )
    return findings


def normalize_attribution_source(
    text: str,
    path: str = "<memory>",
    *,
    current_year: int = CURRENT_YEAR,
    original_start_year: int | None = None,
    contribution_domain: str | None = None,
    attribution_tools: str | None = None,
    required_models: tuple[str, ...] = (),
) -> str:
    """Return one source with its single authoritative header normalized."""

    invalid_models = [
        model for model in required_models if not valid_required_model(model)
    ]
    if invalid_models:
        raise ValueError(f"invalid required model token(s): {invalid_models}")

    block = attribution_block(text)
    if block is None:
        if contribution_domain not in {"development", "documentation", "both"}:
            raise ValueError(
                "missing attribution requires an explicit contribution domain"
            )
        if not required_models:
            raise ValueError(
                "missing attribution requires at least one explicit model"
            )
        if attribution_tools not in {"chatgpt", "codex", "both"}:
            raise ValueError(
                "missing attribution requires an explicit OpenAI tool set"
            )
        if not all(EXPLICIT_MODEL.fullmatch(model) for model in required_models):
            raise ValueError(
                "missing attribution requires explicit model identifiers"
            )
        activity = {
            "development": "development",
            "documentation": "documentation",
            "both": "development and documentation",
        }[contribution_domain]
        tool_name, verb = {
            "chatgpt": ("OpenAI ChatGPT", "was"),
            "codex": ("OpenAI Codex", "was"),
            "both": ("OpenAI ChatGPT and Codex", "were"),
        }[attribution_tools]
        rendered = f"# {tool_name} ({', '.join(required_models)}) {verb} used in {activity}."
    else:
        _, _, rendered = block
        observation = parse_attribution(rendered)
        if observation is None:
            raise ValueError("existing attribution is not an approved OpenAI form")
        if (
            contribution_domain is not None
            and observation.contribution_domain != contribution_domain
        ):
            raise ValueError(
                "existing attribution contribution domain conflicts with applicability"
            )
        if attribution_tools is not None and observation.tools != attribution_tools:
            raise ValueError(
                "existing attribution tool set conflicts with applicability"
            )

    observation = parse_attribution(rendered)
    assert observation is not None
    declared_counts = Counter(observation.models)
    for model in required_models:
        if declared_counts[model] == 0:
            if not EXPLICIT_MODEL.fullmatch(model):
                raise ValueError(
                    "cannot add a missing generic-series requirement "
                    "without an explicit attribution-style decision"
                )
            normalized = " ".join(
                line.removeprefix("#").strip()
                for line in rendered.splitlines()
            )
            match = ATTRIBUTION_FORM.fullmatch(normalized)
            if match is None:
                raise ValueError("recognized attribution has no model list")
            declaration = parse_model_declaration(match.group("models"))
            assert declaration is not None
            models = match.group("models")
            if declaration.style == "explicit_model_list" or "; most recent: " in models:
                replacement = f"{models}, {model}"
            else:
                replacement = f"{models}; most recent: {model}"
            rendered = (
                normalized[: match.start("models")]
                + replacement
                + normalized[match.end("models") :]
            )
            declared_counts[model] += 1

    attribution_text = " ".join(
        line.removeprefix("#").strip() for line in rendered.splitlines()
    )
    attribution_lines = [
        f"# {line}"
        for line in textwrap.wrap(
            attribution_text,
            width=77,
            break_long_words=False,
            break_on_hyphens=False,
        )
    ]

    lines = text.splitlines()
    basename = expected_basename(path, lines) or "source"
    copyright_match = next(
        (COPYRIGHT.fullmatch(line) for line in lines if COPYRIGHT.fullmatch(line)),
        None,
    )
    start_year = (
        int(copyright_match.group("start"))
        if copyright_match
        else original_start_year or current_year
    )
    years = str(start_year) if start_year == current_year else f"{start_year}-{current_year}"
    shebang = lines[0] if lines and lines[0].startswith("#!") else (
        "#!/usr/bin/env python3" if basename.endswith(".py") else "#!/usr/bin/env bash"
    )
    license_index = next((index for index, line in enumerate(lines) if line == LICENSE), None)
    if license_index is not None:
        remainder = lines[license_index + 1 :]
    else:
        remainder = lines[1:] if lines and lines[0].startswith("#!") else lines
    while remainder and not remainder[0].strip():
        remainder.pop(0)
    if remainder and remainder[0] == ENCODING:
        remainder.pop(0)
        while remainder and remainder[0] == "#":
            remainder.pop(0)
        while remainder and not remainder[0].strip():
            remainder.pop(0)
    header = [
        shebang,
        ENCODING,
        "#",
        f"# Script: {basename}",
        "#",
        f"# Copyright {years} by Kris Alavattam",
        EMAIL,
        "#",
        *attribution_lines,
        "#",
        LICENSE,
        "",
        "",
    ]
    return "\n".join(header + remainder) + ("\n" if text.endswith("\n") else "")


def repository_start_year(root: Path, path: str, current_year: int) -> int:
    """Return the preserved header year or earliest tracked commit year."""

    target = root / path
    text = target.read_text(encoding="utf-8")
    match = next(
        (COPYRIGHT.fullmatch(line) for line in text.splitlines() if COPYRIGHT.fullmatch(line)),
        None,
    )
    if match:
        return int(match.group("start"))
    history = run_git(
        root,
        ["log", "--follow", "--format=%ad", "--date=format:%Y", "--", path],
        check=False,
    ).stdout.splitlines()
    return int(history[-1]) if history else current_year


def attribution_paths(root: Path) -> list[str]:
    """Return every current-diff shell and Python source path."""

    tracked_py = run_git(
        root,
        ["diff", "--name-only", "--diff-filter=ACMR", "HEAD", "--", "*.py"],
    ).stdout.splitlines()
    untracked_py = run_git(
        root,
        ["ls-files", "--others", "--exclude-standard", "--", "*.py"],
    ).stdout.splitlines()
    return sorted(set(shell_paths(root)) | set(tracked_py) | set(untracked_py))


def scan_repository(
    root: Path,
    paths: list[str] | None = None,
    *,
    current_year: int = CURRENT_YEAR,
    required_models: tuple[str, ...] = (),
    requirements: dict[str, AttributionRequirement] | None = None,
) -> list[Finding]:
    """Check current-diff source files or explicit source paths."""

    root = root.resolve()
    selected = paths if paths is not None else attribution_paths(root)
    return [
        finding
        for path in sorted(set(selected))
        for finding in check_attribution_source(
            (root / path).read_text(encoding="utf-8"),
            path,
            current_year=current_year,
            required_models=(
                requirements[path].required_models
                if requirements is not None and path in requirements
                else required_models
            ),
            required_contribution_domain=(
                requirements[path].contribution_domain
                if requirements is not None and path in requirements
                else None
            ),
            required_tools=(
                requirements[path].tools
                if requirements is not None and path in requirements
                else None
            ),
        )
        if (root / path).is_file()
    ]


def source_header_inventory(
    root: Path,
    paths: list[str],
    *,
    current_year: int = CURRENT_YEAR,
    required_models: tuple[str, ...] = (),
    requirements: dict[str, AttributionRequirement] | None = None,
) -> list[dict[str, object]]:
    """Return deterministic per-file header facts for enforced sources."""

    rows: list[dict[str, object]] = []
    for path in sorted(set(paths)):
        target = root / path
        if not target.is_file():
            continue
        lines = target.read_text(encoding="utf-8").splitlines()
        if not lines or not lines[0].startswith("#!"):
            continue
        copyright_match = next(
            (COPYRIGHT.fullmatch(line) for line in lines if COPYRIGHT.fullmatch(line)),
            None,
        )
        block = attribution_block("\n".join(lines))
        rendered = block[2] if block else ""
        prior = run_git(root, ["show", f"HEAD:{path}"], check=False)
        prior_match = next(
            (
                COPYRIGHT.fullmatch(line)
                for line in prior.stdout.splitlines()
                if COPYRIGHT.fullmatch(line)
            ),
            None,
        )
        prior_block = attribution_block(prior.stdout) if prior.returncode == 0 else None
        prior_rendered = prior_block[2] if prior_block else ""
        prior_observation = (
            parse_attribution(prior_rendered) if prior_rendered else None
        )
        observation = parse_attribution(rendered) if rendered else None
        prior_style = (
            prior_observation.attribution_style if prior_observation else None
        )
        current_style = observation.attribution_style if observation else None
        observed_counts = Counter(observation.models if observation else ())
        observed_models = sorted(observed_counts)
        models = (
            requirements[path].required_models
            if requirements is not None and path in requirements
            else required_models
        )
        if prior.returncode:
            disposition = "new_source_header"
        elif prior_match is None:
            disposition = "added_header_preserving_repository_start_year"
        elif int(prior_match.group("end") or prior_match.group("start")) != current_year:
            disposition = "extended_copyright_through_current_year"
        else:
            disposition = "normalized_existing_current_year_header"
        rows.append(
            {
                "path": path,
                "shebang": lines[0],
                "script_basename": Path(path).name,
                "copyright_start": (
                    int(copyright_match.group("start")) if copyright_match else None
                ),
                "copyright_end": (
                    int(copyright_match.group("end") or copyright_match.group("start"))
                    if copyright_match
                    else None
                ),
                "prior_copyright_start": (
                    int(prior_match.group("start")) if prior_match else None
                ),
                "prior_copyright_end": (
                    int(prior_match.group("end") or prior_match.group("start"))
                    if prior_match
                    else None
                ),
                "change_disposition": disposition,
                "attribution_style": current_style,
                "prior_attribution_style": prior_style,
                "attribution_style_preserved": (
                    prior_style is None or prior_style == current_style
                ),
                "observed_models": observed_models,
                "required_models": list(models),
                "required_model_counts": {
                    model: observed_counts[model] for model in models
                },
                "required_models_agree": all(
                    observed_counts[model] == 1 for model in models
                ),
                "observed_contribution_domain": (
                    observation.contribution_domain if observation else None
                ),
                "required_contribution_domain": (
                    requirements[path].contribution_domain
                    if requirements is not None and path in requirements
                    else None
                ),
                "contribution_domain_agrees": (
                    requirements is None
                    or path not in requirements
                    or (
                        observation is not None
                        and observation.contribution_domain
                        == requirements[path].contribution_domain
                    )
                ),
                "observed_tools": observation.tools if observation else None,
                "required_tools": (
                    requirements[path].tools
                    if requirements is not None and path in requirements
                    else None
                ),
                "tools_agree": (
                    requirements is None
                    or path not in requirements
                    or (
                        observation is not None
                        and observation.tools == requirements[path].tools
                    )
                ),
            }
        )
    return rows


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--current-year", type=int)
    parser.add_argument("--fix", action="store_true")
    parser.add_argument(
        "--contribution-domain",
        choices=("development", "documentation", "both"),
    )
    parser.add_argument(
        "--attribution-tools",
        choices=("chatgpt", "codex", "both"),
    )
    parser.add_argument(
        "--required-model",
        action="append",
        default=[],
        metavar="MODEL",
    )
    parser.add_argument("--applicability-manifest", type=Path)
    parser.add_argument(
        "--rule",
        action="append",
        choices=(
            RULE_HEADER,
            RULE_BASENAME,
            RULE_WIDTH,
            RULE_YEAR,
            RULE_ATTRIBUTION,
            RULE_ATTRIBUTION_UNIQUE,
        ),
    )
    parser.add_argument("--inventory-output", type=Path)
    parser.add_argument("paths", nargs="*")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """Report bounded attribution findings."""

    args = parse_args(argv)
    root = args.root.resolve()
    requirements: dict[str, AttributionRequirement] | None = None
    manifest_year: int | None = None
    if args.applicability_manifest is not None:
        if args.required_model:
            raise SystemExit(
                "--required-model and --applicability-manifest are mutually exclusive"
            )
        try:
            manifest_year, requirements = load_applicability_manifest(
                root,
                args.applicability_manifest.resolve(),
            )
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            raise SystemExit(f"invalid applicability manifest: {exc}") from exc
        if args.current_year is not None and args.current_year != manifest_year:
            raise SystemExit(
                "--current-year conflicts with applicability manifest assessment_year"
            )
    current_year = args.current_year or manifest_year
    if current_year is None:
        raise SystemExit(
            "an explicit --current-year or applicability-manifest assessment_year is required"
        )
    if args.fix and requirements is None:
        raise SystemExit("--fix requires --applicability-manifest")
    if requirements is not None:
        selected = args.paths or sorted(requirements)
        unknown = sorted(set(selected) - set(requirements))
        if unknown:
            raise SystemExit(
                f"selected paths are absent from applicability manifest: {unknown}"
            )
    else:
        selected = args.paths or attribution_paths(root)
    required_models = tuple(dict.fromkeys(args.required_model))
    invalid_models = [
        model for model in required_models if not valid_required_model(model)
    ]
    if invalid_models:
        raise SystemExit(f"invalid required model token(s): {invalid_models}")
    fix_errors = False
    if args.fix:
        for path in selected:
            target = root / path
            if not target.is_file():
                continue
            text = target.read_text(encoding="utf-8")
            if not text.startswith("#!"):
                continue
            try:
                normalized = normalize_attribution_source(
                    text,
                    path,
                    current_year=current_year,
                    original_start_year=repository_start_year(
                        root,
                        path,
                        current_year,
                    ),
                    contribution_domain=(
                        requirements[path].contribution_domain
                        if requirements is not None
                        else args.contribution_domain
                    ),
                    attribution_tools=(
                        requirements[path].tools
                        if requirements is not None
                        else args.attribution_tools
                    ),
                    required_models=(
                        requirements[path].required_models
                        if requirements is not None
                        else required_models
                    ),
                )
            except ValueError as exc:
                fix_errors = True
                print(f"{RULE_ATTRIBUTION}: {path}:1: {exc}", file=sys.stderr)
                continue
            if normalized != text:
                target.write_text(normalized, encoding="utf-8")
    findings = scan_repository(
        root,
        selected,
        current_year=current_year,
        required_models=required_models,
        requirements=requirements,
    )
    if args.rule:
        selected_rules = set(args.rule)
        findings = [
            finding for finding in findings if finding.rule_id in selected_rules
        ]
    for finding in findings:
        print(finding.format())
    if args.inventory_output:
        inventory = source_header_inventory(
            root,
            selected,
            current_year=current_year,
            required_models=required_models,
            requirements=requirements,
        )
        args.inventory_output.parent.mkdir(parents=True, exist_ok=True)
        args.inventory_output.write_text(
            json.dumps(inventory, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    if findings or fix_errors:
        print(f"Source header: {len(findings)} violation(s)")
        return 1
    print("Source header: pass (ordered header and declared attribution)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
