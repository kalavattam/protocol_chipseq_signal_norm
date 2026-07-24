#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: shell_help_pilot.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


"""Extract bounded download-fastqs shell/help pilot facts from a snapshot."""

from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Any

SECTION_NAMES = {"Usage", "Parameters", "Returns", "Notes", "References", "See Also", "Examples"}
PLACEHOLDERS = {
    "flag", "str", "bool", "int", "flt", "num", "path", "file", "dir", "csv",
    "mode", "method", "format", "engine", "layout", "equation", "aligner",
    "algorithm", "choice", "spec", "time", "size",
}


def line_number(text: str, needle: str) -> int | None:
    """Return the first one-based line containing 'needle'."""

    for number, line in enumerate(text.splitlines(), 1):
        if needle in line:
            return number
    return None


def help_body(text: str) -> str | None:
    """Return the literal EOM help heredoc body when its narrow shape is known."""

    match = re.search(
        r"function\s+help_[A-Za-z0-9_]+\(\)\s*\{\n"
        r"(?:[ \t]*(?:#.*)?\n)*[ \t]*cat\s+<<\s*EOM\n"
        r"(.*?)\nEOM\n\}",
        text,
        re.DOTALL,
    )
    return match.group(1) if match else None


def add_fact(rows: list[dict[str, object]], topic: str, path: str, value: Any, certainty: str = "certain", line: int | None = None) -> None:
    """Append one normalized fact without treating semantic uncertainty as failure."""

    rows.append({"topic": topic, "path": path, "line": line, "value": value, "certainty": certainty})


def function_body(text: str, name: str) -> str | None:
    """Return a narrow top-level Bash function body when its delimiters are clear."""

    match = re.search(rf"^function\s+{re.escape(name)}\(\)\s*\{{\n(.*?)(?=^\}}\s*$)", text, re.MULTILINE | re.DOTALL)
    return match.group(1) if match else None


def complete_assignment(value: str) -> bool:
    """Reject dynamic, partial, and compound assignment expressions."""

    return bool(value) and "$(" not in value and value not in {"(", "{"} and not value.endswith("\\")


def assignment_facts(path: str, text: str, facts: list[dict[str, object]]) -> None:
    """Extract only interface assignments with an accurate category."""

    scopes = (
        ("init_args_hardcoded", "hardcoded_interface_value"),
        ("init_arg_defs", "declared_argument_default"),
        ("init_defs", "declared_argument_default"),
    )
    rows: list[dict[str, object]] = []
    for name, category in scopes:
        body = function_body(text, name)
        if body is None:
            continue
        for number, line in enumerate(body.splitlines(), 1):
            match = re.match(r"\s*([A-Za-z_][A-Za-z0-9_]*)=(\"[^\"]*\"|'[^']*'|true|false|[0-9]+)\s*$", line)
            if match and complete_assignment(match.group(2)):
                rows.append({"name": match.group(1), "value": match.group(2), "assignment_kind": category, "function": name, "function_line": number})
    parse = function_body(text, "parse_args")
    if parse:
        for number, line in enumerate(parse.splitlines(), 1):
            match = re.match(r"\s*([A-Za-z_][A-Za-z0-9_]*)=\"\$\{2\}\"\s*$", line)
            if match:
                rows.append({"name": match.group(1), "value": "${2}", "assignment_kind": "parsed_argument_assignment", "function": "parse_args", "function_line": number})
    if rows:
        add_fact(facts, "interface_assignments", path, rows)


def parser_alias_facts(path: str, text: str, facts: list[dict[str, object]]) -> None:
    """Extract policy-neutral parser case-arm observations."""

    records: list[dict[str, object]] = []
    for function in ("resolve_dir_scr", "parse_args"):
        body = function_body(text, function)
        if body is None:
            continue
        for number, line in enumerate(body.splitlines(), 1):
            match = re.match(r"\s*(-[^)]*)\)\s*$", line)
            if not match:
                continue
            pattern = match.group(1).strip()
            records.append({
                "pattern": pattern,
                "observation_kind": "fallback_case_arm" if pattern == "-*" else "option_case_arm",
                "function": function,
                "line": number,
            })
    if records:
        add_fact(facts, "parser_alias_observations", path, records)
    main = function_body(text, "main") or ""
    startup = []
    if "--h[e]?lp" in main:
        startup = ["-h", "--help", "--hlp"]
    if startup:
        add_fact(facts, "startup_help_alias_observations", path, startup)


def expand_static_alias_pattern(pattern: str) -> list[str] | None:
    """Expand the narrow literal and underscore/hyphen parser syntax."""

    aliases: list[str] = []
    for token in pattern.split("|"):
        token = token.strip()
        if not token.startswith("-") or (any(character in token for character in "*?[]") and "[_-]" not in token):
            return None
        if "[_-]" in token:
            aliases.extend((token.replace("[_-]", "_"), token.replace("[_-]", "-")))
        else:
            aliases.append(token)
    return aliases


def resolve_alias_facts(
    facts: list[dict[str, object]],
    hidden_aliases: dict[str, set[str]] | None = None,
) -> None:
    """Append policy resolutions without changing raw parser observations."""

    declared_hidden = hidden_aliases or {}
    observed: dict[str, dict[str, list[dict[str, object]]]] = defaultdict(lambda: defaultdict(list))
    fallbacks: dict[str, list[dict[str, object]]] = defaultdict(list)
    for item in facts:
        if item["topic"] != "parser_alias_observations":
            continue
        path = str(item["path"])
        for record in item["value"]:
            if record["observation_kind"] == "fallback_case_arm":
                fallbacks[path].append(record)
                continue
            aliases = expand_static_alias_pattern(str(record["pattern"]))
            if aliases is None:
                continue
            for alias in aliases:
                observed[path][alias].append({"function": record["function"], "line": record["line"]})
    for item in facts:
        if item["topic"] == "startup_help_alias_observations":
            for alias in item["value"]:
                observed[str(item["path"])][str(alias)].append({"function": "main", "line": item.get("line")})
    for path in sorted(set(observed) | set(fallbacks)):
        aliases = set(observed[path])
        resolved: list[dict[str, object]] = []
        for alias in sorted(aliases):
            underscore = "--" + alias[2:].replace("-", "_") if alias.startswith("--") else alias
            if alias in declared_hidden.get(path, set()):
                visibility = "hidden_declared_compatibility"
                retention = "indefinite"
            elif alias == "--hlp":
                visibility = "hidden_legacy_compatibility"
                retention = "indefinite"
            elif alias.startswith("--") and "-" in alias[2:] and underscore in aliases:
                visibility = "hidden_systematic_compatibility"
                retention = "indefinite"
            else:
                visibility = "public"
                retention = "canonical_or_supported"
            resolved.append({"alias": alias, "visibility": visibility, "retention": retention, "locations": observed[path][alias]})
        if fallbacks[path]:
            resolved.append({"alias": "-*", "visibility": "unsupported_fallback", "retention": "not_an_option", "locations": [{"function": row["function"], "line": row["line"]} for row in fallbacks[path]]})
        add_fact(facts, "resolved_alias_classifications", path, resolved)


def documented_alias_facts(path: str, text: str, facts: list[dict[str, object]]) -> None:
    """Extract documented aliases from public help parameter rows only."""

    body = help_body(text)
    if body is None:
        return
    aliases = [
        token
        for line in body.splitlines()
        if " : " in line
        for token in re.findall(
            r"--?[A-Za-z0-9][A-Za-z0-9_-]*",
            line.split(" : ", 1)[0],
        )
    ]
    add_fact(facts, "documented_aliases", path, sorted(set(aliases)))
    duplicates = sorted(
        alias for alias in set(aliases) if aliases.count(alias) > 1
    )
    if duplicates:
        add_fact(facts, "duplicate_documented_aliases", path, duplicates)


def validate_interfaces(value: object) -> list[dict[str, object]]:
    """Validate declarative wrapper, documentation-source, and removal policy."""

    if not isinstance(value, list):
        raise ValueError("interface configuration must be a list")
    interfaces: list[dict[str, object]] = []
    owners: set[str] = set()
    for item in value:
        if not isinstance(item, dict) or set(item) != {"path", "documentation_source", "removed_aliases"}:
            raise ValueError("interface entries need path, documentation_source, and removed_aliases")
        path, source, removed = item["path"], item["documentation_source"], item["removed_aliases"]
        if not isinstance(path, str) or not path or not isinstance(source, str) or not source:
            raise ValueError("interface paths must be nonempty strings")
        if path in owners:
            raise ValueError(f"duplicate interface path: {path}")
        if not isinstance(removed, list) or not all(isinstance(alias, str) and re.fullmatch(r"--?[A-Za-z][A-Za-z0-9_-]*", alias) for alias in removed):
            raise ValueError("removed aliases must preserve exact option spelling")
        if len(removed) != len(set(removed)):
            raise ValueError("removed aliases must be unique")
        owners.add(path)
        interfaces.append({"path": path, "documentation_source": source, "removed_aliases": list(removed)})
    return interfaces


def alias_documentation_findings(
    facts: list[dict[str, object]],
    interfaces: list[dict[str, object]],
    findings: list[dict[str, object]],
) -> None:
    """Compare exact accepted spellings with each declaratively bound help source."""

    resolved = {
        str(item["path"]): {str(entry["alias"]): entry for entry in item["value"]}
        for item in facts
        if item["topic"] == "resolved_alias_classifications"
    }
    documented = {
        str(item["path"]): {str(alias) for alias in item["value"]}
        for item in facts
        if item["topic"] == "documented_aliases"
    }
    duplicates = {
        str(item["path"]): {str(alias) for alias in item["value"]}
        for item in facts
        if item["topic"] == "duplicate_documented_aliases"
    }
    configured = {str(item["path"]) for item in interfaces}
    unbound = set(resolved) - configured
    if unbound:
        raise ValueError(f"accepted aliases lack interface documentation binding: {sorted(unbound)}")
    for interface in interfaces:
        path = str(interface["path"])
        source = str(interface["documentation_source"])
        if source not in documented:
            raise ValueError(f"configured documentation source was not extracted: {source}")
        aliases = resolved.get(path, {})
        public = {alias for alias, entry in aliases.items() if entry["visibility"] == "public"}
        docs = documented[source]
        removed = {str(alias) for alias in interface["removed_aliases"]}
        add_fact(facts, "documentation_source_associations", path, {"documentation_source": source})
        add_fact(
            facts,
            "removed_alias_classifications",
            path,
            [{"alias": alias, "acceptance": "rejected", "retention": "removed"} for alias in sorted(removed)],
        )
        for alias in sorted(duplicates.get(source, set())):
            findings.append({"path": path, "line": None, "topic": "aliases", "message": f"documented alias is duplicated: {alias}", "evidence": alias})
        for alias in sorted(public - docs):
            findings.append({"path": path, "line": None, "topic": "aliases", "message": f"accepted public alias is undocumented: {alias}", "evidence": alias})
        for alias in sorted((docs - public) - removed):
            findings.append({"path": path, "line": None, "topic": "aliases", "message": f"documented alias is not accepted by parser extraction: {alias}", "evidence": alias})
        for alias in sorted(removed & (set(aliases) | docs)):
            findings.append({"path": path, "line": None, "topic": "aliases", "message": f"removed alias is accepted or documented: {alias}", "evidence": alias})


def examples_in_help(body: str) -> list[tuple[str, tuple[str, ...]]]:
    """Return numbered example labels and normalized command shapes."""

    blocks: list[tuple[str, tuple[str, ...]]] = []
    lines = body.splitlines()
    index = 0
    while index < len(lines):
        if re.match(r"\s*\d+\.\s+", lines[index]):
            label = lines[index].strip()
            if index + 1 < len(lines) and "'''bash" in lines[index + 1]:
                index += 2
                command: list[str] = []
                while index < len(lines) and "'''" not in lines[index]:
                    command.append(lines[index].strip())
                    index += 1
                shape = tuple(re.sub(r'"[^\"]*"|\$\{[^}]+\}', 'VALUE', line) for line in command if line)
                blocks.append((label, shape))
        index += 1
    return blocks


def line_length_facts(path: str, text: str, facts: list[dict[str, object]]) -> None:
    """Emit actionable candidate records while leaving classification unresolved."""

    candidates: list[dict[str, object]] = []
    lines = text.splitlines()
    heredoc_kind: dict[int, str] = {}
    active: str | None = None
    for number, line in enumerate(lines, 1):
        if re.search(r"<<\s*EOM", line):
            active = "user_facing_help_heredoc" if "function help_" in "\n".join(lines[max(0, number - 4):number]) else "local_maintainer_help_heredoc"
            continue
        if active and line == "EOM":
            active = None
            continue
        if active:
            heredoc_kind[number] = active
    for number, line in enumerate(lines, 1):
        if len(line) <= 80:
            continue
        cues = []
        location = heredoc_kind.get(number, "ordinary_shell_code")
        relationship = "independent_pilot_only_candidate"
        if location.endswith("heredoc"):
            relationship = "excluded_by_documented_source_checker_policy"
        if "http://" in line or "https://" in line:
            cues.append("url")
            location = "url_or_literal_exemption"
            relationship = "excluded_by_documented_source_checker_policy"
        if re.search(r"\[[_^\-]", line):
            cues.append("parser_pattern")
            location = "parser_pattern"
            relationship = "excluded_by_documented_source_checker_policy"
        if line.lstrip().startswith("#"):
            cues.append("comment")
            location = "comment"
        excerpt = "\n".join(f"{idx + 1}: {lines[idx]}" for idx in range(max(0, number - 2), min(len(lines), number + 1)))
        candidates.append({"line": number, "length": len(line), "excerpt": excerpt, "exception_cues": cues, "location_kind": location, "source_checker_relationship": relationship, "classification_status": "semantic/manual"})
    add_fact(facts, "line_length_candidates", path, candidates)


def supporting_alignment(path: str, text: str, registered: bool, facts: list[dict[str, object]], evidence: dict[str, object] | None = None) -> None:
    """Summarize public-interface evidence from one supporting smoke test."""

    if evidence is not None:
        groups = evidence.get("source_assertion_groups") if isinstance(evidence, dict) else None
        if not isinstance(groups, dict):
            raise ValueError("supporting controlled-smoke evidence lacks source_assertion_groups")
        coverage = {
            name: {"required_markers": markers, "present": all(isinstance(marker, str) and marker in text for marker in markers)}
            for name, markers in groups.items()
        }
        add_fact(facts, "controlled_smoke_source_coverage", path, {
            "registration_status": "registered" if registered else "not_registered",
            "evidence_kind": evidence.get("kind"),
            "test_identifier": evidence.get("test_identifier"),
            "coverage_groups": coverage,
            "execution_status": "not_evaluated_by_source_inspection",
        })
        return
    if not path.startswith((
        "tests/integration/local/download_fastqs/test_execute_download_fastqs_",
        "tests/integration/parallel/download_fastqs/test_execute_download_fastqs_",
    )):
        return
    name = next((line.split("=", 1)[1].strip().strip('"') for line in text.splitlines() if line.startswith("TEST_NAME=")), "unknown")
    lines = text.splitlines()
    start = next((index for index, line in enumerate(lines) if "execute_download_fastqs.sh" in line), None)
    command_lines: list[str] = []
    if start is not None:
        for line in lines[start:start + 24]:
            command_lines.append(line)
            if line.strip() == "then":
                break
    options = sorted(set(re.findall(r"--[A-Za-z][A-Za-z0-9_-]*", "\n".join(command_lines))))
    required = {"--fil_in", "--dir_out", "--dir_sym"}
    modes = "parallel" if "parallel" in path else "mixed" if "mixed" in path else "paired_end" if "pe_" in path else "single_end"
    direct_submit_evidence = [
        {"line": number, "source": line.strip()}
        for number, line in enumerate(lines, 1)
        if "submit_download_fastqs.sh" in line and re.search(r'"\$\{TEST_BASH\}"\s+"\$\{ROOT_REPO\}/bin/submit_download_fastqs\.sh"', line)
    ]
    direct_submit = bool(direct_submit_evidence)
    gap = "dynamic invocation boundary" if start is None else f"Direct submit-wrapper invocation at line {direct_submit_evidence[0]['line']}." if direct_submit else "No direct submit-wrapper assertion."
    add_fact(facts, "supporting_test_alignment", path, {"registration_status": "registered" if registered else "not_registered", "tested_mode": modes, "test_name": name, "public_options_invoked": options, "required_arguments_represented": sorted(required & set(options)), "dry_run_tested": "--dry_run" in options, "submit_wrapper_directly_invoked": direct_submit, "direct_submit_evidence": direct_submit_evidence, "relevant_interface_coverage": "execute wrapper invocation" if start is not None else "uncertain", "apparent_gap_or_uncertainty": gap})


def validate_command_registry(value: object) -> tuple[set[str], dict[str, set[str]]]:
    """Validate and index exact callable and conceptual registry spellings."""

    if not isinstance(value, dict) or value.get("schema_version") != 1 or not isinstance(value.get("commands"), list):
        raise ValueError("command registry must be a schema_version 1 object with commands")
    callables: set[str] = set()
    concepts: dict[str, set[str]] = defaultdict(set)
    for index, entry in enumerate(value["commands"]):
        if not isinstance(entry, dict) or set(entry) != {"callable", "conceptual_names"}:
            raise ValueError(f"command registry entry {index} must contain only callable and conceptual_names")
        callable_name = entry["callable"]
        conceptual_names = entry["conceptual_names"]
        if not isinstance(callable_name, str) or not callable_name or callable_name != callable_name.strip() or re.search(r"\s", callable_name):
            raise ValueError(f"command registry entry {index} has an invalid callable")
        if callable_name in callables:
            raise ValueError(f"duplicate command registry callable: {callable_name}")
        if not isinstance(conceptual_names, list) or not all(isinstance(name, str) and name and name == name.strip() for name in conceptual_names):
            raise ValueError(f"command registry entry {index} has invalid conceptual_names")
        if len(conceptual_names) != len(set(conceptual_names)):
            raise ValueError(f"command registry entry {index} repeats a conceptual name")
        callables.add(callable_name)
        for name in conceptual_names:
            concepts[name].add(callable_name)
    return callables, concepts


def external_program_observations(path: str, text: str) -> list[dict[str, object]]:
    """Return references from a rendered help External programs block."""

    body = help_body(text)
    if body is None:
        return []
    lines = body.splitlines()
    observations: list[dict[str, object]] = []
    active_indent: int | None = None
    for number, line in enumerate(lines, 1):
        if line.strip() == "External programs:":
            active_indent = len(line) - len(line.lstrip())
            continue
        if active_indent is None:
            continue
        indent = len(line) - len(line.lstrip())
        if line.strip() and indent <= active_indent:
            active_indent = None
            continue
        match = re.match(r"\s*-\s+(.+?)\s*$", line)
        if not match:
            continue
        reference = match.group(1).split(",", 1)[0].strip()
        candidate = re.sub(r"\s+(?:>=|<=|==|>|<)\s*[^ ]+$", "", reference).strip()
        observations.append({"path": path, "line": number, "scope": "external_programs", "reference": reference, "candidate": candidate})
    return observations


def runtime_version_observations(path: str, text: str, known_names: set[str]) -> list[dict[str, object]]:
    """Return registered or unknown names used in version comparisons."""

    observations: list[dict[str, object]] = []
    comparison = re.compile(r"(?:>=|<=|==|>|<)\s*[0-9]+(?:\.[0-9]+)*")
    ordered = sorted(known_names, key=len, reverse=True)
    for number, line in enumerate(text.splitlines(), 1):
        match = comparison.search(line)
        if not match:
            continue
        prefix = line[:match.start()].rstrip()
        candidate = next((name for name in ordered if re.search(rf"(?<![A-Za-z0-9_]){re.escape(name)}$", prefix)), None)
        if candidate is None:
            continue
        observations.append({"path": path, "line": number, "scope": "runtime_version", "reference": line.strip(), "candidate": candidate})
    return observations


def command_reference_facts(
    path: str,
    text: str,
    callables: set[str],
    concepts: dict[str, set[str]],
    scopes: set[str],
    facts: list[dict[str, object]],
    findings: list[dict[str, object]],
    policy_questions: list[dict[str, object]],
) -> None:
    """Resolve scoped command references without normalizing unknown names."""

    observations: list[dict[str, object]] = []
    if "external_programs" in scopes:
        observations.extend(external_program_observations(path, text))
    if "runtime_version" in scopes:
        version_rows = runtime_version_observations(path, text, callables | set(concepts))
        external_versions = {
            (str(row["candidate"]), str(row["reference"]))
            for row in observations
            if row["scope"] == "external_programs" and re.search(r"(?:>=|<=|==|>|<)\s*[0-9]", str(row["reference"]))
        }
        observations.extend(
            row for row in version_rows
            if not any(candidate == str(row["candidate"]) and reference in str(row["reference"]) for candidate, reference in external_versions)
        )
    unique = {(row["line"], row["scope"], row["candidate"], row["reference"]): row for row in observations}
    observations = [unique[key] for key in sorted(unique)]
    if observations:
        add_fact(facts, "command_reference_observations", path, [{key: value for key, value in row.items() if key != "path"} for row in observations])
    resolutions: list[dict[str, object]] = []
    for row in observations:
        candidate = str(row["candidate"])
        if candidate in callables:
            status = "exact_callable"
            expected: str | None = candidate
        elif len(concepts.get(candidate, set())) == 1:
            status = "conceptual_name_in_callable_scope"
            expected = next(iter(concepts[candidate]))
            findings.append({"path": path, "line": row["line"], "topic": "runtime_requirements", "message": f"callable reference must use exact registry spelling '{expected}' instead of conceptual name '{candidate}'", "evidence": row["reference"]})
        elif candidate in concepts:
            status = "ambiguous_conceptual_name"
            expected = None
            policy_questions.append({"topic": "runtime_requirements", "question": f"Resolve ambiguous runtime-command reference '{candidate}' without guessing among {sorted(concepts[candidate])}.", "paths": [path]})
        else:
            status = "unknown_command_name"
            expected = None
            policy_questions.append({"topic": "runtime_requirements", "question": f"Resolve unknown runtime-command reference '{candidate}' without spelling normalization.", "paths": [path]})
        resolutions.append({"line": row["line"], "scope": row["scope"], "candidate": candidate, "status": status, "exact_callable": expected})
    if resolutions:
        add_fact(facts, "command_reference_resolutions", path, resolutions)


def source_style(path: str, text: str, facts: list[dict[str, object]], findings: list[dict[str, object]]) -> None:
    """Extract source-help vocabulary facts and narrow documented syntax findings."""

    body = help_body(text)
    if body is None:
        return
    lines = body.splitlines()
    headings = [line for line in lines if line in SECTION_NAMES]
    add_fact(facts, "headings", path, headings)
    examples = examples_in_help(body)
    add_fact(facts, "examples", path, {"count": len(examples), "distinct_shapes": len({shape for _, shape in examples})})
    for index, line in enumerate(lines):
        if line in SECTION_NAMES and (index + 1 >= len(lines) or not re.fullmatch(r"-{3,}", lines[index + 1])):
            findings.append({"path": path, "line": index + 1, "topic": "headings", "message": f"section '{line}' lacks an underline", "evidence": line})
    if path.endswith("help_execute_download_fastqs.sh") or path.endswith("help_submit_download_fastqs.sh"):
        if "Examples" not in headings:
            findings.append({"path": path, "line": None, "topic": "examples", "message": "top-level wrapper help is missing Examples", "evidence": "Examples"})
        elif len(examples) < 2 or len({shape for _, shape in examples}) < 2:
            findings.append({"path": path, "line": None, "topic": "examples", "message": "wrapper help requires two materially distinct examples", "evidence": str(examples)})
    parameter_rows = [line.strip() for line in lines if re.match(r"\s*(?:-\S|\d+\s+\S).*\s:\s", line)]
    add_fact(facts, "help_parameter_rows", path, parameter_rows)
    for index, line in enumerate(lines, 1):
        for placeholder in re.findall(r"<([^>]+)>", line):
            if placeholder not in PLACEHOLDERS:
                findings.append({"path": path, "line": index, "topic": "types_and_placeholders", "message": f"noncanonical placeholder <{placeholder}>", "evidence": line})
        if re.match(r"\s*-\S.*:\s*(?:enum:|csv:|spec\s*$)", line):
            findings.append({"path": path, "line": index, "topic": "types_and_placeholders", "message": "compact parameter pseudo-type", "evidence": line})
    inventories = [line for line in lines if re.match(r"\s*(?:Functions|Function scripts|Sourced functions|Shell functions|Helpers|Helper functions|Sourced helpers)\s*:?\s*$", line)]
    add_fact(facts, "user_facing_helper_inventory", path, inventories)
    if inventories:
        findings.append({"path": path, "line": line_number(body, inventories[0]), "topic": "help_ownership", "message": "user-facing helper inventory", "evidence": inventories[0]})
def wrapper_contract(path: str, text: str, facts: list[dict[str, object]], findings: list[dict[str, object]]) -> None:
    """Extract conservative execute/submit topology facts."""

    is_execute = path == "bin/execute_download_fastqs.sh"
    is_submit = path == "bin/submit_download_fastqs.sh"
    if not (is_execute or is_submit):
        return
    expected_help = f"help_{path.rsplit('/', 1)[-1].removesuffix('.sh')}"
    inline_main_help = f"function {expected_help}()" in text
    external_help_reference = expected_help in text
    add_fact(facts, "help_ownership", path, {"external_help_reference": external_help_reference, "inline_main_help": inline_main_help})
    if inline_main_help:
        findings.append({"path": path, "line": None, "topic": "help_ownership", "message": "top-level wrapper defines inline main help", "evidence": expected_help})
    if not external_help_reference:
        findings.append({"path": path, "line": None, "topic": "help_ownership", "message": "top-level wrapper does not reference external main help", "evidence": expected_help})
    if is_execute:
        physical = 'dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}"' in text
        derived_submit = 'scr_sub="${dir_scr}/submit_download_fastqs.sh"' in text
        forwarded = "--dir_scr \"${dir_scr}\"" in text
        add_fact(facts, "execute_submit_dir_scr", path, {"physical_default": physical, "derived_submit": derived_submit, "forwarded": forwarded})
        if not (physical and derived_submit and forwarded):
            findings.append({"path": path, "line": None, "topic": "dir_scr", "message": "execute-to-submit --dir_scr contract is incomplete", "evidence": "physical_default/derived_submit/forwarded"})
    if is_submit:
        bootstrap = "-ds|--dir[_-]scr" in text
        add_fact(facts, "submit_dir_scr_bootstrap", path, {"parser": bootstrap, "help_sourced_before_parse": "help_submit_download_fastqs.sh" in text})
        if not bootstrap:
            findings.append({"path": path, "line": None, "topic": "dir_scr", "message": "submit wrapper lacks -ds|--dir_scr bootstrap parser", "evidence": "parse_args"})


def interface_facts(path: str, text: str, facts: list[dict[str, object]]) -> None:
    """Collect facts that require human interpretation across the selected family."""

    if path.startswith("bin/"):
        parser_alias_facts(path, text, facts)
        assignment_facts(path, text, facts)
    documented_alias_facts(path, text, facts)
    body = help_body(text)
    if body is not None:
        headings = [line for line in body.splitlines() if line in SECTION_NAMES]
        is_wrapper_help = path.endswith("help_execute_download_fastqs.sh") or path.endswith("help_submit_download_fastqs.sh")
        is_submit_help = path.endswith("help_submit_download_fastqs.sh")
        add_fact(facts, "section_assessment", path, {
            "examples": "present" if "Examples" in headings else "required_missing" if is_wrapper_help else "not_assessed",
            "runtime_requirements": "present" if "Runtime requirements:" in body else "required_missing" if is_submit_help else "not_applicable",
        })
    stale = {name: [number for number, line in enumerate(text.splitlines(), 1) if re.search(rf"(?<![A-Za-z_]){name}(?![A-Za-z_])", line)] for name in ("err_out", "infile", "outfile")}
    add_fact(facts, "stale_name_occurrences", path, stale)
    line_length_facts(path, text, facts)


def main(argv: list[str] | None = None) -> int:
    """Emit bounded structured facts and narrow source-style findings."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", required=True, choices=("source_style", "command_references", "wrapper_contract", "interface_facts"))
    parser.add_argument("--snapshot", required=True)
    args = parser.parse_args(argv)
    snapshot = json.loads(Path(args.snapshot).read_text(encoding="utf-8"))
    if snapshot.get("schema_version") != 1 or not isinstance(snapshot.get("targets"), list):
        raise ValueError("snapshot must be a schema_version 1 object with targets")
    facts: list[dict[str, object]] = []
    findings: list[dict[str, object]] = []
    policy_questions: list[dict[str, object]] = []
    callables: set[str] = set()
    concepts: dict[str, set[str]] = {}
    reference_scopes: set[str] = set()
    if args.mode == "command_references":
        config = snapshot.get("adapter_config")
        if not isinstance(config, dict) or not isinstance(config.get("registry_path"), str) or not isinstance(config.get("reference_scopes"), list):
            raise ValueError("command_references mode requires registry_path and reference_scopes adapter configuration")
        registry_content = next((item.get("content") for item in snapshot.get("context", []) if item.get("path") == config["registry_path"]), None)
        if not isinstance(registry_content, str):
            raise ValueError("configured command registry is absent from snapshot context")
        callables, concepts = validate_command_registry(json.loads(registry_content))
        reference_scopes = {str(scope) for scope in config["reference_scopes"]}
        allowed_scopes = {"runtime_version", "external_programs"}
        if not reference_scopes or not reference_scopes <= allowed_scopes:
            raise ValueError(f"unsupported command-reference scopes: {sorted(reference_scopes - allowed_scopes)}")
    for target in snapshot["targets"]:
        path, text = target["path"], target["content"]
        if args.mode == "source_style":
            source_style(path, text, facts, findings)
        elif args.mode == "command_references":
            command_reference_facts(path, text, callables, concepts, reference_scopes, facts, findings, policy_questions)
        elif args.mode == "wrapper_contract":
            wrapper_contract(path, text, facts, findings)
        else:
            interface_facts(path, text, facts)
    if args.mode == "interface_facts":
        config = snapshot.get("adapter_config", {})
        if not isinstance(config, dict):
            raise ValueError("interface_facts adapter configuration must be an object")
        interfaces = validate_interfaces(config.get("interfaces", []))
        documentation_sources = snapshot.get("documentation_sources", [])
        if not isinstance(documentation_sources, list) or not all(isinstance(item, dict) and set(item) == {"path", "content"} and isinstance(item["path"], str) and isinstance(item["content"], str) for item in documentation_sources):
            raise ValueError("documentation_sources must contain path and content strings")
        required_sources = {str(item["documentation_source"]) for item in interfaces}
        supplied_sources = {str(item["path"]) for item in documentation_sources}
        if required_sources != supplied_sources:
            raise ValueError("configured and supplied documentation sources differ")
        existing_documentation = {str(item["path"]) for item in facts if item["topic"] == "documented_aliases"}
        for item in documentation_sources:
            if item["path"] not in existing_documentation:
                documented_alias_facts(str(item["path"]), str(item["content"]), facts)
        context = snapshot.get("context", [])
        registration = next((item["content"] for item in context if item["path"] == "tests/run_tests.sh"), "")
        for target in snapshot["targets"]:
            supporting_alignment(target["path"], target["content"], target["path"].split("/")[-1] in registration, facts, target.get("evidence"))
        resolve_alias_facts(facts)
        alias_documentation_findings(facts, interfaces, findings)
        submit_section = next((item for item in facts if item["topic"] == "section_assessment" and item["path"].endswith("help_submit_download_fastqs.sh")), None)
        if submit_section and submit_section["value"]["runtime_requirements"] != "present":
            findings.append({"path": submit_section["path"], "line": None, "topic": "runtime_requirements", "message": "submit help is missing Runtime requirements for bash >= 4.4, wget, and ln", "evidence": "Runtime requirements"})
    limitations_by_mode = {
        "source_style": [
            ("static_help_ownership", "false_negative", "Static ownership detection cannot establish broader function visibility or dynamic sourcing."),
            ("example_shape_count", "false_positive", "Example counting detects absence/count but cannot determine usefulness or semantic distinctness beyond normalized shape."),
            ("help_syntax", "false_negative", "Vocabulary and heading checks do not establish conceptual documentation quality."),
        ],
        "wrapper_contract": [
            ("static_topology", "false_negative", "Literal topology checks can miss dynamically constructed paths, options, and dispatch behavior."),
            ("behavioral_diff", "semantic_only", "Focused diff evidence cannot prove preserved runtime behavior without independent semantic review."),
        ],
        "command_references": [
            ("registry_scope", "semantic_only", "Only configured executable-reference scopes are resolved; conceptual prose remains out of scope."),
            ("unknown_or_ambiguous_command", "semantic_only", "Unknown and ambiguous command references are never normalized by guesswork."),
        ],
        "interface_facts": [
            ("alias_acceptance", "false_positive", "Parser acceptance does not by itself establish intended public visibility; dynamic aliases may be missed."),
            ("assignment_scope", "false_negative", "Assignment extraction does not establish effective runtime values, requiredness, or environment propagation."),
            ("line_and_test_static_analysis", "false_negative", "Heredoc classification and smoke invocation alignment depend on static delimiters and recognizable command boundaries."),
        ],
    }
    limitations = [{"adapter": args.mode, "limitation_area": area, "risk_type": risk, "limitation": limitation, "review_implication": "Treat this area as bounded evidence and retain semantic review where indicated."} for area, risk, limitation in limitations_by_mode[args.mode]]
    if args.mode == "interface_facts":
        for target in snapshot["targets"]:
            evidence = target.get("evidence")
            if isinstance(evidence, dict) and isinstance(evidence.get("limitation"), str):
                limitations.append({"adapter": args.mode, "limitation_area": "controlled_smoke_provenance", "risk_type": "semantic_only", "limitation": evidence["limitation"], "review_implication": "Treat the selected-environment and binary-provenance boundary as semantic review evidence."})
    print(json.dumps({"schema_version": 1, "facts": facts, "findings": findings, "policy_questions": policy_questions, "limitations": limitations}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
