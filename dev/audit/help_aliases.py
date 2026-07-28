#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: help_aliases.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


"""
Compare bounded shell parser aliases with public Parameter rows.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import re
import shlex
import sys
from pathlib import Path

from dev.audit.help_heredoc_reflow import extract_help_heredocs, shell_paths
from dev.audit.shell_help_pilot import expand_static_alias_pattern

FUNCTION_START = re.compile(
    r"^\s*(?:function\s+)?(?P<name>[A-Za-z_][A-Za-z0-9_]*)\s*\(\s*\)\s*\{",
)
OPTION_ARM = re.compile(r"^\s*(?P<pattern>-[^)]*)\)\s*$")
PARAMETER_ROW = re.compile(
    r"^(?P<indent>\s*)(?P<head>-[^:]+?)\s+:\s+(?P<type>\S.*)$",
)
UNDERLINE = re.compile(r"^-{3,}$")
COMPATIBILITY_EXEC = re.compile(
    r"exec\s+[^\n]*?/(?P<target>(?:execute|submit)_[A-Za-z0-9_]+\.sh)"
    r"[^\n]*?\"?\$@\"?",
)


@dataclasses.dataclass(frozen=True)
class AliasChunk:
    """
    One bounded logical option group from a shell case arm.
    """

    owner: str
    aliases: tuple[str, ...]


@dataclasses.dataclass(frozen=True)
class ParameterRow:
    """
    One documented logical option row.
    """

    owner: str
    line: int
    aliases: tuple[str, ...]
    indent: str
    type_name: str


@dataclasses.dataclass(frozen=True)
class Finding:
    """
    One exact parser/documentation alias mismatch.
    """

    path: str
    line: int
    owner: str
    documented: tuple[str, ...]
    expected: tuple[str, ...]
    hidden: tuple[str, ...]
    message: str

    def format(self) -> str:
        """
        Render one stable diagnostic.
        """

        location = (
            f"HELP.PARAMETER.ALIAS_SET: {self.path}:{self.line}: "
            f"owner={self.owner}"
        )
        documented = f"documented={','.join(self.documented)}"
        expected = f"expected={','.join(self.expected)}"
        hidden = f"hidden={','.join(self.hidden)}"

        return "; ".join(
            (location, self.message, documented, expected, hidden),
        )


@dataclasses.dataclass(frozen=True)
class DelegatedParserBinding:
    """
    One documented owner bound to option chunks parsed elsewhere.
    """

    documentation_path: str
    documented_owner: str
    parser_path: str
    parser_owner: str
    source_call_line: int
    fixed_mode: str
    relation: str
    confidence: str
    applicable_chunks: tuple[AliasChunk, ...]
    rejected_aliases: tuple[str, ...]
    source_evidence: tuple[str, ...]


@dataclasses.dataclass(frozen=True)
class DelegatedParserPolicy:
    """
    Stable identity metadata for one source-reviewed conditional parser.
    """

    documentation_path: str
    documented_owner: str
    parser_path: str
    parser_owner: str
    mode_argument_index: int
    fixed_mode: str
    applicable_logical_options: tuple[str, ...]
    rejected_aliases: tuple[str, ...]
    source_reference: str


FILTER_SHARED_LOGICAL_OPTIONS = (
    "--threads",
    "--fil_in",
    "--fil_out",
    "--mito",
    "--chk_chr",
    "--ref_fa",
)

DELEGATED_PARSER_POLICIES = (
    DelegatedParserPolicy(
        documentation_path="lib/bash/workflows/filter_alignment.sh",
        documented_owner="filter_alignment_sc",
        parser_path="lib/bash/workflows/filter_alignment.sh",
        parser_owner="_parse_args_filter_alignment",
        mode_argument_index=1,
        fixed_mode="sc",
        applicable_logical_options=FILTER_SHARED_LOGICAL_OPTIONS,
        rejected_aliases=("-tg", "--tg", "-mr", "--mr", "--mtr"),
        source_reference=(
            "filter_alignment_sc calls _parse_args_filter_alignment with "
            "fixed organism selector 'sc'"
        ),
    ),
    DelegatedParserPolicy(
        documentation_path="lib/bash/workflows/filter_alignment.sh",
        documented_owner="filter_alignment_sp",
        parser_path="lib/bash/workflows/filter_alignment.sh",
        parser_owner="_parse_args_filter_alignment",
        mode_argument_index=1,
        fixed_mode="sp",
        applicable_logical_options=(
            *FILTER_SHARED_LOGICAL_OPTIONS,
            "--tg",
            "--mtr",
        ),
        rejected_aliases=(),
        source_reference=(
            "filter_alignment_sp calls _parse_args_filter_alignment with "
            "fixed organism selector 'sp'"
        ),
    ),
)


def csv_short_alias_findings(
    inventory: list[dict[str, object]],
) -> list[Finding]:
    """
    Require public shorts for canonical --csv_* options to start with -c.
    """

    findings: list[Finding] = []

    for row in inventory:
        logical = str(row.get("logical_option", ""))
        if not logical.startswith("--csv_"):
            continue

        aliases = tuple(str(alias) for alias in row.get("public_aliases", []))
        shorts = tuple(
            alias
            for alias in aliases
            if alias.startswith("-") and not alias.startswith("--")
        )
        invalid = tuple(
            alias for alias in shorts if not alias.startswith("-c")
        )
        if not invalid:
            continue

        invalid_text = ",".join(invalid)
        findings.append(
            Finding(
                path=str(row.get("path", "<unknown>")),
                line=int(row.get("line", 1)),
                owner=str(row.get("owner", "<unknown>")),
                documented=aliases,
                expected=tuple(
                    alias for alias in aliases if alias not in invalid
                ),
                hidden=(),
                message=(
                    "public short aliases for canonical --csv_* options must "
                    f"begin with -c; invalid={invalid_text}"
                ),
            ),
        )

    return findings


def masked_source_lines(text: str) -> list[str]:
    """
    Mask recognized help heredocs while preserving source line numbers.
    """

    lines = text.splitlines()

    for heredoc in extract_help_heredocs(text):
        for line_number in range(heredoc.start_line, heredoc.end_line + 1):
            lines[line_number - 1] = ""

    return lines


def function_body_spans(text: str) -> dict[str, tuple[int, int, str]]:
    """
    Return one-based source spans and bodies for bounded functions.
    """

    lines = masked_source_lines(text)
    result: dict[str, tuple[int, int, str]] = {}
    index = 0

    while index < len(lines):
        match = FUNCTION_START.match(lines[index])

        if match is None:
            index += 1

            continue

        depth = lines[index].count("{") - lines[index].count("}")
        end = index + 1

        while end < len(lines) and depth > 0:
            depth += lines[end].count("{") - lines[end].count("}")
            end += 1

        if depth == 0:
            result[match.group("name")] = (
                index + 1,
                end,
                "\n".join(lines[index + 1 : end - 1]),
            )
            index = end
        else:
            index += 1

    return result


def function_bodies(text: str) -> dict[str, str]:
    """
    Return bounded function bodies while preserving nested brace blocks.
    """

    return {
        owner: body
        for owner, (_, _, body) in function_body_spans(text).items()
    }


def split_logical_aliases(aliases: list[str]) -> list[tuple[str, ...]]:
    """
    Split combined case arms at a short alias after observed long aliases.
    """

    chunks: list[list[str]] = []
    current: list[str] = []
    saw_long = False

    for alias in aliases:
        is_short = alias.startswith("-") and not alias.startswith("--")

        if is_short and saw_long:
            chunks.append(current)
            current = []
            saw_long = False

        current.append(alias)
        saw_long = saw_long or alias.startswith("--")

    if current:
        chunks.append(current)

    return [tuple(chunk) for chunk in chunks]


def _alias_chunks_from_body(owner: str, body: str) -> list[AliasChunk]:
    """
    Extract bounded alias groups from one parser-bearing source body.
    """

    chunks: list[AliasChunk] = []

    for line in body.splitlines():
        match = OPTION_ARM.match(line)
        if match is None:
            continue

        expanded = expand_static_alias_pattern(match.group("pattern"))
        if expanded is None or expanded == ["-*"]:
            continue

        for aliases in split_logical_aliases(expanded):
            chunks.append(AliasChunk(owner=owner, aliases=aliases))

    regex_groups = (
        (r"-h\|--h\[e\]\?lp", ("-h", "--hlp", "--help")),
        (r"-d\|--details", ("-d", "--details")),
        (
            r"-ah\|--all\[_-\]h\[e\]\?lp",
            ("-ah", "--all_hlp", "--all_help", "--all-hlp", "--all-help"),
        ),
    )
    observed = {chunk.aliases for chunk in chunks}

    for pattern, aliases in regex_groups:
        if re.search(pattern, body) and aliases not in observed:
            chunks.append(AliasChunk(owner=owner, aliases=aliases))

    return chunks


def alias_chunks(
    text: str,
    owners: set[str] | None = None,
) -> list[AliasChunk]:
    """
    Extract literal and underscore/hyphen option groups from functions.
    """

    chunks: list[AliasChunk] = []

    for owner, body in function_bodies(text).items():
        if owners is not None and owner not in owners:
            continue

        chunks.extend(_alias_chunks_from_body(owner, body))

    return chunks


def file_alias_chunks(text: str) -> list[AliasChunk]:
    """
    Extract aliases from a top-level parser when no parser function exists.
    """

    return _alias_chunks_from_body("<file>", text)


def parameter_rows(text: str) -> list[ParameterRow]:
    """
    Extract option rows only from recognized Parameters sections.
    """

    rows: list[ParameterRow] = []

    for heredoc in extract_help_heredocs(text):
        if heredoc.owner == "<file>":
            continue

        lines = list(heredoc.lines)
        in_parameters = False

        for index, (number, line) in enumerate(lines):
            if (
                line.strip() == "Parameters"
                and index + 1 < len(lines)
                and UNDERLINE.fullmatch(lines[index + 1][1].strip())
            ):
                in_parameters = True
                continue

            if (
                in_parameters
                and line.strip()
                and index + 1 < len(lines)
                and UNDERLINE.fullmatch(lines[index + 1][1].strip())
            ):
                in_parameters = False

            if not in_parameters:
                continue

            match = PARAMETER_ROW.match(line)
            if match is None:
                continue

            aliases = tuple(
                re.findall(
                    r"--?[A-Za-z0-9][A-Za-z0-9_-]*",
                    match.group("head"),
                ),
            )
            rows.append(
                ParameterRow(
                    owner=heredoc.owner,
                    line=number,
                    aliases=aliases,
                    indent=match.group("indent"),
                    type_name=match.group("type"),
                ),
            )

    return rows


def usage_aliases_by_owner(text: str) -> dict[str, tuple[str, ...]]:
    """
    Return exact option spellings observed in recognized Usage sections.
    """

    result: dict[str, tuple[str, ...]] = {}
    heredocs = extract_help_heredocs(text)

    for heredoc in heredocs:
        lines = list(heredoc.lines)
        aliases: list[str] = []
        in_usage = False

        for index, (_, line) in enumerate(lines):
            if (
                line.strip() == "Usage"
                and index + 1 < len(lines)
                and UNDERLINE.fullmatch(lines[index + 1][1].strip())
            ):
                in_usage = True
                continue

            if (
                in_usage
                and line.strip()
                and index + 1 < len(lines)
                and UNDERLINE.fullmatch(lines[index + 1][1].strip())
            ):
                break

            if in_usage:
                aliases.extend(
                    re.findall(r"--?[A-Za-z0-9][A-Za-z0-9_-]*", line),
                )

        result[heredoc.owner] = tuple(dict.fromkeys(aliases))

    shared = result.get("<file>", ())

    if shared:
        for heredoc in heredocs:
            if heredoc.owner == "<file>":
                continue

            if any("${usage}" in line for _, line in heredoc.lines):
                result[heredoc.owner] = shared

    return result


def hidden_alias(alias: str, aliases: tuple[str, ...]) -> bool:
    """
    Return whether settled repository convention makes an alias hidden.
    """

    if alias == "--hlp" or alias.endswith(("_hlp", "-hlp")):
        return True

    chunk_size_is_public = "--chunk_size" in aliases
    chunk_size_compatibility_alias = alias in {
        "--chnk_size",
        "--chnk-size",
    }

    if chunk_size_is_public and chunk_size_compatibility_alias:
        return True

    if alias.startswith("--") and "-" in alias[2:]:
        underscore = "--" + alias[2:].replace("-", "_")
        return underscore in aliases

    return False


def expected_aliases(
    row: ParameterRow,
    chunk: AliasChunk,
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """
    Return public parser spellings and deliberately hidden compatibility.
    """

    public: list[str] = []
    hidden: list[str] = []

    for alias in chunk.aliases:
        if hidden_alias(alias, chunk.aliases):
            hidden.append(alias)
        else:
            public.append(alias)

    return tuple(public), tuple(hidden)


def logical_option(chunk: AliasChunk) -> str:
    """
    Return the canonical public identity for one logical option chunk.
    """

    public, _ = expected_aliases(
        ParameterRow(chunk.owner, 1, (), "  ", "option"),
        chunk,
    )

    return next(
        (alias for alias in reversed(public) if alias.startswith("--")),
        public[-1] if public else "",
    )


def logical_shell_statements(
    body: str,
    body_start_line: int,
) -> list[tuple[int, str]]:
    """
    Join explicit shell continuations into line-addressable statements.
    """

    statements: list[tuple[int, str]] = []
    parts: list[str] = []
    start = body_start_line

    for offset, line in enumerate(body.splitlines()):
        number = body_start_line + offset
        stripped = line.strip()

        if not parts:
            start = number

        if stripped.endswith("\\"):
            parts.append(stripped[:-1].rstrip())

            continue

        parts.append(stripped)
        statement = " ".join(part for part in parts if part)

        if statement:
            statements.append((start, statement))

        parts = []

    if parts:
        statements.append((start, " ".join(part for part in parts if part)))

    return statements


def parser_calls(
    text: str,
    caller: str,
    parser_owner: str,
) -> list[tuple[int, tuple[str, ...], str]]:
    """
    Return direct calls from one function to one parser-bearing function.
    """

    span = function_body_spans(text).get(caller)
    if span is None:
        return []

    start_line, _, body = span
    calls: list[tuple[int, tuple[str, ...], str]] = []

    for line, statement in logical_shell_statements(body, start_line + 1):
        try:
            tokens = shlex.split(statement, comments=True, posix=True)
        except ValueError:
            tokens = statement.split()

        indexes = [
            index
            for index, token in enumerate(tokens)
            if token == parser_owner
        ]

        for index in indexes:
            arguments: list[str] = []

            for token in tokens[index + 1 :]:
                if token in {";", "&&", "||", "then", "do"}:
                    break

                arguments.append(token)

            calls.append((line, tuple(arguments), statement))

    return calls


def policy_for_binding(
    documentation_path: str,
    documented_owner: str,
    parser_path: str,
    parser_owner: str,
) -> DelegatedParserPolicy | None:
    """
    Return exact conditional-contract metadata for one stable identity.
    """

    requested_identity = (
        documentation_path,
        documented_owner,
        parser_path,
        parser_owner,
    )

    for policy in DELEGATED_PARSER_POLICIES:
        policy_identity = (
            policy.documentation_path,
            policy.documented_owner,
            policy.parser_path,
            policy.parser_owner,
        )
        if policy_identity == requested_identity:
            return policy

    return None


def _binding_for_call(
    documentation_path: str,
    documented_owner: str,
    parser_path: str,
    parser_owner: str,
    parser_chunks: list[AliasChunk],
    call: tuple[int, tuple[str, ...], str],
) -> DelegatedParserBinding:
    """
    Resolve one direct call through optional fixed-mode contract metadata.
    """

    line, arguments, statement = call
    policy = policy_for_binding(
        documentation_path,
        documented_owner,
        parser_path,
        parser_owner,
    )
    fixed_mode = ""
    applicable = tuple(parser_chunks)
    rejected: tuple[str, ...] = ()
    evidence = [f"{parser_path}:{line}: {statement}"]
    confidence = "high"

    if policy is not None:
        fixed_mode = (
            arguments[policy.mode_argument_index]
            if len(arguments) > policy.mode_argument_index
            else ""
        )
        allowed = set(policy.applicable_logical_options)
        applicable = tuple(
            chunk
            for chunk in parser_chunks
            if logical_option(chunk) in allowed
        )
        rejected = policy.rejected_aliases
        evidence.append(policy.source_reference)

        if fixed_mode != policy.fixed_mode:
            confidence = "invalid_fixed_mode"

    return DelegatedParserBinding(
        documentation_path=documentation_path,
        documented_owner=documented_owner,
        parser_path=parser_path,
        parser_owner=parser_owner,
        source_call_line=line,
        fixed_mode=fixed_mode,
        relation="direct_delegated_parser_call",
        confidence=confidence,
        applicable_chunks=applicable,
        rejected_aliases=rejected,
        source_evidence=tuple(evidence),
    )


def owning_source(root: Path, documentation_path: str) -> Path | None:
    """
    Resolve one external help module to its top-level parser source.
    """

    name = Path(documentation_path).name
    if not name.startswith("help_"):
        return None

    script_name = name.removeprefix("help_")
    candidates = (
        root / "bin" / script_name,
        root / "install" / "scripts" / script_name,
    )

    return next((path for path in candidates if path.is_file()), None)


def delegated_parser_bindings(
    root: Path,
    path: str,
    text: str,
) -> list[DelegatedParserBinding]:
    """
    Discover high-confidence parser/help ownership across function bounds.

    Parameters
    ----------
    root : Path
        Repository root containing help and parser sources.
    path : str
        Repository-relative documentation-source path.
    text : str
        Complete help source to associate with parser calls.

    Returns
    -------
    bindings : list[DelegatedParserBinding]
        Source-derived direct and centralized parser ownership bindings.
    """

    bindings: list[DelegatedParserBinding] = []
    source = owning_source(root, path)

    if source is not None:
        source_path = str(source.relative_to(root))
        source_text = source.read_text(encoding="utf-8")
        parser_owners = {
            "parse_args",
            "resolve_dir_scr",
            "check_args_light",
            "main",
        }
        chunks_by_owner: dict[str, list[AliasChunk]] = {}

        for chunk in alias_chunks(source_text, parser_owners):
            chunks_by_owner.setdefault(chunk.owner, []).append(chunk)

        source_lines = source_text.splitlines()

        for heredoc in extract_help_heredocs(text):
            if heredoc.owner == "<file>":
                continue

            call_line = next(
                (
                    number
                    for number, line in enumerate(source_lines, 1)
                    if re.search(rf"\b{re.escape(heredoc.owner)}\b", line)
                ),
                1,
            )

            for parser_owner, chunks in sorted(chunks_by_owner.items()):
                bindings.append(
                    DelegatedParserBinding(
                        documentation_path=path,
                        documented_owner=heredoc.owner,
                        parser_path=source_path,
                        parser_owner=parser_owner,
                        source_call_line=call_line,
                        fixed_mode="",
                        relation="centralized_script_parser",
                        confidence="high",
                        applicable_chunks=tuple(chunks),
                        rejected_aliases=(),
                        source_evidence=(
                            (
                                f"{source_path}:{call_line}: centralized help "
                                f"owner "
                                f"{heredoc.owner}"
                            ),
                        ),
                    ),
                )

        return sorted(
            bindings,
            key=lambda item: (
                item.documentation_path,
                item.documented_owner,
                item.parser_path,
                item.parser_owner,
                item.source_call_line,
            ),
        )

    chunks_by_owner: dict[str, list[AliasChunk]] = {}

    for chunk in alias_chunks(text):
        chunks_by_owner.setdefault(chunk.owner, []).append(chunk)

    documented_owners = sorted(
        {
            heredoc.owner
            for heredoc in extract_help_heredocs(text)
            if heredoc.owner != "<file>"
        },
    )

    for documented_owner in documented_owners:
        for parser_owner, chunks in sorted(chunks_by_owner.items()):
            if parser_owner == documented_owner:
                continue

            policy = policy_for_binding(
                path,
                documented_owner,
                path,
                parser_owner,
            )
            if policy is None and not re.search(
                r"(?:^|_)(?:parse|parser)(?:_|$)|(?:^|_)check_args(?:_|$)",
                parser_owner,
            ):
                continue

            for call in parser_calls(text, documented_owner, parser_owner):
                bindings.append(
                    _binding_for_call(
                        path,
                        documented_owner,
                        path,
                        parser_owner,
                        chunks,
                        call,
                    ),
                )

    return sorted(
        bindings,
        key=lambda item: (
            item.documentation_path,
            item.documented_owner,
            item.parser_path,
            item.parser_owner,
            item.source_call_line,
        ),
    )


def parser_chunks_for_document(
    root: Path,
    path: str,
    text: str,
) -> tuple[list[AliasChunk], bool]:
    """
    Return parser chunks and whether they come from an external owner.
    """

    source = owning_source(root, path)

    if source is not None:
        owners = {"parse_args", "resolve_dir_scr", "check_args_light", "main"}
        return alias_chunks(source.read_text(encoding="utf-8"), owners), True

    if Path(path).name == "install_envs_entrypoint.sh":
        return alias_chunks(text, {"check_args_light"}), True

    chunks = alias_chunks(text)
    if chunks:
        return chunks, False

    return file_alias_chunks(text), True


def best_chunk(
    row: ParameterRow,
    chunks: list[AliasChunk],
    external: bool,
    delegated: tuple[AliasChunk, ...] = (),
) -> AliasChunk | None:
    """
    Bind one row to the parser group with the strongest exact overlap.
    """

    candidates = (
        chunks
        if external
        else [
            *[chunk for chunk in chunks if chunk.owner == row.owner],
            *delegated,
        ]
    )
    ranked = sorted(
        (
            (len(set(row.aliases) & set(chunk.aliases)), chunk)
            for chunk in candidates
        ),
        key=lambda item: item[0],
        reverse=True,
    )

    return ranked[0][1] if ranked and ranked[0][0] else None


def check_document(
    root: Path,
    path: str,
    text: str,
) -> tuple[list[Finding], dict[int, tuple[str, ...]]]:
    """
    Return mismatches and safe normalized row heads for one source.

    Parameters
    ----------
    root : Path
        Repository root used to resolve maintained paths.
    path : str
        Repository-relative path associated with the source.
    text : str
        Source text to inspect or normalize.

    Returns
    -------
    findings, replacements : tuple[
        list[Finding], dict[int, tuple[str, ...]]
    ]
        Alias-documentation findings for the selected owner.
    """

    chunks, external = parser_chunks_for_document(root, path, text)
    bindings = delegated_parser_bindings(root, path, text)
    delegated_by_owner: dict[str, tuple[AliasChunk, ...]] = {}

    for owner in {binding.documented_owner for binding in bindings}:
        delegated_by_owner[owner] = tuple(
            chunk
            for binding in bindings
            if binding.documented_owner == owner
            for chunk in binding.applicable_chunks
        )

    findings: list[Finding] = []
    replacements: dict[int, tuple[str, ...]] = {}
    usage_by_owner = usage_aliases_by_owner(text)
    observed_binding_keys = {
        (
            binding.documentation_path,
            binding.documented_owner,
            binding.parser_path,
            binding.parser_owner,
        )
        for binding in bindings
    }

    for policy in DELEGATED_PARSER_POLICIES:
        key = (
            policy.documentation_path,
            policy.documented_owner,
            policy.parser_path,
            policy.parser_owner,
        )

        if (
            policy.documentation_path == path
            and key not in observed_binding_keys
        ):
            findings.append(
                Finding(
                    path,
                    1,
                    policy.documented_owner,
                    (),
                    policy.applicable_logical_options,
                    (),
                    (
                        "source-reviewed delegated parser call is missing or "
                        "undiscoverable"
                    ),
                ),
            )

    for row in parameter_rows(text):
        if len(row.aliases) != len(set(row.aliases)):
            findings.append(
                Finding(
                    path,
                    row.line,
                    row.owner,
                    row.aliases,
                    (),
                    (),
                    "duplicate alias",
                ),
            )

            continue

        delegated = delegated_by_owner.get(row.owner, ())
        chunk = best_chunk(row, chunks, external, delegated)

        if chunk is None:
            rejected_bindings = [
                binding
                for binding in bindings
                if binding.documented_owner == row.owner
                and set(row.aliases) & set(binding.rejected_aliases)
            ]

            if rejected_bindings:
                findings.append(
                    Finding(
                        path,
                        row.line,
                        row.owner,
                        row.aliases,
                        (),
                        (),
                        (
                            "documented alias is rejected by delegated "
                            "fixed-mode contract"
                        ),
                    ),
                )

            continue

        expected, hidden = expected_aliases(row, chunk)
        documented = tuple(row.aliases)
        advertised_hidden = tuple(
            alias for alias in documented if alias in hidden
        )
        unaccepted = tuple(
            alias for alias in documented if alias not in chunk.aliases
        )

        if documented != expected or advertised_hidden or unaccepted:
            message = (
                "public alias set differs from parser and visibility evidence"
            )

            if chunk in delegated:
                message = (
                    "delegated public alias set differs from parser and "
                    "visibility evidence"
                )

            if advertised_hidden:
                message = (
                    "delegated hidden compatibility alias is documented"
                    if chunk in delegated
                    else "hidden compatibility alias is documented"
                )
            elif unaccepted:
                message = (
                    (
                        "documented alias is not accepted by its delegated "
                        "parser arm"
                    )
                    if chunk in delegated
                    else "documented alias is not accepted by its parser arm"
                )

            findings.append(
                Finding(
                    path,
                    row.line,
                    row.owner,
                    documented,
                    expected,
                    hidden,
                    message,
                ),
            )

        if (
            not advertised_hidden
            and not unaccepted
            and set(documented) <= set(expected)
        ):
            replacements[row.line] = expected

        usage_aliases = tuple(
            alias
            for alias in usage_by_owner.get(row.owner, ())
            if alias in chunk.aliases
        )
        invalid_usage = tuple(
            alias
            for alias in usage_aliases
            if (alias.startswith("-") and not alias.startswith("--"))
            or alias in hidden
        )
        public_long_usage = tuple(
            alias
            for alias in usage_aliases
            if alias.startswith("--") and alias in expected
        )

        if usage_aliases and (invalid_usage or len(public_long_usage) != 1):
            findings.append(
                Finding(
                    path,
                    row.line,
                    row.owner,
                    usage_aliases,
                    public_long_usage[:1],
                    hidden,
                    (
                        "Usage requires exactly one public canonical long "
                        "spelling and no short or hidden aliases"
                    ),
                ),
            )

    rows_by_owner: dict[str, list[ParameterRow]] = {}

    for row in parameter_rows(text):
        rows_by_owner.setdefault(row.owner, []).append(row)

    for binding in bindings:
        if binding.relation != "direct_delegated_parser_call":
            continue

        owner_rows = rows_by_owner.get(binding.documented_owner, [])

        if binding.confidence != "high":
            findings.append(
                Finding(
                    path,
                    binding.source_call_line,
                    binding.documented_owner,
                    (),
                    (),
                    (),
                    (
                        "delegated parser fixed mode differs from stable "
                        "source-reviewed contract"
                    ),
                ),
            )

            continue

        policy = policy_for_binding(
            binding.documentation_path,
            binding.documented_owner,
            binding.parser_path,
            binding.parser_owner,
        )
        applicable_options = tuple(
            logical_option(chunk) for chunk in binding.applicable_chunks
        )

        if (
            policy is not None
            and applicable_options != policy.applicable_logical_options
        ):
            findings.append(
                Finding(
                    path,
                    binding.source_call_line,
                    binding.documented_owner,
                    (),
                    policy.applicable_logical_options,
                    (),
                    (
                        "delegated parser policy names a logical option "
                        "absent from parser source"
                    ),
                ),
            )

        seen: set[tuple[str, ...]] = set()

        for chunk in binding.applicable_chunks:
            if chunk.aliases in seen:
                continue

            seen.add(chunk.aliases)
            expected, hidden = expected_aliases(
                ParameterRow(
                    binding.documented_owner,
                    binding.source_call_line,
                    (),
                    "  ",
                    "option",
                ),
                chunk,
            )
            if not expected or any(
                set(row.aliases) & set(chunk.aliases) for row in owner_rows
            ):
                continue

            findings.append(
                Finding(
                    path,
                    binding.source_call_line,
                    binding.documented_owner,
                    (),
                    expected,
                    hidden,
                    "delegated public alias row is missing from Parameters",
                ),
            )

    for heredoc in extract_help_heredocs(text):
        if heredoc.owner == "<file>":
            continue

        candidates = (
            chunks
            if external
            else [chunk for chunk in chunks if chunk.owner == heredoc.owner]
        )
        help_chunk = next(
            (chunk for chunk in candidates if "--help" in chunk.aliases),
            None,
        )
        if help_chunk is None:
            continue

        owner_rows = rows_by_owner.get(heredoc.owner, [])
        if any("--help" in row.aliases for row in owner_rows):
            continue

        expected, hidden = expected_aliases(
            ParameterRow(heredoc.owner, heredoc.start_line, (), "  ", "flag"),
            help_chunk,
        )
        findings.append(
            Finding(
                path,
                heredoc.start_line,
                heredoc.owner,
                (),
                expected,
                hidden,
                "accepted public help aliases are missing from Parameters",
            ),
        )

    return findings, replacements


def selected_documents(root: Path, paths: list[str] | None) -> list[str]:
    """
    Return public production shell documents with recognized help.
    """

    selected = paths or shell_paths(root)
    return [
        path
        for path in sorted(set(selected))
        if path.startswith(("bin/", "lib/bash/", "install/scripts/", "tests/"))
        and not path.startswith("artifacts/tests/")
        and (root / path).is_file()
        and extract_help_heredocs((root / path).read_text(encoding="utf-8"))
    ]


def wrapper_documents(paths: list[str]) -> list[str]:
    """
    Return only external help modules owned by execute/submit wrappers.
    """

    return [
        path
        for path in paths
        if re.fullmatch(
            r"lib/bash/help/help_(?:execute|submit)_[A-Za-z0-9_]+\.sh",
            path,
        )
    ]


def normalize_document(
    text: str,
    replacements: dict[int, tuple[str, ...]],
) -> str:
    """
    Apply exact alias-head replacements without changing descriptions.
    """

    lines = text.splitlines()

    for number, aliases in replacements.items():
        match = PARAMETER_ROW.match(lines[number - 1])
        if match is None:
            continue

        indentation = match.group("indent")
        alias_head = ", ".join(aliases)
        value_type = match.group("type")
        lines[number - 1] = f"{indentation}{alias_head} : {value_type}"

    return "\n".join(lines) + ("\n" if text.endswith("\n") else "")


def insert_missing_help_rows(
    root: Path,
    path: str,
    text: str,
) -> str:
    """
    Insert the one settled public help row into documented parser owners.
    """

    chunks, external = parser_chunks_for_document(root, path, text)
    documented = {
        row.owner for row in parameter_rows(text) if "--help" in row.aliases
    }
    insertions: list[tuple[int, tuple[str, ...]]] = []

    for heredoc in extract_help_heredocs(text):
        if heredoc.owner == "<file>":
            continue

        if heredoc.owner in documented:
            continue

        candidates = (
            chunks
            if external
            else [chunk for chunk in chunks if chunk.owner == heredoc.owner]
        )
        help_chunk = next(
            (chunk for chunk in candidates if "--help" in chunk.aliases),
            None,
        )
        if help_chunk is None:
            continue

        expected, _ = expected_aliases(
            ParameterRow(heredoc.owner, heredoc.start_line, (), "  ", "flag"),
            help_chunk,
        )
        body = list(heredoc.lines)

        for index, (_, line) in enumerate(body[:-1]):
            if line.strip() != "Parameters" or not UNDERLINE.fullmatch(
                body[index + 1][1].strip(),
            ):
                continue

            insertions.append((body[index + 1][0], expected))

            break

    lines = text.splitlines()

    for underline_line, aliases in sorted(insertions, reverse=True):
        insertion = [
            f"  {', '.join(aliases)} : flag",
            "    Display this help message and exit.",
        ]

        if underline_line >= len(lines) or lines[underline_line].strip():
            insertion.append("")

        lines[underline_line:underline_line] = insertion

    return "\n".join(lines) + ("\n" if text.endswith("\n") else "")


def delegated_parser_inventory(
    root: Path,
    paths: list[str] | None = None,
) -> list[dict[str, object]]:
    """
    Return deterministic source-derived cross-function parser mappings.

    Parameters
    ----------
    root : Path
        Repository root containing parser and help sources.
    paths : list[str] | None
        Optional bounded document paths.

    Returns
    -------
    records : list[dict[str, object]]
        Delegated parser bindings with source-derived alias evidence.
    """

    inventory: list[dict[str, object]] = []

    for path in selected_documents(root, paths):
        text = (root / path).read_text(encoding="utf-8")
        rows_by_owner: dict[str, list[ParameterRow]] = {}

        for row in parameter_rows(text):
            rows_by_owner.setdefault(row.owner, []).append(row)

        for binding in delegated_parser_bindings(root, path, text):
            accepted = tuple(
                dict.fromkeys(
                    alias
                    for chunk in binding.applicable_chunks
                    for alias in chunk.aliases
                ),
            )

            public: list[str] = []
            hidden: list[str] = []
            logical: list[str] = []

            for chunk in binding.applicable_chunks:
                expected, hidden_chunk = expected_aliases(
                    ParameterRow(
                        binding.documented_owner,
                        binding.source_call_line,
                        (),
                        "  ",
                        "option",
                    ),
                    chunk,
                )
                public.extend(expected)
                hidden.extend(hidden_chunk)
                logical.append(logical_option(chunk))

            public_tuple = tuple(dict.fromkeys(public))
            hidden_tuple = tuple(dict.fromkeys(hidden))
            related = set(accepted) | set(binding.rejected_aliases)

            documented_rows = rows_by_owner.get(binding.documented_owner, [])

            if binding.relation == "centralized_script_parser":
                documented_rows = [
                    row for rows in rows_by_owner.values() for row in rows
                ]

            documented = tuple(
                alias
                for row in documented_rows
                if set(row.aliases) & related
                for alias in row.aliases
            )
            missing = tuple(
                alias for alias in public_tuple if alias not in documented
            )
            overdocumented = tuple(
                alias for alias in documented if alias not in public_tuple
            )

            inventory.append(
                {
                    "documentation_path": binding.documentation_path,
                    "documented_owner": binding.documented_owner,
                    "parser_path": binding.parser_path,
                    "parser_owner": binding.parser_owner,
                    "source_call_location": (
                        f"{binding.parser_path}:{binding.source_call_line}"
                    ),
                    "fixed_mode": binding.fixed_mode,
                    "accepted_logical_options": list(dict.fromkeys(logical)),
                    "accepted_aliases": list(accepted),
                    "public_aliases": list(public_tuple),
                    "documented_aliases": list(documented),
                    "missing_public_aliases": list(missing),
                    "overdocumented_aliases": list(overdocumented),
                    "hidden_aliases": list(hidden_tuple),
                    "rejected_aliases": list(binding.rejected_aliases),
                    "relation": binding.relation,
                    "confidence": binding.confidence,
                    "source_evidence": list(binding.source_evidence),
                    "remediation_performed": bool(
                        policy_for_binding(
                            binding.documentation_path,
                            binding.documented_owner,
                            binding.parser_path,
                            binding.parser_owner,
                        )
                        and not missing
                        and not overdocumented,
                    ),
                },
            )

    if paths is None:
        canonical = list(inventory)

        for shim in sorted((root / "bin").glob("*.sh")):
            text = shim.read_text(encoding="utf-8")
            match = COMPATIBILITY_EXEC.search(text)
            if match is None:
                continue

            target_path = f"bin/{match.group('target')}"
            call_line = text[: match.start()].count("\n") + 1

            for row in canonical:
                if row["parser_path"] != target_path:
                    continue

                copied = dict(row)
                copied.update(
                    {
                        "delegator_path": str(shim.relative_to(root)),
                        "delegator_owner": shim.name,
                        "source_call_location": (
                            f"{shim.relative_to(root)}:{call_line}"
                        ),
                        "relation": "compatibility_script_delegator",
                        "source_evidence": [
                            f"{shim.relative_to(root)}:{call_line}: "
                            f"exec delegates '$@' to {match.group('target')}",
                        ],
                        "remediation_performed": False,
                    },
                )
                inventory.append(copied)

        forwarder = re.compile(
            r"^(?P<target>[A-Za-z_][A-Za-z0-9_]*)\s+\"\$@\"$",
        )

        for relative in shell_paths(root):
            source = root / relative
            if not source.is_file():
                continue

            text = source.read_text(encoding="utf-8")

            for owner, (start, _, body) in function_body_spans(text).items():
                match = forwarder.fullmatch(body.strip())
                if match is None:
                    continue

                target = match.group("target")

                for row in canonical:
                    if row["documented_owner"] != target:
                        continue

                    copied = dict(row)
                    copied.update(
                        {
                            "delegator_path": relative,
                            "delegator_owner": owner,
                            "source_call_location": f"{relative}:{start + 1}",
                            "relation": "compatibility_function_delegator",
                            "source_evidence": [
                                (
                                    f"{relative}:{start + 1}: {owner} "
                                    f"forwards '$@' to "
                                    f"{target}"
                                ),
                            ],
                            "remediation_performed": False,
                        },
                    )
                    inventory.append(copied)

    return sorted(
        inventory,
        key=lambda row: (
            str(row["documentation_path"]),
            str(row["documented_owner"]),
            str(row["parser_path"]),
            str(row["parser_owner"]),
        ),
    )


def scan_repository(
    root: Path,
    paths: list[str] | None = None,
) -> tuple[list[Finding], list[dict[str, object]]]:
    """
    Return findings plus a machine-readable row inventory.

    Parameters
    ----------
    root : Path
        Repository root containing the selected help and parser sources.
    paths : list[str] | None
        Optional bounded help-document paths.

    Returns
    -------
    findings, inventory : tuple[list[Finding], list[dict[str, object]]]
        Alias findings and complete parser/help ownership inventory rows.
    """

    findings: list[Finding] = []
    inventory: list[dict[str, object]] = []

    for path in selected_documents(root, paths):
        text = (root / path).read_text(encoding="utf-8")
        usage_by_owner = usage_aliases_by_owner(text)

        path_findings, _ = check_document(root, path, text)
        findings.extend(path_findings)
        chunks, external = parser_chunks_for_document(root, path, text)

        bindings = delegated_parser_bindings(root, path, text)
        delegated_by_owner: dict[str, tuple[AliasChunk, ...]] = {}

        for owner in {binding.documented_owner for binding in bindings}:
            delegated_by_owner[owner] = tuple(
                chunk
                for binding in bindings
                if binding.documented_owner == owner
                for chunk in binding.applicable_chunks
            )

        for row in parameter_rows(text):
            chunk = best_chunk(
                row,
                chunks,
                external,
                delegated_by_owner.get(row.owner, ()),
            )
            if chunk is None:
                continue

            expected, hidden = expected_aliases(row, chunk)
            logical_option = next(
                (
                    alias
                    for alias in reversed(expected)
                    if alias.startswith("--")
                ),
                expected[-1] if expected else "",
            )
            rejected: list[str] = []

            if (
                path == "lib/bash/help/help_execute_download_fastqs.sh"
                and logical_option == "--fil_in"
            ):
                rejected.extend(["-i", "--infile"])

            if (
                path == "lib/bash/help/help_execute_download_fastqs.sh"
                and logical_option == "--dir_eo"
            ):
                rejected.append("-eo")

            if (
                path == "lib/bash/core/check_args.sh"
                and logical_option == "--asgmt"
            ):
                rejected.append("--asmgt")

            if (
                path == "lib/bash/help/help_submit_compute_signal.sh"
                and logical_option == "--csv_usr_frg"
            ):
                rejected.extend(["-uf", "--usr_frg", "--usr-frg"])

            if (
                path == "lib/bash/help/help_compress_remove_files.sh"
                and logical_option == "--dry_run"
            ):
                rejected.extend(
                    [
                        "-ce",
                        "-cu",
                        "--chk_exc",
                        "--chk-exc",
                        "--chk_exu",
                        "--chk-exu",
                        "--dry",
                        "--dry-run",
                    ],
                )

            if path in {
                "install/scripts/install_envs_entrypoint.sh",
                "lib/bash/help/help_install_envs.sh",
            }:
                if logical_option == "--channels":
                    rejected.extend(
                        ["--channel", "--channel_list", "--channel-list"],
                    )
                elif logical_option == "--override_channels":
                    rejected.extend(
                        ["--override_channel", "--override-channel"],
                    )

            inventory.append(
                {
                    "path": path,
                    "owner": row.owner,
                    "line": row.line,
                    "logical_option": logical_option,
                    "usage_aliases": [
                        alias
                        for alias in usage_by_owner.get(row.owner, ())
                        if alias in chunk.aliases
                    ],
                    "documented_aliases": list(row.aliases),
                    "public_aliases": list(expected),
                    "restored_aliases": [
                        alias
                        for alias in expected
                        if alias.startswith("-") and not alias.startswith("--")
                    ],
                    "hidden_aliases": list(hidden),
                    "rejected_or_retired_aliases": sorted(set(rejected)),
                },
            )

    findings.extend(csv_short_alias_findings(inventory))

    return findings, inventory


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
        Parsed alias-audit and rewrite options.
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--fix", action="store_true")
    parser.add_argument("--inventory-json", action="store_true")
    parser.add_argument("--inventory-output", type=Path)
    parser.add_argument("--delegated-inventory-json", action="store_true")
    parser.add_argument("--delegated-inventory-output", type=Path)
    parser.add_argument("--wrapper-only", action="store_true")
    parser.add_argument("paths", nargs="*")

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """
    Check or normalize current-diff shell Parameter alias sets.

    Parameters
    ----------
    argv : list[str] | None
        Explicit arguments, or None to read the process arguments.

    Returns
    -------
    status : int
        Zero when alias documentation is valid and one otherwise.
    """

    args = parse_args(argv)
    root = args.root.resolve()
    selected = selected_documents(root, args.paths or None)

    if args.wrapper_only:
        selected = wrapper_documents(selected)

    if args.fix:
        for path in selected:
            target = root / path
            text = target.read_text(encoding="utf-8")
            _, replacements = check_document(root, path, text)
            normalized = normalize_document(text, replacements)
            normalized = insert_missing_help_rows(root, path, normalized)

            if normalized != text:
                target.write_text(normalized, encoding="utf-8")

    findings, inventory = scan_repository(root, selected)
    delegated_paths = selected if args.paths or args.wrapper_only else None
    delegated = delegated_parser_inventory(root, delegated_paths)

    if args.delegated_inventory_json:
        print(json.dumps(delegated, indent=2, sort_keys=True))
    elif args.inventory_json:
        print(json.dumps(inventory, indent=2, sort_keys=True))
    else:
        for finding in findings:
            print(finding.format())

    if args.inventory_output is not None:
        args.inventory_output.write_text(
            json.dumps(inventory, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    if args.delegated_inventory_output is not None:
        args.delegated_inventory_output.write_text(
            json.dumps(delegated, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    if findings:
        if not args.inventory_json and not args.delegated_inventory_json:
            print(f"HELP.PARAMETER.ALIAS_SET: {len(findings)} violation(s)")

        return 1

    if not args.inventory_json and not args.delegated_inventory_json:
        print("HELP.PARAMETER.ALIAS_SET: pass (zero unexplained findings)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
