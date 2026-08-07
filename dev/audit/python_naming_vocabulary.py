#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: python_naming_vocabulary.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


"""
Load the governed Python naming vocabulary.
"""

from __future__ import annotations

import dataclasses
import tomllib
from pathlib import Path

VOCABULARY_PATH = (
    Path(__file__).resolve().parents[1]
    / "config"
    / "python_naming_vocabulary.toml"
)
MATCH_KINDS = ("casefold_segment", "exact_identifier")
STATUSES = (
    "ordinary_short_word",
    "allowed_domain_term",
    "protected_external_or_interface_spelling",
    "prohibited_internal_shorthand",
    "review_candidate",
)
CONTEXTS = (
    "conventional_import_alias",
    "direct_adapter",
    "exact_identifier",
    "external_or_public_boundary",
    "identifier_segment",
)


@dataclasses.dataclass(frozen=True, order=True)
class VocabularyEntry:
    """
    Represent one governed Python vocabulary spelling.
    """

    spelling: str
    status: str
    match_kind: str
    contexts: tuple[str, ...]
    evidence_candidate: bool
    replacement: str | None


def load_vocabulary(
    path: Path = VOCABULARY_PATH,
) -> tuple[VocabularyEntry, ...]:
    """
    Load and validate the canonical naming vocabulary.
    """

    data = tomllib.loads(path.read_text(encoding="utf-8"))
    if data.get("schema_version") != 2:
        raise ValueError("Unsupported Python naming-vocabulary schema.")

    raw_entries = data.get("entry")
    if not isinstance(raw_entries, list):
        raise ValueError("Vocabulary entries must be a list.")

    entries: list[VocabularyEntry] = []
    seen: set[tuple[str, str]] = set()

    required = {
        "spelling",
        "status",
        "match_kind",
        "contexts",
        "evidence_candidate",
    }
    allowed = required | {"replacement"}
    expected_evidence_statuses = {
        "allowed_domain_term",
        "prohibited_internal_shorthand",
        "review_candidate",
    }

    for raw_entry in raw_entries:
        if not isinstance(raw_entry, dict):
            raise ValueError("Vocabulary entries must be tables.")
        if set(raw_entry) - allowed or not required <= set(raw_entry):
            raise ValueError(
                "Vocabulary entry fields are incomplete or unknown."
            )

        spelling = raw_entry["spelling"]
        status = raw_entry["status"]
        match_kind = raw_entry["match_kind"]
        contexts = raw_entry["contexts"]
        evidence_candidate = raw_entry["evidence_candidate"]
        replacement = raw_entry.get("replacement")

        if not isinstance(spelling, str) or not spelling:
            raise ValueError("Vocabulary spellings must be nonempty.")
        if status not in STATUSES:
            raise ValueError(f"Unknown vocabulary status: {status}.")
        if match_kind not in MATCH_KINDS:
            raise ValueError(f"Unknown vocabulary match kind: {match_kind}.")
        if (
            not isinstance(contexts, list)
            or not contexts
            or any(context not in CONTEXTS for context in contexts)
            or len(contexts) != len(set(contexts))
        ):
            raise ValueError(f"Invalid vocabulary contexts: {spelling}.")
        if not isinstance(evidence_candidate, bool):
            raise ValueError(
                f"Evidence-candidate membership must be Boolean: {spelling}.",
            )

        expected_evidence = (
            match_kind == "casefold_segment"
            and status in expected_evidence_statuses
        )
        if evidence_candidate != expected_evidence:
            raise ValueError(
                f"Invalid evidence-candidate membership: {spelling}.",
            )

        if status == "prohibited_internal_shorthand":
            if not isinstance(replacement, str) or not replacement:
                raise ValueError(
                    "Prohibited shorthand requires a replacement: "
                    f"{spelling}.",
                )
        elif replacement is not None:
            raise ValueError(
                f"Only prohibited shorthand records a replacement: {spelling}.",
            )

        identity = (
            match_kind,
            (
                spelling.casefold()
                if match_kind == "casefold_segment"
                else spelling
            ),
        )
        if identity in seen:
            raise ValueError(f"Duplicate vocabulary spelling: {identity}.")
        seen.add(identity)

        entries.append(
            VocabularyEntry(
                spelling=spelling,
                status=status,
                match_kind=match_kind,
                contexts=tuple(contexts),
                evidence_candidate=evidence_candidate,
                replacement=replacement,
            ),
        )

    return tuple(sorted(entries))


def spellings_for(
    *,
    match_kind: str,
    statuses: tuple[str, ...],
) -> frozenset[str]:
    """
    Return spellings for one matching mode and status selection.
    """

    unknown_statuses = set(statuses) - set(STATUSES)
    if unknown_statuses:
        raise ValueError(
            f"Unknown vocabulary statuses: {sorted(unknown_statuses)}.",
        )

    return frozenset(
        entry.spelling
        for entry in load_vocabulary()
        if entry.match_kind == match_kind and entry.status in statuses
    )


def ordinary_short_words() -> frozenset[str]:
    """
    Return ordinary short words excluded from abbreviation evidence.
    """

    return spellings_for(
        match_kind="casefold_segment",
        statuses=("ordinary_short_word",),
    )


def evidence_candidate_segments() -> frozenset[str]:
    """
    Return the historical abbreviation-evidence candidate projection.
    """

    return frozenset(
        entry.spelling
        for entry in load_vocabulary()
        if entry.match_kind == "casefold_segment" and entry.evidence_candidate
    )


def prohibited_internal_segments() -> frozenset[str]:
    """
    Return shorthand segments rejected in internal identifiers.
    """

    return spellings_for(
        match_kind="casefold_segment",
        statuses=("prohibited_internal_shorthand",),
    )
