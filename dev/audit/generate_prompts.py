#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: generate_prompts.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


"""
Build bounded, review-only Codex prompts from normalized findings.
"""

from __future__ import annotations

import difflib
import hashlib
import json
import re
import subprocess
from collections import Counter, defaultdict
from collections.abc import Callable
from itertools import pairwise
from pathlib import Path

REQUIRED_BEHAVIORAL_ANCHORS = {
    "bin/execute_download_fastqs.sh": {
        "dir_scr_forwarding": ("--dir_scr",),
        "argument_parsing": ("function parse_args()",),
        "environment_resolution": (
            "function setup_env()",
            'handle_env "${env_nam}"',
        ),
        "final_dispatch": (
            "setup_env               || return 1",
            "run_jobs                || return 1",
        ),
    },
    "bin/submit_download_fastqs.sh": {
        "help_source": (
            (
                '"${_dir_scr_hlp}/../lib/bash/help/help_submit_download_fastq'
                's.sh"'
            ),
        ),
        "early_help": (
            "-h|--hlp|--help)",
            "help_submit_download_fastqs >&2",
            "return 0",
        ),
        "bootstrap_dir_scr": (
            "-ds|--dir[_-]scr)",
            "(( i + 1 >= ${#args[@]} ))",
            '-z "${args[i + 1]:-}"',
            '"${args[i + 1]}" == -*',
            'printf "%s\\n" "${args[i + 1]}"',
            "required option '--dir_scr' was not supplied.",
            "return 1",
        ),
        "positional_validation": (
            "if [[ ${#args_pos[@]} -ne 8 ]]; then",
            'msg="but ${#args_pos[@]} were supplied."',
            'srr="${args_pos[0]}"',
            'url_1="${args_pos[1]}"',
            'url_2="${args_pos[2]}"',
            'dir_out="${args_pos[3]}"',
            'dir_sym="${args_pos[4]}"',
            'nam_cus="${args_pos[5]}"',
            'dir_eo="${args_pos[6]}"',
            'nam_job="${args_pos[7]}"',
        ),
        "worker_after_bootstrap": (
            'source_helpers_submit "${0##*/}" "${dir_scr}" || return 1',
            'parse_args "$@"',
            "validate_args",
            "run_downloads",
        ),
    },
    "lib/bash/help/help_execute_download_fastqs.sh": {
        "renamed_options_defaults": ("--env_nam", "default"),
        "environment_documentation": (
            "Runtime requirements:",
            "A compatible Conda environment providing the listed tools",
        ),
        "examples_context": ("Examples",),
    },
    "lib/bash/help/help_submit_download_fastqs.sh": {
        "complete_user_help": ("Usage", "Parameters", "Notes"),
    },
}

SOURCE_ONLY_ANCHORS = {
    ("bin/submit_download_fastqs.sh", "early_help"),
    ("bin/submit_download_fastqs.sh", "bootstrap_dir_scr"),
    ("bin/submit_download_fastqs.sh", "positional_validation"),
    ("bin/submit_download_fastqs.sh", "worker_after_bootstrap"),
    (
        "lib/bash/help/help_execute_download_fastqs.sh",
        "environment_documentation",
    ),
}


def _shell_source_units(path: str, source: str) -> list[dict[str, object]]:
    """
    Split source into complete top-level functions and intervening regions.

    Parameters
    ----------
    path : str
        Repository-relative path used in unit identities and diagnostics.
    source : str
        Complete Shell source to partition.

    Returns
    -------
    units : list[dict[str, object]]
        Ordered complete function and intervening top-level source units.

    Raises
    ------
    ValueError
        If a recognized function cannot be bounded as valid Shell source.
    """

    lines = source.splitlines()
    functions: list[tuple[str, int, int]] = []
    starts: list[tuple[str, int]] = []

    for number, line in enumerate(lines, 1):
        match = re.match(
            r"^function\s+([A-Za-z_][A-Za-z0-9_]*)\(\)\s*\{\s*$",
            line,
        )

        if match:
            starts.append((match.group(1), number))

    for index, (name, start) in enumerate(starts):
        next_start = (
            starts[index + 1][1] if index + 1 < len(starts) else len(lines) + 1
        )
        end = None

        for number in range(start + 1, next_start):
            if lines[number - 1].strip() != "}":
                continue

            fragment = "\n".join(lines[start - 1 : number]) + "\n"
            parsed = subprocess.run(
                ["bash", "-n"],
                input=fragment,
                text=True,
                capture_output=True,
                check=False,
            )

            if parsed.returncode == 0:
                end = number
                break

        if end is None:
            raise ValueError(
                (
                    f"cannot identify complete top-level shell function "
                    f"{name}: "
                    f"{path}:{start}"
                ),
            )

        functions.append((name, start, end))

    units: list[dict[str, object]] = []
    cursor = 1

    for name, start, end in functions:
        if cursor < start:
            units.append(
                {
                    "unit_kind": "top_level_region",
                    "unit_name": f"top_level_lines_{cursor}_{start - 1}",
                    "start_line": cursor,
                    "end_line": start - 1,
                },
            )

        units.append(
            {
                "unit_kind": "shell_function",
                "unit_name": name,
                "start_line": start,
                "end_line": end,
            },
        )
        cursor = end + 1

    if cursor <= len(lines):
        units.append(
            {
                "unit_kind": "top_level_region",
                "unit_name": f"top_level_lines_{cursor}_{len(lines)}",
                "start_line": cursor,
                "end_line": len(lines),
            },
        )

    if not units and lines:
        units.append(
            {
                "unit_kind": "top_level_region",
                "unit_name": f"top_level_lines_1_{len(lines)}",
                "start_line": 1,
                "end_line": len(lines),
            },
        )

    for unit in units:
        start, end = int(unit["start_line"]), int(unit["end_line"])
        unit.update(
            {"path": path, "source": "\n".join(lines[start - 1 : end]) + "\n"},
        )

    return units


def _relocation_origin(root: Path, path: str) -> str:
    """
    Return the tracked origin for a configured hard relocation.
    """

    config = root / "dev/config/path_relocations.json"
    if not config.is_file():
        return path

    value = json.loads(config.read_text(encoding="utf-8"))

    return str(value.get("exact", {}).get(path, path))


def _historical_source(root: Path, path: str) -> str:
    """
    Read the newest reachable revision that contains 'path'.
    """

    revisions = subprocess.run(
        ["git", "rev-list", "HEAD", "--", path],
        cwd=root,
        text=True,
        capture_output=True,
        check=False,
    )
    candidates = ["HEAD", *revisions.stdout.splitlines()]

    for revision in dict.fromkeys(candidates):
        source = subprocess.run(
            ["git", "show", f"{revision}:{path}"],
            cwd=root,
            text=True,
            capture_output=True,
            check=False,
        )
        if source.returncode == 0:
            return source.stdout

    raise ValueError(f"cannot read relocation origin: {path}")


def _changed_blocks(
    root: Path,
    path: str,
    line_count: int,
    baseline_path: str | None = None,
) -> list[dict[str, object]]:
    """
    Return exact hunks with independent before/current coordinates.
    """

    baseline_path = baseline_path or path

    if baseline_path != path:
        baseline = _historical_source(root, baseline_path)
        current = (root / path).read_text(encoding="utf-8")
        diff_text = "".join(
            difflib.unified_diff(
                baseline.splitlines(keepends=True),
                current.splitlines(keepends=True),
                fromfile=f"a/{baseline_path}",
                tofile=f"b/{path}",
                n=0,
            ),
        )
    else:
        tracked_result = subprocess.run(
            ["git", "ls-files", "--error-unmatch", "--", path],
            cwd=root,
            text=True,
            capture_output=True,
            check=False,
        )
        tracked = tracked_result.returncode == 0
        argv = ["git", "diff", "--no-ext-diff", "--unified=0", "--", path]
        normal_statuses = {0}

        if not tracked:
            argv = [
                "git",
                "diff",
                "--no-index",
                "--no-ext-diff",
                "--unified=0",
                "--",
                "/dev/null",
                str(root / path),
            ]
            normal_statuses.add(1)

        diff = subprocess.run(
            argv,
            cwd=root,
            text=True,
            capture_output=True,
            check=False,
        )
        if diff.returncode not in normal_statuses:
            raise ValueError(
                f"cannot read diff for semantic-unit selection: {path}",
            )

        diff_text = diff.stdout

    pattern = re.compile(
        r"^@@ -(\d+)(?:,(\d+))? \+(\d+)(?:,(\d+))? @@.*$",
        re.MULTILINE,
    )
    matches = list(pattern.finditer(diff_text))
    blocks: list[dict[str, object]] = []

    for index, match in enumerate(matches):
        old_count = 1 if match.group(2) is None else int(match.group(2))
        new_count = 1 if match.group(4) is None else int(match.group(4))
        end = (
            matches[index + 1].start()
            if index + 1 < len(matches)
            else len(diff_text)
        )
        blocks.append(
            {
                "old_start": int(match.group(1)),
                "old_count": old_count,
                "new_start": int(match.group(3)),
                "new_count": new_count,
                "diff_evidence": diff_text[match.start() : end].rstrip()
                + "\n",
            },
        )

    return blocks


def configured_semantic_units(
    root: Path,
    selection: dict[str, object],
) -> tuple[list[dict[str, object]], dict[str, object]]:
    """
    Render configured target changes and named context as complete shell units.

    Parameters
    ----------
    root : Path
        Repository root used to resolve maintained paths.
    selection : dict[str, object]
        Validated target-selection payload.

    Returns
    -------
    units, coverage : tuple[list[dict[str, object]], dict[str, object]]
        Selected semantic units and complete coverage metadata.

    Raises
    ------
    ValueError
        If a supplied value violates the validated contract.
    """

    evidence_selection = selection["evidence_selection"]
    reason = str(evidence_selection["target_selector"]["selection_reason"])

    rendered: list[dict[str, object]] = []
    uncovered: list[dict[str, object]] = []
    changed_blocks: list[dict[str, object]] = []
    overlapping: list[dict[str, object]] = []
    changed_top_level: list[str] = []
    segmentation_failures: list[dict[str, str]] = []

    changed_block_count = 0

    for target in selection["targets"]:
        path = str(target["path"])
        source = (root / path).read_text(encoding="utf-8")

        try:
            units = _shell_source_units(path, source)
        except ValueError as exc:
            segmentation_failures.append({"path": path, "error": str(exc)})

            continue

        for left, right in pairwise(units):
            if int(left["end_line"]) >= int(right["start_line"]):
                overlapping.append(
                    {
                        "path": path,
                        "left": left["unit_name"],
                        "right": right["unit_name"],
                        "source_state": "current",
                    },
                )

        baseline_path = _relocation_origin(root, path)

        if baseline_path != path:
            baseline_source = _historical_source(root, baseline_path)
        else:
            baseline = subprocess.run(
                ["git", "show", f"HEAD:{baseline_path}"],
                cwd=root,
                text=True,
                capture_output=True,
                check=False,
            )
            baseline_source = (
                baseline.stdout if baseline.returncode == 0 else ""
            )

        before_units: list[dict[str, object]] = []

        if baseline_source:
            try:
                before_units = _shell_source_units(path, baseline_source)
            except ValueError as exc:
                segmentation_failures.append(
                    {"path": path, "error": f"HEAD: {exc}"},
                )

                continue

        blocks = _changed_blocks(
            root,
            path,
            len(source.splitlines()),
            baseline_path,
        )
        changed_block_count += len(blocks)
        current_records: dict[int, dict[str, object]] = {}
        before_records: dict[int, dict[str, object]] = {}

        def intersects(
            unit: dict[str, object],
            start: int,
            count: int,
        ) -> bool:
            return (
                count > 0
                and int(unit["start_line"]) <= start + count - 1
                and int(unit["end_line"]) >= start
            )

        def current_record(
            index: int,
            units: list[dict[str, object]] = units,
            records: dict[int, dict[str, object]] = current_records,
        ) -> dict[str, object]:
            unit = units[index]

            if index not in records:
                records[index] = {
                    **unit,
                    "source_state": "current",
                    "evidence_role": "target",
                    "selection_reason": reason,
                    "source_boundary": {
                        "start_line": unit["start_line"],
                        "end_line": unit["end_line"],
                    },
                    "diff_evidence": [],
                }

            return records[index]

        def before_record(
            index: int,
            units: list[dict[str, object]] = before_units,
            records: dict[int, dict[str, object]] = before_records,
        ) -> dict[str, object]:
            unit = units[index]

            if index not in records:
                records[index] = {
                    **unit,
                    "source_state": "before",
                    "evidence_role": "target",
                    "selection_reason": reason,
                    "source_boundary": {
                        "start_line": unit["start_line"],
                        "end_line": unit["end_line"],
                    },
                    "diff_evidence": [],
                }

            return records[index]

        for block_index, block in enumerate(blocks, 1):
            current_indexes = [
                index
                for index, unit in enumerate(units)
                if intersects(
                    unit,
                    int(block["new_start"]),
                    int(block["new_count"]),
                )
            ]
            before_indexes = [
                index
                for index, unit in enumerate(before_units)
                if intersects(
                    unit,
                    int(block["old_start"]),
                    int(block["old_count"]),
                )
            ]

            if int(block["new_count"]) == 0:
                for before_index in before_indexes:
                    before_unit = before_units[before_index]

                    if before_unit["unit_kind"] == "shell_function":
                        same_function = [
                            index
                            for index, unit in enumerate(units)
                            if unit["unit_kind"] == "shell_function"
                            and unit["unit_name"] == before_unit["unit_name"]
                        ]
                        current_indexes.extend(same_function)

                        continue

                    old_end = (
                        int(block["old_start"]) + int(block["old_count"]) - 1
                    )
                    deletes_complete_region = int(block["old_start"]) <= int(
                        before_unit["start_line"],
                    ) and old_end >= int(before_unit["end_line"])

                    if not deletes_complete_region:
                        new_start = int(block["new_start"])

                        for index, unit in enumerate(units):
                            if unit["unit_kind"] != "top_level_region":
                                continue

                            unit_start = int(unit["start_line"]) - 1
                            unit_end = int(unit["end_line"])

                            if unit_start <= new_start <= unit_end:
                                current_indexes.append(index)

            current_indexes = sorted(set(current_indexes))
            before_indexes = sorted(set(before_indexes))
            diff_text = str(block["diff_evidence"])

            for index in current_indexes:
                record = current_record(index)

                if diff_text not in record["diff_evidence"]:
                    record["diff_evidence"].append(diff_text)

            for index in before_indexes:
                before_unit = before_units[index]
                current_identities = {
                    (
                        units[current_index]["unit_kind"],
                        units[current_index]["unit_name"],
                    )
                    for current_index in current_indexes
                }
                before_identity = (
                    before_unit["unit_kind"],
                    before_unit["unit_name"],
                )
                has_current_region = before_unit[
                    "unit_kind"
                ] == "top_level_region" and any(
                    kind == "top_level_region"
                    for kind, _ in current_identities
                )
                has_current_identity = (
                    before_identity in current_identities or has_current_region
                )

                if int(block["new_count"]) == 0 and not has_current_identity:
                    record = before_record(index)

                    if diff_text not in record["diff_evidence"]:
                        record["diff_evidence"].append(diff_text)

            block_record = {
                "block_id": f"{path}:hunk-{block_index}",
                **{
                    key: value
                    for key, value in block.items()
                    if key != "diff_evidence"
                },
                "diff_evidence_fingerprint": _sha256_text(diff_text),
                "exact_diff_retained": bool(diff_text),
                "complete_current_units": [
                    f"current:{path}:{units[index]['unit_name']}"
                    for index in current_indexes
                ],
                "complete_before_units": [
                    f"before:{path}:{before_units[index]['unit_name']}"
                    for index in before_indexes
                ],
            }

            block_record["represented"] = (
                bool(diff_text)
                and (int(block["new_count"]) == 0 or bool(current_indexes))
                and (int(block["old_count"]) == 0 or bool(before_indexes))
            )
            changed_blocks.append(block_record)

            if not block_record["represented"]:
                uncovered.append(block_record)

        target_records = [
            current_records[index] for index in sorted(current_records)
        ] + [before_records[index] for index in sorted(before_records)]

        for record in target_records:
            if record["unit_kind"] == "top_level_region":
                changed_top_level.append(
                    f"{record['source_state']}:{path}:{record['unit_name']}",
                )

        rendered.extend(target_records)

    for selector in evidence_selection["context_selectors"]:
        path = str(selector["path"])
        source = (root / path).read_text(encoding="utf-8")
        matches = [
            unit
            for unit in _shell_source_units(path, source)
            if unit["unit_kind"] == "shell_function"
            and unit["unit_name"] == selector["name"]
        ]
        if len(matches) != 1:
            raise ValueError(
                (
                    f"configured context function must resolve exactly once: "
                    f"{path}:{selector['name']}"
                ),
            )

        unit = matches[0]
        rendered.append(
            {
                **unit,
                "source_state": "current",
                "evidence_role": "context",
                "selection_reason": selector["selection_reason"],
                "source_boundary": {
                    "start_line": unit["start_line"],
                    "end_line": unit["end_line"],
                },
                "diff_evidence": [],
            },
        )

    return (
        rendered,
        {
            "configured_semantic_unit_count": len(rendered),
            "function_unit_count": sum(
                unit["unit_kind"] == "shell_function" for unit in rendered
            ),
            "top_level_region_count": sum(
                unit["unit_kind"] == "top_level_region" for unit in rendered
            ),
            "changed_block_count": changed_block_count,
            "changed_blocks": changed_blocks,
            "uncovered_changed_blocks": uncovered,
            "overlapping_units": overlapping,
            "changed_top_level_regions": changed_top_level,
            "segmentation_failures": segmentation_failures,
            "all_changed_blocks_covered": not uncovered
            and not overlapping
            and not segmentation_failures,
            "evidence_truncated": False,
        },
    )


def _sha256_text(text: str) -> str:
    """
    Return a prefixed SHA-256 fingerprint for text.
    """

    return "sha256:" + hashlib.sha256(text.encode("utf-8")).hexdigest()


def _current_evidence_text(evidence: str, evidence_kind: str) -> str:
    """
    Return evidence text eligible to establish current behavior.
    """

    if evidence_kind != "complete_line_diff_window":
        return evidence

    return "\n".join(
        line
        for line in evidence.splitlines()
        if not line.startswith(("--- ", "-"))
    )


def _table_cells(line: str) -> list[str]:
    """
    Parse one pipe table row while preserving escaped literal pipes.
    """

    if not line.startswith("|") or not line.endswith("|"):
        raise ValueError("Markdown table rows must begin and end with a pipe")

    cells: list[str] = []
    current: list[str] = []
    escaped = False

    for character in line[1:-1]:
        if escaped:
            current.extend(("\\", character))
            escaped = False
        elif character == "\\":
            escaped = True
        elif character == "|":
            cells.append("".join(current).strip())
            current = []
        else:
            current.append(character)

    if escaped:
        current.append("\\")

    cells.append("".join(current).strip())

    return cells


def _table_cell(value: object) -> str:
    """
    Escape one human-facing Markdown table cell.
    """

    return (
        str(value)
        .replace("\\", "\\\\")
        .replace("|", "\\|")
        .replace("\n", "<br>")
    )


def validate_markdown_table(lines: list[str]) -> None:
    """
    Reject malformed pipe tables before a review bundle is written.

    Parameters
    ----------
    lines : list[str]
        Complete table rows beginning with a header and separator.

    Raises
    ------
    ValueError
        If the table is incomplete, ragged, or uses invalid alignment cells.
    """

    if len(lines) < 2:
        raise ValueError("Markdown table requires header and separator rows")

    headers = _table_cells(lines[0])
    separators = _table_cells(lines[1])
    if not headers or len(headers) != len(separators):
        raise ValueError(
            "Markdown table header and separator column counts differ",
        )

    if any(not re.fullmatch(r":?-{3,}:?", cell) for cell in separators):
        raise ValueError("Markdown table has invalid alignment marker")

    for row in lines[2:]:
        if len(_table_cells(row)) != len(headers):
            raise ValueError("Markdown table has a ragged body row")


def markdown_table(
    headers: list[str],
    rows: list[list[object]],
    alignments: list[str] | None = None,
) -> str:
    """
    Render and validate one structurally sound Markdown pipe table.
    """

    if not headers:
        raise ValueError("Markdown table requires at least one header")

    markers = alignments or [":---"] * len(headers)
    lines = [
        "| " + " | ".join(_table_cell(value) for value in headers) + " |",
        "| " + " | ".join(markers) + " |",
    ]
    lines.extend(
        "| " + " | ".join(_table_cell(value) for value in row) + " |"
        for row in rows
    )

    validate_markdown_table(lines)

    return "\n".join(lines)


def validate_bundle_tables(markdown: str) -> int:
    """
    Validate every generated non-fenced pipe table and return its count.
    """

    lines = markdown.splitlines()
    fence_marker: str | None = None
    count = 0
    index = 0

    while index < len(lines):
        marker = re.match(r"^(`{3,})", lines[index])

        if fence_marker is None and marker:
            fence_marker = marker.group(1)
            index += 1

            continue

        if fence_marker is not None and lines[index] == fence_marker:
            fence_marker = None
            index += 1

            continue

        if (
            fence_marker is None
            and lines[index].startswith("|")
            and index + 1 < len(lines)
            and lines[index + 1].startswith("|")
        ):
            table = [lines[index], lines[index + 1]]
            index += 2

            while index < len(lines) and lines[index].startswith("|"):
                table.append(lines[index])
                index += 1

            validate_markdown_table(table)

            count += 1

            continue

        index += 1

    return count


def write_configured_bundle(
    report_dir: Path,
    selection: dict[str, object],
    targets: list[dict[str, object]],
    facts: list[dict[str, object]],
    findings: list[dict[str, object]],
    limitations: list[dict[str, object]],
    semantic_units: list[dict[str, object]],
    unit_coverage: dict[str, object],
    run_id: str,
) -> list[dict[str, object]]:
    """
    Write one declaratively named five-artifact semantic-review Markdown view.

    Parameters
    ----------
    report_dir : Path
        Destination report directory.
    selection : dict[str, object]
        Validated package, linkage, and size-limit configuration.
    targets : list[dict[str, object]]
        Selected source records and their current fingerprints.
    facts : list[dict[str, object]]
        Deterministic facts rendered into the review document.
    findings : list[dict[str, object]]
        Deterministic findings rendered by topic.
    limitations : list[dict[str, object]]
        Adapter limitations retained for reviewer interpretation.
    semantic_units : list[dict[str, object]]
        Complete source and change units supplied for semantic review.
    unit_coverage : dict[str, object]
        Fail-closed coverage proof for the supplied semantic units.
    run_id : str
        Stable run identity attached to the package record.

    Returns
    -------
    records : list[dict[str, object]]
        One semantic-review artifact record for the written Markdown.

    Raises
    ------
    ValueError
        If coverage is incomplete or a linked package exceeds its size limits.
    """

    invalid_coverage = (
        unit_coverage.get("all_changed_blocks_covered") is not True
        or bool(unit_coverage.get("uncovered_changed_blocks"))
        or bool(unit_coverage.get("overlapping_units"))
        or bool(unit_coverage.get("segmentation_failures"))
        or unit_coverage.get("evidence_truncated") is not False
        or any(
            block.get("represented") is not True
            or block.get("exact_diff_retained") is not True
            for block in unit_coverage.get("changed_blocks", [])
        )
    )

    if invalid_coverage:
        raise ValueError(
            "incomplete semantic-unit coverage; refusing bundle generation",
        )

    package = selection["package"]

    def fence(language: str, value: str) -> str:
        longest = max(
            (len(item) for item in re.findall(r"`+", value)),
            default=0,
        )
        marker = "`" * max(3, longest + 1)

        return f"{marker}{language}\n{value.rstrip()}\n{marker}"

    target_table = markdown_table(
        [
            "Path",
            "Target role",
            "Internal subcohort",
            "Git state",
            "Checks",
            "Findings",
        ],
        [
            [
                f"`{target['path']}`",
                target.get("target_role", target["role"]),
                target.get("subcohort", "unclassified"),
                ", ".join(target.get("git_state_labels", [])),
                ", ".join(target.get("checks_run", [])) or "none",
                target.get("findings_count", 0),
            ]
            for target in targets
        ],
    )

    rule_table = markdown_table(
        [
            "Rule",
            "Configured path scope",
            "Enforcement",
            "Exact matched targets",
            "Classification",
        ],
        [
            [
                f"`{scope['rule_id']}`",
                ", ".join(f"`{path}`" for path in scope["paths"]),
                scope["enforcement"],
                ", ".join(
                    f"`{path}`"
                    for path in selection.get("rule_scope_report", {})
                    .get(scope["rule_id"], {})
                    .get("matched_targets", [])
                ),
                scope["classification"],
            ]
            for scope in selection["rule_scopes"]
        ],
    )

    edge_table = markdown_table(
        [
            "From",
            "To",
            "Edge kind",
            "Reachability",
            "Runtime status",
            "Evidence reference",
        ],
        [
            [
                f"`{edge['from']}`",
                f"`{edge['to']}`",
                edge["kind"],
                edge["reachability"],
                edge["runtime_status"],
                edge.get(
                    "evidence_reference",
                    "not applicable or not supplied",
                ),
            ]
            for edge in selection["dependency_edges"]
        ],
    )

    finding_table = markdown_table(
        ["Rule", "Path", "Finding"],
        [
            [f"`{row['rule_id']}`", f"`{row['path']}`", row["message"]]
            for row in findings
        ],
    )

    limitation_table = markdown_table(
        ["Rule", "Adapter", "Limitation"],
        [
            [f"`{row['rule_id']}`", row["adapter"], row["limitation"]]
            for row in limitations
        ],
    )

    questions = "\n".join(
        f"{index}. [{question['topic']}] {question['question']} Paths: "
        + ", ".join(f"`{path}`" for path in question["paths"])
        for index, question in enumerate(selection["semantic_questions"], 1)
    )

    evidence_blocks: list[tuple[str, str]] = []
    target_subcohorts = {
        str(target["path"]): str(target.get("subcohort", "package_overhead"))
        for target in targets
    }

    for unit in semantic_units:
        block = (
            f"### `{unit['path']}` — `{unit['unit_name']}`\n\n"
            f"Evidence role: `{unit['evidence_role']}`  \n"
            f"Source state: `{unit.get('source_state', 'current')}`  \n"
            f"Unit kind: `{unit['unit_kind']}`  \n"
            f"Source boundary: lines {unit['start_line']}–{unit['end_line']}  "
            f"\n"
            f"Selection reason: {unit['selection_reason']}\n\n"
            f"Rendered complete unit ({unit.get('source_state', 'current')} "
            f"state):\n\n" + fence("bash", str(unit["source"]))
        )

        if "before_source" in unit:
            boundary = unit.get("before_source_boundary", {})
            block += (
                "\n\nBefore-state complete unit:"
                f" lines "
                f"{boundary.get('start_line', 'unknown')}–"
                f"{boundary.get('end_line', 'unknown')}\n\n"
                + fence("bash", str(unit["before_source"]))
            )

        if unit.get("diff_evidence"):
            block += "\n\nExact change evidence:\n\n" + "\n\n".join(
                fence("diff", str(diff)) for diff in unit["diff_evidence"]
            )

        category = (
            "context_only_callers"
            if unit["evidence_role"] == "context"
            else target_subcohorts[str(unit["path"])]
        )
        evidence_blocks.append((category, block))

    prefix = f"# {package['title']}\n\nRun ID: `{run_id}`\n\n"
    prefix += "## Package identity and review boundary\n\n"

    linkage = selection.get("linked_package")

    if linkage:
        sibling_edges = linkage["sibling_cross_package_edges"]
        sibling_relationships = (
            "; ".join(
                f"`{edge['from']}` → `{edge['to']}`" for edge in sibling_edges
            )
            if sibling_edges
            else "none"
        )
        prefix += (
            f"This is linked child `{linkage['child_package_id']}`, not the "
            f"entire closure."
            f" Review it together with sibling"
            f" "
            f"`{linkage['sibling_package_id']}`. "
            "The complete graph is the disjoint union of both edge partitions."
            " Shared graph fingerprint:"
            f" "
            f"`{linkage['graph_ownership_fingerprint']}`; "
            f"complete edge count: {linkage['umbrella_edge_count']}.\n\n"
            "Answer the shared umbrella dependency-closure question once for "
            "the linked pair.\n\n"
            f"Cross-package edges rendered by the sibling: "
            f"{sibling_relationships}.\n\n"
        )

    prefix += (
        "\n\n".join(str(statement) for statement in package["statements"])
        + "\n\n"
    )
    prefix += (
        "Semantic-review-only topics: "
        + "; ".join(str(topic) for topic in package["semantic_only_topics"])
        + ".\n\n"
    )

    prefix += "## Target status\n\n" + target_table + "\n\n"
    prefix += "## Scoped rule applicability\n\n" + rule_table + "\n\n"

    prefix += (
        (
            "## Dependency edges\n\nEdge kind, reachability, and runtime "
            "status are independent classifications and must not be collapsed "
            "into a generic dependency-proven label.\n\n"
        )
        + edge_table
        + "\n\n"
    )

    prefix += "## Declarative semantic questions\n\n" + questions + "\n\n"
    prefix += "## Deterministic findings\n\n" + finding_table + "\n\n"
    prefix += "## Adapter limitations\n\n" + limitation_table + "\n\n"

    prefix += "## Complete semantic-unit evidence\n\n"
    prefix += (
        "All configured evidence is rendered without truncation. Omission "
        "markers may appear only outside complete units; this package emits "
        "none.\n\n"
    )
    suffix = (
        "## Supplied review package\n\n"
        + "\n".join(f"- `{path}`" for path in package["supplied_artifacts"])
        + "\n"
    )

    parts: list[tuple[str, str]] = [("package_overhead", prefix)]
    parts.extend(
        (category, block + "\n\n") for category, block in evidence_blocks
    )
    parts.append(("package_overhead", suffix))
    body = "".join(text for _, text in parts)

    contribution_names = (
        "shared_runtime",
        "download_production_help",
        "shared_test_infrastructure",
        "download_tests_registration",
        "context_only_callers",
        "package_overhead",
    )
    contributions = {
        name: {"lines": 0, "bytes": 0} for name in contribution_names
    }

    for category, text in parts:
        contributions[category]["lines"] += text.count("\n")
        contributions[category]["bytes"] += len(text.encode("utf-8"))

    semantic_path = report_dir / str(package["semantic_review_path"])
    semantic_path.parent.mkdir(parents=True, exist_ok=True)
    semantic_path.write_text(body, encoding="utf-8")

    line_count = len(body.splitlines())
    byte_count = len(body.encode("utf-8"))
    limits = package["size_limits"]

    if linkage and (
        line_count > int(limits["semantic_markdown_max_lines"])
        or byte_count > int(limits["semantic_markdown_max_bytes"])
    ):
        raise ValueError(
            "linked child semantic Markdown exceeds configured size limits",
        )

    exceeds_byte_limit = byte_count > int(
        limits["semantic_markdown_max_bytes"],
    )
    exceeds_line_limit = line_count > int(
        limits["semantic_markdown_max_lines"],
    )

    record = {
        "schema_version": 1,
        "run_id": run_id,
        "bundle_id": package["bundle_id"],
        "path": package["semantic_review_path"],
        "primary_paths": [
            target["path"] for target in targets if target["role"] == "primary"
        ],
        "target_fingerprints": {
            target["path"]: target["content_fingerprint"] for target in targets
        },
        "prompt_fingerprint": _sha256_text(body),
        "markdown_lines": line_count,
        "markdown_bytes": byte_count,
        "markdown_tables_validated": validate_bundle_tables(body),
        "markdown_contributions": contributions,
        "semantic_unit_coverage": unit_coverage,
        "size_limits": {
            **limits,
            "semantic_markdown_exceeds_byte_limit": exceeds_byte_limit,
            "semantic_markdown_exceeds_line_limit": exceeds_line_limit,
        },
        "fact_count": len(facts),
    }

    if linkage:
        record.update(
            {
                "child_package_id": linkage["child_package_id"],
                "sibling_package_id": linkage["sibling_package_id"],
                "umbrella_run_id": linkage["umbrella_run_id"],
                "config_fingerprint": linkage["config_fingerprint"],
                "target_ownership_fingerprint": linkage[
                    "target_ownership_fingerprint"
                ],
                "graph_ownership_fingerprint": linkage[
                    "graph_ownership_fingerprint"
                ],
            },
        )

    return [record]


def build_batches(
    findings: list[dict[str, object]],
    ledger: dict[str, dict[str, object]],
    max_paths: int,
    max_findings: int,
) -> list[list[dict[str, object]]]:
    """
    Group findings by subsystem/risk while enforcing hard batch bounds.
    """

    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)

    for finding in findings:
        row = ledger.get(str(finding["path"]), {})
        grouped[
            (
                str(row.get("subsystem", "repository")),
                str(row.get("risk", "medium")),
            )
        ].append(finding)

    batches: list[list[dict[str, object]]] = []

    for group in grouped.values():
        batch: list[dict[str, object]] = []
        paths: set[str] = set()

        for finding in group:
            path = str(finding["path"])

            if batch and (
                len(batch) >= max_findings
                or (path not in paths and len(paths) >= max_paths)
            ):
                batches.append(batch)
                batch, paths = [], set()

            batch.append(finding)
            paths.add(path)

        if batch:
            batches.append(batch)

    return batches


def write_prompt_bundle(
    report_dir: Path,
    batches: list[list[dict[str, object]]],
    ledger: dict[str, dict[str, object]],
    rules: dict[str, dict[str, object]],
    rule_excerpt: Callable[[dict[str, object]], str],
    source_excerpt: Callable[[str, int | None], str],
    max_excerpt_chars: int,
) -> list[dict[str, object]]:
    """
    Write bounded prompt files and return their NDJSON metadata records.

    Parameters
    ----------
    report_dir : Path
        Audit report directory that owns the prompt subdirectory.
    batches : list[list[dict[str, object]]]
        Ordered bounded finding batches.
    ledger : dict[str, dict[str, object]]
        Source-ledger records keyed by maintained path.
    rules : dict[str, dict[str, object]]
        Rule metadata keyed by rule ID.
    rule_excerpt : Callable[[dict[str, object]], str]
        Renderer for the applicable normative provision.
    source_excerpt : Callable[[str, int | None], str]
        Renderer for a bounded source window around one finding.
    max_excerpt_chars : int
        Hard character bound applied to rendered excerpts.

    Returns
    -------
    records : list[dict[str, object]]
        Prompt paths, batch identities, fingerprints, and target metadata.
    """

    prompt_dir = report_dir / "prompts"
    prompt_dir.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, object]] = []

    for index, batch in enumerate(batches, 1):
        batch_id = f"batch-{index:03d}"

        paths = sorted(
            {str(item["path"]) for item in batch if item["path"] != "."},
        )
        rule_ids = sorted({str(item["rule_id"]) for item in batch})

        rule_text = "\n\n".join(
            f"[{rule_id}]\n{rule_excerpt(rules[rule_id])}"
            for rule_id in rule_ids
        )

        details: list[str] = []

        for finding in batch:
            path = str(finding["path"])
            line = finding["line"]
            details.append(
                "- "
                + json.dumps(
                    {
                        "rule_id": finding["rule_id"],
                        "path": path,
                        "line": line,
                        "message": finding["message"],
                        "evidence": finding["evidence"],
                        "content_fingerprint": finding["content_fingerprint"],
                    },
                    sort_keys=True,
                ),
            )

        excerpt_parts = []

        for path in paths:
            first_line = next(
                (item["line"] for item in batch if item["path"] == path),
                None,
            )
            excerpt = source_excerpt(path, first_line)[:max_excerpt_chars]
            fingerprint = ledger[path]["content_fingerprint"]
            header = f"--- {path} ({fingerprint}) ---\n"
            excerpt_parts.append(f"{header}{excerpt}")

        excerpts = "\n\n".join(excerpt_parts)

        verification_lines = []

        for rule_id in rule_ids:
            command = rules[rule_id].get(
                "command",
                rules[rule_id].get("commands"),
            )
            verification_lines.append(f"- {command}")

        verification = "\n".join(verification_lines)

        body = (
            f"# Cleanup audit review: {batch_id}\n\n"
            "Task: review only the allowed target paths below. Do not edit "
            "files or produce a patch.\n\n"
            "Allowed target paths:\n"
            f"{chr(10).join(f'- {path}' for path in paths)}\n\n"
            "Applicable rules:\n"
            f"{', '.join(rule_ids)}\n\n"
            "Relevant rule excerpts:\n"
            f"{rule_text}\n\n"
            "Normalized findings:\n"
            f"{chr(10).join(details)}\n\n"
            "Relevant source or diff excerpts:\n"
            f"{excerpts}\n\n"
            "Required verification commands:\n"
            f"{verification}\n\n"
            "Prohibited actions: source edits, staging, commits, resets, "
            "cleanup, autofixes, or model-generated patches.\n\n"
            "Expected response format: one JSON object per target path with "
            "`path`, `disposition`, `rationale`, `verification`, and "
            "`stale_fingerprint`.\n"
        )

        fingerprint = _sha256_text(body)

        (prompt_dir / f"{batch_id}.md").write_text(body, encoding="utf-8")
        records.append(
            {
                "batch_id": batch_id,
                "paths": paths,
                "target_fingerprints": {
                    path: ledger[path]["content_fingerprint"] for path in paths
                },
                "finding_ids": [item["finding_id"] for item in batch],
                "prompt_fingerprint": fingerprint,
                "path": f"prompts/{batch_id}.md",
            },
        )

    return records


def write_pilot_bundle(
    report_dir: Path,
    targets: list[dict[str, object]],
    context: list[dict[str, object]],
    facts: list[dict[str, object]],
    findings: list[dict[str, object]],
    policy_questions: list[dict[str, object]],
    coverage: dict[str, object],
    source_excerpt: Callable[[str], str],
    max_excerpt_chars: int,
    rule_excerpts: dict[str, str] | None = None,
    diff_excerpts: dict[str, str] | None = None,
    run_id: str | None = None,
    anchor_evidence: dict[str, dict[str, dict[str, str]]] | None = None,
) -> list[dict[str, object]]:
    """
    Write a concise, fenced human review view over machine artifacts.

    Parameters
    ----------
    report_dir : Path
        Directory in which review artifacts are written.
    targets : list[dict[str, object]]
        Selected review targets and their stable metadata.
    context : list[dict[str, object]]
        Supporting review context records included in the bundle.
    facts : list[dict[str, object]]
        Normalized evidence facts rendered in the review bundle.
    findings : list[dict[str, object]]
        Mutable finding collection populated by the check.
    policy_questions : list[dict[str, object]]
        Unresolved policy questions rendered for review.
    coverage : dict[str, object]
        Coverage reconciliation rendered for human review.
    source_excerpt : Callable[[str], str]
        Callable returning bounded source text for one path.
    max_excerpt_chars : int
        Maximum characters retained in each source excerpt.
    rule_excerpts : dict[str, str] | None
        Optional normative excerpts keyed by rule identifier.
    diff_excerpts : dict[str, str] | None
        Optional bounded diff excerpts keyed by source path.
    run_id : str | None
        Optional stable identifier for the evidence run.
    anchor_evidence : dict[str, dict[str, dict[str, str]]] | None
        Optional evidence indexed by review anchor.

    Returns
    -------
    records : list[dict[str, object]]
        Metadata records for the written human-review bundle.

    Raises
    ------
    ValueError
        If a supplied value violates the validated contract.
    """

    def fence(language: str, text: str) -> str:
        longest = max(
            (len(value) for value in __import__("re").findall(r"`+", text)),
            default=0,
        )
        marker = "`" * max(3, longest + 1)

        return f"{marker}{language}\n{text.rstrip()}\n{marker}"

    def rows(topic: str) -> list[dict[str, object]]:
        return [row for row in facts if row.get("topic") == topic]

    def rule(topic: str) -> str:
        return fence(
            "text",
            (rule_excerpts or {}).get(
                topic,
                "No bounded governing provision was extracted.",
            ),
        )

    def result(topic: str, fact_topics: tuple[str, ...]) -> str:
        """
        Render one finding summary with its normalized evidence count.

        Parameters
        ----------
        topic : str
            Finding topic to summarize.
        fact_topics : tuple[str, ...]
            Evidence topics whose records support a zero-finding result.

        Returns
        -------
        result : str
            Finding details or the evidence-backed zero-finding statement.
        """

        matched = [row for row in findings if row.get("topic") == topic]
        fact_count = sum(len(rows(name)) for name in fact_topics)

        if matched:
            grouped: dict[str, Counter[str]] = defaultdict(Counter)

            for finding in matched:
                grouped[str(finding["message"])][str(finding["path"])] += 1

            labels = {
                "runtime/version requirement must use 'bash >= 4.4'": (
                    "runtime/version capitalization findings"
                ),
                (
                    "submit help is missing Runtime requirements for bash >= "
                    "4.4, wget, and ln"
                ): "missing submit Runtime requirements block",
            }
            path_labels = {
                "bin/execute_download_fastqs.sh": "execute wrapper",
                "bin/submit_download_fastqs.sh": "submit wrapper",
                "lib/bash/help/help_execute_download_fastqs.sh": (
                    "execute help"
                ),
                "lib/bash/help/help_submit_download_fastqs.sh": "submit help",
            }

            summary = [f"{len(matched)} deterministic findings:"]

            for message, by_path in grouped.items():
                subtotal = sum(by_path.values())
                summary.append(f"- {subtotal} {labels.get(message, message)}:")
                summary.extend(
                    f"  - {path_labels.get(path, path)}: {count}"
                    for path, count in sorted(by_path.items())
                )

            return "\n".join(summary)

        return (
            f"No deterministic finding in this run; {fact_count} normalized "
            f"evidence record(s) support that result."
        )

    def topic(
        title: str,
        classification: str,
        result_text: str,
        evidence: str,
        uncertainty: str,
        question: str,
        paths: str,
        rule_name: str,
    ) -> str:
        return (
            f"## {title}\n\nCoverage classification: "
            f"`{classification}`.\nDeterministic result: "
            f"{result_text}\nEvidence "
            f"summary: {evidence}\nSemantic uncertainty: "
            f"{uncertainty}\nReviewer "
            f"question: {question}\nRelevant paths: {paths}\nGoverning "
            f"provision:\n{rule(rule_name)}\n"
        )

    def alias_table() -> str:
        """
        Render accepted and documented alias agreement as a Markdown table.

        Returns
        -------
        table : str
            One table covering every accepted or documented wrapper alias.
        """

        accepted: dict[str, dict[str, tuple[str, str]]] = defaultdict(dict)
        documented: dict[str, set[str]] = defaultdict(set)
        owners = {
            str(row["value"]["documentation_source"]): str(row["path"])
            for row in rows("documentation_source_associations")
        }

        for row in rows("resolved_alias_classifications"):
            for item in row["value"]:
                accepted[str(row["path"])][str(item["alias"])] = (
                    str(item["visibility"]),
                    ", ".join(
                        f"{loc['function']}:{loc['line']}"
                        for loc in item.get("locations", [])
                    ),
                )

        for row in rows("documented_aliases"):
            owner = owners.get(str(row["path"]))

            if owner is not None:
                documented[owner].update(str(alias) for alias in row["value"])

        output: list[list[object]] = []
        wrappers = sorted(set(accepted) | set(documented))

        for path in wrappers:
            for alias in sorted(set(accepted[path]) | documented[path]):
                classification, source = accepted[path].get(
                    alias,
                    ("documented_alias", "help parameter row"),
                )
                is_accepted, is_documented = (
                    alias in accepted[path],
                    alias in documented[path],
                )
                resolved_hidden = classification in {
                    "hidden_legacy_compatibility",
                    "hidden_systematic_compatibility",
                }

                if classification == "unsupported_fallback":
                    outcome = "parser fallback; not a supported option"
                elif resolved_hidden:
                    outcome = (
                        "hidden compatibility; intentionally undocumented"
                    )
                elif is_accepted and not is_documented:
                    outcome = "public accepted but undocumented"
                elif is_documented and not is_accepted:
                    outcome = "documented but not accepted"
                else:
                    outcome = "accepted and documented"

                if (
                    is_accepted
                    and is_documented
                    and alias not in {"-h", "--help"}
                ):
                    continue

                displayed_classification = (
                    "resolved_hidden_compatibility_alias"
                    if resolved_hidden
                    else classification
                )
                output.append(
                    [
                        f"`{path}`",
                        f"`{alias}`",
                        "yes" if is_accepted else "no",
                        "yes" if is_documented else "no",
                        f"`{displayed_classification}`",
                        source,
                        outcome,
                    ],
                )

        return markdown_table(
            [
                "Path",
                "Alias",
                "Accepted",
                "Documented",
                "Classification",
                "Source/interface phase",
                "Result",
            ],
            output,
        )

    directory = report_dir / "semantic_review"
    directory.mkdir(parents=True, exist_ok=True)
    rule_excerpts, diff_excerpts = rule_excerpts or {}, diff_excerpts or {}

    primary = [row for row in targets if row["role"] == "primary"]
    assignments = rows("interface_assignments")
    lines = rows("line_length_candidates")
    supporting = rows("supporting_test_alignment")
    controlled_source = rows("controlled_smoke_source_coverage")
    controlled_runtime = rows("controlled_smoke_execution_result")
    stale = rows("stale_name_occurrences")

    target_table = markdown_table(
        ["Path", "Role", "Git state", "Checks", "Findings"],
        [
            [
                f"`{row['path']}`",
                row["role"],
                ", ".join(row["git_state_labels"]),
                len(row.get("checks_run", [])),
                row.get("findings_count", 0),
            ]
            for row in targets
        ],
        [":---", ":---", ":---", "---:", "---:"],
    )

    assignment_count = sum(len(row["value"]) for row in assignments)
    stale_counts = {
        name: sum(len(row["value"].get(name, [])) for row in stale)
        for name in ("err_out", "infile", "outfile")
    }

    topics = [
        topic(
            "Help ownership",
            "fully_deterministic",
            result("help_ownership", ("help_ownership",)),
            "Static wrapper ownership facts were extracted.",
            (
                "The narrow ownership contract does not establish broader "
                "public-function visibility."
            ),
            (
                "Confirm the static ownership boundary matches the intended "
                "interface."
            ),
            "execute and submit wrappers",
            "help ownership",
        ),
        topic(
            "Aliases",
            "deterministic extraction plus semantic review",
            result(
                "aliases",
                (
                    "parser_alias_observations",
                    "documented_aliases",
                    "startup_help_alias_observations",
                    "resolved_alias_classifications",
                ),
            ),
            (
                "Raw parser/startup observations and resolved compatibility "
                "classifications are recorded separately."
            ),
            (
                "Dynamic parser behavior outside the recognized case-arm "
                "syntax remains semantic-only."
            ),
            (
                "Review only aliases not covered by the established public "
                "and hidden-compatibility policies."
            ),
            "wrappers and associated help",
            "aliases",
        ),
        topic(
            "Headings and required sections",
            "fully_deterministic for required wrapper Examples",
            result("examples", ("headings", "examples", "section_assessment")),
            (
                "Examples counts and normalized command shapes are "
                "machine-recorded."
            ),
            (
                "Whether example content is materially useful remains "
                "review-only."
            ),
            (
                "Do the examples reflect accepted interfaces and materially "
                "distinct uses?"
            ),
            "primary help modules",
            "headings and required sections",
        ),
        topic(
            "Types and placeholders",
            (
                "vocabulary/syntax fully deterministic; semantic "
                "appropriateness review-only"
            ),
            result("types_and_placeholders", ("help_parameter_rows",)),
            (
                "Parameter rows remain in facts.ndjson; only "
                "vocabulary/syntax is mechanically evaluated."
            ),
            "File/path/directory/string choice is domain-specific.",
            (
                "Are the selected types and placeholders conceptually "
                "appropriate?"
            ),
            "primary help modules",
            "types and placeholders",
        ),
        topic(
            "Requiredness and defaults",
            "deterministic extraction plus semantic review",
            result("requiredness_defaults", ("interface_assignments",)),
            (
                f"{assignment_count} path-scoped interface assignment "
                f"record(s) were extracted; details remain in"
                f" facts.ndjson."
            ),
            (
                "Conditional requiredness across parser, prose, and examples "
                "is not inferred."
            ),
            "Do documented requiredness and defaults match runtime behavior?",
            "primary wrappers and help",
            "requiredness/defaults",
        ),
        topic(
            "Runtime requirements",
            "registry-backed deterministic checks plus semantic review",
            result(
                "runtime_requirements",
                (
                    "command_reference_observations",
                    "command_reference_resolutions",
                    "section_assessment",
                ),
            ),
            (
                "Executable/version and External programs references are "
                "resolved against exact callable spellings."
            ),
            (
                "Unknown and ambiguous references remain semantic-review "
                "items; conceptual prose is excluded."
            ),
            "Review only unknown or ambiguous command references.",
            "primary scripts and help",
            "runtime requirements",
        ),
        topic(
            "--dir_scr",
            "fully_deterministic",
            result(
                "dir_scr",
                ("execute_submit_dir_scr", "submit_dir_scr_bootstrap"),
            ),
            (
                "Physical default, derived submit path, forwarding, and "
                "submit bootstrap are extracted separately."
            ),
            "This does not establish a general wrapper-topology verdict.",
            "Is this narrow forwarding contract sufficient for the workflow?",
            "execute and submit wrappers",
            "--dir_scr",
        ),
        topic(
            "Environment handling",
            "semantic review only",
            result(
                "environment_handling",
                ("interface_assignments", "supporting_test_alignment"),
            ),
            (
                f"{assignment_count} scoped assignment record(s) and "
                f"supporting-test evidence are"
                f" available."
            ),
            "Resolved environment propagation is behavior-sensitive.",
            (
                "Does the documented environment interface match resolution "
                "and propagation?"
            ),
            "execute wrapper and supporting tests",
            "environment handling",
        ),
        topic(
            "Stale naming",
            "deterministic extraction plus semantic review",
            result("stale_naming", ("stale_name_occurrences",)),
            (
                f"Observed counts: `err_out` {stale_counts['err_out']}; "
                f"`infile` "
                f"{stale_counts['infile']}; `outfile` "
                f"{stale_counts['outfile']}."
            ),
            "Occurrences may be internal or intentional.",
            "Are any occurrences user-facing or behaviorally relevant?",
            "selected targets",
            "stale naming",
        ),
        topic(
            "Behavior-sensitive diff review",
            "semantic review only",
            "Not applicable; this topic is semantic-review only.",
            (
                f"{len(diff_excerpts)} bounded tracked-primary diff "
                f"excerpt(s) are available below."
            ),
            "No behavior conclusion is inferred from a diff excerpt.",
            (
                "Do the focused changes preserve command construction and "
                "startup behavior?"
            ),
            "tracked primary targets",
            "behavior-sensitive diff review",
        ),
    ]

    unresolved = list(policy_questions)
    semantic_decisions = [
        (
            "Review unresolved alias intent not covered by established "
            "hidden-compatibility policy; review type appropriateness, "
            "requiredness/default semantics, ordinary-code line lengths, "
            "unknown or ambiguous command references, and focused "
            "behavior-sensitive diffs."
        ),
    ]

    decisions = [
        str(row.get("question")) for row in unresolved
    ] + semantic_decisions

    support_rows: list[list[object]] = []

    for row in supporting:
        value = row["value"]
        support_rows.append(
            [
                f"`{str(row['path']).split('/')[-1]}`",
                value["registration_status"],
                value["tested_mode"],
                ", ".join(value["public_options_invoked"]) or "none",
                ", ".join(value["required_arguments_represented"]) or "none",
                value["relevant_interface_coverage"],
                value["apparent_gap_or_uncertainty"],
            ],
        )

    support_table = markdown_table(
        [
            "Test",
            "Registered",
            "Mode",
            "Public options",
            "Required options",
            "Coverage",
            "Gap/uncertainty",
        ],
        support_rows,
    )

    source_evidence = (
        fence(
            "json",
            json.dumps(
                [row["value"] for row in controlled_source],
                indent=2,
                sort_keys=True,
            ),
        )
        if controlled_source
        else "No declarative controlled-smoke source coverage was selected."
    )

    runtime_evidence = (
        fence(
            "json",
            json.dumps(
                [row["value"] for row in controlled_runtime],
                indent=2,
                sort_keys=True,
            ),
        )
        if controlled_runtime
        else "No normalized controlled-smoke execution result was supplied."
    )

    actionable, excluded = [], defaultdict(int)

    for row in lines:
        for candidate in row["value"]:
            relationship = candidate.get("source_checker_relationship")

            if relationship == "excluded_by_documented_source_checker_policy":
                excluded[
                    (
                        str(row["path"]),
                        str(candidate.get("location_kind", "unknown")),
                    )
                ] += 1
            else:
                actionable.append((row, candidate))

    grouped_lines: dict[tuple[str, str, str], list[dict[str, object]]] = (
        defaultdict(list)
    )
    individual_lines: list[tuple[dict[str, object], dict[str, object]]] = []

    for row, candidate in actionable:
        excerpt_lines = str(candidate["excerpt"]).splitlines()
        candidate_line = next(
            (
                line
                for line in excerpt_lines
                if line.startswith(f"{candidate['line']}:")
            ),
            "",
        )
        is_log_assignment = str(row["path"]).startswith(
            "tests/smoke/",
        ) and bool(
            re.match(r"\d+:\s+log_(?:out|err)_[A-Za-z0-9_]+=", candidate_line),
        )

        if is_log_assignment:
            grouped_lines[
                (
                    str(row["path"]),
                    str(candidate.get("classification_status")),
                    "log-path assignment",
                )
            ].append(candidate)
        else:
            individual_lines.append((row, candidate))

    line_blocks = [
        (
            f"`{row['path']}`:{candidate['line']} — {candidate['length']} "
            f"characters; "
            f"`{candidate.get('location_kind')}`; "
            f"`{candidate.get('source_checker_relationship')}`; cues: "
            f"{', '.join(candidate.get('exception_cues', [])) or 'none'}; "
            f"status: "
            f"`{candidate.get('classification_status')}`.\n\n"
            f"{fence('bash', candidate['excerpt'])}"
        )
        for row, candidate in individual_lines
    ]

    grouped_summary = []

    for (path, status, pattern), candidates in sorted(grouped_lines.items()):
        numbers = sorted(int(item["line"]) for item in candidates)
        lengths = sorted(int(item["length"]) for item in candidates)
        grouped_summary.append(
            (
                f"- `{path}`: lines {numbers[0]}–{numbers[-1]}; "
                f"{len(candidates)} {pattern} "
                f"candidate(s); {lengths[0]}–{lengths[-1]} characters; "
                f"status: `{status}`. "
                f"Details remain in "
                f"`facts.ndjson`."
            ),
        )

    exclusion_table = markdown_table(
        ["Path", "Location class", "Count", "Status"],
        [
            [f"`{path}`", f"`{kind}`", count, "excluded_not_review_required"]
            for (path, kind), count in sorted(excluded.items())
        ],
        [":---", ":---", "---:", ":---"],
    )

    diffs = [
        f"### `{path}`\n\n{fence('diff', excerpt)}"
        for path, excerpt in diff_excerpts.items()
        if excerpt
    ]

    if anchor_evidence is not None:
        missing_paths = {
            str(target["path"])
            for target in primary
            if str(target["path"]) in REQUIRED_BEHAVIORAL_ANCHORS
        } - set(anchor_evidence)
        if missing_paths:
            raise ValueError(
                "missing required behavioral anchor evidence for: "
                + ", ".join(sorted(missing_paths)),
            )

        for target in primary:
            path = str(target["path"])
            required_groups = set(REQUIRED_BEHAVIORAL_ANCHORS.get(path, {}))
            missing_groups = required_groups - set(
                anchor_evidence.get(path, {}),
            )
            if missing_groups:
                raise ValueError(
                    (
                        f"missing required behavioral anchor group "
                        f"{sorted(missing_groups)[0]}: "
                        f"{path}"
                    ),
                )

    anchor_rows: list[list[object]] = []
    rendered_anchor_evidence: list[str] = []
    anchor_metadata: list[dict[str, object]] = []

    for path, groups in (anchor_evidence or {}).items():
        for group, record in groups.items():
            terms = REQUIRED_BEHAVIORAL_ANCHORS.get(path, {}).get(group)
            if terms is None:
                raise ValueError(
                    f"unknown behavioral anchor group {group}: {path}",
                )

            evidence = record.get("rendered_evidence", "")

            if not isinstance(evidence, str):
                raise ValueError(
                    (
                        f"missing required behavioral anchor group {group}: "
                        f"{path}"
                    ),
                )

            evidence_kind = record.get(
                "evidence_kind",
                "complete_line_source_window",
            )

            if not isinstance(evidence_kind, str):
                raise ValueError(
                    f"invalid behavioral anchor evidence kind {group}: {path}",
                )

            if (
                (path, group) in SOURCE_ONLY_ANCHORS
                and evidence_kind != "complete_line_source_window"
            ):
                raise ValueError(
                    f"anchor requires current source evidence {group}: {path}",
                )

            if not all(
                term in _current_evidence_text(evidence, evidence_kind)
                for term in terms
            ):
                raise ValueError(
                    (
                        f"missing required behavioral anchor group {group}: "
                        f"{path}"
                    ),
                )

            anchor_rows.append([f"`{path}`", group, evidence_kind])
            anchor_metadata.append(
                {
                    "path": path,
                    "anchor_group": group,
                    "expected_source_markers": list(terms),
                    "evidence_kind": evidence_kind,
                    "rendered_evidence": evidence,
                    "rendered_evidence_fingerprint": _sha256_text(evidence),
                },
            )

            evidence_language = (
                "diff"
                if evidence_kind == "complete_line_diff_window"
                else "bash"
            )
            rendered_block = fence(evidence_language, evidence)
            rendered_anchor_evidence.append(
                f"#### `{path}` — {group}\n\n"
                f"evidence_kind: `{evidence_kind}`\n\n"
                f"{rendered_block}",
            )

    anchor_table = markdown_table(
        ["Path", "Anchor group", "Evidence status"],
        anchor_rows,
    )

    untracked = []

    for target in primary:
        if "untracked" in target["git_state_labels"]:
            source = source_excerpt(str(target["path"]))
            match = re.search(
                r"function\s+help_[A-Za-z0-9_]+\(\)\s*\{\s*cat"
                r"\s+<<\s*EOM\n(.*?)\nEOM",
                source,
                re.DOTALL,
            )

            if match:
                evidence_kind, excerpt = (
                    "new_file_help_heredoc",
                    match.group(1),
                )
            else:
                evidence_kind, excerpt = (
                    "new_file_source_fallback",
                    "\n".join(source.splitlines()[:80]),
                )

            untracked.append(
                (
                    f"### `{target['path']}`\n\ngit_state: untracked  "
                    f"\nevidence_kind: "
                    f"{evidence_kind}  \nfingerprint: "
                    f"`{target['content_fingerprint']}`\n\n"
                    f"{fence('bash', excerpt)}"
                ),
            )

    findings_table = markdown_table(
        ["Topic", "Path", "Finding"],
        [
            [
                f"`{row.get('topic', 'other')}`",
                f"`{row['path']}`",
                row["message"],
            ]
            for row in findings
        ],
    )

    body = "# Download FASTQs shell/help pilot: semantic review\n\n"
    body += f"Run ID: `{run_id or 'unspecified'}`\n\n"

    body += (
        "## Executive summary\n\nReview-only bounded pilot. This rendering "
        "summarizes versioned machine artifacts; no patch, semantic "
        "disposition, or production action is to be applied in this review. "
        "Fully deterministic describes the extractor/check path; medium "
        "confidence records that pilot-adapter scope remains subject to "
        "independent false-positive review.\n\n"
    )

    body += (
        "## Immediate reviewer decisions\n\n"
        + "\n".join(
            f"{index}. {item}" for index, item in enumerate(decisions, 1)
        )
        + "\n\n"
    )

    body += "## Target status\n\n" + target_table + "\n\n"
    body += "\n".join(topics) + "\n"
    body += "## Alias comparison\n\n" + alias_table() + "\n\n"
    body += "## Supporting-test alignment\n\n" + support_table + "\n\n"

    body += (
        "## Phase 2 direct-submit evidence\n\n### Source-derived coverage\n\n"
        + source_evidence
        + "\n\n### Normalized execution-result evidence\n\n"
        + runtime_evidence
        + (
            "\n\n### Remaining runtime limitation\n\nThe normalized result "
            "records only directly captured outer-shell facts and configured "
            "assertions. Environment selection is classified separately; "
            "broader host-versus-Conda inheritance and complete binary "
            "provenance remain review-only.\n\n"
        )
    )

    body += "## Deterministic findings\n\n" + findings_table + "\n\n"

    body += (
        "## Semantic-review items\n\n"
        + "\n".join(f"- {item}" for item in semantic_decisions)
        + "\n\n"
    )

    body += (
        "## Policy questions\n\n"
        + (
            fence("json", json.dumps(unresolved, indent=2, sort_keys=True))
            if unresolved
            else "No unresolved policy questions."
        )
        + "\n\n"
    )

    body += (
        "## Line-length candidates\n\n"
        + (
            "\n\n".join(line_blocks)
            or "No individually rendered ordinary-code candidates."
        )
        + "\n\n"
    )

    body += (
        "### Grouped repetitive candidates\n\n"
        + (
            "\n".join(grouped_summary)
            or "No repetitive candidates were grouped."
        )
        + "\n\n### Documented exclusions\n\n"
        + exclusion_table
        + "\n\n"
    )

    body += (
        "## Behavior-sensitive diff excerpts\n\n"
        + (
            "\n\n".join(diffs)
            or "No focused tracked-primary diff hunk was available."
        )
        + "\n\n"
    )

    body += "### Behavioral anchor coverage\n\n" + anchor_table + "\n\n"

    body += (
        "### Behavioral anchor evidence\n\n"
        + (
            "\n\n".join(rendered_anchor_evidence)
            or "No anchor-specific evidence was requested."
        )
        + "\n\n"
    )

    if untracked:
        body += (
            "## Untracked primary evidence\n\n"
            + "\n\n".join(untracked)
            + "\n\n"
        )

    body += (
        "## Adapter limitations and false-positive risks\n\nDynamic shell "
        "construction, conceptual documentation quality, and "
        "behavior-sensitive command construction remain review-only. "
        "False-positive assessment remains pending independent review.\n\n"
    )

    body += (
        (
            "## Structured response template\n\nDisposition must be "
            "`true_positive`, `false_positive`, or `uncertain`.\n\n"
        )
        + fence(
            "json",
            json.dumps(
                {
                    "finding_id": f"{run_id or 'run'}:...",
                    "topic": "aliases",
                    "path": "bin/submit_download_fastqs.sh",
                    "line": None,
                    "disposition": "uncertain",
                    "rationale": "",
                    "content_fingerprint": "sha256:...",
                    "recommended_action": None,
                    "adapter_implication": None,
                },
                indent=2,
            ),
        )
        + "\n\n"
    )

    body += (
        "## Supplied review package\n\n"
        + "\n".join(
            f"- `{name}`"
            for name in (
                (
                    "semantic_review/download-fastqs-shell-help-pilot.md "
                    "(this file)"
                ),
                "findings.ndjson",
                "facts.ndjson",
                "adapter_limitations.ndjson",
                "pilot_report.json",
            )
        )
        + (
            "\n\nOther runner files are generator-side provenance artifacts "
            "and are intentionally not supplied for this review.\n"
        )
    )

    table_count = validate_bundle_tables(body)
    line_count = len(body.splitlines())
    byte_count = len(body.encode("utf-8"))

    within_line_target = 300 <= line_count <= 500
    target_reason = (
        None
        if within_line_target
        else (
            "Reviewer-critical bounded source, diff, and individual "
            "production-code evidence exceeds the preferred line target."
        )
    )

    fingerprint = _sha256_text(body)
    path = directory / "download-fastqs-shell-help-pilot.md"
    path.write_text(body, encoding="utf-8")

    return [
        {
            "schema_version": 1,
            "run_id": run_id,
            "bundle_id": "download-fastqs-shell-help-pilot",
            "path": "semantic_review/download-fastqs-shell-help-pilot.md",
            "primary_paths": [row["path"] for row in primary],
            "target_fingerprints": {
                row["path"]: row["content_fingerprint"] for row in targets
            },
            "prompt_fingerprint": fingerprint,
            "markdown_lines": line_count,
            "markdown_bytes": byte_count,
            "markdown_tables_validated": table_count,
            "behavioral_anchor_coverage": anchor_metadata,
            "line_length_rendering": {
                "detailed_candidate_count": len(actionable),
                "human_rendered_individual_candidate_count": len(
                    individual_lines,
                ),
                "grouped_candidate_count": len(grouped_lines),
            },
            "markdown_target": {
                "minimum_lines": 300,
                "maximum_lines": 500,
                "met": within_line_target,
                "reviewer_critical_reason": target_reason,
            },
        },
    ]
