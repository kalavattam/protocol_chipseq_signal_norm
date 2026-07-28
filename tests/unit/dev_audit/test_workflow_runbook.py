#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: test_workflow_runbook.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


"""
Structural regressions for maintained workflow runbook fences.
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
WORKFLOW = ROOT / "notebooks/workflows/workflow.md"
SECTION = re.compile(r"^### (?P<label>[H-J])\.", re.MULTILINE)
BASH_FENCE = re.compile(r"```bash\n(?P<body>.*?)\n```", re.DOTALL)


def ratio_fences() -> dict[str, str]:
    """
    Return the first Bash fence from each ratio workflow section.
    """

    text = WORKFLOW.read_text(encoding="utf-8")
    matches = list(SECTION.finditer(text))
    output: dict[str, str] = {}

    for index, match in enumerate(matches):
        end = matches[index + 1].start() if index + 1 < len(matches) else None
        fence = BASH_FENCE.search(text[match.end() : end])

        assert fence is not None

        output[match.group("label")] = fence.group("body")

    return output


def workflow_fence(label: str) -> str:
    """
    Return the first Bash fence from one lettered workflow section.
    """

    text = WORKFLOW.read_text(encoding="utf-8")
    heading = (
        rf"^### {re.escape(label)}\. Obtain and organize "
        r"ChIP-seq FASTQ files\.$"
    )

    start = re.search(heading, text, re.MULTILINE)

    assert start is not None

    end = re.search(r"^### [A-Z][0-9]?\.", text[start.end() :], re.MULTILINE)
    section_end = start.end() + end.start() if end is not None else None
    fence = BASH_FENCE.search(text[start.end() : section_end])

    assert fence is not None

    return fence.group("body")


def test_download_fastqs_section_uses_current_interface() -> None:
    """
    Keep Data C on maintained paths, helpers, and one command array.
    """

    body = workflow_fence("C")
    required = [
        'dir_bin="${dir_rep}/bin"',
        'dir_bash="${dir_rep}/lib/bash"',
        'source "${dir_bash}/core/source_helpers.sh"',
        'source_helpers "${dir_bash}" check_env format_outputs handle_env',
        'handle_env "${env_nam}"',
        "check_pgrm_path parallel",
        "check_pgrm_path wget",
        'bash "${dir_bin}/execute_download_fastqs.sh"',
        '--fil_in "${pth_tsv}"',
        '--dir_eo "${dir_log}"',
        'cmd_download+=( --slurm --time "${time}" )',
        "printf '%q ' \"${cmd_download[@]}\"",
        '"${cmd_download[@]}"',
        'bash "${dir_bin}/compress_remove_files.sh"',
    ]
    retired = [
        "${dir_rep}/scripts",
        "check_program_path",
        "echo_warning",
        "--infile",
        "--err_out",
        "${dir_scr}/execute_download_fastqs.sh",
        "${dir_scr}/compress_remove_files.sh",
    ]

    for value in required:
        assert value in body

    for value in retired:
        assert value not in body

    assert body.count('"${cmd_download[@]}"') == 2


def test_download_fastqs_section_fails_closed() -> None:
    """
    Require successful download before independently guarded cleanup.
    """

    body = workflow_fence("C")
    download_guard = body.index('if ! \\\n    "${cmd_download[@]}" \\')
    download_error = body.index(
        'echo "error(workflow): download_fastqs failed; cleanup was not run." '
        ">&2",
        download_guard,
    )
    download_exit = body.index("    exit 1", download_error)
    download_end = body.index("\nfi", download_exit)

    cleanup_guard = body.index(
        'if ! \\\n    bash "${dir_bin}/compress_remove_files.sh"',
        download_end,
    )
    cleanup_error = body.index(
        'echo "error(workflow): log cleanup failed." >&2',
        cleanup_guard,
    )
    cleanup_exit = body.index("    exit 1", cleanup_error)
    cleanup_end = body.index("\nfi", cleanup_exit)

    assert [
        download_guard,
        download_error,
        download_exit,
        download_end,
        cleanup_guard,
        cleanup_error,
        cleanup_exit,
        cleanup_end,
    ] == sorted(
        [
            download_guard,
            download_error,
            download_exit,
            download_end,
            cleanup_guard,
            cleanup_error,
            cleanup_exit,
            cleanup_end,
        ],
    )


def test_download_fastqs_section_is_valid_bash() -> None:
    """
    Keep the final Data C fence valid Bash.
    """

    body = workflow_fence("C")
    subprocess.run(
        ["bash", "-n"],
        input=body,
        text=True,
        check=True,
        capture_output=True,
    )


@pytest.mark.parametrize("label", ["H", "I", "J"])
def test_ratio_environment_precedes_installed_command_use(
    label: str,
) -> None:
    """
    Require one H/I/J activation before table and floor commands.
    """

    body = ratio_fences()[label]
    source = 'source "${dir_bash}/core/source_helpers.sh"'
    tables = 'source_helpers "${dir_bash}" workflows/process_tables'
    helpers = (
        'source_helpers "${dir_bash}" check_env format_outputs handle_env'
    )
    activate = 'handle_env "${env_nam}"'
    verify = "check_pgrm_path compute_input_floor"
    extract = "extract_field_str"
    invoke = re.search(r"(?m)^    compute_input_floor \\$", body)

    assert body.count(source) == 1
    assert body.count(tables) == 1
    assert body.count(helpers) == 1
    assert body.count(activate) == 1
    assert body.count(verify) == 1
    assert invoke is not None

    ordered = [
        body.index(source),
        body.index(tables),
        body.index(helpers),
        body.index(activate),
        body.index(verify),
        body.index(extract),
        invoke.start(),
    ]

    assert ordered == sorted(ordered)
