"""Normalize output from strict cleanup-audit check adapters."""

from __future__ import annotations

import json
import re
from dataclasses import dataclass


@dataclass(frozen=True)
class CommandResult:
    """Captured subprocess result supplied to an output parser."""

    argv: list[str]
    stdout: str
    stderr: str
    returncode: int
    duration_seconds: float = 0.0
    timed_out: bool = False
    launch_error: str | None = None


def _line(value: str | None) -> int | None:
    """Return a positive line number when the checker supplied one."""

    if value is None:
        return None
    number = int(value)
    return number if number > 0 else None


def parse_result(
    parser_name: str,
    result: CommandResult,
    target_path: str,
) -> list[dict[str, object]]:
    """Return normalized parser fragments without audit-run metadata."""

    text = "\n".join(part for part in (result.stdout, result.stderr) if part)
    fragments: list[dict[str, object]] = []
    if result.timed_out or result.launch_error:
        return fragments
    if parser_name == "bash_stderr":
        pattern = re.compile(r"^(.*?): line ([0-9]+): (.+)$", re.MULTILINE)
        for match in pattern.finditer(text):
            fragments.append(
                {
                    "path": match.group(1),
                    "line": _line(match.group(2)),
                    "message": match.group(3).strip(),
                    "evidence": match.group(0),
                }
            )
    elif parser_name == "python_compile_stderr":
        matches = list(re.finditer(r'File "([^"]+)", line ([0-9]+)', text))
        for index, match in enumerate(matches):
            block = text[match.start():matches[index + 1].start() if index + 1 < len(matches) else len(text)]
            message = next((line.strip() for line in reversed(block.splitlines()) if line.strip()), "Python compilation failed")
            fragments.append(
                {
                    "path": match.group(1),
                    "line": _line(match.group(2)),
                    "message": message,
                    "evidence": block.strip(),
                }
            )
        if not fragments:
            fragments.append(
                {
                    "path": target_path,
                    "line": None,
                    "message": next((line.strip() for line in reversed(text.splitlines()) if line.strip()), "Python compilation failed"),
                    "evidence": text.strip(),
                }
            )
    elif parser_name == "git_diff_check":
        for line in text.splitlines():
            match = re.match(r"^(.+?):([0-9]+):\s*(.*)$", line)
            fragments.append(
                {
                    "path": match.group(1) if match else ".",
                    "line": _line(match.group(2)) if match else None,
                    "message": match.group(3) if match else line,
                    "evidence": line,
                }
            )
    if not fragments and result.returncode != 0:
        fragments.append(
            {
                "path": target_path,
                "line": None,
                "message": "Checker failed without parseable output",
                "evidence": text.strip(),
            }
        )
    return fragments


def parse_pilot_payload(result: CommandResult) -> dict[str, list[dict[str, object]]]:
    """Validate the structured output from a snapshot-only pilot adapter."""

    if result.timed_out or result.launch_error:
        return {"facts": [], "findings": [], "policy_questions": [], "limitations": []}
    try:
        payload = json.loads(result.stdout)
    except json.JSONDecodeError as exc:
        raise ValueError(f"pilot adapter emitted invalid JSON: {exc.msg}") from exc
    if not isinstance(payload, dict) or payload.get("schema_version") != 1:
        raise ValueError("pilot adapter payload must be a schema_version 1 object")
    parsed: dict[str, list[dict[str, object]]] = {}
    for key in ("facts", "findings", "policy_questions", "limitations"):
        value = payload.get(key, [])
        if not isinstance(value, list) or not all(isinstance(item, dict) for item in value):
            raise ValueError(f"pilot adapter payload {key} must be a list of objects")
        parsed[key] = value
    return parsed
