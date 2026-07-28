#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: shell_validation.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


"""
Resolve and describe the repository's authoritative shell syntax gates.
"""

from __future__ import annotations

import dataclasses
import re
import subprocess
from pathlib import Path

MINIMUM_BASH = (4, 4)
POSIX_BOOTSTRAP = "install/scripts/install_envs_entrypoint.sh"


@dataclasses.dataclass(frozen=True)
class BashInterpreter:
    """
    One absolute Bash path and parsed version from env_protocol.
    """

    path: Path
    version: str
    version_tuple: tuple[int, int, int]

    @property
    def satisfies_minimum(self) -> bool:
        """
        Return whether the interpreter satisfies the repository minimum.
        """

        return self.version_tuple[:2] >= MINIMUM_BASH

    @property
    def label(self) -> str:
        """
        Return a diagnostic label containing path, version, and status.
        """

        status = "satisfies" if self.satisfies_minimum else "does not satisfy"
        return (
            f"authoritative Bash: {self.path} (GNU Bash {self.version}; "
            f"{status} >= {MINIMUM_BASH[0]}.{MINIMUM_BASH[1]})"
        )


def resolve_env_bash(environment: str = "env_protocol") -> BashInterpreter:
    """
    Resolve Bash inside one Conda environment without trusting caller PATH.
    """

    result = subprocess.run(
        [
            "conda",
            "run",
            "--no-capture-output",
            "-n",
            environment,
            "bash",
            "-c",
            'command -v bash; printf "%s\\n" "${BASH_VERSION}"',
        ],
        text=True,
        capture_output=True,
        check=True,
    )
    rows = [row.strip() for row in result.stdout.splitlines() if row.strip()]
    if len(rows) < 2:
        raise RuntimeError(
            "could not resolve Bash path and version from env_protocol",
        )

    path = Path(rows[0]).resolve()
    match = re.match(
        r"(?P<major>\d+)\.(?P<minor>\d+)\.(?P<patch>\d+)",
        rows[1],
    )
    if not path.is_absolute() or match is None:
        raise RuntimeError(
            "env_protocol returned an invalid Bash path or version",
        )

    return BashInterpreter(
        path=path,
        version=rows[1],
        version_tuple=tuple(
            int(match.group(name)) for name in ("major", "minor", "patch")
        ),
    )


def syntax_command(path: str, interpreter: BashInterpreter) -> list[str]:
    """
    Return the authoritative syntax command for one repository shell source.
    """

    executable = "sh" if path == POSIX_BOOTSTRAP else str(interpreter.path)
    return [executable, "-n", path]
