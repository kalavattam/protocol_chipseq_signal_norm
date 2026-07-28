#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script: coordinator.py
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


"""
Prepare, transfer, run, and verify portable Slurm validation bundles.

The public commands operate on an isolated, run-specific directory. The
'remote-launch' and 'wet-run' commands are internal handoff interfaces used by
the prepared launcher and 'tests/integration/slurm/run_wet_tests.sh'.

Notes
-----
Python >= 3.11 is required. Network and scheduler actions occur only for an
explicit 'push', 'pull', remote 'status', 'remote-launch', or 'wet-run'
command. Preparing and inspecting a bundle are local operations.

Examples
--------
Prepare and inspect a fixed run without contacting a remote host::

    python tests/integration/slurm/coordinator.py prepare \
        --run_id wet-20260716 \
        --ssh_host rhino --partition campus-new --account my-account
    python tests/integration/slurm/coordinator.py status --run_id wet-20260716

Render transfer and remote commands without running them::

    python tests/integration/slurm/coordinator.py push \
        --run_id wet-20260716 --dry_run
    python tests/integration/slurm/coordinator.py instructions \
        --run_id wet-20260716
"""

from __future__ import annotations

import argparse
import dataclasses
import datetime as dt
import gzip
import hashlib
import json
import os
import platform
import re
import shlex
import shutil
import stat
import subprocess
import sys
import tarfile
import tempfile
import time
from collections.abc import Iterable, Sequence
from pathlib import Path, PurePosixPath
from typing import Any

from protocol_chipseq_signal_norm.utilities.utils_cli import (
    CapArgumentParser,
    add_help_cap,
)

RUNNER_VERSION = "2.0.0"
FORMAT_VERSION = 1
DEFAULT_BUNDLE_DIR = Path("artifacts/slurm/sessions")
DEFAULT_RESULT_DIR = Path("artifacts/slurm/results")
DEFAULT_REMOTE_ROOT = PurePosixPath("protocol_chipseq_signal_norm_slurm_runs")
DEFAULT_SCOPE = (
    "submit bootstrap plus tiny Bowtie2 alignment and BAM signal workers"
)
RUN_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$")
REMOTE_TOKEN_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
TERMINAL_STATES = {
    "BOOT_FAIL",
    "CANCELLED",
    "COMPLETED",
    "DEADLINE",
    "FAILED",
    "NODE_FAIL",
    "OUT_OF_MEMORY",
    "PREEMPTED",
    "REVOKED",
    "TIMEOUT",
}
REQUIRED_RESULT_PATHS = (
    "preflight.json",
    "run_manifest.json",
    "jobs.json",
    "commands.log",
    "summary.json",
    "summary.txt",
    "exit_status.txt",
    "checksums.sha256",
    "stdout",
    "stderr",
    "artifacts",
)
REQUIRED_JOB_KEYS = {
    "assertions",
    "cleanup_result",
    "command",
    "exit_code",
    "final_state",
    "finish_timestamp",
    "job_id",
    "job_name",
    "requested_resources",
    "start_timestamp",
    "stderr_path",
    "stdout_path",
    "submit_timestamp",
}
EXCLUDED_PARTS = {"__pycache__", ".DS_Store"}
EXCLUDED_PREFIXES = (PurePosixPath("artifacts"),)
REQUIRED_SOURCE_PREFIXES = (
    PurePosixPath("bin"),
    PurePosixPath("lib/bash"),
    PurePosixPath("src/protocol_chipseq_signal_norm"),
    PurePosixPath("tests/integration/slurm"),
    PurePosixPath("tests/fixtures/slurm"),
)
REQUIRED_SOURCE_PATHS = {
    PurePosixPath("install/envs/env_protocol.yml"),
}


class BundleError(RuntimeError):
    """
    Report a user-facing bundle or result validation failure.
    """


@dataclasses.dataclass(frozen=True)
class CommandResult:
    """
    Capture one subprocess command and its decoded output.
    """

    command: tuple[str, ...]
    returncode: int
    stdout: str
    stderr: str


def utc_now() -> str:
    """
    Return a UTC timestamp in a stable ISO-8601 representation.
    """

    return dt.datetime.now(dt.UTC).isoformat().replace("+00:00", "Z")


def canonical_json(value: Any) -> bytes:
    """
    Serialize JSON with deterministic key order and terminal newline.
    """

    return (
        json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    ).encode("utf-8")


def sha256_bytes(data: bytes) -> str:
    """
    Return the lowercase SHA-256 digest for bytes.
    """

    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    """
    Return a streaming SHA-256 digest for a regular file.
    """

    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


def write_json(path: Path, value: Any) -> None:
    """
    Write canonical JSON, creating the parent directory.
    """

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(canonical_json(value))


def load_json(path: Path) -> Any:
    """
    Load a JSON file and convert parse failures to 'BundleError'.
    """

    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise BundleError(f"malformed JSON file: {path}: {exc}") from exc


def validate_run_id(run_id: str) -> str:
    """
    Validate and return a portable run identifier.
    """

    if not RUN_ID_RE.fullmatch(run_id) or ".." in run_id:
        raise BundleError(
            "run ID must be 1-64 portable characters, start with an "
            "alphanumeric character, and contain no '..'",
        )

    return run_id


def validate_remote_token(
    value: str,
    label: str,
    *,
    empty_ok: bool = False,
) -> str:
    """
    Validate an SSH, partition, account, or environment token.
    """

    if empty_ok and not value:
        return value

    if not REMOTE_TOKEN_RE.fullmatch(value):
        raise BundleError(f"invalid {label}: {value!r}")

    return value


def safe_remote_run_dir(remote_root: str, run_id: str) -> PurePosixPath:
    """
    Construct a confined remote run directory without traversal.

    Parameters
    ----------
    remote_root : str
        Portable non-root POSIX directory selected for remote sessions.
    run_id : str
        Validated portable run identity appended below the root.

    Returns
    -------
    run_dir : PurePosixPath
        Confined remote path for the selected run.

    Raises
    ------
    BundleError
        If either component permits traversal or nonportable path syntax.
    """

    validate_run_id(run_id)

    root = PurePosixPath(remote_root)
    if not re.fullmatch(r"/?[A-Za-z0-9._/-]+", remote_root):
        raise BundleError("remote run root contains nonportable characters")

    if not remote_root or root == PurePosixPath("/") or ".." in root.parts:
        raise BundleError(
            "remote run root must be a non-root path without '..'",
        )

    if any(part in {"", "."} for part in root.parts[1:]):
        raise BundleError("remote run root contains an unsafe component")

    return root / run_id


def repo_root_from_script() -> Path:
    """
    Resolve the repository root from this script's installed location.
    """

    return Path(__file__).resolve().parents[3]


def validate_remote_python(value: str) -> str:
    """
    Validate an explicit absolute POSIX Python executable path.
    """

    path = PurePosixPath(value)
    if not path.is_absolute() or ".." in path.parts:
        raise BundleError(
            "remote Python must be an absolute path without '..'",
        )

    if not re.fullmatch(r"/[A-Za-z0-9._/-]+", value):
        raise BundleError("remote Python contains nonportable characters")

    return value


def run_command(
    command: Sequence[str],
    *,
    cwd: Path | None = None,
    env: dict[str, str] | None = None,
    timeout: float | None = None,
) -> CommandResult:
    """
    Run a command without a shell and capture text output.
    """

    completed = subprocess.run(
        list(command),
        cwd=cwd,
        env=env,
        text=True,
        capture_output=True,
        timeout=timeout,
        check=False,
    )

    return CommandResult(
        tuple(str(item) for item in command),
        completed.returncode,
        completed.stdout,
        completed.stderr,
    )


def git_bytes(root: Path, args: Sequence[str]) -> bytes:
    """
    Return exact bytes from a required Git command.
    """

    completed = subprocess.run(
        ["git", *args],
        cwd=root,
        capture_output=True,
        check=False,
    )

    if completed.returncode != 0:
        message = completed.stderr.decode("utf-8", errors="replace").strip()

        raise BundleError(f"git {' '.join(args)} failed: {message}")

    return completed.stdout


def excluded_source_path(relative: PurePosixPath) -> bool:
    """
    Return whether a repository path is excluded from source bundles.
    """

    if relative.is_absolute() or ".." in relative.parts:
        raise BundleError(f"unsafe repository path: {relative}")

    if any(part in EXCLUDED_PARTS for part in relative.parts):
        return True

    if relative.suffix == ".pyc":
        return True

    return any(
        relative == prefix or relative.is_relative_to(prefix)
        for prefix in EXCLUDED_PREFIXES
    )


def repository_paths(root: Path) -> list[PurePosixPath]:
    """
    Inventory tracked and required untracked, nonignored repository files.
    """

    raw = git_bytes(root, ["ls-files", "-co", "--exclude-standard", "-z"])
    values: list[PurePosixPath] = []

    for item in raw.split(b"\0"):
        if not item:
            continue

        try:
            relative = PurePosixPath(item.decode("utf-8"))
        except UnicodeDecodeError as exc:
            raise BundleError("repository paths must be UTF-8") from exc

        if excluded_source_path(relative):
            continue

        required = relative in REQUIRED_SOURCE_PATHS or any(
            relative == prefix or relative.is_relative_to(prefix)
            for prefix in REQUIRED_SOURCE_PREFIXES
        )
        if not required:
            continue

        absolute = root.joinpath(*relative.parts)

        if absolute.is_file() or absolute.is_symlink():
            values.append(relative)

    return sorted(set(values), key=lambda item: item.as_posix())


def source_entry(
    root: Path,
    relative: PurePosixPath,
    tracked: set[PurePosixPath],
) -> dict[str, Any]:
    """
    Build one manifest entry for a regular file or symbolic link.
    """

    path = root.joinpath(*relative.parts)
    info = path.lstat()
    mode = stat.S_IMODE(info.st_mode)

    if path.is_symlink():
        target = os.readlink(path)
        digest = sha256_bytes(target.encode("utf-8"))
        kind = "symlink"
        size = len(target.encode("utf-8"))
        link_target: str | None = target
    elif path.is_file():
        digest = sha256_file(path)
        kind = "file"
        size = info.st_size
        link_target = None
    else:
        raise BundleError(f"unsupported source entry: {relative}")

    return {
        "path": relative.as_posix(),
        "type": kind,
        "mode": f"{mode:04o}",
        "size": size,
        "sha256": digest,
        "link_target": link_target,
        "git_state": "tracked"
        if relative in tracked
        else "required_untracked",
    }


def build_source_manifest(
    root: Path,
    run_id: str,
    scope: str,
    created_at: str,
) -> dict[str, Any]:
    """
    Build the complete deterministic source manifest for one run.
    """

    validate_run_id(run_id)

    paths = repository_paths(root)
    tracked = {
        PurePosixPath(item.decode("utf-8"))
        for item in git_bytes(root, ["ls-files", "-z"]).split(b"\0")
        if item
    }
    head = git_bytes(root, ["rev-parse", "HEAD"]).decode().strip()

    status = git_bytes(root, ["status", "--porcelain=v1", "-z"])
    staged = git_bytes(
        root,
        ["diff", "--cached", "--binary", "--full-index", "--no-ext-diff"],
    )
    worktree = git_bytes(
        root,
        ["diff", "--binary", "--full-index", "--no-ext-diff"],
    )

    return {
        "format_version": FORMAT_VERSION,
        "runner_version": RUNNER_VERSION,
        "run_id": run_id,
        "created_at": created_at,
        "declared_validation_scope": scope,
        "public_git_head": head,
        "porcelain_status": status.decode("utf-8", errors="surrogateescape"),
        "porcelain_status_sha256": sha256_bytes(status),
        "staged_state_sha256": sha256_bytes(staged),
        "working_tree_diff_sha256": sha256_bytes(worktree),
        "files": [source_entry(root, item, tracked) for item in paths],
    }


def source_checksums_text(manifest: dict[str, Any]) -> bytes:
    """
    Render the source checksum file in lexical path order.
    """

    lines = [
        f"{entry['sha256']}  source/{entry['path']}"
        for entry in manifest["files"]
    ]

    return ("\n".join(lines) + "\n").encode("utf-8")


def add_tar_bytes(
    archive: tarfile.TarFile,
    name: str,
    data: bytes,
    *,
    mode: int = 0o644,
) -> None:
    """
    Add deterministic in-memory bytes to a tar archive.
    """

    info = tarfile.TarInfo(name)
    info.size = len(data)
    info.mode = mode
    info.mtime = 0
    info.uid = 0
    info.gid = 0
    info.uname = "root"
    info.gname = "root"

    archive.addfile(info, fileobj=_BytesReader(data))


class _BytesReader:
    """
    Provide the minimal file interface needed by 'TarFile.addfile'.
    """

    def __init__(self, data: bytes) -> None:
        self.data = data
        self.offset = 0

    def read(self, size: int = -1) -> bytes:
        """
        Read the next byte slice.
        """

        if size < 0:
            size = len(self.data) - self.offset

        chunk = self.data[self.offset : self.offset + size]
        self.offset += len(chunk)

        return chunk


def create_source_archive(
    root: Path,
    archive_path: Path,
    manifest: dict[str, Any],
) -> str:
    """
    Create a deterministic gzip-compressed tar source bundle.
    """

    archive_path.parent.mkdir(parents=True, exist_ok=True)

    with (
        archive_path.open("wb") as raw_handle,
        gzip.GzipFile(
            filename="",
            fileobj=raw_handle,
            mode="wb",
            mtime=0,
        ) as gzip_handle,
        tarfile.open(fileobj=gzip_handle, mode="w") as archive,
    ):
        add_tar_bytes(
            archive,
            "bundle/manifest.json",
            canonical_json(manifest),
        )
        add_tar_bytes(
            archive,
            "bundle/source_checksums.sha256",
            source_checksums_text(manifest),
        )

        for entry in manifest["files"]:
            relative = PurePosixPath(entry["path"])
            source = root.joinpath(*relative.parts)
            name = f"bundle/source/{relative.as_posix()}"

            if entry["type"] == "symlink":
                info = tarfile.TarInfo(name)
                info.type = tarfile.SYMTYPE
                info.linkname = entry["link_target"]
                info.mode = int(entry["mode"], 8)
                info.mtime = 0
                info.uid = 0
                info.gid = 0
                info.uname = "root"
                info.gname = "root"

                archive.addfile(info)
            else:
                add_tar_bytes(
                    archive,
                    name,
                    source.read_bytes(),
                    mode=int(entry["mode"], 8),
                )

    return sha256_file(archive_path)


def validate_archive_member(member: tarfile.TarInfo) -> PurePosixPath:
    """
    Reject archive members that could escape the extraction directory.

    Parameters
    ----------
    member : tarfile.TarInfo
        Archive member whose name, type, and optional target are inspected.

    Returns
    -------
    path : PurePosixPath
        Safe relative member path.

    Raises
    ------
    BundleError
        If the member can escape, has an unsafe link, or uses a denied type.
    """

    path = PurePosixPath(member.name)
    if path.is_absolute() or ".." in path.parts or not path.parts:
        raise BundleError(f"unsafe archive member: {member.name}")

    if member.isdev() or member.isfifo() or member.islnk():
        raise BundleError(f"unsupported archive member: {member.name}")

    if member.issym():
        target = PurePosixPath(member.linkname)
        if target.is_absolute() or ".." in target.parts:
            raise BundleError(f"unsafe archive symlink: {member.name}")

    return path


def safe_extract_archive(archive_path: Path, destination: Path) -> Path:
    """
    Extract a source archive after validating every member and checksum.
    """

    if destination.exists():
        raise BundleError(
            f"extraction destination already exists: {destination}",
        )

    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(
            prefix=f".{destination.name}.",
            dir=destination.parent,
        ),
    )

    try:
        with tarfile.open(archive_path, mode="r:gz") as archive:
            members = archive.getmembers()

            for member in members:
                validate_archive_member(member)

            archive.extractall(temporary, members=members, filter="data")

        bundle = temporary / "bundle"
        manifest = load_json(bundle / "manifest.json")

        verify_source_tree(bundle / "source", manifest)

        bundle.rename(destination)
        temporary.rmdir()
    except Exception:
        shutil.rmtree(temporary, ignore_errors=True)
        raise

    return destination / "source"


def verify_source_tree(source_root: Path, manifest: dict[str, Any]) -> None:
    """
    Verify exact source-tree membership, types, modes, and hashes.

    Parameters
    ----------
    source_root : Path
        Extracted source-tree root.
    manifest : dict[str, Any]
        Expected file, symlink, mode, and digest inventory.

    Raises
    ------
    BundleError
        If membership, entry type, mode, link target, or digest differs.
    """

    if not isinstance(manifest, dict) or not isinstance(
        manifest.get("files"),
        list,
    ):
        raise BundleError("source manifest is malformed")

    expected = {str(item.get("path")): item for item in manifest["files"]}
    if len(expected) != len(manifest["files"]):
        raise BundleError("source manifest contains duplicate paths")

    actual: set[str] = set()

    for path in source_root.rglob("*"):
        if path.is_file() or path.is_symlink():
            actual.add(path.relative_to(source_root).as_posix())

    if actual != set(expected):
        missing = sorted(set(expected) - actual)
        extra = sorted(actual - set(expected))

        raise BundleError(
            f"source membership mismatch: missing={missing}, extra={extra}",
        )

    for relative, entry in expected.items():
        path = source_root / relative

        if entry.get("type") == "symlink":
            if not path.is_symlink():
                raise BundleError(f"expected source symlink: {relative}")

            digest = sha256_bytes(os.readlink(path).encode("utf-8"))
        elif entry.get("type") == "file":
            if not path.is_file() or path.is_symlink():
                raise BundleError(f"expected regular source file: {relative}")

            digest = sha256_file(path)
        else:
            raise BundleError(f"unknown source type: {relative}")

        if digest != entry.get("sha256"):
            raise BundleError(f"source checksum mismatch: {relative}")


def prepared_session(bundle_dir: Path, run_id: str) -> Path:
    """
    Return a validated local session directory.
    """

    validate_run_id(run_id)

    base = bundle_dir.resolve()
    candidate = (base / run_id).resolve()
    if candidate.parent != base:
        raise BundleError(
            "local session path escaped the configured bundle directory",
        )

    return candidate


def session_config(session: Path) -> dict[str, Any]:
    """
    Load and minimally validate a prepared session configuration.
    """

    config = load_json(session / "incoming" / "remote_config.json")
    if config.get("run_id") != session.name:
        raise BundleError(
            "prepared config run ID does not match session directory",
        )

    validate_run_id(str(config.get("run_id", "")))

    return config


def build_remote_config(
    args: argparse.Namespace,
    digest: str,
) -> dict[str, Any]:
    """
    Build the immutable configuration transferred with a source archive.
    """

    run_id = validate_run_id(args.run_id)

    validate_remote_token(args.ssh_host, "SSH host")
    validate_remote_token(args.ssh_user, "SSH user", empty_ok=True)

    validate_remote_token(args.partition, "partition")

    validate_remote_token(args.account, "account")

    validate_remote_token(args.env_name, "environment name")

    validate_remote_python(args.remote_python)

    remote_run_dir = safe_remote_run_dir(args.remote_run_root, run_id)
    resources = {
        "job_count": 2,
        "cpus_per_task": 1,
        "memory": args.memory,
        "time_limit": args.time_limit,
        "poll_seconds": args.poll_seconds,
        "wait_seconds": args.wait_seconds,
    }

    validate_resource_bounds(resources)

    return {
        "format_version": FORMAT_VERSION,
        "runner_version": RUNNER_VERSION,
        "run_id": run_id,
        "declared_validation_scope": args.scope,
        "ssh_host": args.ssh_host,
        "ssh_user": args.ssh_user,
        "remote_run_root": args.remote_run_root,
        "remote_run_dir": remote_run_dir.as_posix(),
        "partition": args.partition,
        "account": args.account,
        "environment_name": args.env_name,
        "remote_python": args.remote_python,
        "result_destination": str(Path(args.result_destination)),
        "resources": resources,
        "source_archive": "source.tar.gz",
        "source_archive_sha256": digest,
        "source_manifest": "source_manifest.json",
    }


def transfer_manifest(incoming: Path) -> dict[str, Any]:
    """
    Build checksums for every prepared transfer file except itself.
    """

    entries = []

    for path in sorted(incoming.rglob("*")):
        if path.is_file() and path.name != "transfer_manifest.json":
            entries.append(
                {
                    "path": path.relative_to(incoming).as_posix(),
                    "sha256": sha256_file(path),
                    "size": path.stat().st_size,
                },
            )

    return {"format_version": FORMAT_VERSION, "files": entries}


def prepare(args: argparse.Namespace) -> int:
    """
    Prepare a source archive and isolated transfer directory.

    Parameters
    ----------
    args : argparse.Namespace
        Validated preparation arguments, including the run ID and bundle root.

    Returns
    -------
    status : int
        Zero after a dry-run inventory or successful bundle preparation.

    Raises
    ------
    BundleError
        If the run ID is unsafe or the isolated session already exists.
    """

    root = Path(args.root).resolve()
    created_at = getattr(args, "created_at", "") or utc_now()
    manifest = build_source_manifest(root, args.run_id, args.scope, created_at)
    session = prepared_session(Path(args.bundle_dir), args.run_id)

    if args.dry_run:
        print(f"run_id={args.run_id}")
        print(f"session={session}")
        print(f"source_files={len(manifest['files'])}")
        print(f"scope={args.scope}")

        for entry in manifest["files"]:
            print(f"include {entry['path']}")

        return 0

    if session.exists():
        raise BundleError(f"session already exists: {session}")

    incoming = session / "incoming"
    incoming.mkdir(parents=True)
    archive = incoming / "source.tar.gz"

    digest = create_source_archive(root, archive, manifest)

    write_json(incoming / "source_manifest.json", manifest)
    (incoming / "source.tar.gz.sha256").write_text(
        f"{digest}  source.tar.gz\n",
        encoding="utf-8",
    )

    launcher_root = incoming / "launcher"
    launcher_script = (
        launcher_root / "tests" / "integration" / "slurm" / "coordinator.py"
    )
    launcher_package = launcher_root / "src" / "protocol_chipseq_signal_norm"
    launcher_script.parent.mkdir(parents=True)
    shutil.copy2(
        root / "tests" / "integration" / "slurm" / "coordinator.py",
        launcher_script,
    )
    shutil.copytree(
        root / "src" / "protocol_chipseq_signal_norm",
        launcher_package,
    )

    config = build_remote_config(args, digest)
    write_json(incoming / "remote_config.json", config)
    write_json(
        incoming / "transfer_manifest.json",
        transfer_manifest(incoming),
    )

    state = {
        "run_id": args.run_id,
        "created_at": created_at,
        "session": str(session),
        "archive_sha256": digest,
        "source_manifest_sha256": sha256_file(
            incoming / "source_manifest.json",
        ),
        "session_created_by": f"slurm_bundle.py/{RUNNER_VERSION}",
    }
    write_json(session / "local_state.json", state)

    print(f"prepared {session}")
    print(f"archive_sha256={digest}")
    print(f"source_files={len(manifest['files'])}")

    return 0


def ssh_target(config: dict[str, Any]) -> str:
    """
    Return the configured SSH target without shell metacharacters.
    """

    host = validate_remote_token(str(config["ssh_host"]), "SSH host")
    user = validate_remote_token(
        str(config.get("ssh_user", "")),
        "SSH user",
        empty_ok=True,
    )

    return f"{user}@{host}" if user else host


def push_commands(session: Path, config: dict[str, Any]) -> list[list[str]]:
    """
    Render isolated remote-directory creation and rsync transfer commands.
    """

    remote = PurePosixPath(str(config["remote_run_dir"])) / "incoming"
    target = ssh_target(config)

    return [
        ["ssh", target, "mkdir", "-p", "--", remote.as_posix()],
        [
            "rsync",
            "-a",
            f"{session / 'incoming'}/",
            f"{target}:{remote.as_posix()}/",
        ],
    ]


def print_commands(commands: Iterable[Sequence[str]]) -> None:
    """
    Print shell-copyable commands.
    """

    for command in commands:
        print(shlex.join(list(command)))


def push(args: argparse.Namespace) -> int:
    """
    Transfer a prepared session into its isolated remote run directory.
    """

    session = prepared_session(Path(args.bundle_dir), args.run_id)
    config = session_config(session)
    commands = push_commands(session, config)

    if args.dry_run:
        print_commands(commands)
        return 0

    for command in commands:
        result = run_command(command)

        if result.returncode:
            raise BundleError(
                f"command failed ({result.returncode}): "
                f"{shlex.join(command)}\n"
                f"{result.stderr.strip()}",
            )

    print(f"pushed {args.run_id} to {ssh_target(config)}")

    return 0


def remote_launch_command(config: dict[str, Any]) -> str:
    """
    Render the explicit three-gate remote launch command.
    """

    incoming = PurePosixPath(str(config["remote_run_dir"])) / "incoming"
    launcher_source = incoming / "launcher" / "src"
    launcher = (
        incoming
        / "launcher"
        / "tests"
        / "integration"
        / "slurm"
        / "coordinator.py"
    )
    config_path = incoming / "remote_config.json"
    archive = incoming / "source.tar.gz"
    command = [
        str(config["remote_python"]),
        launcher.as_posix(),
        "remote-launch",
        "--config",
        config_path.as_posix(),
        "--archive",
        archive.as_posix(),
    ]

    python_path = shlex.quote(launcher_source.as_posix())
    environment = (
        f"PYTHONDONTWRITEBYTECODE=1 PYTHONPATH={python_path} "
        "RUN_SLURM=1 WAIT_SLURM=1 CONFIRM_SLURM_WET=1"
    )

    return f"{environment} {shlex.join(command)}"


def instructions(args: argparse.Namespace) -> int:
    """
    Print exact remote verification, launch, and status commands.
    """

    session = prepared_session(Path(args.bundle_dir), args.run_id)
    config = session_config(session)
    incoming = PurePosixPath(str(config["remote_run_dir"])) / "incoming"
    result = PurePosixPath(str(config["remote_run_dir"])) / "results"

    print(f"# Log in: ssh {shlex.quote(ssh_target(config))}")
    print(
        f"(cd {shlex.quote(incoming.as_posix())} && "
        "sha256sum -c source.tar.gz.sha256)",
    )
    print(remote_launch_command(config))
    print("# Inspect while or after running:")

    launcher_source = incoming / "launcher" / "src"
    status_command = shlex.join(
        [
            str(config["remote_python"]),
            (
                incoming / "launcher/tests/integration/slurm/coordinator.py"
            ).as_posix(),
            "status",
            "--session_dir",
            PurePosixPath(str(config["remote_run_dir"])).as_posix(),
        ],
    )

    print(
        "PYTHONDONTWRITEBYTECODE=1 "
        f"PYTHONPATH={shlex.quote(launcher_source.as_posix())} "
        f"{status_command}",
    )
    print(f"find {shlex.quote(result.as_posix())} -maxdepth 2 -type f -print")

    return 0


def pull_command(
    session: Path,
    config: dict[str, Any],
    destination: Path,
) -> list[str]:
    """
    Render the result-only rsync command.
    """

    remote = PurePosixPath(str(config["remote_run_dir"])) / "results"
    return [
        "rsync",
        "-a",
        f"{ssh_target(config)}:{remote.as_posix()}/",
        f"{destination.resolve()}/",
    ]


def pull(args: argparse.Namespace) -> int:
    """
    Pull the result bundle without modifying remote state.
    """

    session = prepared_session(Path(args.bundle_dir), args.run_id)
    config = session_config(session)
    destination = Path(args.result_destination or config["result_destination"])

    if not destination.is_absolute():
        destination = (repo_root_from_script() / destination).resolve()

    destination = destination / args.run_id
    command = pull_command(session, config, destination)

    if args.dry_run:
        print_commands([command])
        return 0

    if destination.exists():
        raise BundleError(f"result destination already exists: {destination}")

    destination.mkdir(parents=True)

    result = run_command(command)

    if result.returncode:
        raise BundleError(
            (
                f"result pull failed ({result.returncode}): "
                f"{result.stderr.strip()}"
            ),
        )

    print(f"pulled results to {destination}")

    return 0


def parse_checksum_file(path: Path) -> dict[str, str]:
    """
    Parse a strict SHA-256 checksum file.

    Parameters
    ----------
    path : Path
        Checksum file containing two-space-delimited digest rows.

    Returns
    -------
    checksums : dict[str, str]
        Relative paths mapped to lowercase SHA-256 digests.

    Raises
    ------
    BundleError
        If the file cannot be read or a row is malformed or duplicated.
    """

    checksums: dict[str, str] = {}

    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeDecodeError) as exc:
        raise BundleError(f"cannot read checksum file: {path}") from exc

    for number, line in enumerate(lines, 1):
        match = re.fullmatch(r"([0-9a-f]{64})  ([^\0\r\n]+)", line)
        if not match:
            raise BundleError(f"malformed checksum row {number}: {line!r}")

        relative = PurePosixPath(match.group(2))
        if relative.is_absolute() or ".." in relative.parts:
            raise BundleError(f"unsafe checksum path: {relative}")

        name = relative.as_posix()
        if name in checksums:
            raise BundleError(f"duplicate checksum path: {name}")

        checksums[name] = match.group(1)

    return checksums


def verify_result_checksums(result_dir: Path) -> None:
    """
    Verify exact result-file coverage and all declared SHA-256 digests.
    """

    checksums = parse_checksum_file(result_dir / "checksums.sha256")
    actual = {
        path.relative_to(result_dir).as_posix()
        for path in result_dir.rglob("*")
        if path.is_file() and path.name != "checksums.sha256"
    }
    if set(checksums) != actual:
        raise BundleError(
            "result checksum membership mismatch: "
            f"missing={sorted(actual - set(checksums))}, "
            f"extra={sorted(set(checksums) - actual)}",
        )

    for relative, expected in checksums.items():
        if sha256_file(result_dir / relative) != expected:
            raise BundleError(f"result checksum mismatch: {relative}")


def _validated_result_payload(
    result_dir: Path,
    name: str,
    run_id: str,
) -> dict[str, Any]:
    """
    Load one result payload and verify its object and run identity.
    """

    payload = load_json(result_dir / f"{name}.json")

    if not isinstance(payload, dict):
        raise BundleError(f"{name} result must be a JSON object")

    if payload.get("run_id") != run_id:
        raise BundleError(f"{name} run ID mismatch")

    return payload


def _validated_job_inventory(
    jobs_payload: dict[str, Any],
    run_manifest: dict[str, Any],
) -> tuple[dict[Any, dict[str, Any]], list[Any]]:
    """
    Return the complete, uniquely keyed required-job inventory.

    Parameters
    ----------
    jobs_payload : dict[str, Any]
        Result payload containing the observed job records.
    run_manifest : dict[str, Any]
        Prepared manifest containing required job identities.

    Returns
    -------
    by_key, required_ids : tuple[dict[Any, dict[str, Any]], list[Any]]
        Unique job records by identity and the ordered required identities.

    Raises
    ------
    BundleError
        If either inventory is malformed, duplicated, missing, or extraneous.
    """

    jobs = jobs_payload.get("jobs")
    required_ids = run_manifest.get("required_job_keys")
    if not isinstance(jobs, list) or not isinstance(required_ids, list):
        raise BundleError("job or required-job inventory is malformed")

    if len(required_ids) != len(set(required_ids)):
        raise BundleError("required-job inventory contains duplicates")

    by_key = {job.get("job_key"): job for job in jobs if isinstance(job, dict)}
    if len(by_key) != len(jobs) or set(by_key) != set(required_ids):
        raise BundleError(
            "a required job is missing, duplicated, or silently omitted",
        )

    return by_key, required_ids


def _required_jobs_succeeded(
    by_key: dict[Any, dict[str, Any]],
    required_ids: list[Any],
) -> bool:
    """
    Verify every required job and return their combined success.

    Parameters
    ----------
    by_key : dict[Any, dict[str, Any]]
        Complete job records keyed by required job identity.
    required_ids : list[Any]
        Ordered required job identities.

    Returns
    -------
    succeeded : bool
        Whether every required job completed with passing assertions.

    Raises
    ------
    BundleError
        If a job is incomplete, nonterminal, or internally inconsistent.
    """

    computed_success = True

    for key in required_ids:
        job = by_key[key]
        if not set(job) >= REQUIRED_JOB_KEYS:
            raise BundleError(f"job record is incomplete: {key}")

        state = str(job["final_state"]).split("+")[0]
        if state not in TERMINAL_STATES:
            raise BundleError(f"job is not terminal: {key}: {state}")

        assertions = job["assertions"]
        if not isinstance(assertions, list) or not assertions:
            raise BundleError(f"job assertions are missing: {key}")

        passed = all(
            isinstance(item, dict) and item.get("passed") is True
            for item in assertions
        )
        exit_ok = job["exit_code"] == "0:0"
        job_ok = state == "COMPLETED" and exit_ok and passed

        if bool(job.get("success")) != job_ok:
            raise BundleError(
                f"job success field disagrees with evidence: {key}",
            )

        computed_success = computed_success and job_ok

    return computed_success


def _verify_result_summary(
    summary: dict[str, Any],
    job_count: int,
    required_job_count: int,
    computed_success: bool,
) -> None:
    """
    Verify summary counts and success against the job evidence.

    Parameters
    ----------
    summary : dict[str, Any]
        Result summary to compare with verified job evidence.
    job_count : int
        Number of observed job records.
    required_job_count : int
        Number of jobs required by the prepared manifest.
    computed_success : bool
        Success derived from every required job record.

    Raises
    ------
    BundleError
        If any summary count or success value disagrees with job evidence.
    """

    if summary.get("job_count") != job_count:
        raise BundleError("summary job count disagrees with job inventory")

    if summary.get("required_job_count") != required_job_count:
        raise BundleError("summary required-job count disagrees with manifest")

    if bool(summary.get("success")) != computed_success:
        raise BundleError("summary success disagrees with job evidence")


def _read_result_exit_status(result_dir: Path) -> int:
    """
    Read the required integer result exit status.
    """

    try:
        return int((result_dir / "exit_status.txt").read_text().strip())
    except (OSError, ValueError) as exc:
        raise BundleError("exit_status.txt is malformed") from exc


def verify_results(
    result_dir: Path,
    run_id: str,
    source_manifest_sha256: str,
) -> dict[str, Any]:
    """
    Prove result completeness, provenance, checksums, and job agreement.

    Parameters
    ----------
    result_dir : Path
        Pulled result directory for one prepared run.
    run_id : str
        Expected portable run identifier.
    source_manifest_sha256 : str
        Expected prepared source-manifest digest.

    Returns
    -------
    summary : dict[str, Any]
        Verified run identity, job count, success state, and result path.

    Raises
    ------
    BundleError
        If required evidence is absent, malformed, inconsistent, or failed.
    """

    validate_run_id(run_id)

    for relative in REQUIRED_RESULT_PATHS:
        if not (result_dir / relative).exists():
            raise BundleError(f"required result path is missing: {relative}")

    verify_result_checksums(result_dir)

    preflight = _validated_result_payload(result_dir, "preflight", run_id)
    run_manifest = _validated_result_payload(
        result_dir,
        "run_manifest",
        run_id,
    )

    jobs_payload = _validated_result_payload(result_dir, "jobs", run_id)
    summary = _validated_result_payload(result_dir, "summary", run_id)

    if run_manifest.get("source_manifest_sha256") != source_manifest_sha256:
        raise BundleError(
            "result source manifest does not match prepared bundle",
        )

    by_key, required_ids = _validated_job_inventory(
        jobs_payload,
        run_manifest,
    )
    computed_success = _required_jobs_succeeded(by_key, required_ids)
    _verify_result_summary(
        summary,
        len(by_key),
        len(required_ids),
        computed_success,
    )

    exit_status = _read_result_exit_status(result_dir)

    if (exit_status == 0) != computed_success:
        raise BundleError("exit status disagrees with result summary")

    if preflight.get("success") is not True:
        raise BundleError("remote preflight did not succeed")

    if not computed_success:
        raise BundleError("one or more required wet-validation jobs failed")

    return {
        "run_id": run_id,
        "jobs": len(by_key),
        "success": computed_success,
        "result_dir": str(result_dir),
    }


def verify(args: argparse.Namespace) -> int:
    """
    Verify pulled results against the prepared source bundle.
    """

    session = prepared_session(Path(args.bundle_dir), args.run_id)
    state = load_json(session / "local_state.json")
    config = session_config(session)
    destination = Path(args.result_destination or config["result_destination"])

    if not destination.is_absolute():
        destination = (repo_root_from_script() / destination).resolve()

    result = verify_results(
        destination / args.run_id,
        args.run_id,
        str(state["source_manifest_sha256"]),
    )
    print(json.dumps(result, sort_keys=True))

    return 0


def status(args: argparse.Namespace) -> int:
    """
    Inspect a local prepared session or an explicitly selected remote one.

    Parameters
    ----------
    args : argparse.Namespace
        Status arguments selecting a prepared session and local or remote mode.

    Returns
    -------
    status : int
        Zero after printing the selected status document or dry-run command.

    Raises
    ------
    BundleError
        If remote inspection fails or no local status document exists.
    """

    if getattr(args, "session_dir", None):
        session = Path(args.session_dir).resolve()
    else:
        session = prepared_session(Path(args.bundle_dir), args.run_id)

    if args.remote:
        config = session_config(session)
        remote_status = (
            PurePosixPath(str(config["remote_run_dir"]))
            / "results"
            / "live_status.json"
        )
        command = [
            "ssh",
            ssh_target(config),
            "cat",
            "--",
            remote_status.as_posix(),
        ]

        if args.dry_run:
            print_commands([command])
            return 0

        result = run_command(command)

        if result.returncode:
            raise BundleError(f"remote status failed: {result.stderr.strip()}")

        print(result.stdout, end="")

        return 0

    candidates = [
        session / "results" / "live_status.json",
        session / "results" / "summary.json",
        session / "local_state.json",
    ]

    for candidate in candidates:
        if candidate.exists():
            print(candidate.read_text(encoding="utf-8"), end="")

            return 0

    raise BundleError(f"no status is available under {session}")


def clean(args: argparse.Namespace) -> int:
    """
    Remove only a coordinator-created local session directory.
    """

    session = prepared_session(Path(args.bundle_dir), args.run_id)
    marker = session / "local_state.json"
    if not marker.is_file():
        raise BundleError(f"refusing to clean unmarked session: {session}")

    state = load_json(marker)
    if state.get("run_id") != args.run_id or not str(
        state.get("session_created_by", ""),
    ).startswith("slurm_bundle.py/"):
        raise BundleError("refusing to clean a session with an invalid marker")

    if args.dry_run:
        print(f"would remove {session}")

        return 0

    shutil.rmtree(session)
    print(f"removed {session}")

    return 0


def verify_transfer(incoming: Path) -> None:
    """
    Verify exact incoming transfer membership and checksums.

    Parameters
    ----------
    incoming : Path
        Isolated incoming directory containing files and transfer manifest.

    Raises
    ------
    BundleError
        If the manifest, membership, digest form, or file checksum is invalid.
    """

    manifest = load_json(incoming / "transfer_manifest.json")
    entries = manifest.get("files") if isinstance(manifest, dict) else None
    if not isinstance(entries, list):
        raise BundleError("transfer manifest is malformed")

    expected = {item.get("path"): item.get("sha256") for item in entries}
    if len(expected) != len(entries) or not all(
        isinstance(key, str) and isinstance(value, str)
        for key, value in expected.items()
    ):
        raise BundleError("transfer manifest entries are malformed")

    actual = {
        path.relative_to(incoming).as_posix()
        for path in incoming.rglob("*")
        if path.is_file() and path.name != "transfer_manifest.json"
    }
    if actual != set(expected):
        raise BundleError("incoming transfer membership mismatch")

    for relative, digest in expected.items():
        if (
            not SHA256_RE.fullmatch(digest)
            or sha256_file(incoming / relative) != digest
        ):
            raise BundleError(
                f"incoming transfer checksum mismatch: {relative}",
            )


def remote_launch(args: argparse.Namespace) -> int:
    """
    Verify, extract, and launch the bundled wet runner in isolation.

    Parameters
    ----------
    args : argparse.Namespace
        Remote configuration and source-archive paths.

    Returns
    -------
    status : int
        Bundled runner process status.

    Raises
    ------
    BundleError
        If transfer, run identity, archive, configuration, or extraction fails.
    """

    config_path = Path(args.config).resolve()
    archive = Path(args.archive).resolve()
    incoming = config_path.parent

    verify_transfer(incoming)

    config = load_json(config_path)
    run_id = validate_run_id(str(config.get("run_id", "")))
    if archive.name != config.get("source_archive"):
        raise BundleError("remote archive name does not match configuration")

    if sha256_file(archive) != config.get("source_archive_sha256"):
        raise BundleError("remote source archive checksum mismatch")

    session = incoming.parent
    if session.name != run_id:
        raise BundleError("remote session directory does not match run ID")

    source_bundle = session / "source_bundle"
    source_root = safe_extract_archive(archive, source_bundle)
    source_manifest = source_bundle / "manifest.json"
    prepared_manifest = incoming / str(config["source_manifest"])
    if sha256_file(source_manifest) != sha256_file(prepared_manifest):
        raise BundleError(
            "extracted source manifest differs from prepared manifest",
        )

    results = session / "results"
    if results.exists():
        raise BundleError(f"remote result directory already exists: {results}")

    runtime = dict(config)
    runtime.update(
        {
            "session_dir": str(session),
            "source_dir": str(source_root),
            "result_dir": str(results),
            "source_manifest_sha256": sha256_file(source_manifest),
            "runtime_created_at": utc_now(),
        },
    )
    runtime_path = session / "runtime_config.json"
    write_json(runtime_path, runtime)

    runner = (
        source_root / "tests" / "integration" / "slurm" / "run_wet_tests.sh"
    )
    command = ["bash", str(runner), "--config", str(runtime_path)]

    environment = os.environ.copy()
    environment["SLURM_REMOTE_PYTHON"] = validate_remote_python(
        str(config["remote_python"]),
    )
    environment["PROTOCOL_PYTHON"] = environment["SLURM_REMOTE_PYTHON"]
    environment["PY_INVOKE"] = "module"
    bundled_source = str(source_root / "src")
    environment["PYTHONPATH"] = (
        bundled_source + os.pathsep + environment["PYTHONPATH"]
        if environment.get("PYTHONPATH")
        else bundled_source
    )
    environment["PYTHONDONTWRITEBYTECODE"] = "1"

    os.execvpe(command[0], command, environment)
    return 1


def shell_version() -> tuple[str, str]:
    """
    Return the Bash path and first version line.
    """

    bash = shutil.which("bash")
    if not bash:
        return "", ""

    result = run_command([bash, "--version"])

    return bash, result.stdout.splitlines()[0] if result.stdout else ""


def command_probe(name: str) -> dict[str, Any]:
    """
    Record a required remote command's path and version-like output.
    """

    path = shutil.which(name)
    return {"name": name, "path": path, "available": bool(path)}


def validate_resource_bounds(resources: dict[str, Any]) -> None:
    """
    Reject resource requests outside the dedicated tiny wet scope.

    Parameters
    ----------
    resources : dict[str, Any]
        Job count, CPU, memory, wall-time, polling, and wait limits.

    Raises
    ------
    BundleError
        If any resource exceeds the fixed bounded wet-validation scope.
    """

    if resources.get("job_count") != 2 or resources.get("cpus_per_task") != 1:
        raise BundleError("wet scope requires exactly two single-CPU jobs")

    memory = str(resources.get("memory", ""))
    match = re.fullmatch(r"([1-9][0-9]*)([MG])", memory)
    if not match:
        raise BundleError("memory must use a bounded integer M or G value")

    mebibytes = int(match.group(1)) * (1024 if match.group(2) == "G" else 1)
    if mebibytes > 1024:
        raise BundleError("wet memory limit may not exceed 1G")

    limit = str(resources.get("time_limit", ""))
    match_time = re.fullmatch(r"00:0([1-5]):00", limit)
    if not match_time:
        raise BundleError("wet wall time must be between one and five minutes")

    poll = int(resources.get("poll_seconds", 0))
    wait = int(resources.get("wait_seconds", 0))
    if not 1 <= poll <= 30 or not 60 <= wait <= 900:
        raise BundleError(
            "wet polling or wait bound is outside the safe range",
        )


def writable_probe(path: Path) -> bool:
    """
    Prove a directory is writable using a run-local temporary file.
    """

    path.mkdir(parents=True, exist_ok=True)
    probe = path / ".slurm_bundle_write_probe"

    try:
        probe.write_text("probe\n", encoding="utf-8")
        probe.unlink()
    except OSError:
        return False

    return True


def run_preflight(config: dict[str, Any]) -> dict[str, Any]:
    """
    Record and validate the complete remote scheduler preflight.

    Parameters
    ----------
    config : dict[str, Any]
        Verified run configuration for the isolated remote session.

    Returns
    -------
    checks : dict[str, Any]
        Environment, source, scheduler, resource, and directory checks.

    Raises
    ------
    BundleError
        If resources are malformed or outside the bounded wet-test policy.
    """

    resources = config.get("resources")
    if not isinstance(resources, dict):
        raise BundleError("remote resources are malformed")

    validate_resource_bounds(resources)

    session = Path(config["session_dir"])
    results = Path(config["result_dir"])
    source = Path(config["source_dir"])

    bash_path, bash_version = shell_version()
    commands = {
        name: command_probe(name)
        for name in ("sbatch", "squeue", "sacct", "scontrol")
    }

    provider = shutil.which("conda") or shutil.which("mamba")
    env_check = CommandResult((), 1, "", "environment provider unavailable")

    if provider:
        env_check = run_command([provider, "env", "list"])

    partition_check = CommandResult((), 1, "", "scontrol unavailable")

    if commands["scontrol"]["available"]:
        partition_check = run_command(
            ["scontrol", "show", "partition", str(config["partition"]), "-o"],
        )

    source_manifest = source.parent / "manifest.json"
    source_ok = (
        source_manifest.is_file()
        and sha256_file(source_manifest) == config["source_manifest_sha256"]
    )

    archive = session / "incoming" / str(config["source_archive"])
    archive_ok = (
        archive.is_file()
        and sha256_file(archive) == config["source_archive_sha256"]
    )

    try:
        source_payload = load_json(source_manifest)

        verify_source_tree(source, source_payload)
    except BundleError:
        source_ok = False

    python_ok = sys.version_info >= (3, 11)
    bash_match = re.search(r"version ([0-9]+)\.([0-9]+)", bash_version)
    bash_ok = bool(
        bash_match
        and (int(bash_match.group(1)), int(bash_match.group(2))) >= (4, 4),
    )

    env_name = str(config["environment_name"])
    environment_match = re.search(
        rf"(^|\s){re.escape(env_name)}(\s|$)",
        env_check.stdout,
        re.MULTILINE,
    )
    environment_listed = environment_match is not None
    env_ok = env_check.returncode == 0 and environment_listed

    checks = {
        "python_version": python_ok,
        "bash_version": bash_ok,
        "scheduler_commands": all(
            item["available"] for item in commands.values()
        ),
        "environment": env_ok,
        "partition": partition_check.returncode == 0,
        "source_archive_checksum": archive_ok,
        "source_tree_checksum": source_ok,
        "run_directory_writable": writable_probe(session),
        "result_directory_writable": writable_probe(results),
        "partition_selected": bool(config.get("partition")),
        "account_selected": bool(config.get("account")),
        "resource_bounds": True,
    }

    return {
        "run_id": config["run_id"],
        "timestamp": utc_now(),
        "hostname": platform.node(),
        "operating_system": platform.platform(),
        "bash": {"path": bash_path, "version": bash_version},
        "python": {
            "path": sys.executable,
            "version": platform.python_version(),
        },
        "commands": commands,
        "environment": {
            "provider": provider,
            "name": env_name,
            "check_returncode": env_check.returncode,
        },
        "source_archive": {
            "path": str(archive),
            "sha256": config["source_archive_sha256"],
        },
        "directories": {
            "run": str(session),
            "results": str(results),
        },
        "partition": config["partition"],
        "account": config["account"],
        "resources": resources,
        "checks": checks,
        "success": all(checks.values()),
    }


def job_definitions(config: dict[str, Any]) -> list[dict[str, Any]]:
    """
    Define the two bounded fixture-backed submit-worker validations.
    """

    source = Path(config["source_dir"])
    result = Path(config["result_dir"])
    fixtures = result / "artifacts" / "fixtures"

    run_id = str(config["run_id"])
    suffix = hashlib.sha256(run_id.encode()).hexdigest()[:8]
    provider = shutil.which("conda") or shutil.which("mamba") or "conda"

    align_art = result / "artifacts" / "align"
    signal_art = result / "artifacts" / "signal"

    return [
        {
            "job_key": "align_bowtie2_se",
            "job_name": f"wet7i-align-{suffix}",
            "worker": source / "bin" / "submit_align_fastqs.sh",
            "arguments": [
                "--env_nam",
                str(config["environment_name"]),
                "--dir_scr",
                str(source / "bin"),
                "--threads",
                "1",
                "--aligner",
                "bowtie2",
                "--bt2_mode",
                "global",
                "--mapq",
                "0",
                "--index",
                str(fixtures / "bowtie2" / "tiny"),
                "--csv_fil_in",
                str(source / "tests/fixtures/slurm/fastq/tiny_se.fastq"),
                "--dir_out",
                str(align_art),
                "--out_ext",
                "bam",
                "--sfx_se",
                ".fastq",
                "--sfx_pe",
                "_R1.fastq",
                "--dir_eo",
                str(align_art / "logs"),
                "--nam_job",
                f"wet7i-align-worker-{suffix}",
            ],
            "output_directories": [align_art, align_art / "logs"],
            "assertions": [
                {"kind": "nonempty", "path": str(align_art / "tiny_se.bam")},
                {
                    "kind": "nonempty",
                    "path": str(align_art / "tiny_se.bam.bai"),
                },
                {
                    "kind": "command",
                    "command": [
                        provider,
                        "run",
                        "-n",
                        str(config["environment_name"]),
                        "samtools",
                        "quickcheck",
                        str(align_art / "tiny_se.bam"),
                    ],
                    "expected_returncode": 0,
                },
                {
                    "kind": "command_contains",
                    "command": [
                        provider,
                        "run",
                        "-n",
                        str(config["environment_name"]),
                        "samtools",
                        "view",
                        str(align_art / "tiny_se.bam"),
                    ],
                    "text": "tiny_se_read_1\t",
                },
            ],
        },
        {
            "job_key": "compute_signal_bam_se",
            "job_name": f"wet7i-signal-{suffix}",
            "worker": source / "bin" / "submit_compute_signal.sh",
            "arguments": [
                "--env_nam",
                str(config["environment_name"]),
                "--dir_scr",
                str(source / "bin"),
                "--threads",
                "1",
                "--mode",
                "signal",
                "--method",
                "unadj",
                "--csv_fil_in",
                str(fixtures / "tiny_signal.bam"),
                "--chr_sizes",
                str(fixtures / "tiny.fa.fai"),
                "--csv_fil_out",
                str(signal_art / "tiny_se_signal_unadj.bdg"),
                "--csv_usr_frg",
                "NA",
                "--csv_scl_fct",
                "NA",
                "--siz_bin",
                "10",
                "--engine",
                "window",
                "--dp",
                "3",
                "--dir_eo",
                str(signal_art / "logs"),
                "--nam_job",
                f"wet7i-signal-worker-{suffix}",
            ],
            "output_directories": [signal_art, signal_art / "logs"],
            "assertions": [
                {
                    "kind": "nonempty",
                    "path": str(signal_art / "tiny_se_signal_unadj.bdg"),
                },
                {
                    "kind": "contains",
                    "path": str(signal_art / "tiny_se_signal_unadj.bdg"),
                    "text": "I\t0\t10\t10",
                },
                {
                    "kind": "contains",
                    "path": str(signal_art / "tiny_se_signal_unadj.bdg"),
                    "text": "I\t20\t30\t10",
                },
            ],
        },
    ]


def prepare_remote_fixtures(
    config: dict[str, Any],
    commands: list[str],
) -> tuple[bool, str]:
    """
    Materialize tiny derived fixtures inside the isolated result tree.
    """

    source = Path(config["source_dir"])
    result = Path(config["result_dir"])
    fixtures = result / "artifacts" / "fixtures"
    index = fixtures / "bowtie2" / "tiny"

    fixtures.mkdir(parents=True, exist_ok=True)
    index.parent.mkdir(parents=True, exist_ok=True)

    reference = fixtures / "tiny.fa"
    shutil.copy2(
        source / "tests" / "fixtures" / "slurm" / "reference" / "tiny.fa",
        reference,
    )
    provider = shutil.which("conda") or shutil.which("mamba")
    if not provider:
        return False, "Conda or Mamba is unavailable"

    environment = str(config["environment_name"])
    command_list = [
        [
            provider,
            "run",
            "-n",
            environment,
            "bowtie2-build",
            str(reference),
            str(index),
        ],
        [
            provider,
            "run",
            "-n",
            environment,
            "samtools",
            "faidx",
            str(reference),
        ],
        [
            provider,
            "run",
            "-n",
            environment,
            "samtools",
            "view",
            "-b",
            "-o",
            str(fixtures / "tiny_signal.bam"),
            str(source / "tests/fixtures/slurm/sam/tiny_signal.sam"),
        ],
        [
            provider,
            "run",
            "-n",
            environment,
            "samtools",
            "index",
            str(fixtures / "tiny_signal.bam"),
        ],
    ]

    for command in command_list:
        commands.append(shlex.join(command))
        completed = run_command(command)
        if completed.returncode:
            return False, completed.stderr.strip() or completed.stdout.strip()

    expected = [
        reference.with_suffix(".fa.fai"),
        fixtures / "tiny_signal.bam",
        fixtures / "tiny_signal.bam.bai",
        *[
            Path(f"{index}{suffix}")
            for suffix in (
                ".1.bt2",
                ".2.bt2",
                ".3.bt2",
                ".4.bt2",
                ".rev.1.bt2",
                ".rev.2.bt2",
            )
        ],
    ]

    if not all(
        path.is_file() and path.stat().st_size > 0 for path in expected
    ):
        return False, "one or more derived fixture products are missing"

    return True, "tiny committed inputs materialized successfully"


def render_job_script(job: dict[str, Any], destination: Path) -> None:
    """
    Write a small run-local Bash script that execs one submit worker.
    """

    command = ["bash", str(job["worker"]), *job["arguments"]]
    text_value = (
        "#!/usr/bin/env bash\nset -euo pipefail\nexec "
        + shlex.join(command)
        + "\n"
    )
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(text_value, encoding="utf-8")
    destination.chmod(0o700)


def prepare_job_output_directories(
    definitions: Sequence[dict[str, Any]],
) -> tuple[bool, str]:
    """
    Create every worker-configured output directory before submission.
    """

    try:
        for definition in definitions:
            for directory in definition["output_directories"]:
                Path(directory).mkdir(parents=True, exist_ok=True)
    except (KeyError, OSError, TypeError) as exc:
        return False, f"cannot prepare worker output directories: {exc}"

    return True, "all worker output directories prepared"


def sbatch_command(
    config: dict[str, Any],
    job: dict[str, Any],
    script: Path,
) -> list[str]:
    """
    Render one bounded 'sbatch' command.
    """

    resources = config["resources"]
    result = Path(config["result_dir"])

    return [
        "sbatch",
        "--parsable",
        f"--job-name={job['job_name']}",
        "--nodes=1",
        "--ntasks=1",
        "--cpus-per-task=1",
        f"--mem={resources['memory']}",
        f"--time={resources['time_limit']}",
        f"--partition={config['partition']}",
        f"--account={config['account']}",
        f"--output={result / 'stdout' / (job['job_key'] + '.out')}",
        f"--error={result / 'stderr' / (job['job_key'] + '.err')}",
        str(script),
    ]


def parse_sacct(text_value: str, job_id: str) -> dict[str, str] | None:
    """
    Parse the exact top-level job row from pipe-delimited 'sacct' output.
    """

    for line in text_value.splitlines():
        fields = line.strip().split("|")
        if len(fields) < 6 or fields[0] != job_id:
            continue

        return {
            "job_id": fields[0],
            "state": fields[1].split()[0],
            "exit_code": fields[2],
            "submit": fields[3],
            "start": fields[4],
            "end": fields[5],
        }

    return None


def evaluate_assertion(
    spec: dict[str, Any],
    commands: list[str] | None = None,
) -> dict[str, Any]:
    """
    Evaluate a file assertion and retain its evidence.
    """

    path = Path(str(spec.get("path", "")))
    kind = spec.get("kind")

    if kind == "nonempty":
        passed = path.is_file() and path.stat().st_size > 0
        detail = f"size={path.stat().st_size}" if path.is_file() else "missing"
    elif kind == "contains":
        needle = str(spec.get("text", ""))

        try:
            passed = needle in path.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            passed = False

        detail = f"contains={needle!r}"
    elif kind in {"command", "command_contains"}:
        command = spec.get("command")

        if not isinstance(command, list) or not all(
            isinstance(item, str) for item in command
        ):
            passed = False
            detail = "malformed assertion command"
        else:
            if commands is not None:
                commands.append(shlex.join(command))

            completed = run_command(command)
            expected_returncode = int(spec.get("expected_returncode", 0))
            passed = completed.returncode == expected_returncode

            if kind == "command_contains":
                passed = (
                    passed and str(spec.get("text", "")) in completed.stdout
                )

            detail = (
                f"returncode={completed.returncode}; "
                f"stdout_sha256={sha256_bytes(completed.stdout.encode())}; "
                f"stderr_sha256={sha256_bytes(completed.stderr.encode())}"
            )
    else:
        passed = False
        detail = f"unknown assertion kind={kind!r}"

    return {**spec, "passed": passed, "detail": detail}


def write_result_checksums(result: Path) -> None:
    """
    Write sorted checksums for every result file except the list itself.
    """

    lines = []

    for path in sorted(result.rglob("*")):
        if path.is_file() and path.name != "checksums.sha256":
            relative = path.relative_to(result).as_posix()
            lines.append(f"{sha256_file(path)}  {relative}")

    (result / "checksums.sha256").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def finalize_results(
    config: dict[str, Any],
    preflight: dict[str, Any],
    jobs: list[dict[str, Any]],
    commands: list[str],
    required_keys: list[str],
) -> int:
    """
    Write the fixed result schema, summary, exit status, and checksums.
    """

    result = Path(config["result_dir"])
    success = preflight.get("success") is True and len(jobs) == len(
        required_keys,
    )
    success = success and all(job.get("success") is True for job in jobs)

    summary = {
        "run_id": config["run_id"],
        "success": success,
        "job_count": len(jobs),
        "required_job_count": len(required_keys),
        "passed_job_count": sum(job.get("success") is True for job in jobs),
        "failed_job_keys": [
            job["job_key"] for job in jobs if not job.get("success")
        ],
        "missing_job_keys": sorted(
            set(required_keys) - {str(job.get("job_key")) for job in jobs},
        ),
        "finished_at": utc_now(),
    }

    run_manifest = {
        "format_version": FORMAT_VERSION,
        "runner_version": RUNNER_VERSION,
        "run_id": config["run_id"],
        "declared_validation_scope": config["declared_validation_scope"],
        "source_archive_sha256": config["source_archive_sha256"],
        "source_manifest_sha256": config["source_manifest_sha256"],
        "required_job_keys": required_keys,
        "resources": config["resources"],
        "partition": config["partition"],
        "account": config["account"],
        "environment_name": config["environment_name"],
    }

    write_json(result / "preflight.json", preflight)
    write_json(result / "run_manifest.json", run_manifest)
    write_json(
        result / "jobs.json",
        {"run_id": config["run_id"], "jobs": jobs},
    )
    (result / "commands.log").write_text(
        "\n".join(commands) + ("\n" if commands else ""),
        encoding="utf-8",
    )
    write_json(result / "summary.json", summary)
    (result / "summary.txt").write_text(
        f"run_id: {config['run_id']}\n"
        f"success: {'yes' if success else 'no'}\n"
        f"jobs: {len(jobs)}/{len(required_keys)}\n",
        encoding="utf-8",
    )

    status_value = 0 if success else 1
    (result / "exit_status.txt").write_text(
        f"{status_value}\n",
        encoding="utf-8",
    )
    write_json(result / "live_status.json", summary)
    write_result_checksums(result)

    return status_value


def wet_run(args: argparse.Namespace) -> int:
    """
    Run the gated, bounded wet suite and emit its complete result bundle.

    Parameters
    ----------
    args : argparse.Namespace
        Internal wet-run arguments containing the runtime configuration path.

    Returns
    -------
    status : int
        Zero when every required job succeeds; otherwise one.

    Raises
    ------
    BundleError
        If a wet gate is absent or required runtime evidence is malformed.
    """

    for gate in ("RUN_SLURM", "WAIT_SLURM", "CONFIRM_SLURM_WET"):
        if os.environ.get(gate) != "1":
            raise BundleError(
                f"{gate}=1 is required before any wet Slurm action",
            )

    config = load_json(Path(args.config).resolve())
    run_id = validate_run_id(str(config.get("run_id", "")))
    result = Path(config["result_dir"])

    for name in ("stdout", "stderr", "artifacts"):
        (result / name).mkdir(parents=True, exist_ok=False)

    definitions = job_definitions(config)
    required_keys = [job["job_key"] for job in definitions]
    commands: list[str] = []
    jobs: list[dict[str, Any]] = []

    try:
        preflight = run_preflight(config)
    except BundleError as exc:
        preflight = {
            "run_id": run_id,
            "timestamp": utc_now(),
            "success": False,
            "error": str(exc),
        }

    if not preflight.get("success"):
        return finalize_results(
            config,
            preflight,
            jobs,
            commands,
            required_keys,
        )

    fixture_ok, fixture_detail = prepare_remote_fixtures(config, commands)
    preflight["checks"]["fixture_materialization"] = fixture_ok
    preflight["fixture_materialization"] = {
        "success": fixture_ok,
        "detail": fixture_detail,
    }
    preflight["success"] = fixture_ok and all(preflight["checks"].values())
    if not preflight["success"]:
        return finalize_results(
            config,
            preflight,
            jobs,
            commands,
            required_keys,
        )

    output_ok, output_detail = prepare_job_output_directories(definitions)

    preflight["checks"]["worker_output_directories"] = output_ok
    preflight["worker_output_directories"] = {
        "success": output_ok,
        "detail": output_detail,
    }
    preflight["success"] = output_ok and all(preflight["checks"].values())
    if not preflight["success"]:
        return finalize_results(
            config,
            preflight,
            jobs,
            commands,
            required_keys,
        )

    write_json(
        result / "live_status.json",
        {"run_id": run_id, "state": "submitting"},
    )

    for definition in definitions:
        job_key = definition["job_key"]
        artifact = result / "artifacts" / job_key
        script = artifact / "job.sh"
        render_job_script(definition, script)
        command = sbatch_command(config, definition, script)
        test_command = [*command[:1], "--test-only", *command[1:]]

        commands.append(shlex.join(test_command))

        test_result = run_command(test_command)

        if test_result.returncode:
            jobs.append(
                {
                    "job_key": job_key,
                    "job_id": None,
                    "job_name": definition["job_name"],
                    "command": command,
                    "requested_resources": config["resources"],
                    "submit_timestamp": utc_now(),
                    "start_timestamp": None,
                    "finish_timestamp": utc_now(),
                    "final_state": "PREFLIGHT_FAILED",
                    "exit_code": None,
                    "stdout_path": f"stdout/{job_key}.out",
                    "stderr_path": f"stderr/{job_key}.err",
                    "assertions": [
                        {
                            "kind": "sbatch_test_only",
                            "passed": False,
                            "detail": test_result.stderr.strip(),
                        },
                    ],
                    "cleanup_result": "passed",
                    "success": False,
                },
            )

            break

        commands.append(shlex.join(command))
        submitted = run_command(command)
        submit_at = utc_now()
        raw_id = submitted.stdout.strip().split(";", 1)[0]

        if submitted.returncode or not raw_id.isdigit():
            jobs.append(
                {
                    "job_key": job_key,
                    "job_id": None,
                    "job_name": definition["job_name"],
                    "command": command,
                    "requested_resources": config["resources"],
                    "submit_timestamp": submit_at,
                    "start_timestamp": None,
                    "finish_timestamp": utc_now(),
                    "final_state": "SUBMISSION_FAILED",
                    "exit_code": None,
                    "stdout_path": f"stdout/{job_key}.out",
                    "stderr_path": f"stderr/{job_key}.err",
                    "assertions": [
                        {
                            "kind": "sbatch_acceptance",
                            "passed": False,
                            "detail": submitted.stderr.strip(),
                        },
                    ],
                    "cleanup_result": "passed",
                    "success": False,
                },
            )

            break

        record = {
            "job_key": job_key,
            "job_id": raw_id,
            "job_name": definition["job_name"],
            "command": command,
            "requested_resources": config["resources"],
            "submit_timestamp": submit_at,
            "start_timestamp": None,
            "finish_timestamp": None,
            "final_state": "SUBMITTED",
            "exit_code": None,
            "stdout_path": f"stdout/{job_key}.out",
            "stderr_path": f"stderr/{job_key}.err",
            "assertions": [],
            "cleanup_result": "passed",
            "success": False,
        }
        jobs.append(record)

        deadline = time.monotonic() + int(config["resources"]["wait_seconds"])

        while time.monotonic() < deadline:
            queue_command = ["squeue", "-h", "-j", raw_id, "-o", "%T"]
            commands.append(shlex.join(queue_command))
            run_command(queue_command)

            accounting_command = [
                "sacct",
                "-n",
                "-P",
                "-j",
                raw_id,
                "--format=JobIDRaw,State,ExitCode,Submit,Start,End",
            ]
            commands.append(shlex.join(accounting_command))
            accounting = run_command(accounting_command)
            state_row = parse_sacct(accounting.stdout, raw_id)

            if (
                state_row
                and state_row["state"].split("+")[0] in TERMINAL_STATES
            ):
                record["final_state"] = state_row["state"]
                record["exit_code"] = state_row["exit_code"]
                record["start_timestamp"] = state_row["start"] or None
                record["finish_timestamp"] = state_row["end"] or utc_now()

                break

            time.sleep(int(config["resources"]["poll_seconds"]))

        if record["final_state"] == "SUBMITTED":
            record["final_state"] = "UNVERIFIABLE_TIMEOUT"
            record["finish_timestamp"] = utc_now()

        record["assertions"] = [
            evaluate_assertion(spec, commands)
            for spec in definition["assertions"]
        ]
        state = str(record["final_state"]).split("+")[0]
        record["success"] = (
            state == "COMPLETED"
            and record["exit_code"] == "0:0"
            and all(item["passed"] for item in record["assertions"])
        )

        write_json(
            result / "live_status.json",
            {"run_id": run_id, "state": "running", "jobs": jobs},
        )
        if not record["success"]:
            break

    return finalize_results(config, preflight, jobs, commands, required_keys)


def add_common_session_args(parser: CapArgumentParser) -> None:
    """
    Add run-ID and local session-directory arguments.
    """

    parser.add_argument(
        "-ri",
        "--run_id",
        dest="run_id",
        required=True,
        help="Portable prepared run ID.",
    )
    parser.add_argument(
        "-bd",
        "--bundle_dir",
        dest="bundle_dir",
        default=str(DEFAULT_BUNDLE_DIR),
        help=f"Local session base directory (default: {DEFAULT_BUNDLE_DIR}).",
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse coordinator subcommands and their explicit safety options.

    Parameters
    ----------
    argv : list[str] | None
        Explicit arguments, or None to read the process arguments.

    Returns
    -------
    arguments : argparse.Namespace
        Parsed coordinator action, session, and safety options.
    """

    parser = CapArgumentParser(description=__doc__)
    add_help_cap(parser)
    subparsers = parser.add_subparsers(dest="command", required=True)

    prep = subparsers.add_parser("prepare", add_help=False)
    add_help_cap(prep)
    add_common_session_args(prep)

    prep.add_argument(
        "-rt",
        "--root",
        dest="root",
        default=str(repo_root_from_script()),
        help="Public repository root.",
    )
    prep.add_argument(
        "-sc",
        "--scope",
        dest="scope",
        default=DEFAULT_SCOPE,
        help="Declared validation scope.",
    )

    prep.add_argument(
        "-sh",
        "--ssh_host",
        dest="ssh_host",
        required=True,
        help="SSH host alias or name.",
    )
    prep.add_argument(
        "-su",
        "--ssh_user",
        dest="ssh_user",
        default="",
        help="SSH user (default: current SSH configuration).",
    )
    prep.add_argument(
        "-rr",
        "--remote_run_root",
        dest="remote_run_root",
        default=DEFAULT_REMOTE_ROOT.as_posix(),
        help=f"Isolated remote run root (default: {DEFAULT_REMOTE_ROOT}).",
    )

    prep.add_argument(
        "-pt",
        "--partition",
        dest="partition",
        required=True,
        help="Slurm partition.",
    )
    prep.add_argument(
        "-ac",
        "--account",
        dest="account",
        required=True,
        help="Slurm account.",
    )

    prep.add_argument(
        "-en",
        "--env_name",
        dest="env_name",
        default="env_protocol",
        help="Conda environment name (default: env_protocol).",
    )
    prep.add_argument(
        "-rp",
        "--remote_python",
        dest="remote_python",
        required=True,
        help="Absolute remote Python executable path.",
    )

    prep.add_argument(
        "-rd",
        "--result_destination",
        dest="result_destination",
        default=str(DEFAULT_RESULT_DIR),
        help="Local pulled-result base directory.",
    )
    prep.add_argument(
        "-mm",
        "--memory",
        dest="memory",
        default="1G",
        help="Per-job memory, at most 1G (default: 1G).",
    )
    prep.add_argument(
        "-tl",
        "--time_limit",
        dest="time_limit",
        default="00:05:00",
        help="Per-job wall time, at most five minutes.",
    )

    prep.add_argument(
        "-ps",
        "--poll_seconds",
        dest="poll_seconds",
        type=int,
        default=5,
        help="Scheduler poll interval (default: 5).",
    )
    prep.add_argument(
        "-ws",
        "--wait_seconds",
        dest="wait_seconds",
        type=int,
        default=600,
        help="Terminal-state wait bound (default: 600).",
    )
    prep.add_argument(
        "-dr",
        "--dry_run",
        dest="dry_run",
        action="store_true",
        default=False,
        help=(
            "Run script in dry-run mode. Inventory without writing artifacts."
        ),
    )
    prep.set_defaults(handler=prepare)

    for name, handler in (("push", push), ("pull", pull)):
        child = subparsers.add_parser(name, add_help=False)
        add_help_cap(child)
        add_common_session_args(child)
        child.add_argument(
            "-dr",
            "--dry_run",
            dest="dry_run",
            action="store_true",
            default=False,
            help=(
                "Run script in dry-run mode. Print commands without network "
                "access."
            ),
        )

        if name == "pull":
            child.add_argument(
                "-rd",
                "--result_destination",
                dest="result_destination",
                default="",
                help="Override local pulled-result base directory.",
            )

        child.set_defaults(handler=handler)

    inst = subparsers.add_parser("instructions", add_help=False)
    add_help_cap(inst)
    add_common_session_args(inst)
    inst.set_defaults(handler=instructions)

    verify_parser = subparsers.add_parser("verify", add_help=False)
    add_help_cap(verify_parser)
    add_common_session_args(verify_parser)
    verify_parser.add_argument(
        "-rd",
        "--result_destination",
        dest="result_destination",
        default="",
        help="Override local pulled-result base directory.",
    )
    verify_parser.set_defaults(handler=verify)

    status_parser = subparsers.add_parser("status", add_help=False)

    add_help_cap(status_parser)

    status_group = status_parser.add_mutually_exclusive_group(required=True)

    status_group.add_argument(
        "-ri",
        "--run_id",
        dest="run_id",
        default=None,
        help="Prepared run ID.",
    )
    status_group.add_argument(
        "-sd",
        "--session_dir",
        dest="session_dir",
        default=None,
        help="Explicit local or remote session directory.",
    )
    status_parser.add_argument(
        "-bd",
        "--bundle_dir",
        dest="bundle_dir",
        default=str(DEFAULT_BUNDLE_DIR),
        help="Local session base directory.",
    )
    status_parser.add_argument(
        "-rm",
        "--remote",
        dest="remote",
        action="store_true",
        default=False,
        help="Inspect through SSH; never the default.",
    )
    status_parser.add_argument(
        "-dr",
        "--dry_run",
        dest="dry_run",
        action="store_true",
        default=False,
        help=(
            "Run script in dry-run mode. Print the remote status command only."
        ),
    )
    status_parser.set_defaults(handler=status)

    clean_parser = subparsers.add_parser("clean", add_help=False)
    add_help_cap(clean_parser)
    add_common_session_args(clean_parser)
    clean_parser.add_argument(
        "-dr",
        "--dry_run",
        dest="dry_run",
        action="store_true",
        default=False,
        help=(
            "Run script in dry-run mode. Show confined removal without "
            "deleting."
        ),
    )
    clean_parser.set_defaults(handler=clean)

    remote_parser = subparsers.add_parser("remote-launch", add_help=False)
    add_help_cap(remote_parser)
    remote_parser.add_argument(
        "-cf",
        "--config",
        dest="config",
        required=True,
        help="Transferred remote configuration file.",
    )
    remote_parser.add_argument(
        "-ar",
        "--archive",
        dest="archive",
        required=True,
        help="Transferred source archive.",
    )
    remote_parser.set_defaults(handler=remote_launch)

    wet_parser = subparsers.add_parser("wet-run", add_help=False)
    add_help_cap(wet_parser)
    wet_parser.add_argument(
        "-cf",
        "--config",
        dest="config",
        required=True,
        help="Extracted runtime configuration file.",
    )
    wet_parser.set_defaults(handler=wet_run)

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """
    Dispatch one coordinator command and print fail-closed diagnostics.

    Parameters
    ----------
    argv : list[str] | None
        Explicit arguments, or None to read the process arguments.

    Returns
    -------
    status : int
        Stable process status from the selected coordinator action.
    """

    try:
        args = parse_args(argv)
        status_value = int(args.handler(args))

        if status_value == 0:
            return 0

        return status_value
    except BundleError as exc:
        print(f"error(slurm_bundle.py): {exc}", file=sys.stderr)

        return 1


if __name__ == "__main__":
    raise SystemExit(main())
