#!/usr/bin/env python3
"""Validate read-only numeric-output applicability and migration readiness."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from jsonschema import Draft202012Validator

RULE_SCHEMA = "OUTPUT.NUMERIC.APPLICABILITY.SCHEMA"
RULE_REFERENCE = "OUTPUT.NUMERIC.APPLICABILITY.REFERENCE"
RULE_READINESS = "OUTPUT.NUMERIC.MIGRATION.READINESS"
TERMINAL = {
    "applicable_dp_governed_migrated",
    "explicitly_excepted",
    "proved_not_dp_governed",
}


def _sha(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _schema_findings(
    data: dict[str, Any], schema: dict[str, Any]
) -> list[dict[str, str]]:
    findings: list[dict[str, str]] = []
    for error in Draft202012Validator(schema).iter_errors(data):
        locator = "/".join(str(x) for x in error.absolute_path) or "/"
        findings.append(
            {"rule_id": RULE_SCHEMA, "message": f"{locator}: {error.message}"}
        )
    return sorted(
        findings, key=lambda item: (item["rule_id"], item["message"])
    )


def validate(
    config: dict[str, Any], inventory: dict[str, Any], schema: dict[str, Any]
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    """Validate immutable inventory references without remapping a site."""

    findings = _schema_findings(config, schema)
    if findings:
        return findings, {}
    sites = {item["site_id"]: item for item in inventory.get("sites", [])}
    records = config["records"]
    seen: set[str] = set()
    unresolved = 0
    stale = 0
    for record in records:
        site_id = record["site_id"]
        if site_id in seen:
            findings.append(
                {
                    "rule_id": RULE_REFERENCE,
                    "message": f"duplicate stable surface ID: {site_id}",
                }
            )
            continue
        seen.add(site_id)
        site = sites.get(site_id)
        if site is None:
            stale += 1
            findings.append(
                {
                    "rule_id": RULE_REFERENCE,
                    "message": f"missing inventory surface: {site_id}",
                }
            )
            continue
        for key in ("path", "symbol", "fingerprint"):
            if record[key] != site[key]:
                stale += 1
                findings.append(
                    {
                        "rule_id": RULE_REFERENCE,
                        "message": f"stale {key} for surface: {site_id}",
                    }
                )
        if record[
            "classification"
        ] == "explicitly_excepted" and not record.get("exception_owner"):
            findings.append(
                {
                    "rule_id": RULE_SCHEMA,
                    "message": f"exception lacks owner: {site_id}",
                }
            )
        if record["classification"] == "unresolved":
            unresolved += 1

    candidate_ids = {
        item["site_id"]
        for item in inventory.get("sites", [])
        if item["language"] == "python"
        and item["numeric_relevance"] == "candidate"
    }
    unregistered = sorted(candidate_ids - seen)
    unresolved += len(unregistered)
    cohort = config["selected_cohort"]
    cohort_findings = []
    for site_id in cohort["site_ids"]:
        record = next(
            (item for item in records if item["site_id"] == site_id), None
        )
        if record is None or record["classification"] not in TERMINAL:
            cohort_findings.append(site_id)
    if cohort["declared_ready"] and cohort_findings:
        findings.append(
            {
                "rule_id": RULE_READINESS,
                "message": "ready cohort contains no terminal disposition: "
                + ", ".join(cohort_findings),
            }
        )

    artifact = {
        "schema_version": 1,
        "inventory_sha256": _sha(inventory),
        "config_sha256": _sha(config),
        "registered_record_count": len(records),
        "unresolved_count": unresolved,
        "unregistered_candidate_count": len(unregistered),
        "stale_reference_result": (
            "not_exercised_no_registered_records"
            if not records
            else "pass"
            if stale == 0
            else "fail"
        ),
        "selected_cohort_readiness": "ready"
        if cohort["declared_ready"] and not cohort_findings
        else "not_ready",
        "global_s3_mig_002_readiness": "blocked" if unresolved else "ready",
        "global_pre_0_3_0_readiness": "blocked" if unresolved else "ready",
        "records": records,
    }
    return sorted(
        findings, key=lambda item: (item["rule_id"], item["message"])
    ), artifact


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("dev/config/numeric_output_applicability.json"),
    )
    parser.add_argument(
        "--schema",
        type=Path,
        default=Path("dev/schemas/numeric_output_applicability.schema.json"),
    )
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args(argv)
    root = args.root.resolve()

    def resolve(path: Path) -> Path:
        return path if path.is_absolute() else root / path

    findings, artifact = validate(
        json.loads(resolve(args.config).read_text(encoding="utf-8")),
        json.loads(resolve(args.inventory).read_text(encoding="utf-8")),
        json.loads(resolve(args.schema).read_text(encoding="utf-8")),
    )
    if args.json_out is not None:
        target = resolve(args.json_out)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(
            json.dumps(artifact, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    for finding in findings:
        print(f"{finding['rule_id']}: {finding['message']}")
    return 1 if findings else 0


if __name__ == "__main__":
    raise SystemExit(main())
