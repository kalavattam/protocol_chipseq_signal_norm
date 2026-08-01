from __future__ import annotations

import copy
import json
from pathlib import Path

from dev.audit.numeric_output_applicability import validate
from dev.audit.run_numeric_output_applicability import main as adapter_main


def inventory() -> dict:
    return {
        "sites": [
            {
                "site_id": "site-1",
                "path": "src/example.py",
                "symbol": "emit",
                "fingerprint": "a" * 64,
                "language": "python",
                "numeric_relevance": "candidate",
            }
        ]
    }


def config() -> dict:
    return {
        "schema_version": 1,
        "records": [
            {
                "site_id": "site-1",
                "path": "src/example.py",
                "symbol": "emit",
                "emission_role": "user_output",
                "fingerprint": "a" * 64,
                "classification": "unresolved",
                "dp_provenance": "inherited args.dp",
                "owner": "OUTPUT.NUMERIC.RENDERING",
                "rationale": "awaiting migration",
                "review_state": "evidence_only",
                "migration_disposition": "deferred",
            }
        ],
        "selected_cohort": {
            "id": "none",
            "declared_ready": False,
            "site_ids": [],
        },
    }


def schema() -> dict:
    from pathlib import Path

    return json.loads(
        (
            Path(__file__).resolve().parents[3]
            / "dev/schemas/numeric_output_applicability.schema.json"
        ).read_text()
    )


def test_unresolved_evidence_passes_but_blocks_global_readiness() -> None:
    findings, artifact = validate(config(), inventory(), schema())
    assert findings == []
    assert artifact["global_s3_mig_002_readiness"] == "blocked"


def test_ready_cohort_rejects_unresolved_member() -> None:
    data = config()
    data["selected_cohort"] = {
        "id": "pilot",
        "declared_ready": True,
        "site_ids": ["site-1"],
    }
    findings, _ = validate(data, inventory(), schema())
    assert [item["rule_id"] for item in findings] == [
        "OUTPUT.NUMERIC.MIGRATION.READINESS"
    ]


def test_stale_reference_rejects_without_remap() -> None:
    data = copy.deepcopy(config())
    data["records"][0]["fingerprint"] = "b" * 64
    findings, _ = validate(data, inventory(), schema())
    assert any(
        item["rule_id"] == "OUTPUT.NUMERIC.APPLICABILITY.REFERENCE"
        for item in findings
    )


def test_schema_and_each_owned_fault_diagnostic_are_reachable() -> None:
    data = copy.deepcopy(config())
    del data["records"][0]["owner"]
    findings, _ = validate(data, inventory(), schema())
    assert {item["rule_id"] for item in findings} == {
        "OUTPUT.NUMERIC.APPLICABILITY.SCHEMA"
    }

    data = copy.deepcopy(config())
    data["records"][0]["fingerprint"] = "b" * 64
    findings, _ = validate(data, inventory(), schema())
    assert {item["rule_id"] for item in findings} == {
        "OUTPUT.NUMERIC.APPLICABILITY.REFERENCE"
    }

    data = copy.deepcopy(config())
    data["selected_cohort"] = {
        "id": "fault",
        "declared_ready": True,
        "site_ids": ["site-1"],
    }
    findings, _ = validate(data, inventory(), schema())
    assert {item["rule_id"] for item in findings} == {
        "OUTPUT.NUMERIC.MIGRATION.READINESS"
    }


def test_registered_adapter_emits_each_owned_diagnostic(
    tmp_path: Path, capsys
) -> None:
    inventory_path = tmp_path / "inventory.json"
    inventory_path.write_text(json.dumps(inventory()))
    schema_path = tmp_path / "schema.json"
    schema_path.write_text(json.dumps(schema()))
    for expected, mutate in (
        (
            "OUTPUT.NUMERIC.APPLICABILITY.SCHEMA",
            lambda value: value["records"][0].pop("owner"),
        ),
        (
            "OUTPUT.NUMERIC.APPLICABILITY.REFERENCE",
            lambda value: value["records"][0].update(
                {"fingerprint": "b" * 64}
            ),
        ),
        (
            "OUTPUT.NUMERIC.MIGRATION.READINESS",
            lambda value: value.update(
                {
                    "selected_cohort": {
                        "id": "fault",
                        "declared_ready": True,
                        "site_ids": ["site-1"],
                    }
                }
            ),
        ),
    ):
        payload = copy.deepcopy(config())
        mutate(payload)
        config_path = tmp_path / f"{expected}.json"
        config_path.write_text(json.dumps(payload))
        assert (
            adapter_main(
                [
                    "--config",
                    str(config_path),
                    "--schema",
                    str(schema_path),
                    "--inventory",
                    str(inventory_path),
                ]
            )
            == 1
        )
        assert expected in capsys.readouterr().out
