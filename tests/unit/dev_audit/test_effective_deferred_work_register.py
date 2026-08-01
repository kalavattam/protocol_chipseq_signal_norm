"""Negative coverage for the Phase C effective-register validator."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
FOLLOW_UP = (
    ROOT
    / "artifacts/reviews/2026-07-30_post_pilot_standards_phase_c/follow_up"
)


def load_validator():
    spec = importlib.util.spec_from_file_location(
        "phase_c_effective_register_validator",
        FOLLOW_UP / "validate_effective_deferred_work_register.py",
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize(
    "fault",
    [
        "hash",
        "unauthorized",
        "missing_new",
        "parity",
        "mig4_missing",
        "mig4_extra",
        "duplicate_amendment",
        "orphan_amendment",
        "acceptance_state",
        "semantic_old",
        "missing_semantic_amendment",
        "stale_blocker",
        "dormant_blocking",
        "callable_cross_reference",
        "overview_omission",
        "overview_duplication",
        "overview_label_drift",
        "overview_markdown",
        "overview_timing",
    ],
)
def test_effective_register_rejects_negative_cases(
    tmp_path: Path, fault: str
) -> None:
    module = load_validator()
    data = json.loads(
        (FOLLOW_UP / "effective_deferred_work_register.json").read_text()
    )
    markdown = (FOLLOW_UP / "effective_deferred_work_register.md").read_text()
    contracts = json.loads(
        (ROOT / "dev/config/help_contracts.json").read_text()
    )
    if fault == "hash":
        data["base"]["json_sha256"] = "0" * 64
    elif fault == "unauthorized":
        data["records"][0]["decision"] = "unauthorized change"
    elif fault == "missing_new":
        data["semantic_reaudit_amendments"][0]["fields"]["current_status"].pop(
            "new"
        )
    elif fault == "mig4_missing":
        data["effective_amendments"][-1]["fields"].pop("decision")
    elif fault == "mig4_extra":
        data["effective_amendments"][-1]["fields"]["extra"] = {
            "old": None,
            "new": "extra",
        }
    elif fault == "duplicate_amendment":
        data["effective_amendments"].append(data["effective_amendments"][0])
    elif fault == "orphan_amendment":
        data["semantic_reaudit_amendments"].append(
            {"id": "S3-ORPHAN", "provenance": "fault", "fields": {}},
        )
    elif fault == "acceptance_state":
        standard = next(
            item for item in data["records"] if item["id"] == "S3-STD-001"
        )
        standard["current_status"] = "phase_c_foundation_corrected_accepted"
    elif fault == "semantic_old":
        data["semantic_reaudit_amendments"][0]["fields"]["next_action"][
            "old"
        ] = "wrong old value"
    elif fault == "missing_semantic_amendment":
        data["semantic_reaudit_amendments"].pop()
    elif fault == "stale_blocker":
        standard = next(
            item for item in data["records"] if item["id"] == "S3-STD-005"
        )
        standard["blocking_condition"] = "The blockquote owner is absent."
    elif fault == "dormant_blocking":
        dormant = next(
            item for item in data["records"] if item["id"] == "S3-DORMANT-001"
        )
        dormant["blocking_condition"] = "Dormant languages block Phase C."
    elif fault == "callable_cross_reference":
        callable_record = next(
            item
            for item in contracts["examples"]
            if item["surface_id"] == "compute_input_floor_callable"
        )
        callable_record["deferred_record"] = "S3-MIG-003"
    elif fault == "overview_omission":
        data["overview"]["primary_categories"]["roadmaps_0_3_0"].pop()
    elif fault == "overview_duplication":
        data["overview"]["primary_categories"]["dormant_languages"].append(
            {"id": "S3-GATE-001", "label": "duplicate"},
        )
    elif fault == "overview_label_drift":
        data["overview"]["primary_categories"]["completed_transitions"][0][
            "label"
        ] = "drift"
    elif fault == "overview_markdown":
        markdown = markdown.replace(
            "### Dormant languages",
            "### Drifted languages",
            1,
        )
    elif fault == "overview_timing":
        data["overview"]["cross_cutting_residual_notes"][1]["note"] = (
            "Its main roadmap is post-0.3.0."
        )
    else:
        markdown = markdown.replace(
            '"phase_c_foundation_implemented_accepted"',
            '"parity drift"',
            1,
        )
    (tmp_path / "effective_deferred_work_register.json").write_text(
        json.dumps(data)
    )
    (tmp_path / "effective_deferred_work_register.md").write_text(markdown)
    contracts_path = tmp_path / "help_contracts.json"
    contracts_path.write_text(json.dumps(contracts))
    module.HERE = tmp_path
    module.HELP_CONTRACTS = contracts_path
    with pytest.raises((AssertionError, KeyError)):
        module.main()
