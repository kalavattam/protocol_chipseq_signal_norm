
# Semantic-movement test fixtures
These fixtures are synthetic micro-fixtures for fast, deterministic tests of the semantic-movement record validator in `dev/audit/semantic_movement.py`, which decides whether an obligation moved between rule owners completely and under human authority.

They are intentionally small and hand-checkable. Running `make.sh` regenerates the fixture set deterministically.

Regenerate fixtures from the repository root with:
```bash
bash tests/fixtures/semantic_movement/make.sh
```

Generated fixture outputs are ignored by Git. `tests/run_tests.sh` regenerates this fixture set automatically when required inputs are missing.

This fixture set is one record of every verdict at once rather than a set of inputs each carrying one, so it sits at the fixture root with no verdict directories beneath it. The unit test derives its rejected, boundary, and non-applicable cases by perturbing the accepted record, which keeps each rejection one edit away from an accepted record instead of a separately drifting document.

<br />

## Files
Readable provenance:
- `make.sh`

Generated verdict record:
- `cases.json`: one complete accepted movement record, in the canonical form owned by `JSON.SOURCE.FORM`

<br />

## Expected movement behavior
The accepted record validates with no errors. It carries both directions of the completeness claim — old-to-new dispositions and new-to-old provenance — along with source and destination fingerprints, so a movement cannot be recorded as complete on one side alone.

Perturbing the record must fail in specific ways rather than generally: a consequential delta without explicit approval is rejected for that reason, an LLM named as the authorizing role is rejected for that reason, and adding or removing a field is rejected by the schema before either check runs.

<br />

## Current and deferred test coverage
Current coverage in `tests/unit/dev_audit/test_semantic_movement.py`:
- the accepted record validates cleanly;
- a consequential delta lacking approval is rejected;
- an LLM cannot be the authorizing role, only a supporting comparison; and
- the supplied schema rejects both missing and extra fields.

Deferred:
- multi-record reconciliation across a whole roadmap, which is human-owned and has no checker to exercise.
