
# Submit-bootstrap test fixtures
These fixtures are for `tests/contract/interfaces/test_submit_bootstrap.sh`, which checks how the `submit_*.sh` drivers find and load their helper libraries. They need no managed environment.

Regenerate fixtures from the repository root with:
```bash
bash tests/fixtures/submit_bootstrap/make.sh
```

Generated fixture outputs are ignored by Git. `tests/run_tests.sh` regenerates this fixture set automatically when required inputs are missing.

<br />

## Files
Readable provenance:
- `make.sh`

A sourced helper:
- `script/counted_helper.sh`: adds one to `COUNTED_HELPER_LOADS` each time it is sourced

<br />

## Current and deferred test coverage
Current coverage:
- `source_once` sources a helper once, however often it is called.

Deferred:
- none.
