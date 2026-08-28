
# Test suite
Test code is classified by what it proves, not by implementation language.

```text
tests/
├── unit/                 isolated Python behavior
├── contract/             repository and command-interface policy
├── integration/
│   ├── local/            fixture-backed local workflows
│   ├── parallel/         gated GNU Parallel workflows
│   └── slurm/            gated scheduler integration and wet coordinator
├── fixtures/             tracked static inputs/recipes/docs plus generated data
└── support/              shared test-only helpers and cleanup
```

Production entrypoints live in `bin/`, sourced Bash in `lib/bash/`, and Python package code in `src/`. Executable validation coordinators, assertions, fakes, and test harnesses belong under `tests/`.

<br />

## Canonical runner
Run the safe suite from the repository root:
```bash
bash tests/run_tests.sh
```

The default `all-safe` selection runs unit, contract, local integration, and parallel integration groups. Parallel scripts self-skip unless `RUN_PARALLEL=1`. Select groups explicitly when narrowing validation:
```bash
bash tests/run_tests.sh unit contract
bash tests/run_tests.sh integration-local
RUN_PARALLEL=1 bash tests/run_tests.sh integration-parallel
```

Inspect the exact, deduplicated discovery set without running it:
```bash
bash tests/run_tests.sh --list all-safe
```

Logs and temporary products default to `artifacts/tests/`. Set `TEST_ARTIFACT_ROOT` to an absolute path outside the repository for non-mutating validation. In-repository overrides are rejected. The runner discovers every selected shell test at most once and never discovers `tests/integration/slurm/run_wet_tests.sh`.

<br />

## Fixtures and cleanup
Fixture recipes are `tests/fixtures/<feature>/make.sh`. Generated fixture data is ignored; each tracked `README.md` describes provenance and expectations. `make.sh` generates and nothing else — tests and checkers do the validating.

```bash
bash tests/support/clean_artifacts.sh --dry_run --all
bash tests/support/clean_artifacts.sh --fixtures
bash tests/support/clean_artifacts.sh --outputs
```

<br />

## Gates
`RUN_ATRIA=1`, `RUN_DOWNLOAD=1`, `RUN_INSTALL_ENVS=1`, and `RUN_PARALLEL=1` enable their bounded optional integrations. `RUN_INSTALL_ENVS=1` creates and updates a disposable real `env_siqchip` under the test artifact root. Scheduler submission is separate: the ordinary Slurm integration requires `RUN_SLURM=1`, while the checksummed two-job wet validation requires exact `RUN_SLURM=1 WAIT_SLURM=1 CONFIRM_SLURM_WET=1` gates and the workflow in [`integration/slurm/README.md`](integration/slurm/README.md). The wet runner interprets `WAIT_SLURM` as an exact confirmation; ordinary Boolean normalization does not apply to it.

Boolean test gates accept true, t, yes, y, 1 and false, f, no, n, 0 case-insensitively. Unset or empty gates are disabled; surrounding whitespace and other nonempty values are invalid.
