
# Slurm dry-run test fixtures
These fixtures let a driver's Slurm-mode dry run be tested on a machine without Slurm. A Slurm dry run prints the `sbatch` call it would submit, but the driver first checks that `sbatch` is on `PATH`; the stub only lets it pass that check.

The stub fails if it is ever called, so a test using it also proves that the dry run did not submit anything. Tests that Slurm accepts a submission belong to the real-Slurm suite under `tests/integration/slurm/`, whose fixtures are in `tests/fixtures/slurm/`.

Regenerate fixtures from the repository root with:
```bash
bash tests/fixtures/slurm_dry_run/make.sh
```

Generated fixture outputs are ignored by Git. `tests/run_tests.sh` regenerates this fixture set automatically when required inputs are missing.

<br />

## Files
Readable provenance:
- `make.sh`

Generated command stub:
- `tool/sbatch`: prints an error and exits 1 when called

<br />

## Current and deferred test coverage
Current coverage in `tests/contract/interfaces/test_download_fastqs_interface.sh`:
- a Slurm dry run of `execute_download_fastqs.sh` with `--threads 2` succeeds with the stub first on `PATH`, prints one `sbatch` call, and wraps `parallel --jobs 2` around the run's GNU Parallel config.

Deferred:
- Slurm dry-run rows for the other `execute_*.sh` drivers.
