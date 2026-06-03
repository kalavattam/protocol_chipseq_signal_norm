
# Test suite
The repository uses Bash smoke tests to exercise wrapper startup, help output, fixture-backed local workflows, and selected optional execution paths.

<br />

## Run smoke tests
Run the default suite from the repository root:
```bash
bash tests/scripts/run_smoke_tests.sh
```

Logs and temporary products are written under `tests/outputs/`. The default suite exercises local serial workflow paths and default-safe parser, syntax, help, and installation-layout checks.

<br />

## Generated fixtures
Fixture recipes live under workflow-specific directories:
```text
tests/align_fastqs/scripts/make_fixtures.sh
tests/calculate_scaling_factor/scripts/make_fixtures.sh
tests/compute_signal/scripts/make_fixtures.sh
tests/download_fastqs/scripts/make_fixtures.sh
tests/filter_alignments/scripts/make_fixtures.sh
tests/trim_fastqs/scripts/make_fixtures.sh
```

`tests/scripts/run_smoke_tests.sh` detects missing required fixtures and regenerates them automatically. By default, it prepares `align_fastqs`, `calculate_scaling_factor`, `compute_signal`, `download_fastqs`, and `filter_alignments` fixtures. It prepares `trim_fastqs` fixtures only when `RUN_ATRIA=1` is set (see below).

Generated fixture outputs are ignored by Git and are not to be committed. Each `tests/*/fixtures/README.md` file is tracked documentation for its fixture set.

<br />

## Clean generated artifacts
Preview cleanup without deleting files:
```bash
bash tests/scripts/clean_test_outputs.sh --dry_run --all
```

Remove ignored generated fixtures, smoke-test outputs, or both:
```bash
bash tests/scripts/clean_test_outputs.sh --fixtures
bash tests/scripts/clean_test_outputs.sh --outputs
bash tests/scripts/clean_test_outputs.sh --all
```

The cleaner uses scoped `git clean -dX` commands and will not remove tracked fixture files.

<br />

## Optional gates
Optional dependency classes are enabled explicitly:
```bash
RUN_ATRIA=1 bash tests/scripts/run_smoke_tests.sh
RUN_PARALLEL=1 bash tests/scripts/run_smoke_tests.sh
RUN_SLURM=1 bash tests/scripts/run_smoke_tests.sh
RUN_SLURM=1 SLURM_WAIT=1 bash tests/scripts/run_smoke_tests.sh
```
| gate             | function                                                                              | todo                                                                                                |
| :----            | :----                                                                                 | :----                                                                                               |
| `RUN_ATRIA=1`    | Enables Atria-backed `trim_fastqs` tests and fixture bootstrap.                       | -                                                                                                   |
| `RUN_PARALLEL=1` | Enables local GNU Parallel tests.                                                     | -                                                                                                   |
| `RUN_SLURM=1`    | Enables the Slurm-backed `align_fastqs` smoke test. Run it on a Slurm-capable system. | Add Slurm coverage for `trim_fastqs`, `download_fastqs`, `filter_alignments`, and `compute_signal`. |
| `SLURM_WAIT=1`   | Makes the Slurm test poll for expected outputs after submission.                      | Rename to verb-first `WAIT_SLURM=1`.                                                                |

Tests that require multiple dependency classes require all relevant gates, such as `RUN_ATRIA=1 RUN_PARALLEL=1` for parallel trimming.

<br />

## Tools and coverage
Use Bash >= 4.4. The default fixture-backed suite expects `env_protocol` to be active or available through Conda. Depending on the workflow, it uses tools such as Python, Samtools, Bowtie2, BWA, bwa-mem2, `wget`, and `gzip`. Optional suites additionally require Atria with `pigz` and `pbzip2`, GNU Parallel, or a Slurm installation with `sbatch`.

Smoke groups cover:
- shell syntax, Python startup, help output/style, installation layout, and cleanup dry-runs;
- direct spike-in and siQ-ChIP scaling-factor part-file combination and serial
  execute-layer spike-in assembly;
- local `download_fastqs`, `trim_fastqs`, `align_fastqs`, `filter_alignments`, and `compute_signal` wrapper paths;
- selected GNU Parallel paths and one remote Slurm submission path.

Fixture details and expected products are documented in each `tests/*/fixtures/README.md`. Broader `calculate_scaling_factor` coverage for siQ-ChIP, GNU Parallel, and Slurm paths is still pending.
