
# Portable Slurm wet validation
This workflow captures the current public working tree on the Mac, transfers an exact source archive into an isolated FHCC run directory, executes a two-job real-Slurm validation without Codex, and returns a checksummed evidence tree for local verification. It never depends on the remote repository checkout being current or clean.

The coordinator requires Python >= 3.11:
```text
tests/integration/slurm/coordinator.py
```

The scheduler entrypoint is deliberately separate from the ordinary test suite:
```text
tests/integration/slurm/run_wet_tests.sh
```

`tests/run_tests.sh` does not invoke the wet runner.

<br />

## Validation scope
The bounded wet suite submits exactly two single-CPU jobs using dedicated text-only fixtures under `tests/fixtures/slurm/`:
1. `submit_align_fastqs.sh` aligns one 30-base FASTQ read with Bowtie2 and asserts that the BAM and BAI are nonempty.
2. `submit_compute_signal.sh` computes an unadjusted 10-base-bin signal track from two tiny SAM-derived alignments and asserts the expected chromosome-I rows.

Before submission, the runner copies the committed reference into `results/artifacts/fixtures/`, builds its tiny Bowtie2 index, and converts the committed SAM to an indexed BAM in that same run-specific directory. It does not use production data, download data, or write generated products into the extracted source tree.

The fixed scheduler bounds are two jobs, one node, one task, one CPU per task, at most 1 GiB memory, and at most five minutes wall time per job. Preparation rejects wider values.

<br />

## Local preparation
Choose a unique run ID and supply the FHCC values for your account. The remote root may be absolute or relative to the SSH login directory; it must not contain `..` and cannot be `/`.

```bash
run_id=7i-slurm-bundle-20260716

python3 tests/integration/slurm/coordinator.py prepare \
    --run_id "${run_id}" \
    --ssh_host RHINO_HOST_OR_ALIAS \
    --ssh_user FHCC_USER \
    --remote_run_root protocol_chipseq_signal_norm_slurm_runs \
    --partition FHCC_PARTITION \
    --account FHCC_ACCOUNT \
    --env_name env_protocol \
    --remote_python /absolute/path/to/env_protocol/bin/python3 \
    --result_destination artifacts/slurm/results
```

Omit `--ssh_user` if the SSH configuration already selects it. `prepare` refuses to replace an existing run directory.

For a write-free inventory before preparation:
```bash
python3 tests/integration/slurm/coordinator.py prepare \
    --run_id "${run_id}" \
    --ssh_host RHINO_HOST_OR_ALIAS \
    --partition FHCC_PARTITION \
    --account FHCC_ACCOUNT \
    --remote_python /absolute/path/to/env_protocol/bin/python3 \
    --dry_run
```

Inspect local state and render the network commands without executing them:
```bash
python3 tests/integration/slurm/coordinator.py status --run_id "${run_id}"
python3 tests/integration/slurm/coordinator.py push --run_id "${run_id}" --dry_run
python3 tests/integration/slurm/coordinator.py instructions --run_id "${run_id}"
python3 tests/integration/slurm/coordinator.py pull --run_id "${run_id}" --dry_run
```

<br />

## Source bundle format
The local session is `artifacts/slurm/sessions/<run-id>/`. Its `incoming/` directory is the complete transfer unit:
```text
incoming/
├── launcher/
│   ├── src/protocol_chipseq_signal_norm/utilities/utils_cli.py
│   └── tests/integration/slurm/coordinator.py
├── remote_config.json
├── source.tar.gz
├── source.tar.gz.sha256
├── source_manifest.json
└── transfer_manifest.json
```

`source.tar.gz` has deterministic gzip and tar metadata and contains:
```text
bundle/
├── manifest.json
├── source_checksums.sha256
└── source/<exact repository paths>
```

The source inventory is the sorted union of tracked and nonignored untracked files in the declared runtime closure: `bin/`, `lib/bash/`, `src/protocol_chipseq_signal_norm/`, `tests/integration/slurm/`, `tests/fixtures/slurm/`, and `install/envs/env_protocol.yml`. Thus current dirty submit/bootstrap dependencies and required untracked runner files are captured while unrelated data, reports, tests, and documentation are not copied. It excludes `.git`, ignored caches, `__pycache__`, `.pyc`, `.DS_Store`, `artifacts/` and prior bundle sessions, and all private-repository content. The manifest records every source path, type, mode, size, SHA-256 digest, public Git HEAD, NUL-safe porcelain status as JSON text, porcelain digest, staged-state digest, working-tree diff digest, creation timestamp, runner version, run ID, and declared validation scope.

The source archive checksum is recorded both in `source.tar.gz.sha256` and `remote_config.json`. `transfer_manifest.json` checks every transferred file other than itself. The launcher refuses unsafe tar paths, hard links, devices, FIFOs, escaping symlinks, checksum drift, unexpected transfer files, or a run-directory/run-ID mismatch.

<br />

## Transfer and remote launch
Execute the prepared push only when ready to contact FHCC:
```bash
python3 tests/integration/slurm/coordinator.py push --run_id "${run_id}"
```

This creates only `<remote-run-root>/<run-id>/incoming/` and rsyncs the prepared transfer unit there. It does not use `--delete`.

Print the authoritative remote commands:
```bash
python3 tests/integration/slurm/coordinator.py instructions --run_id "${run_id}"
```

Log in manually, then run the printed checksum and launch commands. Their shape is:
```bash
ssh FHCC_USER@RHINO_HOST_OR_ALIAS
(cd protocol_chipseq_signal_norm_slurm_runs/7i-slurm-bundle-20260716/incoming && \
    sha256sum -c source.tar.gz.sha256)
PYTHONDONTWRITEBYTECODE=1 RUN_SLURM=1 WAIT_SLURM=1 CONFIRM_SLURM_WET=1 \
    /absolute/path/to/env_protocol/bin/python3 \
        protocol_chipseq_signal_norm_slurm_runs/7i-slurm-bundle-20260716/incoming/launcher/tests/integration/slurm/coordinator.py \
        remote-launch \
        --config protocol_chipseq_signal_norm_slurm_runs/7i-slurm-bundle-20260716/incoming/remote_config.json \
        --archive protocol_chipseq_signal_norm_slurm_runs/7i-slurm-bundle-20260716/incoming/source.tar.gz
```

All three gate values must be exactly `1`. `RUN_SLURM` retains its established meaning of enabling scheduler execution, `WAIT_SLURM` retains its established meaning of waiting for completion, and `CONFIRM_SLURM_WET` is the additional deliberate real-submission confirmation. The shell runner and Python runner both enforce the gate before any `sbatch` call.

The launcher verifies the transfer, extracts into `<remote-run-root>/<run-id>/source_bundle/`, verifies exact source membership and hashes, writes a run-local runtime configuration, and then invokes the dedicated wet runner. It never modifies the normal remote checkout.

<br />

## Remote preflight and status
Before submission, `preflight.json` records the hostname, operating system, Bash path/version, Python path/version, paths for `sbatch`, `squeue`, `sacct`, and `scontrol`, Conda/Mamba environment availability, archive and exact source-tree checksums, writable run/result directories, partition, account, resource limits, and derived-fixture preparation. The runner fails without submitting when any required check fails. Each job also passes `sbatch --test-only` before its real submission so the selected account, partition, and resource request are scheduler-validated.

While logged in, inspect structured status with the exact command printed by `instructions`, or inspect files directly:
```bash
/absolute/path/to/env_protocol/bin/python3 \
    protocol_chipseq_signal_norm_slurm_runs/7i-slurm-bundle-20260716/incoming/launcher/tests/integration/slurm/coordinator.py \
    status \
    --session_dir protocol_chipseq_signal_norm_slurm_runs/7i-slurm-bundle-20260716

find protocol_chipseq_signal_norm_slurm_runs/7i-slurm-bundle-20260716/results \
    -maxdepth 2 -type f -print
```

The coordinator waits on `sacct` terminal evidence, using `squeue` only for live observation. A numeric job ID from successful `sbatch` acceptance is never treated as job success.

<br />

## Result format
The remote result directory has fixed names:
```text
results/
├── artifacts/
├── checksums.sha256
├── commands.log
├── exit_status.txt
├── jobs.json
├── live_status.json
├── preflight.json
├── run_manifest.json
├── stderr/
├── stdout/
├── summary.json
└── summary.txt
```

Every job record contains the exact submission command, job ID and unique name, requested resources, submit/start/finish timestamps, final scheduler state, exit code, stdout/stderr paths, output assertions, cleanup result, and derived success value. The alignment assertions include `samtools quickcheck` and an expected read-name check in addition to nonempty BAM/BAI checks. `commands.log` retains preflight fixture commands, `sbatch --test-only`, real submission, `squeue`, `sacct`, and command-based assertion calls. `checksums.sha256` covers every result file except itself with exact membership; unexpected or omitted files fail local verification.

<br />

## Pull and local verification
Back on the Mac:
```bash
python3 tests/integration/slurm/coordinator.py pull --run_id "${run_id}"
python3 tests/integration/slurm/coordinator.py verify --run_id "${run_id}"
```

The verifier checks the result run ID against the prepared session, compares the returned source-manifest digest with the prepared bundle, requires every fixed result path, validates exact result checksum membership and digests, requires the full declared job inventory, rejects duplicate or nonterminal jobs, recomputes job success from state/exit/assertions, compares that evidence with `summary.json` and `exit_status.txt`, and returns nonzero for any failed, missing, malformed, or unverifiable job.

For later Codex review, provide the run ID and the pulled directory:
```text
artifacts/slurm/results/<run-id>/
```

Codex needs only the local prepared session and pulled result tree; no Codex process is needed on FHCC.

<br />

## Confined local cleanup
Preview and then remove only a marked coordinator-created local session:
```bash
python3 tests/integration/slurm/coordinator.py clean --run_id "${run_id}" --dry_run
python3 tests/integration/slurm/coordinator.py clean --run_id "${run_id}"
```

Cleanup validates the run ID, proves the target is a direct child of the configured bundle directory, and requires a matching coordinator marker. It does not remove pulled results stored elsewhere and never cleans remote paths.
