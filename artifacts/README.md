
# Generated artifacts
All repository-generated, non-source outputs live under this boundary and are ignored by Git.

- `tests/` contains test logs, temporary products, and preserved legacy test outputs.
- `dev/audit/` contains generated audit reports.
- `slurm/sessions/` contains immutable prepared transfer bundles.
- `slurm/results/` contains checksummed pulls from isolated remote runs.

Do not place executable source, test definitions, fixture recipes, or documentation here. Preserve validation evidence unless its owner explicitly authorizes cleanup.
