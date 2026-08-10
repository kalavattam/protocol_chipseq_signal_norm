
# Generated artifacts
All repository-generated, non-source outputs live under this boundary and are ignored by Git.

- `tests/` contains test logs, temporary products, and preserved legacy test outputs.
- `dev/audit/` contains generated audit reports.
- `slurm/sessions/` contains immutable prepared transfer bundles.
- `slurm/results/` contains checksummed pulls from isolated remote runs.
- `reviews/` contains review evidence. It is **not** regenerable: a unit test imports a validator from `reviews/2026-07-30_post_pilot_standards_phase_c/follow_up`, so this tree must not be cleared with the rest. Pass `T1` removes that dependency, after which it can go.

Everything else here regenerates from the checkers and test runs, so it may be deleted to reclaim space.

Roadmap planning documents no longer live here. They moved to the private repository on 2026-08-10 and are version-controlled there; nothing in this repository may reference them.
