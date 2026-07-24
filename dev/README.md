
# Cleanup audit MVP
Run the audit from the repository root with explicit public and private roots:
```bash
PYTHONDONTWRITEBYTECODE=1 \
    conda run -n env_protocol python dev/audit/run.py \
    --public-root "${PWD}" \
    --private-root ../protocol_chipseq_signal_norm_private
```

Reports are generated under `artifacts/dev/audit/` and ignored by Git. The tool inventories every current dirty public path, preserves the immutable original cleanup cohort in `dev/config/baseline_cleanup_cohort.json`, and does not edit audit targets.

Strict mode runs only declared safe adapters. Existing smoke checks that write `artifacts/tests/` are reported as unavailable rather than run. The generated prompts are review-only and never invoke patches or autofixes.

The run summary separates strict-safe coverage from rules that are unavailable by design or unavailable because their required environment is absent. A completed report with zero findings means only that the strict-safe transactions that executed produced no findings; it is not a repository-wide cleanup verdict.

`--verify` accepts a report as fresh only when its source `run.json` says `completed` and all recorded baseline, manifest, evidence, ledger, finding, and prompt fingerprints still match. Use `--allow-partial` only to inspect an aborted bundle; its verification summary remains explicitly partial.

`--verify-read-only REPORT` performs the same freshness and package-integrity checks in memory without writing verification artifacts or changing package state. A configured package may declare at most two linked children; generate one with `--package-child PACKAGE_ID --umbrella-run-id RUN_GROUP_ID`. Verify the ordered pair with `--verify-linked-pair-read-only CHILD1_REPORT CHILD2_REPORT`; this additionally checks the shared identities, ownership partitions, graph and configuration fingerprints, rule scopes, evidence coverage, and five-artifact package structure without mutating either report.

The two configured download-FASTQs dependency-closure children use a common 16,000-line semantic-Markdown ceiling. The former 12,500-line ceiling became counterproductive during Runtime-requirements remediation: complete source-grounded documentation reduced the remaining capacity to an unsafe margin and encouraged deletion or compression of valid semantic evidence. The 16,000-line ceiling preserves fail-closed behavior while accommodating complete Runtime-requirements evidence through the remaining 7h3c cohorts; it is scoped to these children, not a permanent universal value. The retained final-v3 reports measure 12,377 lines for runtime/production and 13,068 lines for tests/registration, leaving 3,623 lines (22.64% of the ceiling; 29.27% growth over the retained size) and 2,932 lines (18.33% of the ceiling; 22.44% growth over the retained size) of historical headroom, respectively. The current configuration has expanded since those reports, so a newly generated audit is required before describing current rendered headroom. The independent byte ceiling remains unchanged at 1,048,576 bytes for both children.

Target-family pilots use either repeated `--path PATH` arguments or a versioned `--paths-from` JSON selection. The latter records primary/supporting roles and declared context while preserving the full dirty-file ledger. Selected tracked paths use path-scoped staged and worktree `git diff --check`; selected untracked or clean paths use a `git diff --no-index --check` new-file comparison and are reported with that independent coverage label. Pilot reports version their target, fact, policy-question, limitation, semantic-review, and summary artifacts. A generated semantic-review bundle is review-only and includes untracked target source as explicit `new_file_source` evidence.

`dev/config/command_names.json` is the authoritative registry for exact callable spellings. Registry-backed adapters keep callable names distinct from conceptual project or language names, preserve mixed-case commands, and route unknown or ambiguous runtime-command references to semantic review without guessed normalization.

Shell-interface extraction records raw parser observations separately from resolved visibility policy. `Usage` advertises canonical long names only. `Parameters` lists every public short and long alias, including `-h` where settled interface evidence makes it public. Systematic underscore-to-hyphen long-option aliases and the legacy `--hlp` spelling remain indefinitely retained, intentionally undocumented compatibility interfaces.

The pilot Markdown is a concise human review view, not a raw NDJSON dump. It groups evidence by topic, fences every source/diff/JSON/command excerpt, and links the retained machine-readable artifacts for complete detail.

Ruff is available through `env_protocol`, but repository configuration, rule selection, formatting, and enforcement are deferred to a dedicated follow-up session. Audit validation must not treat Ruff as a pilot-closure acceptance gate.
