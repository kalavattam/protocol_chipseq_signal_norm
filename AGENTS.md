
# Repository guidelines
## Read this first
Use the repository as the source of truth. Inspect current files before editing, and do not let old chat context override the working tree. Keep patches small and focused.

Do not touch generated fixture outputs, `artifacts/`, manuscript/workflow notes, or troubleshooting drafts unless the user explicitly asks for that work.

<br />

## Project layout
Maintained shell entrypoints live in `bin/`; sourced Bash lives under `lib/bash/` by responsibility. Installable Python CLIs and utilities live under `src/protocol_chipseq_signal_norm/`.

Tests are organized semantically under `tests/unit/`, `tests/contract/`, and `tests/integration/`; helpers live in `tests/support/`. Fixture recipes and documentation live in `tests/fixtures/<workflow>/`. Generated fixtures and all run evidence stay ignored under the documented fixture paths and `artifacts/`.

<br />

## Common commands
Create the main environment, run the safe test suite, and syntax-check shell scripts:
```bash
sh install/scripts/install_envs_entrypoint.sh --env_nam env_protocol --yes
bash tests/run_tests.sh
bash -n bin/execute_align_fastqs.sh
conda run -n env_protocol ruff check --no-fix src tests dev
conda run -n env_protocol ruff format --check src tests dev
```

<br />

## Standards and tooling
This file is the routing layer, not the full standards manual. Start with the ordinary task router in `docs/standards/README.md`, then load only the owner sections and focused commands named for the task. Load a complete domain standard only when the work spans that complete domain or the router has no applicable task:
- `docs/standards/governance.md` owns rule lifecycle, governed change, anti-accretion, and lossless semantic movement.
- `docs/standards/markdown.md` owns maintained Markdown.
- `docs/standards/source_layout.md` owns cross-language source semantics.
- `docs/standards/help.md` owns shared help and documentation semantics.
- `docs/standards/output.md` owns output contracts.
- Active language owners (`shell.md` and `python.md`) own language realization; `docs/standards/testing.md` and `tests/README.md` own test taxonomy and gates.

Prefer deterministic checks over agent memory. When a standard becomes important, document it in `docs/standards/` and add an advisory or enforced check when practical. Keep Markdown prose natural; do not hard-wrap it only for source line-length preferences.

For governed multi-domain changes, approve the normative owner before changing registries, applicability, automation, fixtures, contracts, or maintained realizations. Route through [`GOV.CHANGE.GOLDEN_FIRST`](docs/standards/governance.md#authoritative-standard-first-changes-govchangegoldenfirst), perform its anti-accretion review, record both old-to-new and new-to-old preservation, validate successive owner-specific passes, and finish with all-owner reconciliation.

<br />

## Testing expectations
Run focused tests for changed paths before broader suites. Use existing gates: `RUN_ATRIA=1`, `RUN_PARALLEL=1`, `RUN_SLURM=1`, and `WAIT_SLURM=1`. Do not add fake dependency shims or broad integration gates.

For shell changes, run `bash -n` on changed scripts and `git diff --check`.

<br />

## Commits
Use focused commits. The first line should follow `verb(subject): description`, for example `fix(trim): validate discovered SE output`. Follow-up message paragraphs should be complete sentences grouped by context or subject. Commit generated fixture outputs only if the user explicitly changes the fixture policy. When Codex contributes materially to a commit, include:
```txt
Co-authored-by: Codex <codex@openai.com>
```
