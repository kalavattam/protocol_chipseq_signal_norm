# Repository Guidelines
## Read This First
Use the repository as the source of truth. Inspect current files before editing, and do not let old chat context override the working tree. Keep patches small and focused.

Do not touch generated fixture outputs, `tests/outputs/`, manuscript/workflow notes, or troubleshooting drafts unless the user explicitly asks for that work.

<br />

## Project Layout
Workflow entrypoints live in `scripts/`: `execute_*.sh` orchestrates work, and `submit_*.sh` handles per-job execution. Shared Bash functions are in `scripts/functions/`; Python utilities are in `scripts/` and `scripts/functions/`.

Smoke tests live in `tests/scripts/smoke/`, with helpers in `tests/scripts/lib/`. Fixture recipes live in `tests/<workflow>/scripts/`. Generated fixture outputs stay ignored; hand-written fixture READMEs stay tracked.

<br />

## Common Commands
Create the main environment, run smoke tests, and syntax-check shell scripts:
```bash
sh install/scripts/install_envs_entrypoint.sh --env_nam env_protocol --yes
bash tests/scripts/run_smoke_tests.sh
bash -n scripts/execute_align_fastqs.sh
```

<br />

## Style and Tooling
This file is the routing layer, not the full style manual. Use:
- `docs/style/shell.md` for Bash structure, wrapper organization, and shell formatting.
- `docs/style/help.md` for shell help text and heredoc-based CLI docs.
- `docs/style/python.md` for Python script structure and output formatting.
- `tests/README.md` for smoke-test structure, gates, fixtures, and generated output policy.

Prefer deterministic checks over agent memory. When a style rule becomes important, document it in `docs/style/` and add an advisory or enforced check when practical. Keep Markdown prose natural; do not hard-wrap it only for source line-length preferences.

<br />

## Testing Expectations
Run focused tests for changed paths before broader suites. Use existing gates: `RUN_ATRIA=1`, `RUN_PARALLEL=1`, `RUN_SLURM=1`, and `SLURM_WAIT=1`. Do not add fake dependency shims or broad integration gates.

For shell changes, run `bash -n` on changed scripts and `git diff --check`.

<br />

## Commits
Use focused commits. The first line should follow `verb(subject): description`, for example `fix(trim): validate discovered SE output`. Follow-up message paragraphs should be complete sentences grouped by context or subject. Commit generated fixture outputs only if the user explicitly changes the fixture policy. When Codex contributes materially to a commit, include:
```txt
Co-authored-by: Codex <codex@openai.com>
```
