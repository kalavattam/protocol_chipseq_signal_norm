
# Repository Guidelines
## Project Structure and Module Organization
Workflow entrypoints live in `scripts/`: use `execute_*.sh` for orchestration and `submit_*.sh` for per-job work. Shared Bash functions are in `scripts/functions/`; Python signal-processing utilities are in `scripts/` and `scripts/functions/`. Installation support is under `install/`, with Conda YAML files in `install/envs/` and bootstrap scripts in `install/scripts/`.

Smoke tests live in `tests/scripts/smoke/`, with shared helpers in `tests/scripts/lib/`. Workflow-specific fixture recipes are in `tests/<workflow>/scripts/make_fixtures.sh`. Generated fixture outputs are not source files; keep only the hand-written `tests/*/fixtures/README.md` docs tracked.

<br />

## Build, Test, and Development Commands
There is no compilation step. Create the main environment with:
```bash
sh install/scripts/install_envs_entrypoint.sh --env_nam env_protocol --yes
```

Run the default smoke suite with:
```bash
bash tests/scripts/run_smoke_tests.sh
```

Syntax-check a changed shell script before committing:
```bash
bash -n scripts/execute_align_fastqs.sh
```

Generate missing fixtures explicitly when needed, for example:
```bash
bash tests/download_fastqs/scripts/make_fixtures.sh
```
<br />

## Coding Style and Naming Conventions
Target Bash >= 4.4. Use four-space indentation, `snake_case` variables and functions, and `set -euo pipefail` where appropriate. Use two empty lines between shell function definitions. Prefer existing short variable idioms where clear, such as `fil_in`, `fil_out`, `dir_out`, `arr_*`, `nam_job`, and `idx`.

Shell help text is the explicit indentation exception: help heredocs follow `docs/style/help.md`, which uses two spaces per indentation level. For wrapped Bash command arrays, indent option lines under the invoked script or command:
```bash
arr_cmd=(
    "${TEST_BASH}" "${scr}"
        --arg_1 "${val_1}"
        --arg_2 "${val_2}"
)
```

Prefer keeping source code, code comments, diagnostics, and command-construction lines at or under 80 characters where practical; this applies across Bash, Python, R, and test code. Markdown prose is not subject to this line-length preference: do not hard-wrap Markdown prose solely for line length, and let the renderer or editor wrap natural paragraphs. Split long shell diagnostics as adjacent quoted string arguments at natural phrase boundaries. Do not split inside command substitutions, variable expansions, compact parser alias patterns, or other syntax where wrapping would reduce clarity. Shell help text is also exempt from this rule and is governed by `docs/style/help.md`. The smoke test `test_shell_line_length.sh` reports long non-help shell lines as advisory warnings.

For multiline smoke-test helper calls, put the helper name on its own line and give each helper positional argument its own line:
```bash
assert_pattern_found \
    "${file}" \
    "${pattern}" \
    "label"
```

Indent `case` patterns one level under the `case ... in` line, and indent each pattern body one additional level:
```bash
case "${mode}" in
    signal)
        mode="signal"
        ;;
    *)
        return 1
        ;;
esac
```

For `run_capture`, split the label and log/output path onto separate lines; the captured command and its ordinary command arguments may stay grouped when that is clearer. Keep final human-readable labels as one quoted argument when reasonably short. Do not split command substitutions or variable expansions inside labels; split long labels only at natural phrase boundaries for helpers that safely collect final label text with `$*`.

In Markdown headings and prose, write `and` rather than `&`, and use the Oxford/serial comma.

Bash wrappers should use explicit entrypoint functions: keep argument parsing in `parse_args()` and final control flow in `main()`. Prefer the sequence: version guard, safe mode, path resolution and helper sourcing, local helper functions, `init_defs`, `parse_args`, `canonicalize_args`, `validate_args`, vector preparation/validation, execution setup, environment activation, then dispatch.

Order shell functions from low-level helpers toward lifecycle orchestration. Put script-local help near `init_defs()` and `parse_args()`, keep runtime helpers after validation/preparation helpers, and keep `main()` last.

Keep `main()` as lifecycle orchestration. When setup, validation, debug printing, environment activation, vector reconstruction, or execution dispatch becomes nontrivial, extract those blocks into named helpers and have `main()` call them in order. Helpers should return nonzero on failure; `main()` should decide whether to return or exit.

Do not place live default/global assignments between function definitions. Prefer setting defaults in `init_defs()` and calling it at the start of `main()` before help output or argument parsing. Define `init_defs()` near the argument lifecycle functions, before `parse_args()` when practical. For very small scripts, a compact default block immediately before `main "$@"` is acceptable.

When a wrapper has both hardcoded execution controls and user-facing argument defaults, split them into `init_args_hardcoded()` and `init_arg_defs()`, with `init_defs()` calling both in that order.

Keep `execute_*.sh` responsible for orchestration: deriving output paths, splitting work, and dispatching serial, GNU Parallel, or Slurm jobs. Keep `submit_*.sh` responsible for worker-level execution: validating direct inputs, constructing Python/tool calls, and writing per-task logs. Build commands as arrays, keep dry-run output explicit, and avoid reading mode-specific globals unless that mode is active under `set -u`.

Keep wrapper changes surgical: preserve existing `submit_` / `execute_` boundaries and source shared helpers instead of duplicating logic. Use section comments for major phases, `#+` for wrapped explanatory comments, and local `show_help` heredocs where a helper needs its own usage text. For detailed shell-function help, prefer `Usage:`, `Description:`, `Arguments:` or `Positional arguments:`, `Expected globals:` or `Expected inputs:` when applicable, `Returns:`, then `Notes:`. Use `Expected globals:` only for helpers that intentionally read or write important global state; group entries by role, such as required scalar globals, optional scalar globals, reconstructed arrays, and generated global outputs. Nontrivial shell-function help should include a `Returns:` section; use hyphen bullets when that section lists multiple distinct items. Very small obvious helpers may stay terse. Use `git diff --check` to catch whitespace errors. Do not commit generated fixtures or `tests/outputs/`.

<br />

## Detailed Style References
Treat this file as the top-level contributor and agent guide. For detailed shell help-text rules, consult `docs/style/help.md` before editing `show_help`, `detail_*`, `help_*`, or heredoc-based CLI documentation.

For smoke-test structure, fixture policy, generated-output handling, and test gates, consult `tests/README.md`. Add new focused style guides under `docs/style/` only when a topic outgrows this file. When a style policy changes, update this guide and the relevant detailed guide in the same patch when practical.

<br />

## Testing Guidelines
Name smoke tests `test_<layer>_<workflow>_<case>.sh`, such as `test_execute_filter_alignments_cram_to_cram.sh`. Run focused tests for changed paths before the full runner. Explicit gates are dependency-specific: `RUN_ATRIA=1`, `RUN_PARALLEL=1`, and `RUN_SLURM=1`. Do not introduce broad integration gates or fake tool shims.

<br />

## Commit and Pull Request Guidelines
Recent commits use concise scoped subjects such as `test(align): add Bowtie2 CRAM smoke coverage` and `fix(execute): skip parallel detection for serial local runs`. Keep commits focused. PRs should describe behavior changes, list validation commands, note required environments or gates, and call out deferred work.
