
# Shell Style Guide
This file owns detailed Bash style rules for this repository. `AGENTS.md` is the top-level repository guide.

For shell help text, `show_help`, `detail_*`, `help_*`, and heredoc-based CLI documentation, use `docs/style/help.md`.

<br />

## Basic Conventions
Target Bash >= 4.4. Use `set -euo pipefail` where appropriate.

Use four-space indentation for shell code. Shell help text is the explicit indentation exception: help heredocs follow `docs/style/help.md`, which uses two spaces per indentation level.

Use `snake_case` variables and functions. Prefer existing short variable idioms where clear, such as `fil_in`, `fil_out`, `dir_out`, `dir_eo`, `arr_*`, `nam_job`, and `idx`.

Use `dir_eo` for directories that collect stderr/stdout logs, paired with the user-facing option `--dir_eo`. The short option may remain `-eo`. Older `err_out` names may still exist in scripts that have not yet been migrated, but new or actively refactored wrappers should prefer `dir_eo`.

Use two empty lines between shell function definitions.

Use section comments for major phases and `#+` for wrapped explanatory comments.

<br />

## Line Length
Prefer keeping source code, code comments, diagnostics, and command-construction lines at or under 80 characters where practical. This applies across production and test code.

Markdown prose is not subject to this line-length preference. Do not hard-wrap Markdown prose solely for line length, and let the renderer or editor wrap natural paragraphs.

User-facing shell help text (e.g., that in heredocs) is also exempt from the source-code line-length preference and is governed by `docs/style/help.md`.

Split long shell diagnostics as adjacent quoted string arguments at natural phrase boundaries. Do not split inside command substitutions, variable expansions, compact parser alias patterns, or other syntax where wrapping would reduce clarity.

The smoke test `test_shell_line_length.sh` reports long non-help shell lines as advisory warnings.

<br />

## Indentation Patterns
For wrapped Bash command arrays, indent option lines under the invoked script or command:
```bash
arr_cmd=(
    "${TEST_BASH}" "${scr}"
        --arg_1 "${val_1}"
        --arg_2 "${val_2}"
)
```

For `sbatch` command arrays, `sbatch` is the invoked command:
```bash
cmd_slurm=(
    sbatch
        --job-name="${nam_job}"
        --nodes=1
        --cpus-per-task="${threads}"
        "${cmd_bld[@]}"
)
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

<br />

## Heredocs
Use `EOM` as the standard delimiter for Bash heredocs:
```bash
cat << EOM
literal text
EOM
```

Quoted forms are fine when expansion should be disabled:
```bash
cat << 'EOM'
literal text
EOM
```

Interpreter stdin heredocs may use a language-specific delimiter when that is clearer. For example, use `PY` for Python code passed to Python stdin.

<br />

## Wrapper Structure
Bash wrappers should use explicit entrypoint functions. Keep argument parsing in `parse_args()` and final control flow in `main()`.

Use the compute-signal wrappers as the current template for complex `execute_*.sh` and `submit_*.sh` scripts. Do not force tiny helper scripts into the full wrapper shape when a compact structure is clearer.

Prefer this sequence:
```txt
version guard
safe mode
path resolution and helper sourcing
local helper functions
init_defs
parse_args
canonicalize_args
validate_args
vector preparation and validation
execution setup
environment activation
dispatch
main
```

Order shell functions from low-level helpers toward lifecycle orchestration. Put script-local help near `init_defs()` and `parse_args()`, keep runtime helpers after validation and preparation helpers, and keep `main()` last.

Keep `main()` as lifecycle orchestration. When setup, validation, debug printing, environment activation, vector reconstruction, or execution dispatch becomes nontrivial, extract those blocks into named helpers and have `main()` call them in order.

Inside functions, including `main()`, prefer `return 0` and `return 1` over `exit`. Reserve `exit` for top-level code outside functions, especially Bash-version guards and fatal bootstrap checks before `main()` is available. Let the final top-level `main "$@"` call propagate the returned status under `set -euo pipefail`.

Submit wrappers may need a small bootstrap exception before the ordinary lifecycle when parser helpers depend on a user-supplied scripts directory. In that case, use narrowly named bootstrap helpers, such as `resolve_dir_scr()` and `source_submit_helpers()`, before normal argument parsing.

Do not place live default/global assignments between function definitions. Prefer setting defaults in `init_defs()` and calling it at the start of `main()` before help output or argument parsing. Define `init_defs()` near the argument lifecycle functions, before `parse_args()` when practical. For very small scripts, a compact default block immediately before `main "$@"` is acceptable.

When a wrapper has both hardcoded execution controls and user-facing argument defaults, split them into `init_args_hardcoded()` and `init_arg_defs()`, with `init_defs()` calling both in that order.

Top-level indexed arrays declared with `declare -a` must become `declare -ga` when moved into a function (lifecycle or otherwise) if later helpers need to read or mutate them. Use plain `declare -a` only for arrays that are intentionally function-local.

Use stable lifecycle helper names when the behavior exists:
- `config_exec` configures serial, GNU Parallel, or Slurm execution.
- `setup_env` activates the requested environment and updates environment variables such as `PYTHONPATH`.
- `check_tools` checks executable dependencies required by the selected execution path.
- `print_state_debug` and `print_vecs_debug` print scalar and vector state without mutating it.
- `prepare_vecs` reconstructs or derives arrays from scalar inputs.
- `validate_vecs` validates reconstructed arrays and cross-array constraints.
- `run_jobs` is the final dispatch step.

Keep public API renames separate from structure-only lifecycle refactors. For example, do not migrate `--err_out` to `--dir_eo` in the same patch unless the user explicitly asks for an option rename. *(`#TODO`: come back to adjust this rule later.)*

<br />

## Execute and Submit Boundaries
Keep `execute_*.sh` responsible for orchestration: deriving output paths, splitting work, and dispatching serial, GNU Parallel, or Slurm jobs.

Keep `submit_*.sh` responsible for worker-level execution: validating direct inputs, constructing Python/tool calls, and writing per-task logs.

Execute wrappers should prepare and validate vectors before dispatch, then serialize arrays to comma-delimited values only when constructing downstream submit commands.

Submit wrappers should reconstruct comma-delimited inputs into arrays, validate direct worker inputs, and handle worker-specific context such as `SLURM_ARRAY_TASK_ID`.

Submit wrappers that support both Slurm-array task execution and local iteration should keep worker dispatch in `run_jobs`; `main()` should call that helper after parsing, validation, setup, and vector reconstruction.

Build commands as arrays, keep dry-run output explicit, and avoid reading mode-specific globals unless that mode is active under `set -u`.

Keep wrapper changes surgical: preserve existing `submit_` / `execute_` boundaries and source shared helpers instead of duplicating logic.

<br />

## State and Command Helpers
Keep debug printing, state mutation, command construction, and job dispatch in separate helpers once each block becomes more than a few lines.

Debug-print helpers should be read-only and should return 0 when debugging or verbosity is disabled.

State-conversion helpers, such as vector serialization or reconstruction, should do one transformation and leave validation to a separate validation helper.

Command builders should construct Bash arrays, not command strings. Use a stable argument order: executable, execution controls, mode, inputs, outputs, mode-specific options, output formatting, logging, and job controls.

If a helper intentionally writes a global, document that contract in the helper help. For example, a command builder that writes `cmd_bld` should say that `cmd_bld` is a generated global indexed array.

<br />

## Smoke-Test Helper Calls
For multiline smoke-test helper calls, put the helper name on its own line and give each helper positional argument its own line:
```bash
assert_pattern_found \
    "${file}" \
    "${pattern}" \
    "label"
```

For `run_capture`, split the label and log/output path onto separate lines; the captured command and its ordinary command arguments may stay grouped when that is clearer.

Keep final human-readable labels as one quoted argument when reasonably short. Do not split command substitutions or variable expansions inside labels. Split long labels only at natural phrase boundaries for helpers that safely collect final label text with `$*`.

<br />

## Function Help
Local `show_help` heredocs are appropriate when a helper needs its own usage text. Detailed shell-function help follows `docs/style/help.md`.

For nontrivial shell-function help, prefer:
```txt
Usage:
Description:
Keyword arguments: / Positional arguments:
Expected globals: / Expected inputs:
Returns:
Notes:
```

Use `Expected globals:` only for helpers that intentionally read or write important global state. Group entries by role, such as required scalar globals, optional scalar globals, reconstructed arrays, and generated global outputs.

Nontrivial shell-function help should include a `Returns:` section. Use hyphen bullets when that section lists multiple distinct items. Very small obvious helpers may stay terse.
