
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
Bash wrappers and nontrivial utility scripts should use explicit entrypoint functions. Keep argument parsing in `parse_args()` and final control flow in `main()`.

Use the compute-signal wrappers as the current template for complex `execute_*.sh` and `submit_*.sh` scripts. Do not force tiny helper scripts into the full wrapper shape when a compact structure is clearer.

Prefer this bootstrap and lifecycle sequence:
```txt
Bash version guard
safe mode
path resolution
minimal bootstrap parsing only if needed
helper sourcing
local helper functions
default initialization
argument parsing
canonicalization
validation
preparation
execution setup
environment activation
dispatch
main
```

For `execute_*.sh` wrappers, the usual `main()` order is:
```txt
init_defs
source_helpers_execute
help handling
parse_args
canonicalize_args, if needed
validate_args
prepare_vecs and validate_vecs, if needed
config_exec, if needed
setup_env
check_tools
debug/state printing, if enabled
run_jobs
```

For `submit_*.sh` wrappers that need user-supplied `--dir_scr` before helpers can be sourced, keep the bootstrap path narrow:
```txt
init_defs
early help handling
resolve_dir_scr
source_helpers_submit
parse_args
canonicalize_args, if needed
validate_args
prepare_vecs and validate_vecs, if needed
debug/state printing, if enabled
setup_env
check_tools
run_jobs
```

Use only the lifecycle helpers that correspond to real work in the script. Do not add empty `canonicalize_args()`, `config_exec()`, `setup_env()`, or similar functions just to satisfy the template.

Order shell functions from low-level helpers toward lifecycle orchestration. Put script-local help near `init_defs()` and `parse_args()`, keep runtime helpers after validation and preparation helpers, and keep `main()` last.

Keep `main()` as lifecycle orchestration. When setup, validation, debug printing, environment activation, vector reconstruction, or execution dispatch becomes nontrivial, extract those blocks into named helpers and have `main()` call them in order.

Inside functions, including `main()`, prefer `return 0` and `return 1` over `exit`. Reserve `exit` for top-level code outside functions, especially Bash-version guards and fatal bootstrap checks before `main()` is available. Let the final top-level `main "$@"` call propagate the returned status under `set -euo pipefail`.

Put helper sourcing in a named helper for executable scripts:
- use `source_helpers_execute()` in `execute_*.sh` wrappers,
- use `source_helpers_submit()` in `submit_*.sh` wrappers, and
- use `source_helpers_script()` in utility scripts.

Submit wrappers may need a small bootstrap exception before the ordinary lifecycle when parser helpers depend on a user-supplied scripts directory. In that case, use narrowly named bootstrap helpers, such as `resolve_dir_scr()` and `source_helpers_submit()`, before normal argument parsing.

Execute wrappers usually source help functions from `scripts/functions/help/`, so call `source_helpers_execute()` near the start of `main()` before help handling, normal argument parsing, and validation. Utility scripts may call `source_helpers_script()` before help handling when help or parser helpers are sourced.

Keep bootstrap parsing minimal. It should only collect information required to source helpers or establish the script environment, such as `--dir_scr`; ordinary options belong in `parse_args()`.

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

Keep `submit_*.sh` responsible for worker-level execution: validating direct inputs, selecting a processing primitive, and writing per-task logs.

Execute wrappers should prepare and validate vectors before dispatch, then serialize arrays to comma-delimited values only when constructing downstream submit commands.

Submit wrappers should reconstruct comma-delimited inputs into arrays, validate direct worker inputs, and handle worker-specific context such as `SLURM_ARRAY_TASK_ID`.

Submit wrappers that support both Slurm-array task execution and local iteration should keep worker dispatch in `run_jobs`; `main()` should call that helper after parsing, validation, setup, and vector reconstruction.

Wrappers should check user-visible tools required by the selected execution path, even when the actual external-tool command is assembled and run by a sourced processing helper.

Build commands as arrays, keep dry-run output explicit, and avoid reading mode-specific globals unless that mode is active under `set -u`.

Keep wrapper changes surgical: preserve existing `submit_` / `execute_` boundaries and source shared helpers instead of duplicating logic.

<br />

## Wrapper and Processing Boundaries
Wrappers orchestrate work. Sourced functions perform substantive processing. This is the practical version of "no processing logic in `execute` or `submit`."

Wrappers may:
- parse, canonicalize, and validate arguments;
- derive output paths, sample names, log paths, and task indexes;
- reconstruct or serialize vectors;
- activate environments and check tools;
- choose which named processing primitive to call;
- call Python entrypoints, external tools, or helper functions through a sourced function;
- handle wrapper-level success, failure, dry-run, and dispatch behavior.

Wrappers should not embed:
- substantive domain algorithms;
- long pipelines such as `samtools | awk | sort | gzip`;
- inline AWK, Python, R, or Perl programs;
- multi-step file transformations with temporary files;
- BAM/CRAM, FASTQ, bedGraph, TSV, or metadata-processing logic that can be named and tested as a helper.

Use `scripts/functions/` for processing primitives. Good candidates include alignment sorting, filtering, format conversion, counting, table transformation, and domain-specific calculations. Prefer explicit helper arguments over hidden globals. Pass paths such as `fil_in`, `fil_out`, `ref_fa`, `log_out`, and `log_err` directly when the helper owns command execution or redirection.

When a processing helper owns an external-tool call, it should own command-array construction and redirection. Have the wrapper pass explicit input, output, and log paths; the wrapper should not rely on the helper to reconstruct wrapper-local log names from globals. Append optional command arguments conditionally with arrays instead of command substitution, so paths with spaces remain safe.

A tiny wrapper-local command is acceptable only when it is pure dispatch glue and has no domain behavior beyond invoking a named entrypoint. If the command needs mode branching, format-specific BAM/CRAM handling, temporary files, pipes, an embedded program, or nontrivial cleanup, move it to a sourced function.

For example, do not keep a paired-fragment conversion pipeline inside `submit_*::run_job()`:
```bash
samtools view ... | awk '...' | sort ... | gzip > "${fil_out}"
```

Prefer a named processing primitive:
```bash
convert_alignment_to_bed_awk \
    "${threads}" \
    "${fil_in}" \
    "${fil_out}" \
    "${ref_fa}" \
    "${log_out}" \
    "${log_err}"
```

The wrapper remains responsible for deciding that the AWK branch is active, deriving the output and log paths, and reporting wrapper-level failure context.

<br />

## Shell Test Structure
Shell smoke tests are usually executable recipes, not user-facing command-line programs. Do not force ordinary `test_*.sh` scripts into `main()`, `parse_args()`, or `validate_args()` solely for style.

Every smoke-test script must enable safe mode with `set -euo pipefail`. Use
this order for smoke tests unless a local script has a clear reason to differ:
```txt
shebang, attribution/header comments
set -euo pipefail
TEST_NAME and required test metadata
source test helpers
local helper functions
print the test section/banner
handle gates or clean skips
static fixture/reference paths and arrays
temporary, output, log, and input paths
case-specific deterministic paths/logs/status files/job names/expected outputs
deterministic expected values, such as rows, headers, regex tails, and fragments
static arrays grouped by purpose, including command arrays and case arguments
cleanup and mkdir
environment and fixture validation
runtime-dependent declarations immediately after the step that creates their required values, including arrays
copied/runtime fixture preparation
test cases
assertions near the cases they validate
finish
```

Deterministic top-level declarations should be grouped in the earliest valid declaration section before test execution begins. Like declarations belong with like, and distinct declaration groups should be separated by one blank line. Do not define deterministic top-level variables near first use in a test case when their values are knowable during setup. Deterministic case-specific assertion artifacts should also be declared in this section, including primary outputs, retained part-file paths, derived index files such as `bai_*` and `crai_*`, captured logs such as `log_*`, stats/count files, status files, and expected output paths. Expected headers, rows, regex tails, and output fragments are deterministic assertion artifacts. Do not place them between static arrays and cleanup or validation code. Fixture inputs are not assertion artifacts and do not need extra aliases solely for this rule.

Static arrays are dependency consumers. Place scalar/path declarations and deterministic expected values before static arrays when the arrays and expected values use only earlier declarations. Static arrays should be grouped by purpose, such as fixture manifests, command arrays, case argument arrays, expected-file arrays, and skip lists. Like arrays belong with like arrays, and distinct static-array groups should be separated by one blank line.

Use descriptive case-specific names instead of repeatedly reassigning generic names across cases. For example, prefer `fil_out_se_signal` and `log_se_signal` over reusing `fil_out` and `log` below a case comment. Case comments should remain near the test action, but deterministic variables should not be introduced under those comments unless they are genuinely runtime-derived. Keep assertion calls near the cases they validate. Do not separate an explanatory comment from the command or assertion it directly introduces. Blank lines separate declaration groups and major logical sections, not a case comment from its associated case action.

Place declarations in the earliest section where all values they use have already been defined. Do not move a declaration above the runtime step that creates one of its values. For example, `require_env_project env_nam` defines `env_nam`, so submit-layer `arr_cmd_filter` arrays that include `--env_nam "${env_nam}"` must stay after that call.

Function-local variables, loop/per-iteration variables, and values derived from runtime output may stay near the code that creates them. Within functions, define repeated or clarity-important derived values as `local` variables near the top of the function. This includes paths, labels, expected rows, log names, count files, and command fragments derived from function arguments. Do not repeatedly construct the same derived value inline throughout a function, and do not hoist function-argument-dependent values into top-level test sections.

Local helper functions belong immediately after helper sourcing. Keep a helper local when it is used by one script and reduces repetition or clarifies assertions. If the same helper behavior appears in two or more scripts, promote it to a shared helper and source it from the appropriate helper library.

Use `tests/scripts/lib/test_helpers.sh` for shared smoke-test behavior and `tests/scripts/lib/fixture_helpers.sh` for shared fixture-generation behavior. Do not create shared helpers for superficially similar snippets that have different meaning or different failure semantics.

Use `main()` and phase functions for user-facing or complex test utilities, such as `run_smoke_tests.sh` and `clean_test_outputs.sh`. Consider `main()` for unusually large smoke tests with multiple modes, complex setup, or repeated phases, but keep simple smoke tests as top-to-bottom recipes.

Fixture-generation scripts should follow a similar recipe: Bash guard, safe mode, path resolution, helper sourcing, local helpers, static fixture paths/content declarations, cleanup and directory creation, generation steps, validation, and final report.

<br />

## State and Command Helpers
Keep debug printing, state mutation, command construction, and job dispatch in separate helpers once each block becomes more than a few lines.

Debug-print helpers should be read-only and should return 0 when debugging or verbosity is disabled.

State-conversion helpers, such as vector serialization or reconstruction, should do one transformation and leave validation to a separate validation helper.

Command builders should construct Bash arrays, not command strings. Use a stable argument order: executable, execution controls, mode, inputs, outputs, mode-specific options, output formatting, logging, and job controls.

If a helper intentionally writes a global, document that contract in the helper help. For example, a command builder that writes `cmd_bld` should say that `cmd_bld` is a generated global indexed array.

<br />

## Test Helper Calls
For multiline test helper calls, put the helper name on its own line and give each helper positional argument its own line:
```bash
assert_pattern_found \
    "${file}" \
    "${pattern}" \
    "label"
```

For `run_capture`, split the label and log/output path onto separate lines; the captured command and its ordinary command arguments may stay grouped when that is clearer.

Keep final human-readable labels as one quoted argument when reasonably short. Do not split command substitutions or variable expansions inside labels. Split long labels only at natural phrase boundaries for helpers that safely collect final label text with `$*`. Keep explanatory comments directly attached to the helper call or assertion they introduce.

Treat multiline helper calls as block-like. Place one blank line between adjacent multiline helper-call blocks, but do not add blank lines inside shell control-flow conditions.

<br />

## Function Help
Local `show_help` heredocs are appropriate when a helper needs its own usage text. Detailed shell-function help follows `docs/style/help.md`.

For nontrivial shell-function help, prefer:
```txt
Usage:
Description:
Arguments: / Keyword arguments: / Positional arguments:
Expected globals: / Expected inputs:
Returns:
Notes:
Examples:
```

Use `Expected globals:` only for helpers that intentionally read or write important global state. Group entries by role, such as required scalar globals, optional scalar globals, reconstructed arrays, and generated global outputs.

Nontrivial shell-function help should include a `Returns:` section. Use hyphen bullets when that section lists multiple distinct items. Very small obvious helpers may stay terse.
