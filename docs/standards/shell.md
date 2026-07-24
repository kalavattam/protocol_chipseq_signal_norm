
# Shell standard
This document owns maintained Bash runtime, source form, static analysis, execution semantics, sourcing, wrapper and processing topology, submit bootstrap, parser diagnostics, and shell-test source form. `help.md` owns help semantics, `source_headers.md` owns literal headers and shebang profiles, `testing.md` owns general test operations, and `repository_layout.md` owns placement.

<br />

## Shell syntax and runtime (`SHELL.SYNTAX`)
**Classification:** `advisory` with deterministic parseability and declared-runtime evidence.

**Scope:** Maintained shell entrypoints, sourced libraries, installers, fixture recipes, and shell tests.

Use Bash for maintained shell except the explicit POSIX environment bootstrap. Bash code requires Bash 4.4 or newer and may use arrays, `[[ ]]`, arithmetic evaluation, `mapfile`, process substitution, and other supported Bash features. Parse Bash files with the resolved project Bash and parse `install/scripts/install_envs_entrypoint.sh` with `sh -n`.

The literal shebang row belongs to `source_headers.md`; this rule owns whether the chosen interpreter and syntax match runtime behavior. Do not claim POSIX compatibility for Bash source.

**Automation:** `tests/contract/repository/test_shell_syntax.sh` parses maintained shell and checks the Bash floor and POSIX exception. Registry coverage is `subset` because parsing does not prove runtime portability.

**Semantic remainder:** Review external shell-version behavior and constructs not exercised on all supported hosts.

**Exceptions:** Only the named POSIX bootstrap currently uses `/bin/sh` semantics.

<br />

## Safe-mode boundaries (`SHELL.RUNTIME.SAFE_MODE`)
**Classification:** `advisory` with deterministic declaration and exception portions.

**Scope:** Directly executed shell entrypoints, sourced libraries, and test orchestration.

Direct workflow and maintenance scripts normally use `set -euo pipefail` after runtime guards and before substantive work. Sourced libraries must not unilaterally change the caller's shell options. Handle expected nonzero statuses explicitly rather than disabling safety for a broad region.

The canonical test runner uses `set -u` and aggregates individual failures intentionally; it is the explicit runner exception. Individual shell tests use `set -euo pipefail` unless the shared harness requires a narrower, documented control-flow exception.

**Automation:** Source-form contracts inventory recognized safe-mode declarations, ordering, and named exceptions with `subset` coverage.

**Semantic remainder:** Review commands whose status is data, cleanup that must continue after failure, and sourced code that inherits caller options.

**Exceptions:** The canonical aggregate runner may omit `-e` and `pipefail` only while its tests prove complete failure aggregation and final nonzero status.

<br />

## Shell source form (`SHELL.SOURCE.FORM`)
**Classification:** `advisory` with deterministic simple-form portions.

**Scope:** Maintained Bash outside literal help and other heredoc bodies.

Use four-space code indentation and spaces rather than tabs. Put lifecycle and public functions before their invocation. Prefer `function name() {` for maintained Bash functions. Keep conditionals and continuations readable, quote expansions unless intentional splitting or globbing is documented, and use arrays for command construction rather than assembled command strings.

Maintained function and variable names use `snake_case`. Arrays use the established `arr_*` convention where it clarifies type. Controlled project abbreviations remain acceptable; avoid introducing opaque abbreviations whose meaning is not local and evident.

Indent `case` patterns one level inside `case`, commands one level inside the pattern, and terminators with their pattern. Use one empty line between adjacent nontrivial arms when it improves scanning; bounded case-arm spacing remains advisory.

Use `EOM` for ordinary shell-text heredocs and quote it when interpolation is not intended. Use a language-specific delimiter such as `PY` for embedded interpreter input. Another meaningful delimiter is permitted when required to avoid a collision or express literal semantics. Heredoc bodies follow their content owner: help bodies follow `help.md`, embedded Python follows Python source expectations, and literal fixtures remain unchanged. Do not reindent or reflow a heredoc without reviewing rendered behavior.

Diagnostics identify the failing script or function and direct errors to stderr. Use `printf` when format control matters; `echo` remains acceptable for simple fixed diagnostics. Do not print success before the operation succeeds.

**Automation:** `dev/audit/shell_source_form.py` checks only fixture-defined simple naming declarations, selected indentation forms, and recognized heredoc delimiters. Shell syntax and related audits provide additional evidence. Coverage is `subset`; no checker claims to parse arbitrary Bash layout, decide readability, or infer quoting intent.

**Semantic remainder:** Review readability, quoting intent, complex continuations, heredoc semantics, and whether command construction preserves argument boundaries.

**Exceptions:** Generated command fixtures and literal expected output preserve their owned representation.

<br />

## Source directives and roles (`SHELL.SOURCE.DIRECTIVES`)
**Classification:** `advisory` with deterministic recognized-policy portions.

**Scope:** Bash `source` directives, source-only libraries, dual-role helpers, and direct-execution guards.

Resolve repository libraries from a path anchored to the current source file, not the caller's working directory. Use the shared source helper when repeated sourcing must be idempotent. A source-only library returns to its caller and must not execute workflow work while loaded. A dual-role file separates definitions from an explicit direct-execution guard.

Place ShellCheck source directives immediately above the source operation they describe and keep their target accurate. Do not use a directive to conceal a dynamically unresolved path when the repository can compute a stable source path.

**Automation:** `dev/audit/source_policy.py` discovers current maintained shell paths and checks recognized source anchoring, direct guards, help sourcing, and role contracts with `subset` coverage.

**Semantic remainder:** Classify source-only versus dual-role intent and assess dynamic paths that static analysis cannot resolve.

**Exceptions:** Environment activation hooks may source an externally supplied path only when the public input contract validates it before use.

<br />

## Shell static analysis (`SHELL.STATIC_ANALYSIS`)
**Classification:** `advisory` with deterministic inventory and provenance evidence.

**Scope:** Maintained Bash entrypoints, sourced libraries, installers, fixture recipes, and shell tests discovered by the repository ShellCheck runner.

Run ShellCheck through `dev/audit/run_shellcheck.sh` with external-source following and the repository source path. The runner uses `${CONDA_PREFIX}/bin/shellcheck` directly, requires that path to be executable, and verifies version 0.10.0; it never selects ShellCheck through `PATH`. Check maintained Bash with `--shell=bash`, exclude `install/scripts/install_envs_entrypoint.sh` from that set, and check that POSIX bootstrap separately with `--shell=sh`. An ambient Homebrew or system executable is not repository evidence.

Keep statically resolvable `source=` directives exact. Use `SC1090` or `SC1091` suppressions only for genuinely dynamic or externally supplied sources. Every suppression has a narrow owner, rationale, and review trigger. A clean inventory is evidence, not automatic promotion to blocking enforcement.

Treat ShellCheck status 0 as a completed clean inventory and status 1 as a completed inventory with findings. Any other status is an infrastructure failure. Retain raw structured findings plus the executable, version, language-specific and total file counts, finding counts, raw statuses, and final status with each validation run.

**Automation:** `dev/audit/run_shellcheck.sh`, fixture-backed status and discovery contracts, and the registered command check explicit executable provenance in activated and non-activated caller contexts, Bash/POSIX splitting, invocation flags, structured evidence, and status classification with `subset` coverage. ShellCheck diagnostics remain upstream facets beneath this owner.

**Semantic remainder:** Review false positives, dynamic-source intent, suppression scope, portability, and whether findings should become blocking repository contracts.

**Exceptions:** A suppression does not authorize use of an executable outside `env_protocol` or conceal a statically resolvable repository source.

<br />

## Shell line length (`SHELL.LINE_LENGTH`)
**Classification:** `advisory` with deterministic bounded-pilot portions.

**Scope:** Maintained shell source outside exempt heredoc, URL, checksum, regex, parser-label, and fixture-literal contexts.

Keep ordinary shell source within 79 columns when a readable split preserves behavior. Prefer semantic command-array and condition boundaries. Do not split indivisible URLs, checksums, regular expressions, diagnostics, or literal interface text merely to satisfy a count.

**Automation:** `tests/contract/repository/test_shell_line_length.sh` reports repository-wide candidates and strict bounded pilots. Coverage is `subset`, not `exact`.

**Semantic remainder:** Decide whether a line is indivisible and whether a proposed split improves readability without changing quoting or execution.

**Exceptions:** Owned literal and indivisible forms remain documented checker exclusions.

<br />

## Unknown-option diagnostics (`SHELL.PARSER.UNKNOWN_OPTION`)
**Classification:** `deterministic` for recognized parser roles.

**Scope:** Maintained top-level shell parsers and local function parsers.

An unknown option emits the repository-standard error through the role-appropriate helper, preserves the exact offending token, emits applicable help, and returns or exits with status 1. Function diagnostics identify `${FUNCNAME[0]}`; script diagnostics identify the script basename. Do not silently shift or reinterpret an unknown option as a positional value.

**Automation:** `dev/audit/unknown_option_helpers.py` emits `SHELL.UNKNOWN_OPTION` summary evidence and maps its message, helper, function-name, help, and status facets to this owner. Focused fixtures provide `subset` coverage for recognized parser forms.

**Semantic remainder:** Identify forwarding and parser roles not represented by the recognized source forms.

**Exceptions:** A command that intentionally forwards an option must document the forwarding boundary and stop local parsing before that token.

<br />

## Wrapper topology (`SHELL.WRAPPER_TOPOLOGY`)
**Classification:** `advisory` with deterministic lifecycle and interface portions.

**Scope:** Maintained `execute_*.sh`, `submit_*.sh`, related user-facing wrappers, and their sourced libraries.

A wrapper owns user-interface orchestration: runtime guard, safe mode, source resolution, help and argument parsing, validation, environment handling, command construction, dispatch, and final status. Keep full help in one centralized owner and route help semantics to `help.md`.

Execute wrappers own user-facing orchestration and select serial, GNU Parallel, or Slurm dispatch. They construct the submit-worker command and may invoke it directly, through Parallel, or through `sbatch`. Submit wrappers are bounded worker entrypoints: they validate one worker payload, resolve task-specific selection such as `SLURM_ARRAY_TASK_ID`, and invoke the processing primitive without resubmitting work. Shared processing functions own reusable per-item computation and must not acquire wrapper parsing, environment activation, engine selection, or scheduler policy. Build commands as arrays and preserve dry-run output from the same argument vector used for execution.

Separate debug printing, state mutation, command construction, and dispatch once any block becomes more than a few lines. Debug helpers are read-only and return success when debugging is disabled. State-conversion helpers perform one transformation and leave validation separate. A helper that writes a global documents that generated global and its type through `help.md`.

Do not encode brittle repository-wide script counts. Discover maintained wrappers from owned paths and classify explicit compatibility delegates separately.

**Automation:** Wrapper-style, interface, and source-policy contracts check recognized lifecycle, centralized help, command arrays, dry-run parity, and directory contracts with `subset` coverage.

**Semantic remainder:** Review responsibility placement, environment lifetime, dispatch equivalence, and behavior-sensitive command changes.

**Exceptions:** A compatibility delegate may omit ordinary lifecycle stages only when it performs no independent work and has an explicit migration owner.

<br />

## Processing boundaries (`SHELL.PROCESSING.BOUNDARY`)
**Classification:** `advisory` with deterministic dependency-direction portions.

**Scope:** Execute, submit, and processing functions and the calls among them.

Execute wrappers dispatch submit-worker commands serially, through Parallel, or through Slurm. Submit wrappers validate and adapt one dispatched unit, then call the reusable processor; they do not submit a second scheduler job or duplicate the processor. Processing functions receive validated values, perform bounded workflow work, and return status; they do not parse top-level CLI options, emit top-level help, activate the project environment, select an execution engine, or submit scheduler jobs.

Keep local, Parallel, and Slurm paths behaviorally aligned at the processing boundary. Mode-specific setup that changes the processor's required inputs is explicit and tested.

**Automation:** Dependency and wrapper contracts inspect recognized calls and retired topology with `subset` coverage.

**Semantic remainder:** Decide whether setup belongs to orchestration or reusable processing and whether execution engines remain equivalent.

**Exceptions:** A scheduler-only packaging helper may live beside submit logic when it has no local processing responsibility.

<br />

## Submit bootstrap (`SHELL.SUBMIT.BOOTSTRAP`)
**Classification:** `advisory` with deterministic recognized-interface portions.

**Scope:** Maintained submit wrappers, compatibility delegates, and source helpers used before submission.

A submit wrapper retains its approved Bash declaration even when intentionally non-executable. Before scheduler or environment work, it establishes the Bash version contract, distinguishes sourced from direct use, resolves repository sources, and runs `main` only under the direct-execution guard. A worker-only wrapper rejects an unsupported Bash with top-level status 1 before newer syntax or side effects; functions return status rather than exiting the caller.

Direct submit interfaces must not source another submit wrapper as implementation. A compatibility submit interface delegates explicitly to its canonical owner, preserves arguments and status, and performs no competing orchestration. Guards and compatibility decisions precede work with side effects. When helper discovery genuinely depends on user-supplied `--dir_scr`, bootstrap parsing collects only that value, resolves the helper root, and delegates all ordinary parsing and orchestration to the owned lifecycle.

**Automation:** `dev/audit/source_policy.py --submit-bootstrap` and runtime interface contracts emit shebang, mode, source-contract, compatibility, guard-order, guard-contract, and direct-contract facets beneath this owner. Coverage is `subset`.

**Semantic remainder:** Review environment ordering, scheduler semantics, and whether a delegate remains justified.

**Exceptions:** Only registered compatibility delegates may omit independent `main` orchestration.

<br />

## Startup sourcing (`STARTUP.SOURCING`)
**Classification:** `advisory` with deterministic isolation-load evidence.

**Scope:** Maintained sourceable files below `lib/bash/`.

Every sourceable library loads in isolation after its declared prerequisites are available, defines behavior without starting a workflow, and leaves the caller in control. Load-time code may establish constants and definitions but must not parse CLI arguments, invoke production commands, or exit the caller.

**Automation:** `tests/contract/repository/test_startup_sources.sh` sources every maintained library in isolation and reports load failures or immediate exits. Coverage is `subset`.

**Semantic remainder:** Review global mutation, traps, option changes, expensive initialization, and side effects not detected by isolation loading.

**Exceptions:** None permit production workflow execution at source time.

<br />

## Shell-test source form (`SHELL.TEST.SOURCE_FORM`)
**Classification:** `advisory` with deterministic recognized-structure evidence.

**Scope:** Maintained shell tests below `tests/contract/` and `tests/integration/`.

Individual shell tests declare Bash, enable ordinary safe mode, set `TEST_NAME`, source `tests/support/test_helpers.sh`, report through the shared harness, and call `finish` when they reach normal completion. Tests use harness-provided scratch and log roots rather than repository-relative ad hoc paths.

For a multiline assertion or capture helper call, put the helper name on its own line and each helper positional argument on its own continuation line. Keep one short human-readable label as one quoted argument, attach explanatory comments to the call they introduce, and separate adjacent multiline helper-call blocks with one empty line except inside shell control-flow conditions.

`testing.md` owns taxonomy, gates, fixtures, evidence, cleanup, and result semantics. The aggregate runner is not an individual test and follows its explicit safe-mode and reporting exceptions.

**Automation:** The test-layout contract scans contract, local, Parallel, and ordinary Slurm shell tests and maps source-form findings to this owner. Coverage is `subset` for conditional and generated forms.

**Semantic remainder:** Review clarity, appropriate helper reuse, and intentional early-exit paths.

**Exceptions:** A test that cannot reach `finish` after a fatal prerequisite must record that disposition through the shared helper.
