
# Python standard
This document owns maintained Python package, import, API, command-line, dependency, installation, docstring, Ruff, and Python-specific testing obligations. Repository placement belongs to `repository_layout.md`, shared documentation semantics to `help.md`, test operations to `testing.md`, and source-header profiles to `source_headers.md`.

<br />

## Package discovery and names (`PY.PACKAGE.SRC`, `PY.PACKAGE.NAMES`)
**Classification:** `deterministic` for current metadata and discovery, with advisory extensibility review.

**Scope:** `pyproject.toml`, the installed distribution, and maintained imports below `src/`, `tests/`, and `dev/`.

Use setuptools `src` discovery. Maintained import packages live below `src/`; the current import package is `protocol_chipseq_signal_norm` and the distribution is `protocol-chipseq-signal-norm`. Preserve the intentional underscore and hyphen distinction in metadata, installation evidence, and user documentation.

Tests and production code import the installed package. Repository-relative aliases and the retired top-level `scripts` package are prohibited. A later domain subpackage may be added below the maintained package without changing this rule; a new top-level import package requires a governed package-boundary decision.

**Automation:** Package and layout contracts check current discovery, import names, retired paths, and installed metadata with `subset` coverage.

**Semantic remainder:** Review whether a new package is a subdomain of the maintained API or a new governed distribution boundary.

**Exceptions:** None preserve retired import aliases or top-level compatibility packages.

<br />

## Module layering (`PY.MODULE.LAYERS`)
**Classification:** `advisory` with deterministic import-boundary portions.

**Scope:** Maintained modules below `src/protocol_chipseq_signal_norm/`.

Command modules live in `protocol_chipseq_signal_norm.cli`. Reusable methods, classes, parsers, formatting, and I/O behavior live in `utilities` or a later governed domain package. Utilities do not import CLI modules. CLI modules may import utilities; peer CLI imports require a documented shared-dispatch exception. Move shared behavior downward rather than create circular imports.

**Automation:** `tests/contract/repository/test_python_boundaries.sh` rejects retired imports and recognized inverted dependencies with `subset` coverage.

**Semantic remainder:** Decide whether behavior is orchestration, reusable domain logic, or an intentional shared dispatcher.

**Exceptions:** A peer CLI dispatcher must be documented, acyclic, and covered by direct tests.

<br />

## Import safety (`PY.IMPORT.SIDE_EFFECTS`)
**Classification:** `advisory` with deterministic runtime probes.

**Scope:** Importing maintained package modules in a clean process.

Importing a module must not parse arguments, perform workflow I/O, start subprocesses, mutate process-global state, or emit user output. Constants and definitions may initialize at import time when deterministic and inexpensive.

**Automation:** Import and package contracts probe selected modules for import success and captured output with `subset` coverage.

**Semantic remainder:** Review expensive initialization, environment reads, registration side effects, and mutations not visible to bounded probes.

**Exceptions:** None permit workflow execution or argument parsing during import.

<br />

## Public APIs (`PY.API.PUBLIC`)
**Classification:** `semantic-only` with deterministic naming evidence.

**Scope:** Maintained package modules, documented utilities, re-exports, and callers.

Public APIs are deliberate. A leading underscore marks private implementation. Define `__all__` only when a module or package intentionally publishes a bounded API; do not add it mechanically or use re-exports to preserve retired paths. Migrate callers and tests atomically when a documented public utility changes.

**Automation:** Documentation and import audits inventory names and re-exports but do not decide public usefulness or compatibility.

**Semantic remainder:** Classify public versus private behavior and assess compatibility impact.

**Exceptions:** A temporary compatibility export requires an approved migration record and removal condition.

<br />

## CLI registration (`PY.CLI.ENTRYPOINT`)
**Classification:** `deterministic`.

**Scope:** Maintained Python CLIs, `[project.scripts]`, and explicit module execution.

Register every maintained CLI in `[project.scripts]` using an underscore-named command. The console command and `python -m protocol_chipseq_signal_norm.cli.<module>` invoke the same implementation and must not fork behavior.

**Automation:** Package and CLI contracts compare entry points, modules, command names, selected help and status behavior, and console/module parity with `subset` coverage.

**Semantic remainder:** `None` for registered parity.

**Exceptions:** Internal modules that are not public commands remain unregistered.

<br />

## Typed CLI main (`PY.CLI.MAIN`)
**Classification:** `deterministic`.

**Scope:** Every maintained CLI module.

Provide this contract:
```python
def main(argv: list[str] | None = None) -> int:
    ...


if __name__ == "__main__":
    raise SystemExit(main())
```
`main()` owns orchestration and returns an explicit integer, including zero on success. Parser construction and reusable domain work remain in testable helpers.

**Automation:** Python CLI contracts parse maintained modules and execute selected status paths with `subset` coverage.

**Semantic remainder:** Review whether orchestration has leaked into parser or reusable helpers.

**Exceptions:** None for maintained public CLIs.

<br />

## CLI parsers (`PY.CLI.PARSER`)
**Classification:** `advisory` with deterministic registered-interface portions.

**Scope:** Parser construction, positional names, optional arguments, aliases, help, and destinations for maintained CLIs.

Use `CapArgumentParser` and `add_help_cap()`. Every visible optional argument has help text, an explicit canonical `dest`, and canonical aliases before compatibility aliases. Positional names use canonical `snake_case`. Hidden hyphenated aliases may map to canonical underscore options, but retired semantic names are not aliases. Preserve registered canonical aliases and `--dp`.

**Automation:** Parser, alias, and help contracts check registered facts and recognized source forms with `subset` coverage.

**Semantic remainder:** Decide public visibility, alias compatibility, and whether a parser abstraction remains appropriate.

**Exceptions:** A hidden compatibility alias requires an active migration reason and must remain absent from public help.

<br />

## Error and status ownership (`PY.ERROR.EXIT`)
**Classification:** `advisory` with deterministic selected-path evidence.

**Scope:** Anticipated CLI failures, reusable exceptions, diagnostics, and exit status.

The CLI layer translates anticipated user or input failures into diagnostics and stable nonzero statuses. Reusable functions raise specific exceptions instead of exiting. Do not catch exceptions merely to discard context; unexpected programmer errors remain visible during tests.

**Automation:** CLI tests cover selected parsing, validation, diagnostic, and status paths with `subset` coverage.

**Semantic remainder:** Classify anticipated user errors versus programmer defects and choose stable diagnostic context.

**Exceptions:** None permit a reusable utility to terminate the process for an ordinary validation failure.

<br />

## Dependency categories (`PY.PACKAGE.DEPENDENCIES`)
**Classification:** `advisory` with deterministic metadata-parity portions.

**Scope:** Package runtime dependencies, build requirements, environment-management tools, and external workflow programs.

Declare import-time and runtime Python dependencies in `[project.dependencies]` and build requirements in `[build-system]`. Keep environment-management tools in the environment YAML. Aligners, samtools, GNU Parallel, Slurm, R, and other workflow executables are not Python package dependencies; their owners document and validate them.

**Automation:** Package/environment parity contracts compare registered dependency categories with `subset` coverage.

**Semantic remainder:** Decide whether a dependency is imported by the installed package, required only to build, or invoked by a workflow.

**Exceptions:** An optional Python feature requires explicit extra or plugin ownership before its dependency is omitted from the core set.

<br />

## Package metadata (`PY.PACKAGE.METADATA`)
**Classification:** `deterministic` for recorded metadata with semantic release decisions.

**Scope:** Project name, version, license, authorship, README, and package metadata.

Keep one canonical metadata value for each field and require the installed distribution to agree with it. `pyproject.toml`, rather than this standard, owns the current version data. Changing that value is an explicit release decision.

**Automation:** Package contracts inspect `pyproject.toml`, built metadata, and installed metadata with `subset` coverage.

**Semantic remainder:** Current parity is deterministic. Approval of a future release version or descriptive metadata change is a prospective decision and does not change present conformance.

**Exceptions:** None permit divergent source and installed metadata.

<br />

## Supported Python floor (`PY.VERSION.FLOOR`)
**Classification:** `deterministic` for recorded and executed version evidence.

**Scope:** Package metadata, environment selection, Ruff target, executed interpreter, maintained source syntax, local guards, and documentation.

Maintained Python requires `Python >= 3.11`. `requires-python`, the authoritative environment, Ruff `target-version`, documentation, local guards, and executed validation interpreter must express that same floor. Source syntax across maintained `src/`, `tests/`, and `dev/` must parse at the floor.

**Automation:** `dev/audit/python_version_policy.py` emits `PY.VERSION.FLOOR` with structured facets for interpreter, documentation, environment, metadata, Ruff target, guard, and source evidence. It runs through `env_protocol` and parses maintained source with the Python 3.11 grammar. Coverage is `subset` because grammar acceptance cannot prove runtime behavior on every supported interpreter.

**Semantic remainder:** Review dependency-level interpreter constraints and syntax not exercised by the configured floor interpreter.

**Exceptions:** A tool requiring a newer interpreter must run in a separately owned environment and remain outside the maintained floor scope.

<br />

## Floor-compatible types and syntax (`PY.TYPE.FLOOR`)
**Classification:** `advisory` with deterministic compilation evidence.

**Scope:** Public and reusable function boundaries, CLI signatures, annotations, and maintained source syntax.

Type public and reusable boundaries and every CLI `main(argv) -> int`. Use annotations and syntax valid at Python 3.11. Compilation is evidence for syntax compatibility, not proof of type correctness.

**Automation:** External-scratch compilation and focused tests provide `subset` evidence. The retired `PYTHON.STARTUP_COMPILE` adapter is registered under this owner rather than as a separate rule.

**Semantic remainder:** Review annotation usefulness, runtime typing effects, and compatibility not covered by compilation.

**Exceptions:** Dynamic boundaries may use broader types when a narrower contract would be false; explain the accepted shape.

<br />

## Editable installation (`PY.INSTALL.EDITABLE`)
**Classification:** `advisory` with deterministic managed-installation evidence.

**Scope:** Environment creation, reuse, managed YAML update, dry runs, and editable package installation.

Development uses a managed editable install into `env_protocol` with `--no-deps --no-build-isolation`. Creation, reuse, and managed update refresh the install. Update reconciles the environment YAML with installed packages frozen, without pruning or opportunistically changing unrelated packages. Repeatable `--update_package` selections are the safer deliberate path for a reviewed incremental update such as this migration's ShellCheck installation. Preserve the selected manager—Mamba-driven invocation uses `mamba run`, Conda-driven invocation uses `conda run`—and never install into an unrelated environment. Dry runs display environment and package actions but perform neither. Capturing and retaining a before/after package transaction is a validation responsibility; the installer does not automatically preserve that report.

**Automation:** Installation and package contracts exercise command construction, reuse, update, dry-run noninstallation, and installed metadata with `subset` coverage.

**Semantic remainder:** Review manager provenance and environment ownership.

**Exceptions:** None permit an implicit install into the caller's active unrelated environment.

<br />

## Checkout-independent runtime (`PY.RUNTIME.CWD`, `PY.RESOURCE.ACCESS`)
**Classification:** `advisory` with deterministic registered-command evidence.

**Scope:** Installed imports, console commands, package resources, explicit runtime inputs, and the Slurm bundle.

Installed imports and console commands work outside the checkout. Resolve package resources through package-aware APIs and runtime inputs through explicit paths; do not derive resources from the caller's current directory. The Slurm bundle is the provenance exception: its launcher prepends bundled `src/` and uses the required remote Python so another checkout's editable install cannot win.

**Automation:** Package and Slurm bundle contracts execute registered commands from external directories and inspect bundled provenance with `subset` coverage.

**Semantic remainder:** Identify resources that must be packaged and review provenance boundaries.

**Exceptions:** The documented Slurm bundle behavior is the only current checkout-isolation exception.

<br />

## Python docstrings (`PY.DOCSTRING.NUMPY`)
**Classification:** `advisory` with deterministic coverage and syntax portions.

**Scope:** Public modules, functions, classes, parser builders, CLI `main()` functions, and nontrivial private helpers.

Give public objects concise useful docstrings. Use NumPy/SciPy-style sections when they materially clarify parameters, returns, yields, or exceptions. Do not add empty sections or docstrings that repeat only a name. Tiny private helpers may remain undocumented; nontrivial private helpers should state assumptions or side effects. `help.md` owns shared section vocabulary; this rule owns Python applicability, recursive coverage, NumPy semantics, and usefulness.

**Automation:** Repository doc-coverage checks provide `subset` evidence. All Ruff D rules are deferred to the later docstring phase; the reserved convention is D211/D212 with D203/D213 ignored. Neither syntax checking nor coverage proves usefulness.

**Semantic remainder:** Decide public applicability, useful content, and whether a section materially clarifies the contract.

**Exceptions:** Generated modules and tiny private helpers may be excluded by documented file or role.

<br />

## Python line length (`PY.FORMAT.LINE_LENGTH`)
**Classification:** `advisory` with deterministic configured evidence.

**Scope:** Maintained Python source below `src/`, `tests/`, and `dev/`.

Keep ordinary source within 79 columns where Ruff can do so without degrading URLs, command strings, diagnostics, regexes, or other indivisible content. Markdown prose and shell heredocs have separate owners.

**Automation:** An explicit Ruff `E501` inventory may report candidates but remains nonblocking and unapproved for enforcement.

**Semantic remainder:** Review whether an overlength construct is indivisible and clearer than an available rewrite.

**Exceptions:** Stable URLs and other indivisible literals may remain overlength with an owned disposition.

<br />

## Ruff governance (`PY.RUFF.POLICY`)
**Classification:** `deterministic` for the approved lint core.

**Scope:** Ruff configuration, inventories, rollout decisions, and maintained Python roots.

`python.md` is normative and `pyproject.toml` is Ruff's sole machine configuration. Target Python 3.11 and list `src`, `tests`, and `dev` plus generated/private exclusions explicitly. The approved 7j lint core selects the explicit E4, E7, E9, F, I001, UP, B, C4, SIM, and RUF codes recorded in `pyproject.toml`; it excludes UP009, E501, and every D rule. Run lint and format checks separately. Configuration changes do not promote a proposed rule without standards review.

**Automation:** `ruff check --no-fix src tests dev` is blocking and provides `subset` coverage for the approved core. Safe-fix idempotence is part of rollout evidence; unsafe fixes remain prohibited.

**Semantic remainder:** Review any future selection change and the separately deferred E501, docstring, and formatter phases.

**Exceptions:** E402 is ignored only for `dev/**/*.py` and `tests/**/*.py`, whose executable-tool and fixture roles may establish repository imports after setup. No exception allows configuration to silently supersede an owner rule.

<br />

## Ruff selection (`PY.RUFF.SELECTION`)
**Classification:** `semantic-only`.

**Scope:** Enabled Ruff rule families and individual codes.

Enable a family only after findings are inspected and mapped to an owner. Ruff's `E*`, `F*`, `I*`, `UP*`, `B*`, `C4*`, `SIM*`, `RUF*`, and `D*` codes remain upstream diagnostic facets, not repository rule IDs.

**Automation:** Inventories count findings by upstream code. The selected core is enforced only through `PY.RUFF.POLICY`; deferred codes remain nonblocking evidence.

**Semantic remainder:** Review future selection changes. E501 and all D rules remain deferred; the later docstring phase uses D211/D212 and retains the mutually exclusive D203/D213 ignores.

**Exceptions:** UP009 is excluded because the source-header standard requires the encoding row. No deferred rule is implicitly promoted.

<br />

## Ruff boundary (`PY.RUFF.BOUNDARY`)
**Classification:** `advisory`.

**Scope:** Division between upstream Ruff behavior and repository-owned checkers.

Ruff owns configured generic syntax, imports, modernization, bug patterns, simplifications, and consistency. Repository checkers retain source headers, CLI aliases, package layers, NumPy-section semantics, and other behavior Ruff cannot prove. When Ruff conflicts with a repository owner, the owner wins until policy is explicitly changed.

**Automation:** Registry reconciliation maps Ruff coverage to owner sections and reports overreach.

**Semantic remainder:** Decide whether Ruff has semantic parity with a repository contract.

**Exceptions:** UP009 is not selected and must not remove the required source-header encoding row.

<br />

## Ruff suppressions (`PY.RUFF.SUPPRESSIONS`)
**Classification:** `advisory` with deterministic record portions.

**Scope:** Inline `noqa`, per-file ignores, and lasting Ruff exclusions.

Prefer correcting code. Inline suppression is limited to a genuine local exception. Per-file ignores require a documented file-role distinction. Every lasting suppression has an owner, rationale, and review condition; unexplained blanket ignores are prohibited.

**Automation:** Ruff and registry inventories collect suppressions and configured ignores with `subset` coverage.

**Semantic remainder:** Assess whether the finding is false, the exception is narrow, and the review trigger is adequate.

**Exceptions:** Only an owned suppression record authorizes deviation.

<br />

## Ruff formatting (`PY.RUFF.FORMAT`)
**Classification:** `advisory` with deterministic inventory and idempotence portions.

**Scope:** `ruff format --check`, formatter approval, and any future formatting write.

Keep formatting separate from lint. Do not run the formatter in write mode until formatter policy and representative changes are approved. A formatting inventory proves only difference from configured output; it does not prove approval, readability, or semantic safety. Any approved write must be idempotent.

**Automation:** `ruff format --check` has `subset` nonblocking coverage. Repository-wide formatter adoption remains deferred by the 7j approval.

**Semantic remainder:** Approve formatter adoption and review representative changes.

**Exceptions:** None permit repository-wide formatting from an inventory-only configuration.

<br />

## Ruff upgrades (`PY.RUFF.UPGRADE`)
**Classification:** `advisory` with deterministic version evidence.

**Scope:** Ruff pin changes in `env_protocol.yml` and resulting behavior.

Pin Ruff. An upgrade requires before-and-after findings, standards and configuration review, and passing focused plus repository contracts.

**Automation:** Environment and inventory reports record the installed and configured version with `subset` coverage.

**Semantic remainder:** Review rule changes, changed fixes, formatter output, and rollout impact.

**Exceptions:** Security remediation may accelerate the upgrade only with recorded residual comparison work.

<br />

## Python-specific tests (`PY.TEST.LAYERS`)
**Classification:** `advisory` with deterministic registered contracts.

**Scope:** Python unit, CLI, package, and workflow-integration obligations.

Unit-test pure parsing, calculation, formatting, and stream transformations. CLI tests cover parsing, help, diagnostics, status ownership, and console/module parity. Package tests cover discovery, build, editable installation, imports outside the checkout, registered commands, reuse refresh, and dry-run noninstallation. Workflow integration covers shell-to-console dispatch and relevant external boundaries. `testing.md` owns the repository taxonomy and runner.

**Automation:** Focused Python unit and contract groups provide `subset` coverage.

**Semantic remainder:** Decide representative values, failure cases, and proportional integration depth.

**Exceptions:** None replace required package or CLI boundary evidence with unit-only tests.

<br />

## Python-generated state (`PY.ARTIFACTS.LOCATION`)
**Classification:** `deterministic` for path cleanliness.

**Scope:** Placement of Python bytecode, caches, wheels, build trees, installation reports, and Ruff inventories.

Keep generated Python state outside maintained production roots, normally below the test artifact root owned by `TEST.EVIDENCE`. Maintained roots remain free of `__pycache__`, `.pyc`, wheels, build trees, and generated reports after validation. `TEST.EVIDENCE` owns retention decisions.

**Automation:** Cleanliness contracts and final repository inventories check known generated paths with `subset` coverage.

**Semantic remainder:** `None` for path cleanliness.

**Exceptions:** None permit caches or build products in maintained source roots. A maintained release artifact follows `TEST.EVIDENCE` and its release owner.
