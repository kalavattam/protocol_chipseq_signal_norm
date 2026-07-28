
# Testing standard
This document owns repository test taxonomy, runner behavior, gates, fixtures, generated evidence, shared reporting, cleanup safety, and proportional proof. Python- and shell-specific source obligations remain in their language standards; detailed wet Slurm operation remains in [`tests/integration/slurm/README.md`](../../tests/integration/slurm/README.md).

<br />

## Test layers (`TEST.LAYERS`)
**Classification:** `advisory` with deterministic placement evidence.

**Scope:** Maintained tests below `tests/unit/`, `tests/contract/`, and `tests/integration/`, plus shared test-only code and fixtures.

Classify a test by what it proves:
- `unit` isolates pure parsing, calculation, formatting, or state transitions without invoking a production workflow.
- `contract` proves repository structure, a command interface, or a policy-to-implementation connection.
- `integration/local` executes a bounded fixture-backed workflow on the local host.
- `integration/parallel` executes a bounded GNU Parallel path and remains optional.
- `integration/slurm` contains ordinary scheduler integration and the separately coordinated wet workflow.

Tests, assertions, fakes, coordinators, and harnesses are test-only dependencies. Production code must not import from `tests/`, depend on generated evidence, or require fixture roots at runtime. [`repository_layout.md`](repository_layout.md) owns directory placement and dependency direction.

**Automation:** Repository-layout contracts provide deterministic placement and dependency evidence with `subset` coverage. They do not decide whether a test's assertions prove the claimed layer.

**Semantic remainder:** Review the behavior proved, external boundaries crossed, and whether a lower-cost layer would provide equivalent confidence.

**Exceptions:** A deliberately cross-layer contract must state why its combined boundary is necessary and keep external effects behind the applicable gate.

<br />

## Safe runner (`TEST.RUNNER.SAFE`)
**Classification:** `deterministic`.

**Scope:** Group selection and discovery performed by `tests/run_tests.sh`.

With no group, the canonical runner selects `all-safe`: unit, contract, local integration, and Parallel integration. Parallel tests remain selected but self-skip unless their gate is enabled. Explicit groups may narrow the run.

Resolve every selected test to one canonical path and execute it at most once, even when group arguments overlap or repeat. `--list` prints the deduplicated selection without creating artifacts, preparing fixtures, or executing tests. Validate every requested group before discovery, gate evaluation, fixture preparation, artifact creation, or test execution. Unknown groups fail without performing any of those actions.

The safe runner must never discover or execute `tests/integration/slurm/run_wet_tests.sh`. It may select the ordinary Slurm integration only through the explicit `integration-slurm` group; that test must still enforce its own gate before scheduler work.

**Automation:** A runner contract exercises default, explicit, repeated, overlapping, and invalid group selections; verifies unique paths; and proves categorical wet-runner exclusion. Coverage is `subset` because the contract does not prove every future discovery implementation or execution failure mode.

**Semantic remainder:** `None` for discovery and deduplication.

**Exceptions:** None permit the wet coordinator to enter safe-runner discovery.

<br />

## Optional gates (`TEST.GATES.OPTIONAL`)
**Classification:** `deterministic`.

**Scope:** Ordinary optional test gates, including `RUN_ATRIA`, `RUN_DOWNLOAD`, `RUN_PARALLEL`, `RUN_SLURM`, and ordinary non-wet `WAIT_SLURM` handling.

Ordinary Boolean gates accept `true`, `t`, `yes`, `y`, and `1` as enabled and `false`, `f`, `no`, `n`, and `0` as disabled, case-insensitively. Unset or empty values are disabled. Surrounding whitespace and every other nonempty value are invalid and fail before the gated action.

Do not prepare fixtures solely for a disabled optional group. Every gated test checks its own gate before workflow-specific dependency discovery, fixture preparation, network access, parallel execution, or scheduler submission. Shared harness bootstrap may resolve the maintained test interpreter and canonical Boolean normalizer before the test checks its gate; it must not probe or prepare the gated workflow.

**Automation:** `dev/audit/boolean_contracts.py` and focused positive and negative fixtures inventory ordinary gates and normalized semantics with `subset` coverage for recognized source forms.

**Semantic remainder:** `None` for normalized values; deciding whether an integration should be optional remains test design.

**Exceptions:** Wet Slurm confirmation follows `TEST.SLURM.WET`, not this normalization contract.

<br />

## Wet Slurm workflow (`TEST.SLURM.WET`)
**Classification:** `advisory` with deterministic confirmation and discovery evidence.

**Scope:** `tests/integration/slurm/run_wet_tests.sh`, its coordinator, and commands that can submit or wait for the checksummed two-job wet workflow.

Wet execution is separate from the safe runner and requires exact values for all three gates before any scheduler action:
```text
RUN_SLURM=1 WAIT_SLURM=1 CONFIRM_SLURM_WET=1
```
True-like variants are insufficient. The coordinator must verify the gates, inputs, destination, and run identity before submission. Detailed preparation, monitoring, collection, and recovery procedures remain in the operational README.

**Automation:** Slurm coordinator units and source contracts cover gate order, exact values, and safe-runner exclusion. They do not perform a wet submission in the safe suite.

**Semantic remainder:** Confirm scheduler availability, account and partition selection, input provenance, expected cost, and authorization for each wet run.

**Exceptions:** None bypass the three exact gates or permit safe-runner discovery.

<br />

## Fixtures (`TEST.FIXTURES`)
**Classification:** `advisory` with deterministic inventory and ignore-policy portions.

**Scope:** Tracked static fixtures, generated workflow fixtures, fixture recipes, provenance documentation, and fixture consumers.

Every registered generated workflow fixture root contains a tracked `README.md` and `make.sh`. The source-derived fixture inventory, rather than a durable fixed count, owns the current set. The README records provenance, generation prerequisites, deterministic expectations, and ignored outputs. Generated data remains ignored; recipes and documentation remain tracked. Generated workflow roots are distinct from documented tracked static fixtures, including Slurm fixtures and the text-only `python_source_policy` checker inputs, which are not cleanup targets.

Retries must not hide nondeterministic failures. Control time, randomness, ordering, and concurrency inputs where practical. An intentionally observational test records its instability boundary and preserves the first failure rather than converting a later retry into unconditional success.

Tests may generate a missing fixture only after the applicable layer and gate are enabled. Production commands must not depend on test fixtures. A fixture refresh is a separately reviewable change and must not silently replace generated outputs during standards validation.

**Automation:** Fixture and ignore contracts inventory required recipes, documentation, sentinels, tracked-file boundaries, and production independence with `subset` coverage.

**Semantic remainder:** Review provenance, biological or domain representativeness, deterministic sufficiency, and whether a fixture remains minimal.

**Exceptions:** Upstream immutable inputs may use a documented checksum-based recipe when deterministic local generation is impossible.

<br />

## Generated evidence (`TEST.EVIDENCE`)
**Classification:** `advisory` with deterministic path and isolation portions.

**Scope:** Test logs, temporary products, reports, caches, bytecode, build products, and inventory output.

Generated test evidence is artifact-rooted. The exact validated `TEST_ARTIFACT_ROOT` is its owner: use the canonical ignored `artifacts/tests/` default or an approved absolute external scratch root that is nonempty, outside the repository, and not `/`. Arbitrary in-repository overrides are invalid. Derive `logs/`, `tmp/`, pytest cache, and Python bytecode beneath that exact root.

Evidence from the current run may support review, but production behavior must not read it as authority. Generated evidence is reproducible or explicitly labeled observational, ignored when repository-local, and disposable. Historical checkpoint evidence is immutable input to later assessment and must not be overwritten by validation.

**Automation:** The shared harness validates the canonical default or an approved external root before creating directories. Repository contracts check exact cache derivation, rejection of arbitrary repository roots, ignored paths, production independence, and the absence of caches in maintained roots with `subset` coverage.

**Semantic remainder:** Decide whether evidence is sufficient, reproducible, and appropriately retained for the decision it supports.

**Exceptions:** A maintained release artifact requires its own owner and is not test evidence merely because a test produced it.

<br />

## Shared harness reporting (`TEST.HARNESS.REPORTING`)
**Classification:** `deterministic`.

**Scope:** Maintained shell tests that source `tests/support/test_helpers.sh`.

Set `TEST_NAME`, source the shared helper, and use `print_section`, `record_pass`, `record_fail`, `record_warn`, `record_skip`, and `finish` for repository-standard reporting. `record_fail` increments the failure count; warnings and skips remain visible but nonfatal. `finish` prints the complete summary and returns nonzero exactly when failures were recorded.

Tests must not replace shared counters or emit a second competing final summary. A test may exit early only through a shared skip or fatal-precondition path that records the disposition.

**Automation:** The harness unit and source-form contract exercise counter transitions, output, exit status, required calls, and applicable exclusions with `subset` coverage for recognized shell tests.

**Semantic remainder:** `None` for reporting mechanics.

**Exceptions:** Pure Python unit tests use pytest reporting. The canonical aggregate runner owns its own group summary and is not an individual harness test.

<br />

## Scoped cleanup (`TEST.CLEANUP.SCOPED`)
**Classification:** `deterministic`.

**Scope:** `tests/support/clean_artifacts.sh` and any operation that deletes generated fixture or test-output paths.

Cleanup requires an explicit target class and supports a dry run. Before deletion, resolve the exact repository root, reject a mismatched checkout, verify that output roots contain no tracked paths, and verify that generated fixture roots contain only their permitted tracked README and recipe. Delete ignored content only through explicitly enumerated roots.

Cleanup is destructive maintenance, not ordinary validation. Never infer authority to clean from a test request, a failed run, or the presence of ignored files.

**Automation:** Focused cleanup contracts use temporary repositories to prove target selection, dry-run behavior, tracked-file refusal, root verification, and ignored-only deletion with `subset` coverage. The contract does not authorize or exhaustively model every future cleanup implementation.

**Semantic remainder:** `None` after targets and authorization are explicit.

**Exceptions:** None permit broad `git clean`, an unresolved path, or deletion of tracked or historical evidence.

<br />

## Proportional proof (`TEST.PROOF.PROPORTIONAL`)
**Classification:** `advisory`.

**Scope:** Validation planning and handoff for repository changes.

Run focused checks for changed owners and implementations before broader suites. Run write-free checks before scratch-backed tests, safe local tests before optional external integrations, and source-form checks before behavior that depends on that source. Proof must cover the changed boundary and its likely regressions without enabling unrelated network, Parallel, scheduler, or wet work.

Record unavailable proof, skipped optional gates, warnings, generated-output locations, and residual semantic review. A clean checker does not replace the semantic remainder owned by a rule.

**Automation:** No checker can decide proportionality. Test selection inventories and coverage reports provide bounded evidence for review.

**Semantic remainder:** Decide whether the selected evidence is sufficient for the risk and scope of the change.

**Exceptions:** Emergency validation may be narrower only when the omitted proof, risk, owner, and required follow-up are recorded.
