
# Repository layout standard
This document owns repository-root responsibilities, cross-root dependency direction, and the placement boundary between maintained source, support material, generated state, and private-only work. Domain standards own behavior within those boundaries and must not define competing locations.

<br />

## Production architecture (`LAYOUT.PRODUCTION.ROOTS`)
**Classification:** `deterministic`.

**Scope:** Required production and governance directories, prohibited retired directories, and the discovered Python package set.

Maintained production code and its authoritative standards use these required roots:
```text
bin/                                         maintained shell entrypoints
lib/bash/core/                               validation, bootstrap, state, output
lib/bash/dispatch/                           local, Parallel, Slurm, and Python dispatch
lib/bash/workflows/                          domain and workflow processing
lib/bash/help/                               executable shell help owners
src/protocol_chipseq_signal_norm/            sole distribution import tree
src/protocol_chipseq_signal_norm/cli/         Python command-line modules
src/protocol_chipseq_signal_norm/utilities/   reusable Python utilities
docs/standards/                               normative repository standards
```

The following retired roots are prohibited:
```text
scripts/
analysis/
scraps/
protocol_chipseq_signal_norm
lib/bash/blog/
docs/style/
```

**Automation:** `tests/contract/repository/test_repository_layout.sh` checks the current enumerated directory presence and absence and discovered package set with `subset` coverage. It does not provide a complete accepted/rejected boundary-fixture matrix for every future root shape. Content responsibilities and dependency direction are governed by the sections below.

**Semantic remainder:** `None` for the enumerated path-shape checks.

**Exceptions:** Only an approved staged migration under `LAYOUT.MIGRATION.NO_SHIMS` may temporarily alter the required or prohibited set.

<br />

## Shell ownership (`LAYOUT.SHELL.BOUNDARY`)
**Classification:** `advisory`.

**Scope:** Maintained shell source under `bin/` and `lib/bash/`.

`bin/` contains maintained user-facing shell entrypoints and approved delegates that preserve a supported command interface. Sourced implementation files do not live in `bin/`. `lib/bash/core/` owns generic validation, bootstrap, state, and output helpers; `dispatch/` owns execution backends; `workflows/` owns domain processing; and `help/` owns centralized executable help functions.

Cross-directory sourcing must resolve from the entrypoint or source file location, never from the caller's current working directory. A basename-only helper lookup is permitted only when the loader proves the basename is unique across its declared responsibility roots.

Supported interface delegates stay in `bin/`; they must not recreate a retired source tree. Detailed shell organization and sourcing behavior belong to [`shell.md`](shell.md), and help ownership belongs to [`help.md`](help.md).

**Automation:** No execution-registry entry is dedicated to this rule ID. The repository-layout, shell-syntax, wrapper, source-policy, help, and interface contracts provide related evidence without establishing registered coverage for this rule.

**Semantic remainder:** Review determines whether a helper's responsibility belongs in `core/`, `dispatch/`, or `workflows/`, and whether a compatibility delegate still represents a supported interface.

**Exceptions:** Use the exception process in [`governance.md`](governance.md); an exception may not establish an alternate shell architecture.

<br />

## Python ownership (`LAYOUT.PYTHON.SRC`)
**Classification:** `advisory`.

**Scope:** Maintained importable Python under `src/protocol_chipseq_signal_norm/`.

`src/` is the package-discovery root, and `src/protocol_chipseq_signal_norm/` is the only distribution import tree. Command-line modules live in `cli/`; importable reusable logic lives in `utilities/` or a later explicitly governed domain package. Tests and developer tools import the installed package name rather than a repository-directory alias.

Repository scripts, test helpers, and developer tooling do not become import namespaces merely because they contain Python. Detailed module roles, import direction, entry-point registration, and package behavior belong to [`python.md`](python.md).

**Automation:** No execution-registry entry is dedicated to this rule ID. The repository-layout, Python-package, and Python-boundary contracts provide related evidence without establishing registered coverage for this rule.

**Semantic remainder:** Review determines whether new reusable logic belongs in an existing package or justifies a new governed domain package.

**Exceptions:** A new package directory requires an approved standards change plus synchronized package discovery, registry, contract, and documentation updates.

<br />

## Installation ownership (`LAYOUT.INSTALL.BOUNDARY`)
**Classification:** `advisory`.

**Scope:** Maintained installation support under `install/`.

`install/scripts/` contains bootstrap and installation-specific entrypoints. `install/envs/` contains versioned environment definitions consumed by those entrypoints. `install/README.md` documents the supported installation surface. Workflow execution belongs in `bin/` and `lib/bash/`, not in `install/`.

Installation code may consume repository packaging metadata and environment definitions, but normal installed runtime behavior must not depend on the checkout's `install/` tree. Generated installer logs, package builds, and environment evidence belong under `artifacts/`.

**Automation:** No execution-registry entry is dedicated to this rule ID. Installation-layout, shell-syntax, help, and Python-package contracts provide related evidence without establishing registered coverage for this rule.

**Semantic remainder:** Review distinguishes installation lifecycle behavior from ordinary workflow execution and reusable runtime logic.

**Exceptions:** A bootstrap with interpreter constraints may use an explicitly documented profile, but it remains under `install/scripts/` and must retain its dedicated contract.

<br />

## Tests, development tools, and generated evidence (`LAYOUT.SUPPORT.BOUNDARY`)
**Classification:** `advisory`.

**Scope:** `tests/`, `dev/`, and `artifacts/`.

Tests follow the semantic taxonomy in [`testing.md`](testing.md): unit tests in `tests/unit/`, repository and interface contracts in `tests/contract/`, environment-facing tests in `tests/integration/`, reusable test-only helpers in `tests/support/`, and fixture recipes and documentation in `tests/fixtures/`. A fixture root tracks its `make.sh` and `README.md`; its generated inputs stay ignored.

Repository audit implementations live in `dev/audit/`, explicitly invoked mutation tools in `dev/tools/`, machine registry and configuration in `dev/config/`, and schemas in `dev/schemas/`. Developer tooling may inspect production source, but production source must not import or source `dev/` or `tests/`.

All repository-generated test, build, installation, audit, checkpoint, and scheduler evidence lives under the ignored `artifacts/` boundary documented in [`artifacts/README.md`](../../artifacts/README.md). Maintained tests may inspect artifacts created during the same run; production behavior must not require pre-existing artifacts. Tool-native caches outside `artifacts/` must remain ignored, disposable, and non-authoritative.

Production roots must not contain test coordinators, fake dependencies, generated fixtures, bytecode, caches, wheels, build metadata, audit reports, or preserved run output.

**Automation:** No execution-registry entry is dedicated to this rule ID. Test-layout, package, registry, and generated-artifact contracts provide related evidence without establishing registered coverage for this rule.

**Semantic remainder:** Review distinguishes reusable production behavior from test coordination, audit logic, and run evidence.

**Exceptions:** None for fixtures. A fixture is generated by its recipe and ignored; committing one is not an approvable exception, because a tracked fixture forfeits the discovery isolation the ignore rule provides and reintroduces the maintained-input burden the recipe exists to remove. An upstream immutable input that cannot be generated locally remains [`TEST.FIXTURES`](testing.md#fixtures-testfixtures)' documented checksum-recipe case, not a committed fixture.

<br />

## Documentation and executable narratives (`LAYOUT.DOCS.BOUNDARY`)
**Classification:** `advisory`.

**Scope:** Top-level and subsystem `README.md` files, `docs/`, and `notebooks/`.

The top-level `README.md` owns the public project entry point. A subsystem `README.md` documents the maintained surface beside its owner. `docs/standards/` contains normative standards; `docs/design/` contains proposals, decisions, and architecture rationale that remain non-normative unless a maintained standard explicitly adopts them.

`notebooks/workflows/` contains reproducible runbooks, and `notebooks/validation/` contains validation narratives and recorded interpretation, as indexed by [`notebooks/README.md`](../../notebooks/README.md). Automated assertions, test coordinators, and generated notebook products do not live in documentation roots.

Tracked static documentation assets may live under a documented `docs/` asset boundary. Generated documents, rendered notebooks, plots, and transient exports belong under `artifacts/` or another explicitly documented ignored output boundary.

**Automation:** No execution-registry entry is dedicated to this rule ID. Documentation-consistency, Markdown, link, and test-layout checks provide related evidence without establishing registered coverage for this rule.

**Semantic remainder:** Review determines whether material is normative, explanatory, a reproducible runbook, a validation narrative, or generated evidence.

**Exceptions:** Historical evidence may retain obsolete path text when it is clearly identified as immutable history and cannot be mistaken for current instructions.

<br />

## Data boundaries (`LAYOUT.DATA.BOUNDARY`)
**Classification:** `advisory`.

**Scope:** Repository-local project inputs and results under `data/`.

`data/` contains project data, reference inputs, and associated manifests; it does not contain executable production source, package resources, test fixtures, or developer tools. Small tracked inputs must have documented provenance and a maintained use. Large, local, downloaded, symlinked, or generated data remains ignored under a documented data boundary.

Reusable test data belongs in `tests/fixtures/`, not `data/`. Repository-generated test and audit evidence belongs in `artifacts/`. A file required at installed runtime is either an explicit runtime input or an intentionally packaged resource under `src/protocol_chipseq_signal_norm/`; code must not infer such a resource from a developer's local `data/` tree.

**Automation:** No execution-registry entry is dedicated to this rule ID. Ignore-policy and repository-layout checks may provide related placement and tracked-state evidence; data meaning and provenance require review.

**Semantic remainder:** Review determines whether a data product is a maintained source input, a test fixture, generated evidence, or private analysis output.

**Exceptions:** A tracked generated data product requires explicit ownership, provenance, regeneration instructions, and a reason it cannot remain in an ignored output boundary.

<br />

## Cross-root dependency direction (`LAYOUT.DEPENDENCY.DIRECTION`)
**Classification:** `advisory`.

**Scope:** Dependencies among maintained repository roots.

Allowed dependencies point toward production implementation and configuration:
```text
bin/ --------------------------> lib/bash/
bin/ and lib/bash/ ------------> installed Python entry points or modules
install/ ----------------------> packaging metadata and install/envs/
tests/ and dev/ ---------------> production roots
docs/ and notebooks/ ----------> documented public interfaces
```

Reverse dependencies are prohibited. In particular, `src/`, `bin/`, and `lib/bash/` must not depend on `tests/`, `dev/`, `notebooks/`, `artifacts/`, ignored local data, or the sibling private repository. Documentation and notebooks may demonstrate public interfaces but must not supply required runtime implementation. Tests and developer tools must not become undeclared production dependencies.

**Automation:** No execution-registry entry is dedicated to this rule ID. Python-boundary, package, shell-source, dependency-closure, and scheduler-bundle contracts provide related evidence without establishing registered coverage for this rule.

**Semantic remainder:** Review is required for dynamic shell dispatch, optional developer-only integrations, and dependencies that cross more than one responsibility boundary.

**Exceptions:** An optional audit may accept the sibling private repository as explicit evidence input, but public build, installation, tests, and runtime remain complete without it, and tracked public output must not disclose private content.

<br />

## Public and private material (`LAYOUT.PRIVATE.FORBIDDEN`)
**Classification:** `advisory`.

**Scope:** Public-repository analysis, scratch, blog, and local working material.

Maintained public analysis belongs in an appropriate public documentation, notebook, test, or data boundary according to its role. Private analysis, long-lived scratch material, private blog drafts, and troubleshooting records belong in the sibling private repository.

An ignored `tmp/` directory may hold short-lived local work or preserved migration evidence, but it is not maintained source and must not be imported, sourced, packaged, tested as an authoritative input, or referenced by current user instructions. Material that becomes durable must move to an owned public boundary or the private repository.

**Automation:** No execution-registry entry is dedicated to this rule ID. The repository-layout contract rejects prohibited public roots on behalf of `LAYOUT.PRODUCTION.ROOTS`; ignore and dependency audits provide related evidence without establishing registered coverage for this rule.

**Semantic remainder:** Review determines whether material is durable public documentation, private work, or disposable local state.

**Exceptions:** None may make private-only material or an ignored local path a public runtime dependency.

<br />

## Layout changes and retired paths (`LAYOUT.MIGRATION.NO_SHIMS`)
**Classification:** `advisory`.

**Scope:** Moves, new roots, new packages, compatibility delegates, and references to retired paths.

A layout change updates the normative owner, execution registry when applicable, packaging metadata, source and import references, help, current documentation, tests, contracts, and ignore rules in one reconciled change or an approved staged migration. A new maintained top-level root must have a documented owner and must be added to the relevant standards, indexes, and contracts before source is placed there.

Do not preserve retired architecture through compatibility directories, filesystem symlinks, duplicate wrappers, Python re-export modules, import aliases, or obsolete module paths. A supported user-interface delegate is permitted only in its canonical current root, with explicit interface ownership and contract coverage; it must not recreate the retired source location.

Historical immutable evidence is not rewritten to imitate current paths. Current manifests, standards, tests, help, runbooks, and live documentation must use the current architecture.

**Automation:** No execution-registry entry is dedicated to this rule ID. Repository-layout, Python-package, Python-boundary, source-policy, documentation-consistency, and interface contracts provide related retired-path evidence without establishing registered coverage for this rule.

**Semantic remainder:** Review distinguishes a supported interface delegate from an architectural shim and determines whether historical text is immutable evidence or current instruction.

**Exceptions:** A staged migration must satisfy [`governance.md`](governance.md), identify its temporary paths and compatibility behavior, and define an objective removal condition.
