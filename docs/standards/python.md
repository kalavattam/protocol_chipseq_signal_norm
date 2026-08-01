
# Python standard
This document owns maintained Python package, import, API, command-line, dependency, installation, source-layout realization, comments, naming, docstrings, Ruff, and Python-specific testing obligations. Shared source-layout semantics belong to [`source_layout.md`](source_layout.md), repository placement to [`repository_layout.md`](repository_layout.md), shared documentation semantics to [`help.md`](help.md), test operations to [`testing.md`](testing.md), and source-header profiles to [`source_headers.md`](source_headers.md).

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

Within a maintained module, use this topology when the corresponding parts exist: module docstring; future imports; standard-library, third-party, and local import groups; module constants and narrowly scoped initialization; public classes and functions; private helpers next to the public behavior they support; parser construction; `main()` orchestration; and the `__main__` guard last. Keep a small cohesive module intact. Extract a helper or domain module only when it creates a reusable responsibility boundary; topology is not authority to repartition working source mechanically.

**Automation:** `tests/contract/repository/test_python_boundaries.sh` rejects retired imports and recognized inverted dependencies with `subset` coverage.

**Semantic remainder:** Decide whether behavior is orchestration, reusable domain logic, or an intentional shared dispatcher. Also decide whether the behavior is sufficiently reusable to justify a module boundary.

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

## CLI help-literal layout (`PY.CLI.HELP.LAYOUT`)
**Classification:** `deterministic` for recognized static-literal source form.

**Scope:** Constant `help=` values built from one ordinary string literal or adjacent ordinary string literals inside maintained Python `parse_args()` functions.

[`HELP.AUDIENCE`](help.md) owns which help surface and level of detail serves the user. [`SOURCE.PROSE.WRAP`](source_layout.md) owns greedy source wrapping; this rule owns Python recognition and exact rendered-value preservation.

Keep a help value on one source line when it fits. When prose requires adjacent literals, wrap greedily at whole-word boundaries: fill each physical source line through the last complete word that fits within 79 columns, and break only when the next complete prose word and its required separator would exceed that limit. Preserve explicit paragraph escapes such as `\n\n`, literal newlines, indivisible tokens, and deliberate semantic boundaries.

Reflowing source literals must not change the rendered help value. Use one-line double-quoted literals under `PY.STRING.QUOTES`; do not use manual padding or trailing source whitespace to approach the limit.

**Automation:** `dev/audit/python_source_policy.py` checks physical width and premature word-boundary breaks for recognized static literals. `dev/tools/python_help_format.py` provides preview-by-default formatting for the proven adjacent-literal subset, measures the same Unicode-preserving one-line literal serialization that it emits, and rejects write operations outside it. Focused fixtures prove exact rendered values and idempotence. Coverage is `subset`.

**Semantic remainder:** Decide whether a token is indivisible, an explicit break is meaningful, or clearer help prose should replace mechanical reflow.

**Exceptions:** Dynamic help values, translated text, generated literals, and explicit non-prose formatting require an owned disposition before this source form is applied.

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

## Callable annotations (`PY.TYPE.ANNOTATIONS`)
**Classification:** `advisory` with deterministic presence portions.

**Scope:** Every maintained `def` and `async def` below `src/`, `tests/`, and `dev/`, including methods, tests, fixtures, nested helpers, callbacks, generators, and context-manager helpers.

Annotate every parameter and return value. Exclude only the actual first `self` or `cls` receiver of an instance or class method; a free-function parameter with the same spelling is not excluded. Annotate `*args` and `**kwargs` according to the values they collect, and annotate test functions with `-> None`. Do not require annotations on every local assignment.

Use an honest broad type when a narrow type would be false. Prefer a meaningful `Mapping`, `Sequence`, protocol, union, or `object` boundary over reflexive `Any`; use `Any` when the boundary is intentionally dynamic and record a non-obvious reason. `PY.TYPE.FLOOR` continues to own Python 3.11-compatible syntax.

**Automation:** `dev/audit/python_source_policy.py` checks parameter and return annotation presence across maintained Python. `dev/audit/python_source_evidence.py` inventories presence separately from correctness. No static type checker is selected or installed by this migration. Coverage is `subset`.

**Semantic remainder:** Review whether annotations are true, appropriately broad, useful to callers, and consistent with runtime behavior.

**Exceptions:** No callable role is exempt from presence. A deliberately dynamic boundary may use a broad annotation but not omit it.

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

## Exact docstring layout (`PY.DOCSTRING.LAYOUT`)
**Classification:** `deterministic` for recognized docstring source form.

**Scope:** Every docstring that exists in maintained Python below `src/`, `tests/`, and `dev/`.

Every docstring is one triple-double-quoted token; implicit concatenation is prohibited. Put the opening delimiter, summary, and closing delimiter on separate physical lines, including for a summary-only docstring. Align the opener exactly one indentation level inside a definition, and align the summary and closing delimiter exactly with the opener. Module docstrings begin at column zero. Use spaces for indentation, leave blank docstring rows empty, and do not dedent nonblank content before the opener's base indentation. Use lowercase `r"""` only when literal backslashes materially require raw-string semantics. Triple-single-quoted docstrings and `f`, `b`, `u`, or combined prefixes are prohibited.

Put no blank line between a definition signature and its docstring. Put exactly one blank line after a function or method docstring before the first executable statement and after a class docstring before the first member.

The summary is an imperative or descriptive sentence ending with a full stop. For a full docstring, put one blank line between the summary and the unheaded extended summary. Applicability and semantic content belong to `PY.DOCSTRING.NUMPY`.

The exact summary-only source form is:
```python
def load_config(path: Path) -> Config:
    """
    Load and validate one configuration file.
    """

    text = path.read_text(encoding="utf-8")
```

The source form follows [PEP 257](https://peps.python.org/pep-0257/) where compatible and intentionally adds the repository-required blank line after function and method docstrings.

**Automation:** `dev/audit/python_source_policy.py` checks the one-token form, delimiter selection, permitted prefix, opener, summary line and full stop, exact base indentation, empty blank rows, closer, signature adjacency, and the post-docstring blank line across maintained Python. Focused fixtures provide `subset` coverage. Ruff D200, D202, D203, and D212 conflict with the project form; D209, D211, D213, D300, and D204 can supply only partial evidence and are not selected by this migration.

**Semantic remainder:** Decide whether raw semantics are materially required and whether the summary is meaningful rather than merely punctuated.

**Exceptions:** A generated or literal fixture requires an owned applicability exclusion. Triple-single-quoted docstrings and prohibited prefixes have no stylistic exception.

<br />

## NumPy docstrings and applicability (`PY.DOCSTRING.NUMPY`)
**Classification:** `advisory` with deterministic coverage and section portions.

**Scope:** Maintained modules and applicable public, private, test, fixture, nested, abstract, overload, override, and protocol objects below `src/`, `tests/`, and `dev/`.

A summary-only docstring is permitted only when the annotated signature and summary communicate the complete contract and the object has no non-obvious side effects, failure behavior, units, formats, state invariants, ordering requirements, or caller obligations. A full docstring contains a short summary, an unheaded extended summary, and every applicable NumPy-style section in canonical order. Render every `Parameters` entry as `name : type`, with a type that agrees with the annotated and observed contract. Render every `Returns` or `Yields` entry as `name : type`; use meaningful result identities, including comma-separated identities for a tuple returned as one contractual unit. A multiline textual type follows [`SOURCE.DELIMITED.MULTILINE`](source_layout.md): indent inner rows once and align its closing-only delimiter with the `name : type` entry line. Do not add empty, name-repeating, or signature-repeating sections.

Every maintained module requires a useful docstring unless an intentionally empty package marker has an owned exception. Give public interfaces useful docstrings. Nontrivial private helpers and shared fixtures document assumptions, side effects, state, or error behavior; tiny private adapters and predicates may remain undocumented when their contract is evident. Test functions and classes may remain undocumented when their names and bodies fully expose intent. Overrides document material deviations; authoritative abstract and overload declarations document the public contract.

[`help.md`](help.md) owns shared section vocabulary and straight-single-quote prose delimiters, including the prohibition on Markdown backtick prose outside literal examples. `PY.DOCSTRING.LAYOUT` owns exact source form, and this rule owns applicability, recursive coverage, NumPy rendering, summary-only versus full selection, and usefulness. The [numpydoc style guide](https://numpydoc.readthedocs.io/en/latest/format.html) supplies the section-rendering basis.

**Automation:** `dev/audit/python_source_policy.py` deterministically checks canonical section names, exact underlines, canonical order, blank section boundaries, section-body indentation, malformed or incomplete boundaries, colon pseudo-headings, empty sections, complete annotated-parameter membership, typed `Parameters` entries, named and typed `Returns` and `Yields` entries, multiline textual-type closer alignment, agreement between documented and declared annotation spellings, agreement with stable direct returned or yielded tuple-component names, and overlength summary lines throughout maintained Python. `PY.DOCSTRING.NUMPY` is the sole diagnostic owner for the textual-type closer facet. The checker does not replace a meaningful scalar result identity with an incidental local name or guess an identity from a computed expression or branch-dependent transfer. `dev/audit/python_source_evidence.py` inventories object roles, source forms, recognized sections, and explicit applicability and content dispositions. Neither tool infers applicability, prose usefulness, or whether matching declared and documented types agree with runtime behavior. Coverage is `subset`.

**Semantic remainder:** Decide documentation applicability, useful content, summary-only versus full form, applicable sections, meaningful result identities, and whether types, side effects, failures, invariants, and examples are accurate.

**Exceptions:** Intentionally empty package markers, generated modules, tiny private helpers, and self-explanatory tests may be excluded only by documented file or role.

<br />

## Ordinary string quotes (`PY.STRING.QUOTES`)
**Classification:** `advisory` with deterministic recognized-literal portions.

**Scope:** Ordinary non-docstring string literals in maintained Python below `src/`, `tests/`, and `dev/`.

Use double quotes for an ordinary one-line string. A one-line single-quoted string is permitted when it materially avoids escaping embedded double quotes. Use triple double quotes for an ordinary multiline string; triple single quotes are prohibited by default and require a bounded literal-content exception when the value itself contains triple double quotes.

String prefixes remain available when runtime semantics require them. This rule does not prohibit raw strings, f-strings, or bytes literals in executable code; `PY.DOCSTRING.LAYOUT` owns the stricter docstring-prefix boundary.

Constructed prose literal wrapping and multiline triple-quoted content boundaries follow [`SOURCE.DELIMITED.MULTILINE`](source_layout.md). Put the first content character on the line after the opener and the last content character on the line before the closer. Use explicit `removeprefix("\n")` or `removesuffix("\n")` normalization when separating a delimiter must not change the runtime value.

**Automation:** `dev/audit/python_source_policy.py` checks recognized ordinary string quote selection and records embedded-quote and f-string-expression exceptions across maintained Python. Focused fixtures cover one-line, multiline, prefix, and literal-content boundaries. Coverage is `subset`.

**Semantic remainder:** Decide whether avoiding an escape is material, whether literal content requires the exception, and whether a prefix expresses required runtime semantics.

**Exceptions:** Literal fixtures and values containing the otherwise canonical delimiter may retain the narrowest quote form that preserves content clearly.

<br />

## Python source layout (`PY.SOURCE.LAYOUT`)
**Classification:** `advisory` with deterministic topology and candidate-evidence portions.

**Scope:** Maintained Python statements, expressions, calls, assignments, literals, comprehensions, parser builders, control flow, definitions, and result or transfer regions below `src/`, `tests/`, and `dev/`.

Realize [`SOURCE.LAYOUT.PARAGRAPHS`](source_layout.md) with four-space indentation and Python syntax. Keep related statements contiguous; put one blank line between distinct setup, acquisition, validation, transformation, side-effect, cleanup, and result phases. Put a visible blank-line boundary before and after each noncompact sibling `if`, `for`, `while`, `with`, `try`, or `match` statement so the end of one control-flow phase cannot be mistaken for the start of another. A direct single-transfer `if` guard keeps its transfer attached inside the guard, then uses one blank line before the following sibling statement. This following boundary also separates back-to-back compact guards. Semantic paragraphs may occur inside any suite without separating an `if` chain, `try` chain, decorator stack, block header, or other syntactically connected structure. These sibling, guard, phase, and syntactically connected-structure semantics apply in every Python suite.

Place module imports before executable module statements. Dependency availability checks do not create an exception: required project dependencies import normally, while a genuinely optional dependency uses an import architecture that does not leave later static imports below a call.

Keep related short repeated calls or assignments together by topic. Review long independent calls, assignments, literals, comprehensions, and parser-builder regions for both hyper-density and mechanical fragmentation. Use phase comments or helper extraction when blank lines alone do not expose the source organization.

Multiline Python calls and analogous supported delimited structures follow [`SOURCE.DELIMITED.MULTILINE`](source_layout.md), use one item per block-indented line, and include a trailing comma. The Python realization uses one item per block-indented line and a trailing comma in expanded calls and analogous supported structures. The analogous supported structures in the restored sentence are expanded forms; a bounded cohesive compact call row remains the narrow formatter-preserving exception and omits the trailing comma that would make the formatter expand it. PEP 8 supplies the upstream [indentation](https://peps.python.org/pep-0008/#indentation) and [trailing-comma](https://peps.python.org/pep-0008/#when-to-use-trailing-commas) basis.

Do not let a formatter-expanded membership collection displace the logical relationship in a compound condition. Extract truthful predicates, then compose the `if` or `while` condition from their names.

Apply compact-guard and transfer semantics to `return`, `raise`, `yield`, `yield from`, `break`, `continue`, and statically named process-exit calls such as `sys.exit()`. Keep a direct transfer compact inside its guard and separate the completed guard from its following sibling. Separate the transfer itself after substantive preparation or a distinct diagnostic, mutation, cleanup, control-flow, or result-building phase.

`X`, `Y`, and `Z` remain evidence concepts under [`SOURCE.LAYOUT.CANDIDATES`](source_layout.md). This standard establishes no numeric values before the representative Python pilot.

**Automation:** `dev/audit/python_source_policy.py` checks top-level definition and class-method blank-line topology, visible sibling boundaries before and after every noncompact control-flow statement, the following boundary after every standalone compact transfer guard, transitions in every suite nested within a callable involving local imports, explicit result-bearing action assignments, validation or verification after a substantive setup paragraph, completed action or validation followed by a new phase, assertions, and multiline data or case inventories immediately consumed by a loop, recognized transfers after an attached substantive phase, and static module imports that follow executable module statements across maintained Python. The import check deliberately supplements Ruff's E402 exemptions for calls such as `pytest.importorskip()`. The control-flow check treats one direct guarded action as compact and does not pretend to identify every semantic phase change. The transfer check recognizes Python transfer statements and statically named process-exit calls; it treats multiple attached statements, noncompact control flow, direct mutation or output, and multiline preparation as substantive. The semantic-phase check recognizes result assignments through an explicit vocabulary, keeps consecutive captures together, applies strict arrange–act separation in tests, requires substantive preparation before an initial action boundary in non-test helpers, recognizes validation and verification calls by static callable name, permits one directly attached mutation or observation that the call proves, recognizes direct `assert`, unittest-style assertion calls, `raises` or `warns` context managers, and recognizes a multiline assignment whose bound name is directly consumed by a following loop; X candidates, dynamically dispatched behavior, and same-syntax subject transitions remain review-owned. `dev/audit/python_source_evidence.py` produces X, two separately unit-bearing Y metrics, Z candidates, complete repeated-run members, every recognized suite transfer with preceding-phase membership, and separate syntactic blank-line-region density. Exact source fingerprints invalidate stale decisions but do not supply the semantic decision. The policy unit suite also discovers and checks every Python tool in the corrected-evidence successor so migration tooling cannot exempt itself through path placement. Ruff formatting remains separately governed by `PY.RUFF.FORMAT`.

**Semantic remainder:** Review paragraph boundaries not established by the deterministic subsets, repeated-construction topics, dynamically dispatched exit behavior, and whether wrapping, comments, intermediate names, or helper extraction best resolve a candidate.

**Exceptions:** Python-required connected syntax remains attached. Generated source and literal fixtures require owned applicability exclusions; other deviations use the governed exception process.

<br />

## Python comments (`PY.COMMENT.FORM`)
**Classification:** `advisory` with deterministic marker, spacing, and width portions.

**Scope:** Ordinary full-line, continuation, paragraph-separator, and inline comments in maintained Python below `src/`, `tests/`, and `dev/`.

[`SOURCE.COMMENT.ATTACHMENT`](source_layout.md) is the sole owner for attachment, prose, role, and usefulness. Python realizes marker and separator form by beginning each nonempty ordinary full-line or continuation comment with `# ` and using exactly `#` for an empty separator inside a multiline comment. Inline comments use exactly two spaces before `#` and one space after it.

Keep the complete physical source line, including indentation and marker, within 79 columns. Wrap at word boundaries and do not hyphenate an ordinary word merely to satisfy the limit. Rewrite or disposition an indivisible overlength token through `PY.FORMAT.LINE_LENGTH`.

Shebangs, source headers, encoding rows, tool directives, coverage pragmas, type-checker directives, generated markers, and literal fixtures retain their separately owned syntax. The superseded `#  ` and `#+ ` proposal is not a permitted ordinary-comment form.

These forms are consistent with [PEP 8 block and inline comment guidance](https://peps.python.org/pep-0008/#comments).

**Automation:** `dev/audit/python_source_policy.py` checks recognized ordinary markers, separator context, inline spacing, trailing whitespace, header and directive exclusions, 79-column width, and safely identifiable whole-paragraph capitalization, terminal punctuation, and sentence spacing across maintained Python. Wrapped prose is evaluated as one paragraph. Literal labels and statically recognizable code or operator fragments remain outside the deterministic prose subset; attachment and usefulness remain semantic review.

**Semantic remainder:** Review comment role, attachment, usefulness, prose, indivisible content, and whether the explanation belongs in a docstring or maintained documentation. Apply that review through the shared owner.

**Exceptions:** Separately owned directives, headers, generated content, and literal fixtures are applicability exclusions rather than ordinary-comment exceptions.

<br />

## Python identifiers (`PY.NAMING.IDENTIFIERS`)
**Classification:** `advisory` with deterministic AST and casing opportunities.

**Scope:** Maintained Python packages, modules, functions, methods, parameters, locals, attributes, properties, types, constants, tests, files, visibility markers, and exact external boundaries.

Apply [`SOURCE.NAMING.SEMANTICS`](source_layout.md) using Python-native forms:
- packages and modules use lowercase names, with underscores where needed for readability;
- functions, methods, parameters, locals, attributes, and properties use `snake_case`;
- classes, exceptions, protocols, enums, and other user-defined types use `UpperCamelCase`;
- module- or class-level semantic constants and constant-like enum members use `SCREAMING_SNAKE_CASE`;
- one leading underscore marks a non-public interface;
- one trailing underscore resolves a keyword collision; and
- interpreter-owned documented `__dunder__` names are used only for their defined protocols.

Do not use `lowerCamelCase` for new project-defined names. Preserve it only at an externally owned interface, override, generated schema, or serialized field that requires exact spelling, and adapt it at the boundary when practical. Do not use `UpperCamelCase` as general emphasis or promote mutable module state to constant spelling.

Type variables may use conventional short capitals or concise `UpperCamelCase`. Two leading underscores without trailing underscores require a demonstrated class name-mangling need. Controlled scientific, file-format, and project abbreviations remain subject to the shared abbreviation and evidence rules.

This owner governs five vocabulary statuses:
- `ordinary_short_word` is an unambiguous ordinary word and is not abbreviation evidence;
- `allowed_domain_term` is an accepted scientific, format, ecosystem, or project term in its recorded matching context;
- `protected_external_or_interface_spelling` must remain exact only at the recorded external or public boundary and in a direct adapter;
- `prohibited_internal_shorthand` is rejected in new or changed internal names and records a descriptive replacement; and
- `review_candidate` remains evidence for contextual review, not an automatic defect.

The reviewed Python vocabulary permits unambiguous domain and ecosystem terms such as `BAM`, `CRAM`, `IP`, `SIQ`, `bedGraph`, `CSV`, `TSV`, `JSON`, `gzip`, `CLI`, `SSH`, `cwd`, and the conventional imported NumPy alias `np`. Existing names that mirror protected CLI destinations may retain project spellings such as `fil_in`, `fil_out`, `siz_bin`, `skp_pfx`, `scl_fct`, `frg`, `psc`, and `qntl` only at that interface or in code directly adapting it. New implementation locals use descriptive words. In particular, do not introduce private shorthand segments such as `cfg`, `cmb`, `col`, `cvg`, `dp`, `ext`, `fh`, `fmt`, `py`, `rc`, `sb`, `src`, `str`, `xs`, or `ys`; use names such as `configuration`, `combined`, `column`, `coverage`, `decimal_places`, `extension`, `handle`, `format`, `python`, `return_code`, `start_bins`, `source`, `text`, `values`, or another truthful domain name.

[`dev/config/python_naming_vocabulary.toml`](../../dev/config/python_naming_vocabulary.toml) is the sole maintained entry list. Every entry records its spelling, exactly one approved status, matching kind, one or more approved contexts, explicit evidence-candidate membership, and a replacement only when prohibited internal shorthand requires one. `dev/audit/python_naming_vocabulary.py` validates the list and supplies exact projections to the checker and evidence producer; neither derived surface creates or promotes an entry. Adding or revising an entry must pass the anti-accretion procedure in [`GOV.CHANGE.GOLDEN_FIRST`](governance.md).

No role-specific too-short or too-long threshold is normative before the Python naming inventory and representative pilot. Approved internal renames are direct atomic changes without compatibility shims; production wrapper and production Python command names remain protected by [`SOURCE.NAMING.SEMANTICS`](source_layout.md).

The casing and visibility baseline follows [PEP 8 naming conventions](https://peps.python.org/pep-0008/#naming-conventions).

**Automation:** `dev/audit/python_source_policy.py` checks module, callable, parameter, class, local, and stored-attribute casing or prohibited opaque-shorthand forms across maintained Python. `dev/audit/python_source_evidence.py` produces role-aware length, abbreviation-candidate, maintained-Python reference-surface, and rename-status evidence. It does not prove grammar, external references, or rename safety.

**Semantic remainder:** Review grammatical role, abbreviation clarity, visibility, external ownership, candidate length, dynamic references, and whether a proposed rename improves meaning.

**Exceptions:** Exact external contracts and documented language protocols may retain required spelling at a bounded interface. Other deviations use the governed exception process.

<br />

## Naming-length evidence (`PY.NAMING.LENGTH`)
**Classification:** `advisory` evidence only.

**Scope:** Role-aware Python identifier length, segmentation, abbreviation candidates, scope use, maintained-Python reference surfaces, and pilot rename proposals.

Name length identifies review candidates rather than violations. Compare adjacent thresholds by role and inspect short, long, public, private, test, local, attribute, type, abbreviation, and repeated-context examples before selecting a value. Keep casing and visibility under `PY.NAMING.IDENTIFIERS` and grammatical meaning and migration under [`SOURCE.NAMING.SEMANTICS`](source_layout.md).

A direct-rename status is evidence, not authorization. Before an approved rename, reconcile definitions, imports, calls, tests, fixtures, documentation, Shell or configuration references, serialization, reflection, generated commands, and external ownership. Pilot proposals remain unresolved until human semantic review.

**Automation:** `dev/audit/python_source_evidence.py` reports length histograms, adjacent threshold counts, abbreviation candidates, scope use, maintained-Python path references, string mentions, and conservative rename statuses. Coverage is `subset`; non-Python and dynamic surfaces remain unresolved.

**Semantic remainder:** Select thresholds, judge abbreviations and grammar, inspect every reference surface, and decide whether a rename improves meaning.

**Exceptions:** Conventional short names and clear long names may be retained with their role-specific rationale; no length alone requires a rename.

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

Repository audit tools form the `dev.audit` package. Imports use package-qualified names, and registry, documentation, fixture, and contract invocations launch tools from the repository root with `python -m dev.audit.<module>`. Direct path execution is unsupported unless a tool documents a narrow compatibility entry point.

**Automation:** `ruff check --no-fix src tests dev` is blocking and provides `subset` coverage for the approved core. Safe-fix idempotence is part of rollout evidence; unsafe fixes remain prohibited.

**Semantic remainder:** Review any future selection change and the separately deferred E501, docstring, and formatter phases.

**Exceptions:** No maintained path has an E402 exception. No exception allows configuration to silently supersede an owner rule.

<br />

## Ruff selection (`PY.RUFF.SELECTION`)
**Classification:** `semantic-only`.

**Scope:** Enabled Ruff rule families and individual codes.

Enable a family only after findings are inspected and mapped to an owner. Ruff's `E*`, `F*`, `I*`, `UP*`, `B*`, `C4*`, `SIM*`, `RUF*`, and `D*` codes remain upstream diagnostic facets, not repository rule IDs.

**Automation:** Inventories count findings by upstream code. The selected core is enforced only through `PY.RUFF.POLICY`; deferred codes remain nonblocking evidence.

**Semantic remainder:** Review future selection changes. E501 and all D rules remain deferred. D200, D202, D203, and D212 conflict with the approved docstring form; any nonconflicting D-rule subset requires a later owner, registry, fixture, and pilot decision.

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
