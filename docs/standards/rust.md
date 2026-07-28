
# Rust standard
This document owns Rust-specific source syntax, naming, documentation markers, stable-language boundaries, layout realization, and future toolchain decisions. The requirements become applicable when maintained Rust source is introduced. [`source_layout.md`](source_layout.md) owns shared semantics, [`help.md`](help.md) owns shared documentation concepts, [`repository_layout.md`](repository_layout.md) owns placement, and [`testing.md`](testing.md) owns repository test operations.

<br />

## Applicability and package boundary (`RUST.PACKAGE.APPLICABILITY`)
**Classification:** `semantic-only`.

**Scope:** Any future maintained Rust crate, package, binary, library, workspace, build script, example, benchmark, or test source.

Before maintained Rust source is added, approve its repository root, crate and package role, dependency direction, public interface, build outputs, runtime integration, and proportional proof. Use Cargo's conventional layout unless repository evidence justifies an explicitly governed alternative. A workspace, multiple package roots, build script, or generated bindings require separate ownership because they expand execution and artifact boundaries.

The official [Cargo package-layout guidance](https://doc.rust-lang.org/cargo/guide/project-layout.html) is the primary basis for the future package decision.

**Automation:** No Cargo metadata, Rust source, package contract, fixture, registry entry, or maintained Rust artifact currently implements this rule.

**Semantic remainder:** Classify the crate or package role, choose its public boundary, and determine whether one package or a workspace is justified.

**Exceptions:** Vendored, generated, or externally copied Rust is not maintained project source unless an approved owner adopts it and records provenance and applicability.

<br />

## Rust identifiers (`RUST.NAMING.IDENTIFIERS`)
**Classification:** `advisory` with deterministic compiler-lint opportunities.

**Scope:** Project-defined Rust crates, modules, functions, methods, locals, fields, lifetimes, types, variants, constants, statics, macros, tests, files, and external boundaries.

Apply [`SOURCE.NAMING.SEMANTICS`](source_layout.md) using Rust-native casing:
- modules, functions, methods, local variables, fields, and lifetimes use `snake_case`;
- structs, enums, traits, type aliases, and enum variants use `UpperCamelCase`;
- constants and statics use `SCREAMING_SNAKE_CASE`;
- macros ordinarily use `snake_case`; and
- package and crate names follow Cargo and Rust conventions rather than Python package rules.

Adapt externally owned names at serialization, FFI, protocol, or generated-code boundaries when the interface permits it. Do not corrupt an external contract merely to satisfy project casing, and do not propagate an external spelling through internal code when a bounded adapter can contain it.

Rust's standard nonstandard-style lint group includes casing checks for types, names, and globals; see the [rustc lint documentation](https://doc.rust-lang.org/rustc/lints/listing/warn-by-default.html#nonstandard-style). Compiler availability does not by itself select lint levels or prove grammatical clarity.

**Automation:** No dedicated checker or registry entry exists for Rust naming. No compiler invocation, naming inventory, proposal ledger, fixture, or threshold exists; candidate thresholds remain unset until maintained Rust symbols provide role-aware evidence.

**Semantic remainder:** Review grammatical role, abbreviation clarity, external ownership, generic and lifetime scope, and too-short or too-long candidates beyond compiler casing.

**Exceptions:** Language-defined conventional names, generated bindings, and exact external interfaces may retain required spelling through owned applicability exclusions.

<br />

## Rust comments and documentation (`RUST.COMMENT.DOCUMENTATION`)
**Classification:** `advisory` with deterministic marker opportunities.

**Scope:** Ordinary Rust comments, outer item documentation, inner module or crate documentation, attributes that carry documentation, and generated rustdoc output.

Ordinary full-line, continuation, and inline comments use `// `; use exactly `//` as an empty separator inside a multiline ordinary comment. Exact spacing before an inline marker remains formatter-owned once a Rust toolchain is selected. Follow [`SOURCE.COMMENT.ATTACHMENT`](source_layout.md). Use `///` for outer documentation attached to the following item and `//!` for inner documentation attached to the containing module or crate. Do not use documentation markers or attributes as phase-comment decoration.

Document public and non-obvious interfaces with applicable purpose, parameters or fields, returned values, errors, panics, safety requirements, side effects, runtime requirements, examples, and references. Empty boilerplate and prose that repeats only an item name or signature are prohibited. Shared concepts remain with [`HELP.SECTION.SCHEMA`](help.md); Rust source renders them through rustdoc-native Markdown and conventions.

The syntax and attachment meaning of `//`, `///`, and `//!` follow the [Rust Reference](https://doc.rust-lang.org/reference/comments.html). Future content and rendered-output policy should be assessed against the [rustdoc book](https://doc.rust-lang.org/rustdoc/how-to-write-documentation.html).

**Automation:** No dedicated checker or registry entry exists for Rust documentation. No rustdoc invocation, coverage inventory, or fixture exists.

**Semantic remainder:** Decide documentation applicability, useful content, item versus module attachment, and which error, panic, safety, and example obligations materially affect callers.

**Exceptions:** Generated bindings and external source require owned applicability exclusions. Exact literal examples retain the syntax they demonstrate.

<br />

## Rust source layout (`RUST.SOURCE.LAYOUT`)
**Classification:** `semantic-only`.

**Scope:** Calls, builders, literals, iterator chains, macro invocations, `match` expressions, guards, loops, tail expressions, `?`, explicit transfers, and item organization in maintained Rust source.

Realize [`SOURCE.LAYOUT.PARAGRAPHS`](source_layout.md) in agreement with the official Rust style and selected `rustfmt` output. Keep related construction and transformation steps together; separate setup, acquisition, validation, transformation, side effects, cleanup, and result production when they form distinct phases.

Put a blank-line boundary before and after each noncompact sibling `if`, `for`, `while`, `loop`, or `match` control-flow unit when `rustfmt` preserves that boundary. Keep `else` and connected match-arm syntax attached to the construct they complete, and retain a direct single-transfer guard as a compact exception.

Review calls, builders, struct and collection literals, iterator chains, macro invocations, and `match` arms for both hyper-density and mechanical fragmentation. A formatter-owned wrap does not decide whether a builder or iterator chain should be named, split, or extracted. Macro token trees and `match` arms require construct-aware evidence rather than raw punctuation counting.

Multiline calls use block indentation, one argument per line, an aligned closing delimiter, and a trailing comma:
```rust
let rows = load_table(
    args.tbl_met,
    skp_pfx,
    args.verbose,
    other_value,
);
```

Treat tail expressions as Rust-native result forms. Do not add an explicit `return` merely to imitate Python or R. Apply semantic transfer separation to explicit `return`, `break`, and `continue`, to diverging forms such as `panic!`, and to a terminal tail expression after substantive preparation. The `?` operator propagates control on the residual path but often remains attached to the expression it qualifies; separate the surrounding phase when acquisition, validation, cleanup, or result construction changes.

The multiline form follows the official [Rust Style Guide for expressions](https://doc.rust-lang.org/style-guide/expressions.html) and [items](https://doc.rust-lang.org/style-guide/items.html). The `?` behavior follows the [Rust Reference](https://doc.rust-lang.org/reference/expressions/operator-expr.html#the-question-mark-operator).

**Automation:** No selected `rustfmt` version or configuration, parser-aware checker, density producer, candidate threshold, fixture, or registry entry exists. `X`, `Y`, and `Z` values remain unset pending representative Rust source and a focused pilot.

**Semantic remainder:** Review paragraph boundaries, chain and builder responsibilities, macro and `match` structure, tail-result clarity, `?` attachment, guard compactness, and whether wrapping or extraction is clearer.

**Exceptions:** Formatter output may refine syntax only after the tool relationship, version, configuration, exclusions, and representative evidence are recorded. Generated macro output is outside maintained source unless explicitly adopted.

<br />

## Stable-language boundary (`RUST.LANGUAGE.STABILITY`)
**Classification:** `semantic-only`.

**Scope:** Rust language features, compiler channels, feature gates, standard-library APIs, examples, and normative source claims.

Maintained Rust policy is written against stable Rust unless a separately approved toolchain boundary states otherwise. Do not present an unstable or nightly-only feature as available stable syntax. Any proposal that mentions an unstable feature must label it explicitly, identify its feature gate and channel, isolate it from stable requirements, and state why stable alternatives are insufficient.

Stable Rust currently has no general project rule for generator-style `yield`; do not import Python's `yield` layout into Rust. If a future stable release changes the relevant language boundary, update this owner from the then-current official documentation and repository compiler floor before using the feature.

The Rust release-channel and feature-gate model is described in [The Rust Programming Language](https://doc.rust-lang.org/book/appendix-07-nightly-rust.html).

**Automation:** No Rust compiler floor, channel check, feature-gate inventory, fixture, or registry entry exists.

**Semantic remainder:** Determine the supported stable compiler floor, identify unstable features, and evaluate whether any nightly-only boundary is justified and containable.

**Exceptions:** A nightly-only tool or experiment requires a separately owned environment, explicit feature gates, isolated paths, and proof that public build and runtime remain complete without it unless the standards owner approves a broader contract.

<br />

## Rust toolchain promotion (`RUST.TOOLCHAIN.PROMOTION`)
**Classification:** `semantic-only`.

**Scope:** Formatter, compiler, lint, documentation, test, package-layout, fixture, registry, and migration decisions for future maintained Rust source.

Do not claim implementation merely because Rust ships integrated tools. Before promotion, decide and record:
- stable compiler floor, Cargo package or workspace layout, dependency policy, lockfile policy, and generated-output boundary;
- `rustfmt` version, configuration, check and write modes, and representative output;
- compiler lint levels and whether [Clippy](https://doc.rust-lang.org/clippy/configuration.html) adds mapped, nonduplicative value;
- rustdoc applicability, warnings, links, examples, and rendered-output checks;
- unit, integration, documentation, and interface tests using the native test framework described in [The Rust Programming Language](https://doc.rust-lang.org/book/ch11-01-writing-tests.html); and
- source-layout candidate evidence for calls, builders, literals, chains, macros, `match`, tail expressions, `?`, and explicit transfers.

Promotion requires pinned or otherwise reproducibly selected tools, invocation options, positive and negative fixtures for repository-added behavior, current-source findings, representative formatter previews, threshold evidence, applicability exclusions, registry coverage metadata, focused contracts, idempotence, and a bounded migration plan. Upstream defaults do not silently create repository owners.

**Automation:** No Rust tool is selected, installed, configured, registered, or claimed as conformance evidence by this standard.

**Semantic remainder:** Choose the smallest stable toolchain and proof surface that support the maintained Rust role while preserving semantic review.

**Exceptions:** An exploratory compiler or tool invocation is not repository conformance and must not mutate maintained source or establish policy by precedent.
