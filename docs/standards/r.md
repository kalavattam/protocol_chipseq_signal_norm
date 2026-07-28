
# R standard
This document owns R-specific source syntax, naming, documentation markers, layout realization, applicability, and future toolchain decisions. The requirements become applicable when maintained R source is introduced. [`source_layout.md`](source_layout.md) owns shared semantics, [`help.md`](help.md) owns shared documentation concepts, [`repository_layout.md`](repository_layout.md) owns placement, and [`testing.md`](testing.md) owns repository test operations.

<br />

## Applicability and project form (`R.PROJECT.APPLICABILITY`)
**Classification:** `semantic-only`.

**Scope:** Proposed standalone R scripts, reusable R modules, and any future R package introduced as maintained repository source.

Before maintained R source is added, classify it as a standalone script or package-owned source and update the repository layout, dependency direction, runtime documentation, and proportional proof for that role. A standalone script follows the source-layout, naming, comment, and runtime requirements that apply without package metadata. Package-only obligations such as `DESCRIPTION`, `NAMESPACE`, package documentation, exported API review, and package checks apply only after a package boundary is approved.

Do not introduce package infrastructure merely to host one analysis script, and do not let a growing reusable R surface remain an undeclared collection of scripts when package ownership would materially improve installation, documentation, testing, or dependency control.

The [Writing R Extensions manual](https://cran.r-project.org/doc/manuals/r-release/R-exts.html) is the primary basis for future package structure and checking decisions.

**Automation:** No checker, fixture, registry entry, package metadata, or maintained R source currently implements this rule.

**Semantic remainder:** Classify the intended R surface, choose its public and reusable boundaries, and determine when package ownership becomes justified.

**Exceptions:** Exploratory or generated R content is not maintained source unless an approved repository owner adopts it and defines its applicability.

<br />

## R identifiers (`R.NAMING.IDENTIFIERS`)
**Classification:** `advisory` with deterministic casing and dispatch opportunities.

**Scope:** Project-defined R functions, methods, arguments, objects, columns, fields, properties, class names, constants, files, and exact external boundaries.

Apply the grammatical roles and evidence obligations in [`SOURCE.NAMING.SEMANTICS`](source_layout.md) using R-native forms:
- ordinary functions, methods, arguments, objects, project-created columns, fields, properties, and instances use `snake_case`;
- S3 class strings use `snake_case`, while S3 method functions retain the required `generic.class_name` dispatch punctuation;
- S4, R6, and S7 class definitions or generators use `UpperCamelCase`;
- package-level semantic constants use `SCREAMING_SNAKE_CASE`, recognizing that spelling alone does not make an R binding immutable;
- a leading dot is reserved for intentionally hidden names or required R conventions such as `.onLoad`; and
- `lowerCamelCase` is not used for new project-defined names.

Periods are not ordinary project word separators because they carry S3 method-dispatch meaning. Preserve noncanonical spellings only at exact external data, API, serialization, or interoperability boundaries, and convert to project names internally when practical.

R's `generic.class` lookup is defined by the [R Language Definition](https://stat.ethz.ch/R-manual/R-devel/doc/manual/R-lang.html#Method-dispatching). The ordinary `snake_case` choice and reservation of dots for S3 follow the [tidyverse naming guidance](https://style.tidyverse.org/syntax.html#object-names). S4 class construction is defined by [`setClass()`](https://stat.ethz.ch/R-manual/R-devel/library/methods/html/setClass.html); R6 and S7 class forms are described by the primary [R6](https://r6.r-lib.org/articles/Introduction.html) and [S7](https://rconsortium.github.io/S7/reference/new_class.html) documentation.

**Automation:** No dedicated checker or registry entry exists for R naming. No inventory, proposal ledger, or fixture exists; candidate thresholds remain unset until maintained R symbols provide role-aware evidence.

**Semantic remainder:** Review grammatical role, S3 dispatch intent, class-system role, external ownership, abbreviation clarity, and too-short or too-long candidates.

**Exceptions:** Required external names and established language hooks retain exact spelling at their boundaries. Other deviations require the governed exception process.

<br />

## R comments and documentation (`R.COMMENT.DOCUMENTATION`)
**Classification:** `advisory` with deterministic marker opportunities.

**Scope:** Ordinary R implementation comments, inline comments, documentation blocks, standalone-script documentation, and future package documentation.

Ordinary R comments use `# `, with exactly `#` as an empty separator inside a multiline comment. An inline ordinary comment uses `# ` after the code and remains short; exact pre-marker spacing remains review-owned until representative maintained R source supports a formatter decision. Attach comments according to [`SOURCE.COMMENT.ATTACHMENT`](source_layout.md), use complete prose when the comment forms a sentence, and do not use documentation markers as visual decoration.

Reserve `#'` for roxygen2 documentation. A future R package should use roxygen2 when generated `.Rd` help, namespace declarations, package collation, or object-system documentation is selected as part of the package architecture. A standalone script may use ordinary attached comments and an owned script-level help surface instead; it does not acquire roxygen2 merely for visual consistency.

Documentation must communicate applicable purpose, inputs, outputs or returned values, side effects, failures, runtime requirements, notes, references, and examples using the concepts in [`HELP.SECTION.SCHEMA`](help.md). Empty boilerplate and name-repeating documentation are prohibited.

Roxygen2 block and generated-document behavior follows the [roxygen2 documentation](https://roxygen2.r-lib.org/articles/roxygen2.html).

**Automation:** No dedicated checker or registry entry exists for R comments or documentation. No roxygen2, documentation-coverage, rendered-documentation, or fixture workflow is selected.

**Semantic remainder:** Decide when a comment should become maintained documentation, whether roxygen2 is appropriate to the project form, and which documentation concepts apply to each object.

**Exceptions:** Generated, vendored, literal-fixture, and specialized directive syntax requires an owned applicability exclusion before maintained R source is introduced.

<br />

## R source layout (`R.SOURCE.LAYOUT`)
**Classification:** `semantic-only`.

**Scope:** Assignments, calls, pipelines, constructed vectors, lists, data frames, formulas, control flow, guards, transfers, and top-level organization in maintained R source.

Realize [`SOURCE.LAYOUT.PARAGRAPHS`](source_layout.md) using two-space R indentation unless an approved formatter decision establishes another compatible source form. Keep related assignments, vectorized operations, and short pipeline stages together; separate setup, acquisition, validation, transformation, constructed-object assembly, side effects, and results when they are distinct phases.

Put a blank-line boundary before and after each noncompact sibling `if`, `for`, `while`, `repeat`, or handler-oriented control-flow unit. Keep `else` attached to its owning `if`, and retain a direct single-transfer guard as a compact exception.

Review calls, assignments, pipelines, and constructed vectors, lists, and data frames for both hyper-density and mechanical fragmentation. A long pipeline is not made readable merely by placing every stage on a separate line; expose phase changes with intermediate names or helpers when the stages represent distinct transformations.

Keep a call on one line when it fits and remains readable. Once broken, place the opening call on its own line, one argument on each indented line, and the closing `)` on its own aligned line. Multiline R calls do not use the Python/Rust trailing comma:
```r
rows <- load_table(
  args$tbl_met,
  skp_pfx = skp_pfx,
  verbose = args$verbose,
  other_arg = args$other_val
)
```

Use explicit `return()` when it makes a maintained function's terminal result or an early transfer clear. External code may use R's implicit final-value semantics; do not rewrite it without ownership. Apply compact-guard and transfer separation semantically to `return()`, `stop()`, `break`, and `next`.

The multiline call form and ordinary naming baseline are grounded in the [tidyverse syntax guide](https://style.tidyverse.org/syntax.html#long-lines). R evaluation and function behavior remain governed by the [R Language Definition](https://stat.ethz.ch/R-manual/R-devel/doc/manual/R-lang.html).

**Automation:** No formatter, parser-aware checker, density producer, candidate threshold, fixture, or registry entry exists for R. `X`, `Y`, and `Z` values remain unset pending representative R source and a focused pilot.

**Semantic remainder:** Review paragraph boundaries, pipeline stages, constructed-object density, guard compactness, result clarity, and whether wrapping or extraction is the clearer repair.

**Exceptions:** Language-required syntax and approved formatter output may refine source form only after the formatter relationship and representative evidence are recorded here.

<br />

## R toolchain promotion (`R.TOOLCHAIN.PROMOTION`)
**Classification:** `semantic-only`.

**Scope:** Formatter, documentation, type-analysis, lint, test, package-check, fixture, registry, and migration decisions for future maintained R source.

Do not select tools before the maintained R role and representative source are known. The promotion decision must evaluate:
- formatter behavior and configuration, including whether [styler](https://styler.r-lib.org/) preserves the approved multiline and semantic-grouping policy;
- documentation generation and validation, including roxygen2 only when package or object documentation needs it;
- type-contract needs and whether annotations, runtime validation, type analysis, or no additional type tool best matches the source;
- lint scope and configuration, including [lintr](https://lintr.r-lib.org/) only after findings are mapped to normative owners;
- unit and interface testing, including [testthat](https://testthat.r-lib.org/) only after test roles and repository-runner integration are defined; and
- package layout and `R CMD check` only when a package boundary exists.

Promotion requires selected versions and invocation options, positive and negative fixtures for repository-added behavior, current-source findings, representative formatter previews, threshold evidence, applicability exclusions, registry coverage metadata, focused contracts, idempotence, and a bounded migration plan. Tool defaults do not silently become repository policy.

**Automation:** No R tool is selected, installed, configured, registered, or claimed as conformance evidence by this standard.

**Semantic remainder:** Choose the smallest toolchain that proves the maintained R boundary without duplicating upstream behavior or obscuring semantic review.

**Exceptions:** A temporary exploratory tool invocation is not repository conformance and must not mutate maintained source or establish policy by precedent.
