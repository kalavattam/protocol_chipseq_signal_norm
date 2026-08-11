
# Markdown policy test fixtures
These fixtures are synthetic micro-fixtures for fast, deterministic tests of the Markdown rules owned by [`markdown.md`](../../../docs/standards/markdown.md) and of the formatter in `dev/tools/markdown_format.py`.

They are intentionally small and hand-checkable. Running `make.sh` regenerates the fixture set deterministically.

Regenerate fixtures from the repository root with:
```bash
bash tests/fixtures/markdown/make.sh
```

Generated fixture outputs are ignored by Git. `tests/run_tests.sh` regenerates this fixture set automatically when required inputs are missing.

Each directory names the verdict the checker must return for the documents inside it, so a verdict with no fixture is visible as an absent directory rather than as a filename nobody happened to write. `format/` is the one directory that does not name a verdict, because its pairs claim a rewrite rather than a reading.

Because the outputs are ignored, they are also invisible to the Markdown checker itself, whose discovery lists tracked and non-ignored files. A document under `rejected/` therefore cannot be mistaken for maintained documentation violating the rule it exists to exercise, and the formatter needs no exclusion list of its own.

<br />

## Files
Readable provenance:
- `make.sh`

Generated documents the checker must accept:
- `accepted/anchors.md`: canonical anchor placement and casing
- `accepted/basic.md`: the ordinary document skeleton
- `accepted/blockquotes.md`: blockquote spacing and nesting
- `accepted/nested_fences.md`: a fenced block containing fences
- `accepted/section_boundaries.md`: canonical section breaks
- `accepted/standard_section.md`: the standard-document section form

Generated documents the checker must reject:
- `rejected/anchors.md`: a noncanonical anchor
- `rejected/blockquotes.md`: blockquote spacing faults
- `rejected/nested_fences.md`: a mishandled inner fence
- `rejected/section_boundaries.md`: six missing section breaks
- `rejected/spacing.md`: file-boundary, heading, colon, and list spacing faults
- `rejected/standard_section.md`: a malformed standard section
- `rejected/unclosed_fence.md`: a fence that never closes

Generated documents on a decision edge, which must still pass:
- `boundary/blockquotes.md`
- `boundary/nested_fences.md`

Generated documents the rule does not reach:
- `non_applicable/blockquotes.md`
- `non_applicable/nested_fences.md`

Generated formatter pairs:
- `format/input.md` and `format/expected.md`
- `format/blockquotes_input.md` and `format/blockquotes_expected.md`
- `format/nested_fences_input.md` and `format/nested_fences_expected.md`

<br />

## Expected Markdown behavior
Every document under `accepted/`, `boundary/`, and `non_applicable/` reports no deterministic finding. Each document under `rejected/` reports the rules named for it in the unit test and no others, so a regression that broadens one rule is separable from a regression that narrows another.

Each `format/` pair pins one rewrite and its fixpoint: formatting the input yields the expected document, and formatting the expected document yields itself.

<br />

## Current and deferred test coverage
Current coverage in `tests/unit/dev_audit/test_markdown_policy.py`:
- accepted documents report no deterministic finding;
- rejected documents report their own owned rules;
- boundary documents pass and survive formatting unchanged; and
- each formatter pair matches its expected output and is idempotent.

Deferred:
- table and link rendering, which the maintained corpus exercises directly and which no fixture would exercise more sharply; and
- multi-file discovery behavior, which belongs to the repository gate rather than to a unit fixture.
