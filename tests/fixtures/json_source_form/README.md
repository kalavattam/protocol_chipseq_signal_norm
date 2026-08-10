
# JSON source-form test fixtures
These fixtures are synthetic micro-fixtures for fast, deterministic tests of the canonical JSON rendering owned by [`JSON.SOURCE.FORM`](../../../docs/standards/json.md).

They are intentionally small and hand-checkable. Running `make.sh` regenerates the fixture set deterministically.

Regenerate fixtures from the repository root with:
```bash
bash tests/fixtures/json_source_form/make.sh
```

Generated fixture outputs are ignored by Git. `tests/run_tests.sh` regenerates this fixture set automatically when required inputs are missing.

Each negative fixture departs from the canonical rendering along exactly one axis. That isolation is deliberate: a fixture that violated two rules at once would let a checker regression hide behind whichever finding still fired.

Because the outputs are ignored, they are also invisible to the checker itself, whose discovery lists tracked and non-ignored files. A fixture that is deliberately unreadable, wrongly indented, or tab-indented therefore cannot be mistaken for maintained source violating the rule it exists to exercise. This is the same isolation the AI-attribution fixtures rely on, and it is what allows a negative fixture to be pure text with no tooling dependency.

<br />

## Files
Readable provenance:
- `make.sh`

Generated source fixtures:
- `source/canonical.json`: the canonical rendering itself, exercising both treatments the budget selects
- `source/inline_overflow.json`: one array packed onto its key's line, past the budget
- `source/expanded_fits.json`: one structure broken across lines that fits the budget inline
- `source/hybrid_delimiter.json`: one object opened inline and then continued vertically
- `source/wrong_indent.json`: one expanded structure indented by three spaces instead of two
- `source/tab_indent.json`: one expanded structure indented with tabs
- `source/no_trailing_newline.json`: one canonical document with its final newline removed
- `source/duplicate_key.json`: one object declaring the same key twice
- `source/unreadable.json`: one document that is not JSON at all

<br />

## Expected canonical-form behavior
`canonical.json` is the only conforming input. It carries a short structure that stays inline, a long structure that expands, a record-per-line array whose rows fit and are therefore preserved, and a scalar longer than the budget that the rule deliberately leaves alone. Every other fixture reports at least one finding.

`duplicate_key.json` is reported as unreadable rather than reformatted. `json.loads` keeps only the last occurrence of a repeated key, so rewriting the file would silently discard the first value; refusing to read it is what keeps `--fix` from destroying data.

<br />

## Current and deferred test coverage
Current coverage in `tests/unit/dev_audit/test_json_source_form.py`:
- the canonical rendering reports nothing, and re-rendering it is a fixpoint;
- each negative fixture reports its own class and no other;
- rewriting every fixture that parses preserves the parsed value exactly; and
- a rewritten fixture reports nothing on a second pass.

Deferred:
- non-ASCII and escape-sequence rendering, which the maintained corpus exercises directly and which no fixture would exercise more sharply; and
- multi-file discovery behavior, which belongs to the repository gate rather than to a unit fixture.
