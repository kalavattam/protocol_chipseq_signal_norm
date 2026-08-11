
# Python source-policy test fixtures
These fixtures are synthetic micro-fixtures for fast, deterministic tests of the bounded Python source-policy checker in `dev/audit/python_source_policy.py`, which owns the rules recorded in [`python.md`](../../../docs/standards/python.md) and [`source_layout.md`](../../../docs/standards/source_layout.md).

They are intentionally small and hand-checkable. Running `make.sh` regenerates the fixture set deterministically.

Regenerate fixtures from the repository root with:
```bash
bash tests/fixtures/python_source_policy/make.sh
```

Generated fixture outputs are ignored by Git. `tests/run_tests.sh` regenerates this fixture set automatically when required inputs are missing.

Each directory names the verdict the checker must return for the module inside it. `format/` is the one directory that does not name a verdict, because its pairs claim a rewrite by the help-literal producer rather than a reading by the checker.

The recipe writes and does nothing else. An earlier revision ran the checker over its own output and required eleven rule owners before it would report success, which inverted the separation the standard states: `make.sh` generates, and tests and checkers validate. Those claims now live in `tests/unit/dev_audit/test_python_source_policy.py`, where deciding which owners a fixture must provoke is a property of the checker rather than a precondition of having a fixture at all.

Because the outputs are ignored, they are also invisible to the checker's own discovery. `maintained_python_paths()` is the exception that proves the point: it globs the filesystem rather than asking Git, so it cannot see that these modules are ignored and excludes `tests/fixtures` explicitly.

<br />

## Files
Readable provenance:
- `make.sh`

Generated modules the checker must accept:
- `accepted/canonical.py`: canonical docstrings, annotations, comments, strings, definition topology, multiline parameter and construction delimiters, a multiline NumPy type, compact nested tuples, and a greedily wrapped `parse_args()` help literal
- `accepted/exceptions.py`: governed source-header and directive exclusions, a necessary raw docstring, embedded-double-quote selection, and the multiline triple-single literal-content exception

Generated modules on the exact source-form edge, which must also pass:
- `boundary/source_form.py`: parameter markers, generated 79- and 80-column cases, and greedy docstring prose whose structural boundaries — entry headers, bullet continuation indentation, multiline textual types, and doctest rows — are not prose breaks

Generated modules the checker must reject:
- `rejected/deterministic_owners.py`: bounded near misses for every deterministic owner, including docstring alignment, a misaligned multiline NumPy type closer, multiline parameters, and an overlong single static help literal

Generated formatter pairs:
- `format/help_input.py` and `format/help_expected.py`
- `format/help_unicode_input.py` and `format/help_unicode_expected.py`

<br />

## Expected source-policy behavior
The accepted and boundary modules produce zero deterministic findings. The rejected module produces at least one finding for each of `HELP.PROSE.SENTENCES`, `PY.CLI.HELP.LAYOUT`, `PY.COMMENT.FORM`, `PY.DOCSTRING.LAYOUT`, `PY.DOCSTRING.NUMPY`, `PY.NAMING.IDENTIFIERS`, `PY.SOURCE.LAYOUT`, `PY.STRING.QUOTES`, `PY.TYPE.ANNOTATIONS`, `SOURCE.DELIMITED.MULTILINE`, and `SOURCE.PROSE.WRAP`.

Semantic readability, docstring usefulness, naming quality, abbreviation clarity, paragraph boundaries, and candidate dispositions are deliberately outside what these fixtures decide.

<br />

## Current and deferred test coverage
Current coverage in `tests/unit/dev_audit/test_python_source_policy.py` and `tests/unit/dev_audit/test_python_help_format.py`:
- accepted and boundary modules report no deterministic finding;
- the rejected module provokes every expected owner;
- analysis of one module is stable across repeated runs; and
- each formatter pair is greedy, value-preserving, and idempotent, including the width-sensitive non-ASCII pair.

The focused unit tests add generated in-memory cases for exact 79-column help and comment boundaries, rendered help equivalence, `self` receiver scope, generator-expression exclusions, X/Y/Z thresholds, candidate-schema stability, and checker or producer idempotence.

Deferred:
- semantic review of the migration cohort, which is human-owned and recorded separately; and
- repository-wide inspection, which belongs to `tests/contract/repository/test_python_source_policy_gate.sh`.
