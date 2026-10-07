
# Shell source-form test fixtures
These fixtures are synthetic micro-fixtures for fast, deterministic tests of the comment-wrap facet of the bounded shell source-form checker in `dev/audit/shell_source_form.py`. The facet realizes [`SOURCE.PROSE.WRAP`](../../../docs/standards/source_layout.md) for ordinary Shell comments under [`SHELL.COMMENT.FORM`](../../../docs/standards/shell.md).

They are intentionally small and hand-checkable. Running `make.sh` regenerates the fixture set deterministically.

Regenerate fixtures from the repository root with:
```bash
bash tests/fixtures/shell_source_form/make.sh
```

Generated fixture outputs are ignored by Git. `tests/run_tests.sh` regenerates this fixture set automatically when required inputs are missing.

Each directory names the verdict the checker must return for the script inside it, except `format/`, whose pair claims a rewrite by `dev/tools/comment_wrap_format.py` rather than a verdict. The recipe writes and does nothing else; which lines the checker must report is asserted in `tests/unit/dev_audit/test_shell_source_form.py`.

Because the outputs are ignored, they are invisible to the checker's own discovery, which asks Git for maintained shell paths, so the rejected script never reaches a repository-wide run.

<br />

## Files
Readable provenance:
- `make.sh`

Generated scripts the checker must accept:
- `accepted/comment_wrap.sh`: greedily wrapped comment prose, a blank-comment paragraph separator, a list with an indented continuation, a `TODO:` item, a function inventory, a ShellCheck directive, an indented comment inside a function, and comment-like lines in a heredoc body

Generated scripts on the exact edge of the fit test, which must also pass:
- `boundary/comment_wrap.sh`: a next word that would end at column 80, a bracketed unit whose first word fits but whole unit does not, a hyphenated join one column too long, the line that ends a list, an indented display, and a fenced block

Generated scripts the checker must reject:
- `rejected/comment_wrap.sh`: five blocks, each with one premature break: plain prose whose next word ends at column 79, a list-item continuation, a hyphenated join that fits only without a separator, a bracketed unit that fits whole, and a comment indented inside a function

Generated formatter pair:
- `format/comment_wrap_input.sh` and `format/comment_wrap_expected.sh`: premature breaks in prose, in a list item, and before a bracketed unit, refilled by hand in the expected file; a paragraph with a line ending in a hyphen and a heredoc body stay as written

Every case lies below the 16 header rows the checker skips.

<br />

## Expected checker behavior
The accepted and boundary scripts produce no findings. The rejected script produces exactly one `SHELL.COMMENT.FORM` wrap finding in each of its five blank-line-separated blocks.

Formulas and other multiword units the facet cannot recognize, comment usefulness, and whether a break is semantically deliberate are deliberately outside what these fixtures decide.

<br />

## Current and deferred test coverage
Current coverage in `tests/unit/dev_audit/test_shell_source_form.py`:
- accepted and boundary scripts report no finding; and
- each rejected block reports exactly one wrap finding.

Current coverage in `tests/unit/dev_audit/test_comment_wrap_format.py`:
- the formatter turns the input into the expected file, is idempotent, and leaves for review exactly the lines the checker still reports.

The checker imports its wrap rules from the formatter, so the two cannot disagree.

Deferred:
- repository-wide inspection, which reports the current wrap findings as an advisory inventory rather than a gate.
