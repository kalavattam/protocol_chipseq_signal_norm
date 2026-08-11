
# Rendered-help test fixtures
This fixture is the expected Examples section of `bin/combine_parts_scaling_factor.sh --help`, held byte-for-byte so a change to the rendered surface has to be made deliberately.

It is intentionally small and hand-checkable. Running `make.sh` regenerates the fixture set deterministically.

Regenerate fixtures from the repository root with:
```bash
bash tests/fixtures/help/make.sh
```

Generated fixture outputs are ignored by Git. `tests/run_tests.sh` regenerates this fixture set automatically when required inputs are missing.

`expected/` names the role the file plays rather than a verdict, because the file is not an input a rule reads: it is the result the rule's extractor must produce.

**The text is authored, not derived, and that is the whole point.** "Author text and derive everything else" stops here. Deriving this file by running the same extraction the contract runs would leave the contract comparing the extractor's output against the extractor's output, which passes for any output at all — including an empty one. A golden file is evidence only while a person wrote it, so the Examples section is typed into `make.sh` and changes only when someone decides the rendered help should change.

<br />

## Files
Readable provenance:
- `make.sh`

Generated expected output:
- `expected/combine_parts_scaling_factor.examples.txt`: the two-example Examples section, including the backslash continuations and the triple-single-quote fence rows

<br />

## Expected rendering behavior
`dev/audit/help_heredoc_reflow --extract-rendered-examples` applied to the captured `--help` output of `bin/combine_parts_scaling_factor.sh` must reproduce this file exactly. A mismatch is reported as a unified diff rather than as a bare inequality, because the differences that matter here are whitespace and line breaks.

The section is the repository's representative multiline Examples block: two numbered examples, each a fenced Bash invocation wrapped across continuation lines.

<br />

## Current and deferred test coverage
Current coverage in `tests/contract/interfaces/test_help_output.sh`:
- the extracted Examples section matches this fixture byte-for-byte.

Deferred:
- equivalent golden sections for the remaining shell entry points, which would multiply maintenance without exercising a distinct rule; the structural checks in the same contract already cover every surface.
