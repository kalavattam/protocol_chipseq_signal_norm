
# Embedded-program checker fixtures
These fixtures are for `tests/unit/dev_audit/test_embedded_programs.py`, which tests `dev/audit/embedded_programs.py`, the checker for `SOURCE.EMBED.PROGRAMS`, its Shell and Python realizations, and the recorded-input and inline-stub parts of `TEST.FIXTURES`. They need no managed environment.

Regenerate fixtures from the repository root with:
```bash
bash tests/fixtures/embedded_programs/make.sh
```

Generated fixture outputs are ignored by Git. `tests/run_tests.sh` regenerates this fixture set automatically when required inputs are missing.

<br />

## Files
Readable provenance:
- `make.sh`

Each verdict directory is a small repository root, since the checker judges a file by its path from the root (a test is a `test_*` file under `tests/unit/`, `tests/contract/`, or `tests/integration/`).

`accepted/`: programs within their limits and writes the rule allows:
- `bin/short_programs.sh`: one-line `python -c`, `bash -c`, `awk`, and `sed` programs, and longer programs shown only inside help text in a heredoc
- `tests/contract/test_marked.sh`: copies marked with a reason, two of them continued commands whose marker is on the last line, and redirections that capture output
- `tests/unit/test_short.py`: a long docstring, a three-line literal, and a copy marked on its statement's last line
- `tests/fixtures/demo/make.sh`: a recipe that copies a stub and makes that copy executable

`rejected/`: one case of each finding:
- `bin/long_programs.sh`: Python through a heredoc, a two-line `python -c`, a two-line `bash -c`, a four-line `awk`, a two-line `sed`, `sed -i`, and `perl`
- `tests/contract/test_writes.sh`: unmarked `printf`, `cp`, `: >`, and heredoc writes, and a marker with no reason
- `tests/unit/test_literal.py`: a four-line literal that the test writes to disk, and an f-string with a short and a long piece between its fields, found once
- `tests/fixtures/demo/make.sh`: a recipe that writes a stub inline and makes it executable, and makes files executable in an indented loop

`boundary/`: cases at a limit, or writer text where nothing is written:
- `bin/limits.sh`: a three-line `awk`, and a one-line `bash -c` inside a command substitution
- `tests/contract/test_strings.sh`: `printf … >` inside a multi-line string and inside an `awk` program
- `tests/unit/test_three.py`: a literal and an f-string of exactly three lines
- `tests/fixtures/demo/make.sh`: a recipe whose `chmod` names two copied stubs across a continued line

`nonapplicable/`: files the rule does not reach:
- `tests/fixtures/demo/make.sh`: a recipe whose `chmod` sits inside a script it writes as data
- `tests/fixtures/demo/generated/test_generated.sh`: a generated fixture that holds a Python heredoc
- `tests/integration/helper.py`: a helper, not a test, that writes a file; written last, it is the runner's sentinel

<br />

## Current and deferred test coverage
Current coverage:
- every rejected case is found once, at its line, and nothing is found in the other verdicts;
- the checker exits 1 only on a blocking finding, and its JSON report counts each diagnostic.

Deferred:
- judging whether an `awk` call in `lib/` or `bin/` reads sequencing or alignment records, which `SHELL.EMBED.PROGRAMS` routes to `lib/awk/`; the input can follow the program's closing quote, so the checker leaves it to review.
