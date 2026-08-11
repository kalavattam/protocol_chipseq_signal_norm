
# ShellCheck runner fixtures
These fixtures are synthetic micro-fixtures for fast, deterministic tests of the ShellCheck runner in `dev/audit/run_shellcheck.sh`.

They are intentionally small and hand-checkable. Running `make.sh` regenerates the fixture set deterministically.

Regenerate fixtures from the repository root with:
```bash
bash tests/fixtures/shellcheck/make.sh
```

Generated fixture outputs are ignored by Git. `tests/run_tests.sh` regenerates this fixture set automatically when required inputs are missing.

Neither directory names a verdict, because neither fixture carries one. `script/` holds inputs whose declared language is the whole subject, and `tool/` holds an executable standing in for ShellCheck itself. Naming the role rather than a verdict keeps a stub from being filed as though a rule had accepted it.

The stub exists so the runner's behavior can be exercised without applying repository findings, which would make this test's result depend on the state of the repository it is testing.

<br />

## Files
Readable provenance:
- `make.sh`

Generated scripts, distinguished only by their shebang:
- `script/bash.sh`: a Bash script the runner must classify as Bash
- `script/posix.sh`: a POSIX script the runner must classify as POSIX

Generated tooling:
- `tool/fake_shellcheck.sh`: a ShellCheck stub reporting a fixed version, logging its invocation when `FAKE_SHELLCHECK_LOG` is set, and emitting an empty or one-comment result according to `FAKE_SHELLCHECK_STATUS`

<br />

## Expected runner behavior
The runner must select the stub through the configured prefix rather than through `PATH`, split the two scripts by declared language, parse the stub's structured JSON output, and propagate the stub's exit status. Status `0` yields no comments, status `1` yields exactly one warning-level comment, and any other status is an infrastructure failure that must not be reported as a clean run.

The stub's findings payload is a single line of JSON. It is assembled from fragments inside the stub rather than typed as one long line, so the recipe that writes it stays inside the shell line-length budget while the emitted output is unchanged.

<br />

## Current and deferred test coverage
Current coverage in `tests/contract/repository/test_shellcheck_inventory.sh`:
- explicit executable selection through the managed prefix;
- Bash and POSIX splitting across a synthetic repository, including a tracked file deleted from the working tree;
- structured-finding parsing and inventory accounting; and
- exit-status handling for clean, finding, and infrastructure-failure results.

Deferred:
- real ShellCheck invocation, which the repository inventory performs directly and which a stub cannot stand in for.
