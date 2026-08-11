
# Help-contract test fixtures
These fixtures are synthetic micro-fixtures for fast, deterministic tests of the help-surface contracts owned by [`help.md`](../../../docs/standards/help.md).

They are intentionally small and hand-checkable. Running `make.sh` regenerates the fixture set deterministically.

Regenerate fixtures from the repository root with:
```bash
bash tests/fixtures/help_contracts/make.sh
```

Generated fixture outputs are ignored by Git. `tests/run_tests.sh` regenerates this fixture set automatically when required inputs are missing.

`accepted/` holds the surfaces the contracts must read without complaint. `cases.json` sits at the fixture root rather than under a verdict directory because it is the record of every verdict at once — accepted, rejected, boundary, and non-applicable — rather than an input carrying one of them.

Because the outputs are ignored, they are invisible to the checkers' own discovery, so a fixture help surface cannot be mistaken for a maintained tool's help surface.

<br />

## Files
Readable provenance:
- `make.sh`

Generated verdict record:
- `cases.json`: the accepted, rejected, boundary, and non-applicable readings the help-contract rules must produce, in the canonical form owned by `JSON.SOURCE.FORM`

Generated help surfaces the contracts must accept:
- `accepted/shell.sh`: one shell help function whose body is a heredoc carrying the standard Usage and Parameters sections
- `accepted/python.py`: one `add_argument` help literal using the accepted inline-token form
- `accepted/python_callable.py`: one callable whose Examples section is source-language rather than a rendered CLI invocation

<br />

## Expected help-contract behavior
All three generated surfaces are conforming inputs. `accepted/shell.sh` exercises the structured Usage and Parameters rows, and `accepted/python_callable.py` exercises the callable example form, whose fingerprint covers both the source and the displayed expected result.

`cases.json` is a record rather than an input. It states what each verdict means for the audience, token, parameter, example, and Usage-row obligations, so a later change to a rule can be checked against a written expectation instead of against memory.

<br />

## Current and deferred test coverage
Current coverage in `tests/unit/dev_audit/test_help_contracts.py` and `tests/unit/dev_audit/test_help_examples.py`:
- the shell surface parses into its declared Usage and Parameters rows; and
- the callable surface yields two example blocks and two signatures, and a mutation of either is detected.

Deferred:
- `cases.json` and `accepted/python.py` have no consuming assertion yet. They are generated so the record and the surface stay together and stay canonical, and wiring them to assertions is the remaining work.
