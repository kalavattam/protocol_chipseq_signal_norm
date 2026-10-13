
# Diagnostic-context test fixtures
These fixtures are synthetic micro-fixtures for fast, deterministic tests of the diagnostic context that the helpers in `lib/bash/core/format_outputs.sh` print, as prescribed by [`SHELL.DIAGNOSTIC.FORM`](../../../docs/standards/shell.md): `<file>::<function>`, with the run script and `>` first when the defining file is not the run script.

They are intentionally small and hand-checkable. Running `make.sh` regenerates the fixture set deterministically.

Regenerate fixtures from the repository root with:
```bash
bash tests/fixtures/diagnostics/make.sh
```

Generated fixture outputs are ignored by Git. `tests/run_tests.sh` regenerates this fixture set automatically when required inputs are missing.

The recipe writes and does nothing else; which prefix each call must print is asserted in `tests/contract/interfaces/test_diagnostic_context.sh`.

<br />

## Files
Readable provenance:
- `make.sh`

Generated scripts:
- `script/run.sh`: a run script that sources the repository helpers and `lib.sh`, then makes a top-level call or calls the function its second argument names; it defines a function, a function reported on by a library helper, and a user function named `main`
- `script/lib.sh`: a library with one function per helper form (`echo_err`, `echo_err_func` with and without a name, `echo_warn`, `echo_warn_func` and `echo_note_func` for functions off the call stack, `echo_note`) and a top-level note

<br />

## Current and deferred test coverage
Current coverage in `tests/contract/interfaces/test_diagnostic_context.sh`:
- each run-script case prints exactly the expected stderr, writes nothing to stdout, and succeeds;
- under `bash -c`, with `$0` set to `bash` or to a login shell's `-bash`, a library function's diagnostic names the shell, the library, and the function; and
- a real driver's library refusal names the run script, the library, and the function.

Deferred:
- the direct `echo` diagnostics in `source_helpers_*` bootstrap functions, which write their prefix by hand until the sweep tracked in [#37](https://github.com/kalavattam/protocol_chipseq_signal_norm/issues/37).
