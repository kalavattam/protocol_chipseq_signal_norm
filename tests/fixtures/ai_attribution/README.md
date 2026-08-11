
# AI-attribution test fixtures
These fixtures are synthetic micro-fixtures for fast, deterministic tests of the bounded source-header attribution rules owned by [`SOURCE.HEADER.AI_ATTRIBUTION`](../../../docs/standards/source_headers.md).

They are intentionally small and hand-checkable. Running `make.sh` regenerates the fixture set deterministically.

Regenerate fixtures from the repository root with:
```bash
bash tests/fixtures/ai_attribution/make.sh
```

Generated fixture outputs are ignored by Git. `tests/run_tests.sh` regenerates this fixture set automatically when required inputs are missing.

Each directory names the verdict the checker must return for the header inside it: the two credited forms are accepted, and a coherent no-AI profile is outside what the rule reaches rather than a header it accepts. Every fixture shares one bounded header and differs only in its attribution block, which is the subject under test. That isolation is deliberate: a difference in any other header row would confound a finding about attribution with a finding about header structure.

Because the outputs are ignored, they are also invisible to the attribution checker itself, whose discovery lists tracked and non-ignored files. A fixture whose header is deliberately unattributed therefore cannot be mistaken for a maintained source in violation of the rule it exists to exercise.

<br />

## Files
Readable provenance:
- `make.sh`

Generated headers the checker must accept:
- `accepted/single_vendor.sh`: one vendor credited in the bounded prose form
- `accepted/multi_vendor.sh`: two vendors credited in the lead-in and semicolon-list form, in first-use order

Generated headers the rule does not reach:
- `non_applicable/no_attribution.sh`: a coherent no-AI profile that declares no attribution at all

<br />

## Expected attribution behavior
All three fixtures are conforming inputs. `accepted/single_vendor.sh` and `accepted/multi_vendor.sh` exercise the two accepted representations, and `non_applicable/no_attribution.sh` exercises the profile that must report null observed attribution rather than a finding.

Trailer-agreement behavior is exercised against injected evidence rather than against these files, because a fixture has no commit history of its own.

<br />

## Current and deferred test coverage
Current coverage in `tests/unit/dev_audit/test_ai_attribution.py`:
- both accepted attribution representations parse and report no finding;
- a no-AI profile reports null observed attribution;
- a focused commit crediting a vendor the header omits reports one finding; and
- a broad checkpoint crediting that vendor reports nothing.

Deferred:
- malformed and duplicate-token cases, which remain inline in the unit test because each is a one-line perturbation of a conforming block rather than an independent input; and
- Python-language fixtures, which would duplicate the shell cases without exercising a distinct rule.
