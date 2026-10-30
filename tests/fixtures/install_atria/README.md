
# Atria installer fixtures
This generated fake toolchain exercises installer reference resolution and rebuild dispatch without network access or an upstream Julia/Atria build.

`make.sh` creates a fake Git client that offers stable and prerelease tags, a minimal Julia extractor, and required runtime commands. The contract records their invocations in test artifacts and asserts the resolved stable tag.

<br />

## Files
- `make.sh`: Generates the local fake toolchain.
- `tool/`: Ignored generated command stubs.
- `snippet/legacy_path.sh`: a PATH snippet holding the exports of an earlier installation that carried no markers.
- `snippet/malformed_path.sh`: a PATH snippet whose installer markers are in the wrong order.
- `snippet/sentinel_path.sh` and `snippet/sentinel_atria.txt`: contents that a failed installation must leave unchanged.
- `snippet/stale_manifest.toml`: generated dependency state that an Atria tag update must clear.
- `snippet/julia_0.0.0`: a selected Julia that reports a version no release has, which a failed Julia update must keep.

The installer rewrites a PATH snippet, so the test copies each one before a run.

<br />

## Current and deferred test coverage
The repository contract covers dry-run modes, stable-tag selection, checkout, version verification, Julia-aware current-build promotion, retained-build reporting, managed PATH snippets, non-mutating fail/reuse behavior, and staged invalid-Julia replacement. A real upstream installation remains manual confirmation.
