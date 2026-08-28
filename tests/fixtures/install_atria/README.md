
# Atria installer fixtures
This generated fake toolchain exercises installer reference resolution and rebuild dispatch without network access or an upstream Julia/Atria build.

`make.sh` creates a fake Git client that offers stable and prerelease tags, a minimal Julia extractor, and required runtime commands. The contract records their invocations in test artifacts and asserts the resolved stable tag.

<br />

## Files
- `make.sh`: Generates the local fake toolchain.
- `tool/`: Ignored generated command stubs.

<br />

## Current and deferred test coverage
The repository contract covers dry-run modes, stable-tag selection, checkout, version verification, Julia-aware current-build promotion, retained-build reporting, and managed PATH snippets. A real upstream installation remains manual confirmation.
