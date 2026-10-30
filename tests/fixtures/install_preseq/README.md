
# preseq installer test fixtures
These fixtures let `tests/contract/repository/test_install_preseq.sh` run `install/scripts/install_preseq.sh` against a stand-in environment prefix, with no compiler, `make`, or `htslib` that works. They need no managed environment.

Regenerate fixtures from the repository root with:
```bash
bash tests/fixtures/install_preseq/make.sh
```

Generated fixture outputs are ignored by Git. `tests/run_tests.sh` regenerates this fixture set automatically when required inputs are missing.

<br />

## Files
Readable provenance:
- `make.sh`

A file named like a tarball:
- `tarball/not_a_tarball.tar.gz`: text, so its SHA-256 never matches the pinned one

Prefix states, each with `bin/make`, `include/htslib/sam.h`, and `lib/libhts.so.3`. Every stub exits 1, since the installer must stop at its checks before it would use them:
- `prefix/bare/`: also has `bin/clang++`
- `prefix/no_compiler/`: has no `bin/clang++`
- `prefix/installed/`: also has `bin/clang++`, and a `bin/preseq` that reports version 3.2.0

The installer writes into a prefix, so the test copies a state before each case.

<br />

## Current and deferred test coverage
Current coverage:
- a tarball whose checksum does not match is refused before any build, and nothing is installed;
- a prefix without a compiler is reported;
- an existing installation is left alone under `--if_exists fail`, accepted under `reuse` when its version is the pinned one, and refused under `reuse` otherwise.

Deferred:
- a real build, which needs a working compiler and `htslib`.
