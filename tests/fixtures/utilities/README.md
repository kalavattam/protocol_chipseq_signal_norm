
# Utility-script test fixtures
These fixtures are inputs for the tests of the small utility scripts in `bin/`: `symlink_files.sh`, `compress_remove_files.sh`, and `find_files.sh`. They need no managed environment.

Regenerate fixtures from the repository root with:
```bash
bash tests/fixtures/utilities/make.sh
```

Generated fixture outputs are ignored by Git. `tests/run_tests.sh` regenerates this fixture set automatically when required inputs are missing.

<br />

## Files
Readable provenance:
- `make.sh`

Inputs for dry runs, in `tests/contract/interfaces/test_dry_run_commands.sh`:
- `text/sample.txt`: one line, for `symlink_files.sh`
- `log/large.log`: 2,048 bytes, which `compress_remove_files.sh` would compress if its dry run acted
- `log/empty.log`: empty, which it would delete

The test copies `log/` before each run, since a dry run that acted would change it.

FASTQ names, in `tests/integration/local/utilities/test_find_files_opt_applicability.sh`:
- `fastq/find/a_R1.fastq.gz` and `fastq/find/a_R2.fastq.gz`: one record each, since `find_files.sh` matches names only

<br />

## Current and deferred test coverage
Current coverage:
- dry runs of `symlink_files.sh` and `compress_remove_files.sh` change no file;
- `find_files.sh` warns that `--shw_str` and `--fastqs` are ignored when it prints no one-line string.

Deferred:
- inputs that other utility-script tests still build inline, until they are converted.
