
# Trim-fastqs test fixtures
These fixtures are synthetic micro-fixtures for fast, deterministic tests of the trim-fastqs workflow.

They are intentionally small and hand-checkable. Running `scripts/make_fixtures.sh` regenerates the fixture set deterministically.

Regenerate fixtures from the repository root with:
```bash
bash tests/trim_fastqs/scripts/make_fixtures.sh
```

Generated fixture outputs are ignored by Git. With `RUN_ATRIA=1`, `tests/scripts/run_smoke_tests.sh` regenerates this fixture set automatically when required inputs are missing.

The fixture set currently focuses on single-end and paired-end Atria trimming smoke tests. The FASTQ records are encoded directly in `scripts/make_fixtures.sh` so that the expected read names and sequences are easy to inspect.

<br />

## Files
Readable provenance:
- `scripts/make_fixtures.sh`

The FASTQ records are generated deterministically by `scripts/make_fixtures.sh` as compressed workflow inputs that can be inspected with `gzip -cd`.

Compressed FASTQ input:
- `fastq/tiny_se.fastq.gz`
- `fastq/tiny_pe_R1.fastq.gz`
- `fastq/tiny_pe_R2.fastq.gz`

The compressed FASTQ files are generated with `gzip -n -c` so gzip headers do not encode timestamps or source filenames.

<br />

## Expected trimming behavior
The single-end FASTQ fixture contains one clean 64-bp read:

| Read name             | Sequence                                                           |
|:---                   |:---                                                                |
| `tiny_trim_se_read_1` | `ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT` |

`submit_trim_fastqs.sh` currently calls Atria with `--length-range 35:500`, so this read is intentionally longer than the minimum length filter. Submit- and execute-layer smoke tests should verify that Atria-backed trimming produces a non-empty `gzip` FASTQ output and preserves the expected read name.

The paired-end FASTQ fixtures contain one clean 64-bp read pair:

| Read name               | Mate | Sequence                                                           |
|:---                     |:---  |:---                                                                |
| `tiny_trim_pe_pair_1/1` | R1   | `ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT` |
| `tiny_trim_pe_pair_1/2` | R2   | `TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA` |

The paired-end submit- and execute-layer smoke tests should pass the two compressed inputs as one comma-delimited `--csv_infile` entry, verify that Atria-backed trimming produces non-empty gzip FASTQ outputs for R1 and R2, and confirm that each output contains the expected paired read name prefix.

<br />

## Current and deferred smoke-test coverage
The smoke suite covers submit- and execute-layer single-end and paired-end Atria trimming when gated by `RUN_ATRIA=1`. It also covers local GNU Parallel dispatch when gated by `RUN_ATRIA=1 RUN_PARALLEL=1`.

`#TODO`: Add remote Slurm coverage gated by `RUN_SLURM=1`.
