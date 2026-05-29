# Compute-signal test fixtures
These fixtures are synthetic micro-fixtures for fast, deterministic tests of the compute-signal workflow.

They are intentionally small and hand-checkable. Running `scripts/make_fixtures.sh` from `env_protocol` regenerates the fixture set deterministically.

Regenerate fixtures from the repository root with:
```bash
source activate env_protocol
bash tests/compute_signal/scripts/make_fixtures.sh
```

For now, tests assume these fixtures are already present. A later test-runner refactor will allow `tests/scripts/run_smoke_tests.sh` to generate missing fixtures automatically when required inputs are not detected (`#TODO`).

The fixture set covers ratio-mode tests with plain text bedGraph files and signal/coordinate-mode tests with tiny BAM and CRAM files. These fixtures are not derived from real sequencing data. They were written by hand so that expected behavior is easy to inspect and reason about.

<br />

## Files
Readable provenance:
- `scripts/make_fixtures.sh`
- `reference/tiny.fa`
- `sam/tiny_se.sam`
- `sam/tiny_pe.sam`

Ratio bedGraph:
- `bedgraph/ratio_A.bdg`
- `bedgraph/ratio_B.bdg`
- `bedgraph/ratio_headers_A.bdg`
- `bedgraph/ratio_headers_B.bdg`

Alignment fixtures:
- `reference/tiny.fa`
- `reference/tiny.fa.fai`
- `sam/tiny_se.sam`
- `sam/tiny_pe.sam`
- `bam/tiny_se.bam`
- `bam/tiny_se.bam.bai`
- `bam/tiny_pe.bam`
- `bam/tiny_pe.bam.bai`
- `cram/tiny_se.cram`
- `cram/tiny_se.cram.crai`
- `cram/tiny_pe.cram`
- `cram/tiny_pe.cram.crai`

The `ratio_A.bdg` and `ratio_B.bdg` files are plain four-column bedGraph files with matching bins and no header lines.

The `ratio_headers_A.bdg` and `ratio_headers_B.bdg` files preserve the same data rows, but include header-like lines for `--skp_pfx` coverage. They include default skipped prefixes (`track`, `browser`, and `#`) plus the custom prefix `customHeader`.

<br />

## Row-level expected behavior
| Bin       | A  | B    | Purpose |
|:---       |:---|:---  |:---     |
| `I:0-10`  | 4  | 2    | Simple unadjusted ratio: `A / B = 2`.                                               |
| `I:10-20` | 0  | 2    | Zero numerator: unadjusted ratio is `0`; with pseudocounts `1:1`, ratio is `1 / 3`. |
| `I:20-30` | 5  | 0    | Zero denominator: test non-finite, stabilized, or filtered behavior.                |
| `I:30-40` | 0  | 0    | Zero-zero bin for `skip_00` behavior before or after scaling.                       |
| `I:40-50` | 2  | 0.5  | Scaling-sensitive finite bin: baseline ratio is `4`.                                |
| `I:50-60` | 1  | 0.04 | Denominator-floor case: with `dep_min = 0.1`, ratio is `10`.                        |
| `I:60-70` | 1  | 3    | Decimal-rounding case: with `dp = 3`, ratio is `0.333`.                             |
| `I:70-80` | 1  | 1    | Stable finite row: ratio is `1`.                                                    |

<br />

## Alignment fixtures
`signal`-mode and `coord`-mode tests use tiny SAM/BAM/CRAM fixtures generated from synthetic SAM and FASTA provenance files.

The SAM and FASTA files serve as readable provenance. The FASTA index, BAM files, BAM indexes, CRAM files, and CRAM indexes are generated with `samtools` by `scripts/make_fixtures.sh`.

CRAM support should not be treated as a distant optional feature. If the current compute-signal path cannot pass a reference FASTA cleanly to every layer that needs it, that gap should be exposed by the tiny CRAM fixture tests and addressed as part of the compute-signal testing/refactor work.

<br />

## Expected smoke-test behavior
Smoke tests should use these fixtures to verify ratio outputs, BAM-backed signal output, CRAM-backed signal output with an explicit reference FASTA, and local GNU Parallel behavior where gated by `RUN_PARALLEL=1`.

<br />

## Deferred fixture batches
`#TODO`: Later compute-signal batches should add...
