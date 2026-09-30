
# Compute-signal test fixtures
These fixtures are synthetic micro-fixtures for fast, deterministic tests of the compute-signal workflow.

They are intentionally small and hand-checkable. Running `make.sh` from `env_protocol` regenerates the fixture set deterministically.

Regenerate fixtures from the repository root with:
```bash
source activate env_protocol
bash tests/fixtures/compute_signal/make.sh
```

Generated fixture outputs are ignored by Git. When required inputs are missing, `tests/run_tests.sh` regenerates this fixture set automatically.

The fixture set covers ratio-mode tests with plain text bedGraph files and signal/coordinate-mode tests with tiny BAM and CRAM files. These fixtures are not derived from real sequencing data. They were written by hand so that expected behavior is easy to inspect and reason about.

<br />

## Files
Readable provenance:
- `make.sh`
- `reference/tiny.fa`
- `sam/se/tiny_se.sam`
- `sam/pe/tiny_pe.sam`

Ratio bedGraph:
- `bedgraph/ratio_A.bdg`
- `bedgraph/ratio_B.bdg`
- `bedgraph/ratio_headers_A.bdg`
- `bedgraph/ratio_headers_B.bdg`
- `bedgraph/ratio_A.bdg.gz`
- `bedgraph/ratio_B.bdg.gz`

Count bedGraph and per-sample counts:
- `bedgraph/count/tiny_se.bedGraph`
- `bedgraph/count/tiny_se.n_frg.txt`
- `bedgraph/count/tiny_se.n_ovlp.txt`
- `bedgraph/count/tiny_pe.bedGraph`
- `bedgraph/count/tiny_pe.n_frg.txt`
- `bedgraph/count/tiny_pe.n_ovlp.txt`

Alignment fixtures:
- `reference/tiny.fa`
- `reference/tiny.fa.fai`
- `sam/se/tiny_se.sam`
- `sam/pe/tiny_pe.sam`
- `bam/se/tiny_se.bam`
- `bam/se/tiny_se.bam.bai`
- `bam/pe/tiny_pe.bam`
- `bam/pe/tiny_pe.bam.bai`
- `cram/se/tiny_se.cram`
- `cram/se/tiny_se.cram.crai`
- `cram/pe/tiny_pe.cram`
- `cram/pe/tiny_pe.cram.crai`

The `ratio_A.bdg` and `ratio_B.bdg` files are plain four-column bedGraph files with matching bins and no header lines.

The `ratio_headers_A.bdg` and `ratio_headers_B.bdg` files preserve the same data rows, but include header-like lines for `--skp_pfx` coverage. They include default skipped prefixes (`track`, `browser`, and `#`) plus the custom prefix `customHeader`.

The `ratio_A.bdg.gz` and `ratio_B.bdg.gz` files are `gzip -n` copies of the plain pair, for gzip input coverage. `-n` omits the name and timestamp, so regeneration is byte-identical.

The `count/` tracks are what `--method count --siz_bin 10` writes for `sam/se/tiny_se.sam` and `sam/pe/tiny_pe.sam`, written literally rather than by running the tool. Beside each are the two counts a signal run writes and `--csv_pseudo edger` reads back: `n_frg`, the fragment count `N`, and `n_ovlp`, the number of bins those fragments touch, `L`.

| Track     | Fragments              | Bins with value 1 | `N`  | `L`  |
| :---      | :---                   | :---              | :--- | :--- |
| `tiny_se` | `[0, 10)`, `[20, 30)`  | 0, 2              | 2    | 2    |
| `tiny_pe` | `[10, 40)`, `[40, 60)` | 1, 2, 3, 4, 5     | 2    | 5    |

The two tracks differ in `N` and `L`, so an A/B swap changes the derived pseudocount and the pairing test can see one.

<br />

## Row-level expected behavior
| Bin       | A    | B    | Purpose                                                                             |
| :---      | :--- | :--- | :---                                                                                |
| `I:0-10`  | 4    | 2    | Simple unadjusted ratio: `A / B = 2`.                                               |
| `I:10-20` | 0    | 2    | Zero numerator: unadjusted ratio is `0`; with pseudocounts `1:1`, ratio is `1 / 3`. |
| `I:20-30` | 5    | 0    | Zero denominator: test non-finite, stabilized, or filtered behavior.                |
| `I:30-40` | 0    | 0    | Zero-zero bin for `skip_00` behavior before or after scaling.                       |
| `I:40-50` | 2    | 0.5  | Scaling-sensitive finite bin: baseline ratio is `4`.                                |
| `I:50-60` | 1    | 0.04 | Denominator-floor case: with `dep_min = 0.1`, ratio is `10`.                        |
| `I:60-70` | 1    | 3    | Decimal-rounding case: with `dp = 3`, ratio is `0.333`.                             |
| `I:70-80` | 1    | 1    | Stable finite row: ratio is `1`.                                                    |

<br />

## Alignment fixtures
`signal`-mode and `coord`-mode tests use tiny SAM/BAM/CRAM fixtures generated from synthetic SAM and FASTA provenance files.

The SAM and FASTA files serve as readable provenance. The FASTA index, BAM files, BAM indexes, CRAM files, and CRAM indexes are generated with `samtools` by `make.sh`.

<br />

## Expected test behavior
Tests should use these fixtures to verify ratio outputs, BAM-backed signal output, CRAM-backed signal output with an explicit reference FASTA, and local GNU Parallel behavior where gated by `RUN_PARALLEL=1`.

<br />

## Current and deferred test coverage
The test suite covers submit- and execute-layer ratio mode, BAM-backed signal and coordinate output, and CRAM-backed signal and coordinate output with an explicit reference FASTA. It also covers local GNU Parallel ratio dispatch when gated by `RUN_PARALLEL=1`.

`#TODO`: Add Slurm coverage. Consider extending GNU Parallel coverage to BAM- and CRAM-backed signal paths where that broader matrix provides useful regression protection. Add direct shell- and Python-script coverage where applicable.
