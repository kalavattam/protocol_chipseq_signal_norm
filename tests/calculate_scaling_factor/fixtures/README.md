# Calculate-scaling-factor test fixtures
These fixtures are realistic micro-fixtures for fast, deterministic tests of the calculate-scaling-factor part-file combiner.

They are intentionally small and hand-checkable. Running `scripts/make_fixtures.sh`
from `env_protocol` regenerates the fixture set deterministically.

Regenerate fixtures from the repository root with:

```bash
source activate env_protocol
bash tests/calculate_scaling_factor/scripts/make_fixtures.sh
```

Generated fixture outputs are ignored by Git. When required inputs are missing, `tests/scripts/run_smoke_tests.sh` regenerates this fixture set automatically.

The fixture rows are adapted from tracked example scaling-factor tables under
`data/processed/compute_signal/`. They preserve realistic file naming, core
table shapes, scaling values, and alignment depths. Synthetic role-specific
SAM, BAM, and CRAM fixtures exercise Samtools-backed alignment counting.

Spike-in rows also record the canonical coefficient method, such as
`fractional`, alongside the computed value.

<br />

## Files
Readable provenance:
- `scripts/make_fixtures.sh`

Role-specific SE alignment fixtures:
- `sam/IP_A.sc.sam`, `bam/IP_A.sc.bam`
- `sam/IP_A.sp.sam`, `bam/IP_A.sp.bam`
- `sam/in_A.sc.sam`, `bam/in_A.sc.bam`
- `sam/in_A.sp.sam`, `bam/in_A.sp.bam`
- `sam/IP_B.sc.sam`, `bam/IP_B.sc.bam`
- `sam/IP_B.sp.sam`, `bam/IP_B.sp.bam`
- `sam/in_B.sc.sam`, `bam/in_B.sc.bam`
- `sam/in_B.sp.sam`, `bam/in_B.sp.bam`

Reference FASTA:
- `reference/tiny.fa`
- `reference/tiny.fa.fai`

Role-specific PE BAM fixtures:
- `bam/pe/IP_A.sc.bam`, `bam/pe/IP_A.sp.bam`
- `bam/pe/in_A.sc.bam`, `bam/pe/in_A.sp.bam`
- `bam/pe/IP_B.sc.bam`, `bam/pe/IP_B.sp.bam`
- `bam/pe/in_B.sc.bam`, `bam/pe/in_B.sp.bam`

Reference-backed CRAM fixtures:
- `cram/se/*.cram`, `cram/se/*.cram.crai`
- `cram/pe/*.cram`, `cram/pe/*.cram.crai`

Spike-in scaling-factor part files:
- `parts/example_scaling_factors.spike.tsv.part.000000`
- `parts/example_scaling_factors.spike.tsv.part.000002`

siQ-ChIP scaling-factor part files:
- `parts/example_scaling_factors.siq.tsv.part.000000`
- `parts/example_scaling_factors.siq.tsv.part.000002`

Malformed negative-test part file:
- `parts/malformed_scaling_factors.spike.tsv.part.000003`

The nonconsecutive numeric suffixes exercise deterministic numeric ordering when input paths are supplied in reverse order. The malformed part contains too few fields and exercises combiner validation.

<br />

## Alignment-backed spike-in expected behavior
| sample | main IP | spike IP | main input | spike input | `chiprx_alpha_ratio` |
| :----  | ------: | -------: | ---------: | ----------: | -------------------: |
| `A`    | 3       | 1        | 2          | 2           | `2`                  |
| `B`    | 2       | 2        | 3          | 1           | `0.5`                |

The serial execute-layer smoke test uses these generated alignment files
directly. Its baseline case does not provide depth overrides, so the workflow
counts alignments with Samtools and automatically detects the SE layout.

<br />

## Current and deferred smoke-test coverage
The smoke suite covers direct spike-in and siQ-ChIP part-file combination,
numeric ordering, dry-run behavior, malformed input rejection, overwrite
protection, explicit part-file cleanup, and serial execute-layer spike-in
assembly. Serial spike-in coverage exercises every canonical coefficient,
accepted alias normalization, alignment-derived counts, automatic SE/PE
detection, BAM/CRAM input, mixed input lists, and broadcast depth overrides.

The scaling-factor rows intentionally contain only core workflow values. Compute downstream denominator floors separately with `python -m scripts.compute_input_floor`.

`#TODO`: Add siQ-ChIP wrapper coverage, direct submit-worker coverage, GNU Parallel coverage, and Slurm coverage.
