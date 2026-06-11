
# Calculate-scaling-factor test fixtures
These fixtures are realistic micro-fixtures for fast, deterministic tests of the calculate-scaling-factor workflow.

They are intentionally small and hand-checkable. Running `scripts/make_fixtures.sh` from `env_protocol` regenerates the fixture set deterministically.

Regenerate fixtures from the repository root with:
```bash
source activate env_protocol
bash tests/calculate_scaling_factor/scripts/make_fixtures.sh
```

Generated fixture outputs are ignored by Git. When required inputs are missing, `tests/scripts/run_smoke_tests.sh` regenerates this fixture set automatically.

The fixture rows are adapted from tracked example scaling-factor tables under `data/processed/compute_signal/`. They preserve realistic file naming, core table shapes, scaling values, and alignment depths. Synthetic role-specific SAM, BAM, and CRAM fixtures exercise Samtools-backed alignment counting.

Alignment fixture basenames are biologically minded because siQ-ChIP tests exercise production filename parsing. PE fixtures use Hho1-like names, and SE fixtures use Brn1-like names.

Spike-in rows also record the canonical coefficient method used alongside the computed value. The default fixture examples use `chiprx_alpha_ratio`.

Minimal siQ-ChIP metadata tables exercise `parse_metadata_siqchip.py` with the production parsing YAML in `data/raw/docs/parse_metadata_siqchip.yml`. These tables are intentionally small and cover canonical column names, compressed input, skip-prefix handling, curated aliases, ambiguous matches, treatment-aware matching, and malformed negative cases.

<br />

## Files
Readable provenance:
- `scripts/make_fixtures.sh`

Role-specific SE alignment fixtures:
- `sam/se/IP_WT_log_Brn1_rep1.sc.sam`, `bam/se/IP_WT_log_Brn1_rep1.sc.bam`, `bam/se/IP_WT_log_Brn1_rep1.sc.bam.bai`
- `sam/se/IP_WT_log_Brn1_rep1.sp.sam`, `bam/se/IP_WT_log_Brn1_rep1.sp.bam`, `bam/se/IP_WT_log_Brn1_rep1.sp.bam.bai`
- `sam/se/in_WT_log_Brn1_rep1.sc.sam`, `bam/se/in_WT_log_Brn1_rep1.sc.bam`, `bam/se/in_WT_log_Brn1_rep1.sc.bam.bai`
- `sam/se/in_WT_log_Brn1_rep1.sp.sam`, `bam/se/in_WT_log_Brn1_rep1.sp.bam`, `bam/se/in_WT_log_Brn1_rep1.sp.bam.bai`
- `sam/se/IP_WT_log_Brn1_rep2.sc.sam`, `bam/se/IP_WT_log_Brn1_rep2.sc.bam`, `bam/se/IP_WT_log_Brn1_rep2.sc.bam.bai`
- `sam/se/IP_WT_log_Brn1_rep2.sp.sam`, `bam/se/IP_WT_log_Brn1_rep2.sp.bam`, `bam/se/IP_WT_log_Brn1_rep2.sp.bam.bai`
- `sam/se/in_WT_log_Brn1_rep2.sc.sam`, `bam/se/in_WT_log_Brn1_rep2.sc.bam`, `bam/se/in_WT_log_Brn1_rep2.sc.bam.bai`
- `sam/se/in_WT_log_Brn1_rep2.sp.sam`, `bam/se/in_WT_log_Brn1_rep2.sp.bam`, `bam/se/in_WT_log_Brn1_rep2.sp.bam.bai`

Role-specific PE alignment fixtures:
- `sam/pe/IP_WT_G1_Hho1_6336.sc.sam`, `bam/pe/IP_WT_G1_Hho1_6336.sc.bam`, `bam/pe/IP_WT_G1_Hho1_6336.sc.bam.bai`
- `sam/pe/IP_WT_G1_Hho1_6336.sp.sam`, `bam/pe/IP_WT_G1_Hho1_6336.sp.bam`, `bam/pe/IP_WT_G1_Hho1_6336.sp.bam.bai`
- `sam/pe/in_WT_G1_Hho1_6336.sc.sam`, `bam/pe/in_WT_G1_Hho1_6336.sc.bam`, `bam/pe/in_WT_G1_Hho1_6336.sc.bam.bai`
- `sam/pe/in_WT_G1_Hho1_6336.sp.sam`, `bam/pe/in_WT_G1_Hho1_6336.sp.bam`, `bam/pe/in_WT_G1_Hho1_6336.sp.bam.bai`
- `sam/pe/IP_WT_G1_Hho1_6337.sc.sam`, `bam/pe/IP_WT_G1_Hho1_6337.sc.bam`, `bam/pe/IP_WT_G1_Hho1_6337.sc.bam.bai`
- `sam/pe/IP_WT_G1_Hho1_6337.sp.sam`, `bam/pe/IP_WT_G1_Hho1_6337.sp.bam`, `bam/pe/IP_WT_G1_Hho1_6337.sp.bam.bai`
- `sam/pe/in_WT_G1_Hho1_6337.sc.sam`, `bam/pe/in_WT_G1_Hho1_6337.sc.bam`, `bam/pe/in_WT_G1_Hho1_6337.sc.bam.bai`
- `sam/pe/in_WT_G1_Hho1_6337.sp.sam`, `bam/pe/in_WT_G1_Hho1_6337.sp.bam`, `bam/pe/in_WT_G1_Hho1_6337.sp.bam.bai`

Treatment-name PE alignment variants:
- `bam/pe/IP_WT_G1_HU_Hho1_6336.sc.bam`, `bam/pe/IP_WT_G1_HU_Hho1_6336.sc.bam.bai`
- `bam/pe/in_WT_G1_HU_Hho1_6336.sc.bam`, `bam/pe/in_WT_G1_HU_Hho1_6336.sc.bam.bai`
- `bam/pe/IP_WT_G1_HU_Hho1_6337.sc.bam`, `bam/pe/IP_WT_G1_HU_Hho1_6337.sc.bam.bai`
- `bam/pe/in_WT_G1_HU_Hho1_6337.sc.bam`, `bam/pe/in_WT_G1_HU_Hho1_6337.sc.bam.bai`

The treatment-name variants copy the generated PE main-alignment fixtures under filenames that include a treatment token. They are used only to exercise treatment-aware metadata parsing through the submit and execute wrappers.

Reference FASTA:
- `reference/tiny.fa`
- `reference/tiny.fa.fai`

Reference-backed CRAM fixtures:
- `cram/se/*.cram`, `cram/se/*.cram.crai`
- `cram/pe/*.cram`, `cram/pe/*.cram.crai`

Spike-in scaling-factor part files:
- `parts/example_scaling_factors.spike.tsv.part.000000`
- `parts/example_scaling_factors.spike.tsv.part.000002`

siQ-ChIP scaling-factor part files:
- `parts/example_scaling_factors.siq.tsv.part.000000`
- `parts/example_scaling_factors.siq.tsv.part.000002`

Spike-in malformed negative-test part file:
- `parts/malformed_scaling_factors.spike.tsv.part.000003`
- `parts/header_scaling_factors.spike.tsv.part.000004`
- `parts/duplicate_index_A.spike.tsv.part.000005`
- `parts/duplicate_index_B.spike.tsv.part.000005`

siQ-ChIP metadata parser fixtures:
- `metadata/measurements_siqchip.tsv`
- `metadata/measurements_siqchip.tsv.gz`
- `metadata/measurements_siqchip_aliases.tsv`
- `metadata/measurements_siqchip_missing_required.tsv`
- `metadata/measurements_siqchip_unsupported_alias.tsv`
- `metadata/measurements_siqchip_skip_prefixes.tsv`
- `metadata/measurements_siqchip_alias_collision.tsv`
- `metadata/measurements_siqchip_duplicate_match.tsv`
- `metadata/measurements_siqchip_treatment.tsv`

siQ-ChIP metadata parser config fixtures:
- `config/parse_metadata_siqchip_regex.yml`
- `config/parse_metadata_siqchip_match_treatment.yml`

The nonconsecutive numeric suffixes exercise deterministic numeric ordering when input paths are supplied in reverse order. The malformed part contains too few fields and exercises combiner validation. The header-looking part and duplicate-index pair exercise additional combiner validation without creating those static inputs inside smoke tests.

<br />

## Alignment-backed spike-in expected behavior
| sample | main IP | spike IP | main input | spike input | `chiprx_alpha_ratio` |
| :----  | ------: | -------: | ---------: | ----------: | -------------------: |
| `A`    | 3       | 1        | 2          | 2           | `2`                  |
| `B`    | 2       | 2        | 3          | 1           | `0.5`                |

The serial execute-layer smoke test uses these generated alignment files directly. Its baseline case does not provide depth overrides, so the workflow counts alignments with Samtools and automatically detects the SE layout.

<br />

## Current and deferred smoke-test coverage
The smoke suite covers direct spike-in and siQ-ChIP part-file combination, numeric ordering, dry-run behavior, malformed input rejection, overwrite protection, explicit part-file cleanup, direct spike-in and siQ-ChIP Python calculator behavior, siQ-ChIP metadata parsing, serial submit/execute-layer siQ-ChIP behavior, and serial execute-layer spike-in assembly. Metadata parser coverage includes compressed tables, skip prefixes, alias normalization and collision rejection, ambiguous-match rejection, regex-mode parsing, and treatment-aware matching. Serial spike-in coverage exercises every canonical coefficient, accepted alias normalization, alignment-derived counts, automatic SE/PE detection, BAM/CRAM input, mixed input lists, and broadcast depth overrides.

The scaling-factor rows intentionally contain only core workflow values. Compute downstream denominator floors separately with `python -m scripts.compute_input_floor`.

`#TODO`: Add calculate-scaling-factor GNU Parallel and Slurm coverage.
