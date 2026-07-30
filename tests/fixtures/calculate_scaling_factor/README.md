
# Calculate-scaling-factor test fixtures
These fixtures are realistic micro-fixtures for fast, deterministic tests of the calculate-scaling-factor workflow.

They are intentionally small and hand-checkable. Running `make.sh` from `env_protocol` regenerates the fixture set deterministically.

Regenerate fixtures from the repository root with:
```bash
source activate env_protocol
bash tests/fixtures/calculate_scaling_factor/make.sh
```

Generated fixture outputs are ignored by Git. When required inputs are missing, `tests/run_tests.sh` regenerates this fixture set automatically.

The fixture rows are adapted from tracked example scaling-factor tables under `data/processed/compute_signal/`. They preserve realistic file naming, core table shapes, scaling values, and alignment depths. Synthetic role-specific SAM, BAM, and CRAM fixtures exercise Samtools-backed alignment counting.

Alignment fixture basenames are biologically minded because siQ-ChIP tests exercise production filename parsing. PE fixtures use Hho1-like names, and SE fixtures use Brn1-like names.

Spike-in rows also record the canonical coefficient method used alongside the computed value. The default fixture examples use `chiprx_alpha_ratio`.

Minimal siQ-ChIP metadata tables exercise `parse_metadata_siqchip.py` with the deterministic production parsing YAML in `data/raw/docs/parse_metadata_siqchip.yml`. These tables are intentionally small and cover
- exact filename-token parsing,
- configured `matching.fields`,
- explicit `field_to_column` mapping,
- gzip-compressed metadata input,
- duplicate-match rejection,
- missing required calculator inputs,
- optional precomputed `dep_*` and `len_*` values, and
- optional paired library-loading volume correction.

The parser fixture config `config/parse_metadata_siqchip_field_to_column.yml` proves that filename field names do not need to match metadata table column names when an explicit mapping is supplied.

The production reference config `data/raw/docs/parse_metadata_siqchip_PRJNA857063.yml` and table `data/raw/docs/measurements_siqchip_PRJNA857063.tsv` model factor-tagged logical input names from the human siQ-ChIP paper metadata for future reanalysis readiness without downloading reads.

<br />

## Files
Readable provenance:
- `make.sh`

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

The treatment-name variants copy the generated PE main-alignment fixtures under filenames that include a treatment token. They remain available as alignment fixtures, but the deterministic parser tests now use the configured `matching.fields` contract rather than a treatment-specific parser mode.

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

siQ-ChIP deterministic metadata parser fixtures:
- `metadata/measurements_siqchip.tsv`
- `metadata/measurements_siqchip.tsv.gz`
- `metadata/measurements_siqchip_lib_volume.tsv`
- `metadata/measurements_siqchip_lib_volume_one_sided.tsv`
- `metadata/measurements_siqchip_lib_volume_zero.tsv`
- `metadata/measurements_siqchip_missing_required.tsv`
- `metadata/measurements_siqchip_duplicate_match.tsv`
- `metadata/measurements_siqchip_precomputed.tsv`

siQ-ChIP metadata parser config fixtures:
- `config/parse_metadata_siqchip_field_to_column.yml`

The siQ metadata fixtures cover
- the production deterministic parser config,
- explicit field-to-column mapping,
- gzip input,
- duplicate matches,
- optional library-loading volume correction,
- precomputed depth and fragment-length values, and
- malformed metadata tables.

Library-volume fixtures exercise
- the valid paired correction path and
- malformed one-sided, zero, and nonnumeric values.

The nonconsecutive numeric suffixes exercise deterministic numeric ordering when input paths are supplied in reverse order.

The malformed part contains too few fields and exercises combiner validation.

The header-looking part and duplicate-index pair exercise additional combiner validation without creating those static inputs inside tests.

siQ-ChIP part files may include optional `lib_vol_ip` and `lib_vol_in` fields. Output headers include those columns only when at least one combined part row supplies the paired library-loading volume correction.

<br />

## Alignment-backed spike-in expected behavior
| sample | main IP | spike IP | main input | spike input | `chiprx_alpha_ratio` |
| :---   | ---:    | ---:     | ---:       | ---:        | ---:                 |
| `A`    | 3       | 1        | 2          | 2           | `2`                  |
| `B`    | 2       | 2        | 3          | 1           | `0.5`                |

The serial execute-layer test uses these generated alignment files directly. Its baseline case does not provide depth overrides, so the workflow counts alignments with Samtools and automatically detects the SE layout.

<br />

## Current and deferred test coverage
The test suite covers
- direct spike-in and siQ-ChIP part-file combination,
- numeric ordering,
- dry-run behavior,
- malformed input rejection,
- overwrite protection,
- explicit part-file cleanup,
- direct spike-in and siQ-ChIP Python calculator behavior,
- siQ-ChIP metadata parsing,
- serial submit/execute wrapper behavior for both scaling-factor branches, and
- gated local GNU Parallel wet execution.

Metadata parser coverage includes
- compressed tables,
- deterministic filename token parsing,
- configured matching fields,
- explicit field-to-column mapping,
- ambiguous-match rejection,
- zero-match rejection,
- missing required metadata rejection,
- precomputed `dep_*` and `len_*` values, and
- paired library-volume validation.

Wrapper coverage exercises
- expected failures,
- CRAM/reference behavior,
- mixed SE/PE layouts,
- mixed BAM/CRAM inputs,
- siQ equation validation,
- spike-in method validation,
- alignment-derived counts,
- automatic SE/PE detection,
- every canonical spike-in coefficient, and
- broadcast depth overrides.

The scaling-factor rows intentionally contain only core workflow values. Compute downstream denominator floors separately with `python -m protocol_chipseq_signal_norm.cli.compute_input_floor`.

`#TODO`: Add calculate-scaling-factor Slurm coverage.
