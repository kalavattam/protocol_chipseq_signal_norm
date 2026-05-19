# filter_bams fixtures

These fixtures exercise `submit_filter_bams.sh` and `execute_filter_bams.sh`
with tiny deterministic alignment inputs.

Regenerate them from the repository root with:

```bash
bash tests/filter_bams/scripts/make_fixtures.sh
```

## Files

- `sam/filter_sc_sp.sam`
  - Coordinate-sorted SAM with one read each on `I`, `Mito`, `SP_I`,
    `SP_II_TG`, `SP_MTR`, `SP_Mito`, and `chrUn`.
  - Used by BAM smoke tests to build BAM input at test runtime.
  - Used by CRAM smoke tests to build CRAM input at test runtime.

- `reference/filter_sc_sp.fa`
  - Tiny reference FASTA with all contigs present in the SAM header.
  - Every contig is 100 bp.
  - Required for deterministic CRAM input generation, CRAM reading, and CRAM
    output writing.

- `reference/filter_sc_sp.fa.fai`
  - Deterministic FASTA index for `reference/filter_sc_sp.fa`.
  - Generated directly by this script because the reference uses fixed
    one-line contigs.

## Expected smoke-test behavior

- `retain=sc` keeps chromosome `I` and drops `Mito` unless `--mito` is used.
- `retain=sc` drops S. pombe contigs such as `SP_I`.
- `retain=sp --tg --mtr --mito` keeps `SP_I`, `SP_II_TG`, `SP_MTR`, and
  `SP_Mito`.
- `retain=sp` drops S. cerevisiae contigs such as `I` and unrelated contigs
  such as `chrUn`.
- CRAM input tests must pass `--ref_fa reference/filter_sc_sp.fa` and still
  produce BAM output.
- CRAM output tests must pass `--out_ext cram --ref_fa reference/filter_sc_sp.fa`
  and produce indexed CRAM output.
