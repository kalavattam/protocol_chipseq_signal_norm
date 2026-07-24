
# Filter-alignments test fixtures
These fixtures exercise `submit_filter_alignments.sh` and `execute_filter_alignments.sh` with tiny deterministic alignment inputs.

They are intentionally small and hand-checkable. Running `make.sh` regenerates the fixture set deterministically.

Regenerate them from the repository root with:
```bash
bash tests/fixtures/filter_alignments/make.sh
```

Generated fixture outputs are ignored by Git. When required inputs are missing, `tests/run_tests.sh` regenerates this fixture set automatically.

<br />

## Files
- `sam/filter_sc_sp.sam`
  + Coordinate-sorted SAM with one read each on `I`, `Mito`, `SP_I`, `SP_II_TG`, `SP_MTR`, `SP_Mito`, and `chrUn`.
  + Used by BAM tests to build BAM input at test runtime.
  + Used by CRAM tests to build CRAM input at test runtime.
- `reference/filter_sc_sp.fa`
  + Tiny reference FASTA with all contigs present in the SAM header.
  + Every contig is 100 bp.
  + Required for deterministic CRAM input generation, CRAM reading, and CRAM output writing.
- `reference/filter_sc_sp.fa.fai`
  + Deterministic FASTA index for `reference/filter_sc_sp.fa`.
  + Generated directly by this script because the reference uses fixed one-line contigs.

<br />

## Expected test behavior
- `retain=sc` keeps chromosome `I` and drops `Mito` unless `--mito` is used.
- `retain=sc` drops *S. pombe* contigs such as `SP_I`.
- `retain=sp --tg --mtr --mito` keeps `SP_I`, `SP_II_TG`, `SP_MTR`, and `SP_Mito`.
- `retain=sp` drops *S. cerevisiae* contigs such as `I` and unrelated contigs such as `chrUn`.
- CRAM input tests must pass `--ref_fa reference/filter_sc_sp.fa` and still produce BAM output.
- CRAM output tests must pass `--out_ext cram --ref_fa reference/filter_sc_sp.fa` and produce indexed CRAM output.

<br />

## Current and deferred test coverage
The test suite covers submit- and execute-layer BAM-to-BAM, CRAM-to-BAM, BAM-to-CRAM, and CRAM-to-CRAM paths. It also covers CRAM-to-CRAM GNU Parallel dispatch for `retain=sc` and `retain=sp` when gated by `RUN_PARALLEL=1`.

`#TODO`: Add Slurm coverage. Consider extending GNU Parallel coverage to additional input/output format combinations where that broader matrix provides useful regression protection.
