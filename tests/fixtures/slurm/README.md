
# Slurm wet-validation fixtures
These fixtures are the dedicated inputs for `tests/integration/slurm/run_wet_tests.sh`. They are intentionally tiny and independent of the workflow fixtures elsewhere under `tests/fixtures/`.

They are intentionally small and hand-checkable. Running `make.sh` regenerates the fixture set deterministically.

Regenerate fixtures from the repository root with:
```bash
bash tests/fixtures/slurm/make.sh
```

Generated fixture outputs are ignored by Git. `tests/run_tests.sh` regenerates this fixture set when the Slurm gate selects the wet suite, and not otherwise, because preparing fixtures for a disabled optional group is what `TEST.GATES.OPTIONAL` forbids.

These are workflow inputs rather than checker inputs, so the directories name the kind of data they hold, matching `compute_signal` and the other workflow fixture sets. No verdict directory would mean anything here: a reference sequence is not accepted or rejected by a rule.

No fixture preparation writes to the extracted source tree or to the ordinary remote repository checkout.

<br />

## Files
Readable provenance:
- `make.sh`

Generated workflow inputs:
- `reference/tiny.fa`: one 108-base sequence on contig `I`
- `fastq/tiny_se.fastq`: one 30-base single-end read
- `sam/tiny_signal.sam`: two 10-base single-end alignments, one forward and one reverse, against the 108-base reference

<br />

## Expected wet-validation behavior
The remote runner builds a Bowtie 2 index from `reference/tiny.fa` and converts `sam/tiny_signal.sam` into a sorted, indexed BAM beneath the isolated run's `results/artifacts/fixtures/` directory before submitting jobs. The bounded wet suite then submits exactly two single-CPU jobs.

The SAM rows are written through `write_sam_line` rather than typed into a heredoc, so no tab is committed as an invisible control character that an editor could convert to spaces. A converted tab would be accepted by the recipe and rejected by Samtools on the remote host, which is the worst place to discover it.

<br />

## Current and deferred test coverage
Current coverage in `tests/integration/slurm/`:
- the coordinator includes this fixture root in the bundled runtime closure and records a digest for every file it copies; and
- the gated wet suite consumes the reference, read, and alignments in its two submitted jobs.

Deferred:
- paired-end wet coverage, which needs a second read set and a scheduler allocation that the bounded suite deliberately avoids.
