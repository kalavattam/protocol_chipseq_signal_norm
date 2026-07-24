
# Slurm wet-validation fixtures
These text-only fixtures are the dedicated inputs for `tests/integration/slurm/run_wet_tests.sh`. They are intentionally tiny and independent of the generated, Git-ignored workflow fixtures elsewhere under `tests/`.

The remote runner builds a Bowtie2 index from `reference/tiny.fa` and converts `sam/tiny_signal.sam` to a sorted, indexed BAM beneath the isolated run's `results/artifacts/fixtures/` directory before submitting jobs. The source FASTQ contains one 30-base read, and the signal SAM contains two 10-base alignments on the 108-base reference.

No fixture preparation writes to the extracted source tree or the normal remote repository checkout.
