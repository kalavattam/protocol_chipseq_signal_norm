
# Compute-pseudocount test fixtures
These fixtures are synthetic micro-fixtures for fast, deterministic tests of `compute_pseudo`, whose edgeR pseudocount arithmetic is owned by `utilities.utils_stabilizer.compute_pseudo_edger`.

They are intentionally small and hand-checkable. Running `make.sh` regenerates the fixture set deterministically.

Regenerate fixtures from the repository root with:
```bash
bash tests/fixtures/compute_pseudo/make.sh
```

Generated fixture outputs are ignored by Git. `tests/run_tests.sh` regenerates this fixture set automatically when required inputs are missing.

The directory names the role each file plays in the assertion. These are bedGraph inputs, which is the role a workflow input plays, rather than a checker verdict, so the root has no `accepted/` or `rejected/` split.

Only the CLI surface needs these files. The estimator itself takes library sizes as floats and is tested directly in `tests/unit/utilities/test_stabilizer.py` with no input on disk at all; a fixture would add a file read to arithmetic that has no file in it.

<br />

## Files
Readable provenance:
- `make.sh`

Generated bedGraph inputs:
- `bedgraph/pair_A.bdg`: three rows on a uniform 10 bp grid, summing to a library size of 6
- `bedgraph/pair_B.bdg`: the same grid, summing to a library size of 18

<br />

## Expected estimator behavior
The pair is deliberately imbalanced 1:3 rather than near-equal, so a test asserting the per-sample depth correction cannot pass by returning the nominal `prior.count` for both tracks.

Every quantity the estimator derives is a round number:

| quantity | value |
| :--- | ---: |
| `L_A` | 6 |
| `L_B` | 18 |
| `L_bar` | 12 |
| `prior_scaled_A` = `prior.count * L_A / L_bar` | 1.0 |
| `prior_scaled_B` = `prior.count * L_B / L_bar` | 3.0 |

The two priors sum to `2 * prior.count` and their ratio is `L_A / L_B`, which are the invariants the consuming tests assert. A uniform 10 bp grid also lets the bin width be inferred rather than supplied, so the same pair exercises `--siz_bin` omission.

<br />

## Current and deferred test coverage
Current coverage in `tests/unit/compute_signal/test_pseudo.py`:
- the `--prt_jsn` payload reports `prior_scaled`, and its ratio and sum hold; and
- under `--normalization nc` the prior tracks the fragment counts rather than the library sizes, and `pseudo_i / scale_i` does not recover it.

Deferred:
- a track carrying `track` or `browser` header lines, which `--skp_pfx` handling needs and which the inline constructions in `test_pseudo.py` still cover; and
- a nonuniform bin grid, which belongs with `sum_counts_bdg` in `tests/unit/utilities/test_bdg.py` rather than with the pseudocount CLI.
