
# bedGraph unit-test fixtures
These fixtures are the small bedGraph, BED, and chromosome-size inputs that the Python unit tests of the bedGraph and signal code read. Each file is named for the case it serves and sits in a directory named for its consumer. They need no managed environment.

Regenerate fixtures from the repository root with:
```bash
bash tests/fixtures/bedgraph/make.sh
```

Generated fixture outputs are ignored by Git. `tests/run_tests.sh` regenerates this fixture set automatically when required inputs are missing.

<br />

## Files
Readable provenance:
- `make.sh`, which writes every file literally, one argument per line

Inputs by consumer:
- `bdg/`: `utils_bdg`, in `tests/unit/utilities/test_bdg.py`: bin-size and grid checks, chromosome-size tables, bounds, and bin counting
- `input_floor/`: `compute_input_floor`, in `tests/unit/compute_signal/test_input_floor.py` and `test_input_floor_opt_applicability.py`
- `merge_bins/`: `merge_bins`, in `tests/unit/bedgraph/test_merge_bins.py`, space-delimited as the script accepts
- `pseudo/`: `compute_pseudo`, in `tests/unit/compute_signal/test_pseudo.py`
- `ratio/`: `comp_sig_rat`, in `tests/unit/compute_signal/test_ratio.py`, as `A` and `B` pairs, space-delimited
- `signal/`: a chromosome-size table, in `tests/unit/compute_signal/test_signal.py`
- `sum/`: `sum_bdg`, in `tests/unit/bedgraph/test_sum.py`, space-delimited

<br />

## Current and deferred test coverage
Current coverage: each file serves the unit test named in its case, and the tests read the files in place.

Deferred:
- `test_pseudo.py` still derives normalized copies of the `compute_pseudo` pair at run time (`_norm_like`).
