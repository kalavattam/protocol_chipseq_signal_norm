
# Relocated algorithm-testing tests
`algotest_variants_test.py` holds the tests that moved out of
`tests/unit/compute_signal/` together with the experimental machinery they
prove: the dict-kernel and array-kernel assertions, the variant-format
equivalences against the public `B1` baseline, and the `B3`/array arms of the
non-negativity property. The file uses the `*_test.py` suffix deliberately:
repository contracts reserve the `test_*.py` prefix for the maintained suites
under `tests/`, and this suite is not part of them.

Run from the repository root with the project environment active:
```bash
python -m pytest dev/algorithm_testing/tests
```

The maintained safe suite must stay green without these tests.
