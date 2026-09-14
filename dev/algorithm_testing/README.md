
# Algorithm-testing quarantine for compute_signal
This directory is the temporary quarantine home for the experimental signal machinery removed from `src/protocol_chipseq_signal_norm/cli/compute_signal.py` when the production file was reduced to its single public arithmetic (`direct_sparse_np`, private-roadmap shorthand "B1"). A future archive repository is expected to adopt this directory wholesale; until then it lives here so the July 2026 benchmark and correctness comparisons stay runnable from source.

`compute_signal_algotest.py` carries, unchanged: the full-featured sparse accumulator with its `B3` (`use_bincount`), `touched`, `local`, and `idx` variants; the reference dict kernel and the array kernel; the dense and event accumulators; the seven-strategy merge; and the digest/profile-only write helpers. The S3-MIG-009 masking fix is retained at every accumulation site. The module imports shared helpers from the installed package; production code must never import from it.

`tests/` holds the relocated test module. It is deliberately outside the maintained safe suite (`tests/run_tests.sh` never discovers it); run it directly:
```bash
python -m pytest dev/algorithm_testing/tests
```
