## [Unreleased]
- Placeholder for changes in progress.

<br />

## 2025-09-14
### Continued Development
- *Organized and committed all previously uncommitted code from the last four weeks.*
- Signal computation (Python)
    + `compute_signal.py`: Correct PE/SE fragment interval construction (leftmost 99/163; TLEN > 0), strand-aware SE 5’ extension using `query_alignment_length`, guards for `usr_frg <= 0` and rare `reference_end=None`, clamping to contig bounds, and explicit half-open [start,end) invariant. Expanded `parse_bam` docstring.
    + `compute_signal_ratio.py`: Hardened ratio pipeline&mdash;bin-size parity check, optional denominator clamp via `--dep_min`, consistent division-by-zero handling, optional `--log2` path with clear -inf/NaN policy, optional `--track` output (filters -inf/NaN), clearer CLI help and argument validation.
- Executors & submit wrappers (Shell)
    + `execute_compute_signal.sh`: Mode-aware command construction, safer dry-run (`eval "$(build_cmd)"`), bin-size validation, and cleanup of echo/serialization; default `interactive=false`.
    - `execute_calculate_scaling_factor.sh`: Pruned stale blocks and prepped for explicit mode handling.
    - `execute_filter_bams.sh`, `submit_filter_bams.sh`: Standardized logging, CSV-style stdout for easier parsing, aligned flag echoes.
    - `submit_compute_signal.sh`, `submit_calculate_scaling_factor.sh`: Debug to stderr; parse `set_logs_slurm` via CSV (`IFS=','`); small usage/docs tidy-ups.
- Shared Shell functions
    - `scripts/functions/submit.sh`: `debug_var` prints to stderr; `set_logs_slurm` uses hard-linked logs and clearer errors.
    - `scripts/functions/calculate_scaling_factor.sh`: Safer `bc` inputs (leading zeros), pass `--verbose` through, dynamic field assembly, tabbed output; mirrored for spike-in path.
    - `scripts/functions/populate_array_empty.sh`: New nounset-safe helper&mdash;creates target if missing, validates integers, parameterized fill value, avoids unbound vars.
    - `scripts/functions/check_arrays_lengths.sh`: Generalized array length checker using indirect lookups; fixed second-array size bug; clearer failures and non-zero exits.
- Headers and outputs
    - `write_header.sh`: Idempotent header writing (prepend only if absent), validate output directory (not file), create file when missing; in spike mode, rename header column from “sf” to “spike.”
- Documentation
   - `README.md`, `workflow.md`, `validate_siq_chip.md`, `download_process_fasta_gff3.md`, `CHANGELOG.md`: Refreshed quickstart/requirements and terminology, clarified FASTA/GFF3 prep and sanity checks, added validation tips (track integrals, bin width, IGV), and tidied `CHANGELOG` wording.
- *The repository remains under active, heavy development. Ongoing refactors and test coverage are in progress, and some interfaces will continue to evolve as we stabilize the workflow.*

<br />

## 2025-08-18
### Returning to Development
- Organized and committed all previously uncommitted code after a break from active work, bringing the repository back into a clean, trackable state.
- Refactored multiple Bash driver and submission scripts for consistency and reduced redundancy (to be continued).
- Updated core functions for robustness and naming clarity (e.g., renamed `extract_field_str.sh`, `summarize_sig_norm.sh`).
- Added new utility functions (`format_print_cmd.sh`, `make_lower.sh`) for string handling and command formatting.
- Refined signal computation scripts (`compute_signal.py`, `compute_signal_ratio.py`) with in-progress integration testing.
- Updated documentation (`workflow.md`) to reflect the current flow.
- Corrected help messages and replaced the case variable `type` with `mode` in write_header.sh.
- Began to resume active work following publication in *Bio-protocol*. Note that some features are still in progress, not fully tested, or likely to change.
