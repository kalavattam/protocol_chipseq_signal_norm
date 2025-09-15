## [Unreleased]
- Placeholder for changes in progress.

## [2025-08-18]
### Returning to Development
- Organized and committed all previously uncommitted code after a break from active work, bringing the repository back into a clean, trackable state.
- Refactored multiple Bash driver and submission scripts for consistency and reduced redundancy (to be continued).
- Updated core functions for robustness and naming clarity (e.g., renamed `extract_field_str.sh`, `summarize_sig_norm.sh`).
- Added new utility functions (`format_print_cmd.sh`, `make_lower.sh`) for string handling and command formatting.
- Refined signal computation scripts (`compute_signal.py`, `compute_signal_ratio.py`) with in-progress integration testing.
- Updated documentation (`workflow.md`) to reflect the current flow.
- Corrected help messages and replaced the case variable `type` with `mode` in write_header.sh.
- Began to resume active work following publication in *Bio-protocol*. Note that some features are still in progress, not fully tested, or likely to change.