# Changelog

## [Unreleased]
- Placeholder for changes in progress.

<br />

## [v0.2.0-rc.1] - 2026-05-02
### Release Status
- Provisional pre-release milestone for the large helper, wrapper, and Python utility refactor.
- This version is still pre-testing and should not be treated as a fully validated stable release.
- Needs review: smoke/integration testing coverage, remaining help-doc consistency passes, and manuscript/blog support material should be reviewed before promoting this milestone to a final release.

### Helper Consolidation
- Added shared Bash helper-loading infrastructure with source-once behavior, source-only guards, and centralized helper dependency sourcing.
- Added shared command-printing, Python-runner, and wrapper log helper utilities for refactored shell wrapper chains.
- Normalized foundational helper libraries for argument checks, numeric validation, input validation, environment handling, output formatting, parallel execution, Slurm log handling, and array population.
- Continued consolidating standalone helper scripts into broader helper modules, including argument, numeric, input, environment, output-formatting, sequence, table, region, BAM-filter, parallelization, and interactive-exit helpers.
- Renamed `scripts/functions/submit.sh` to `scripts/functions/manage_slurm.sh` to clarify its Slurm-specific role.

### Domain Helper Refactors
- Refactored alignment helpers to support current option naming, BWA/BWA-MEM2-related handling, BAM/CRAM output paths, reference FASTA validation, and consistent error reporting.
- Refactored filtering helpers for shared parsing, validation, chromosome checks, and organism-specific BAM filtering paths.
- Refactored region, sequence, table, and find-command helpers around current naming conventions, stricter validation, and safer command or array construction.
- Expanded scaling-factor helper logic for siQ-ChIP and spike-in workflows, including metadata parsing, alignment-count handling, depth-factor handling, per-sample output parts, and Python runner integration.

### Wrapper API and Wiring
- Refactored active `execute_*` and `submit_*` wrapper chains around shared helper sourcing, `require_optarg`-based parsing, command-array construction, and consistent error/help behavior.
- Aligned FASTQ download, trimming, alignment, BAM filtering, compute-signal, and calculate-scaling-factor wrapper paths across Slurm, GNU Parallel, serial, dry-run, and debug execution modes.
- Propagated newer alignment, filtering, signal, ratio, metadata, and scaling-factor options through execute and submit layers.
- Standardized canonical option naming across wrapper APIs, including `--dp` as the preferred decimal-precision option while preserving compatibility aliases such as `--rnd`, `--round`, `--decimals`, and `--digits`.
- Preserved compatibility aliases where practical while moving wrapper call sites toward clearer long-option forms.
- Added `scripts/symlink_files.sh` as a utility wrapper following the newer serialized input/output and dry-run conventions.

### Python Signal and Scaling Utilities
- Added shared Python utility modules for CLI formatting, interactive argument handling, I/O, validation, chromosome sorting, bedGraph parsing, and robust stabilizer selection.
- Modernized signal, ratio, and input-floor utilities around shared Python helpers, stronger CLI validation, and clearer bedGraph handling.
- Added standalone utilities for pseudocount estimation, bedGraph interval merging, bedGraph summation, and namespaced coefficient-table augmentation.
- Expanded spike-in scaling-factor tooling with additional coefficient outputs, output formats, rounding aliases, and documentation of related spike-in scaling frameworks.
- Modernized siQ-ChIP scaling-factor and metadata parsing utilities with configurable metadata handling, shell-safe output, stricter validation, and updated CLI behavior.

### Help Documentation
- Extracted and normalized help documentation for active wrappers and utility entrypoints into sourced help modules.
- Standardized usage blocks, option descriptions, examples, notes, dependency sections, and compatibility-alias documentation across the newer shell entrypoints.
- Updated environment-installation help to describe dry-run behavior, Conda/Mamba fallback behavior, current package lists, and hidden legacy environments.

### Setup, Documentation, and Repository Cleanup
- Added and refined a POSIX-compatible `install_envs` entrypoint that provides setup guidance before Bash >= 5 is available.
- Updated `install_envs.sh` to use shared helpers, support dry-run mode, fall back from Mamba to Conda, and refresh environment package lists.
- Updated README setup guidance for the install-env entrypoint and clarified workflow documentation language.
- Clarified genome-processing documentation around `bowtie2` alignment and index use.
- Expanded `.gitignore` coverage for generated analysis files, scratch data, compressed outputs, PDFs, caches, tests, and temporary directories.
- Removed the root `__init__.py` package marker while retaining package markers under `scripts/` and `scripts/functions/` for module-style imports.

### Known Follow-Ups / Not Yet Validated
- Needs review: run smoke and integration tests across Slurm, GNU Parallel, serial, dry-run, and debug execution paths.
- Needs review: add or formalize unit, regression, and fixture-based tests for key Bash helper functions, wrapper command construction, and Python signal/ratio/scaling utilities.
- Needs review: continue validating wrapper API consistency, help-text dependency claims, canonical aliases, and hidden/transitional aliases after the helper-layer refactor.
- Needs review: normalize remaining help dependency headings and audit stale helper names such as legacy `submit.sh`, `echo_error`, and `echo_warning` references.
- Needs review: decide whether legacy one-off Slurm scripts should remain explicitly legacy or be migrated to the current helper/parser/help conventions.
- Needs review: defer Bash 5-specific cleanup until after correctness and smoke testing are complete.
- Needs review: decide how manuscript/blog support material should be tracked and documented relative to production workflow scripts.
- Needs review: continue testing Python signal, ratio, pseudocount, metadata, and scaling-factor utilities against representative data.

<br />

## [v0.1.1] - 2025-09-14
### Release Status
- Continued-development cleanup milestone after several weeks of uncommitted work.
- The repository remained under active development; interfaces and tests were still being stabilized.

### Python Signal and Ratio Updates
- Corrected PE/SE fragment interval construction in `compute_signal.py`, including strand-aware SE 5' extension, TLEN handling, contig-bound clamping, and half-open interval invariants.
- Hardened `compute_signal_ratio.py` with bin-size parity checks, optional denominator clamping via `--dep_min`, clearer division-by-zero handling, optional `--log2` behavior, optional cleaned track output, and improved argument validation.

### Wrapper API and Wiring
- Updated `execute_compute_signal.sh` with mode-aware command construction, safer dry-run behavior, bin-size validation, and cleaner serialization handling.
- Pruned stale logic in `execute_calculate_scaling_factor.sh` in preparation for more explicit mode handling.
- Standardized logging and CSV-style parsing behavior in filtering and scaling submit wrappers.
- Routed debug output to stderr in compute-signal and calculate-scaling-factor submit wrappers.

### Shared Shell Helpers
- Improved Slurm log handling in the then-current submit helper, including stderr debug output and hard-linked descriptive log files.
- Improved scaling-factor helper behavior around safer `bc` inputs, verbose propagation, dynamic field assembly, and tabular output.
- Added nounset-safe array population behavior and generalized array-length checking.

### Headers, Outputs, and Documentation
- Made `write_header.sh` idempotent, validated output directories more carefully, and adjusted spike-mode header naming.
- Refreshed README, workflow, validation, genome-processing, and changelog documentation for the active development state.

<br />

## [v0.1.0] - 2025-08-18
### Release Status
- Baseline return-to-development milestone after publication-related downtime.
- The repository was brought back into a cleaner, trackable state, but features remained in progress and not fully validated.

### Wrapper and Function Refactors
- Refactored multiple Bash driver and submission scripts for consistency and reduced redundancy.
- Updated core function scripts for robustness and naming clarity.
- Added utility helpers for string formatting and command formatting.
- Corrected help messages and replaced the `type` case variable with `mode` in `write_header.sh`.

### Python Signal Work
- Refined `compute_signal.py` and `compute_signal_ratio.py` as part of ongoing integration testing.

### Documentation
- Updated workflow documentation to reflect the then-current pipeline flow.
- Recorded that continued refactoring and test coverage were still in progress.
