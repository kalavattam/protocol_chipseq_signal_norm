
# Output standard
This document owns shared output contracts that are not scientific computation rules.

<br />

## Numeric rendering (`OUTPUT.NUMERIC.RENDERING`)
**Classification:** `semantic-only` with an approved post-round invariant and deterministic evidence inventory.

**Scope:** All inventoried maintained Python user-facing numeric emission sites with an interface-owned or inherited `dp` are presumed applicable. A site may be excluded only through explicit reviewed classification as not `dp`-governed or as an owned protocol exception.

Only after explicit site enrollment, an at-most contract retains no more than `dp` decimal places, removes non-informative trailing zeros and a trailing decimal point, normalizes emitted negative zero, keeps finite and nonfinite handling distinct, and records approved protocol exceptions. These requirements apply to an already selected post-round value and do not choose how rounding occurs.

No choice is made here for rounding, binary-floating semantics, ties, `dp=0`, very small values, nonfinite spelling or rejection, a Python implementation, or a Shell mechanism. The separately authorized Python migration must resolve those choices before `0.3.0`; Shell implementation remains deferred with inventory and intended parity preserved.

**Automation:** `dev/audit/numeric_emission_inventory.py` completely enumerates only its versioned, recognized static Python and Shell sink classes. `dev/audit/numeric_output_applicability.py` separately validates stable applicability references and migration readiness. The inventory records dynamic, indirect, generated, sourced, eval-like, embedded-language, dynamic-stream, unrecognized-writer, and otherwise unsupported constructs as explicit review candidates. The result is not a complete semantic or runtime emission-site inventory. It changes no output. Coverage is `independent evidence`.

**Semantic remainder:** Global `S3-MIG-002` and pre-0.3.0 numeric-migration readiness may close only when every inventoried maintained Python user-facing numeric candidate has exactly one terminal disposition: applicable `dp`-governed and migrated; explicitly excepted by an authoritative owner with rationale; or proved not to be `dp`-governed. Any unresolved candidate blocks global closure. A selected bounded cohort may be ready when every member has a valid terminal disposition, but never implies global readiness while unresolved candidates remain elsewhere.

**Exceptions:** Fixed-width machine protocols and significant-digit diagnostics are candidates for review, not implicit exceptions.

<br />

## Alignment-file header record (`OUTPUT.ALIGNMENT.RECORD`)
**Classification:** `advisory` with deterministic fixture evidence.

**Scope:** Every SAM, BAM, or CRAM file that a maintained chain writes as an output, including additional outputs such as the `--qname` BAM, in every chain that writes alignment files (today `align_fastqs` and `filter_alignments`).

A reader must be able to learn from the header alone everything the pipeline did to the file, whether the settings that shaped it are public or hidden. Keep every tool's own `@PG` line: never pass `--no-PG` when writing an output, and never strip `@PG` lines from one.

Each chain also adds a pipeline record to every output: `@CO` lines, one fact per line, in two groups.

- Content lines, each starting `pipeline:`, hold the facts that determine what the file contains: the chain; the pipeline commit, marked when the working tree has uncommitted changes; every option in effect that can change the file's content, with its value, defaults and hidden settings included; every step the chain skipped; every setting that changes what the file means, such as a hidden measurement setting whose output is unfiltered; and the inputs, named without their directories. These lines hold nothing that varies between runs that make the same file, so the same inputs, options, and commit give the same content lines.
- Run lines, each starting `pipeline-run:`, hold the context of the run: thread counts and other resource options, the start time, the host name, and the full input and output paths.

A chain that reads an alignment file with a record keeps it and appends its own, so the records accumulate in the order the chains ran.

The record survives every later step and conversion within the chain, including sorting, duplicate marking, and conversion to CRAM. A chain step that rewrites the header carries the record forward.

**Automation:** No dedicated checker, and no execution-registry entry is dedicated to this rule ID. Each chain's fixture tests read the headers of all its outputs, with `samtools view -H --no-PG` so that reading adds no `@PG` line, prove that both groups of the record are present and complete, and fail when a step drops it. A test running one chain twice with different thread counts proves that the content lines are identical and the run lines differ.

**Semantic remainder:** Decide which options can change a file's content and so belong in the content lines, and how an option's value is spelled in the record. Whether downstream tools such as `compute_signal` read the record, for example to refuse a file made with a hidden measurement setting, is a separate decision.

**Exceptions:** A tool's own `@PG` lines keep what the tool writes, including paths and thread counts. Header comparisons across runs, such as determinism checks, therefore compare the content lines and the `@PG` program sequence, not the run lines or the tools' command lines.
