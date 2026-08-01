
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
