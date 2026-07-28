
# Standards governance
This document defines how repository standards are authored, implemented, verified, changed, and excepted. The maintained standards listed in [`README.md`](README.md) are the normative source; prose summaries, machine configuration, checkers, evidence, fixtures, contracts, and implementation must remain traceable to them.

<br />

## Requirement language and terms (`GOV.TERMS`)
**Classification:** `advisory` with deterministic controlled-vocabulary portions.

**Scope:** Maintained standards and derived artifacts that use repository requirement, ownership, coverage, and exception terminology.

The keywords **must**, **must not**, **should**, **should not**, and **may** express requirements, recommendations, and permission. A deviation from **must** or **must not** requires an approved exception. A deviation from **should** or **should not** requires a documented reason when it is material to review.

- **Normative owner:** the single authoritative standards section that defines a rule, its classification, and its scope.
- **Maintained standard:** a normative document listed in [`README.md`](README.md). A backup, generated copy, proposal, or other unlisted document is not a maintained standard.
- **Derived artifact:** registry data, checker configuration or implementation, an evidence producer, fixture, test, repository contract, suppression, or generated reproduction that represents or applies a normative rule.
- **Checker:** automation that decides conformance for a stated deterministic scope without human judgment.
- **Evidence producer:** automation that collects bounded facts or review candidates without making the final semantic decision.
- **Semantic remainder:** the part of a rule whose final application still requires contextual human judgment after automation has run.
- **Repository contract:** an executable assertion that verifies a rule's checker or review procedure is connected to the repository's required workflow with the registered configuration.
- **Standards owner:** the person or review group authorized to approve a maintained standard and its changes. A standard or maintained ownership record may name a more specific owner; otherwise the repository maintainer identified in the top-level `README.md` is the standards owner.
- **Exception owner:** the person or review group accountable for the rationale, narrow scope, and timely review of an approved exception. The applicable standards owner approves the exception.
- **Suppression:** the narrow technical mechanism that prevents or changes a finding for an approved exception; a suppression is not itself approval.

**Automation:** No dedicated registry entry or checker owns the controlled vocabulary. The standards registry and domain checkers consume these terms but do not prove that prose applies them with the intended force.

**Semantic remainder:** Determine whether wording creates a requirement, recommendation, permission, ownership assignment, coverage claim, or exception obligation in context.

**Exceptions:** None alter the meanings defined here; deviations from a governed rule use that rule's approved exception mechanism.

<br />

## Authoritative ownership (`GOV.OWNER.UNIQUE`)
**Classification:** `advisory` with deterministic identity and registry portions.

**Scope:** Normative rules, rule IDs, owner sections, maintained standards, and derived artifacts that represent or apply those rules.

Every normative rule has one stable rule ID and one authoritative section in one maintained standard. Another document may link to the owner or provide a clearly marked non-normative summary or generated reproduction, but it must not present an independent normative version.

Rule IDs use stable uppercase dot-separated names scoped by domain, such as `PY.CLI.MAIN`, `MD.FENCE.CLOSED`, and `LAYOUT.PYTHON.SRC`. Renaming, retiring, merging, or splitting an ID is a standards change and requires an explicit migration of owner anchors, registry entries, checker output, fixtures or equivalent evidence, contracts, suppressions, and historical references that remain active.

The authoritative rule inventory is the set of stable rule IDs owned by the maintained standards listed in [`README.md`](README.md). Coverage tooling must enumerate that inventory from the owner documents. `dev/config/rules.toml` is the execution registry: it registers deterministic and advisory coverage portions and must reference their normative owner, but it does not create rules. A semantic-only rule with no executable coverage still remains visible in coverage reports through the authoritative inventory.

Each rule owner must make these facts explicit and machine-locatable using the canonical labeled-field representation in [`markdown.md`](markdown.md) under `MD.STANDARD.SECTION`:
- rule ID and classification;
- exact applicability scope, including normal exclusions;
- normative requirement;
- checker or evidence producer and its coverage relation, or an explicit statement that none exists;
- semantic remainder and required review evidence, or `None`;
- permitted exception mechanism and approval requirements.

**Automation:** `dev/audit/standards_registry.py` is registered directly under `GOV.OWNER.UNIQUE`. It enumerates every H2 owner section outside fenced examples, fails on a missing or duplicate index target and an unindexed non-backup standard, and reports owner-only, registry-only, duplicate, and mismatched identities. This is `subset` coverage because it cannot prove that two differently worded requirements are semantically unique.

**Semantic remainder:** Determine whether prose is normative, whether two provisions compete for ownership, and whether an ID should remain stable, split, merge, or retire.

**Exceptions:** Duplicate or ownerless rules are prohibited. An approved ID migration follows `GOV.CHANGE.GOLDEN_FIRST` and must reconcile all active references.

<br />

## Rule classifications (`GOV.CLASSIFICATION`)
**Classification:** `advisory` with deterministic vocabulary and presence portions.

**Scope:** The classification assigned to every normative rule and to independently registered coverage portions derived from that rule.

Each rule has exactly one classification describing its final decision boundary, regardless of whether automation has been implemented:
- `deterministic`: a checker can decide full conformance for the stated scope without human judgment. Positive and negative fixtures, or equivalent executable conformance tests, demonstrate the boundary defined by the normative owner.
- `advisory`: automation can identify bounded review candidates or decide a proper subset, but contextual human judgment determines final conformance.
- `semantic-only`: final conformance requires contextual human review and automation, if present, only collects supporting evidence.

A rule may combine deterministic coverage with a semantic remainder. The owner must state the overall classification, the independently registered deterministic or advisory coverage portion, and the remaining review obligation. If an independently enforceable provision needs a distinct lifecycle, scope, or exception policy, give it a separate rule ID. A deterministic rule that temporarily lacks a checker remains deterministic and is reported as an implementation gap; it does not become semantic-only.

Repeated semantic findings must be considered for deterministic codification when a stable boundary and representative conformance evidence can be written.

**Automation:** `dev/audit/standards_registry.py` compares each executable registry entry's `owner_classification` with its canonical owner field and requires an explicit review procedure for unregistered semantic-only owners. Coverage is `subset` because classification choice remains semantic.

**Semantic remainder:** Choose the classification that matches the rule's final decision boundary and decide when a deterministic subset needs independent registration or a separate rule ID.

**Exceptions:** A rule may not omit or hold multiple classifications. A temporary owner/registry mismatch requires an approved staged migration with an objective completion condition.

<br />

## Authority and precedence (`GOV.PRECEDENCE`)
**Classification:** `semantic-only`.

**Scope:** Conflicts among task authorization, applicable `AGENTS.md` instructions, normative owners, derived artifacts, and implementation behavior.

Task authorization and normative truth are separate. The current user request defines the requested outcome, permitted scope, and stopping boundaries. It may authorize a standards change, but it does not by itself amend a repository standard. All applicable `AGENTS.md` files provide repository-operation constraints and route work to normative owners. A more specific file refines broader instructions for its scope without discarding non-conflicting ancestor instructions; it changes a normative rule only when it is that rule's declared owner or has explicit delegated authority. An unresolved operational conflict requires clarification before affected work continues.

When determining the repository's current standard, use this precedence:
1. The normative owner defines the rule and its scope.
2. Registry data, checkers, evidence producers, fixtures, tests, manifests, suppressions, and repository contracts are derived artifacts that must agree with the owner.
3. Implementation is evidence of current behavior; it cannot silently create, broaden, or narrow a standard.

A conflict between a derived artifact or implementation and its owner is a defect in the lower-precedence layer. Derived artifacts cannot resolve conflicts among themselves by majority, age, or convention; reconcile them to the owner.

A checker must not enforce outside the rule's normative text or scope. It may cover less than the full rule only when the execution registry identifies partial coverage and the owner states the remaining advisory or semantic review obligation.

**Automation:** No checker can determine full authority or resolve semantic conflicts. Registry and scope audits may collect evidence about owner references and checker reach without replacing review.

**Semantic remainder:** Determine which instructions apply, whether provisions conflict, and whether the requested work is a conformance repair, policy change, or unauthorized scope expansion.

**Exceptions:** None reverse normative precedence. An unresolved conflict requires clarification or an approved owner change rather than an implicit exception.

<br />

## Authoritative-standard-first changes (`GOV.CHANGE.GOLDEN_FIRST`)
**Classification:** `advisory` with deterministic recorded-order portions.

**Scope:** Policy changes, conformance repairs, rule-ID migrations, staged migrations, formatter-policy changes, and repository-wide or multi-package transformations.

For a policy change, update and approve the normative owner before changing checker strictness, formatter policy, or maintained source to embody the new policy. This requirement does not apply when repairing a derived artifact or implementation that already conflicts with the current owner.

Before proposing a new obligation, governed vocabulary entry, checker facet, or rule ID, search the existing owners and record why clarification, language specialization, or routing cannot express the concern without duplication. Classify its shared or language-specific ownership and its deterministic, advisory, or semantic-only decision boundary. Show reusable benefit rather than a one-off reaction, and estimate the added ordinary instruction-loading and execution burden. Name the canonical owner, semantic remainder, exceptions, and focused proof. Default a concern already covered by an owner to clarification or specialization, and default a one-off concern to no new rule.

A standards proposal must identify motivation and source evidence, reusable benefit, existing violations, classification, deterministic or advisory coverage, semantic remainder, migration impact, and representative before-and-after examples or an explicit statement that examples are not applicable. A proposal that passes the anti-accretion screen must also estimate ordinary instruction-loading and execution burden. A local pattern becomes normative only when evidence shows that reuse reduces ambiguity or drift.

The applicable standards owner grants approval. Record the decision in a pull-request review, decision record, or other maintained change record that identifies the approving owner, affected rule IDs, and approved scope. An ID migration requires explicit approval; silence, checker output, or implementation precedent is not approval.

Related artifacts may be updated in one changeset. Their logical dependency order is normative owner → execution registry when applicable → checker or evidence producer → conformance or review evidence → repository contract → implementation. The merged repository state must reconcile. A staged migration must record temporary compatibility behavior, responsible owner, tracking record, and objective completion condition.

Do not perform a repository-wide or multi-package automated transformation while its governing rule or formatter policy is awaiting approval.

**Automation:** No dedicated checker establishes whether a change is a policy change or conformance repair. Diffs, registry audits, focused fixtures, contracts, and migration records provide bounded evidence for review.

**Semantic remainder:** Classify the change, evaluate proposal evidence, identify the approving owner, and decide whether a staged migration is narrow and objectively completable.

**Exceptions:** Repairs that reconcile derived artifacts to an already approved owner are not policy changes. Any other departure requires an approved staged migration containing the owner, tracking record, temporary behavior, and completion condition.

<br />

## Traceability (`GOV.TRACE.EXACT`)
**Classification:** `advisory` with a deterministic registry-validation subset.

**Scope:** Active paths from normative rule owners through registry coverage, checkers or evidence producers, fixtures or equivalent evidence, review obligations, repository contracts, and approved exceptions.

Traceability depends on classification:
```text
deterministic:
rule ID -> normative owner and scope -> execution-registry coverage
        -> checker -> positive/negative conformance evidence
        -> repository contract

advisory:
rule ID -> normative owner and scope -> execution-registry coverage
        -> checker or evidence producer -> review criteria
        -> semantic remainder and required evidence -> repository contract

semantic-only:
rule ID -> normative owner and scope -> review procedure
        -> required evidence and recorded disposition -> repository contract
```

For a direct interpreter, compiler, formatter, or schema-validator invocation whose baseline accepted and rejected behavior is owned upstream, repository-owned duplicate fixtures are optional only when the execution registry constrains the tool version, records invocation options, and points to a contract that runs the actual invocation. Repository-owned fixtures remain required for behavior added or changed by a local wrapper.

An owner anchor must be stable and associated with the rule ID. A moved or renamed heading must update active registry and contract references in the same reconciled change.

**Automation:** `dev/audit/standards_registry.py` and `tests/unit/dev_audit/test_standards_registry.py` validate the complete static owner inventory against executable registry entries, owner classifications, execution roles, coverage relations, covered and remaining scopes, checker and evidence paths, diagnostic families, migrations, fixtures, and allowlists. The registry coverage relation remains `subset` because the audit does not execute every registered command or decide semantic sufficiency.

**Semantic remainder:** Determine whether coverage matches the owner's real scope, whether upstream evidence is sufficient, and whether advisory or semantic review evidence supports a recorded disposition.

**Exceptions:** The direct-tool baseline and approved partial-coverage cases described above are bounded alternatives, not exemptions from traceability. Any staged gap must identify its remaining obligation and completion condition.

<br />

## Exceptions (`GOV.EXCEPTION.OWNED`)
**Classification:** `advisory` with deterministic record-schema portions.

**Scope:** Applicability exclusions, approved deviations, exception records, suppressions, owners, approvals, expiry conditions, and review triggers.

An applicability exclusion is part of a rule's normal scope and belongs in the normative owner. An exception is an approved deviation inside that normal scope. Generated, vendored, or role-specific files should normally be modeled as applicability exclusions rather than accumulated exceptions.

Every exception committed to the repository must be registered in the execution registry or another maintained exception record named by the owner. Each record must include:
- rule ID;
- narrowest affected path, symbol, or documented file role;
- exception owner and approving standards owner;
- rationale and creation date;
- expiry date or objective review trigger;
- tracking issue or decision record;
- narrowest available suppression mechanism, or `None`.

Prefer correcting source. Use an inline suppression only for a genuine local exception and a per-file suppression only for a documented file-role distinction. Unexplained blanket ignores, unregistered exceptions, and expired exceptions without renewed approval are prohibited.

**Automation:** `dev/audit/standards_registry.py` is registered under this owner with `subset` coverage. Positive and negative fixtures require each `current_exclusions_or_allowlists` record to contain an `owner=` marker. No dedicated checker validates the complete exception-record schema, approval, expiry, tracking record, or suppression narrowness.

**Semantic remainder:** Distinguish normal applicability exclusions from deviations, assess rationale and scope, choose the narrowest suppression, and decide whether renewal remains justified.

**Exceptions:** An exception to another rule may not waive the ownership, approval, registration, or review requirements defined here.

<br />

## Coverage and reconciliation (`GOV.COVERAGE.REPORT`)
**Classification:** `advisory` with deterministic inventory and reconciliation portions.

**Scope:** Full and change-triggered reconciliation of owner rules, registry entries, checker execution, findings, evidence, review dispositions, contracts, exclusions, exception records, and generated reports.

Coverage reports inventory every deterministic, advisory, and semantic-only rule from its normative owner, then reconcile that inventory with the execution registry, checker execution, findings, fixtures or equivalent evidence, review evidence, contracts, exclusions, exception records, and generated evidence. A passing checker does not erase advisory or semantic review obligations.

Reconciliation must run when a change affects any of the following:
- a maintained standard or normative owner anchor;
- execution-registry data;
- checker or evidence-producer configuration or implementation;
- fixtures or equivalent conformance evidence;
- repository contracts;
- suppressions, applicability exclusions, or exception records.

A periodic full-repository reconciliation should also detect drift that incremental checks may miss.

The report must detect at least owner rules absent from the coverage inventory; registry entries without valid owner anchors; duplicate ownership; classification mismatches; scope mismatches in either direction; checker-only or stale rule IDs; registered checkers that were not executed; deterministic rules without required positive and negative evidence; unowned, unregistered, or expired exceptions; and missing advisory or semantic review evidence.

An ownerless rule, invalid owner anchor, duplicate owner, checker broader than its owner, unrecorded classification or scope mismatch, stale emitted ID, expired exception, or missing required review evidence is a reconciliation error. The merged state is not conforming until the error is repaired or covered by an approved staged migration. A checker narrower than its owner is permitted only as an explicitly registered partial-coverage warning with the remaining review obligation identified. Availability and improvement notices may be warnings only when they do not conceal a required conformance decision.

Reconciliation verifies that required review evidence exists and records its disposition; it does not replace the human judgment required by an advisory or semantic-only rule.

**Automation:** `dev/audit/standards_registry.py` produces the complete static owner-to-execution inventory and reports owner-only automated rules, registry-only rules, duplicate owners, stale or unowned emitted diagnostics, classification mismatches, and incomplete coverage metadata. It does not execute every checker, validate the full exception lifecycle, or decide advisory and semantic dispositions, so registered coverage remains `subset`.

**Semantic remainder:** Evaluate classification and scope parity, decide whether review evidence is sufficient, and determine whether warnings represent safe partial coverage or a concealed conformance decision.

**Exceptions:** Only an approved staged migration may leave a temporary reconciliation gap, and it must identify the affected rule IDs, responsible owner, warning state, tracking record, and objective completion condition.
