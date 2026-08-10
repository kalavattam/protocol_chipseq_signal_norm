
# JSON standard
This document owns the source form of maintained JSON in this repository. It governs how a JSON document is laid out, not what it may contain: schema validity, field semantics, and applicability belong to whichever owner consumes the document.

JSON is data rather than implementation source, so this document deliberately adopts fewer obligations than a language owner. It realizes [`SOURCE.DELIMITED.MULTILINE`](source_layout.md#multiline-delimited-structures-sourcedelimitedmultiline) for JSON and adds nothing else. Key naming is not restated here; [`SOURCE.NAMING.SEMANTICS`](source_layout.md#grammatical-naming-and-migration-sourcenamingsemantics) already permits an exact external field or protocol spelling to remain at its bounded interface, which is what governs JSON Schema keywords such as `additionalProperties`, `patternProperties`, and `minLength`.

<br />

## JSON source form (`JSON.SOURCE.FORM`)
**Classification:** `advisory` with a deterministic canonical rendering.

**Scope:** Maintained authored JSON under `dev/config/`, `dev/schemas/`, and `tests/fixtures/`. Serialized inventories emitted by a maintained producer are excluded, because their form is chosen by the producer's serializer rather than by an author.

One canonical rendering defines the form. A maintained JSON document is conformant when it is byte-identical to that rendering of its own value:
- the document parses, and no object declares the same key twice;
- indentation is two spaces per nesting level, with no tab anywhere in the file;
- a key is followed by `": "`, and inline members are separated by `", "`;
- a structure is written inline when its complete physical line — indentation, key prefix, the structure itself, and any trailing separator — fits within 79 columns, and is otherwise expanded one member per line;
- an expanded structure puts its opening delimiter last on its line and its closing delimiter first on its line, aligned with the line that began the construction; and
- the file ends with exactly one newline, and carries no blank line and no trailing whitespace.

Member order is whatever the author wrote. The rendering never reorders anything.

<br />

### The budget decides structure, not line length
The 79-column figure is a structural trigger, not a limit. It chooses between the inline and expanded treatments of a *structure*; it makes no claim about the width of any resulting line. An indivisible scalar longer than the budget stays exactly as it is, because JSON has no string continuation and a rule promising a column maximum could not be honored. This is the one place where a number shared with [`SHELL.LINE_LENGTH`](shell.md#shell-line-length-shelllinelength) means something different: there it bounds a line, here it selects a form. The figure is shared for consistency across the corpus, not because the two obligations are the same.

A consequence worth stating plainly, because it looks like an inconsistency and is not: the same key may be inline in one place and expanded in another, within a single file. That is the fit test working. Two occurrences of one key differ in treatment exactly when they differ in width, and forcing them to agree would expand short structures for no reason.

<br />

### Three populations, one of which this rule cannot reach
Applicability is the hard part of this rule, not style.

| Population                                                                | Who owns the form                                                                                                                               |
| :---                                                                      | :---                                                                                                                                            |
| **Authored configuration** — `dev/config/`, `tests/fixtures/*/cases.json` | This rule.                                                                                                                                      |
| **JSON Schema** — `dev/schemas/*.schema.json`                             | This rule owns the *form*. The vocabulary is externally owned and delegates to `SOURCE.NAMING.SEMANTICS`; no naming obligation is imposed here. |
| **Serialized inventories** — JSON emitted by a maintained producer        | The producer's serializer. An author-facing rule cannot reach a file no author writes, so changing its form means changing the emitting call.   |

A record-per-line array — an expanded array whose elements are short inline objects — is a deliberate and good treatment for tabular data, and the canonical rendering preserves it wherever the rows fit the budget. It is the JSON analogue of the bounded cohesive continuation row that `SOURCE.DELIMITED.MULTILINE` already permits. Where the rows do not fit, they expand like anything else; a table is not a licence to run past the budget.

<br />

### Explicit non-obligations
These are recorded so that a later pass does not adopt them as improvements.

**Keys are not sorted.** No maintained JSON object has ever been key-sorted, and none will be. Alphabetical order would be a large migration that destroys meaningful grouping in exchange for nothing. A serialized inventory may well be sorted, because its serializer sorts it; that is the producer's choice and says nothing about authored configuration.

**Mixed treatment of one key within a file is not a defect.** See the fit test above.

**No maximum line length is imposed.** A line is as long as its content requires.

**Automation:** `dev/audit/json_source_form.py` reports every departure from the canonical rendering, and `--fix` rewrites a document into it. Every finding names one recognized class — budget, delimiter placement, indentation, tabs, whitespace, trailing newline, or unreadable input — and the checker reports an explicit unexplained-difference finding rather than passing a file it cannot account for. A document with a duplicated key is reported as unreadable and never rewritten, because rewriting it would silently discard a value. `tests/contract/repository/test_json_source_form_gate.sh` runs the checker over the repository and asserts the inspected-path count. Coverage is `subset`, not `exact`.

**Semantic remainder:** Decide whether a record-per-line array is genuinely tabular and worth keeping inline where the budget permits either form, whether a new JSON file belongs to the authored or the serialized population, and whether an inline structure reads as one unit. The checker decides form; it does not decide that a document should exist, or which population owns it.

**Exceptions:** Serialized inventories under `dev/audit/` are excluded by applicability rather than by allowlist, and are owned by the `--inventory-output` paths of `dev/audit/ai_attribution.py` and `dev/audit/help_aliases.py`. Generated fixtures are ignored by Git and therefore never discovered.
