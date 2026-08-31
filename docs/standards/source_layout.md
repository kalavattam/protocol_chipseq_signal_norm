
# Cross-language source-layout standard
This document owns source-layout semantics shared by maintained Shell, Python, R, and Rust. Language standards own the syntax that realizes these semantics, numeric candidate thresholds, formatter and checker selection, applicability exclusions, and language-specific exceptions. [`help.md`](help.md) retains shared help and documentation vocabulary, [`source_headers.md`](source_headers.md) retains source-header profiles, and [`markdown.md`](markdown.md) retains Markdown source form.

<br />

## Shared and language-specific ownership (`SOURCE.LAYOUT.OWNERSHIP`)
**Classification:** `semantic-only`.

**Scope:** Source-layout provisions intended to apply across maintained implementation languages and every language-specific standard that realizes them.

This standard is the unique normative owner for common semantic paragraphs, density governance, candidate-evidence concepts, comment attachment, grammatical naming roles, multiline-structure semantics, and comparison-symbol context. A language standard must link to these owners instead of restating a competing shared rule.

Each language standard owns:
- recognized syntax and construct families;
- ordinary and documentation comment markers;
- indentation, source width, and trailing-comma behavior;
- numeric values for candidate parameters;
- formatter, linter, parser, compiler, and checker relationships;
- applicability exclusions and approved exceptions; and
- language-specific migration and proof.

Common meaning does not require identical source spelling or numeric policy. A formatter may decide a language-native representation without deciding semantic grouping, naming clarity, comment usefulness, or whether a measured candidate requires change.

The Shell and Python realizations apply to their current maintained source. The R and Rust realizations become applicable when maintained source in those languages is introduced. Introducing such source requires selecting its owned repository boundary and proportional proof; it does not require duplicating this shared standard.

**Automation:** No checker can decide this ownership boundary. The standards-registry reconciliation inventories unique owner IDs and indexed standards but does not determine whether two differently worded provisions compete semantically.

**Semantic remainder:** Review new or changed language provisions to determine whether they realize a shared concept or establish a genuinely language-specific requirement.

**Exceptions:** A language standard may refine syntax or applicability but must not redefine a shared semantic owner. A change to shared meaning requires an approved change here.

<br />

## Semantic paragraphs and density (`SOURCE.LAYOUT.PARAGRAPHS`)
**Classification:** `advisory` with deterministic Python subsets.

**Scope:** Statements, expressions, declarations, calls, commands, pipelines, constructed objects, branches, handlers, and result or control-transfer regions in maintained implementation source.

Organize source into visible semantic paragraphs. Keep statements that perform one tight operation contiguous, and put one blank line between distinct phases such as setup, acquisition, validation, transformation, side effects, cleanup, and result production.

Put a visible blank-line boundary before and after a noncompact sibling control-flow unit so a reader can distinguish where one branch, loop, match or case dispatch, handler, or managed-resource phase ends and the next sibling statement begins. Keep language-connected clauses such as `elif`/`else`, `else if`/`else`, `elsif`/`else`, handlers, and closing delimiters attached to their owning construct. A direct single-transfer guard may remain compact with a related validation sequence.

Do not insert a blank line after every assignment, call, command, or short statement. Related short repeated constructions remain contiguous when they form one topic. Separate independently meaningful constructions when their size or responsibility makes each a clearer paragraph, and use a brief phase comment or a helper boundary when blank lines alone do not expose the organization.

Prevent both failure modes:
- hyper-density, where distinct phases or topics are compressed into one undivided region; and
- mechanical fragmentation, where related statements become a row of isolated one-statement paragraphs.

Review every region at or above either current p90 density measure, every recorded largest region, and every additional region identified as difficult to scan by human review. Inspect the actual statement subjects and record why the members form one topic or where phase boundaries, phase comments, or helper extraction were introduced. A generic source hash, construct label, file role, or claim that the region is “coherent” is not a semantic decision. New or changed code must undergo the same review before its language migration can be called complete.

Keep syntactically connected structures attached according to the language. Decorators or attributes remain attached to declarations; branch and handler clauses remain attached to their chains; a block header remains attached to its first governed statement; and delimiters, heredocs, macro bodies, and pipeline syntax retain their language-required continuity.

A compact guard remains compact internally when it performs a direct result or control transfer with no substantive preparatory phase. Once that standalone guard completes, put one blank line before its following sibling statement; this makes the transferred path end visibly and separates back-to-back short guards. Separate the transfer within the guard when preparation, diagnostics, mutation, cleanup, or result construction forms a distinct phase. After a substantive result-building phase, visually separate the terminal result or explicit transfer. Language standards identify applicable forms such as `return`, `raise`, `yield`, `exit`, `stop()`, `break`, `continue`, tail expressions, diverging macros, and error propagation.

Layout must not imply false independence. A blank line may clarify phases, but it must not conceal a shared invariant, split one atomic operation, or substitute for a needed refactor.

**Automation:** `dev/audit/python_source_policy.py` deterministically checks the safely provable subsets comprising the sibling boundary before and after every noncompact Python control-flow unit; the following-sibling boundary after every standalone compact transfer guard; transitions in every suite nested within a callable involving local imports, result-bearing action assignments, validation or verification after a substantive setup paragraph, completed validation followed by a new phase, assertions, and multiline data or case inventories immediately consumed by a loop; and recognized result or control transfers after an attached substantive phase. Result-bearing action recognition uses explicit semantic result names and keeps consecutive result captures together; tests receive strict arrange–act separation, while non-test callables require substantive attached preparation before an initial action boundary. One mutation or observation may remain attached to the validation that proves it. Transfer recognition covers `return`, `raise`, `yield`, `yield from`, `break`, `continue`, and statically named process-exit calls. A transfer remains compact only after at most one single-line, non-mutating, non-diagnostic, non-control-flow preparatory statement. `dev/audit/python_source_evidence.py` measures blank-line-delimited Python regions, physical and logical span, statement count, construct mix, nesting, statement subjects, comment presence, and inspectable transfer-phase membership across maintained Python. X-exceeding constructions and changes of subject within one syntactic phase remain semantic review rather than mechanical blank-line mandates. Shell, R, and Rust have no corresponding semantic-layout adapter yet; this absence is recorded scope, not evidence that the shared rule is inapplicable. New or changed code in every maintained language is review-owned by this shared rule now; an absent adapter never licenses hyper-density or mechanical fragmentation, and language-specific automation requires its own approved migration.

**Semantic remainder:** Review topic boundaries, phase changes, guard compactness, transfer separation, every p90-or-larger or largest density record, human-identified scan failures, and whether dense or fragmented code should be regrouped, commented, or refactored.

**Exceptions:** A language-required connected form remains attached. Any other retained dense or fragmented region requires a specific semantic rationale or an approved exception.

<br />

## Candidate evidence and dispositions (`SOURCE.LAYOUT.CANDIDATES`)
**Classification:** `advisory` with deterministic candidate-record opportunities.

**Scope:** Candidate generation, threshold selection, evidence records, semantic dispositions, pilots, migrations, and reconciliation for shared source-layout and naming rules.

`X`, `Y`, and `Z` are shared evidence concepts, not universal numeric laws:
- `X` measures when one independently meaningful multiline construction becomes a review candidate because of its physical span or structure;
- `Y` has two distinct, unsummed measurements: physical lines in one multiline assignment and statement count in one blank-line-delimited simple-assignment run; and
- `Z` measures tightly related preparatory physical lines inside a compact guard before its result or control-transfer form.

No numeric value is established by this standard. Each language must select, retain, adjust, or reject numeric values from representative repository evidence after accounting for its indentation, formatter, construct grammar, line width, and source roles. A candidate threshold identifies evidence for review; it does not make the candidate a violation.

A deterministic candidate record must contain enough stable facts for independent review:
- language, path, enclosing source object, source location, and recognized construct;
- reason, configured parameter or threshold, physical and logical spans, statement or item count, and nesting information relevant to the finding;
- a source fingerprint or equivalent stable identity; and
- the evidence-producer version and configuration when automation exists.

Every in-scope candidate requires one durable semantic disposition:
- changed to conform;
- retained as one genuinely tight semantic unit, with a specific rationale;
- refactored because spacing or renaming alone could not make the source clear; or
- excepted through the governed exception process.

A pilot must compare adjacent candidate values, inspect representative short, boundary, long, nested, repeated, and formatter-sensitive cases, and report false positives as well as useful findings. It must measure paragraph-level density separately rather than treating `X` or `Y` as a paragraph-size limit. A migration reconciles only when every required candidate has a disposition and no unresolved candidate is hidden by a clean deterministic checker.

Promotion from review-only architecture requires maintained source, representative pilot evidence, explicit language thresholds or a recorded decision not to use them, positive and negative fixtures for deterministic subsets, an evidence schema, formatter and checker relationships, idempotence proof, and an approved migration boundary.

**Automation:** `dev/audit/python_source_evidence.py` emits schema-versioned X, both separately named and unit-bearing Y series, and Z candidates plus complete repeated-run members, all-suite transfer records, and separate syntactic blank-line-region density for the bounded Python pilot. An explicit decision names exact stable record keys and its complete inspectable membership. A source or membership fingerprint may activate or invalidate that already recorded decision; neither a matching cohort hash nor a global completion switch may create a disposition. Focused fixtures validate strict threshold boundaries, complete records, signatures, schema shape, idempotence, and stale-membership invalidation. Candidate dispositions and threshold approval remain review-owned.

**Semantic remainder:** Select representative cohorts, judge threshold usefulness, disposition every candidate, and decide whether a measured subset is reliable enough for deterministic enforcement.

**Exceptions:** A migration may not treat missing evidence or unresolved candidates as implicit conformance. An approved exception must satisfy [`governance.md`](governance.md).

<br />

## Comment attachment (`SOURCE.COMMENT.ATTACHMENT`)
**Classification:** `advisory` with deterministic marker and adjacency opportunities.

**Scope:** Ordinary full-line, continuation, paragraph-separator, inline, directive, and language-native documentation comments in maintained implementation source.

An ordinary comment explains purpose, rationale, constraints, or non-obvious behavior and stays accurate with the source it governs. Put one blank line before a phase comment when it begins a new semantic paragraph, and put no blank line between that comment and the paragraph it introduces. Keep inline comments short; promote a long explanation to an attached preceding block or maintained documentation.

Ordinary comment prose uses sentence capitalization, complete sentences, terminal punctuation, and one space between sentences when it forms prose. Punctuate a parenthetical or an abrupt break with a comma, colon, semicolon, or parentheses; a double hyphen is not a prose dash, because `--` reads here as an option prefix or an end-of-options delimiter. A multiline ordinary comment uses the language's ordinary continuation marker and its language-native empty separator. Do not mechanically rewrite shebangs, source headers, formatter or linter directives, coverage pragmas, generated markers, literal fixtures, or documentation comments as ordinary prose.

Language owners define ordinary, continuation, empty-separator, inline, directive, and documentation-comment spelling. That syntax realizes these shared attachment, prose, and role semantics without creating a second shared owner.

The shared categories use language-native spelling:
- Python, Shell, and R ordinary comments use `# ` and use exactly `#` as an empty paragraph separator;
- Python and Shell inline comments use at least two spaces before `#` and one space after it where inline comments are permitted, and each language owner decides whether a bounded shared alignment column may exceed that minimum;
- R reserves `#'` for roxygen2 documentation rather than ordinary implementation comments; and
- Rust ordinary line comments use `// `, outer item documentation uses `///`, and inner module or crate documentation uses `//!`.

Language documentation syntax carries a documentation contract, not merely visual emphasis. The language owner decides when a docstring, roxygen2 block, Rust doc comment, Shell help heredoc, or other documentation surface applies. Shared section vocabulary and content semantics remain with [`help.md`](help.md).

The Python forms follow [PEP 8 comments](https://peps.python.org/pep-0008/#comments). Rust ordinary and documentation comment distinctions follow the [Rust Reference](https://doc.rust-lang.org/reference/comments.html). Roxygen2 documentation blocks follow the [roxygen2 documentation](https://roxygen2.r-lib.org/articles/roxygen2.html).

**Automation:** No dedicated registry entry exists for this shared owner. The registered `PY.COMMENT.FORM` execution delegates Python's deterministic realization to the Python owner rather than duplicating a shared execution row. `dev/audit/python_source_policy.py` checks ordinary Python comment markers, separators, inline spacing, trailing whitespace, width, and safely identifiable prose-paragraph capitalization, terminal punctuation, and sentence spacing across maintained Python. It evaluates a contiguous wrapped comment as one paragraph so continuation lines do not receive sentence-ending punctuation mechanically. Directives, headers, literal labels, and statically recognizable code or operator fragments remain excluded from prose checks. No Shell, R, or Rust prose adapter exists yet; those languages remain review-owned by this shared rule.

**Semantic remainder:** Classify comment role, determine the governed source paragraph, review prose usefulness, and distinguish documentation from implementation commentary and directives.

**Exceptions:** Specialized markers retain their separately owned syntax. A detached ordinary comment requires a documented reason when adjacency would misstate its scope. A literal `--` inside a command, option name, example, or fixture is delimiter syntax rather than prose punctuation and remains unchanged.

<br />

## Greedy ordinary source prose (`SOURCE.PROSE.WRAP`)
**Classification:** `advisory` with deterministic recognized Python portions.

**Scope:** Eligible ordinary source prose whose language owner can distinguish prose words and indivisible units from structure, literals, formulas, URLs, directives, and generated content.

Wrap at the last complete word or indivisible unit fitting through physical column 79. Break only when the next complete unit and required separator would exceed the boundary. A width-only check does not prove greedy wrapping. Language owners recognize syntax and exclusions without restating this shared meaning.

Python adjacent constant CLI `help=` literals form the only approved automatic formatter subset. The formatter preserves the exact evaluated value and keeps terminal escapes attached to the final text-bearing literal. Python comments and docstrings are checker/evidence-first. Shell comments remain review-owned; R and Rust remain dormant.

A structural boundary is not a premature break. A section header, a documented entry header, a multiline textual-type row, a dedent that ends a block, a fenced verbatim block, and a literal example row each end a prose line for reasons the width boundary does not govern; joining them would destroy meaning rather than restore greedy wrapping. A language realization must recognize its own structural boundaries before it reports any break.

Embedded verbatim content is not prose. A block fenced by a line holding only `'''`, ```` ``` ````, or `~~~` inside a documentation comment holds pseudocode, sample output, or another literal whose line breaks are part of its meaning. Neither the fences nor any line between them is a wrap candidate, and no automated repair may rewrite them.

A line may also break inside a word. A trailing hyphen touching the character before it continues one word such as `whitespace-separated`, so rejoining the two lines restores the word and inserts no space; the fit test must not reserve a separator that will not exist. A hyphen with a space before it is a minus sign, dash, or option prefix and keeps its separator. A suspended hyphen — `zero- or negative-length` — reads identically to a broken word but is two words sharing one hyphen, and rejoining it would corrupt the prose; width cannot distinguish the two, so it belongs to review.

Indivisible means the unit cannot be split, not that it cannot be moved. A complete quoted token, a bracketed literal, or any other unit that already occupies one whole word moves to the preceding line intact and is therefore an ordinary candidate. Exclude a unit only when filling forward would move part of it and leave the rest behind, as when a quotation, formula, or bracketed expression spans several whitespace-separated words.

Distinguish where a continuation sits from where the line broke. A list item may continue at the marker's own continuation column or align under the item text; both are indentation choices, and neither exempts the item from the ordinary fit test. Exclude the column, never the wrap point. This distinction is narrow on purpose: a deeper following line is a continuation only when a list marker introduced it, since an entry header followed by its indented description is a structural boundary and joining the two would destroy the entry.

**Automation:** `dev/audit/python_source_policy.py` checks the recognized Python subset. It emits `SOURCE.PROSE.WRAP` for premature breaks in docstring prose paragraphs, excluding section headers and underlines, NumPy `name : type` entry headers, multiline textual-type rows, block-ending dedents, a following list marker, doctest rows, `Examples` sections, fenced verbatim blocks, colon-introduced indented display blocks, and indivisible units. A list item whose continuation is aligned under the item text keeps its wider column but is still fit-tested; a deeper following line without a list marker remains excluded. Recognized `parse_args()` help literals emit under `PY.CLI.HELP.LAYOUT` and adjacent constructed prose under `SOURCE.DELIMITED.MULTILINE`, so no two owners report one break. Ordinary Python comment prose remains review-owned. `dev/tools/python_help_format.py` formats only eligible adjacent constant help literals with explicit `--write`. Coverage is `subset`.

**Semantic remainder:** Decide whether a unit is indivisible, whether a break is semantically deliberate, and whether prose is eligible.

**Exceptions:** Structural, literal, formula, URL, directive, generated, dynamic, and ambiguous content is never blindly rewritten.

<br />

## Grammatical naming and migration (`SOURCE.NAMING.SEMANTICS`)
**Classification:** `advisory` with deterministic inventory opportunities.

**Scope:** Maintained identifiers, test names, source filenames, non-production script names, public command names, abbreviations, and rename migrations across implementation languages.

Names communicate role and current meaning:
- action functions and commands use concise verb-object phrases;
- predicates use a truth-revealing verb when it improves clarity;
- variables, fields, attributes, and properties use noun phrases describing their values;
- classes, exceptions, protocols, traits, and other types use noun phrases describing the modeled concept;
- tests use a subject-behavior-outcome form that represents one coherent contract; and
- source filenames and scripts describe their maintained role without repeating unnecessary directory context.

Names, together with their owned interface and context, must communicate cardinality and value shape truthfully and intelligibly. This shared owner does not require every implementation language or every identifier to encode collection-versus-scalar status through a prefix or other lexical marker. Each language owner decides whether a lexical distinction applies and owns any language-native prefix, container form, or spelling rule.

Controlled scientific, file-format, ecosystem, and project abbreviations are permitted when unambiguous in context. Avoid private shorthand, opaque truncations, redundant type suffixes, copied historical abbreviations, numeric suffixes that conceal a collection or semantic distinction, and names that require repository history to decode. A rename must improve meaning rather than merely increase length.

When a migration establishes a descriptive replacement for private shorthand, record that vocabulary in the language owner and enforce it for new and changed implementation names. An exact external field, option destination, protocol spelling, scientific symbol, file format, or conventional ecosystem alias may remain only at its bounded interface. This requirement applies across implementation languages even when only one language currently has a deterministic adapter.

Name length is evidence, not a universal conformance limit. A role-aware inventory should record language, path, enclosing scope, role, visibility, location, spelling, character and segment counts, abbreviation tokens, scope span, use count, shadowing, repeated context, and public, dynamic, serialized, reflection, generated-command, fixture, and documentation references.

Too-short candidates include unexplained single-letter names outside tiny mathematical, coordinate, index, generic-parameter, or conventional argument scopes; opaque abbreviations; generic placeholders whose value is unclear; and booleans without evident truth meaning. Too-long candidates include repeated namespace context, responsibility-heavy qualifier chains, redundant type words, tests containing multiple contracts, and names that encode implementation history. Language owners select role-specific candidate thresholds only after evidence.

Approved internal and non-production renames are direct atomic changes. Update definitions, imports, calls, tests, fixtures, documentation, help, entry-point metadata, reflection, serialization, and generated-command expectations together. Do not leave aliases, forwarding wrappers, deprecated spellings, or duplicate exports merely to preserve internal history.

Production wrapper names and production Python command names are protected. They change only through a separately authorized public-interface migration with compatibility and proof obligations. Historical immutable evidence may retain retired spellings when it cannot be mistaken for a live interface.

Language standards own casing, visibility markers, dispatch punctuation, constant forms, lifetime syntax, and externally required spellings. Adapt externally owned names at a boundary when the language and interface permit it.

**Automation:** No shared naming inventory, proposal ledger, reference reconciler, fixture set, or execution-registry entry currently exists. The Python adapter rejects the bounded opaque-shorthand vocabulary established by its completed migration; other language and interface checks cover only their registered naming subsets.

**Semantic remainder:** Judge grammar, abbreviation clarity, length candidates, compatibility risk, external ownership, dynamic references, and whether a proposed rename improves the represented concept.

**Exceptions:** Conventional short names and exact external spellings may remain when their scope and ownership are clear. A public interface change follows its existing owner and [`governance.md`](governance.md).

<br />

## Multiline delimited structures (`SOURCE.DELIMITED.MULTILINE`)
**Classification:** `advisory` with deterministic language-specific formatting opportunities.

**Scope:** Calls, argument and parameter lists, literals, builders, collection construction, arrays, macro invocations, and analogous delimited structures in maintained implementation source.

Keep a delimited structure on one line when it fits its language's approved width and remains readable. Once broken:
- break immediately after the opening delimiter and before the closing delimiter;
- put one argument or item on each line unless every call argument forms one bounded, cohesive continuation row or the language owner defines another narrow native exception;
- use one continuation indentation level rather than fragile visual alignment;
- align the closing delimiter with the line that begins the construction; and
- do not place the first argument on the opening line while vertically listing later arguments.

Keep a comparison operator on the final physical line of its left operand. When the left operand is a delimited call or construction, wrap that operand internally so its closing delimiter and the comparison remain together. Do not use an enclosing group merely to strand the comparison on the following line. If the comparison still does not fit, extract an honestly named intermediate value.

Keep the boolean relationship in a compound condition visually available. When formatting a membership collection would push the collection's items into the visual position normally occupied by the condition's logical operands, extract truthfully named predicates and compose the final condition from those names. Apply the same principle to other formatter-owned expression shapes that conceal `and`, `or`, negation, comparison, or guard intent; do not fight a required formatter with unstable manual wrapping.

Do not introduce a multiline wrapper around one short trailing type or result clause merely to shorten a definition header. Keep that simple clause with the definition when it fits; wrap the meaningful parameter or result structure internally when the type itself is genuinely complex.

Put every multiline triple-quoted or language-equivalent block-literal opening and closing delimiter at its own content boundary, regardless of assignment, call, fixture, or generated-text context. Begin content on the line after the opener and end content on the line before the closer; do not attach literal content to either delimiter. These boundaries make the block's first and last content lines as visible as every interior line. When the runtime value must omit either boundary line feed, preserve that value explicitly with a readable boundary-normalization operation such as `removeprefix()` or `removesuffix()`; do not attach content to a delimiter.

Apply the same construction-line alignment to multiline textual type expressions in language documentation. Indent inner type rows as continuations, but align a closing-only delimiter row with the logical entry line that introduced the type construction.

Trailing-comma behavior is language-specific. Python compact call rows omit the trailing comma that would make the formatter expand them; other Python and Rust multiline calls and analogous supported structures use trailing commas. R multiline calls do not use a trailing comma. Shell command arrays and multiline invocations follow Shell syntax and do not acquire Python punctuation. Each language owner defines its other formatter interactions and syntax-required exclusions.

When one logical prose value is built from adjacent or continued source fragments, fill each physical source line through the last complete word that fits within the language's approved width. Do not leave a word, punctuation-only fragment, or other short tail in its own fragment when it can move to the preceding fragment without changing the runtime value. Preserve explicit newline escapes, structured generated text, indivisible tokens, and deliberate semantic boundaries.

Wrapping is not always the clearest repair. Extract a meaningful intermediate value, helper, or data structure when nesting or expression density remains difficult to scan after canonical wrapping.

Python hanging indentation and trailing-comma guidance is grounded in [PEP 8](https://peps.python.org/pep-0008/#indentation). R call layout follows the [tidyverse syntax guide](https://style.tidyverse.org/syntax.html#long-lines) as the provisional ecosystem basis. Rust block indentation and multiline trailing commas follow the official [Rust Style Guide](https://doc.rust-lang.org/style-guide/expressions.html).

**Automation:** `dev/audit/python_source_policy.py` checks recognized multiline Python function and method parameter lists, simple return-annotation source form, calls, and list, tuple, dict, and set displays for bounded opener, closer, continuation, canonical expanded or compact-call item rows, alignment, and trailing-comma forms across maintained Python. It rejects a comparison operator placed on a later physical line than its left operand, a multiline parenthesized wrapper around one simple return annotation, a compound `if` or `while` condition whose expanded collection obscures its boolean operands instead of using named predicates, content attached to either delimiter of any recognized multiline triple-quoted literal, and safely provable premature word- or punctuation-boundary breaks in adjacent static prose fragments while excluding dynamic interpolation, explicit newlines, structured text, and separately owned `parse_args()` help literals. The `PY.DOCSTRING.NUMPY` adapter exclusively reports multiline NumPy textual-type closer alignment so the shared owner does not create duplicate diagnostics. A compact nested value that fits on one line remains one valid outer item. Sole-generator-expression calls, comprehensions, unsupported construct families, and Shell, R, and Rust language adapters remain explicit exclusions. No cross-language formatter exists, so the shared constructed-prose, comparison-attachment, block-boundary, and density requirements remain mandatory review obligations in those languages until their own adapters are approved.

**Semantic remainder:** Decide whether a structure should remain inline, wrap, be extracted, or use a documented language-native exception; decide whether a non-prose or explicit fragment boundary is deliberate. Shell, R, and Rust require their own language adapters and migrations before this shared rule can be claimed as enforced there.

**Exceptions:** Syntax-required or formatter-owned forms may differ only through the applicable language owner. Literal fixtures and generated source follow their owned applicability exclusions.

<br />

## Comparison symbols by context (`SOURCE.COMPARISON.CONTEXT`)
**Classification:** `advisory` with deterministic character-inventory opportunities.

**Scope:** Comparison relations in executable source, identifiers, filenames, configuration, tests, diagnostics, comments, docstrings, help, and human-facing scientific or mathematical prose.

Use ASCII operators such as `>=` and `<=` in executable and machine-oriented contexts, including source expressions, tests, parsers, configuration, command examples, exact diagnostics, identifiers, filenames, option names, and machine-readable output. Prose that quotes executable syntax uses backticked ASCII operators.

Use `≥` and `≤` in human-facing scientific or mathematical prose when the symbols express a mathematical relation rather than quote code. In source comments and language-native documentation, use ASCII when the text directly explains a code predicate; reserve mathematical symbols for genuine mathematical exposition.

A character inventory can identify source regions and symbols but cannot decide the intended audience or whether prose is mathematical. Rewriting must preserve exact diagnostics, machine-readable output, fixtures, identifiers, and externally owned text.

**Automation:** No dedicated checker or registry entry currently implements this contextual distinction. A future evidence producer must classify source regions conservatively and leave ambiguous prose for review.

**Semantic remainder:** Determine whether each relation is executable syntax, machine-oriented text, quoted code, or mathematical prose.

**Exceptions:** Exact upstream quotations, immutable historical evidence, and externally owned interface text may retain their original symbols with a recorded disposition.
