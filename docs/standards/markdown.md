
# Markdown standard
This document owns canonical source form for every maintained Markdown file in the repository. Domain standards may add content requirements, but they must not define competing Markdown syntax or spacing. The rendering baseline is CommonMark plus GitHub Flavored Markdown pipe tables; the exact XHTML `<br />` section separator defined here is an intentional repository convention.

<br />

## File boundaries (`MD.FILE.BOUNDARY`)
**Classification:** `deterministic`.

**Scope:** The beginning and end of every maintained Markdown file.

Begin the file with exactly one empty line before its first heading unit. A heading unit is an optional canonical explicit HTML anchor plus the heading on its immediately following source line. In byte terms, an ordinary UTF-8 file without a byte-order mark begins with one line-feed character followed immediately by the optional anchor when present, or by the heading when absent. No prose, front matter, comment, or additional empty line precedes that unit.

End the file with exactly one terminal line feed after its final content and no additional empty line. In byte terms, the final nonempty source line is followed by one line-feed character and then end of file. Two or more terminal line feeds create one or more extra empty lines.

**Automation:** `dev/audit/markdown_policy.py` and its accepted and rejected fixtures check the exact leading and terminal line-feed boundaries and recognize an optional canonical first-heading anchor without treating it as pre-heading content. The standing registry entry has `exact` coverage for maintained Markdown paths.

**Semantic remainder:** `None` for source boundaries.

**Exceptions:** Literal fixture content may demonstrate invalid boundaries, but the fixture file containing that literal example still follows this rule unless it is an explicitly rejected fixture.

<br />

## Headings (`MD.HEADING.SPACING`)
**Classification:** `deterministic`.

**Scope:** Formal H1-H6 headings and structural informal H7/H8 headings outside fenced code.

Use ATX syntax with one space after the marker for formal H1-H6 headings. Put no blank line between any formal or informal heading and its first content block. Put exactly one blank line before a GFM table or the `<br />` row of an immediately following ordinary `MD.SECTION.BREAK` boundary.

For spacing checks, an optional canonical anchor belongs to the following heading unit; `MD.ANCHOR.CANONICAL` governs the anchor-to-heading edge. A contentless direct-child heading unit follows its parent heading with no intervening blank line. Anchor presence does not change spacing after the heading itself.

CommonMark recognizes ATX headings only through H6. Do not write H7/H8 as lines beginning with seven or eight `#` markers; use the repository-defined informal forms in `MD.HEADING.INFORMAL`.

Canonical examples:
````markdown

# Document title
Opening content begins immediately.

## Table section

| Name | Value |
| :--- | :---  |
| A    | 1     |
````

**Automation:** `dev/audit/markdown_policy.py` checks canonical formal-heading spacing, recognizes an optional anchor as part of the following heading unit, and applies the table and section-boundary spacing cases. The registered relation remains `subset` because informal-heading intent outside safely recognized structural contexts remains review-owned.

**Semantic remainder:** `None` for spacing once a line is classified as a heading. Choosing the correct formal or informal hierarchy remains semantic document design.

**Exceptions:** None for heading adjacency; table, section-boundary, and contentless direct-child forms are classifications, not exceptions.

<br />

## Informal headings (`MD.HEADING.INFORMAL`)
**Classification:** `advisory` with deterministic canonical forms.

**Scope:** Standalone bold or italic lines used as structural headings outside fenced code.

The repository defines two informal heading levels below H6:
```markdown
**H7 section title**
*H8 section title*
```

An H7 line wraps its entire visible title in Markdown bold, `<strong>...</strong>`, or `<b>...</b>`. An H8 line wraps its entire visible title in Markdown italics, `<em>...</em>`, or `<i>...</i>`. Prefer the Markdown forms unless raw HTML is required by the surrounding source.

A standalone emphasized line is safely recognized as structural H7/H8 when its role is established without depending on an already-conforming section boundary: it is the first headed entry; it has an immediately attached heading-anchor candidate; it follows a recognized H6 or H7 parent as a direct child with no intervening real content, whether the candidate boundary is conforming or malformed; or, outside an open `<details>` body, it follows an ordinary or details-close boundary whose surrounding structure independently identifies it as a section heading. H7 follows H6 and H8 follows H7; do not skip directly from H6 to H8.

Every recognized H7/H8 heading follows `MD.HEADING.SPACING`, `MD.SECTION.BREAK`, `MD.HEADING.CASE`, and `MD.ANCHOR.CANONICAL` exactly as a formal heading does. A malformed boundary does not cause an otherwise safely recognized informal heading to be reclassified as ordinary emphasis.

An emphasized-only line whose structural intent cannot be established safely remains an advisory review candidate rather than being silently accepted as ordinary emphasis or deterministically rewritten. Human review decides whether it was intended as H7/H8 and, if so, applies the same boundary rules.

**Automation:** No independent deterministic `MD.HEADING.INFORMAL` checker or registry row decides author intent. The shared structural parser recognizes the safe contexts above without requiring the boundary to conform, and the Markdown policy audit records otherwise ambiguous emphasized-only candidates for semantic review.

**Semantic remainder:** Decide whether an otherwise ambiguous emphasized-only line was intended as H7/H8 and whether its hierarchy is semantically correct.

**Exceptions:** None permit a recognized H7/H8 heading to use different spacing, boundary, or anchor rules, or seven- or eight-marker ATX-like source as a canonical substitute.

<br />

## Heading text (`MD.HEADING.CASE`)
**Classification:** `advisory` with deterministic terminal-punctuation evidence.

**Scope:** The visible title text of formal H1-H6 and informal H7/H8 headings.

Prefer sentence-like capitalization rather than Title-Like Capitalization. Capitalize the first word and preserve required capitalization for proper nouns, acronyms, scientific names, product names, and literal code identifiers; otherwise use ordinary sentence case. Do not end a heading with a period or full stop. Other terminal punctuation is permitted only when it carries real meaning.

Do not mechanically recase heading text. A checker may report title-case candidates and terminal full stops, but a human reviews words whose capitalization may be semantically required.

**Automation:** No registry entry or checker currently owns this rule. Terminal full stops are a bounded deterministic finding; capitalization remains advisory evidence.

**Semantic remainder:** Determine which words require capitalization and whether question marks, exclamation marks, colons, or other punctuation are semantically meaningful.

**Exceptions:** Literal titles quoted from an upstream source require an owned reason if preserved verbatim as a heading.

<br />

## Standards rule sections (`MD.STANDARD.SECTION`)
**Classification:** `advisory` with a deterministic canonical-structure subset.

**Scope:** Every H2 section in maintained standards listed by [`README.md`](README.md), outside fenced examples and generated or backup files.

Render every in-scope rule-owner section using this canonical source order:
````markdown
## Rule title (`DOMAIN.RULE.ID`)
**Classification:** `deterministic`, `advisory`, or `semantic-only`, with bounded qualifiers when needed.

**Scope:** Exact applicability and normal exclusions.

Normative requirements and necessary examples.

**Automation:** Checker or evidence producer, coverage relation, and implementation status, or an explicit statement that none exists.

**Semantic remainder:** Required human judgment and evidence, or `None`.

**Exceptions:** Permitted deviation and approval mechanism, or `None`.
````

The rule heading contains each owned rule ID in inline code. `**Classification:**` follows the heading directly. Each required field appears exactly once in the displayed order; one empty line separates fields and the normative body. `**Automation:**`, `**Semantic remainder:**`, and `**Exceptions:**` are the final three labeled fields. The normative body may contain multiple paragraphs, tight lists, tables, fences, or lower-level subsections when the rule needs them.

[`governance.md`](governance.md) owns what the fields mean and when a provision needs an independent rule ID. This rule owns only their canonical Markdown representation. Maintained owner documents reserve H2 for rule-owner sections. Use the document introduction or an H3-or-lower subsection for non-normative explanation; indexes and literal fenced examples do not become rule-owner sections merely because they use headings.

**Automation:** `dev/audit/markdown_policy.py` examines every H2 in maintained `docs/standards/*.md` files, rejects an H2 without a rule ID, and reports missing, duplicate, or out-of-order fields plus a nonadjacent classification field. `dev/audit/standards_registry.py` independently enumerates the same owner surface. Focused accepted and rejected fixtures exercise the deterministic structure. Both checks remain nonblocking and do not decide whether field contents are semantically accurate.

**Semantic remainder:** Determine whether a section is normative, whether it owns the correct rule IDs, and whether each field accurately describes the rule's real boundary and implementation state.

**Exceptions:** Literal fenced examples are outside scope. An H2 rule-owner section may not omit a required field; unavailable automation or exceptions are recorded explicitly rather than represented by omission.

<br />

## Explicit HTML heading anchors (`MD.ANCHOR.CANONICAL`)
**Classification:** `deterministic` for anchor source, binding, uniqueness, and identified explicit-ID link resolution, with a semantic inbound-link review remainder.

**Scope:** Standalone explicit HTML heading-anchor candidates associated with formal H1-H6 or recognized structural informal H7/H8 headings, plus maintained local fragment links identified as targeting those explicit IDs, outside fenced code.

Write one explicit heading anchor in the exact standalone source form `<a id="ID"></a>`. After standard HTML-entity decoding, the `ID` value is nonempty and contains no ASCII whitespace. Put the anchor on the source line immediately before exactly one heading, with no intervening blank line, `<br />`, anchor, comment, or other content. A heading has at most one explicit anchor. Consecutive, dangling, malformed, and redundant anchors are invalid.

Decoded explicit IDs are unique within a file. A repeated ID is invalid whether it appears before equivalent or distinct heading text. When one repeated ID belongs to distinct headings, inventory every maintained inbound link, choose the reviewed stable ID for each heading, and update affected links in the same change. Tooling must not guess or synthesize that replacement.

A maintained local fragment link identified as targeting an explicit ID must resolve to exactly one decoded explicit ID in its target file. Resolve a relative target path from the source file, URI-decode the fragment once as valid UTF-8, HTML-decode the explicit ID, and compare them exactly and case-sensitively. Do not case-fold, slugify, or silently repair either value.

A fragment that targets a renderer-generated heading identifier is outside this rule even when the target file also contains explicit anchors. The presence of any explicit anchor in a file does not convert every fragment link into an explicit-ID link, and the checker must not report an unmatched fragment merely on that basis.

Before renaming or removing an explicit ID, perform a reviewed repository-wide inbound-link inventory. Update every link identified as targeting that explicit ID in the same change and retain the inventory as evidence. Preserve every otherwise valid stable ID and its links; do not automatically regenerate IDs from heading text.

For `MD.SECTION.BREAK`, a standalone, immediately adjacent anchor candidate remains structurally transparent even when this rule finds its source noncanonical. The anchor defect is reported here without a consequential section-boundary finding solely because of the defect.

**Automation:** `dev/audit/markdown_policy.py` excludes fenced literals; checks canonical anchor source, immediate heading binding, one-anchor binding, and decoded per-file uniqueness; and requires exactly one match for fragments that decode to a present explicit ID. It does not treat every unmatched fragment as an anchor error or implement a broad generated-heading link checker. Reviewed inbound-link inventories supply the semantic evidence when explicit IDs are renamed, removed, or reused by distinct headings.

`dev/tools/markdown_format.py` may remove an earlier exact duplicate in the same heading-boundary cluster and place the single retained anchor immediately before its heading. It never invents or renames an ID or rewrites a link.

**Semantic remainder:** Identify inbound links whose intended explicit target cannot be inferred solely from current exact fragment equality, especially during ID renames or removal, and distinguish them from renderer-generated heading fragments. Record the disposition for every explicit ID changed by a migration.

**Exceptions:** Fenced literal examples are outside scope. There is no exception for duplicate, dangling, consecutive, detached, or malformed maintained anchors.

<br />

## Section separators (`MD.SECTION.BREAK`)
**Classification:** `deterministic` for a recognized heading unit, with a semantic remainder for ambiguous informal-heading intent.

**Scope:** Boundaries before formal H1-H6 and recognized structural informal H7/H8 sections outside fenced code, including headings inside and immediately after a `<details>` element.

For section-boundary classification, an optional explicit HTML heading-anchor candidate and the heading on its immediately following source line form one structural heading unit. Evaluate the boundary exactly as if the heading appeared at the anchor candidate's position. Adding, removing, or making that immediately adjacent candidate noncanonical neither creates nor removes a required separator or an exemption. `MD.ANCHOR.CANONICAL` independently decides whether the candidate is valid.

A contiguous cluster of standalone anchor candidates immediately attached to one heading is likewise transparent to boundary classification. Consecutive, redundant, malformed, or duplicate candidates are reported under `MD.ANCHOR.CANONICAL`; they do not cause a consequential `MD.SECTION.BREAK` finding solely because they are invalid. A line containing genuine non-anchor content is not transparent and may change the boundary classification.

Do not put a separator before the first headed entry in a file; `MD.FILE.BOUNDARY` supplies its leading source boundary. Every later recognized heading unit has exactly one of these canonical boundaries:
1. **Ordinary later section.** Put the exact XHTML break `<br />` on its own line with exactly one empty line before it and one empty line between it and the heading unit.
2. **Contentless direct parent and child.** When one heading is followed by a heading exactly one level deeper, with no real source content between the parent heading line and the child heading unit, put no separator and no blank line between them.
3. **Section after a details close.** When `</details>` is immediately followed by `<br />`, followed by exactly one empty line and the next heading unit, that direct break is the complete boundary. Do not add a second separator before the heading unit.

The contentless-parent form applies uniformly, including H1→H2, and only to an immediate level increase of one. An empty sibling, a shallower heading, or a heading that skips a level uses the ordinary later-section boundary. If prose, a list, a table, a fence, an HTML block, or any other real content intervenes after a parent heading, the following heading uses the ordinary boundary regardless of its level or anchor presence.

Do not substitute `<br>`, `<br/>`, multiple breaks, a horizontal rule, or extra empty lines. A separator is invalid where the contentless direct-parent/child form requires none, and a second separator is invalid after the direct `</details>` plus `<br />` form.

Every formal H1-H6 and safely recognized structural H7/H8 heading follows these same deterministic boundary classifications. Boundary conformance must not be used as the prerequisite for recognizing a safely identifiable informal heading.

**Automation:** `dev/audit/markdown_policy.py` provides exact deterministic boundary coverage for formal H1-H6 and safely recognized structural H7/H8 heading units. It uses a fence-aware heading-unit model, keeps boundary results invariant under canonical or immediately adjacent noncanonical anchor candidates, and records ambiguous emphasized-only lines for semantic review. The registered coverage relation is `subset`.

**Semantic remainder:** Decide whether an otherwise ambiguous emphasized-only line was intended as H7/H8. If it was, review it against the same deterministic boundary forms; ambiguity about intent does not weaken the normative boundary.

**Exceptions:** Fenced literal examples are outside scope. The three canonical forms are classifications, not exceptions, and recognized H7/H8 headings receive no boundary exemption.

<br />

## Natural prose (`MD.PROSE.UNWRAP`)
**Classification:** `advisory` with a bounded deterministic formatter subset.

**Scope:** Ordinary prose outside fences, tables, HTML blocks, block quotes, list structure, link definitions, and explicit hard breaks.

Keep each natural prose paragraph on one source line and let the renderer or editor wrap it. Do not hard-wrap prose to a fixed column. Preserve paragraph boundaries, explicit hard breaks, block quotes, list continuation prefixes, link definitions, raw HTML, tables, and literal code.

Automatic unwrapping is advisory until a parser proves that joining source lines preserves block and inline meaning. Never join across a block boundary or an explicit hard break.

**Automation:** `MD.PROSE.UNWRAP` is registered with subset coverage. `format_document()` and `--mode proposed` identify candidate paragraphs in preview-only mode; the CLI rejects combining that mode with `--write`.

**Semantic remainder:** Review ambiguous inline HTML, reference definitions, container boundaries, and author-intended hard breaks.

**Exceptions:** Literal source demonstrations and immutable historical evidence remain unchanged when rewriting them would alter their evidentiary meaning.

<br />

## Colon-introduced structures (`MD.COLON.STRUCTURE`)
**Classification:** `deterministic`.

**Scope:** A prose line ending in a colon whose immediately associated block is a list, fenced code block, or GFM table.

Put no blank line between the colon-introduced prose and an immediately associated list or fenced code block. Put exactly one blank line between the colon-introduced prose and a GFM table.

Canonical examples:
`````markdown
Items:
- First item.
- Second item.

Example:
```text
literal content
```

Values:

| Name | Value |
| :--- | ---:  |
| A    | 1     |
`````

The distinction is source-canonical and parser-aware: lists and fenced code can interrupt the introducing paragraph, while the repository requires an explicit table boundary.

**Automation:** `dev/audit/markdown_policy.py` checks no blank before immediately associated lists and fences and one blank before immediately associated tables. Focused fixtures establish `subset` coverage for recognized block forms.

**Semantic remainder:** `None` once the associated block is identified. Determining whether prose actually introduces a later, nonadjacent block requires review.

**Exceptions:** None for an immediately associated block.

<br />

## Structural indentation (`MD.INDENT.SPACES`)
**Classification:** `deterministic`.

**Scope:** Leading structural indentation outside fenced code.

Use spaces, never tabs, for structural indentation. Preserve tabs that are literal fenced content. Container-aware tooling must distinguish indentation from literal text before offering a mutation.

**Automation:** `dev/audit/markdown_policy.py` and its registry entry check leading tabs outside fenced literal content with `subset` coverage.

**Semantic remainder:** `None` after fenced and literal contexts are excluded.

**Exceptions:** None for maintained Markdown structure.

<br />

## Lists (`MD.LIST.SPACING`, `MD.LIST.INDENT`)
**Classification:** `advisory` with deterministic simple-list portions.

**Scope:** Ordered and unordered Markdown lists outside fenced code.

Use two-space repository nesting levels when that indentation produces the intended CommonMark container structure. Keep ordinary single-paragraph lists tight: put no blank lines between peer items at any nesting depth or between a parent item and its nested list. Use blank lines only when an item intentionally contains multiple block elements and the list must therefore be loose.

Ordered lists retain decimal markers. List parsing ultimately follows marker width and content indentation under CommonMark. The two-space convention must yield to the indentation needed for multi-digit ordered markers, embedded blocks, and other cases where fixed indentation would change the container structure.

**Automation:** `MD.LIST.SPACING` and `MD.LIST.INDENT` are registered with `subset` coverage for recognized simple lists. Focused fixtures check tight peers, tight nested lists, two-space levels, tabs, and noncanonical simple indentation without deciding complex CommonMark containers.

**Semantic remainder:** Review intentional loose lists, multi-digit ordered containers, embedded blocks, and ambiguous indentation.

**Exceptions:** A literal example of noncanonical Markdown belongs in a fence or rejected fixture rather than an inline suppression.

<br />

## Unordered-list markers (`MD.LIST.MARKER`)
**Classification:** `advisory` with deterministic simple-list portions.

**Scope:** Unordered-list markers outside fenced code.

Use `-` at the top level and alternate markers by depth: `+` for children, `-` for grandchildren, `+` for great-grandchildren, and so on. Do not use `*` as an unordered-list marker. Marker alternation follows structural depth even when an unordered list is nested inside an ordered list:
```markdown
- Top-level item.
  + Child item.
    - Grandchild item.
      + Great-grandchild item.
- Next top-level item.
```

**Automation:** `MD.LIST.MARKER` is new and has no registry entry or checker. The current list parser recognizes all three unordered markers but does not enforce alternating depth.

**Semantic remainder:** Review marker depth inside complex ordered containers and embedded block structures.

**Exceptions:** A literal example of `*` or nonalternating markers belongs in a fence or rejected fixture.

<br />

## Details hierarchy (`MD.DETAILS.HIERARCHY`)
**Classification:** `advisory`; rebasing remains deferred.

**Scope:** Formal and informal heading structure inside HTML `<details>` elements that contain a `<summary>`.

The shallowest heading in a `<details>` body is one level deeper than the nearest preceding parent heading. Preserve relative internal hierarchy and ignore fenced code and nested `<details>` bodies while evaluating the outer body. Apply `MD.SECTION.BREAK` between headed sections inside the body; `<summary>` is not a heading.

When a valid relative hierarchy extends below H6, represent level seven with the H7 form from `MD.HEADING.INFORMAL` and level eight with the H8 form. Do not emit seven- or eight-marker ATX-like lines. If rebasing would require a level below H8, report a semantic hierarchy decision rather than inventing another rendering convention.

**Automation:** `MD.DETAILS.HIERARCHY` is registered to the formatter unit suite. Proposed H1-H6 rebasing exists only in preview mode, while the deterministic formatter test proves that heading levels inside `<details>` remain unchanged. Informal H7/H8 output, nesting, and overflow behavior remain unproved and nonblocking.

**Semantic remainder:** Choose whether a hierarchy should use informal H7/H8, flatten, split, move, or otherwise be redesigned, especially when it would overflow H8.

**Exceptions:** None permit seven- or eight-marker ATX-like source or an invented level below H8.

<br />

## Fenced code (`MD.FENCE.CLOSED`)
**Classification:** `deterministic`.

**Scope:** Backtick and tilde fenced code blocks.

Every opening fence has a closing fence made from the same character and at least the opening run length. Fence content is literal and formatters do not rewrite it. Use a fence longer than any same-character fence that must appear literally inside it.

**Automation:** `MD.FENCE.CLOSED` is registered with subset coverage. `fence_errors()` reports unclosed same-character fences while preserving literal bodies.

**Semantic remainder:** `None` for closure once container context is parsed.

**Exceptions:** Rejected fixtures intentionally contain invalid fences and remain outside repository conformance scope.

<br />

## Fence style (`MD.FENCE.STYLE`)
**Classification:** `advisory` with deterministic marker and info-string opportunities.

**Scope:** Fenced code blocks that are valid under `MD.FENCE.CLOSED`.

Use backtick fences by default. A tilde fence is permitted only when the subject is tilde-fence behavior or when literal backtick runs make a tilde fence materially clearer and safer than a longer backtick fence. Prefer a language info string whenever the content has a known language; use `text` for plain literal output and omit the info string only when no accurate label exists.

**Automation:** No registry entry or checker currently owns this style rule. A future checker can decide marker spelling and missing known info strings for bounded contexts, while unclear language classification remains advisory.

**Semantic remainder:** Determine whether a tilde fence materially improves safety and which info string accurately describes mixed or unusual content.

**Exceptions:** Markdown fixtures may use the fence form they are specifically testing.

<br />

## GFM tables (`MD.TABLE.CANONICAL`)
**Classification:** `deterministic` for supported non-emoji tables.

**Scope:** GFM pipe tables with a header row, delimiter row, and equal logical column counts.

Use leading and trailing pipes. Use `:---`, `---:`, and `:---:` for left, right, and center alignment, then pad cells for readable canonical source alignment. Escaped pipes and pipes inside matched inline-code spans are cell content rather than separators. Put exactly one blank line before a table when it follows a heading or colon-introduced sentence.

Canonical alignment uses a deterministic display-width approximation for combining marks and East Asian wide or full-width code points. Emoji are prohibited in table cells because editor and renderer grapheme widths are not stable enough for canonical source padding. This prohibition does not apply to ordinary prose. If literal tabular source data contains emoji, preserve it in a fenced representation or obtain an explicit exception rather than deleting or silently replacing data.

**Automation:** `MD.TABLE.CANONICAL` is registered with subset coverage. The current parser and formatter canonicalize supported pipe tables but do not reject emoji or all malformed table candidates.

**Semantic remainder:** `None` for supported table source. Identifying whether malformed pipe-like text was intended as a table may require review.

**Exceptions:** A table that must preserve literal emoji data requires an owned exception and is not eligible for automatic source padding.

<br />

## Delimited-text conversion (`MD.TABLE.CONVERT`)
**Classification:** `deterministic`.

**Scope:** Explicit comma-separated or tab-separated input requested for conversion to a GFM table.

Conversion requires an explicit comma or tab delimiter, an explicit header row, and equal column counts. Do not guess between comma, tab, or runs of spaces. Conversion produces a left-aligned canonical GFM table and never mutates a file unless the formatter is invoked with an explicit write option. Input containing emoji must be reported rather than converted into a table that violates `MD.TABLE.CANONICAL`.

**Automation:** `MD.TABLE.CONVERT` is registered with deterministic unit-test coverage for explicit delimiters and equal-width rows. Emoji rejection is not yet implemented.

**Semantic remainder:** `None` for validated input.

**Exceptions:** None permit delimiter guessing or silent data loss.

<br />

## Tooling and scope (`MD.TOOLING.EXPLICIT`)
**Classification:** `deterministic` for tool mutation boundaries; policy rollout follows governance review.

**Scope:** `dev/audit/markdown_policy.py`, `dev/tools/markdown_format.py`, Markdown fixtures, registry entries, and repository-wide Markdown operations.

`python -m dev.audit.markdown_policy` is non-mutating and fails by default only on deterministic findings. `python -m dev.tools.markdown_format` provides pure parser and formatter functions. Its default `deterministic` mode writes only with explicit `--write`; `--mode proposed` previews deferred prose-unwrapping and details-rebasing behavior and rejects `--write`. Formatter and checker primitives remain independent of editor integrations. Generated reports live under `artifacts/tests/`.

Deterministic formatting derives separator placement from the same fence-aware heading-unit and boundary classification used by the checker. It may collapse a redundant second separator, remove a separator prohibited by the contentless direct-parent/child form, remove an earlier exact duplicate anchor in the same heading-boundary cluster, and place the single retained anchor immediately before its heading.

Boundary formatting treats an immediately adjacent anchor candidate as structurally transparent even when the anchor is noncanonical. The formatter may perform independently safe boundary normalization, but it does not silently repair malformed anchor syntax, invent or rename an ID, rewrite a link, or decide whether an ambiguous emphasized-only line is an informal heading.

A duplicate ID attached to distinct headings remains a checker finding until reviewed IDs and every affected maintained inbound link are supplied.

Every CLI mutation requires `--write`. Explicit paths constrain the target set; with `--root` and no paths, the formatter discovers the same maintained scope as the checker. Preview mode prints paths that would change and returns nonzero when changes are pending; it does not modify source.

After an authoritative Markdown policy changes, update the execution registry, checker, formatter, accepted and rejected fixtures, formatter fixtures, and repository contract before performing a repository-wide write. Until that migration reconciles, checker results that encode the superseded policy are drift evidence rather than authority for the new rule.

Temporary `bak.*.md` files may exist during local review but are not maintained standards. Move them to the private repository or otherwise remove them from the public checker scope before repository-wide conformance or commit review.

The ignored `tmp/` scripts remain design evidence until every behavior has a recorded disposition, the adopted subset has idempotence evidence, and deletion is separately authorized. Noncanonical rejected fixtures and formatter inputs remain excluded from repository conformance and are exercised by unit tests.

**Automation:** `MD.TOOLING.EXPLICIT` is registered with `subset` coverage to the focused unit suite; repository enforcement uses `python -m dev.audit.markdown_policy --root .`. The checker and deterministic formatter share fence-aware heading-unit and boundary classification, implement the approved formal and safely recognized informal boundary forms, keep immediately adjacent anchor candidates boundary-transparent, check the deterministic anchor subset and identified explicit-ID links, and retain approved file-boundary, colon-structure, supported-table, structural-indentation, and simple-list behavior. Unit and contract tests prove repository no-write preview, deterministic idempotence, literal-fence preservation, unchanged details heading levels, no-guess ID behavior, and the preview-only write guard. Ambiguous informal-heading intent, inbound-link intent during ID migration, heading case, list markers, fence style, table emoji policy, complex CommonMark containers, prose unwrapping, and details rebasing remain nonblocking.

**Semantic remainder:** Review ambiguous informal-heading and explicit-ID inbound-link intent and separately decide the deferred prose, details, marker, emoji, and complex-container behaviors.

**Exceptions:** No tool may mutate maintained Markdown without an explicit write operation and reviewed scope.
