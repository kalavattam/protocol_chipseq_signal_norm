
# Source-header standard
This document owns the logical source-header contract and approved Bash and Python profiles. Language standards own runtime and source semantics, while this owner defines applicability, field order, canonical rows, widths, years, attribution, and the boundary between header and body.

<br />

## Header structure and profiles (`SOURCE.HEADER.STRUCTURE`)
**Classification:** `advisory` for applicability with a deterministic structural subset.

**Scope:** Maintained changed Bash and Python sources that already bear a shebang, have a repository-owned direct-execution contract, or are explicitly adopted as repository script sources. Generated, vendored, copied third-party, fixture payload, and historical-evidence files are excluded unless explicitly adopted as maintained source.

All in-scope files use this exact base order:
1. approved shebang;
2. `# -*- coding: utf-8 -*-`;
3. comment-only separator `#`;
4. `Script:` basename row;
5. comment-only separator;
6. copyright row;
7. contact email;
8. comment-only separator;
9. the selected profile remainder;
10. two ordinary blank lines before the first source-body construct.

The no-AI profile puts the license row directly after the email separator. It is coherent when no approved applicability evidence establishes AI participation:
```text
<approved shebang>
# -*- coding: utf-8 -*-
#
# Script: exact_basename.ext
#
# Copyright START[-END] by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Distributed under the MIT license.


<source body>
```

The AI-assisted profile inserts one truthful attribution block and one exact separator before the license:
```text
<approved shebang>
# -*- coding: utf-8 -*-
#
# Script: exact_basename.ext
#
# Copyright START[-END] by Kris Alavattam
# Email: kalavattam@gmail.com
#
# <approved OpenAI attribution form>
#
# Distributed under the MIT license.


<source body>
```

Changed-file discovery does not choose between these profiles. Use the AI-assisted profile only for an existing bounded attribution block or explicit applicability evidence; otherwise preserve or establish the no-AI profile. Never fabricate attribution to satisfy source structure.

Python uses `#!/usr/bin/env python3`. Bash uses `#!/usr/bin/env bash`. The POSIX bootstrap `install/scripts/install_envs_entrypoint.sh` uses `#!/bin/sh`. A source-only library or intentionally non-executable submit script may retain the Bash shebang as a language declaration; the shebang does not by itself authorize direct execution. Do not add a shebang to an imported or sourced file merely for visual uniformity.

Parse header fields only inside the bounded opening header. Body comments or literals that resemble header fields do not satisfy or alter the contract. Unknown languages and ambiguous execution roles fail closed for review and are never auto-fixed.

**Automation:** `dev/audit/ai_attribution.py` discovers maintained changed sources and validates the selected no-AI or AI-assisted profile, ordered base fields, exact separators, and the body boundary. Standing discovery validates an existing block and reports attribution-like malformed header content but does not treat an absent block as participation evidence. Explicit applicability supplies the deterministic AI-assisted subset. Applicability and direct-execution intent retain semantic review.

**Semantic remainder:** Determine whether a shebang-less file has a direct-execution contract or an explicit repository script-source role, and whether copied or generated material has been adopted as maintained source.

**Exceptions:** A language or bootstrap exception requires an approved profile naming its paths, shebang, comment syntax, and runtime owner before automated normalization.

<br />

## Script basename (`SOURCE.HEADER.BASENAME`)
**Classification:** `deterministic`.

**Scope:** The `Script:` row of every source in `SOURCE.HEADER.STRUCTURE` scope.

Write `# Script: ` followed by the exact file basename, including its extension. Directory components, aliases, module names, and retired basenames are invalid.

**Automation:** The source-header checker compares the single bounded `Script:` row with `Path.name`. Registry coverage is `subset` because applicability is established by the broader header-profile rule before this comparison runs.

**Semantic remainder:** `None`.

**Exceptions:** None.

<br />

## Header width (`SOURCE.HEADER.WIDTH`)
**Classification:** `deterministic`.

**Scope:** Every physical header line from the shebang through the license row.

Limit each complete physical line, including its comment prefix, to 79 characters. Wrap attribution prose without altering model names or meaning. Do not wrap canonical single-row fields; report an unrepresentable basename or identity for review.

**Automation:** The checker rejects `len(line) > 79`; the normalizer wraps attribution content to fit the same complete-line limit. Positive 79-character and negative 80-character fixtures establish the boundary.

**Semantic remainder:** `None` for measured width.

**Exceptions:** An unavoidably long exact basename requires manual disposition; it is not silently truncated.

<br />

## Copyright years (`SOURCE.HEADER.YEAR`)
**Classification:** `advisory` with deterministic represented-year evidence.

**Scope:** The copyright row of an in-scope source changed in the assessed work.

Preserve the earliest established start year. When a later-year change updates the file, use `START-END` with the explicit assessment year as `END`; when both are the same, use one year. The checker receives an authorized assessment year and does not derive policy from the wall clock.

For an existing header, preserve its start year. For a maintained file receiving its first header, use the earliest tracked commit year when available, otherwise the explicit assessment year.

**Automation:** The checker and normalizer accept `--current-year`, inspect the bounded row, and use Git history only to recover a missing start year. Fixtures cover single years, ranges, preserved starts, and explicit-year overrides.

**Semantic remainder:** Review whether history predates the current path because of a copy, import, or provenance boundary.

**Exceptions:** Third-party copyright is outside automatic normalization and follows the applicability exclusions in `SOURCE.HEADER.STRUCTURE`.

<br />

## AI attribution (`SOURCE.HEADER.AI_ATTRIBUTION`)
**Classification:** `advisory` for participation with a deterministic representation subset.

**Scope:** In-scope maintained sources for which an AI model contributed to development, documentation, or both, plus explicitly declared migration cohorts.

Use one bounded attribution block naming the applicable approved OpenAI tool and the declared model identifiers. Accepted tool forms name OpenAI ChatGPT alone, OpenAI Codex alone, or OpenAI ChatGPT and Codex together. Use a single-tool form when only that tool contributed; do not name both by default. State truthfully whether the contribution concerned development, documentation, or both.

Inside the parentheses, use either an explicit model-token list separated by comma and one space, one `GPT-X-series models` phrase, or the established `GPT-X- and GPT-Y-series models` phrase. A series phrase may be followed by `; most recent: ` and an explicit model-token list. Do not admit arbitrary words, empty elements, stray punctuation, malformed separators, or a tokenless clause. Preserve an established generic-series or explicit-model-list style unless changing that style is part of the authorized work.

Changed-file discovery is evidence that attribution review is needed; it cannot establish model participation by itself. A normalizer may update only explicitly selected maintained files and must be idempotent. When an attribution block is missing, the normalizer requires an explicit `development`, `documentation`, or `both` contribution-domain choice; it reports and leaves the file unchanged when that semantic input is absent. A model becomes required only through an explicit command argument or applicability manifest; no model identifier is a permanent repository-wide requirement.

**Automation:** `dev/audit/ai_attribution.py` recognizes only approved bounded OpenAI forms and model-declaration grammar, parses exact model tokens, observes style, tool set, and contribution domain, and compares them with repeatable `--required-model` values or a path-scoped applicability manifest. Inventory fields distinguish observed from required tools and domains and report agreement; a coherent no-AI profile reports no attribution style or observed attribution data. The standing registry command supplies no required model, tool, or domain. Creating a missing block requires explicit model, tool, and domain inputs. Coverage is `subset`; participation remains semantic.

**Semantic remainder:** Establish which tools and models contributed and whether the activity wording and retained attribution style are accurate.

**Exceptions:** Do not add attribution when the model did not contribute. Record that disposition in migration evidence rather than fabricating a header claim.

<br />

## Unique AI attribution (`SOURCE.HEADER.AI_ATTRIBUTION.UNIQUE`)
**Classification:** `deterministic`.

**Scope:** The bounded attribution block of an in-scope source that declares or requires AI attribution.

Each declared or explicitly required model token appears exactly once. Repeated model tokens, repeated attribution blocks, or a second attribution row inside the bounded opening header are invalid. Attribution-like body comments do not alter header conformance.

**Automation:** The checker uses parsed-token counters only inside the bounded attribution block, so `GPT-5` is not counted inside `GPT-5.6`. It recognizes explicit, generic-series, repeatable required, and future model identifiers. Normalization uses the same counters when adding a required token; an existing duplicate remains a finding rather than being hidden by substring replacement. Positive, exact-duplicate, generic-series, future-token, duplicate-block, and body-lookalike fixtures provide `subset` coverage.

**Semantic remainder:** `None` after applicability is established.

**Exceptions:** None.
