#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: make.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


# Require Bash >= 4.4 before doing any work.
if [[ -z "${BASH_VERSION:-}" ]]; then
    echo "error(make.sh):" \
        "this script must be run under Bash >= 4.4." >&2
    exit 1
elif ((
    BASH_VERSINFO[0] < 4 || ( BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] < 4 )
)); then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "this script requires Bash >= 4.4; current version is" \
        "'${BASH_VERSION}'." >&2
    exit 1
fi

# Run in safe mode, exiting on errors, unset variables, and pipe failures.
set -euo pipefail


# Resolve paths relative to 'tests/fixtures'.
dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
dir_fix="${dir_scr}"

# Source shared fixture-generation helpers.
# shellcheck source=tests/support/fixture_helpers.sh
source "${dir_scr}/../../support/fixture_helpers.sh"

# Declare every generated path up front. The directory names the verdict the
# checker must return for the documents inside it, so a verdict with no
# fixture is visible as an absent directory rather than as a filename nobody
# wrote. 'format/' holds input and expected pairs for the formatter, which
# makes a claim about a rewrite rather than about a verdict.
dir_acc="${dir_fix}/accepted"
dir_bnd="${dir_fix}/boundary"
dir_nap="${dir_fix}/non_applicable"
dir_rej="${dir_fix}/rejected"
dir_fmt="${dir_fix}/format"

fil_acc_anchors="${dir_acc}/anchors.md"
fil_acc_basic="${dir_acc}/basic.md"
fil_acc_quotes="${dir_acc}/blockquotes.md"
fil_acc_fences="${dir_acc}/nested_fences.md"
fil_acc_sections="${dir_acc}/section_boundaries.md"
fil_acc_standard="${dir_acc}/standard_section.md"

fil_bnd_quotes="${dir_bnd}/blockquotes.md"
fil_bnd_fences="${dir_bnd}/nested_fences.md"

fil_nap_quotes="${dir_nap}/blockquotes.md"
fil_nap_fences="${dir_nap}/nested_fences.md"

fil_rej_anchors="${dir_rej}/anchors.md"
fil_rej_quotes="${dir_rej}/blockquotes.md"
fil_rej_fences="${dir_rej}/nested_fences.md"
fil_rej_sections="${dir_rej}/section_boundaries.md"
fil_rej_spacing="${dir_rej}/spacing.md"
fil_rej_standard="${dir_rej}/standard_section.md"
fil_rej_unclosed="${dir_rej}/unclosed_fence.md"

fil_fmt_input="${dir_fmt}/input.md"
fil_fmt_expected="${dir_fmt}/expected.md"
fil_fmt_quotes_input="${dir_fmt}/blockquotes_input.md"
fil_fmt_quotes_expected="${dir_fmt}/blockquotes_expected.md"
fil_fmt_fences_input="${dir_fmt}/nested_fences_input.md"
fil_fmt_fences_expected="${dir_fmt}/nested_fences_expected.md"


# Remove stale outputs so regeneration is idempotent.
rm_files "${dir_fix}" \
    "${fil_acc_anchors}" \
    "${fil_acc_basic}" \
    "${fil_acc_quotes}" \
    "${fil_acc_fences}" \
    "${fil_acc_sections}" \
    "${fil_acc_standard}" \
    "${fil_bnd_quotes}" \
    "${fil_bnd_fences}" \
    "${fil_nap_quotes}" \
    "${fil_nap_fences}" \
    "${fil_rej_anchors}" \
    "${fil_rej_quotes}" \
    "${fil_rej_fences}" \
    "${fil_rej_sections}" \
    "${fil_rej_spacing}" \
    "${fil_rej_standard}" \
    "${fil_rej_unclosed}" \
    "${fil_fmt_input}" \
    "${fil_fmt_expected}" \
    "${fil_fmt_quotes_input}" \
    "${fil_fmt_quotes_expected}" \
    "${fil_fmt_fences_input}" \
    "${fil_fmt_fences_expected}"

mkdirs "${dir_acc}" "${dir_bnd}" "${dir_nap}" "${dir_rej}" "${dir_fmt}"


# Author every document literally. The delimiter is quoted, so no '$' or
# backtick in a fixture body reaches the shell, and a fixture that is
# deliberately malformed Markdown stays malformed exactly as written.

# Accepted: documents the checker must report nothing about
cat << 'EOM' > "${fil_acc_anchors}"

# Accepted explicit-anchor fixture
The [explicit target](#explicit-target) resolves exactly once, while the [generated heading target](#generated-heading-target) remains outside `MD.ANCHOR.CANONICAL`.

<br />

<a id="explicit-target"></a>
## Explicit target
The stable explicit ID remains attached to its heading.

<br />

## Generated heading target
This heading intentionally has no explicit anchor.
EOM

cat << 'EOM' > "${fil_acc_basic}"

# Accepted fixture
Introductory prose remains one natural source line.

Items:
- First top-level item.
- Second top-level item.
  - Nested item.
  - Nested peer.

| Name | Value |
| :--- | :---  |
| A    | 界    |

<br />

## Heading-attached table
| Name | Value |
| :--- | :---  |
| B    | C     |

Colon-introduced table:

| Name | Value |
| :--- | :---  |
| D    | E     |

Quoted note:
> The blockquote remains attached.
EOM

cat << 'EOM' > "${fil_acc_quotes}"

# Blockquote acceptance
> First quoted paragraph.
>
> Second quoted paragraph.
>
> > Nested quoted paragraph.
EOM

cat << 'EOM' > "${fil_acc_fences}"

# Nested fence acceptance
````````````markdown
Outer literal content.

`````````text
Shorter same-character run stays literal.
`````````

~~~
Other-character run stays literal.
~~~
````````````
EOM

cat << 'EOM' > "${fil_acc_sections}"

# Accepted section-boundary matrix
The first section contains content.

<br />

## Ordinary section without an anchor
This section contains content.

<br />

<a id="ordinary-with-anchor"></a>
## Ordinary section with an anchor
This section contains content.

<br />

## Contentless parent without anchors
### Direct child without anchors
This child contains content.

<br />

<a id="parent-with-anchor"></a>
## Contentless parent with an anchor
<a id="child-with-anchor"></a>
### Direct child with an anchor
This child contains content.

<br />

## Parent containing content
Real parent content preserves the ordinary separator before its child.

<br />

<a id="child-after-content"></a>
### Child after content
This child contains content.

<br />

## Empty sibling

<br />

<a id="empty-sibling-with-anchor"></a>
## Empty sibling with an anchor

<br />

#### Skipped deeper section after an ordinary boundary
This section contains content.

<br />

## Shallower section
<details>
<summary>Details boundary</summary>
Content.
</details>
<br />

<a id="after-details"></a>
## Section after a details close
The direct break after the close is the complete boundary.
EOM

cat << 'EOM' > "${fil_acc_standard}"

# Fixture standard
This file demonstrates one canonical rule-owner section.

<br />

## Example rule (`EXAMPLE.RULE`)
**Classification:** `deterministic`.

**Scope:** Maintained example input.

The example must retain its canonical structure.

**Automation:** A focused fixture proves the source shape.

**Semantic remainder:** `None`.

**Exceptions:** None.
EOM

# Boundary: documents that sit on a decision edge and must still pass, and
# non-applicable: documents the rule does not reach at all.
cat << 'EOM' > "${fil_bnd_quotes}"

# Attached heading
> First quoted paragraph.

Introduced:
> Second quoted paragraph.
EOM

cat << 'EOM' > "${fil_bnd_fences}"

# Nested fence boundary
~~~~~~~~~~~~text
````````````
An other-character run of equal length remains literal.
````````````
~~~~~~~~~~~~
EOM

cat << 'EOM' > "${fil_nap_quotes}"

# Literal blockquote
```text
>Literal fenced content.
>>Literal nested content.
```
EOM

cat << 'EOM' > "${fil_nap_fences}"

# Non-applicable fence-like prose
An inline ````````` run is not a fenced block.
EOM

# Rejected: each document violates the rules named for it in
# 'tests/unit/dev_audit/test_markdown_policy.py' and no others.
cat << 'EOM' > "${fil_rej_anchors}"

# Rejected explicit-anchor fixture
This fixture contains independent anchor defects.

<br />

<a name="malformed"></a>
## Malformed candidate
Content.

<br />

<a id="blank-gap"></a>

## Blank gap
Content.

<br />

<a id="break-gap"></a>
<br />

## Break gap
Content.

<br />

<a id="consecutive"></a>
<a id="consecutive"></a>
## Consecutive and duplicate candidates
Content.

<br />

<a id="distinct-heading-duplicate"></a>
## First distinct heading
Content.

<br />

<a id="distinct-heading-duplicate"></a>
## Second distinct heading
Content.

<br />

<a id="dangling"></a>
EOM

cat << 'EOM' > "${fil_rej_quotes}"

# Blockquote rejection
>First quoted paragraph.
> 
>>Nested quoted paragraph.
EOM

cat << 'EOM' > "${fil_rej_fences}"

# Nested fence rejection
````````````markdown
The outer fence is not closed.

`````````text
A shorter same-character run cannot close it.
`````````
EOM

cat << 'EOM' > "${fil_rej_sections}"

# Rejected section-boundary matrix
Content before a missing ordinary separator.
## Ordinary section without an anchor
Content before another missing ordinary separator.
<a id="ordinary-with-anchor"></a>
## Ordinary section with an anchor

<br />

## Contentless parent without an anchor

<br />

### Direct child with a prohibited separator

<br />

<a id="parent-with-anchor"></a>
## Contentless parent with an anchor

<br />

<a id="child-with-anchor"></a>
### Anchored direct child with a prohibited separator
Content before a doubled ordinary separator.

<br />

<br />

## Ordinary section after doubled breaks
<details>
<summary>Details boundary</summary>
Content.
</details>
<br />

<br />

<a id="after-details"></a>
## Section after a prohibited second details break
EOM

cat << 'EOM' > "${fil_rej_spacing}"
# Rejected spacing fixture

Text after an invalid heading gap.

Items:

- First item.

- Second item.

<br />

## Table with an invalid heading gap

| Name | Value |
| :--- | :---  |
| A    | B     |

Table without its required colon gap:
| Name | Value |
| :--- | :---  |
| C    | D     |

Quoted note:

> Invalid detached blockquote.
EOM

cat << 'EOM' > "${fil_rej_standard}"

# Rejected fixture standard
This file demonstrates noncanonical rule-owner sections.

<br />

## Out-of-order rule (`EXAMPLE.ORDER`)
**Scope:** Maintained example input.

**Classification:** `deterministic`.

The example has every field in the wrong order.

**Automation:** A focused fixture reports the source shape.

**Semantic remainder:** `None`.

**Exceptions:** None.

<br />

## Missing-field rule (`EXAMPLE.MISSING`)
**Classification:** `deterministic`.

The example omits its scope.

**Automation:** A focused fixture reports the source shape.

**Semantic remainder:** `None`.

**Exceptions:** None.
EOM

cat << 'EOM' > "${fil_rej_unclosed}"
# Rejected fence fixture

Example:

```text
missing close
EOM

# Format: input and expected pairs pinning the formatter rewrite. Each expected
# document is also a fixpoint, so the pair proves idempotence as well as the
# rewrite itself.
cat << 'EOM' > "${fil_fmt_input}"
# Formatting fixture
This paragraph is hard
wrapped in source.

Items:
- First.
- Second.

Name | Value
--- | ---:
A | 界

<br />

## Heading table

Name | Value
--- | ---
B | C

Introduced table:
Name | Value
--- | ---
D | E

Quoted note:

> Attached after formatting.
EOM

cat << 'EOM' > "${fil_fmt_expected}"

# Formatting fixture
This paragraph is hard wrapped in source.

Items:
- First.
- Second.

| Name | Value |
| :--- | ---:  |
| A    | 界    |

<br />

## Heading table
| Name | Value |
| :--- | :---  |
| B    | C     |

Introduced table:

| Name | Value |
| :--- | :---  |
| D    | E     |

Quoted note:
> Attached after formatting.
EOM

cat << 'EOM' > "${fil_fmt_quotes_input}"

# Blockquote format
>First.
> 
>>Nested.
EOM

cat << 'EOM' > "${fil_fmt_quotes_expected}"

# Blockquote format
> First.
>
> > Nested.
EOM

cat << 'EOM' > "${fil_fmt_fences_input}"

# Nested format proof

````````````markdown
Literal paragraph one.


Literal paragraph two.
`````````text
## Literal heading

Literal prose
stays split.
`````````
````````````
EOM

cat << 'EOM' > "${fil_fmt_fences_expected}"

# Nested format proof
````````````markdown
Literal paragraph one.


Literal paragraph two.
`````````text
## Literal heading

Literal prose
stays split.
`````````
````````````
EOM


succeed "generated Markdown policy fixtures under ${dir_fix}"
