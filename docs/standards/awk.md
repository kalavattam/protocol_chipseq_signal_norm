
# AWK standard
This document owns maintained `.awk` files: their form, the `gawk` they are written for, their documentation block, and how they are checked and run. [`source_layout.md`](source_layout.md) owns when a program leaves its host for a file of its own (`SOURCE.EMBED.PROGRAMS`), [`shell.md`](shell.md) owns the limits on `awk` programs that Shell carries inline (`SHELL.EMBED.PROGRAMS`), [`help.md`](help.md) owns the shared section vocabulary, and [`source_headers.md`](source_headers.md) owns the header rows.

<br />

## AWK file form (`AWK.FILE.FORM`)
**Classification:** `advisory` with deterministic structural portions.

**Scope:** Every `.awk` file below `lib/awk/`, and any later maintained `.awk` file.

An `.awk` file holds one function's program. It has no shebang, since it is never executed, and carries the source header rows as `#` comments. It is written for the `gawk` that `env_protocol.yml` pins (5.4.1) and is checked with `gawk --lint`. The calling function runs it with `gawk -f`, never `awk -f`, so a missing `gawk` fails loudly.

Comments begin with `# ` and wrap at column 79, and blocks indent four spaces.

**Automation:** No dedicated checker or registry entry implements this rule yet.

**Semantic remainder:** Judge whether a `gawk --lint` warning marks a defect or an accepted AWK idiom, such as an uninitialized element read as 0.

**Exceptions:** None.

<br />

## AWK documentation block (`AWK.HELP.FORM`)
**Classification:** `advisory` with deterministic section-order portions.

**Scope:** The comment block that follows the header of every file in scope of `AWK.FILE.FORM`.

AWK has no heredoc or docstring, so an `.awk` file documents itself in `#` comments, laid out as function help under [`HELP.SECTION.SCHEMA`](help.md#shared-section-schema-helpsectionschema): each section name on its own comment line with its underline row beneath it, and the section's text indented two spaces. The block uses the applicable subsequence of `Usage`, `Parameters`, `Returns`, `Output`, `Notes`, and `Examples`:
- `Usage` gives the `gawk -f` call with placeholders, then what the program does and what it reads.
- `Parameters` gives each variable the caller sets with `-v`, as `name : type` with its description indented beneath.
- `Returns` gives the exit status and the cases that exit nonzero.
- `Output` gives what the program writes to stdout and stderr.
- `Notes` gives what the input must already satisfy and the `Runtime requirements` list.
- `Examples` follows [`HELP.EXAMPLES`](help.md#examples-helpexamples): numbered entries, each with one `'''bash` block.

A fact that a design page or another owner already states is cited, not copied, so the two cannot drift. The block wraps at column 79 like any other comment.

**Automation:** No dedicated checker or registry entry implements this rule yet.

**Semantic remainder:** Decide which optional sections help the caller, and which facts belong in a design page instead.

**Exceptions:** None.
