
# Repository standards
This directory is the normative source for maintained repository structure, source, documentation, and testing rules. [`governance.md`](governance.md) owns how rules are identified, reviewed, enforced, excepted, and changed.

| Standard                                       | Ownership                                                                                                               |
| :---                                           | :---                                                                                                                    |
| [`governance.md`](governance.md)               | Rule lifecycle, IDs, precedence, traceability, exceptions, and coverage                                                 |
| [`repository_layout.md`](repository_layout.md) | Directory ownership, package boundaries, and prohibited locations                                                       |
| [`markdown.md`](markdown.md)                   | Golden Markdown source for every maintained Markdown file                                                               |
| [`shell.md`](shell.md)                         | Bash runtime, source form, static analysis, sourcing, wrapper topology, submit bootstrap, and shell-test form           |
| [`help.md`](help.md)                           | Shared help schema, rendered/source semantics, audiences, aliases, examples, Runtime requirements, and parameters       |
| [`python.md`](python.md)                       | Python package, import, API, CLI, dependency, installation, version-floor, docstring, Ruff, and Python-test obligations |
| [`testing.md`](testing.md)                     | Test layers, safe runner, gates, fixtures, evidence, reporting, cleanup, and proportional proof                         |
| [`source_headers.md`](source_headers.md)       | Header applicability, profiles, field order, basenames, widths, years, attribution, and body boundary                   |

Domain standards add only domain-specific semantics. Every listed standard remains subject to [`markdown.md`](markdown.md). Temporary `bak.*.md` snapshots are review evidence, not maintained standards, and are excluded from conformance inventories.
