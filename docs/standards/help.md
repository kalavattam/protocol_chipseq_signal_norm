
# Help and documentation standard
This document owns shared help and documentation meaning across maintained implementation languages. Python and Shell are the active automated realizations. R and Rust inherit the shared architecture but remain dormant, unenforced, and unmigrated until their existing language-owner applicability triggers are satisfied. Language owners retain source syntax and language-native rendering.

<br />

## Shared section schema (`HELP.SECTION.SCHEMA`)
**Classification:** `advisory` with deterministic vocabulary and ordering portions.

**Scope:** Maintained command, helper, and callable documentation when the same documentation concepts apply.

Use NumPy/SciPy-style section names and underline rows when a rendered format uses headed sections. The shared vocabulary is `Description`, `Usage`, `Parameters`, `Expected globals`, `Generated globals`, `Returns`, `Yields`, `Receives`, `Raises`, `Warns`, `Warnings`, `Output`, `See Also`, `Notes`, `References`, `Examples`, `Attributes`, and `Methods`.

Use `Usage` for invocation shape, `Parameters` for positional arguments, options, and function inputs, `Returns` for shell status or a callable return value, and `Output` for stdout, stderr, or generated files. `Output` and `Returns` are not synonyms.

Function help uses the applicable subsequence of:
```text
Usage
-----
Parameters
----------
Expected globals
----------------
Generated globals
-----------------
Returns
-------
Output
------
Notes
-----
Examples
--------
```

Every recognized function-help document requires `Usage` and `Returns`. `Parameters` is required when the function accepts arguments or options. A zero-argument, non-CLI helper may omit other inapplicable sections. Script help does not acquire `Returns` merely because function help uses it. `Examples` is final when present. Do not repeat a top-level section in one rendered document.

Do not create user-facing `Testing`, `TODO`, `Dependency`, `Dependencies`, `Argument`, or `Arguments` sections. Put maintainer work in code or developer documentation, dependency facts under `Notes`, and all inputs in one `Parameters` section.

**Automation:** Help structure and rendered-output contracts check recognized headings, order, uniqueness, and bounded applicability with `subset` coverage.

**Semantic remainder:** Decide which optional sections materially help the audience and whether prose belongs in a shared or domain-specific section.

**Exceptions:** A domain owner may add a section only when the shared vocabulary cannot express a stable concept and the new name is documented with its audience and ordering.

<br />

## Source help structure (`HELP.SOURCE_STYLE`)
**Classification:** `advisory` with deterministic indentation and entry-shape portions.

**Scope:** Heredoc-based shell help and local shell-function documentation.

Top-level headings and underline rows begin at column zero. Content begins two spaces below the section, descriptions two spaces below their entry, and nested lists add two spaces per level. Avoid hard tabs. Preserve physical newlines deliberately because each heredoc newline affects rendered output.

In `Usage`, put the command or function on its own first invocation line and indent argument continuations by two additional spaces:
```text
Usage
-----
  filter_alignment_sc
    [--help]
    [--threads <int>]
    --fil_in <file> [--ref_fa <file>]
    --fil_out <file>
    [--mito] [--chk_chr]
```

Use square brackets as the primary optionality signal. Unconditional required arguments are unbracketed; optional, defaulted, conditional, and mode-specific arguments are bracketed and their descriptions state defaults or conditions directly. Parenthesize required mutually exclusive branches and bracket the whole group when the choice itself is optional.

Use the same argument order in `Usage` and `Parameters`. [`HELP.OPTION.ORDER`](#semantic-option-order-helpoptionorder) owns shared semantic ordering and relationships. This owner retains Shell heredoc indentation and source realization only.

Parameter and global entries use `name : type` followed by an indented description. Positional entries use `ordinal  name : type`; zero-pad ordinals when there are ten or more arguments. Separate logical entries with one blank line. Keep aliases for one option on one row, separated by comma and one space, followed by one shared type.

Inside `Notes`, a pseudo-heading such as `Runtime requirements:` begins at normal section indentation and its content begins two spaces deeper. `Runtime requirements:` is invalid outside the active top-level `Notes` section, even when its indentation and entries are otherwise correct. Structured commands, rows, lists, and examples retain their required lines. Ordinary prose is not reflowed merely to meet a source-width target.

**Automation:** `dev/audit/help_style.py` and `dev/audit/help_heredoc_reflow.py` report source-shape facets and advisory reflow candidates. Registry coverage is `subset`.

**Semantic remainder:** Review readability, intentional rendered line breaks, appropriate grouping, and whether a row communicates the interface accurately.

**Exceptions:** Literal output and command examples preserve the structure required to demonstrate behavior.

<br />

## Complete help prose (`HELP.PROSE.SENTENCES`)
**Classification:** `advisory` with deterministic recognized-literal portions.

**Scope:** User-facing option and parameter descriptions in Shell, Python, R, Rust, and other maintained command-line interfaces.

Write help prose as complete sentences. Begin each prose paragraph with sentence capitalization and end it with terminal punctuation. Preserve exact option names, paths, identifiers, case-sensitive values, and literal syntax within the sentence; a leading quoted literal does not change the capitalization requirement for the prose that follows it.

This rule governs descriptions rather than usage synopses, metavariables, choice displays, tables, code blocks, machine-readable fragments, or exact external text. Language owners govern source quoting and wrapping while this shared owner governs the rendered prose contract.

**Automation:** `dev/audit/python_source_policy.py` checks capitalization and terminal punctuation for recognized constant `help=` prose in the bounded four-file Python pilot. Shell help structure checks provide description boundaries but do not yet enforce this complete cross-language rule repository-wide. Coverage is `subset`.

**Semantic remainder:** Decide whether text is prose, whether punctuation accurately ends the thought, and whether multiple sentences communicate one coherent option contract.

**Exceptions:** Syntax-only fragments, literal choice displays, generated help, and externally owned exact text require a bounded applicability disposition rather than silent prose normalization.

<br />

## Rendered help (`HELP.RENDERED`)
**Classification:** `advisory` with deterministic interface portions.

**Scope:** Rendered default, detailed, and local help for maintained commands and helpers when applicable.

Rendered help agrees with the accepted public interface, shows required inputs and important modes, and avoids source-only artifacts. Default help is concise but documents every public argument. Detailed help may add mode behavior, references, edge cases, and examples. Local helper help may document globals, shared-state mutation, output, and recoverable status.

A rendered document has one copy of each section. Merge repeated notes or runtime facts. Never expose hidden compatibility aliases, sourced-helper inventories, internal file paths, or maintainer-only call graphs.

**Automation:** `tests/contract/interfaces/test_help_output.sh` renders maintained interfaces and checks selected sections, aliases, and source/render parity with `subset` coverage.

**Semantic remainder:** Review completeness, clarity, audience fit, and whether implementation details are user-actionable.

**Exceptions:** An unavailable optional `--details` surface may be a bounded skip when default help is complete and the interface does not promise details.

<br />

## Help audiences and ownership (`HELP.AUDIENCE`)
**Classification:** `semantic-only`.

**Scope:** Installed command help, callable documentation, public design documentation, local helper help, and maintainer-facing documentation.

Installed help remains independently useful and contains every fact required to run the interface correctly. Callable documentation remains caller-focused and preserves the direct-call contract. Public design documentation owns extended derivation, rationale, compatibility, provenance, and historical context without becoming a runtime dependency. Absence of a public design document from an installed wheel does not permit installed help to omit required behavior.

Write script help for people running the command and local function help for maintainers when a helper parses arguments, reads or creates globals, mutates state, emits structured output, or has non-obvious status. One centralized help function owns full help for each maintained top-level shell command; wrappers call it instead of embedding competing copies. `repository_layout.md` owns placement.

For a command with `--details`, that surface owns the required full document and concise `--help` may omit `Examples` only when it advertises `--details` accurately. Without `--details`, default `--help` owns the full document. `--all_help` may combine concise and full output but does not replace an otherwise valid full owner. A compatibility wrapper may delegate full-help ownership to one canonical wrapper only when argument forwarding and the rendered owner and interface are unambiguous. Do not add `--details` merely to relocate required content. A recognized function-help heredoc is that function's full document.

**Automation:** Documentation-coverage tools inventory recognized owners and candidates but cannot decide usefulness.

**Semantic remainder:** Choose the audience, surface, and appropriate detail.

**Exceptions:** Tiny private pass-through helpers may omit local help when their name, signature, and body make behavior obvious. Compatibility delegation follows the single-owner and unambiguous-forwarding boundary above.

<br />

## Public aliases (`HELP.ALIAS.PUBLIC`)
**Classification:** `advisory` with deterministic documented-fact portions.

**Scope:** Canonical option names, accepted public aliases, hidden compatibility aliases, and their appearance in help and examples.

Parser acceptance alone does not make an alias public. A public alias is intentionally documented and covered as supported. `Usage` advertises canonical long options. `Parameters` lists each public short and long alias exactly once, with the public short before the canonical long spelling. A public dry-run interface provides `-dr` and canonical `--dry_run`; narrower compatibility decisions may retire other historical dry-run spellings but do not remove that public pair.

Hidden compatibility aliases may remain accepted but do not appear in `Usage`, `Parameters`, examples, or ordinary prose. Retired semantic names are not compatibility aliases. Preserve established underscore options, the documented underscore/hyphen boundary, public CSV short aliases beginning with `-c`, and canonical `--dp`. Alias changes migrate parser, help, examples, and interface evidence together.

**Automation:** `dev/audit/help_aliases.py`, parser contracts, and rendered checks compare registered facts and recognized forms with `subset` coverage.

**Semantic remainder:** Decide whether an accepted spelling is public, compatibility-only, or retired.

**Exceptions:** None permit a hidden or retired spelling to appear as public documentation.

<br />

## Exact lexical-token quoting (`HELP.TOKEN.QUOTING`)
**Classification:** `advisory` with deterministic recognized-context portions.

**Scope:** An exact option, identifier, sentinel, value, or complete inline interface expression presented as one lexical object in ordinary explanatory help prose for an applicable active realization. Python and Shell are currently automated; R and Rust remain dormant without enforcement or migration.

Delimit the complete lexical object with straight single quotes. Do not split one complete inline expression into separately quoted fragments, and do not use Markdown backticks as ordinary help-prose delimiters. Structured `Usage` rows, fenced command blocks, executable examples, formulas, and grammar whose structure carries invocation meaning are outside this owner. Ambiguous fragments remain review-owned and are not automatically rewritten.

`dev/audit/help_style.py` retains source recognition. Only its approved exact lexical-object subset moves here; its remaining vocabulary diagnostics retain their existing owners and coverage.

**Automation:** `dev/audit/help_style.py` emits `HELP.TOKEN.QUOTING.SHELL` and `HELP.TOKEN.QUOTING.PYTHON` for safely recognized ordinary-prose lexical objects. `dev/audit/help_contracts.py` validates applicability and checker assignments but does not re-emit source defects. Coverage is `subset`.

**Semantic remainder:** Classify ambiguous fragments and decide whether a complete expression is being discussed as one lexical object.

**Exceptions:** Structured syntax, examples, fences, formulas, and exact external text retain their separately owned literal form.

<br />

## Examples (`HELP.EXAMPLES`)
**Classification:** `advisory` with deterministic structure and ownership portions.

**Scope:** Examples in maintained command, helper, and callable documentation when applicable.

This owner alone decides whether examples apply and whether zero, one, or two are required. Every maintained nontrivial top-level command's full help normally provides at least two materially distinct examples. Callable examples apply when caller obligations are not self-evident. Exactly one example is permitted when exactly one meaningful invocation exists. Omission is permitted only for a genuinely trivial, context-obvious callable. Number examples consecutively, describe their purpose, and use a complete invocation owned by the documented interface.

Use this rendered form:
```text
Examples
--------
  1. Brief description.
    '''bash
    bash script_name.sh \
        --argument_1 "value" \
        --argument_2 "value"
    '''
```

Provide exactly one final `Examples` section. Begin numbering at 1 without skips or duplicates. Each rendered CLI entry has exactly one nonempty `'''bash` block, with no blank line between its description and opening delimiter, and invokes the documented owner directly with public accepted spellings.

A rendered CLI example uses the existing `'''bash` block because it shows invocation through a shell. The invoked implementation may be Bash, Python, R, Rust, or another language. A callable example uses the source-language-native representation governed by the applicable language documentation owner. `HELP.EXAMPLES` alone decides applicability and required count; the language owner governs only representation, source recognition, and language-specific documentation form. Markdown backtick fences are not rendered shell-help delimiters.

For an ordinary multiline CLI invocation, every continuation line begins four spaces beyond the first command line, tabs are forbidden, nonfinal lines render exactly one terminal backslash with no trailing whitespace, and the final line has none. An unquoted help heredoc therefore needs two literal source backslashes for each one intended in rendered output; a quoted heredoc needs one. Control flow, arrays, nested heredocs, structured pipelines, and other syntax with independent indentation are classified as complex snippets rather than forced into this simple-command shape.

A synopsis example uses placeholders rather than concrete values to show an interface's operand shape. It takes the multiline invocation form above; this rule adds only the order of its operands. Put required operands first, then the optional operands most closely related to them under the applicable bracketing rules, and last the generic `[options]` placeholder, which stands for the remaining options detailed elsewhere in the same documentation. A synopsis never places `[options]` before an operand it does not govern. This ordering is a distinct obligation from [`HELP.OPTION.ORDER`](#semantic-option-order-helpoptionorder), which governs rendered `Usage`, option or parameter sections, and reporting; a synopsis is a compressed operand shape, not a total option projection.

Examples are materially distinct when they vary a meaningful mode, input shape, execution path, or output. Preserve safe quoting and explicit continuation structure. Unsafe commands, indirect owners, duplicate signatures, hidden aliases, and undispositioned review candidates are incomplete evidence. Concise default help may omit examples only when its documented detailed owner provides them.

If implementation, callers, tests, and existing documentation do not establish two accurate materially distinct examples, stop for an interface or documentation decision. Do not invent a mode, global, setup requirement, behavior, alias, or filler example.

**Automation:** `dev/audit/help_examples.py` remains the sole owner-specific checker for presence, counts, numbering, command structure, ownership, aliases, duplicate signatures, distinctness, and bounded safety facets. `dev/audit/help_contracts.py` validates only applicability-record structure and references. No checker decides synopsis operand order, because requiredness, governing relationships, and placeholder meaning are not recoverable from example source alone. Coverage is `subset`.

**Semantic remainder:** Review usefulness, safety, material distinctness, representativeness, and synopsis operand order.

**Exceptions:** A trivial zero-argument helper may provide one example when no second meaningful invocation exists and the owner records why.

<br />

## Semantic option order (`HELP.OPTION.ORDER`)
**Classification:** `advisory` with deterministic registered-realization portions.

**Scope:** Explicitly reviewed semantic option roles and their applicable language realizations across parser registration, rendered `Usage`, option or parameter sections, and reporting.

The default category precedence is help/reporting; preview/no-write; execution environment/resources; operation discriminants; inputs; output targets; output lifecycle/cleanup; mode-specific options; formatting; logging/job names; and scheduler or Parallel controls. A modifier stays immediately beside the object or effect it governs. Related semantic groups remain contiguous, and the same logical groups remain together across every applicable surface.

Roles are explicit reviewed data; never infer them from option names. A language realization may use source-appropriate syntax, but it must provide a total projection for every applicable surface. Registered rendered `Usage` projections declare ordered semantic `usage_rows`: each configured group occupies one continuation row, remains intact, and preserves its internal order. A group may contain closely related roles only when its nonempty rationale records the intentional approval. The accepted pilot order is evidence for that interface, not a repository-wide migration.

**Automation:** `dev/audit/help_option_order.py` checks registered Python and Shell realizations from `dev/config/help_contracts.json`. `dev/audit/help_contracts.py` validates record structure and references but emits no order defect. Coverage is `subset`.

**Semantic remainder:** Assign roles, approve deviations, and decide whether two options form one semantic group.

**Exceptions:** An unreviewed option remains a review candidate; its name does not create an implicit role or exemption.

<br />

## Runtime requirements (`HELP.RUNTIME.REQUIREMENTS`)
**Classification:** `advisory` with deterministic grammar and source-evidence portions.

**Scope:** User-facing runtime requirements rendered as a pseudo-heading inside the active top-level `Notes` section.

Use one flat `Runtime requirements:` list. A singleton is unbulleted; two or more peers use bullets sorted by complete displayed text using case-folded comparison and exact text as a tie-breaker. Do not add category subheadings, continuation rows, helper inventories, or repository-local paths.

List environments, interpreters, external executables, and user-provided resources people must install, provide, or understand. Repository-internal scripts and sourced functions are implementation details unless users directly provide or invoke them. Inspect the reachable top-level workflow before changing requirements; grep alone is insufficient.

Use exact callable spelling and versions, including `bash >= 4.4`, `python >= 3.11`, `python3 >= 3.11`, `samtools`, `bamCompare`, and `multiBigwigSummary`. Use `A compatible Conda environment` for the resource; list `conda` or `mamba` only when reachable code invokes it. Attach conditions as `requirement (when <source-grounded trigger>)`, quoting complete option expressions with straight single quotes. Do not combine a tool with a redundant environment alternative; keep genuine interchangeable providers in one entry.

**Automation:** `dev/audit/help_runtime_requirements.py` is the authoritative Runtime parser. It checks recognized grammar, cardinality, order, callables, versions, conditions, environment wording, and overclaims with `subset` coverage. Other audits consume it rather than implement competing grammars.

**Semantic remainder:** Determine reachability, user relevance, conditional sufficiency, and provider equivalence.

**Exceptions:** Maintainer-facing inventories may include internal scripts and functions, but they are not rendered user requirements.

<br />

## Parameter vocabulary (`HELP.PARAMETER.VOCABULARY`)
**Classification:** `advisory` with deterministic controlled-vocabulary portions.

**Scope:** Usage placeholders, parameter types, grouped globals, Boolean documentation, and displayed defaults.

Use `<flag>`, `<str>`, `<bool>`, `<int>`, `<flt>`, `<num>`, `<path>`, `<file>`, `<dir>`, `<csv>`, `<mode>`, `<method>`, `<format>`, `<engine>`, `<layout>`, `<equation>`, `<aligner>`, `<algorithm>`, `<choice>`, `<spec>`, `<time>`, and `<size>`. Choose the narrowest familiar noun; explain CSV element types and structured syntax in `Parameters`.

Use readable parameter types such as `flag`, `int`, `list of file`, `structured string`, or a displayed choice set. Do not use compact `enum:`, `csv:`, or bare `spec` mini-grammar. Presence-only flags consume no value. Explicit Boolean parameters document their domain contract.

Explicit Boolean values accept `true`, `t`, `yes`, `y`, and `1` as true and `false`, `f`, `no`, `n`, and `0` as false, case-insensitively. Surrounding whitespace and empty required values are invalid; successful normalization emits exactly `true` or `false`. A presence-only flag does not consume or normalize an explicit Boolean value.

Write defaults as `(default: ...)`. Numeric defaults are unquoted. String-like defaults, paths, sentinels, enum values, and time strings use straight single quotes.

File-format displays distinguish public choices, accepted aliases, internal canonical values, suffix recognition, and hint applicability. A value may be accepted case-insensitively without appearing in the displayed choice set. Named-path inference may override or ignore a hint only where the interface contract says so. [`HELP.ALIAS.PUBLIC`](#public-aliases-helpaliaspublic) remains authoritative for public and hidden option spellings.

**Automation:** `dev/audit/help_style.py` retains its existing placeholder, type, grouped-entry, Boolean, default, and remaining prose-vocabulary diagnostics. `dev/audit/help_contracts.py` validates the structure and references of registered Python display, canonicalization, suffix, and hint facts without re-emitting a Shell alias defect. Existing Python interface tests prove the registered realization-specific facts. Coverage is `subset`.

**Semantic remainder:** Choose an accurate type, placeholder, and amount of syntax detail.

**Exceptions:** A domain-specific placeholder is permitted when the shared vocabulary would obscure a stable public concept and its description defines it.

<br />

## Canonical parameter descriptions (`PARAMETER.DESCRIPTIONS`)
**Classification:** `advisory` with deterministic registered-core comparisons.

**Scope:** Shared public parameter names and concepts in maintained command, helper, callable, and interface documentation when applicable.

The table supplies a shared semantic core only for an explicitly registered concept family. Identical names or spellings do not establish shared meaning, applicability, or membership. A shared core is an owned proposition: literal common wording is required only when the registered family says so; otherwise each approved natural realization is recorded and preserved. Same-name nonmembers are review candidates, not semantic drift. `mode`, `method`, and `fil_in` do not acquire one repository-wide meaning from their spellings.

The `Type` field names the user-facing logical value form described or accepted by the interface. It is not a Python, Shell, R, or Rust source type, does not prescribe a language-native container or annotation, and is not by itself a serialization contract. Language owners govern source representation, and an independently owned interface or protocol governs serialization when one exists.

| Parameter     | Type              | Canonical description                                                      |
| :---          | :---              | :---                                                                       |
| `csv_fil_in`  | list of file      | Comma-separated list of input file paths.                                  |
| `fil_out`     | file              | Output file path.                                                          |
| `csv_fil_out` | list of file      | Comma-separated list of output file paths.                                 |
| `dir_out`     | dir               | Output directory.                                                          |
| `dir_eo`      | dir               | Directory for stderr and stdout log files.                                 |
| `dir_scr`     | dir               | Maintained entrypoint directory used to resolve adjacent shared functions. |
| `env_nam`     | str               | Conda environment to activate.                                             |
| `threads`     | int               | Number of threads to use.                                                  |
| `ref_fa`      | file              | Reference FASTA file.                                                      |
| `dp`          | int               | Maximum number of decimal places retained for finite emitted values.       |
| `verbose`     | flag              | Run script in verbose mode.                                                |
| `dry_run`     | flag              | Run script in dry-run mode.                                                |
| `nam_job`     | str               | Job name.                                                                  |
| `max_job`     | int               | Maximum number of jobs to run concurrently.                                |
| `slurm`       | flag              | Submit jobs to the Slurm scheduler.                                        |
| `time`        | time              | Slurm job time limit.                                                      |
| `out_ext`     | choice            | Final output extension.                                                    |
| `chr_sizes`   | file              | Chromosome sizes file.                                                     |
| `engine`      | choice            | Processing engine.                                                         |
| `chunk_size`  | int               | Number of records to process per chunk.                                    |
| `siz_bin`     | int               | Bin size in base pairs.                                                    |
| `scl_fct`     | num               | Scaling factor.                                                            |
| `track`       | flag              | Write a companion track file.                                              |
| `strict_bins` | flag              | Require strict bin compatibility.                                          |
| `skip_00`     | choice            | Skip rows where both compared values are zero.                             |
| `skp_pfx`     | list of str       | Comma-separated list of header prefixes to skip.                           |
| `csv_fil_A`   | list of file      | Comma-separated list of file A paths.                                      |
| `csv_fil_B`   | list of file      | Comma-separated list of file B paths.                                      |
| `csv_scl_fct` | list of num       | Comma-separated list of scaling factors.                                   |
| `csv_dep_min` | list of num       | Comma-separated list of minimum-depth values.                              |
| `csv_pseudo`  | list of num       | Comma-separated list of pseudocount values.                                |
| `csv_usr_frg` | list of int       | Comma-separated list of fixed fragment-length values.                      |
| `aligner`     | choice            | Alignment program to use.                                                  |
| `bt2_mode`    | choice            | Bowtie 2 alignment type.                                                   |
| `bwa_alg`     | choice            | BWA algorithm.                                                             |
| `mapq`        | int               | MAPQ threshold.                                                            |
| `req_flg`     | flag              | Require SAM flag bit 2 for properly paired alignments.                     |
| `index`       | path              | Path to the aligner index/reference.                                       |
| `qname`       | flag              | Retain queryname-sorted intermediate alignment files.                      |
| `align_typ`   | choice            | Alignment layout type for input alignment files.                           |
| `aln_typ`     | choice            | Alignment layout type for input alignment files.                           |
| `csv_mip`     | list of file      | Comma-separated list of main IP alignment files.                           |
| `csv_min`     | list of file      | Comma-separated list of main input alignment files.                        |
| `csv_sip`     | list of file      | Comma-separated list of spike-in IP alignment files.                       |
| `csv_sin`     | list of file      | Comma-separated list of spike-in input alignment files.                    |
| `tbl_met`     | file              | siQ-ChIP metadata table.                                                   |
| `cfg_met`     | file              | YAML configuration file for metadata parsing.                              |
| `eqn`         | choice            | siQ-ChIP alpha equation.                                                   |
| `len_def`     | int               | Default fragment length for single-end libraries.                          |
| `len_mip`     | list of number    | Fragment length value(s) for main IP alignment files.                      |
| `len_min`     | list of number    | Fragment length value(s) for main input alignment files.                   |
| `dep_mip`     | list of int       | Sequencing/alignment depth value(s) for main IP alignment files.           |
| `dep_sip`     | list of int       | Sequencing/alignment depth value(s) for spike-in IP alignment files.       |
| `dep_sin`     | list of int       | Sequencing/alignment depth value(s) for spike-in input alignment files.    |
| `fil_A`       | file              | First bedGraph input file, file A.                                         |
| `fil_B`       | file              | Second bedGraph input file, file B.                                        |
| `pseudo`      | structured string | Per-file pseudocount spec 'A[:B]'.                                         |
| `drp_nan`     | flag              | Drop non-finite values from main output.                                   |
| `typ_out`     | choice            | Output file format.                                                        |
| `usr_frg`     | int               | Fixed fragment length.                                                     |
| `chk_chr`     | flag              | Check chromosomes in output alignment files.                               |
| `retain`      | choice            | Species chromosomes to retain.                                             |
| `mito`        | flag              | Retain mitochondrial chromosome.                                           |
| `tg`          | flag              | Retain SP_II_TG chromosome.                                                |
| `mtr`         | flag              | Retain SP_MTR chromosome.                                                  |
| `fil_aln`     | file              | Input BAM or CRAM alignment file.                                          |
| `fq_1`        | file              | First FASTQ input file.                                                    |
| `fq_2`        | file              | Second FASTQ input file.                                                   |
| `sfx_se`      | str               | Suffix to strip from single-end FASTQ filenames.                           |
| `sfx_pe`      | str               | Suffix to strip from paired-end FASTQ read-1 filenames.                    |
| `suffix_se`   | str               | Suffix to strip from single-end FASTQ filenames.                           |
| `suffix_pe`   | str               | Suffix to strip from paired-end FASTQ read-1 filenames.                    |
| `log_out`     | file              | Stdout log file.                                                           |
| `log_err`     | file              | Stderr log file.                                                           |
| `dir_fnc`     | dir               | Base function directory.                                                   |
| `pth_scr_py`  | file              | Python converter script.                                                   |
| `coef`        | num               | Coefficient.                                                               |
| `floor`       | num               | Lower bound.                                                               |
| `eps`         | num               | Zero tolerance epsilon.                                                    |
| `qntl_nz`     | num               | Quantile in percent.                                                       |
| `mode_nz`     | choice            | Epsilon/zero-handling mode.                                                |
| `siz_gen`     | int               | Effective genome size.                                                     |

Use `csv_fil_in`, `fil_out`, and `csv_fil_out` for their registered shared concept families. The spelling `fil_in` remains available to local interfaces, but its signal/alignment, metadata-table, and generic input-file meanings are explicitly non-applicable to one shared canonical core. Use `dp` and `--dp` for decimal precision; do not accept or document `--rnd`, `--round`, `--decimals`, or `--digits` in maintained public CLIs.

**Automation:** `tests/contract/repository/test_parameter_docs_consistency.sh` remains the semantic description-consistency checker. `dev/audit/help_contracts.py` validates concept-family structure, applicability, and schema completeness but does not duplicate its prose-consistency diagnostic. Coverage is `subset`.

**Semantic remainder:** Review whether shared prose is integrated naturally and whether a same-name parameter truly has the shared meaning.

**Exceptions:** A same spelling with a deliberately different meaning requires a narrower name or an explicit owned exclusion.

<br />

### Shared help-contract checker ownership
`dev/config/help_contracts.json` shares applicability and accepted realization facts; it does not make every consumer an owner of every defect. The standards registry maps each diagnostic to one authoritative owner. The contract data assigns each configured diagnostic to one permitted emitter. Seeded tests separately prove actual cross-checker non-overlap.

| Governed concern                            | Authoritative source checker and diagnostics                                                                        | `help_contracts.py` boundary                                                                                                     |
| :---                                        | :---                                                                                                                | :---                                                                                                                             |
| `HELP.ALIAS.PUBLIC` Shell alias sets        | `dev/audit/help_aliases.py`: `HELP.PARAMETER.ALIAS_SET`, `HELP.PARAMETER.ALIAS_DUPLICATE`                           | May cross-reference applicability; does not parse or re-emit Shell alias defects.                                                |
| `HELP.ALIAS.PUBLIC` Python display evidence | Existing Python interface tests where the Shell checker is non-applicable                                           | Validates registered structure and references; any future source diagnostic requires a distinct realization-specific assignment. |
| `PARAMETER.DESCRIPTIONS`                    | `tests/contract/repository/test_parameter_docs_consistency.sh`: semantic prose consistency                          | Validates concept-family structure and applicability only.                                                                       |
| `HELP.EXAMPLES`                             | `dev/audit/help_examples.py`: presence, count, structure, ownership, alias visibility, and distinctness diagnostics | Validates applicability structure and references only.                                                                           |
| `HELP.PARAMETER.VOCABULARY`                 | `dev/audit/help_style.py`: every existing vocabulary diagnostic outside the approved movement                       | Does not re-emit source vocabulary defects.                                                                                      |
| `HELP.TOKEN.QUOTING`                        | `dev/audit/help_style.py`: `HELP.TOKEN.QUOTING.SHELL`, `HELP.TOKEN.QUOTING.PYTHON`                                  | Validates applicability and permitted-emitter assignments only.                                                                  |
| `HELP.OPTION.ORDER`                         | `dev/audit/help_option_order.py`: category, adjacency, contiguity, group, role, and surface-parity diagnostics      | Validates record structure, references, and unique permitted-emitter assignments; does not re-emit order defects.                |
| Shared contract integrity                   | `dev/audit/help_contracts.py`: `HELP.CONTRACT.SCHEMA`, `HELP.CONTRACT.REFERENCE`, `HELP.CONTRACT.APPLICABILITY`     | Own checker output only; never authorizes a source or public-interface change.                                                   |

`dev/audit/help_alias_inventory.json` remains protected interface evidence. Shared applicability never authorizes parallel alias enforcement, and no checker may infer diagnostic ownership merely because it reads the shared data.

<br />

## Shared Python-docstring schema (`HELP.DOCSTRING.SCHEMA`)
**Classification:** `advisory` with deterministic shared-vocabulary portions.

**Scope:** Cross-language section names and prose conventions used by Python docstrings.

Python docstrings reuse shared concepts for parameters, returns, yields, raises, notes, references, examples, attributes, and methods. Delimit option names, paths, sentinels, API names, and other literal prose tokens with straight single quotes, such as `'--write'`; do not use Markdown single- or double-backtick prose delimiters. Backticks remain literal only inside `Examples`, doctest rows, or fenced literal content. Avoid empty or name-repeating sections. Python applicability, public API stability, recursive coverage, and NumPy-specific semantics belong to `PY.DOCSTRING.NUMPY`; exact Python delimiter, summary-line, prefix, closing-line, and post-docstring source form belongs to `PY.DOCSTRING.LAYOUT`.

**Automation:** `dev/audit/help_style.py` checks recognized Python docstrings for single- and double-backtick prose delimiters while excluding `Examples`, doctest rows, and fenced literal content. Other Python-specific checks remain registered to the Python owner.

**Semantic remainder:** Review whether a shared concept maps cleanly to the Python object and audience.

**Exceptions:** Literal backtick syntax inside `Examples`, doctest rows, and fenced literal content remains unchanged. Python syntax, signatures, and types remain literal.

<br />

## Documentation coverage (`DOC.COVERAGE`)
**Classification:** `advisory` with deterministic owner-presence portions.

**Scope:** Centralized help for maintained top-level shell commands and local documentation for reusable shell helpers.

Every maintained top-level shell command has one centralized full-help owner. Reusable public helpers should have local help when behavior, inputs, globals, output, side effects, or status are non-obvious. Lifecycle functions and tiny private pass-through helpers may remain undocumented when their role is mechanically recognizable or obvious. Python docstring applicability belongs to `PY.DOCSTRING.NUMPY`.

**Automation:** `tests/contract/repository/test_doc_coverage.sh` inventories shell help ownership and candidates. Coverage is `advisory`.

**Semantic remainder:** Decide whether a helper is public or nontrivial enough to require local documentation and whether its help is useful.

**Exceptions:** Generated, compatibility-only, and tiny private helper roles may be excluded through an owned disposition.
