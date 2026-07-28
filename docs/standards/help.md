
# Help and documentation standard
This document owns shared help and documentation semantics for maintained shell interfaces and the bounded cross-language schema used by Python docstrings. Shell source form, Python applicability, Markdown source form, and path placement remain with their domain owners.

<br />

## Shared section schema (`HELP.SECTION.SCHEMA`)
**Classification:** `advisory` with deterministic vocabulary and ordering portions.

**Scope:** User-facing shell help, local shell-function help, and Python docstrings where the same documentation concepts apply.

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
    [--help] [--threads <int>]
    --fil_in <file> --fil_out <file> [--ref_fa <file>]
```

Use square brackets as the primary optionality signal. Unconditional required arguments are unbracketed; optional, defaulted, conditional, and mode-specific arguments are bracketed and their descriptions state defaults or conditions directly. Parenthesize required mutually exclusive branches and bracket the whole group when the choice itself is optional.

Use the same argument order in `Usage` and `Parameters`. The default workflow order is help and reporting, execution controls, analysis mode, inputs, output location and naming, mode-specific options, output formatting, logging and job naming, then scheduler and Parallel controls. Place dry-run immediately after help and other reporting controls such as verbose, and before environment, thread, workflow-input, and output options. Adapt group names to the workflow without reordering the same interface differently across help surfaces.

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

**Scope:** Rendered default, detailed, and local help for maintained shell commands and helpers.

Rendered help agrees with the accepted public interface, shows required inputs and important modes, and avoids source-only artifacts. Default help is concise but documents every public argument. Detailed help may add mode behavior, references, edge cases, and examples. Local helper help may document globals, shared-state mutation, output, and recoverable status.

A rendered document has one copy of each section. Merge repeated notes or runtime facts. Never expose hidden compatibility aliases, sourced-helper inventories, internal file paths, or maintainer-only call graphs.

**Automation:** `tests/contract/interfaces/test_help_output.sh` renders maintained interfaces and checks selected sections, aliases, and source/render parity with `subset` coverage.

**Semantic remainder:** Review completeness, clarity, audience fit, and whether implementation details are user-actionable.

**Exceptions:** An unavailable optional `--details` surface may be a bounded skip when default help is complete and the interface does not promise details.

<br />

## Help audiences and ownership (`HELP.AUDIENCE`)
**Classification:** `semantic-only`.

**Scope:** Selection of short, detailed, local, and maintainer-facing documentation surfaces.

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

## Examples (`HELP.EXAMPLES`)
**Classification:** `advisory` with deterministic structure and ownership portions.

**Scope:** Examples in full command help and recognized nontrivial function help.

Every maintained top-level command's full help includes `Examples`. A recognized nontrivial function-help document provides at least two examples. Number examples consecutively, describe their purpose, and use a complete invocation owned by the documented interface.

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

Provide exactly one final `Examples` section. Begin numbering at 1 without skips or duplicates. Each numbered entry has exactly one nonempty `'''bash` block, with no blank line between its description and opening delimiter, and invokes the documented owner directly with public accepted spellings. Markdown backtick fences are not rendered shell-help delimiters.

For an ordinary multiline CLI invocation, every continuation line begins four spaces beyond the first command line, tabs are forbidden, nonfinal lines render exactly one terminal backslash with no trailing whitespace, and the final line has none. An unquoted help heredoc therefore needs two literal source backslashes for each one intended in rendered output; a quoted heredoc needs one. Control flow, arrays, nested heredocs, structured pipelines, and other syntax with independent indentation are classified as complex snippets rather than forced into this simple-command shape.

Examples are materially distinct when they vary a meaningful mode, input shape, execution path, or output. Preserve safe quoting and explicit continuation structure. Unsafe commands, indirect owners, duplicate signatures, hidden aliases, and undispositioned review candidates are incomplete evidence. Concise default help may omit examples only when its documented detailed owner provides them.

If implementation, callers, tests, and existing documentation do not establish two accurate materially distinct examples, stop for an interface or documentation decision. Do not invent a mode, global, setup requirement, behavior, alias, or filler example.

**Automation:** `dev/audit/help_examples.py` and `dev/audit/help_heredoc_reflow.py` check presence, counts, numbering, command structure, ownership, aliases, duplicate signatures, and bounded safety facets. Coverage is `subset`.

**Semantic remainder:** Review usefulness, safety, material distinctness, and representativeness.

**Exceptions:** A trivial zero-argument helper may provide one example when no second meaningful invocation exists and the owner records why.

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

**Automation:** Help-style and parameter-consistency audits check recognized placeholders, types, grouped entries, Boolean vocabulary, and defaults with `subset` coverage.

**Semantic remainder:** Choose an accurate type, placeholder, and amount of syntax detail.

**Exceptions:** A domain-specific placeholder is permitted when the shared vocabulary would obscure a stable public concept and its description defines it.

<br />

## Canonical parameter descriptions (`PARAMETER.DESCRIPTIONS`)
**Classification:** `advisory` with deterministic registered-core comparisons.

**Scope:** Shared public parameter names in shell help, Python CLI help, and maintained interface documentation.

The table supplies the shared semantic core. Integrate workflow constraints after that core without repeating it. Same-name drift is enforced only for registered shared families; generic local names are not repository-global contracts.

| Parameter     | Type              | Canonical description                                                      |
| :---          | :---              | :---                                                                       |
| `fil_in`      | file              | Input file path.                                                           |
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
| `mode`        | choice            | Workflow mode.                                                             |
| `method`      | choice            | Workflow method.                                                           |
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

Use `fil_in`, `csv_fil_in`, `fil_out`, and `csv_fil_out` for canonical input and output file arguments. Use `dp` and `--dp` for decimal precision; do not accept or document `--rnd`, `--round`, `--decimals`, or `--digits` in maintained public CLIs.

**Automation:** `tests/contract/repository/test_parameter_docs_consistency.sh` extracts this table only from the `PARAMETER.DESCRIPTIONS` owner section and compares maintained interface roots. Coverage is `subset`.

**Semantic remainder:** Review whether shared prose is integrated naturally and whether a same-name parameter truly has the shared meaning.

**Exceptions:** A same spelling with a deliberately different meaning requires a narrower name or an explicit owned exclusion.

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
