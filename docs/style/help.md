
# Shell Help Style Guide
This file is the detailed style guide for Bash script help, shell-function help, `show_help`, `detail_*`, `help_*`, and heredoc-based CLI documentation.

`AGENTS.md` is the top-level repository guide. This file owns detailed shell help-text rules; if the two diverge, update both in the same patch when practical.

<br />

## Recommended heading order
Use this top-level order everywhere, substituting exactly one appropriate argument-section heading at the argument-section position:
```txt
Usage:
Description:
Arguments: / Positional arguments: / Keyword arguments:
Dependencies:
  Recommended environment:   # if applicable
  External programs:
  Shell scripts:             # if applicable
  Python scripts:            # if applicable
  Configuration files:       # if applicable
  Sourced function scripts:  # only in maintainer-facing docs
Returns:                     # if applicable
Notes:
Examples:
#TODO:                       # if user-facing
```

<br />

## Formatting rules
Use 2 spaces per indentation level in help text.

Markdown prose in this file and other style docs is not subject to the source-code 80-character preference. Do not hard-wrap Markdown prose solely for line length; use natural paragraphs and let the renderer or editor wrap text. Preserve intentional line breaks in code fences, help-output examples, tables, lists, and other structured Markdown.

Shell help text is also not subject to the source-code 80-character wrapping preference. Wrap help text for readability, structure, and user comprehension, not to satisfy code line-length limits. Preserve clear argument tables, examples, and prose when a longer line is easier to read.

Separate major help sections with one blank line. That is, place one empty line before each new top-level section heading such as `Description:`, an argument-section heading, `Dependencies:`, `Returns:`, `Notes:`, `Examples:`, and `#TODO:`.

In argument sections, use this format:
```txt
  -x, --example  <type>
    Description text.
```

If an argument consumes no value, use `<flag>` and still keep two spaces before `<flag>`, e.g.:
```txt
  -h, --help  <flag>
    Display this help message and exit.
```

For positional-argument sections, put the argument number, name, and type on one line. Put the description on the next line, indented one additional level, and separate entries with one blank line:
```txt
  01  arg_name  <type>
    Description text.

  02  next_arg  <type>
    Description text.
```

If a function has 10 or more positional arguments, zero-pad argument numbers 1-9. For example, use `01` through `09`, then `10`, `11`, and so on.

Nontrivial shell-function help docs should include a `Returns:` section.

Within nested sections such as `Dependencies:`, use this format:
```txt
Dependencies:
  Recommended environment:
    - env_protocol

  External programs:
    - program

  Shell scripts:
    - script_name.sh
```

If a `Returns:` section has multiple distinct items, format it as a hyphen list:
```txt
Returns:
  - Prints optional arguments as a single comma-delimited line to stdout.
  - Returns 0 when optional arguments are built successfully; 1 otherwise.
```

If a `Returns:` section is prose rather than a list, keep it as prose and wrap continuation lines with two-space indentation:
```txt
Returns:
  Prints optional arguments as a single comma-delimited line to stdout. Returns
  0 when optional arguments are built successfully; 1 otherwise.
```

Avoid hard tab characters in help text.

Argument-section headings may be one of:
  - `Arguments:`
  - `Positional arguments:`
  - `Keyword arguments:`

Use `Keyword arguments:` when documenting flags/options such as `--help`, `--mode`, or `--fil_out`.

Use `Positional arguments:` when documenting positional `$1`, `$2`, etc. arguments.

Use `Arguments:` only when a help block combines positional and keyword arguments in one section, or when the script/function has a simple mixed interface where splitting the section would be unnecessary.

<br />

## Help audiences
Write user-facing script help for users running the command. It should explain accepted options, required inputs, generated outputs, user-actionable dependencies, important modes, and practical examples.

Write local shell-function help for maintainers reading the implementation. It may describe expected globals, generated globals, stdout/stderr behavior, mutation of shared state, and return status.

Keep exhaustive sourced-helper inventories out of user-facing help. If those inventories are useful, put them in maintainer-facing documentation or deterministic audit tooling.

<br />

### Short help versus detailed help
#### Short help
For default `--help`, use this structure:
```txt
Usage:
Description:
Arguments: / Positional arguments: / Keyword arguments:
Returns:    # if applicable
Notes:
```

The short help should include every accepted argument, but only one concise sentence per argument (although some exceptions apply). For long scripts, it should avoid deep method lists, paper references, implementation details, and long examples. It should tell the user what the script does, what is required, and where to get more.

<br />

#### Detailed help
For `--details`, use this structure:
```txt
Usage:
Description:
Arguments: / Positional arguments: / Keyword arguments:
Dependencies:
Returns:    # if applicable
Notes:
Examples:
#TODO:
```

Detailed help can include synonym lists, mode-specific behavior, edge cases, references, and examples. There is no cap on description length.

<br />

### Local shell-function help
Use local function help for nontrivial helpers whose behavior is not obvious from the function name and body. This is especially useful for helpers that parse their own arguments, read several globals, write generated globals, print structured output, or return nonzero for recoverable validation failures.

For helpers that read important shared state, add `Expected globals:`. For helpers that write important shared state, add `Generated globals:` or describe the mutation in `Notes:`.

For helpers that print to stdout or stderr, say so in `Returns:` or `Notes:`. For helpers that write a global array, such as `cmd_bld`, prefer `Generated globals:` plus a short note about the array shape.

<br />

### Examples format
Use this format for examples:
```txt
Examples:
  1. Brief description.
    '''bash
    bash script_name.sh
      --argument_1 "value"
      --argument_2 "value"
    '''
```

<br />

### Rules:
- Use `'''bash` and `'''`, not Markdown tick fences.
- Do not include a blank line between the numbered example description and `'''bash`.
- Do not use shell line-continuation backslashes in examples.
- Examples do not need to be directly copy-paste runnable. They should be easy to read and sufficient for the user to understand the command structure.
- Use 2-space indentation steps inside examples.

<br />

### Argument ordering rules
Use the same argument order in `Usage:`, short argument sections, and detailed argument sections.

In `Usage:` blocks, group continuation lines by workflow role rather than by
source-code line length. User-facing help is exempt from the 80-character
source-code preference, so keep a semantic group on one line when that is
clearer.

In `Usage:` blocks, show required options without square brackets and
optional, defaulted, or conditional options with square brackets:
```txt
--dir_out <dir>       # Required
[--threads <int>]    # Optional or has a default
[--ref_fa <file>]    # Conditional; explain condition in the argument body
```

Use parentheses for required mutually exclusive groups, and keep optional
members bracketed inside the group:
```txt
(--csv_infile <csv:file> [--ref_fa <file>] | --csv_fil_A <csv:file> --csv_fil_B <csv:file>)
```

Order arguments by workflow role:
1. Help and reporting flags.
2. Execution controls.
3. Analysis mode controls.
4. Input files.
5. Output location and naming.
6. Signal-mode options.
7. Ratio-mode options.
8. Output formatting options.
9. Logging and job naming.
10. Scheduler and parallelization options.

Within each group, keep closely related arguments together, and keep flags near the values they modify when this improves readability.

For example:
```txt
Usage:
  example.sh
    [--help] [--verbose] [--dry_run]
    [--threads <int>]
    [--mode <enum:a|b>] [--method <enum:x|y>]
    (--csv_infile <csv:file> [--ref_fa <file>] | --csv_fil_A <csv:file> --csv_fil_B <csv:file>)
    --dir_out <dir> [--prefix <str>]
    [--csv_scl_fct <csv:spec>]
    [--dp <int>]
    --dir_eo <dir> [--nam_job <str>]
    [--max_job <int>] [--slurm] [--time <time>]
```

<br />

### Canonical option names and synonyms
In `Usage:`, list one clear preferred spelling for each argument. Prefer the established spelling already used in similar scripts and help blocks.

In argument sections, document accepted user-facing aliases together on the same argument line. Do not move accepted aliases into a separate compatibility sentence.

Use this style:
```txt
  -short, -other_short, --long_name, --other_long_name  <type>
```
For example:
```txt
  -i, -fi, --infile, --infiles, --fil_in, --csv_infile, --csv_infiles  <csv:file>
```

Do not use separate compatibility-alias prose such as:

> Compatibility aliases include '--infile', '--infiles', '--fil_in', and '--csv_infiles'.

If separate compatibility-alias prose already exists, treat it as something to fix rather than a pattern to copy.

For alias ordering, prefer the order already established in similar scripts and help blocks. If the best order is unclear, preserve the current order rather than inventing a new global canonical order.

When both underscore and hyphen forms are accepted by the parser, document only the underscore form. For example, document --dry_run, not both --dry_run and --dry-run.

Prefer --dp as the canonical option for decimal precision. Continue accepting and documenting legacy aliases on the same argument line:
```txt
  -dp, --dp, --rnd, --round, --decimals, --digits  <int>
    Maximum number of decimal places retained for finite emitted values (default: ${dp_or_rnd}).
```

Use the script's current internal default variable in this text. Use `${dp}`
in scripts that have migrated to `dp`, and use `${rnd}` in scripts that still
use `rnd`.

<br />

### Canonical variable names
Prefer the canonical long option name as the internal variable name when practical.

For decimal precision, use `dp` as the canonical variable and `--dp` as the preferred long option. Continue accepting legacy aliases such as `--rnd`, `--round`, `--decimals`, and `--digits`, but assign them to `dp`.

<br />

### Default-value formatting
When documenting default values, use `(default: ...)`.

Use no quotes for numeric defaults:
```txt
(default: 4)
(default: 24)
```

Use single quotes for string-like defaults, including enums, paths, filenames, sentinels, empty strings, and time strings:
```txt
(default: 'signal')
(default: 'bedGraph.gz')
(default: '${dir_out}/logs')
(default: '0:30:00')
(default: '')
(default: 'NA')
```

Avoid double quotes in prose defaults unless the double quotes are part of the value being documented.

<br />

## Recommended `<...>` argument identifiers
Use a small fixed vocabulary and avoid inventing a new identifier for each argument. The purpose is quick recognition, not perfect type theory.

Recommended set:
```txt
<flag>      Boolean flag; no value consumed.
<str>       General string.
<bool>      Boolean-like string. Accepted true values: true, t, yes, y, 1. Accepted false values: false, f, no, n, 0. Values are case-insensitive and should be standardized internally to true or false.
<int>       Integer.
<flt>       Floating-point number.
<num>       Integer or floating-point number.
<path>      File or directory path. (Subset of <str>.)
<file>      File path. (Subset of <str>.)
<dir>       Directory path. (Subset of <str>.)
<csv:str>   Comma-separated strings.
<csv:int>   Comma-separated integers.
<csv:flt>   Comma-separated floats.
<csv:num>   Comma-separated integers/floats.
<csv:file>  Comma-separated file paths.
<csv:dir>   Comma-separated directory paths.
<csv:spec>  Comma-separated structured mini-syntax entries, with the element syntax explained in the argument body.
<enum:...>  One value from a closed set.
<spec>      Structured mini-syntax, explained in the argument body.
<time>      Time string, e.g., h:mm:ss. (Subset of <str>.)
<time>      Time string, e.g., h:mm:ss. (Subset of <str>.)
<size>      Size threshold string, e.g., 1, 1k, 2M, 0.5G, depending on script support.
```

Avoid previously used <mlt> because it is too vague. Previously, it seemed to mean “multiple possible values,” “comma-separated list,” “sentinel-enabled value,” or “enum,” depending on context. But those should be distinguished.

Avoid adding secondary element annotations such as `<element:str>` after the main argument identifier. If an argument takes a comma-separated list of structured entries, use `<csv:spec>` and explain the per-element syntax in the argument body.

<br />

## Dependency documentation rules
Before editing dependency lists, inspect the relevant top-level executable directly. Do not infer dependency lists from grep output alone. For each top-level executable, inspect directly invoked commands, project scripts passed downstream, and code called by the script when needed to understand the user-facing workflow.

Use the ordering principle “recommended environment first; dependency details second.” When `env_protocol` or another repository environment is the recommended installation path, list that environment first, then list external programs as transparency and troubleshooting information.

User-facing help should list dependencies users may need to install, provide, or understand. Do not maintain exhaustive sourced-helper inventories in user-facing CLI help. Those inventories are maintainer-facing and should live in developer docs or deterministic audit tooling if needed.

<br />

### Recommended environment:
- Include this subsection when a script is normally run from a project Conda/Mamba environment.
- List the recommended environment by name, e.g., `env_protocol`.
- Do not use the recommended environment as a reason to hide tool-level dependencies. Users still need external-program lists for debugging missing commands, module-based systems, manual environment recreation, and incomplete environment checks.

<br />

### External programs:
- List every external program directly invoked by the top-level executable workflow.
- Include programs used only in conditional branches, but mark the condition.
- Include optional execution tools such as GNU Parallel and Slurm when they affect a user-visible execution path.
- Include tools required for format-specific paths such as CRAM handling, alignment, trimming, compression, or scheduler submission.
- Do not list shell keywords, shell builtins, helper functions, or programs used only internally by sourced helper functions.
  - Do not list shell keywords or builtins such as `case`, `if`, `for`, `while`, `echo`, `printf`, `source`, `declare`, `typeset`, `local`, `shift`, `return`, `exit`, or `unset` as external dependencies.
  - `cat` is an external program, not a Bash builtin, but avoid listing ubiquitous shell-plumbing commands such as `cat` unless the top-level executable meaningfully depends on their external behavior.
- List external program names as the command checked or invoked, using lowercase command names except when the conventional project name is clearer for users.
  - Use `Bash >= 4.4`, not `bash`.
  - Use `GNU Parallel`, not `parallel`.
  - Use lowercase command names such as `awk`, `python`, `samtools`, `sbatch`, `rm`, `gzip`, `sort`, and `realpath`.
- List conditional dependencies with the condition, e.g., `GNU Parallel, when '--threads > 1'` or `rm, when '--threads > 1'`.

<br />

### Shell scripts:
- List project shell scripts that are directly executed, submitted, or passed to an executor by the top-level executable workflow.
- Include downstream entrypoint scripts when the wrapper dispatches them in a user-visible way.
- Do not list sourced function scripts under `Shell scripts:`.

<br />

### Sourced function scripts:
- Use this section only in maintainer-facing documentation.
- Do not include sourced helper scripts or helper-function inventories in normal user-facing CLI help.
- If maintainer-facing dependency documentation lists sourced helpers, list only directly sourced scripts and directly called functions. Do not list transitive helper calls made internally by sourced helper functions.
