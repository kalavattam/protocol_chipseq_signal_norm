
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
  External programs:
  Shell scripts:             # if applicable
  Python scripts:            # if applicable
  Configuration files:       # if applicable
  Sourced function scripts:
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

Separate major help sections with two blank lines. That is, place two empty lines before each new top-level section heading such as `Description:`, an argument-section heading, `Dependencies:`, `Returns:`, `Notes:`, `Examples:`, and `#TODO:`.

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

Within nested sections such as `Dependencies:`, use this format:
```txt
Dependencies:
  External programs:
    - program
  Sourced function scripts:
    - script_name.sh
      + function_name
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
```

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
<flag>       Boolean flag; no value consumed.
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
Before editing dependency lists, inspect the relevant top-level executable directly. Do not infer dependency lists from grep output alone. For each top-level executable, inspect directly invoked commands, project scripts passed downstream, directly called helper functions, helper scripts listed in `source_helpers`, and code sourced and/or called by the script when needed to understand the top-level workflow.

<br />

### External programs:
- List every external program directly invoked by the top-level executable workflow.
- Include programs used only in conditional branches, but mark the condition.
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
- Do not list sourced function scripts under `Shell scripts:`; list those under `Sourced function scripts:`.

<br />

### Sourced function scripts:
- List each sourced function script required by the top-level executable workflow.
- Under each script, list only the helper functions from that script that are directly called by the top-level executable.
- Do not list transitive helper calls made internally by sourced helper functions.
- If the top-level executable sources helper scripts through `source_helpers.sh`, list `source_helpers.sh` and `source_helpers`.
- If a sourced function script is loaded but none of its functions are directly called by the top-level executable, treat that as a cleanup finding. Remove that script from the sourcing list unless there is a documented reason to keep it.

<br />

## Functions within shell function scripts
The inventory below is reference material for looking up which helper functions are provided by each script; it does not prescribe emitted help-text indentation.
```txt
align_fastqs.sh
  - _check_path_whitespace
  - align_fastqs

calculate_scaling_factor.sh
  - _get_len_idx
  - _get_dep_idx
  - _detect_typ_bam
  - _resolve_typ_fil
  - _get_expr_filter
  - _count_align_bam
  - _calculate_frag_avg
  - _compute_scl_fct
  - _import_shell_asgmt
  - _parse_metadata
  - _calculate_dep_fct
  - _calculate_dep_arr
  - _compute_dep_all
  - _generate_fmt_str
  - _get_fil_out_part
  - process_samp_siq
  - process_samp_spike

check_args.sh
  - require_optarg
  - check_arg_supplied
  - check_args_mut_excl
  - check_flags_mut_excl
  - check_match
  - check_str_delim

check_env.sh
  - check_env_installed
  - check_pgrm_path

check_inputs.sh
  - validate_var
  - validate_file
  - validate_dir
  - validate_var_file
  - validate_var_dir
  - debug_var
  - check_arr_files
  - check_arr_lengths
  - check_file_dir_exists
  - debug_arr_contents
  - check_arr_nonempty
  - check_arr_len_bcst

check_numbers.sh
  - check_flt_nonneg
  - check_flt_pos
  - check_format_time
  - check_int_nonneg
  - check_int_pos
  - check_arr_int_pos
  - check_arr_num_pos
  - check_scl_fct

check_source.sh
  - err_source_only

check_unity.sh
  - check_unity

construct_find.sh
  - construct_find

filter_bam.sh
  - _parse_args_filter_bam
  - _validate_args_filter_bam
  - _check_chr_bam
  - _finalize_bam_filter
  - _filter_sam_chr
  - _cleanup_filter_bam_tmp
  - filter_bam_sc
  - filter_bam_sp

format_outputs.sh
  - echo_err
  - echo_err_func
  - echo_warn
  - echo_warn_func
  - format_print_cmd
  - print_banner_pretty
  - print_cmd_array
  - print_cmd_pretty
  - summarize_sig_norm

handle_env.sh
  - activate_env
  - _current_errexit_nounset
  - _restore_errexit_nounset
  - _handle_env_deactivate
  - _handle_env_activate_success
  - _handle_env_activate
  - handle_env

handle_exit_interactive.sh
  - exit_0
  - exit_1

manage_parallel.sh
  - determine_cores
  - print_parallel_info
  - reset_max_job
  - set_params_parallel

manage_slurm.sh
  - set_logs_slurm

populate_array_empty.sh
  - populate_array_empty

process_region.sh
  - check_region
  - check_region_bam
  - check_region_bdg

process_sequences.sh
  - check_seq_type
  - check_string_fastqs
  - get_paired_suffix
  - parse_fastq_entry
  - pair_fastqs
  - pair_fqs

process_tables.sh
  - check_table
  - check_table_column
  - check_table_scaling_factor
  - extract_field_str
  - _validate_arg_csl
  - _validate_args_table
  - _parse_table_core
  - parse_table
  - parse_table_simple

run_python.sh
  - _resolve_dir_rep_run_py
  - to_module
  - run_py

source_helpers.sh
  - _source_helper_err
  - _source_helper_resolve
  - source_once
  - source_helpers

wrap_cmd.sh
  - get_submit_logs
  - print_built_cmd
```
