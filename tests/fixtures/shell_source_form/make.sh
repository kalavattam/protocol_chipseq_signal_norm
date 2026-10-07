#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: make.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
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
# checker must return for the script inside it; 'format/' holds the input and
# expected pair for the comment-wrap formatter, which claims a rewrite rather
# than a verdict.
dir_acc="${dir_fix}/accepted"
dir_bnd="${dir_fix}/boundary"
dir_rej="${dir_fix}/rejected"
dir_fmt="${dir_fix}/format"

fil_acc_wrap="${dir_acc}/comment_wrap.sh"
fil_bnd_wrap="${dir_bnd}/comment_wrap.sh"
fil_rej_wrap="${dir_rej}/comment_wrap.sh"
fil_fmt_input="${dir_fmt}/comment_wrap_input.sh"
fil_fmt_expected="${dir_fmt}/comment_wrap_expected.sh"


# Remove stale outputs so regeneration is idempotent.
rm_files "${dir_fix}" \
    "${fil_acc_wrap}" \
    "${fil_bnd_wrap}" \
    "${fil_rej_wrap}" \
    "${fil_fmt_input}" \
    "${fil_fmt_expected}"

mkdirs "${dir_acc}" "${dir_bnd}" "${dir_rej}" "${dir_fmt}"


# Author every script literally. The delimiter is quoted, so no '$' in a
# fixture body reaches the shell. This recipe only writes: which lines the
# checker must report is asserted in
# 'tests/unit/dev_audit/test_shell_source_form.py'.

# Accepted: greedily wrapped comment prose and every structural boundary the
# wrap facet recognizes, none of which it may report.
cat << 'EOM' > "${fil_acc_wrap}"
#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: comment_wrap.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail


# Read every record of a template before deciding whether to keep the whole
# template, because a decision made one record at a time splits pairs apart.
# Another sentence in the same paragraph fills the remaining width as well.
#
# A blank comment line separates paragraphs, so this short line ends one.
#
# Each step below is short:
# - Sort the records by name so that every template sits on adjacent lines of
#   the file.
# - Count the primary records of each template.
# TODO: decide whether the filter belongs in its own module.
value=1

# _parse_args
# _filter_records
# _write_output

# shellcheck disable=SC2034  # The value is read by the sourcing script.
value=2

function print_value() {
    # Print the value that the caller passed, followed by one newline, so that
    # callers can read it back with a plain command substitution.
    printf '%s\n' "${1}"
}

cat << 'EOS'
# A heredoc body is literal text,
# and its line breaks are not prose.
EOS
EOM

# Boundary: breaks on the exact edge of the fit test, none of which it may
# report: a word ending at column 80, a bracketed unit that fits only in part,
# a hyphenated join one column too long, the line that ends a list, an indented
# display, and a fenced block.
cat << 'EOM' > "${fil_bnd_wrap}"
#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: comment_wrap.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail


# A word that would end at column 80, one past the boundary, stays on the
# second line, so this break is correct.
value=1

# Keep a pair only when both mates pass, as the help states for each
# ('pe' and 'se') layout, and drop the pair otherwise.
value=2

# Pass the reference only to entries that read or write CRAM, as a commas-
# joined list of paths that the caller builds for each entry in turn.
value=3

# Each step below is short:
# - Count the primary records of each template.
# A sentence after a list ends that list, and its break is not tested.
value=4

# Run the filter on one file:
#     filter_records sample.bam sample.kept.bam 30
# The indented line above is a display that no list item introduces.
value=5

# '''
# A fenced block keeps its line breaks,
# because it holds literal text.
# '''
value=6
EOM

# Rejected: five blocks, each with one premature break: plain prose whose next
# word ends at column 79, a list-item continuation, a hyphenated join that fits
# without a separator, a bracketed unit that fits whole, and a comment indented
# inside a function.
cat << 'EOM' > "${fil_rej_wrap}"
#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: comment_wrap.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail


# Read every record of a template before deciding whether to keep the
# template, because a decision made one record at a time splits pairs apart.
value=1

# Each step below is short:
# - Sort the records by name so that every template sits on adjacent lines
#   of the file, which the filter needs before it can count primary records.
value=2

# Pass the reference only to entries that read or write CRAM, as a comma-
# joined list of paths that the caller builds for each entry in turn.
value=3

# Keep a pair only when both mates pass, in both layouts
# ('pe' and 'se'), and drop the pair otherwise.
value=4

function print_value() {
    # Print the value that the caller passed, followed by one newline,
    # so that callers can read it back with a plain command substitution.
    printf '%s\n' "${1}"
}
EOM

# Format: an input with premature breaks in prose, a list item, and before a
# bracketed unit, and the expected output, written by hand. The hyphenated
# paragraph and the heredoc body stay as written.
cat << 'EOM' > "${fil_fmt_input}"
#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: comment_wrap.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail


# Read every record of a template before deciding whether to keep the
# template, because a decision made one record at a time splits pairs
# apart. A second sentence in the same paragraph is refilled with the
# first.
value=1

# Each step below is short:
# - Sort the records by name so that every template sits on adjacent
#   lines of the file, which the filter needs before it can count primary
#   records.
value=2

# Choose GNU Parallel when more than one job runs and otherwise run serially
# ('par_job == 1'), keeping the bracketed expression whole.
value=3

# A line that ends in a hyphen may break one word or share a suspended
# hyphen, as in single-
# and paired-end data, so this paragraph is left for review.
value=4

cat << 'EOS'
# A heredoc body is literal text, and the
# formatter leaves its line breaks alone.
EOS
EOM

cat << 'EOM' > "${fil_fmt_expected}"
#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: comment_wrap.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail


# Read every record of a template before deciding whether to keep the template,
# because a decision made one record at a time splits pairs apart. A second
# sentence in the same paragraph is refilled with the first.
value=1

# Each step below is short:
# - Sort the records by name so that every template sits on adjacent lines of
#   the file, which the filter needs before it can count primary records.
value=2

# Choose GNU Parallel when more than one job runs and otherwise run serially
# ('par_job == 1'), keeping the bracketed expression whole.
value=3

# A line that ends in a hyphen may break one word or share a suspended
# hyphen, as in single-
# and paired-end data, so this paragraph is left for review.
value=4

cat << 'EOS'
# A heredoc body is literal text, and the
# formatter leaves its line breaks alone.
EOS
EOM


succeed "generated shell source-form fixtures under ${dir_fix}"
