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

# Declare every generated path up front. These are workflow inputs rather than
# checker inputs, so the directories name the kind of data they hold, matching
# 'compute_signal' and the other workflow fixture sets.
dir_ref="${dir_fix}/reference"
dir_fq="${dir_fix}/fastq"
dir_sam="${dir_fix}/sam"
fil_ref="${dir_ref}/tiny.fa"
fil_fastq="${dir_fq}/tiny_se.fastq"
fil_sam="${dir_sam}/tiny_signal.sam"


# Remove stale outputs so regeneration is idempotent.
rm_files "${dir_fix}" "${fil_ref}" "${fil_fastq}" "${fil_sam}"

mkdirs "${dir_ref}" "${dir_fq}" "${dir_sam}"


# Write the 108-base single-sequence reference the remote runner indexes.
cat << 'EOM' > "${fil_ref}"
>I
GATCGTACCTAGGCTAACGTTGACCGTTAACGATCGTAGCTAGGATCCGTTACGATCGATGCTAGCTTACCGGATCAAGCTTAGGCTAATCGGCTAAGGTTCCGATTA
EOM

# Write the single 30-base read the alignment job consumes.
cat << 'EOM' > "${fil_fastq}"
@tiny_se_read_1
ACGTTGACCGTTAACGATCGTAGCTAGGAT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOM

# Write two 10-base single-end alignments for the signal job. The rows go
# through 'write_sam_line' rather than a heredoc, so no tab is typed as an
# invisible control character that an editor could convert to spaces.
{
    write_sam_line '@HD' 'VN:1.6' 'SO:coordinate'
    write_sam_line '@SQ' 'SN:I' 'LN:108'

    write_sam_line \
        'se_fwd' '0' 'I' '1' '60' '10M' \
        '*' '0' '0' 'ACGTACGTAC' 'FFFFFFFFFF'

    write_sam_line \
        'se_rev' '16' 'I' '21' '60' '10M' \
        '*' '0' '0' 'ACGTACGTAC' 'FFFFFFFFFF'
} > "${fil_sam}"


succeed "generated Slurm wet-validation fixtures under ${dir_fix}"
