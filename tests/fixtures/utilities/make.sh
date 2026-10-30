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

# Declare every generated path up front.
dir_txt="${dir_fix}/text"
dir_log="${dir_fix}/log"
dir_fq_fnd="${dir_fix}/fastq/find"

fil_smp="${dir_txt}/sample.txt"
fil_lrg="${dir_log}/large.log"
fil_emp="${dir_log}/empty.log"
fq_fnd_r1="${dir_fq_fnd}/a_R1.fastq.gz"
fq_fnd_r2="${dir_fq_fnd}/a_R2.fastq.gz"


# Require gzip for the compressed FASTQ pair.
require_cmd gzip "to generate utilities fixtures."

# Remove stale outputs so regeneration is idempotent.
rm_files \
    "${dir_fix}" \
    "${fil_smp}" \
    "${fil_lrg}" \
    "${fil_emp}" \
    "${fq_fnd_r1}" \
    "${fq_fnd_r2}"

mkdirs "${dir_txt}" "${dir_log}" "${dir_fq_fnd}"


# Write a one-line file for the 'symlink_files.sh' dry run.
printf 'smoke\n' > "${fil_smp}"

# Write a 2048-byte file and an empty file for the 'compress_remove_files.sh'
# dry runs, which would compress the first and delete the second if they acted.
printf '%02048d' 0 > "${fil_lrg}"
: > "${fil_emp}"

# Write a paired-end name pair for 'find_files.sh', one record each, since the
# script matches names only.
printf '@a\nACGT\n+\nIIII\n' | gzip -n -c > "${fq_fnd_r1}"
printf '@a\nACGT\n+\nIIII\n' | gzip -n -c > "${fq_fnd_r2}"


succeed "generated utilities fixtures under ${dir_fix}"
