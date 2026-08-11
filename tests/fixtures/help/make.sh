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


# Resolve paths relative to 'tests/fixtures'
dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
dir_fix="${dir_scr}"

# Source shared fixture-generation helpers.
# shellcheck source=tests/support/fixture_helpers.sh
source "${dir_scr}/../../support/fixture_helpers.sh"

# Declare every generated path up front.
dir_exp="${dir_fix}/expected"
fil_examples="${dir_exp}/combine_parts_scaling_factor.examples.txt"


# Remove stale outputs so regeneration is idempotent.
rm_file "${dir_fix}" "${fil_examples}"

mkdirs "${dir_exp}"


# Author the expected Examples section literally.
#
# "Author text and derive everything else" stops here, and deliberately. This
# file is the expected result of extracting the Examples section from
# 'bin/combine_parts_scaling_factor.sh --help'. Deriving it by running that
# same extraction would leave the contract comparing the extractor's output
# against the extractor's output, which passes for any output at all. A golden
# file is only evidence while a person wrote it, so it is typed here and only
# changes when someone decides the rendered help should change.
#
# The delimiter is quoted, so the backslash continuations and the triple single
# quotes below survive verbatim.
cat << 'EOM' > "${fil_examples}"
Examples
--------
  1. Combine spike-in scaling-factor parts.
    '''bash
    bash combine_parts_scaling_factor.sh \
        --mode spike \
        --csv_fil_in results.spike.tsv.part.000001,results.spike.tsv.part.000000 \
        --fil_out results.tsv
    '''

  2. Preview siQ-ChIP part combination.
    '''bash
    bash combine_parts_scaling_factor.sh \
        --dry_run \
        --mode siq \
        --csv_fil_in results.siq.tsv.part.000001,results.siq.tsv.part.000000 \
        --fil_out results.tsv
    '''
EOM

succeed "generated rendered-help fixtures under ${dir_fix}"
