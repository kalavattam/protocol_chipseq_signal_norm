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
    echo "error(shell):" \
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

# Declare every generated path up front. The directory names the role the file
# plays in the assertion: these are bedGraph inputs, which is what a workflow
# input is, rather than a checker verdict.
dir_bdg="${dir_fix}/bedgraph"
fil_bdg_A="${dir_bdg}/pair_A.bdg"
fil_bdg_B="${dir_bdg}/pair_B.bdg"

# Remove stale outputs so regeneration is idempotent.
rm_files "${dir_fix}" \
    "${fil_bdg_A}" \
    "${fil_bdg_B}"

mkdirs "${dir_bdg}"


# Author both tracks literally. The values are chosen so every quantity the
# edgeR estimator derives from them is a round number a reader can check by
# hand, which is what makes a failure legible rather than merely red:
#
#     L_A = 1 + 2 + 3 =  6        L_B = 4 + 6 + 8 = 18
#     L_bar = (6 + 18) / 2 = 12
#     prior_scaled_A = 2 * 6 / 12 = 1.0
#     prior_scaled_B = 2 * 18 / 12 = 3.0
#
# The pair is deliberately imbalanced 1:3 rather than near-equal, so a test
# that asserts the depth correction cannot pass by returning the nominal
# 'prior.count' for both tracks. Three rows on a uniform 10 bp grid also let
# the bin width be inferred rather than supplied.

cat << 'EOM' > "${fil_bdg_A}"
chrI	0	10	1
chrI	10	20	2
chrI	20	30	3
EOM

cat << 'EOM' > "${fil_bdg_B}"
chrI	0	10	4
chrI	10	20	6
chrI	20	30	8
EOM


succeed "generated compute-pseudo fixtures under ${dir_fix}"
