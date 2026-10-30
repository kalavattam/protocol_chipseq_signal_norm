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
dir_scp="${dir_fix}/script"
fil_cnt="${dir_scp}/counted_helper.sh"


# Remove stale outputs so regeneration is idempotent.
rm_files "${dir_fix}" "${fil_cnt}"

mkdirs "${dir_scp}"


# Write a helper that counts how many times it is sourced, so a test can show
# that 'source_once' loads it once.
# shellcheck disable=SC2016  # The line is written literally, for later use.
printf '%s\n' 'COUNTED_HELPER_LOADS=$(( ${COUNTED_HELPER_LOADS:-0} + 1 ))' \
    > "${fil_cnt}"


succeed "generated submit-bootstrap fixtures under ${dir_fix}"
