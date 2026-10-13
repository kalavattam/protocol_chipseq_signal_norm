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

# Declare every generated path up front. 'script/' holds a run script and the
# library it sources, so a diagnostic can come from either file.
dir_spt="${dir_fix}/script"

fil_run="${dir_spt}/run.sh"
fil_lib="${dir_spt}/lib.sh"


# Remove stale outputs so regeneration is idempotent.
rm_files "${dir_fix}" "${fil_run}" "${fil_lib}"

mkdirs "${dir_spt}"


# Author every script literally. The delimiter is quoted, so no '$' in a
# fixture body reaches the shell. This recipe only writes: which prefix each
# call must print is asserted in
# 'tests/contract/interfaces/test_diagnostic_context.sh'.

# The run script: argument 1 is the repository root, and argument 2 names the
# case, either a top-level call or a function to call.
cat << 'EOM' > "${fil_run}"
#!/usr/bin/env bash

source "${1}/lib/bash/core/format_outputs.sh"

case_diag="${2}"
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"

function fn_scr() {
    echo_err "script function"
}


function fn_scr_for() {
    fn_lib_for "${FUNCNAME[0]}"
}


function main() {
    echo_err "user main"
}


if [[ "${case_diag}" == "top" ]]; then
    echo_err "script top level"
elif [[ "${case_diag}" == "top_func" ]]; then
    echo_err_func "" "script top level, empty name"
elif [[ "${case_diag}" != "lib_top" ]]; then
    "${case_diag}"
fi
EOM

# The library: one function per helper form, and a top-level note.
cat << 'EOM' > "${fil_lib}"
function fn_lib() {
    echo_err "lib function"
}


function fn_lib_func() {
    echo_err_func "${FUNCNAME[0]}" "lib function, named"
}


function fn_lib_for() {
    echo_err_func "${1}" "lib helper for its caller"
}


function fn_lib_off() {
    echo_warn_func "fn_lib_never" "lib function off the stack"
    echo_note_func "fn_scr" "script function off the stack"
}


function fn_lib_never() {
    :
}


function fn_lib_note() {
    echo_note "lib note"
    echo_note_func "" "lib note, empty name"
    echo_warn "lib warning"
}


if [[ "${case_diag:-}" == "lib_top" ]]; then
    echo_note "lib top level"
fi
EOM


succeed "generated diagnostic-context fixtures under ${dir_fix}"
