#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_diagnostic_context.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="diagnostic context"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Define fixture and output paths.
dir_fx="${ROOT_REPO}/tests/fixtures/diagnostics"
fil_run="${dir_fx}/script/run.sh"
fil_lib="${dir_fx}/script/lib.sh"

tmp="${TEST_DIR_TMP}/diagnostic_context"


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${tmp}"

require_files_nonempty "${fil_run}" "${fil_lib}" || {
    finish
    exit $?
}


# Compare one call's stderr with the expected lines; stdout must be empty and
# the call must succeed.
function check_case() {
    local nam="${1}"
    local exp="${2}"
    local rc=0
    shift 2

    "$@" > "${tmp}/${nam}.out" 2> "${tmp}/${nam}.err" || rc=$?

    if (( rc == 0 )) && [[ ! -s "${tmp}/${nam}.out" ]]; then
        record_pass "${nam} writes only to stderr and succeeds"
    else
        record_fail "${nam} exited ${rc} or wrote to stdout"
    fi

    if [[ "$(< "${tmp}/${nam}.err")" == "${exp}" ]]; then
        record_pass "${nam} prints '${exp%%$'\n'*}'"
    else
        record_fail \
            "${nam} printed '$(< "${tmp}/${nam}.err")', expected '${exp}'"
    fi
}


# The expected stderr for each case the run script takes as its argument.
declare -A arr_exp=(
    [top]="error(run.sh): script top level"
    [top_func]="error(run.sh): script top level, empty name"
    [fn_scr]="error(run.sh::fn_scr): script function"
    [main]="error(run.sh::main): user main"
    [fn_lib]="error(run.sh>lib.sh::fn_lib): lib function"
    [fn_lib_func]="error(run.sh>lib.sh::fn_lib_func): lib function, named"
    [fn_scr_for]="error(run.sh::fn_scr_for): lib helper for its caller"
    [lib_top]="note(run.sh>lib.sh): lib top level"
)

arr_exp[fn_lib_off]="warning(run.sh>lib.sh::fn_lib_never):"
arr_exp[fn_lib_off]+=" lib function off the stack"$'\n'
arr_exp[fn_lib_off]+="note(run.sh::fn_scr): script function off the stack"

arr_exp[fn_lib_note]="note(run.sh>lib.sh::fn_lib_note): lib note"$'\n'
arr_exp[fn_lib_note]+="note(run.sh>lib.sh::fn_lib_note): lib note, empty name"$'\n'
arr_exp[fn_lib_note]+="warning(run.sh>lib.sh::fn_lib_note): lib warning"

for arg in \
    top top_func fn_scr main fn_lib fn_lib_func fn_scr_for fn_lib_off \
    fn_lib_note lib_top
do
    check_case \
        "${arg}" \
        "${arr_exp[${arg}]}" \
        "${TEST_BASH}" "${fil_run}" "${ROOT_REPO}" "${arg}"
done


# With no script run, '$0' names the shell, and a login shell's leading '-' is
# dropped.
# shellcheck disable=SC2016  # Expand in the child shell, not this one.
src_lib='source "${1}/lib/bash/core/format_outputs.sh"; source "${2}"; fn_lib'
# shellcheck disable=SC2016  # Expand in the child shell, not this one.
src_top='source "${1}/lib/bash/core/format_outputs.sh"; echo_err "shell top"'

check_case \
    "shell_fn_lib" \
    "error(bash>lib.sh::fn_lib): lib function" \
    "${TEST_BASH}" -c "${src_lib}" bash "${ROOT_REPO}" "${fil_lib}"

check_case \
    "login_fn_lib" \
    "error(bash>lib.sh::fn_lib): lib function" \
    "${TEST_BASH}" -c "${src_lib}" -bash "${ROOT_REPO}" "${fil_lib}"

check_case \
    "shell_top" \
    "error(bash): shell top" \
    "${TEST_BASH}" -c "${src_top}" bash "${ROOT_REPO}"


# A real driver names the run script and the library that refuses.
rc=0
"${TEST_BASH}" "${ROOT_REPO}/bin/submit_download_fastqs.sh" --srr SRR_X \
    > "${tmp}/driver.out" 2> "${tmp}/driver.err" || rc=$?

if (( rc != 0 )); then
    record_pass "submit_download_fastqs.sh without '--url_1' fails"
else
    record_fail "submit_download_fastqs.sh without '--url_1' succeeded"
fi

msg_drv="error(submit_download_fastqs.sh>check_inputs.sh::validate_var):"
msg_drv+=" 'url_1' is empty or unset."

assert_file_exact_line \
    "${tmp}/driver.err" \
    "${msg_drv}" \
    "a library refusal names the run script and the library"


finish
