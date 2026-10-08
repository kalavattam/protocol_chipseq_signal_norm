#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_download_fastqs_opt_applicability.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="download-fastqs option applicability"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Define fixture and output paths. Per 'HELP.PARAMETER.APPLICABILITY', '--time'
# is ignored with a warning without '--slurm'. Dry runs download nothing, so
# the metadata URLs point at an unused loopback port.
dir_fx="${ROOT_REPO}/tests/fixtures/download_fastqs"
mta_tpl="${dir_fx}/metadata/local_se.template.tsv"

tmp="${TEST_DIR_TMP}/download_fastqs_opt_applicability"
mta="${tmp}/local_se.tsv"


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${tmp}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty "${mta_tpl}" || {
    finish
    exit $?
}

sed "s#__BASE_URL__#http://127.0.0.1:9#g" "${mta_tpl}" > "${mta}"


# Dry-run 'execute' locally, writing to its own output directory.
function run_execute() {
    local dir_out="${1}"
    shift 1

    mkdir -p "${dir_out}/out" "${dir_out}/links" "${dir_out}/logs"
    "${TEST_BASH}" "${ROOT_REPO}/bin/execute_download_fastqs.sh" \
        --dry_run \
        --env_nam "${env_nam}" \
        --threads 1 \
        --fil_in "${mta}" \
        --dir_out "${dir_out}/out" \
        --dir_sym "${dir_out}/links" \
        --dir_eo "${dir_out}/logs" \
        "$@" 2>&1
}


# A time limit is for Slurm jobs only, so a local run warns about it.
rc=0
out="$(run_execute "${tmp}/exe_time" --time 1:00:00)" || rc=$?

if [[
    "${rc}" -eq 0
    && "${out}" == *"warning("*"'--time' has no effect without '--slurm'"*
]]; then
    record_pass "execute warns about '--time' without '--slurm'"
else
    record_fail "execute did not warn about '--time' without '--slurm'"
fi

# Without '--time', the default limit draws no warning.
rc=0
out="$(run_execute "${tmp}/exe_default")" || rc=$?

if [[ "${rc}" -eq 0 && "${out}" != *"has no effect"* ]]; then
    record_pass "execute does not warn about the default time limit"
else
    record_fail "execute warned about the default time limit (exit ${rc})"
fi

finish
