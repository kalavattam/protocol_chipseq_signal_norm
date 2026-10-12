#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_download_fastqs_verify_local.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="execute download-fastqs verification local"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


if ! \
    is_download_enabled
then
    record_skip \
        "download-fastqs verification local check disabled; set" \
        "RUN_DOWNLOAD=1 to enable"
    finish
    exit $?
fi


# Define fixture and output paths for loopback-only HTTP downloads.
dir_fx="${ROOT_REPO}/tests/fixtures/download_fastqs"
dir_src="${dir_fx}/source"
dir_mtd="${dir_fx}/metadata"

tmp="${TEST_DIR_TMP}/execute_download_fastqs_verify_local"
dir_log="${TEST_DIR_LOG}/download_fastqs"
log_svr="${dir_log}/execute_download_fastqs_verify_local_http.log"
log_env="${dir_log}/execute_download_fastqs_verify_local_env.log"


print_section "${TEST_NAME}"

trap 'cleanup_server_http "${svr_pid:-}"' EXIT

rm -rf "${tmp}"
mkdir -p "${tmp}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty \
    "${dir_src}/se/tiny_download_se.fastq.gz" \
    "${dir_src}/se/tiny_download_se_truncated.fastq.gz" \
    "${dir_mtd}/local_se_verified.template.tsv" \
    "${dir_mtd}/local_se_missing.template.tsv" || {
    finish
    exit $?
}

# shellcheck disable=SC2154
if ! \
    require_env_download \
        "${env_nam}" \
        "${log_env}"
then
    finish
    exit $?
fi

if ! \
    py="$(find_python_loopback)"
then
    record_fail \
        "local HTTP fixture server requires a Python that can bind loopback" \
        "sockets"
    finish
    exit $?
fi

port="$(find_port_free "${py}")"
url_bas="http://127.0.0.1:${port}"

(
    cd "${dir_src}" || exit 1
    "${py}" -m http.server "${port}" --bind 127.0.0.1
) > "${log_svr}" 2>&1 &
svr_pid=$!

if \
    wait_http_local "${py}" "${url_bas}/se/tiny_download_se.fastq.gz"
then
    record_pass "loopback HTTP fixture server is reachable"
else
    record_fail \
        "loopback HTTP fixture server did not become reachable; see" \
        "$(print_relpath "${log_svr}")"
    finish
    exit $?
fi


# Run one template; print the run log path.
function run_template() {
    local nam="${1}"
    local tpl="${2}"
    local dir="${tmp}/${nam}"
    shift 2

    mkdir -p "${dir}/out" "${dir}/links" "${dir}/logs"
    sed "s#__BASE_URL__#${url_bas}#g" \
        "${dir_mtd}/${tpl}.template.tsv" > "${dir}/input.tsv"

    run_capture \
        "${nam}" \
        "${dir_log}/execute_download_fastqs_verify_${nam}.log" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/execute_download_fastqs.sh" \
            --env_nam "${env_nam}" \
            --threads 1 \
            --fil_in "${dir}/input.tsv" \
            --dir_out "${dir}/out" \
            --dir_sym "${dir}/links" \
            --dir_eo "${dir}/logs" \
            --nam_job "verify_${nam}" \
            "$@"
}


# Each row: case name, template, expected exit (0 or 1), the files expected in
# the output directory, and a pattern the run log must contain ('-' for none).
se="SRR_LOCAL_SE.fastq.gz"
pe="SRR_LOCAL_PE_R1.fastq.gz SRR_LOCAL_PE_R2.fastq.gz"
vrf="failed verification 3 times:"
cnt="column 'fastq_md5' has 1 value(s) but there are 2 URL(s)"
rows=(
    "se_verified|local_se_verified|0|${se}|-"
    "pe_verified|local_pe_verified|0|${pe}|-"
    "pe_md5_only|local_pe_md5_only|0|${pe}|-"
    "se_no_md5|local_se|0|${se}|1 of 1 FASTQ files have no MD5"
    "se_wrong_md5|local_se_wrong_md5|1||${vrf} MD5 is"
    "se_wrong_size|local_se_wrong_size|1||${vrf} size is"
    "se_truncated|local_se_truncated|1||${vrf} file fails 'gzip -t'"
    "se_missing|local_se_missing|1||failed to download"
    "pe_md5_count|local_pe_md5_count|1||${cnt}"
    "pe_bad_md5_r2|local_pe_wrong_md5_r2|1|${pe%% *}|${vrf} MD5 is"
    "pe_conflicting_md5|local_pe_conflicting_md5|1||conflicting FASTQ URL, MD5"
)

for row in "${rows[@]}"; do
    IFS='|' read -r nam tpl rc_exp fil_exp pat <<< "${row}"
    log="${dir_log}/execute_download_fastqs_verify_${nam}.log"

    rc=0
    run_template "${nam}" "${tpl}" || rc=$?

    if (( (rc == 0) == (rc_exp == 0) )); then
        record_pass "${nam} exits as expected (${rc})"
    else
        record_fail \
            "${nam} exited ${rc}, expected ${rc_exp}; see" \
            "$(print_relpath "${log}")"
    fi

    fil_out="$(
        find "${tmp}/${nam}/out" \
            -mindepth 1 \
            -maxdepth 1 \
            -type f \
            -exec basename {} \; \
            | LC_ALL=C sort \
            | tr '\n' ' '
    )"
    fil_out="${fil_out% }"

    if [[ "${fil_out}" == "${fil_exp}" ]]; then
        record_pass "${nam} leaves exactly the expected files"
    else
        record_fail \
            "${nam} left '${fil_out}' in the output directory, expected" \
            "'${fil_exp}'"
    fi

    if [[ "${pat}" != "-" ]]; then
        assert_pattern_found \
            "${log}" \
            "${pat}" \
            "${nam} reports: ${pat}"
    fi
done


# Verified downloads keep the source bytes.
if \
    cmp -s \
        "${tmp}/se_verified/out/SRR_LOCAL_SE.fastq.gz" \
        "${dir_src}/se/tiny_download_se.fastq.gz"
then
    record_pass "verified SE download matches its source"
else
    record_fail "verified SE download differs from its source"
fi

# A file with an MD5 is not reported as unverified.
assert_pattern_absent \
    "${dir_log}/execute_download_fastqs_verify_se_verified.log" \
    'have no MD5' \
    "se_verified gives no unverified-download warning"

# 'wget' retries refused connections and is limited to 5 tries.
assert_pattern_found \
    "${dir_log}/execute_download_fastqs_verify_se_no_md5.log" \
    '^wget .*--tries=5 --retry-connrefused --waitretry=10 ' \
    "wget retries refused connections, up to 5 tries"


finish
