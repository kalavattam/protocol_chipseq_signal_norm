#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_download_fastqs_se_local.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="execute download-fastqs SE local"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

rec_section "${TEST_NAME}"


trap 'cleanup_http_server "${server_pid:-}"' EXIT


#  Define fixture and output paths for a no-network loopback HTTP download
dir_fx="${ROOT_REPO}/tests/download_fastqs/fixtures"
dir_src="${dir_fx}/source"
metadata_template="${dir_fx}/metadata/local_se.template.tsv"
source_fastq="${dir_src}/tiny_download_se.fastq.gz"

tmp="${TEST_DIR_TMP}/execute_download_fastqs_se_local"
dir_out="${tmp}/out"
dir_sym="${tmp}/links"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/download_fastqs"
metadata="${tmp}/local_se.tsv"
view_fastq="${tmp}/SRR_LOCAL_SE.fastq"
count_reads="${tmp}/SRR_LOCAL_SE.read_count.txt"
outfile="${dir_out}/SRR_LOCAL_SE.fastq.gz"
symlink="${dir_sym}/tiny_download_se.fastq.gz"
log_server="${dir_log}/execute_download_fastqs_se_local_http.log"
log_run="${dir_log}/execute_download_fastqs_se_local.log"

rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_sym}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_exist "${source_fastq}" "${metadata_template}" || {
    finish
    exit $?
}

if ! \
    require_download_env \
        "${env_nam}" \
        "${dir_log}/execute_download_fastqs_se_local_env.log"
then
    finish
    exit $?
fi

if ! py="$(find_python)"; then
    rec_fail "local HTTP fixture server requires python or python3 on PATH"
    finish
    exit $?
fi

port="$(find_free_port "${py}")"
url_base="http://127.0.0.1:${port}"
url_source="${url_base}/tiny_download_se.fastq.gz"

(
    cd "${dir_src}" || exit 1
    "${py}" -m http.server "${port}" --bind 127.0.0.1
) > "${log_server}" 2>&1 &
server_pid=$!

if wait_for_local_http "${py}" "${url_source}"; then
    rec_pass "loopback HTTP fixture server is reachable"
else
    rec_fail \
        "loopback HTTP fixture server did not become reachable; see" \
        "$(rec_relpath "${log_server}")"
    finish
    exit $?
fi

sed "s#__BASE_URL__#${url_base}#g" \
    "${metadata_template}" > "${metadata}"

require_files_exist "${metadata}" || {
    finish
    exit $?
}


#  Run execute_download_fastqs.sh serially on one loopback HTTP FASTQ URL
if \
    run_capture \
        "execute download-fastqs local HTTP wet run" "${log_run}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/execute_download_fastqs.sh" \
            --threads 1 \
            --infile "${metadata}" \
            --dir_out "${dir_out}" \
            --dir_sym "${dir_sym}" \
            --err_out "${dir_err}" \
            --nam_job "test_execute_download_se_local"
then
    rec_pass "execute_download_fastqs.sh SE local HTTP wet run exits 0"
else
    rec_fail \
        "execute_download_fastqs.sh SE local HTTP wet run failed; see" \
        "$(rec_relpath "${log_run}")"
fi

assert_downloaded_fastq \
    "${source_fastq}" "${outfile}" \
    "execute download-fastqs local FASTQ output" \
    '^@tiny_download_se_read_1$' "${view_fastq}" "${count_reads}"

assert_custom_symlink \
    "${symlink}" "execute download-fastqs local custom FASTQ"

assert_file_exists \
    "${dir_err}/test_execute_download_se_local.SRR_LOCAL_SE.stdout.txt" \
    "execute download-fastqs local wget stdout log exists"
assert_file_exists \
    "${dir_err}/test_execute_download_se_local.SRR_LOCAL_SE.stderr.txt" \
    "execute download-fastqs local wget stderr log exists"

finish
