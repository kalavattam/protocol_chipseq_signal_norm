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

print_section "${TEST_NAME}"


trap 'cleanup_server_http "${svr_pid:-}"' EXIT


#  Define fixture and output paths for a no-network loopback HTTP download
dir_fx="${ROOT_REPO}/tests/download_fastqs/fixtures"
dir_src="${dir_fx}/source"
mta_tpl="${dir_fx}/metadata/local_se.template.tsv"
src_fq="${dir_src}/se/tiny_download_se.fastq.gz"

tmp="${TEST_DIR_TMP}/execute_download_fastqs_se_local"
dir_out="${tmp}/out"
dir_sym="${tmp}/links"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/download_fastqs"
mta="${tmp}/local_se.tsv"
vw_fq="${tmp}/SRR_LOCAL_SE.fastq"
cnt="${tmp}/SRR_LOCAL_SE.read_count.txt"
fil_out="${dir_out}/SRR_LOCAL_SE.fastq.gz"
sym="${dir_sym}/tiny_download_se.fastq.gz"
log_svr="${dir_log}/execute_download_fastqs_se_local_http.log"
log_run="${dir_log}/execute_download_fastqs_se_local.log"

rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_sym}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty \
    "${src_fq}" \
    "${mta_tpl}" || {
    finish
    exit $?
}

# shellcheck disable=SC2154
if ! \
    require_env_download \
        "${env_nam}" \
        "${dir_log}/execute_download_fastqs_se_local_env.log"
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
url_src="${url_bas}/se/tiny_download_se.fastq.gz"

(
    cd "${dir_src}" || exit 1
    "${py}" -m http.server "${port}" --bind 127.0.0.1
) > "${log_svr}" 2>&1 &
svr_pid=$!

if \
    wait_http_local "${py}" "${url_src}"
then
    record_pass "loopback HTTP fixture server is reachable"
else
    record_fail \
        "loopback HTTP fixture server did not become reachable; see" \
        "$(print_relpath "${log_svr}")"
    finish
    exit $?
fi

sed "s#__BASE_URL__#${url_bas}#g" \
    "${mta_tpl}" > "${mta}"

require_files_nonempty \
    "${mta}" || {
    finish
    exit $?
}


#  Run execute_download_fastqs.sh serially on one loopback HTTP FASTQ URL
if \
    run_capture \
        "execute download-fastqs local HTTP wet run" \
        "${log_run}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/execute_download_fastqs.sh" \
            --threads 1 \
            --infile "${mta}" \
            --dir_out "${dir_out}" \
            --dir_sym "${dir_sym}" \
            --err_out "${dir_err}" \
            --nam_job "test_execute_download_se_local"
then
    record_pass "execute_download_fastqs.sh SE local HTTP wet run exits 0"
else
    record_fail \
        "execute_download_fastqs.sh SE local HTTP wet run failed; see" \
        "$(print_relpath "${log_run}")"
fi

assert_fastq_gzip \
    "${fil_out}" \
    '^@tiny_download_se_read_1$' \
    "${vw_fq}" \
    "${cnt}" \
    "${src_fq}" \
    "execute download-fastqs local FASTQ output"

assert_custom_symlink \
    "${sym}" \
    "execute download-fastqs local custom FASTQ"

assert_file_exists \
    "${dir_err}/test_execute_download_se_local.SRR_LOCAL_SE.stdout.txt" \
    "execute download-fastqs local wget stdout log exists"
assert_file_exists \
    "${dir_err}/test_execute_download_se_local.SRR_LOCAL_SE.stderr.txt" \
    "execute download-fastqs local wget stderr log exists"

finish
