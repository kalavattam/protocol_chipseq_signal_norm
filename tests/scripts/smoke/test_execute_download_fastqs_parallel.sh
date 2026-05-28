#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_download_fastqs_parallel.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="execute download-fastqs GNU Parallel"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

print_section "${TEST_NAME}"

if ! \
    is_parallel_enabled
then
    record_skip \
        "GNU Parallel download-fastqs check disabled;" \
        "set RUN_PARALLEL=1 to enable"
    finish
    exit $?
fi


trap 'cleanup_server_http "${svr_pid:-}"' EXIT


#  Define fixture and output paths for a no-network GNU Parallel HTTP download
dir_fx="${ROOT_REPO}/tests/download_fastqs/fixtures"
dir_src="${dir_fx}/source"
mta_tpl="${dir_fx}/metadata/local_mixed.template.tsv"
src_se="${dir_src}/tiny_download_se.fastq.gz"
src_r1="${dir_src}/tiny_download_pe_R1.fastq.gz"
src_r2="${dir_src}/tiny_download_pe_R2.fastq.gz"

tmp="${TEST_DIR_TMP}/execute_download_fastqs_parallel"
dir_out="${tmp}/out"
dir_sym="${tmp}/links"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/download_fastqs"
mta="${tmp}/local_mixed.tsv"
vw_se="${tmp}/SRR_LOCAL_SE.fastq"
vw_r1="${tmp}/SRR_LOCAL_PE_R1.fastq"
vw_r2="${tmp}/SRR_LOCAL_PE_R2.fastq"
cnt_se="${tmp}/SRR_LOCAL_SE.read_count.txt"
cnt_r1="${tmp}/SRR_LOCAL_PE_R1.read_count.txt"
cnt_r2="${tmp}/SRR_LOCAL_PE_R2.read_count.txt"
out_se="${dir_out}/SRR_LOCAL_SE.fastq.gz"
out_r1="${dir_out}/SRR_LOCAL_PE_R1.fastq.gz"
out_r2="${dir_out}/SRR_LOCAL_PE_R2.fastq.gz"
sym_se="${dir_sym}/tiny_download_se.fastq.gz"
sym_r1="${dir_sym}/tiny_download_pe_R1.fastq.gz"
sym_r2="${dir_sym}/tiny_download_pe_R2.fastq.gz"
cfg="${dir_err}/test_execute_download_parallel.config_parallel.txt"
log_svr="${dir_log}/execute_download_fastqs_parallel_http.log"
log_run="${dir_log}/execute_download_fastqs_parallel.log"

rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_sym}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty \
    "${src_se}" \
    "${src_r1}" \
    "${src_r2}" \
    "${mta_tpl}" || {
        finish
        exit $?
    }

# shellcheck disable=SC2154
if ! \
    require_env_download \
        "${env_nam}" \
        "${dir_log}/execute_download_fastqs_parallel_download_env.log"
then
    finish
    exit $?
fi

if ! \
    require_env_parallel \
        "${env_nam}" \
        "${dir_log}/execute_download_fastqs_parallel_parallel_env.log"
then
    finish
    exit $?
fi

if ! \
    py="$(find_python)"
then
    record_fail "local HTTP fixture server requires python or python3 on PATH"
    finish
    exit $?
fi

port="$(find_port_free "${py}")"
url_bas="http://127.0.0.1:${port}"
url_src="${url_bas}/tiny_download_se.fastq.gz"

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


#  Run execute_download_fastqs.sh through local GNU Parallel on one mixed
#+ SE/PE metadata TSV; with two rows and '--threads 2', the execute wrapper
#+ writes a GNU Parallel config and dispatches one submit command per row
if \
    run_capture \
        "execute download-fastqs GNU Parallel local HTTP wet run" \
        "${log_run}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/execute_download_fastqs.sh" \
            --threads 2 \
            --infile "${mta}" \
            --dir_out "${dir_out}" \
            --dir_sym "${dir_sym}" \
            --err_out "${dir_err}" \
            --nam_job "test_execute_download_parallel"
then
    record_pass \
        "execute_download_fastqs.sh GNU Parallel local HTTP wet run exits 0"
else
    record_fail \
        "execute_download_fastqs.sh GNU Parallel local HTTP wet run failed;" \
        "see $(print_relpath "${log_run}")"
fi

assert_file_nonempty \
    "${cfg}" \
    "execute download-fastqs GNU Parallel config"

if [[ -s "${cfg}" ]]; then
    assert_pattern_found \
        "${cfg}" \
        "${TEST_BASH} ${ROOT_REPO}/scripts/submit_download_fastqs.sh" \
        "execute download-fastqs GNU Parallel config uses Bash-prefixed" \
        "submit command"
fi

assert_fastq_gzip \
    "${out_se}" \
    '^@tiny_download_se_read_1$' \
    "${vw_se}" \
    "${cnt_se}" \
    "${src_se}" \
    "execute download-fastqs GNU Parallel SE accession FASTQ output"
assert_fastq_gzip \
    "${out_r1}" \
    '^@tiny_download_pe_pair_1/1$' \
    "${vw_r1}" \
    "${cnt_r1}" \
    "${src_r1}" \
    "execute download-fastqs GNU Parallel PE R1 accession FASTQ output"
assert_fastq_gzip \
    "${out_r2}" \
    '^@tiny_download_pe_pair_1/2$' \
    "${vw_r2}" \
    "${cnt_r2}" \
    "${src_r2}" \
    "execute download-fastqs GNU Parallel PE R2 accession FASTQ output"

assert_custom_symlink \
    "${sym_se}" \
    "execute download-fastqs GNU Parallel SE custom FASTQ"
assert_custom_symlink \
    "${sym_r1}" \
    "execute download-fastqs GNU Parallel PE R1 custom FASTQ"
assert_custom_symlink \
    "${sym_r2}" \
    "execute download-fastqs GNU Parallel PE R2 custom FASTQ"

assert_file_exists \
    "${dir_err}/test_execute_download_parallel.SRR_LOCAL_SE.stdout.txt" \
    "execute download-fastqs GNU Parallel SE wget stdout log exists"
assert_file_exists \
    "${dir_err}/test_execute_download_parallel.SRR_LOCAL_SE.stderr.txt" \
    "execute download-fastqs GNU Parallel SE wget stderr log exists"
assert_file_exists \
    "${dir_err}/test_execute_download_parallel.SRR_LOCAL_PE_R1.stdout.txt" \
    "execute download-fastqs GNU Parallel PE R1 wget stdout log exists"
assert_file_exists \
    "${dir_err}/test_execute_download_parallel.SRR_LOCAL_PE_R1.stderr.txt" \
    "execute download-fastqs GNU Parallel PE R1 wget stderr log exists"
assert_file_exists \
    "${dir_err}/test_execute_download_parallel.SRR_LOCAL_PE_R2.stdout.txt" \
    "execute download-fastqs GNU Parallel PE R2 wget stdout log exists"
assert_file_exists \
    "${dir_err}/test_execute_download_parallel.SRR_LOCAL_PE_R2.stderr.txt" \
    "execute download-fastqs GNU Parallel PE R2 wget stderr log exists"

finish
