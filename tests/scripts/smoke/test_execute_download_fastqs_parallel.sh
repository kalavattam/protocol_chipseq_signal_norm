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

rec_section "${TEST_NAME}"

if ! is_parallel_enabled; then
    rec_skip \
        "GNU Parallel download-fastqs check disabled;" \
        "set RUN_PARALLEL=1 to enable"
    finish
    exit $?
fi


trap 'cleanup_http_server "${server_pid:-}"' EXIT


#  Define fixture and output paths for a no-network GNU Parallel HTTP download
dir_fx="${ROOT_REPO}/tests/download_fastqs/fixtures"
dir_src="${dir_fx}/source"
metadata_template="${dir_fx}/metadata/local_mixed.template.tsv"
source_se="${dir_src}/tiny_download_se.fastq.gz"
source_r1="${dir_src}/tiny_download_pe_R1.fastq.gz"
source_r2="${dir_src}/tiny_download_pe_R2.fastq.gz"

tmp="${TEST_DIR_TMP}/execute_download_fastqs_parallel"
dir_out="${tmp}/out"
dir_sym="${tmp}/links"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/download_fastqs"
metadata="${tmp}/local_mixed.tsv"
view_fastq_se="${tmp}/SRR_LOCAL_SE.fastq"
view_fastq_r1="${tmp}/SRR_LOCAL_PE_R1.fastq"
view_fastq_r2="${tmp}/SRR_LOCAL_PE_R2.fastq"
count_reads_se="${tmp}/SRR_LOCAL_SE.read_count.txt"
count_reads_r1="${tmp}/SRR_LOCAL_PE_R1.read_count.txt"
count_reads_r2="${tmp}/SRR_LOCAL_PE_R2.read_count.txt"
outfile_se="${dir_out}/SRR_LOCAL_SE.fastq.gz"
outfile_r1="${dir_out}/SRR_LOCAL_PE_R1.fastq.gz"
outfile_r2="${dir_out}/SRR_LOCAL_PE_R2.fastq.gz"
symlink_se="${dir_sym}/tiny_download_se.fastq.gz"
symlink_r1="${dir_sym}/tiny_download_pe_R1.fastq.gz"
symlink_r2="${dir_sym}/tiny_download_pe_R2.fastq.gz"
config="${dir_err}/test_execute_download_parallel.config_parallel.txt"
log_server="${dir_log}/execute_download_fastqs_parallel_http.log"
log_run="${dir_log}/execute_download_fastqs_parallel.log"

rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_sym}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_exist \
    "${source_se}" "${source_r1}" "${source_r2}" "${metadata_template}" || {
        finish
        exit $?
    }

if ! \
    require_download_env \
        "${env_nam}" \
        "${dir_log}/execute_download_fastqs_parallel_download_env.log"
then
    finish
    exit $?
fi

if ! \
    require_parallel_env \
        "${env_nam}" \
        "${dir_log}/execute_download_fastqs_parallel_parallel_env.log"
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


#  Run execute_download_fastqs.sh through local GNU Parallel on one mixed
#+ SE/PE metadata TSV. With two rows and '--threads 2', the execute wrapper
#+ writes a GNU Parallel config and dispatches one submit command per row.
if \
    run_capture \
        "execute download-fastqs GNU Parallel local HTTP wet run" "${log_run}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/execute_download_fastqs.sh" \
            --threads 2 \
            --infile "${metadata}" \
            --dir_out "${dir_out}" \
            --dir_sym "${dir_sym}" \
            --err_out "${dir_err}" \
            --nam_job "test_execute_download_parallel"
then
    rec_pass "execute_download_fastqs.sh GNU Parallel local HTTP wet run exits 0"
else
    rec_fail \
        "execute_download_fastqs.sh GNU Parallel local HTTP wet run failed; see" \
        "$(rec_relpath "${log_run}")"
fi

assert_file_nonempty "${config}" "execute download-fastqs GNU Parallel config"

if [[ -s "${config}" ]]; then
    assert_grep_pattern \
        "${config}" \
        "${TEST_BASH} ${ROOT_REPO}/scripts/submit_download_fastqs.sh" \
        "execute download-fastqs GNU Parallel config uses Bash-prefixed submit command"
fi

assert_downloaded_fastq \
    "${source_se}" "${outfile_se}" \
    "execute download-fastqs GNU Parallel SE accession FASTQ output" \
    '^@tiny_download_se_read_1$' "${view_fastq_se}" "${count_reads_se}"
assert_downloaded_fastq \
    "${source_r1}" "${outfile_r1}" \
    "execute download-fastqs GNU Parallel PE R1 accession FASTQ output" \
    '^@tiny_download_pe_pair_1/1$' "${view_fastq_r1}" "${count_reads_r1}"
assert_downloaded_fastq \
    "${source_r2}" "${outfile_r2}" \
    "execute download-fastqs GNU Parallel PE R2 accession FASTQ output" \
    '^@tiny_download_pe_pair_1/2$' "${view_fastq_r2}" "${count_reads_r2}"

assert_custom_symlink \
    "${symlink_se}" "execute download-fastqs GNU Parallel SE custom FASTQ"
assert_custom_symlink \
    "${symlink_r1}" "execute download-fastqs GNU Parallel PE R1 custom FASTQ"
assert_custom_symlink \
    "${symlink_r2}" "execute download-fastqs GNU Parallel PE R2 custom FASTQ"

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
