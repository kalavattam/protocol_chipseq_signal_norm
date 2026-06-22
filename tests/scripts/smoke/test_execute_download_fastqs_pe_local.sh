#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_download_fastqs_pe_local.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="execute download-fastqs PE local"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

print_section "${TEST_NAME}"

if ! \
    is_download_enabled
then
    record_skip \
        "download-fastqs PE local check disabled; set RUN_DOWNLOAD=1 to" \
        "enable"
    finish
    exit $?
fi


#  Define fixture and output paths for a no-network loopback HTTP PE download
dir_fx="${ROOT_REPO}/tests/download_fastqs/fixtures"
dir_src="${dir_fx}/source"

mta_tpl="${dir_fx}/metadata/local_pe.template.tsv"
mta_tpl_dup="${dir_fx}/metadata/local_pe_duplicate.template.tsv"
mta_tpl_dup_cus="${dir_fx}/metadata/local_pe_duplicate_custom.template.tsv"
mta_tpl_cnfl="${dir_fx}/metadata/local_pe_conflicting_accession.template.tsv"

src_r1="${dir_src}/pe/tiny_download_pe_R1.fastq.gz"
src_r2="${dir_src}/pe/tiny_download_pe_R2.fastq.gz"

tmp="${TEST_DIR_TMP}/execute_download_fastqs_pe_local"

dir_out="${tmp}/out"
dir_sym="${tmp}/links"
dir_err="${tmp}/logs"
dir_out_dup="${tmp}/out_duplicate"
dir_sym_dup="${tmp}/links_duplicate"
dir_err_dup="${tmp}/logs_duplicate"
dir_out_fail="${tmp}/out_failure"
dir_sym_fail="${tmp}/links_failure"
dir_err_fail="${tmp}/logs_failure"
dir_log="${TEST_DIR_LOG}/download_fastqs"

mta="${tmp}/local_pe.tsv"
mta_dup="${tmp}/local_pe_duplicate.tsv"
mta_dup_cus="${tmp}/local_pe_duplicate_custom.tsv"
mta_cnfl="${tmp}/local_pe_conflicting_accession.tsv"

vw_r1="${tmp}/SRR_LOCAL_PE_R1.fastq"
vw_r2="${tmp}/SRR_LOCAL_PE_R2.fastq"
vw_dup_r1="${tmp}/SRR_LOCAL_PE_duplicate_R1.fastq"
vw_dup_r2="${tmp}/SRR_LOCAL_PE_duplicate_R2.fastq"

cnt_r1="${tmp}/SRR_LOCAL_PE_R1.read_count.txt"
cnt_r2="${tmp}/SRR_LOCAL_PE_R2.read_count.txt"
cnt_dup_r1="${tmp}/SRR_LOCAL_PE_duplicate_R1.read_count.txt"
cnt_dup_r2="${tmp}/SRR_LOCAL_PE_duplicate_R2.read_count.txt"

out_r1="${dir_out}/SRR_LOCAL_PE_R1.fastq.gz"
out_r2="${dir_out}/SRR_LOCAL_PE_R2.fastq.gz"
out_dup_r1="${dir_out_dup}/SRR_LOCAL_PE_R1.fastq.gz"
out_dup_r2="${dir_out_dup}/SRR_LOCAL_PE_R2.fastq.gz"

sym_r1="${dir_sym}/tiny_download_pe_R1.fastq.gz"
sym_r2="${dir_sym}/tiny_download_pe_R2.fastq.gz"
sym_dup_h18_r1="${dir_sym_dup}/tiny_download_pe_h3k18_R1.fastq.gz"
sym_dup_h18_r2="${dir_sym_dup}/tiny_download_pe_h3k18_R2.fastq.gz"
sym_dup_h27_r1="${dir_sym_dup}/tiny_download_pe_h3k27_R1.fastq.gz"
sym_dup_h27_r2="${dir_sym_dup}/tiny_download_pe_h3k27_R2.fastq.gz"

log_srvr="${dir_log}/execute_download_fastqs_pe_local_http.log"
log_run="${dir_log}/execute_download_fastqs_pe_local.log"

trap 'cleanup_server_http "${pid_srvr:-}"' EXIT
log_dup="${dir_log}/execute_download_fastqs_pe_local_duplicate.log"
log_dup_cus="${dir_log}/execute_download_fastqs_pe_local_duplicate_custom.log"
log_cnfl="${dir_log}/execute_download_fastqs_pe_local_conflicting_accession.log"
log_env="${dir_log}/execute_download_fastqs_pe_local_env.log"
log_out_r1="${dir_err}/test_execute_download_pe_local.SRR_LOCAL_PE_R1.stdout.txt"
log_err_r1="${dir_err}/test_execute_download_pe_local.SRR_LOCAL_PE_R1.stderr.txt"
log_out_r2="${dir_err}/test_execute_download_pe_local.SRR_LOCAL_PE_R2.stdout.txt"
log_err_r2="${dir_err}/test_execute_download_pe_local.SRR_LOCAL_PE_R2.stderr.txt"

rm -rf "${tmp}"
mkdir -p \
    "${dir_out}"      "${dir_sym}"      "${dir_err}" \
    "${dir_out_dup}"  "${dir_sym_dup}"  "${dir_err_dup}" \
    "${dir_out_fail}" "${dir_sym_fail}" "${dir_err_fail}" \
    "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty \
    "${src_r1}" \
    "${src_r2}" \
    "${mta_tpl}" \
    "${mta_tpl_dup}" \
    "${mta_tpl_dup_cus}" \
    "${mta_tpl_cnfl}" || {
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
url_src="${url_bas}/pe/tiny_download_pe_R1.fastq.gz"

(
    cd "${dir_src}" || exit 1
    "${py}" -m http.server "${port}" --bind 127.0.0.1
) > "${log_srvr}" 2>&1 &
pid_srvr=$!

if \
    wait_http_local "${py}" "${url_src}"
then
    record_pass "loopback HTTP fixture server is reachable"
else
    record_fail \
        "loopback HTTP fixture server did not become reachable; see" \
        "$(print_relpath "${log_srvr}")"
    finish
    exit $?
fi

sed "s#__BASE_URL__#${url_bas}#g" \
    "${mta_tpl}" > "${mta}"

sed "s#__BASE_URL__#${url_bas}#g" \
    "${mta_tpl_dup}" > "${mta_dup}"

sed "s#__BASE_URL__#${url_bas}#g" \
    "${mta_tpl_dup_cus}" > "${mta_dup_cus}"

sed "s#__BASE_URL__#${url_bas}#g" \
    "${mta_tpl_cnfl}" > "${mta_cnfl}"

require_files_nonempty \
    "${mta}" \
    "${mta_dup}" \
    "${mta_dup_cus}" \
    "${mta_cnfl}" || {
    finish
    exit $?
}


#  Run execute_download_fastqs.sh serially on one loopback HTTP PE URL pair
if \
    run_capture \
        "execute download-fastqs PE local HTTP wet run" \
        "${log_run}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/execute_download_fastqs.sh" \
            --threads 1 \
            --infile "${mta}" \
            --dir_out "${dir_out}" \
            --dir_sym "${dir_sym}" \
            --dir_eo "${dir_err}" \
            --nam_job "test_execute_download_pe_local"
then
    record_pass "execute_download_fastqs.sh PE local HTTP wet run exits 0"
else
    record_fail \
        "execute_download_fastqs.sh PE local HTTP wet run failed; see" \
        "$(print_relpath "${log_run}")"
fi

assert_fastq_gzip \
    "${out_r1}" \
    '^@tiny_download_pe_pair_1/1$' \
    "${vw_r1}" \
    "${cnt_r1}" \
    "${src_r1}" \
    "execute download-fastqs PE R1 accession FASTQ output"

assert_fastq_gzip \
    "${out_r2}" \
    '^@tiny_download_pe_pair_1/2$' \
    "${vw_r2}" \
    "${cnt_r2}" \
    "${src_r2}" \
    "execute download-fastqs PE R2 accession FASTQ output"

assert_custom_symlink \
    "${sym_r1}" \
    "execute download-fastqs PE R1 custom FASTQ"

assert_custom_symlink \
    "${sym_r2}" \
    "execute download-fastqs PE R2 custom FASTQ"

assert_file_exists \
    "${log_out_r1}" \
    "execute download-fastqs PE R1 wget stdout log exists"

assert_file_exists \
    "${log_err_r1}" \
    "execute download-fastqs PE R1 wget stderr log exists"

assert_file_exists \
    "${log_out_r2}" \
    "execute download-fastqs PE R2 wget stdout log exists"

assert_file_exists \
    "${log_err_r2}" \
    "execute download-fastqs PE R2 wget stderr log exists"

#  Duplicate accession rows should download once and create all custom links
if \
    run_capture \
        "execute download-fastqs PE duplicate accession local HTTP wet run" \
        "${log_dup}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/execute_download_fastqs.sh" \
            --threads 1 \
            --infile "${mta_dup}" \
            --dir_out "${dir_out_dup}" \
            --dir_sym "${dir_sym_dup}" \
            --dir_eo "${dir_err_dup}" \
            --nam_job "test_execute_download_pe_duplicate"
then
    record_pass "execute_download_fastqs.sh PE duplicate accession exits 0"
else
    record_fail \
        "execute_download_fastqs.sh PE duplicate accession failed; see" \
        "$(print_relpath "${log_dup}")"
fi

assert_fastq_gzip \
    "${out_dup_r1}" \
    '^@tiny_download_pe_pair_1/1$' \
    "${vw_dup_r1}" \
    "${cnt_dup_r1}" \
    "${src_r1}" \
    "execute download-fastqs duplicate PE R1 accession FASTQ output"

assert_fastq_gzip \
    "${out_dup_r2}" \
    '^@tiny_download_pe_pair_1/2$' \
    "${vw_dup_r2}" \
    "${cnt_dup_r2}" \
    "${src_r2}" \
    "execute download-fastqs duplicate PE R2 accession FASTQ output"

assert_custom_symlink \
    "${sym_dup_h18_r1}" \
    "execute download-fastqs duplicate PE H3K18ac R1 custom FASTQ"

assert_custom_symlink \
    "${sym_dup_h18_r2}" \
    "execute download-fastqs duplicate PE H3K18ac R2 custom FASTQ"

assert_custom_symlink \
    "${sym_dup_h27_r1}" \
    "execute download-fastqs duplicate PE H3K27ac R1 custom FASTQ"

assert_custom_symlink \
    "${sym_dup_h27_r2}" \
    "execute download-fastqs duplicate PE H3K27ac R2 custom FASTQ"

if [[ "$(readlink "${sym_dup_h18_r1}")" == "${out_dup_r1}" \
    && "$(readlink "${sym_dup_h27_r1}")" == "${out_dup_r1}" ]]
then
    record_pass "duplicate PE R1 custom names point to one raw download"
else
    record_fail "duplicate PE R1 custom names do not point to one raw download"
fi

num_logs="$(
    find "${dir_err_dup}" -name '*.stdout.txt' -type f | wc -l | tr -d ' '
)"
if [[ "${num_logs}" == "2" ]]; then
    record_pass "duplicate PE accession creates one wget stdout log per mate"
else
    record_fail "duplicate PE accession wrote unexpected wget stdout log count"
fi

if \
    run_capture \
        "execute download-fastqs PE duplicate custom-name expected failure" \
        "${log_dup_cus}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/execute_download_fastqs.sh" \
            --threads 1 \
            --infile "${mta_dup_cus}" \
            --dir_out "${dir_out_fail}" \
            --dir_sym "${dir_sym_fail}" \
            --dir_eo "${dir_err_fail}" \
            --nam_job "test_execute_download_pe_duplicate_custom"
then
    record_fail "execute_download_fastqs.sh accepted duplicate custom_name"
else
    assert_pattern_found \
        "${log_dup_cus}" \
        "duplicate custom_name" \
        "execute_download_fastqs.sh rejects duplicate custom_name"
fi

if \
    run_capture \
        "execute download-fastqs PE conflicting accession expected failure" \
        "${log_cnfl}" \
        "${TEST_BASH}" "${ROOT_REPO}/scripts/execute_download_fastqs.sh" \
            --threads 1 \
            --infile "${mta_cnfl}" \
            --dir_out "${dir_out_fail}" \
            --dir_sym "${dir_sym_fail}" \
            --dir_eo "${dir_err_fail}" \
            --nam_job "test_execute_download_pe_conflict"
then
    record_fail \
        "execute_download_fastqs.sh accepted conflicting accession URLs"
else
    assert_pattern_found \
        "${log_cnfl}" \
        "conflicting.*FASTQ URL" \
        "execute_download_fastqs.sh rejects conflicting accession URLs"
fi

finish
