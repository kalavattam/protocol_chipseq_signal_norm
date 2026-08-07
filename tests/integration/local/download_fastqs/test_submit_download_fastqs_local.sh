#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_submit_download_fastqs_local.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="submit download-fastqs local"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


if ! \
    is_download_enabled
then
    record_skip \
        "submit download-fastqs local check disabled; set RUN_DOWNLOAD=1 to" \
        "enable"
    finish
    exit $?
fi


#  Define fixture and output paths for direct loopback HTTP downloads
dir_fx="${ROOT_REPO}/tests/fixtures/download_fastqs"
dir_src="${dir_fx}/source"
src_se="${dir_src}/se/tiny_download_se.fastq.gz"
src_r1="${dir_src}/pe/tiny_download_pe_R1.fastq.gz"
src_r2="${dir_src}/pe/tiny_download_pe_R2.fastq.gz"

tmp="${TEST_DIR_TMP}/submit_download_fastqs_local"
dir_out_se="${tmp}/se/out"
dir_sym_se="${tmp}/se/links"
dir_eo_se="${tmp}/se/logs"
dir_out_pe="${tmp}/pe/out"
dir_sym_pe="${tmp}/pe/links"
dir_eo_pe="${tmp}/pe/logs"

srr_se="SRR_SUBMIT_LOCAL_SE"
nam_cus_se="submit_local_canonical"
nam_job_se="test_submit_download_se_local"
out_se="${dir_out_se}/${srr_se}.fastq.gz"
sym_se="${dir_sym_se}/${nam_cus_se}.fastq.gz"
log_out_se="${dir_eo_se}/${nam_job_se}.${srr_se}.stdout.txt"
log_err_se="${dir_eo_se}/${nam_job_se}.${srr_se}.stderr.txt"

srr_pe="SRR_SUBMIT_LOCAL_PE"
nam_cus_pe="submit_local_hidden"
nam_job_pe="test_submit_download_pe_local"
out_r1="${dir_out_pe}/${srr_pe}_R1.fastq.gz"
out_r2="${dir_out_pe}/${srr_pe}_R2.fastq.gz"
sym_r1="${dir_sym_pe}/${nam_cus_pe}_R1.fastq.gz"
sym_r2="${dir_sym_pe}/${nam_cus_pe}_R2.fastq.gz"
log_pfx_pe="${dir_eo_pe}/${nam_job_pe}.${srr_pe}"
log_out_r1="${log_pfx_pe}_R1.stdout.txt"
log_err_r1="${log_pfx_pe}_R1.stderr.txt"
log_out_r2="${log_pfx_pe}_R2.stdout.txt"
log_err_r2="${log_pfx_pe}_R2.stderr.txt"

log_svr="${tmp}/loopback_http.log"
log_run_se="${tmp}/submit_canonical_se.log"
log_run_pe="${tmp}/submit_hidden_pe.log"
log_env="${tmp}/submit_download_fastqs_local_env.log"


print_section "${TEST_NAME}"

trap \
    'cleanup_server_http "${svr_pid:-}"' \
    EXIT

rm -rf "${tmp}"
mkdir -p \
    "${dir_out_se}" "${dir_sym_se}" "${dir_eo_se}" \
    "${dir_out_pe}" "${dir_sym_pe}" "${dir_eo_pe}"

require_files_nonempty \
    "${src_se}" \
    "${src_r1}" \
    "${src_r2}" || {
    finish
    exit $?
}

require_env_project env_nam || {
    finish
    exit $?
}

# shellcheck disable=2154
if ! \
    require_env_download \
        "${env_nam}" \
        "${log_env}"
then
    finish
    exit $?
fi

arr_env_cmd=()
if [[ "${CONDA_DEFAULT_ENV:-}" != "${env_nam}" ]]; then
    arr_env_cmd=( conda run -n "${env_nam}" )
fi

if \
    py="$(find_python_loopback)"
then
    record_pass "host Python is usable for loopback fixture server"
    record_pass "env_protocol Python fallback not selected for fixture server"
else
    record_pass "host Python is not usable for loopback fixture server"

    if ! \
        py="$(
            conda run -n env_protocol python -c \
                'import sys; print(sys.executable)' \
                2>> "${log_env}"
        )" || [[
            "${py}" != /* || ! -x "${py}"
        ]]
    then
        record_fail \
            "failed to resolve fixture-server Python from env_protocol; see" \
            "$(print_relpath "${log_env}")"
        finish
        exit $?
    fi

    record_pass "env_protocol Python fallback selected for fixture server"
fi

if ! \
    port="$(find_port_free "${py}")"
then
    record_fail "selected fixture-server Python cannot bind 127.0.0.1"
    finish
    exit $?
fi

record_pass "selected fixture-server Python can bind 127.0.0.1"

url_bas="http://127.0.0.1:${port}"
url_se="${url_bas}/se/tiny_download_se.fastq.gz"
url_r1="${url_bas}/pe/tiny_download_pe_R1.fastq.gz"
url_r2="${url_bas}/pe/tiny_download_pe_R2.fastq.gz"

(
    cd "${dir_src}" || exit 1
    exec env PYTHONDONTWRITEBYTECODE=1 \
        "${py}" -m http.server "${port}" --bind 127.0.0.1
) > "${log_svr}" 2>&1 &
svr_pid=$!

if \
    wait_http_local "${py}" "${url_se}"
then
    record_pass "loopback-only HTTP fixture server is reachable"
    record_pass "fixture HTTP server is bound only to 127.0.0.1"
else
    record_fail \
        "loopback HTTP fixture server did not become reachable; see" \
        "$(print_relpath "${log_svr}")"
    finish
    exit $?
fi

if \
    kill -0 "${svr_pid}" > /dev/null 2>&1
then
    record_pass "loopback HTTP fixture server remains running"
else
    record_fail "loopback HTTP fixture server exited before worker execution"
    finish
    exit $?
fi


#  Run direct submit with canonical --dir_scr on one SE loopback URL
if \
    run_capture \
        "submit download-fastqs canonical SE local HTTP wet run" \
        "${log_run_se}" \
        "${arr_env_cmd[@]}" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/submit_download_fastqs.sh" \
            --dir_scr "${ROOT_REPO}/bin" \
            "${srr_se}" \
            "${url_se}" \
            NA \
            "${dir_out_se}" \
            "${dir_sym_se}" \
            "${nam_cus_se}" \
            "${dir_eo_se}" \
            "${nam_job_se}"
then
    record_pass "submit canonical --dir_scr SE worker exits 0"
else
    record_fail \
        "submit canonical --dir_scr SE worker failed; see" \
        "$(print_relpath "${log_run_se}")"
fi

assert_files_equal \
    "${out_se}" \
    "${src_se}" \
    "submit canonical SE downloaded FASTQ"

assert_custom_symlink \
    "${sym_se}" \
    "submit canonical SE custom FASTQ"

if [[ "${sym_se}" -ef "${out_se}" ]]; then
    record_pass "submit canonical SE symlink resolves to accession output"
else
    record_fail "submit canonical SE symlink resolves to the wrong target"
fi

assert_file_exists \
    "${log_out_se}" \
    "submit canonical SE wget stdout log exists"

if [[ -f "${log_out_se}" && ! -s "${log_out_se}" ]]; then
    record_pass "submit canonical SE wget stdout log is empty"
else
    record_fail "submit canonical SE wget stdout log is unexpectedly nonempty"
fi

assert_file_nonempty \
    "${log_err_se}" \
    "submit canonical SE wget stderr log"

assert_pattern_found \
    "${log_run_se}" \
    "Downloading ${srr_se} from ${url_se}." \
    "submit canonical SE worker reports the mapped accession and URL"

assert_pattern_found \
    "${log_run_se}" \
    "Symlinking ${srr_se} to ${nam_cus_se}." \
    "submit canonical SE worker reports the mapped custom name"

if [[ \
        ! -e "${dir_out_se}/${srr_se}_R2.fastq.gz" \
    && ! -e "${dir_sym_se}/${nam_cus_se}_R2.fastq.gz" \
    && ! -e "${dir_eo_se}/${nam_job_se}.${srr_se}_R2.stdout.txt" \
    && ! -e "${dir_eo_se}/${nam_job_se}.${srr_se}_R2.stderr.txt" \
]]; then
    record_pass "submit canonical SE worker creates no second mate"
else
    record_fail "submit canonical SE worker unexpectedly creates a second mate"
fi

shopt -s nullglob
arr_out_se=( "${dir_out_se}"/* )
arr_sym_se=( "${dir_sym_se}"/* )
arr_eo_se=( "${dir_eo_se}"/* )
shopt -u nullglob

if ((
        ${#arr_out_se[@]} == 1
    && ${#arr_sym_se[@]} == 1
    && ${#arr_eo_se[@]} == 2
)); then
    record_pass "submit canonical SE worker creates only expected files"
else
    record_fail "submit canonical SE worker creates unexpected files"
fi


#  Run direct submit with hidden --dir-scr on one PE loopback URL pair
if \
    run_capture \
        "submit download-fastqs hidden PE local HTTP wet run" \
        "${log_run_pe}" \
        "${arr_env_cmd[@]}" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/submit_download_fastqs.sh" \
            --dir-scr "${ROOT_REPO}/bin" \
            "${srr_pe}" \
            "${url_r1}" \
            "${url_r2}" \
            "${dir_out_pe}" \
            "${dir_sym_pe}" \
            "${nam_cus_pe}" \
            "${dir_eo_pe}" \
            "${nam_job_pe}"
then
    record_pass "submit hidden --dir-scr PE worker exits 0"
else
    record_fail \
        "submit hidden --dir-scr PE worker failed; see" \
        "$(print_relpath "${log_run_pe}")"
fi

assert_files_equal \
    "${out_r1}" \
    "${src_r1}" \
    "submit hidden PE R1 downloaded FASTQ"

assert_files_equal \
    "${out_r2}" \
    "${src_r2}" \
    "submit hidden PE R2 downloaded FASTQ"

assert_custom_symlink \
    "${sym_r1}" \
    "submit hidden PE R1 custom FASTQ"

assert_custom_symlink \
    "${sym_r2}" \
    "submit hidden PE R2 custom FASTQ"

if [[ "${sym_r1}" -ef "${out_r1}" ]]; then
    record_pass "submit hidden PE R1 symlink resolves to accession output"
else
    record_fail "submit hidden PE R1 symlink resolves to the wrong target"
fi

if [[ "${sym_r2}" -ef "${out_r2}" ]]; then
    record_pass "submit hidden PE R2 symlink resolves to accession output"
else
    record_fail "submit hidden PE R2 symlink resolves to the wrong target"
fi

for log_out in "${log_out_r1}" "${log_out_r2}"; do
    assert_file_exists \
        "${log_out}" \
        "submit hidden PE wget stdout log exists"

    if [[ -f "${log_out}" && ! -s "${log_out}" ]]; then
        record_pass "submit hidden PE wget stdout log is empty"
    else
        record_fail "submit hidden PE wget stdout log is unexpectedly nonempty"
    fi
done

assert_file_nonempty \
    "${log_err_r1}" \
    "submit hidden PE R1 wget stderr log"

assert_file_nonempty \
    "${log_err_r2}" \
    "submit hidden PE R2 wget stderr log"

assert_pattern_found \
    "${log_run_pe}" \
    "Downloading ${srr_pe} from ${url_r1}." \
    "submit hidden PE worker reports the mapped R1 URL"

assert_pattern_found \
    "${log_run_pe}" \
    "Downloading ${srr_pe} from ${url_r2}." \
    "submit hidden PE worker reports the mapped R2 URL"

assert_pattern_found \
    "${log_run_pe}" \
    "Symlinking ${srr_pe} to ${nam_cus_pe}." \
    "submit hidden PE worker reports the mapped custom name"

if [[ \
    ! -e "${dir_out_pe}/${srr_pe}.fastq.gz" \
    && ! -e "${dir_sym_pe}/${nam_cus_pe}.fastq.gz" \
]]; then
    record_pass "submit hidden PE worker creates no unsuffixed SE output"
else
    record_fail "submit hidden PE worker unexpectedly creates SE output"
fi

shopt -s nullglob
arr_out_pe=( "${dir_out_pe}"/* )
arr_sym_pe=( "${dir_sym_pe}"/* )
arr_eo_pe=( "${dir_eo_pe}"/* )
shopt -u nullglob

if ((
        ${#arr_out_pe[@]} == 2
    && ${#arr_sym_pe[@]} == 2
    && ${#arr_eo_pe[@]} == 4
)); then
    record_pass "submit hidden PE worker creates only expected files"
else
    record_fail "submit hidden PE worker creates unexpected files"
fi

assert_file_nonempty \
    "${log_svr}" \
    "loopback HTTP server request log"

finish
