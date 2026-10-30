#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_download_fastqs_names_local.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="execute download-fastqs read names local"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


if ! \
    is_download_enabled
then
    record_skip \
        "download-fastqs read-name local check disabled; set" \
        "RUN_DOWNLOAD=1 to enable"
    finish
    exit $?
fi


# Define fixture and output paths for loopback-only HTTP downloads.
dir_fx="${ROOT_REPO}/tests/fixtures/download_fastqs"
dir_src="${dir_fx}/source"
dir_nam="${dir_src}/names"
dir_mtd="${dir_fx}/metadata"
dir_exp="${dir_fx}/expected"

tmp="${TEST_DIR_TMP}/execute_download_fastqs_names_local"
dir_log="${TEST_DIR_LOG}/download_fastqs"
log_svr="${dir_log}/execute_download_fastqs_names_local_http.log"
log_env="${dir_log}/execute_download_fastqs_names_local_env.log"


print_section "${TEST_NAME}"

trap 'cleanup_server_http "${svr_pid:-}"' EXIT

rm -rf "${tmp}"
mkdir -p "${tmp}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty \
    "${dir_nam}/ena_R1.fastq.gz" \
    "${dir_nam}/ena_R2.fastq.gz" \
    "${dir_exp}/names_pe.fastq" \
    "${dir_exp}/names_se.fastq" \
    "${dir_mtd}/local_names_pe_missing.template.tsv" || {
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
    wait_http_local "${py}" "${url_bas}/names/ena_R1.fastq.gz"
then
    record_pass "loopback HTTP fixture server is reachable"
else
    record_fail \
        "loopback HTTP fixture server did not become reachable; see" \
        "$(print_relpath "${log_svr}")"
    finish
    exit $?
fi


# Run one template into the case directory 'dir', which may already hold files;
# the log is named for 'nam'.
function run_template() {
    local nam="${1}"
    local dir="${2}"
    local tpl="${3}"
    shift 3

    mkdir -p "${dir}/out" "${dir}/links" "${dir}/logs"
    sed "s#__BASE_URL__#${url_bas}#g" \
        "${dir_mtd}/${tpl}.template.tsv" > "${dir}/input.tsv"

    run_capture \
        "${nam}" \
        "${dir_log}/execute_download_fastqs_names_${nam}.log" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/execute_download_fastqs.sh" \
            --env_nam "${env_nam}" \
            --fil_in "${dir}/input.tsv" \
            --dir_out "${dir}/out" \
            --dir_sym "${dir}/links" \
            --dir_eo "${dir}/logs" \
            --nam_job "names_${nam}" \
            "$@"
}


# Record whether a case exited as expected.
function check_rc() {
    local nam="${1}"
    local rc="${2}"
    local rc_exp="${3}"
    local log="${dir_log}/execute_download_fastqs_names_${nam}.log"

    if (( rc == rc_exp )); then
        record_pass "${nam} exits as expected (${rc})"
    else
        record_fail \
            "${nam} exited ${rc}, expected ${rc_exp}; see" \
            "$(print_relpath "${log}")"
    fi
}


# Record whether a gzipped FASTQ decompresses to the expected records.
function check_restored() {
    local nam="${1}"
    local fil="${2}"
    local exp="${3}"

    if \
        gzip -cd < "${fil}" | cmp -s - "${exp}"
    then
        record_pass "${nam} restores $(basename "${fil}") as expected"
    else
        record_fail \
            "${nam} did not restore $(print_relpath "${fil}") to" \
            "$(print_relpath "${exp}")"
    fi
}


# Record whether a file is byte-identical to its source.
function check_unchanged() {
    local nam="${1}"
    local fil="${2}"
    local src="${3}"

    if \
        cmp -s "${fil}" "${src}"
    then
        record_pass "${nam} leaves $(basename "${fil}") as downloaded"
    else
        record_fail \
            "${nam} changed $(print_relpath "${fil}"); expected the bytes" \
            "of $(print_relpath "${src}")"
    fi
}


# Record whether a case left no restoration temporaries.
function check_no_tmp() {
    local nam="${1}"
    local dir="${2}"
    local fil_tmp

    fil_tmp="$(find "${dir}" -name '*.names.tmp.gz' -print)"

    if [[ -z "${fil_tmp}" ]]; then
        record_pass "${nam} leaves no restoration temporaries"
    else
        record_fail "${nam} left restoration temporaries: ${fil_tmp}"
    fi
}


# Restore names for good pairs and single-end files, with and without '/1'.
pe_1="SRR_NAMES_PE_R1.fastq.gz"
pe_2="SRR_NAMES_PE_R2.fastq.gz"
se="SRR_NAMES_SE.fastq.gz"

rc=0
run_template pe "${tmp}/pe" local_names_pe || rc=$?
check_rc pe "${rc}" 0
check_restored pe "${tmp}/pe/out/${pe_1}" "${dir_exp}/names_pe.fastq"
check_restored pe "${tmp}/pe/out/${pe_2}" "${dir_exp}/names_pe.fastq"
check_no_tmp pe "${tmp}/pe/out"
assert_pattern_absent \
    "${dir_log}/execute_download_fastqs_names_pe.log" \
    "read-name form" \
    "pe gives no read-name warning"

for nam in se se_no_mate; do
    if [[ "${nam}" == "se" ]]; then
        tpl=local_names_se
    else
        tpl=local_names_se_no_mate
    fi

    rc=0
    run_template "${nam}" "${tmp}/${nam}" "${tpl}" || rc=$?
    check_rc "${nam}" "${rc}" 0
    check_restored \
        "${nam}" \
        "${tmp}/${nam}/out/${se}" \
        "${dir_exp}/names_se.fastq"
done

# '--keep_names' leaves both mates as downloaded.
rc=0
run_template pe_keep "${tmp}/pe_keep" local_names_pe --keep_names || rc=$?
check_rc pe_keep "${rc}" 0
check_unchanged \
    pe_keep \
    "${tmp}/pe_keep/out/${pe_1}" \
    "${dir_nam}/ena_R1.fastq.gz"
check_unchanged \
    pe_keep \
    "${tmp}/pe_keep/out/${pe_2}" \
    "${dir_nam}/ena_R2.fastq.gz"
assert_pattern_absent \
    "${dir_log}/execute_download_fastqs_names_pe_keep.log" \
    'Restoring the read names' \
    "pe_keep does not restore names"

# Each bad pair breaks one check, so both mates stay as downloaded and one
# warning counts the failed record.
wrn="1 of 2 records do not have the archive's read-name form"
arr_rows=(
    "pe_r2_name|local_names_pe_r2_name|ena_R1|r2_name_R2"
    "pe_order|local_names_pe_order|order_R1|order_R2"
    "pe_mate|local_names_pe_mate|mate_R1|ena_R2"
)

for row in "${arr_rows[@]}"; do
    IFS='|' read -r nam tpl src_1 src_2 <<< "${row}"

    rc=0
    run_template "${nam}" "${tmp}/${nam}" "${tpl}" || rc=$?
    check_rc "${nam}" "${rc}" 0
    check_unchanged \
        "${nam}" \
        "${tmp}/${nam}/out/${pe_1}" \
        "${dir_nam}/${src_1}.fastq.gz"
    check_unchanged \
        "${nam}" \
        "${tmp}/${nam}/out/${pe_2}" \
        "${dir_nam}/${src_2}.fastq.gz"
    check_no_tmp "${nam}" "${tmp}/${nam}/out"
    assert_pattern_found \
        "${dir_log}/execute_download_fastqs_names_${nam}.log" \
        "${wrn}" \
        "${nam} warns once with the count"
done

# Names the archive did not store stay as downloaded, with one note.
rc=0
run_template se_lost "${tmp}/se_lost" local_names_se_lost || rc=$?
check_rc se_lost "${rc}" 0
check_unchanged \
    se_lost \
    "${tmp}/se_lost/out/SRR_NAMES_LOST.fastq.gz" \
    "${dir_nam}/lost_se.fastq.gz"
assert_pattern_found \
    "${dir_log}/execute_download_fastqs_names_se_lost.log" \
    '^note(.*read names were not stored by the archive' \
    "se_lost notes that the archive did not store the names"

# Names straight from the sequencer need nothing: the files stay as downloaded,
# with one note, and are not read past their first header.
for nam in seq_pe seq_se; do
    if [[ "${nam}" == "seq_pe" ]]; then
        tpl=local_names_sequencer_pe
        arr_fil=( SEQ_PE_R1 SEQ_PE_R2 )
        arr_src=( sequencer_R1 sequencer_R2 )
    else
        tpl=local_names_sequencer_se
        arr_fil=( SEQ_SE )
        arr_src=( sequencer_se )
    fi

    rc=0
    run_template "${nam}" "${tmp}/${nam}" "${tpl}" || rc=$?
    check_rc "${nam}" "${rc}" 0

    for idx in "${!arr_fil[@]}"; do
        check_unchanged \
            "${nam}" \
            "${tmp}/${nam}/out/${arr_fil[idx]}.fastq.gz" \
            "${dir_nam}/${arr_src[idx]}.fastq.gz"
    done

    assert_pattern_found \
        "${dir_log}/execute_download_fastqs_names_${nam}.log" \
        "^note(.*read names are not in the archive's '<run>.<n>' form" \
        "${nam} notes that the names are left as is"
    assert_pattern_absent \
        "${dir_log}/execute_download_fastqs_names_${nam}.log" \
        'Restoring the read names' \
        "${nam} does not read the files to restore names"
done

# Sequencer files in place keep their MD5 check under '--use_existing': a pair
# whose mates are swapped fails it without a download.
dir="${tmp}/use_existing_seq"
mkdir -p "${dir}/out"
cp \
    "${dir_nam}/sequencer_R2.fastq.gz" \
    "${dir}/out/SEQ_PE_R1.fastq.gz"  # fixture-perturbation: mates swapped
cp \
    "${dir_nam}/sequencer_R1.fastq.gz" \
    "${dir}/out/SEQ_PE_R2.fastq.gz"  # fixture-perturbation: mates swapped

rc=0
run_template ue_seq "${dir}" local_names_sequencer_pe --use_existing || rc=$?
check_rc ue_seq "${rc}" 1
assert_pattern_found \
    "${dir_log}/execute_download_fastqs_names_ue_seq.log" \
    "existing file '.*SEQ_PE_R1.fastq.gz' fails its check: MD5 is" \
    "ue_seq checks the MD5 of sequencer files in place"
assert_pattern_absent \
    "${dir_log}/execute_download_fastqs_names_ue_seq.log" \
    '^Downloading ' \
    "ue_seq downloads nothing"

# '--use_existing' checks files already in place against the TSV's MD5 values
# and sizes, whose URLs the server lacks, and restores them; a rerun notes that
# the names are already restored and changes nothing.
dir="${tmp}/use_existing"
mkdir -p "${dir}/out"
cp \
    "${dir_nam}/ena_R1.fastq.gz" \
    "${dir}/out/${pe_1}"  # fixture-perturbation: restored
cp \
    "${dir_nam}/ena_R2.fastq.gz" \
    "${dir}/out/${pe_2}"  # fixture-perturbation: restored

rc=0
run_template ue "${dir}" local_names_pe_missing --use_existing || rc=$?
check_rc ue "${rc}" 0
check_restored ue "${dir}/out/${pe_1}" "${dir_exp}/names_pe.fastq"
check_restored ue "${dir}/out/${pe_2}" "${dir_exp}/names_pe.fastq"
assert_pattern_absent \
    "${dir_log}/execute_download_fastqs_names_ue.log" \
    '^Downloading ' \
    "ue downloads nothing"

md5_pre="$(cat "${dir}/out/${pe_1}" "${dir}/out/${pe_2}" | md5sum)"

rc=0
run_template ue_rerun "${dir}" local_names_pe_missing --use_existing || rc=$?
check_rc ue_rerun "${rc}" 0
assert_pattern_found \
    "${dir_log}/execute_download_fastqs_names_ue_rerun.log" \
    '^note(.*read names are already restored' \
    "ue_rerun notes that the names are already restored"

if [[
    "$(cat "${dir}/out/${pe_1}" "${dir}/out/${pe_2}" | md5sum)" == "${md5_pre}"
]]; then
    record_pass "ue_rerun changes no file"
else
    record_fail "ue_rerun changed a restored file"
fi

# An existing file that fails its check fails the entry without a download.
dir="${tmp}/use_existing_bad"
mkdir -p "${dir}/out"
cp \
    "${dir_nam}/ena_R1.fastq.gz" \
    "${dir}/out/${pe_1}"  # fixture-perturbation: restored
cp \
    "${dir_nam}/r2_name_R2.fastq.gz" \
    "${dir}/out/${pe_2}"  # fixture-perturbation: an R2 name unlike R1's

rc=0
run_template ue_bad "${dir}" local_names_pe_missing --use_existing || rc=$?
check_rc ue_bad "${rc}" 1
check_unchanged ue_bad "${dir}/out/${pe_2}" "${dir_nam}/r2_name_R2.fastq.gz"
assert_pattern_found \
    "${dir_log}/execute_download_fastqs_names_ue_bad.log" \
    "existing file '.*${pe_2}' fails its check: MD5 is" \
    "ue_bad reports the failed check"
assert_pattern_absent \
    "${dir_log}/execute_download_fastqs_names_ue_bad.log" \
    '^Downloading ' \
    "ue_bad downloads nothing"

# A missing file is downloaded as usual under '--use_existing'.
rc=0
run_template ue_missing "${tmp}/ue_missing" local_names_pe --use_existing \
    || rc=$?
check_rc ue_missing "${rc}" 0
check_restored \
    ue_missing \
    "${tmp}/ue_missing/out/${pe_1}" \
    "${dir_exp}/names_pe.fastq"
assert_pattern_found \
    "${dir_log}/execute_download_fastqs_names_ue_missing.log" \
    '^Downloading SRR_NAMES_PE ' \
    "ue_missing downloads the missing files"

# A dry run forwards both flags to the worker and writes nothing.
rc=0
run_template dry "${tmp}/dry" local_names_pe \
    --dry_run \
    --use_existing \
    --keep_names \
    || rc=$?
check_rc dry "${rc}" 0
assert_pattern_found \
    "${dir_log}/execute_download_fastqs_names_dry.log" \
    'submit_download_fastqs\.sh .*--use_existing --keep_names ' \
    "dry run passes --use_existing and --keep_names to the worker"

if [[ -z "$(find "${tmp}/dry/out" -name '*.fastq.gz' -print)" ]]; then
    record_pass "dry run writes no FASTQ"
else
    record_fail "dry run wrote FASTQ files under $(print_relpath "${tmp}/dry")"
fi

# An empty file fails restoration instead of being taken for lost names.
rc=0
# shellcheck disable=SC2016  # Expand in the child shell, not this one.
run_capture \
    empty \
    "${dir_log}/execute_download_fastqs_names_empty.log" \
    "${TEST_BASH}" -c '
        # shellcheck source=lib/bash/core/source_helpers.sh
        source "${1}/lib/bash/core/source_helpers.sh"

        source_helpers "${1}/lib/bash" \
            core/check_args \
            core/check_inputs \
            core/format_outputs \
            workflows/download_fastqs
        restore_fastq_names SRR_EMPTY "${2}" NA
    ' _ "${ROOT_REPO}" "${dir_nam}/empty.fastq.gz" \
    || rc=$?
check_rc empty "${rc}" 1
assert_pattern_found \
    "${dir_log}/execute_download_fastqs_names_empty.log" \
    'FASTQ file is missing or empty' \
    "an empty file fails restoration"

# GNU Parallel workers restore names too.
if \
    is_parallel_enabled
then
    rc=0
    run_template parallel "${tmp}/parallel" local_names_pe --threads 2 \
        || rc=$?
    check_rc parallel "${rc}" 0
    check_restored \
        parallel \
        "${tmp}/parallel/out/${pe_1}" \
        "${dir_exp}/names_pe.fastq"
    check_restored \
        parallel \
        "${tmp}/parallel/out/${pe_2}" \
        "${dir_exp}/names_pe.fastq"
else
    record_skip \
        "GNU Parallel read-name row disabled; set RUN_PARALLEL=1 to enable"
fi


finish
