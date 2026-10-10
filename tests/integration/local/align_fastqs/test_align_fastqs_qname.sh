#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_align_fastqs_qname.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="align-fastqs queryname copy"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Define fixture and output paths. '--qname' name-sorts the final
# duplicate-marked alignments, so the copy holds the same records and flags.
dir_fx="${ROOT_REPO}/tests/fixtures/align_fastqs"
ref_rpt="${dir_fx}/reference/repeat.fa"
idx_bt2="${dir_fx}/bowtie2/repeat"
idx_bwa="${dir_fx}/bwa/repeat.fa"

tmp="${TEST_DIR_TMP}/align_fastqs_qname"
ref_fa="${tmp}/repeat.fa"
fq_1="${tmp}/qn_R1.fastq"
fq_2="${tmp}/qn_R2.fastq"
len_rd=50


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${tmp}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty \
    "${ref_rpt}" \
    "${idx_bt2}.1.bt2" \
    "${idx_bwa}.bwt" || {
    finish
    exit $?
}

# CRAM writing indexes the reference, so use a copy outside the fixtures.
cp "${ref_rpt}" "${ref_fa}"


# Print one FASTQ record.
function print_fq() {
    printf '@%s\n%s\n+\n%s\n' "${1}" "${2}" "${2//?/I}"
}


# Write one pair: R1 starts at a 0-based offset on the forward strand, and R2
# is the reverse complement of the fragment's last bases.
function print_pair() {
    local nam="${1}"
    local beg="${2}"
    local frg="${3}"
    local seq_2

    seq_2="$(
        rev <<< "${seq_ref:beg + frg - len_rd:len_rd}" | tr ACGT TGCA
    )"
    print_fq "${nam}" "${seq_ref:beg:len_rd}" >> "${fq_1}"
    print_fq "${nam}" "${seq_2}" >> "${fq_2}"
}


# Three 300 bp fragments from the unique blocks; 'dup_a' and 'dup_b' are the
# same fragment under two names, so Step 4 marks one of them.
seq_ref="$(grep -v '^>' "${ref_rpt}" | tr -d '\n')"

print_pair uniq 50 300
print_pair dup_a 1200 300
print_pair dup_b 1200 300


# Run the helper.
function run_helper() {
    # shellcheck disable=SC2016  # Expand in the child shell, not this one.
    "${TEST_BASH}" -c '
        # shellcheck source=lib/bash/core/source_helpers.sh
        source "${1}/lib/bash/core/source_helpers.sh"

        source_helpers "${1}/lib/bash" \
            core/check_args \
            core/check_inputs \
            core/format_outputs \
            workflows/align_fastqs
        shift 1
        align_fastqs "$@"
    ' _ "${ROOT_REPO}" --threads 1 "$@" 2>&1
}


# Print an alignment file's records, sorted as text.
function print_recs() {
    run_samtools view -T "${ref_fa}" "${1}" | LC_ALL=C sort
}


# label|arguments|output format
a_bt2="--aligner bowtie2 --index ${idx_bt2} --fq_2 ${fq_2}"
rows=(
    "Bowtie 2 PE|${a_bt2}|bam"
    "Bowtie 2 PE, CRAM|${a_bt2}|cram"
    "BWA-MEM SE|--aligner bwa --index ${idx_bwa}|bam"
)

for idx in "${!rows[@]}"; do
    IFS='|' read -r lbl arg ext <<< "${rows[idx]}"
    read -r -a args <<< "${arg}"
    dir_out="${tmp}/out_${idx}"
    fil_out="${dir_out}/x.${ext}"
    fil_qnm="${dir_out}/x.qnam.${ext}"
    ref=()
    if [[ "${ext}" == "cram" ]]; then ref=( --ref_fa "${ref_fa}" ); fi

    mkdir -p "${dir_out}"
    rc=0
    out="$(
        run_helper \
            "${args[@]}" \
            --fq_1 "${fq_1}" \
            "${ref[@]}" \
            --qname \
            --fil_out "${fil_out}"
    )" || rc=$?

    if [[ "${rc}" -ne 0 ]]; then
        record_fail "${lbl}: exit ${rc}"
        printf '%s\n' "${out}" >&2
        continue
    fi

    # Only the final file, its index, the duplicate statistics, and the copy
    # remain.
    got="$(
        cd "${dir_out}" && find . -type f | LC_ALL=C sort | paste -sd ' ' -
    )"
    idx_ext="bai"
    if [[ "${ext}" == "cram" ]]; then idx_ext="crai"; fi
    exp="./x.${ext} ./x.${ext}.${idx_ext} ./x.markdup.txt.gz ./x.qnam.${ext}"

    if [[ "${got}" == "${exp}" ]]; then
        record_pass \
            "${lbl}: writes the final file, its index, the statistics, and" \
            "the copy"
    else
        record_fail "${lbl}: wrote '${got}', expected '${exp}'"
        continue
    fi

    if
        run_samtools view -H -T "${ref_fa}" "${fil_qnm}" \
            | grep -q $'^@HD\t.*SO:queryname'
    then
        record_pass "${lbl}: the copy is sorted by queryname"
    else
        record_fail "${lbl}: the copy is not sorted by queryname"
    fi

    if cmp -s <(print_recs "${fil_out}") <(print_recs "${fil_qnm}"); then
        record_pass "${lbl}: the copy holds the final file's records"
    else
        record_fail "${lbl}: the copy's records differ from the final file's"
    fi

    n_dup_out="$(run_samtools view -c -f 1024 -T "${ref_fa}" "${fil_out}")"
    n_dup_qnm="$(run_samtools view -c -f 1024 -T "${ref_fa}" "${fil_qnm}")"

    if [[ "${n_dup_out}" -gt 0 && "${n_dup_qnm}" -eq "${n_dup_out}" ]]; then
        record_pass \
            "${lbl}: the copy carries every duplicate flag (${n_dup_out})"
    else
        record_fail \
            "${lbl}: duplicate flags, final ${n_dup_out}, copy ${n_dup_qnm}"
    fi
done


finish
