#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_align_fastqs_bt2_fragments.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="align-fastqs Bowtie 2 fragments"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Define fixture and output paths. Pairs cut from the repeat reference's unique
# blocks show which fragments each Bowtie 2 paired-end call keeps.
dir_fx="${ROOT_REPO}/tests/fixtures/align_fastqs"
ref_rpt="${dir_fx}/reference/repeat.fa"
idx_rpt="${dir_fx}/bowtie2/repeat"

tmp="${TEST_DIR_TMP}/align_fastqs_bt2_fragments"
fq_1="${tmp}/frg_R1.fastq"
fq_2="${tmp}/frg_R2.fastq"
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
    "${idx_rpt}.1.bt2" || {
    finish
    exit $?
}


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


# The unique blocks are 0-399, 550-949, 1100-1499, and 1650-2049. 'mid' pairs
# are ordinary 300 bp fragments; 'shrt' pairs are 80 bp, so their mates
# overlap; 'long' pairs are 800 bp, longer than Bowtie 2's default '-X' of 500.
seq_ref="$(grep -v '^>' "${ref_rpt}" | tr -d '\n')"

if [[ "${#seq_ref}" -ne 2050 ]]; then
    record_fail "the repeat reference is not the expected 2050 bases"
    finish
    exit $?
fi

print_pair mid_a 50 300
print_pair mid_b 1700 300
print_pair shrt_a 600 80
print_pair shrt_b 1200 80
print_pair long_a 100 800
print_pair long_b 1150 800


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


# label|arguments|pairs kept, in name order
rows=(
    "default||long_a long_b mid_a mid_b shrt_a shrt_b"
    "'--bt2_X 500'|--bt2_X 500|mid_a mid_b shrt_a shrt_b"
    "'--bt2_inherited'|--bt2_inherited|mid_a mid_b"
)

for idx in "${!rows[@]}"; do
    IFS='|' read -r lbl arg exp <<< "${rows[idx]}"
    read -r -a args <<< "${arg}"
    fil_out="${tmp}/out_${idx}.bam"

    rc=0
    out="$(
        run_helper \
            --aligner bowtie2 \
            --index "${idx_rpt}" \
            --fq_1 "${fq_1}" \
            --fq_2 "${fq_2}" \
            "${args[@]}" \
            --fil_out "${fil_out}"
    )" || rc=$?

    if [[ "${rc}" -ne 0 ]]; then
        record_fail "${lbl}: exit ${rc}"
        printf '%s\n' "${out}" >&2
        continue
    fi

    got="$(
        run_samtools view "${fil_out}" | cut -f 1 | LC_ALL=C sort -u \
            | paste -sd ' ' -
    )"

    if [[ "${got}" == "${exp}" ]]; then
        record_pass "${lbl} keeps: ${exp}"
    else
        record_fail "${lbl} kept '${got}', expected '${exp}'"
    fi
done


finish
