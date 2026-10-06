#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_align_fastqs_unaligned.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="align-fastqs unaligned reads"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Define fixture and output paths. BWA and BWA-MEM2 write unaligned reads, so
# the helper must drop them, and with paired-end data the mates of unaligned
# reads too, at any '--mapq'.
dir_fx="${ROOT_REPO}/tests/fixtures/align_fastqs"
in_se="${dir_fx}/fastq/se/tiny_se.atria.fastq.gz"
in_pe_1="${dir_fx}/fastq/pe/tiny_pe_R1.atria.fastq.gz"
in_pe_2="${dir_fx}/fastq/pe/tiny_pe_R2.atria.fastq.gz"
idx_bt2="${dir_fx}/bowtie2/tiny"
idx_bwa="${dir_fx}/bwa/tiny.fa"
idx_bm2="${dir_fx}/bwa-mem2/tiny.fa"

tmp="${TEST_DIR_TMP}/align_fastqs_unaligned"
dir_log="${TEST_DIR_LOG}/align_fastqs"
fq_se="${tmp}/se.fastq"
fq_1="${tmp}/pe_R1.fastq"
fq_2="${tmp}/pe_R2.fastq"

# Sequences absent from the tiny reference.
seq_jnk_1="GGGGCCCCAAAATTTTGGGGCCCCAAAATT"
seq_jnk_2="CACACACAGTGTGTGTCACACACAGTGTGT"


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${tmp}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty \
    "${in_se}" \
    "${in_pe_1}" \
    "${in_pe_2}" \
    "${idx_bt2}.1.bt2" \
    "${idx_bwa}" \
    "${idx_bm2}" || {
    finish
    exit $?
}


# Print one FASTQ record.
function print_fq() {
    printf '@%s\n%s\n+\n%s\n' "${1}" "${2}" "${2//?/I}"
}

# Add an unaligned read to the single-end input. Add an unaligned pair and a
# pair whose second mate is unaligned to the paired-end input.
seq_se="$(gzip -cd "${in_se}" | sed -n 2p)"
seq_r1="$(gzip -cd "${in_pe_1}" | sed -n 2p)"

{ gzip -cd "${in_se}"; print_fq jnk "${seq_jnk_1}"; } > "${fq_se}"
{
    gzip -cd "${in_pe_1}"
    print_fq jnk "${seq_jnk_1}"
    print_fq hlf "${seq_r1}"
} > "${fq_1}"
{
    gzip -cd "${in_pe_2}"
    print_fq jnk "${seq_jnk_2}"
    print_fq hlf "${seq_jnk_2}"
} > "${fq_2}"

if [[ -z "${seq_se}" || -z "${seq_r1}" ]]; then
    record_fail "could not read the fixture FASTQ sequences"
    finish
    exit $?
fi


# Run the helper.
function run_helper() {
    # shellcheck disable=SC2016  # Expand in the child shell, not this one.
    "${TEST_BASH}" -c '
        # shellcheck disable=SC1090
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


# label|input ('se' or 'pe')|arguments|expected primary records
rows=(
    "bwa mem|pe|--aligner bwa --index ${idx_bwa} --mapq 0|2"
    "bwa mem|pe|--aligner bwa --index ${idx_bwa} --mapq 1|2"
    "bwa aln|pe|--aligner bwa --bwa_alg aln --index ${idx_bwa} --mapq 0|2"
    "bwa-mem2|pe|--aligner bwa-mem2 --index ${idx_bm2} --mapq 0|2"
    "bowtie2|pe|--aligner bowtie2 --index ${idx_bt2} --mapq 0|2"
    "bwa mem|se|--aligner bwa --index ${idx_bwa} --mapq 0|1"
    "bwa-mem2|se|--aligner bwa-mem2 --index ${idx_bm2} --mapq 0|1"
)

for idx in "${!rows[@]}"; do
    IFS='|' read -r lbl typ arg n_exp <<< "${rows[idx]}"
    read -r -a args <<< "${arg}"
    fil_out="${tmp}/out_${idx}.bam"
    mapq="${arg##*--mapq }"
    lbl="${lbl} ${typ^^} with '--mapq ${mapq}'"

    if [[ "${typ}" == "pe" ]]; then
        inputs=( --fq_1 "${fq_1}" --fq_2 "${fq_2}" )
    else
        inputs=( --fq_1 "${fq_se}" )
    fi

    rc=0
    out="$(
        run_helper "${args[@]}" "${inputs[@]}" --fil_out "${fil_out}"
    )" || rc=$?

    if [[ "${rc}" -ne 0 || ! -s "${fil_out}" ]]; then
        printf '%s\n' "${out}" > "${dir_log}/unaligned_${idx}.log"
        record_fail "${lbl} failed (exit ${rc})"
        continue
    fi

    # Primary records only; no read is unaligned, and only the fixture's own
    # read or pair is kept, so no mate is left alone.
    names="$(
        run_samtools view -F 0x900 "${fil_out}" | cut -f 1 | sort | uniq -c \
            | awk '{ printf "%s:%s ", $2, $1 }'
    )"
    n_unal="$(run_samtools view -c -f 4 "${fil_out}")"

    if [[ "${n_unal}" -eq 0 && "${names}" == "tiny_${typ}_"*":${n_exp} " ]]
    then
        record_pass "${lbl} keeps no unaligned read and no lone mate"
    else
        record_fail \
            "${lbl} kept unaligned reads or lone mates (unaligned" \
            "${n_unal}; records ${names})"
    fi
done


finish
