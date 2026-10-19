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
# reads too, at any '--mapq'. The 'bridge' fixtures hold a pair whose second
# mate runs off the end of one contig into the next.
dir_fx="${ROOT_REPO}/tests/fixtures/align_fastqs"
in_se="${dir_fx}/fastq/se/tiny_se.atria.fastq.gz"
in_pe_1="${dir_fx}/fastq/pe/tiny_pe_R1.atria.fastq.gz"
in_pe_2="${dir_fx}/fastq/pe/tiny_pe_R2.atria.fastq.gz"
idx_bt2="${dir_fx}/bowtie2/tiny"
idx_bwa="${dir_fx}/bwa/tiny.fa"
idx_bm2="${dir_fx}/bwa-mem2/tiny.fa"
brg_1="${dir_fx}/fastq/pe/bridge_pe_R1.atria.fastq.gz"
brg_2="${dir_fx}/fastq/pe/bridge_pe_R2.atria.fastq.gz"
idx_brg_bt2="${dir_fx}/bowtie2/bridge"
idx_brg_bwa="${dir_fx}/bwa/bridge.fa"
idx_brg_bm2="${dir_fx}/bwa-mem2/bridge.fa"

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
    "${idx_bm2}" \
    "${brg_1}" \
    "${brg_2}" \
    "${idx_brg_bt2}.1.bt2" \
    "${idx_brg_bwa}.bwt" \
    "${idx_brg_bm2}.0123" || {
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


# BWA ALN can flag one mate unaligned and leave its partner's mate-unaligned
# flag clear, so a per-record '-F 12' keeps the partner. Check that the raw
# 'bwa sampe' output of the bridge pair still has that shape: one 'brdg' record
# carries 0x4 and the other lacks 0x8. Otherwise the rows below test nothing.
rc=0
# shellcheck disable=SC2016  # Expand in the child shell, not this one.
flg_raw="$(
    "${TEST_BASH}" -c '
        bwa aln "${1}" "${2}" > "${4}/brg_1.sai" \
            && bwa aln "${1}" "${3}" > "${4}/brg_2.sai" \
            && bwa sampe \
                "${1}" "${4}/brg_1.sai" "${4}/brg_2.sai" "${2}" "${3}"
    ' _ "${idx_brg_bwa}" "${brg_1}" "${brg_2}" "${tmp}" \
        2> "${dir_log}/unaligned_bridge_raw.log" \
        | awk -F '\t' '$1 == "brdg" { print $2 }' \
        | tr '\n' ' '
)" || rc=$?
read -r -a arr_flg_raw <<< "${flg_raw}"

if [[
    "${rc}" -eq 0
    && "${#arr_flg_raw[@]}" -eq 2
    && $(( (arr_flg_raw[0] ^ arr_flg_raw[1]) & 4 )) -ne 0
    && $(( (arr_flg_raw[0] | arr_flg_raw[1]) & 8 )) -eq 0
]]; then
    record_pass "bwa aln bridge pair flags one mate unaligned, not its partner"
else
    record_fail \
        "bwa aln bridge pair no longer has an unaligned mate with an" \
        "unflagged partner (exit ${rc}; flags '${flg_raw}')"
fi

# Each aligner on the bridge pair keeps only whole templates, with or without
# '--req_flg'. BWA ALN must drop 'brdg' whole. Bowtie 2 end-to-end leaves
# 'brdg' R2 unaligned and flags its partner, so it drops 'brdg' whole too. The
# other aligners align R2 in part and keep 'brdg'.
bwa_aln="--aligner bwa --bwa_alg aln --index ${idx_brg_bwa}"
bwa_mem="--aligner bwa --index ${idx_brg_bwa}"
bm2="--aligner bwa-mem2 --index ${idx_brg_bm2}"
bt2="--aligner bowtie2 --index ${idx_brg_bt2}"

# label|arguments|expected primary records, as 'name:count'
brg=(
    "bwa aln|${bwa_aln} --mapq 0|ctrl:2"
    "bwa aln|${bwa_aln} --req_flg --mapq 0|ctrl:2"
    "bwa aln|${bwa_aln} --mapq 1|ctrl:2"
    "bwa aln|${bwa_aln} --req_flg --mapq 1|ctrl:2"
    "bwa mem|${bwa_mem} --mapq 0|brdg:2 ctrl:2"
    "bwa-mem2|${bm2} --mapq 0|brdg:2 ctrl:2"
    "bowtie2 end-to-end|${bt2} --mapq 0|ctrl:2"
    "bowtie2 local|${bt2} --bt2_mode local --mapq 0|brdg:2 ctrl:2"
)

for idx in "${!brg[@]}"; do
    IFS='|' read -r lbl arg exp <<< "${brg[idx]}"
    read -r -a args <<< "${arg}"
    fil_out="${tmp}/out_brg_${idx}.bam"
    mapq="${arg##*--mapq }"
    lbl="${lbl} PE bridge pair with '--mapq ${mapq}'"
    if [[ "${arg}" == *--req_flg* ]]; then lbl+=" and '--req_flg'"; fi

    rc=0
    out="$(
        run_helper "${args[@]}" \
            --fq_1 "${brg_1}" \
            --fq_2 "${brg_2}" \
            --fil_out "${fil_out}"
    )" || rc=$?

    if [[ "${rc}" -ne 0 || ! -s "${fil_out}" ]]; then
        printf '%s\n' "${out}" > "${dir_log}/unaligned_brg_${idx}.log"
        record_fail "${lbl} failed (exit ${rc})"
        continue
    fi

    names="$(
        run_samtools view -F 0x900 "${fil_out}" | cut -f 1 | sort | uniq -c \
            | awk '{ printf "%s%s:%s", (NR > 1 ? " " : ""), $2, $1 }'
    )"
    n_unal="$(run_samtools view -c -f 4 "${fil_out}")"

    if [[ "${n_unal}" -eq 0 && "${names}" == "${exp}" ]]; then
        record_pass "${lbl} keeps no unaligned read and no lone mate"
    else
        record_fail \
            "${lbl} kept unaligned reads or lone mates (unaligned" \
            "${n_unal}; records '${names}', expected '${exp}')"
    fi
done


finish
