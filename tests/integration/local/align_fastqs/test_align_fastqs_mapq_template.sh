#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_align_fastqs_mapq_template.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="align-fastqs MAPQ by template"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Define fixture and output paths. BWA and BWA-MEM2 can give two mates
# different MAPQ values and write supplementary records, so '--mapq' must keep
# or drop each template whole; the 'repeat' fixtures exercise this.
dir_fx="${ROOT_REPO}/tests/fixtures/align_fastqs"
idx_bt2="${dir_fx}/bowtie2/repeat"
idx_bwa="${dir_fx}/bwa/repeat.fa"
idx_bm2="${dir_fx}/bwa-mem2/repeat.fa"
fq_se="${dir_fx}/fastq/se/repeat_se.atria.fastq.gz"
fq_1="${dir_fx}/fastq/pe/repeat_pe_R1.atria.fastq.gz"
fq_2="${dir_fx}/fastq/pe/repeat_pe_R2.atria.fastq.gz"
sam_sec_pe="${dir_fx}/sam/repeat_secondary_pe.sam"
sam_sec_se="${dir_fx}/sam/repeat_secondary_se.sam"

tmp="${TEST_DIR_TMP}/align_fastqs_mapq_template"
dir_log="${TEST_DIR_LOG}/align_fastqs"


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${tmp}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty \
    "${idx_bt2}.1.bt2" \
    "${idx_bwa}.bwt" \
    "${idx_bm2}.0123" \
    "${fq_se}" \
    "${fq_1}" \
    "${fq_2}" \
    "${sam_sec_pe}" \
    "${sam_sec_se}" || {
    finish
    exit $?
}


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


# Print each read name with its record count, as 'name:count' in name order.
function print_names() {
    run_samtools view "${1}" | cut -f 1 | LC_ALL=C sort | uniq -c \
        | awk '{ printf "%s%s:%s", (NR > 1 ? " " : ""), $2, $1 }'
}


# Print, for each read name, its primary MAPQs and then its supplementary and
# secondary MAPQs ('p' and 's'), as 'name:p0,p60,s0' in name order.
function print_mapqs() {
    run_samtools view "${1}" | awk -F '\t' '{
            c = (int($2 / 256) % 2 || int($2 / 2048) % 2) ? "s" : "p"
            print $1, c $5
        }' | LC_ALL=C sort -k 1,1 -k 2,2 \
        | awk '
            $1 != n {
                if (n != "") printf "%s ", s
                n = $1
                s = n ":" $2
                next
            }
            { s = s "," $2 }
            END { if (n != "") print s }
        '
}


# Align and report the output path, or record a failure and print nothing.
function align_case() {
    local idx="${1}" typ="${2}" arg="${3}"
    local fil_out="${tmp}/out_${idx}.bam" rc=0 out
    local -a args inputs

    read -r -a args <<< "${arg}"
    if [[ "${typ}" == "pe" ]]; then
        inputs=( --fq_1 "${fq_1}" --fq_2 "${fq_2}" )
    else
        inputs=( --fq_1 "${fq_se}" )
    fi

    out="$(
        run_helper "${args[@]}" "${inputs[@]}" --fil_out "${fil_out}"
    )" || rc=$?

    if [[ "${rc}" -ne 0 || ! -s "${fil_out}" ]]; then
        printf '%s\n' "${out}" > "${dir_log}/mapq_template_${idx}.log"
        return 1
    fi
    printf '%s\n' "${fil_out}"
}


# Check that the fixture still exercises the cases: unfiltered, each aligner
# writes these MAPQs. If an aligner update changes them, the rows below no
# longer test what they claim.
bwa_mem="--aligner bwa --index ${idx_bwa}"
bwa_aln="--aligner bwa --bwa_alg aln --index ${idx_bwa}"
bm2="--aligner bwa-mem2 --index ${idx_bm2}"
bt2="--aligner bowtie2 --index ${idx_bt2}"
bt2_loc="${bt2} --bt2_mode local"

mq_mem_pe="chim:p0,p60,s0 ctrl:p60,p60 half:p0,p60 rchm:p60,p60,s0"
mq_mem_se="chim:p0,s0 ctrl:p60 rchm:p60,s0"

# label|input|arguments|expected MAPQs
pre=(
    "bwa mem|pe|${bwa_mem}|${mq_mem_pe}"
    "bwa-mem2|pe|${bm2}|${mq_mem_pe}"
    "bwa aln|pe|${bwa_aln}|ctrl:p60,p60 half:p29,p37"
    "bwa mem|se|${bwa_mem}|${mq_mem_se}"
    "bwa-mem2|se|${bm2}|${mq_mem_se}"
)

for idx in "${!pre[@]}"; do
    IFS='|' read -r lbl typ arg exp <<< "${pre[idx]}"
    lbl="${lbl} ${typ^^} unfiltered"

    if ! fil_out="$(align_case "pre_${idx}" "${typ}" "${arg} --mapq 0")"; then
        record_fail "${lbl} failed"
        continue
    fi

    got="$(print_mapqs "${fil_out}")"
    if [[ "${got}" == "${exp}" ]]; then
        record_pass "${lbl} writes the MAPQs the cases need"
    else
        record_fail "${lbl} MAPQs changed (expected '${exp}'; got '${got}')"
    fi
done


# Bowtie 2 keeps the per-record filter, which is whole-template only while both
# mates get one MAPQ and no supplementary or secondary records are written.
# Assert both, including on the chimeric reads in local mode.
guard=(
    "bowtie2 end-to-end|pe|${bt2}"
    "bowtie2 local|pe|${bt2_loc}"
    "bowtie2 local|se|${bt2_loc}"
)

for idx in "${!guard[@]}"; do
    IFS='|' read -r lbl typ arg <<< "${guard[idx]}"
    lbl="${lbl} ${typ^^}"

    if ! fil_out="$(align_case "bt2_${idx}" "${typ}" "${arg} --mapq 0")"; then
        record_fail "${lbl} failed"
        continue
    fi

    # Count records that are supplementary or secondary, or that are primary
    # with a MAPQ unlike their mate's.
    n_bad="$(
        print_mapqs "${fil_out}" | tr ' ' '\n' | awk -F '[:,]' '
            /,s|:s/ { n++; next }
            NF == 3 && $2 != $3 { n++ }
            END { print n + 0 }
        '
    )"

    if [[ "${n_bad}" -eq 0 ]]; then
        record_pass "${lbl} gives mates one MAPQ and writes no split records"
    else
        record_fail \
            "${lbl} wrote ${n_bad} template(s) a per-record '-q' would" \
            "split ($(print_mapqs "${fil_out}"))"
    fi
done


# label|input|arguments|expected 'name:records'
rows=(
    "bwa mem|pe|${bwa_mem} --mapq 0|chim:3 ctrl:2 half:2 rchm:3"
    "bwa mem|pe|${bwa_mem} --mapq 1|ctrl:2 rchm:3"
    "bwa mem|pe|${bwa_mem} --mapq 30|ctrl:2 rchm:3"
    "bwa aln|pe|${bwa_aln} --mapq 1|ctrl:2 half:2"
    "bwa aln|pe|${bwa_aln} --mapq 30|ctrl:2"
    "bwa-mem2|pe|${bm2} --mapq 1|ctrl:2 rchm:3"
    "bwa mem|se|${bwa_mem} --mapq 1|ctrl:1 rchm:2"
    "bwa-mem2|se|${bm2} --mapq 1|ctrl:1 rchm:2"
    "bowtie2|pe|${bt2_loc} --mapq 20|chim:2 ctrl:2 half:2"
    "bowtie2|se|${bt2_loc} --mapq 20|ctrl:1"
)

for idx in "${!rows[@]}"; do
    IFS='|' read -r lbl typ arg exp <<< "${rows[idx]}"
    mapq="${arg##*--mapq }"
    lbl="${lbl} ${typ^^} with '--mapq ${mapq}'"

    if ! fil_out="$(align_case "${idx}" "${typ}" "${arg}")"; then
        record_fail "${lbl} failed"
        continue
    fi

    got="$(print_names "${fil_out}")"
    if [[ "${got}" == "${exp}" ]]; then
        record_pass "${lbl} keeps or drops each template whole"
    else
        record_fail "${lbl} kept '${got}', expected '${exp}'"
    fi
done


# The aligners as called write no secondary records, so check the template
# filter on its own with fixtures that have them: a template is kept or dropped
# by its primary records, never by a secondary record's MAPQ.
function run_filter() {
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
        _filter_bam_mapq_template "$@"
    ' _ "${ROOT_REPO}" "$@" 2>&1
}


# label|SAM file|primary records per template|expected 'name:records'
sec=(
    "paired-end|${sam_sec_pe}|2|sec_keep:3"
    "single-end|${sam_sec_se}|1|sec_keep:2"
)

for idx in "${!sec[@]}"; do
    IFS='|' read -r lbl sam n_pri exp <<< "${sec[idx]}"
    lbl="template filter, ${lbl} with secondary records, MAPQ 30"
    fil_out="${tmp}/sec_${idx}.bam"

    rc=0
    out="$(run_filter "${sam}" "${fil_out}" 30 "${n_pri}" 1)" || rc=$?

    if [[ "${rc}" -ne 0 || ! -s "${fil_out}" ]]; then
        printf '%s\n' "${out}" > "${dir_log}/mapq_template_sec_${idx}.log"
        record_fail "${lbl} failed (exit ${rc})"
        continue
    fi

    got="$(print_names "${fil_out}")"
    if [[ "${got}" == "${exp}" ]]; then
        record_pass "${lbl} follows the primary records"
    else
        record_fail "${lbl} kept '${got}', expected '${exp}'"
    fi
done

# Step 2 reports a failure of the template filter in the caller's name. A
# directory where the filter writes its name list makes the filter fail, even
# for root; the SAM fixture stands in for a single-end BAM work file.
dir_fail="${tmp}/filter_fail"
mkdir -p "${dir_fail}/x.mapq.names.txt"
cp "${sam_sec_se}" "${dir_fail}/x.bam"

rc=0
# shellcheck disable=SC2016  # Expand in the child shell, not this one.
out="$(
    "${TEST_BASH}" -c '
        # shellcheck source=lib/bash/core/source_helpers.sh
        source "${1}/lib/bash/core/source_helpers.sh"

        source_helpers "${1}/lib/bash" \
            core/check_args \
            core/check_inputs \
            core/format_outputs \
            workflows/align_fastqs
        shift 1
        _sort_qname_bam "$@"
    ' _ "${ROOT_REPO}" align_fastqs false 1 bwa 30 "" \
        "${dir_fail}/x.bam" "${dir_fail}/x.qnam.bam" 2>&1
)" || rc=$?
n_own="$(grep -c '^error([^)]*::align_fastqs)' <<< "${out}" || true)"
n_err="$(grep -c '^error(' <<< "${out}" || true)"

if [[
    "${rc}" -ne 0
    && "${n_own}" -eq 2
    && "${n_err}" -eq 2
    && "${out}" == *"failed to list read names"*
    && "${out}" == *"Step #2: failed to filter"*
]]; then
    record_pass "template filter failure is reported in the caller's name"
else
    printf '%s\n' "${out}" > "${dir_log}/mapq_template_filter_fail.log"
    record_fail \
        "template filter failure was not reported in the caller's name" \
        "(exit ${rc}; ${n_own} of ${n_err} errors by name)"
fi


finish
