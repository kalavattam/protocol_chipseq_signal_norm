#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_align_fastqs_determinism.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="align-fastqs determinism (Bowtie 2)"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Define fixture and output paths. The tie fixtures hold reads with several
# equally good placements, so Bowtie 2's seed and output order show in them.
dir_fx="${ROOT_REPO}/tests/fixtures/align_fastqs"
ref_tie="${dir_fx}/reference/tie.fa"
idx_tie="${dir_fx}/bowtie2/tie"
fq_pe_1="${dir_fx}/fastq/pe/tie_pe_R1.atria.fastq.gz"
fq_pe_2="${dir_fx}/fastq/pe/tie_pe_R2.atria.fastq.gz"
fq_cmt_1="${dir_fx}/fastq/pe/tie_cmt_pe_R1.atria.fastq.gz"
fq_cmt_2="${dir_fx}/fastq/pe/tie_cmt_pe_R2.atria.fastq.gz"
fq_idx_1="${dir_fx}/fastq/pe/tie_idx_pe_R1.atria.fastq.gz"
fq_idx_2="${dir_fx}/fastq/pe/tie_idx_pe_R2.atria.fastq.gz"
fq_se="${dir_fx}/fastq/se/tie_se.atria.fastq.gz"
fq_cmt_se="${dir_fx}/fastq/se/tie_cmt_se.atria.fastq.gz"

tmp="${TEST_DIR_TMP}/align_fastqs_determinism"


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${tmp}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty \
    "${ref_tie}" \
    "${idx_tie}.1.bt2" \
    "${fq_pe_1}" \
    "${fq_pe_2}" \
    "${fq_cmt_1}" \
    "${fq_cmt_2}" \
    "${fq_idx_1}" \
    "${fq_idx_2}" \
    "${fq_se}" \
    "${fq_cmt_se}" || {
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
    ' _ "${ROOT_REPO}" "$@"
}


# Print the Bowtie 2 command of the helper's dry run, one argument per line,
# for the given inputs and options.
function print_bt2_cmd() {
    run_helper \
        --threads 1 \
        --aligner bowtie2 \
        --index "${idx_tie}" \
        "$@" \
        --fil_out "${tmp}/dry.bam" \
        --dry_run \
        2> /dev/null \
        | sed -n '/^bowtie2 \\$/,/^ *| samtools/p' \
        | sed -e '$d' -e 's/^ *//' -e 's/ \\$//' \
        | tr ' ' '\n'
}


# Load the Bowtie 2 command for the given helper arguments into 'arr_cmd'.
function load_bt2_cmd() {
    arr_cmd=()
    mapfile -t arr_cmd < <(print_bt2_cmd "$@")

    if (( ${#arr_cmd[@]} == 0 )); then
        record_fail "no Bowtie 2 command in the dry run of: $*"
        return 1
    fi
}


# Print the arguments with one option's value replaced, or with the option and
# its value dropped when the new value is empty.
function print_with_opt() {
    local opt="${1}"
    local val="${2}"
    local -a arr_out=()
    shift 2

    while (( $# > 0 )); do
        if [[ "${1}" == "${opt}" ]]; then
            [[ -n "${val}" ]] && arr_out+=( "${opt}" "${val}" )
            shift 2
            continue
        fi

        arr_out+=( "${1}" )
        shift 1
    done

    printf '%s\n' "${arr_out[@]}"
}


# Run a Bowtie 2 command and keep its alignment records, without the header.
function run_bt2() {
    local out="${1}"
    shift 1

    "$@" 2> /dev/null | grep -v '^@' > "${out}"
}


# Print how many records differ in flag, position, MAPQ, CIGAR, or mate fields
# between two outputs written in input order; read names are left out, since
# they differ between header forms.
function count_moved() {
    if (( $(wc -l < "${1}") != $(wc -l < "${2}") )); then
        echo "mismatch"
        return 0
    fi

    paste <(cut -f 2-9 "${1}") <(cut -f 2-9 "${2}") \
        | awk -F '\t' '{
            for (i = 1; i <= 8; i++) {
                if ($i != $(i + 8)) { n++; break }
            }
        } END { print n + 0 }'
}


# Print a file's alignment records.
function print_records() {
    run_samtools view "${1}"
}


# Print a file's '@PG' program names, in header order.
function print_pg_ids() {
    run_samtools view -H --no-PG "${1}" \
        | grep '^@PG' \
        | tr '\t' '\n' \
        | grep '^PN:' \
        || true
}


# Print a file's pipeline content lines.
function print_cnt() {
    run_samtools view -H --no-PG "${1}" | grep $'^@CO\tpipeline: ' || true
}


# label|helper input arguments
arr_calls=(
    "paired-end|--fq_1 ${fq_pe_1} --fq_2 ${fq_pe_2}"
    "paired-end, '--bt2_inherited'|--fq_1 ${fq_pe_1} --fq_2 ${fq_pe_2} --bt2_inherited"
    "single-end|--fq_1 ${fq_se}"
)


# The command carries '--reorder' and '--seed 0' in every Bowtie 2 call.
for row in "${arr_calls[@]}"; do
    IFS='|' read -r lbl arg <<< "${row}"
    read -r -a arr_args <<< "${arg}"
    load_bt2_cmd "${arr_args[@]}" || continue

    cmd=" ${arr_cmd[*]} "
    if [[ "${cmd}" == *" --reorder "* && "${cmd}" == *" --seed 0 "* ]]; then
        record_pass "${lbl}: the Bowtie 2 call passes '--reorder --seed 0'"
    else
        record_fail "${lbl}: the Bowtie 2 call lacks '--reorder' or '--seed 0'"
    fi
done


# The call's output is the same at one thread and at eight, record for record
# and in order, and '--seed 0' gives the output of no '--seed'.
for idx in "${!arr_calls[@]}"; do
    IFS='|' read -r lbl arg <<< "${arr_calls[idx]}"
    read -r -a arr_args <<< "${arg}"
    load_bt2_cmd "${arr_args[@]}" || continue

    mapfile -t arr_p1 < <(print_with_opt -p 1 "${arr_cmd[@]}")
    mapfile -t arr_p8 < <(print_with_opt -p 8 "${arr_cmd[@]}")
    mapfile -t arr_ns < <(print_with_opt --seed "" "${arr_cmd[@]}")

    run_bt2 "${tmp}/p1_${idx}.sam" "${arr_p1[@]}"
    run_bt2 "${tmp}/p8_${idx}.sam" "${arr_p8[@]}"
    run_bt2 "${tmp}/ns_${idx}.sam" "${arr_ns[@]}"

    if [[ ! -s "${tmp}/p1_${idx}.sam" ]]; then
        record_fail "${lbl}: Bowtie 2 wrote no records"
    elif cmp -s "${tmp}/p1_${idx}.sam" "${tmp}/p8_${idx}.sam"; then
        record_pass "${lbl}: '-p 8' output equals '-p 1', in order"
    else
        record_fail "${lbl}: '-p 8' output differs from '-p 1'"
    fi

    if cmp -s "${tmp}/p1_${idx}.sam" "${tmp}/ns_${idx}.sam"; then
        record_pass "${lbl}: '--seed 0' output equals no '--seed'"
    else
        record_fail "${lbl}: '--seed 0' output differs from no '--seed'"
    fi
done


# Facts about Bowtie 2 that the help states, checked on the fixture: tied
# single-end reads move with the seed and with the header; in a pair, the seed
# and header text both mates share cancel, and only a difference between the
# mates moves the pair.
load_bt2_cmd --fq_1 "${fq_se}" || true
mapfile -t arr_se_7 < <(print_with_opt --seed 7 "${arr_cmd[@]}")
run_bt2 "${tmp}/se_s7.sam" "${arr_se_7[@]}"

load_bt2_cmd --fq_1 "${fq_cmt_se}" || true
run_bt2 "${tmp}/se_cmt.sam" "${arr_cmd[@]}"

load_bt2_cmd --fq_1 "${fq_pe_1}" --fq_2 "${fq_pe_2}" || true
mapfile -t arr_pe_7 < <(print_with_opt --seed 7 "${arr_cmd[@]}")
run_bt2 "${tmp}/pe_s7.sam" "${arr_pe_7[@]}"

load_bt2_cmd --fq_1 "${fq_cmt_1}" --fq_2 "${fq_cmt_2}" || true
run_bt2 "${tmp}/pe_cmt.sam" "${arr_cmd[@]}"

load_bt2_cmd --fq_1 "${fq_idx_1}" --fq_2 "${fq_idx_2}" || true
run_bt2 "${tmp}/pe_idx.sam" "${arr_cmd[@]}"

# label|output A|output B|expectation: 'moved' or 'same'
arr_facts=(
    "single-end, '--seed 7'|p1_2.sam|se_s7.sam|moved"
    "single-end, a FASTQ comment|p1_2.sam|se_cmt.sam|moved"
    "paired-end, '--seed 7'|p1_0.sam|pe_s7.sam|same"
    "paired-end, comments differing between mates|p1_0.sam|pe_cmt.sam|moved"
    "paired-end, an index both mates share|pe_cmt.sam|pe_idx.sam|same"
)

for row in "${arr_facts[@]}"; do
    IFS='|' read -r lbl fil_a fil_b exp <<< "${row}"

    if [[ ! -s "${tmp}/${fil_a}" || ! -s "${tmp}/${fil_b}" ]]; then
        record_fail "${lbl}: '${fil_a}' or '${fil_b}' has no records"
        continue
    fi

    n="$(count_moved "${tmp}/${fil_a}" "${tmp}/${fil_b}")"

    if [[ "${n}" == "mismatch" ]]; then
        record_fail "${lbl}: the two runs wrote different record counts"
    elif [[ "${exp}" == "moved" ]] && (( n > 0 )); then
        record_pass "${lbl}: ${n} tied records moved"
    elif [[ "${exp}" == "same" ]] && (( n == 0 )); then
        record_pass "${lbl}: no record moved"
    else
        record_fail "${lbl}: ${n} records moved, expected '${exp}'"
    fi
done


# The pipeline record names the seed.
fil_rec="${tmp}/rec.bam"
if
    run_helper \
        --threads 1 \
        --aligner bowtie2 \
        --index "${idx_tie}" \
        --fq_1 "${fq_se}" \
        --fil_out "${fil_rec}" \
        > "${tmp}/rec.log" 2>&1 \
        && print_cnt "${fil_rec}" | grep -qx $'@CO\tpipeline: bt2_seed=0'
then
    record_pass "the record gives 'bt2_seed=0'"
else
    record_fail "the record lacks 'bt2_seed=0'; see '${tmp}/rec.log'"
fi


# The chain run twice, at one thread and at eight, writes the same records,
# content lines, and '@PG' program sequence, in the final and '--qname' files.
for row in "pe|--fq_1 ${fq_pe_1} --fq_2 ${fq_pe_2}" "se|--fq_1 ${fq_se}"; do
    IFS='|' read -r lay arg <<< "${row}"
    read -r -a arr_args <<< "${arg}"

    for thr in 1 8; do
        mkdir -p "${tmp}/chain_${lay}_${thr}"
        run_helper \
            --threads "${thr}" \
            --aligner bowtie2 \
            --index "${idx_tie}" \
            "${arr_args[@]}" \
            --qname \
            --fil_out "${tmp}/chain_${lay}_${thr}/x.bam" \
            > "${tmp}/chain_${lay}_${thr}/run.log" 2>&1 || true
    done

    for fil in x.bam x.qnam.bam; do
        fil_1="${tmp}/chain_${lay}_1/${fil}"
        fil_8="${tmp}/chain_${lay}_8/${fil}"

        if [[ ! -s "${fil_1}" || ! -s "${fil_8}" ]]; then
            record_fail "${lay} chain: '${fil}' is missing"
        elif [[
            "$(print_records "${fil_1}")" == "$(print_records "${fil_8}")"
            && "$(print_cnt "${fil_1}")" == "$(print_cnt "${fil_8}")"
            && "$(print_pg_ids "${fil_1}")" == "$(print_pg_ids "${fil_8}")"
        ]]; then
            record_pass \
                "${lay} chain, '${fil}': one and eight threads give the same" \
                "records, content lines, and programs"
        else
            record_fail \
                "${lay} chain, '${fil}': one and eight threads differ"
        fi
    done
done


finish
