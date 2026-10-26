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

TEST_NAME="align-fastqs determinism"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Define fixture and output paths. The tie fixtures hold reads with several
# equally good placements, so they show the effects of Bowtie 2's seed, its
# output order, and BWA's batch size.
dir_fx="${ROOT_REPO}/tests/fixtures/align_fastqs"
ref_tie="${dir_fx}/reference/tie.fa"
idx_tie="${dir_fx}/bowtie2/tie"
idx_bwa="${dir_fx}/bwa/tie.fa"
idx_bm2="${dir_fx}/bwa-mem2/tie.fa"
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
    "${idx_bwa}.bwt" \
    "${idx_bm2}.0123" \
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


# Print the Step 1 aligner command from the helper's dry run with one argument
# per line. The first argument is the line that starts the command ('bowtie2',
# 'bwa mem', or 'bwa-mem2 mem'); the rest are the helper's arguments.
function print_aln_cmd() {
    local head="${1}"
    shift 1

    run_helper \
        --threads 1 \
        "$@" \
        --fil_out "${tmp}/dry.bam" \
        --dry_run \
        2> /dev/null \
        | sed -n "/^${head} \\\\\$/,/^ *| samtools/p" \
        | sed -e '$d' -e 's/^ *//' -e 's/ \\$//' \
        | tr ' ' '\n'
}


# Print the Bowtie 2 command of the helper's dry run for the given inputs and
# options.
function print_bt2_cmd() {
    print_aln_cmd bowtie2 --aligner bowtie2 --index "${idx_tie}" "$@"
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
# between two outputs written in input order; read names are left out because
# they differ between header forms.
function count_moved() {
    if (( $(wc -l < "${1}") != $(wc -l < "${2}") )); then
        echo "mismatch"
        return 0
    fi

    paste <(cut -f 2-9 "${1}") <(cut -f 2-9 "${2}") \
        | awk -F '\t' '
            {
                for (i = 1; i <= 8; i++) {
                    if ($i != $(i + 8)) { n++; break }
                }
            } END {
                print n + 0
            }
        '
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


# label|helper input arguments. Rows below read the first and third by index
# ('p1_0.sam', 'p1_2.sam'), so new rows go at the end.
arg_pe="--fq_1 ${fq_pe_1} --fq_2 ${fq_pe_2}"
arr_calls=(
    "paired-end|${arg_pe}"
    "paired-end, '--bt2_inherited'|${arg_pe} --bt2_inherited"
    "single-end|--fq_1 ${fq_se}"
    "paired-end, '--bt2_relaxed'|${arg_pe} --bt2_relaxed"
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


# At one thread and at eight, the call writes the same records in the same
# order, and '--seed 0' gives the same output as no '--seed'.
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


# Check the help's facts about Bowtie 2 on the fixture. Tied single-end reads
# move with the seed and with the header. In a pair, the seed and any header
# text both mates share cancel out, and only a difference between the mates
# moves the pair.
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


# BWA-MEM and BWA-MEM2 pass a fixed '-K' for either layout.
arr_mem=(
    "BWA-MEM|bwa mem|--aligner bwa --index ${idx_bwa}"
    "BWA-MEM2|bwa-mem2 mem|--aligner bwa-mem2 --index ${idx_bm2}"
)

for row in "${arr_mem[@]}"; do
    IFS='|' read -r lbl head arg <<< "${row}"
    read -r -a arr_args <<< "${arg}"

    for lay in pe se; do
        if [[ "${lay}" == "pe" ]]; then
            arr_in=( --fq_1 "${fq_pe_1}" --fq_2 "${fq_pe_2}" )
        else
            arr_in=( --fq_1 "${fq_se}" )
        fi

        cmd=" $(
            print_aln_cmd "${head}" "${arr_args[@]}" "${arr_in[@]}" \
                | paste -sd ' ' -
        ) "
        if [[ "${cmd}" == *" -K 100000000 "* ]]; then
            record_pass "${lbl} ${lay}: the call passes '-K 100000000'"
        else
            record_fail "${lbl} ${lay}: the call lacks '-K 100000000'"
        fi
    done
done


# Why '-K' is fixed: BWA estimates paired-end insert sizes per batch, so a
# different batch size moves records. At a fixed '-K', the thread count moves
# no records. This part would pass without '-K' too, because the fixture is
# smaller than one default batch; the command check above is the one that fails
# without it.
for idx in "${!arr_mem[@]}"; do
    IFS='|' read -r lbl head arg <<< "${arr_mem[idx]}"
    read -r -a arr_args <<< "${arg}"
    arr_cmd=()
    mapfile -t arr_cmd < <(
        print_aln_cmd "${head}" "${arr_args[@]}" \
            --fq_1 "${fq_pe_1}" --fq_2 "${fq_pe_2}"
    )

    if (( ${#arr_cmd[@]} == 0 )); then
        record_fail "${lbl}: no command in the dry run"
        continue
    fi

    mapfile -t arr_k10 < <(print_with_opt -K 10000 "${arr_cmd[@]}")
    mapfile -t arr_t4 < <(print_with_opt -t 4 "${arr_cmd[@]}")

    run_bt2 "${tmp}/mem_${idx}.sam" "${arr_cmd[@]}"
    run_bt2 "${tmp}/mem_${idx}_k10.sam" "${arr_k10[@]}"
    run_bt2 "${tmp}/mem_${idx}_t4.sam" "${arr_t4[@]}"

    if [[ ! -s "${tmp}/mem_${idx}.sam" ]]; then
        record_fail "${lbl}: the aligner wrote no records"
        continue
    fi

    n="$(count_moved "${tmp}/mem_${idx}.sam" "${tmp}/mem_${idx}_k10.sam")"
    if [[ "${n}" != "mismatch" ]] && (( n > 0 )); then
        record_pass "${lbl}: '-K 10000' moves ${n} records"
    else
        record_fail "${lbl}: '-K 10000' moved '${n}' records, expected some"
    fi

    if cmp -s "${tmp}/mem_${idx}.sam" "${tmp}/mem_${idx}_t4.sam"; then
        record_pass "${lbl}: '-t 4' output equals '-t 1'"
    else
        record_fail "${lbl}: '-t 4' output differs from '-t 1'"
    fi
done


# The pipeline record names Bowtie 2's seed and BWA-MEM's batch size.
# BWA-backtrack has no '-K', so its record has no batch line.
arr_rec=(
    "Bowtie 2|--aligner bowtie2 --index ${idx_tie}|bt2_seed=0|yes"
    "BWA-MEM|--aligner bwa --index ${idx_bwa}|bwa_K=100000000|yes"
    "BWA-MEM2|--aligner bwa-mem2 --index ${idx_bm2}|bwa_K=100000000|yes"
    "BWA-backtrack|--aligner bwa --bwa_alg aln --index ${idx_bwa}|bwa_K=|no"
)

for idx in "${!arr_rec[@]}"; do
    IFS='|' read -r lbl arg lin exp <<< "${arr_rec[idx]}"
    read -r -a arr_args <<< "${arg}"
    fil_rec="${tmp}/rec_${idx}.bam"

    if ! \
        run_helper \
            --threads 1 \
            "${arr_args[@]}" \
            --fq_1 "${fq_pe_1}" \
            --fq_2 "${fq_pe_2}" \
            --fil_out "${fil_rec}" \
            > "${tmp}/rec_${idx}.log" \
            2>&1
    then
        record_fail "${lbl}: the run failed; see '${tmp}/rec_${idx}.log'"
        continue
    fi

    got=no
    if print_cnt "${fil_rec}" | grep -q $'^@CO\tpipeline: '"${lin}"; then
        got=yes
    fi

    if [[ "${got}" == "${exp}" ]]; then
        record_pass "${lbl}: the record's '${lin}' line is present: ${got}"
    else
        record_fail "${lbl}: the record's '${lin}' line is present: ${got}"
    fi
done


# For every aligner and layout, running the chain at one thread and again at
# eight writes the same records, content lines, and '@PG' program sequence to
# both the final and '--qname' files.
arr_chain=(
    "bt2|--aligner bowtie2 --index ${idx_tie}"
    "mem|--aligner bwa --index ${idx_bwa}"
    "aln|--aligner bwa --bwa_alg aln --index ${idx_bwa}"
    "bm2|--aligner bwa-mem2 --index ${idx_bm2}"
)

for row in "${arr_chain[@]}"; do
    IFS='|' read -r aln arg <<< "${row}"
    read -r -a arr_args <<< "${arg}"

    for lay in pe se; do
        if [[ "${lay}" == "pe" ]]; then
            arr_in=( --fq_1 "${fq_pe_1}" --fq_2 "${fq_pe_2}" )
        else
            arr_in=( --fq_1 "${fq_se}" )
        fi

        dir_run="${tmp}/chain_${aln}_${lay}"
        for thr in 1 8; do
            mkdir -p "${dir_run}_${thr}"
            run_helper \
                --threads "${thr}" \
                "${arr_args[@]}" \
                "${arr_in[@]}" \
                --qname \
                --fil_out "${dir_run}_${thr}/x.bam" \
                > "${dir_run}_${thr}/run.log" \
                2>&1 \
                || true
        done

        for fil in x.bam x.qnam.bam; do
            fil_1="${dir_run}_1/${fil}"
            fil_8="${dir_run}_8/${fil}"
            lbl="${aln} ${lay} chain, '${fil}'"

            if [[ ! -s "${fil_1}" || ! -s "${fil_8}" ]]; then
                record_fail "${lbl}: missing"
            elif [[
                "$(print_records "${fil_1}")" == "$(print_records "${fil_8}")"
                    && "$(print_cnt "${fil_1}")" == "$(print_cnt "${fil_8}")"
                    && "$(print_pg_ids "${fil_1}")" \
                        == "$(print_pg_ids "${fil_8}")"
            ]]; then
                record_pass \
                    "${lbl}: one and eight threads give the same records," \
                    "content lines, and programs"
            else
                record_fail "${lbl}: one and eight threads differ"
            fi
        done
    done
done


finish
