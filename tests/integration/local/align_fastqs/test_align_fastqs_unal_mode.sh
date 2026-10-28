#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_align_fastqs_unal_mode.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="align-fastqs unaligned-read mode"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Define fixture and output paths.
#
# '--unal_mode remove' drops unaligned reads and their mates, improper pairs
# under '--req_flg', and low-MAPQ templates. '--unal_mode separate' writes each
# to its own diagnostic file and leaves the alignment file as it is, so every
# record lands in one file.
#
# Every aligner call runs on four fixtures, with and without '--req_flg':
# - 'pairs': one pair per pair class.
# - 'pairs_fill': the pair classes among 300 concordant pairs, so BWA can flag
#   pairs proper.
# - 'bridge': the BWA-backtrack lone mate of #75.
# - 'dup': duplicates among the side file's aligned mates.
#
# Each call's records are checked against the aligner's own output. The 'thr'
# fixture checks that the files do not depend on the thread count.
dir_fx="${ROOT_REPO}/tests/fixtures/align_fastqs"
ref_par="${dir_fx}/reference/pairs.fa"
sam_spl="${dir_fx}/sam/split_template_pe.sam"

tmp="${TEST_DIR_TMP}/align_fastqs_unal_mode"
dir_log="${TEST_DIR_LOG}/align_fastqs"


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${tmp}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty \
    "${ref_par}" \
    "${sam_spl}" \
    "${dir_fx}/bowtie2/pairs.rev.2.bt2" \
    "${dir_fx}/bwa/pairs.fa.bwt" \
    "${dir_fx}/bwa-mem2/pairs.fa.0123" \
    "${dir_fx}/bowtie2/bridge.rev.2.bt2" \
    "${dir_fx}/bwa/bridge.fa.bwt" \
    "${dir_fx}/bwa-mem2/bridge.fa.0123" \
    "${dir_fx}/fastq/pe/pairs_pe_R1.atria.fastq.gz" \
    "${dir_fx}/fastq/pe/pairs_pe_R2.atria.fastq.gz" \
    "${dir_fx}/fastq/pe/pairs_fill_pe_R1.atria.fastq.gz" \
    "${dir_fx}/fastq/pe/pairs_fill_pe_R2.atria.fastq.gz" \
    "${dir_fx}/fastq/se/pairs_se.atria.fastq.gz" \
    "${dir_fx}/fastq/pe/dup_pe_R1.atria.fastq.gz" \
    "${dir_fx}/fastq/pe/dup_pe_R2.atria.fastq.gz" \
    "${dir_fx}/fastq/pe/bridge_pe_R1.atria.fastq.gz" \
    "${dir_fx}/fastq/pe/bridge_pe_R2.atria.fastq.gz" \
    "${dir_fx}/sam/diag/improp.sam" \
    "${dir_fx}/sam/diag/all.markdup.txt.gz" \
    "${dir_fx}/sam/diag/all_neg.markdup.txt.gz" \
    "${dir_fx}/sam/diag/all_odd.markdup.txt.gz" \
    "${dir_fx}/bwa/thr.fa.bwt" \
    "${dir_fx}/fastq/pe/thr_pe_R1.atria.fastq.gz" \
    "${dir_fx}/fastq/pe/thr_pe_R2.atria.fastq.gz" || {
    finish
    exit $?
}


# Run a function from the sourced workflows.
function run_func() {
    # shellcheck disable=SC2016  # Expand in the child shell, not this one.
    "${TEST_BASH}" -c '
        # shellcheck source=lib/bash/core/source_helpers.sh
        source "${1}/lib/bash/core/source_helpers.sh"

        source_helpers "${1}/lib/bash" \
            core/check_args \
            core/check_inputs \
            core/format_outputs \
            workflows/align_fastqs \
            workflows/process_alignments
        shift 1
        "$@"
    ' _ "${ROOT_REPO}" "$@" 2>&1
}


# Run the helper.
function run_helper() {
    run_func align_fastqs "$@"
}


# Run an aligner call as Step 1 does, but write every record it makes, so the
# helper's outputs can be checked against it. Arguments:
# - the call's label,
# - the fixture set,
# - the second FASTQ (or empty), and
# - the output SAM.
function run_raw() {
    local cal="${1}"
    local set="${2}"
    local fq_1="${3}"
    local fq_2="${4}"
    local sam="${5}"
    local idx_bt2="${dir_fx}/bowtie2/${set}"
    local idx_bwa="${dir_fx}/bwa/${set}.fa"
    local idx_bm2="${dir_fx}/bwa-mem2/${set}.fa"
    local -a arr_bt2 arr_in

    arr_bt2=( bowtie2 -p 1 --reorder --seed 0 -x "${idx_bt2}" )
    if [[ -n "${fq_2}" ]]; then
        arr_in=( -1 "${fq_1}" -2 "${fq_2}" )
    else
        arr_in=( -U "${fq_1}" )
    fi

    case "${cal}" in
        bt2)
            if [[ -n "${fq_2}" ]]; then
                arr_in=(
                    --no-mixed --no-discordant -X 1000 "${arr_in[@]}"
                )
            fi
            "${arr_bt2[@]}" --very-sensitive --phred33 "${arr_in[@]}"
            ;;

        bt2_local)
            if [[ -n "${fq_2}" ]]; then
                arr_in=(
                    --no-mixed --no-discordant -X 1000 "${arr_in[@]}"
                )
            fi
            "${arr_bt2[@]}" --very-sensitive-local --phred33 "${arr_in[@]}"
            ;;

        bt2_relaxed)
            "${arr_bt2[@]}" --very-sensitive --phred33 -I 0 -X 1000 \
                "${arr_in[@]}"
            ;;

        mem)
            # shellcheck disable=SC2086  # An empty 'fq_2' adds no word.
            bwa mem -t 1 -K 100000000 "${idx_bwa}" "${fq_1}" ${fq_2}
            ;;

        bm2)
            # shellcheck disable=SC2086  # An empty 'fq_2' adds no word.
            bwa-mem2 mem -t 1 -K 100000000 "${idx_bm2}" "${fq_1}" ${fq_2}
            ;;

        aln)
            bwa aln -t 1 "${idx_bwa}" "${fq_1}" > "${sam%.sam}.R1.sai"
            if [[ -n "${fq_2}" ]]; then
                bwa aln -t 1 "${idx_bwa}" "${fq_2}" > "${sam%.sam}.R2.sai"
                bwa sampe "${idx_bwa}" \
                    "${sam%.sam}.R1.sai" "${sam%.sam}.R2.sai" \
                    "${fq_1}" "${fq_2}"
            else
                bwa samse "${idx_bwa}" "${sam%.sam}.R1.sai" "${fq_1}"
            fi
            ;;
    esac \
        2> "${sam%.sam}.raw.log" \
        1> "${sam}"
}


# Print one key per record: its name, which read it is, whether it is secondary
# or supplementary, and where it aligned. 'fixmate' and 'markdup' change none
# of these for an aligned record. An unaligned record has no alignment, and
# 'fixmate' moves it to its mate's position, so it is keyed by name and read
# alone.
function print_keys() {
    run_samtools view "${1}" \
        | awk -F '\t' '
            function bit(x, b) { return int(x / b) % 2 * b }
            {
                f = bit($2, 64) + bit($2, 128) + bit($2, 256) + bit($2, 2048)
                if (bit($2, 4)) { $3 = "*"; $4 = 0 }
                print $1 "\t" f "\t" $3 "\t" $4
            }
        ' \
        | sort
}


# Print each '@PG' line's ID and parent, in header order.
function print_pg_ids() {
    run_samtools view -H --no-PG "${1}" \
        | awk -F '\t' '
            $1 == "@PG" {
                id = ""
                pp = ""
                for (i = 2; i <= NF; i++) {
                    if ($i ~ /^ID:/) id = $i
                    if ($i ~ /^PP:/) pp = $i
                }
                print id, pp
            }
        '
}


# Print the records sorted, whole.
function print_recs() {
    run_samtools view "${1}" | sort
}


# Print the records sorted, without the duplicate flag and the 'dt' and 'do'
# tags 'markdup' writes with it.
function print_recs_nodup() {
    run_samtools view "${1}" \
        | awk -F '\t' -v OFS='\t' '
            {
                $2 = $2 - int($2 / 1024) % 2 * 1024
                out = $1
                for (i = 2; i <= NF; i++) {
                    if ($i !~ /^d[to]:Z:/) { out = out OFS $i }
                }
                print out
            }
        ' \
        | sort
}


# Print the read names a file holds, sorted, once each.
function print_names() {
    run_samtools view "${1}" | cut -f 1 | sort -u | tr '\n' ' '
}


# Check that each template of the raw output lands whole in one place: one of
# the files given or, unless 'all' is 'yes', none (dropped by '--mapq' or
# '--req_flg'). Print nothing when it holds or the offending names.
function print_split_templates() {
    local sam_raw="${1}"
    local all="${2}"
    local idx=0
    local fil
    shift 2

    {
        print_keys "${sam_raw}" | sed 's/^/R\t/'
        for fil in "$@"; do
            print_keys "${fil}" | sed "s/^/${idx}\t/"
            idx=$(( idx + 1 ))
        done
    } \
    | awk -F '\t' -v all="${all}" '
        {
            key = $2 SUBSEP $3 SUBSEP $4 SUBSEP $5
            if ($1 == "R") { n_raw[$2]++; raw[key] = 1 }
            else {
                if (!(($2, $1) in in_fil)) { in_fil[$2, $1] = 1; n_fil[$2]++ }
                n[$2]++
                if (!(key in raw)) { bad[$2] = 1 }
                if (key in seen) { bad[$2] = 1 }
                seen[key] = 1
            }
        }
        END {
            for (nam in n_raw) {
                if (n_fil[nam] > 1) { bad[nam] = 1 }
                if (n_fil[nam] == 1 && n[nam] != n_raw[nam]) { bad[nam] = 1 }
                if (all == "yes" && n_fil[nam] + 0 == 0) { bad[nam] = 1 }
            }
            for (nam in n) { if (!(nam in n_raw)) { bad[nam] = 1 } }
            for (nam in bad) { print nam }
        }
    ' \
    | sort \
    | tr '\n' ' '
}


# Print each read name a file holds with its number of records, sorted.
function print_counts() {
    run_samtools view "${1}" 2> /dev/null \
        | cut -f 1 \
        | sort \
        | uniq -c \
        | awk '{ printf "%s:%s ", $2, $1 }'
}


# Print a file's content lines of the pipeline record, in order.
function print_rec_lines() {
    run_samtools view -H --no-PG "${1}" | grep $'^@CO\tpipeline: '
}


# Print the records sorted, without the 'MD' and 'NM' tags, which a CRAM
# decoder writes again in its own order. Arguments are for 'samtools view'.
function print_recs_nomd() {
    run_samtools view "$@" \
        | awk -F '\t' -v OFS='\t' '
            {
                out = $1
                for (i = 2; i <= NF; i++) {
                    if ($i !~ /^(MD:Z|NM:i):/) { out = out OFS $i }
                }
                print out
            }
        ' \
        | sort
}


# Calls: label|helper arguments for the fixture set '@SET@'|layouts.
i_bt2="--index ${dir_fx}/bowtie2/@SET@"
i_bwa="--index ${dir_fx}/bwa/@SET@.fa"
i_bm2="--index ${dir_fx}/bwa-mem2/@SET@.fa"
cals=(
    "bt2|--aligner bowtie2 ${i_bt2}|pe se"
    "bt2_local|--aligner bowtie2 --bt2_mode local ${i_bt2}|pe se"
    "bt2_relaxed|--aligner bowtie2 --bt2_relaxed ${i_bt2}|pe"
    "mem|--aligner bwa ${i_bwa}|pe se"
    "bm2|--aligner bwa-mem2 ${i_bm2}|pe se"
    "aln|--aligner bwa --bwa_alg aln ${i_bwa}|pe se"
)

# Inputs: label|fixture set (index)|layout|FASTQ stem.
ins=(
    "pairs|pairs|pe|pairs"
    "pairs_fill|pairs|pe|pairs_fill"
    "bridge|bridge|pe|bridge"
    "dup|pairs|pe|dup"
    "pairs_se|pairs|se|pairs"
)

# The side file's read names, by input and call, as the aligners write them:
# with the default Bowtie 2 call, '--no-mixed' writes every pair that is not
# concordant with both mates unaligned; the relaxed call and BWA keep pairs
# with both mates aligned. BWA-backtrack's 'bridge' partner joins its mate.
declare -A side_exp=(
    ["pairs bt2"]="pair_chrs pair_dove pair_far pair_ff pair_none pair_one\
 pair_rf "
    ["pairs bt2_relaxed"]="pair_none pair_one "
    ["pairs mem"]="pair_none pair_one "
    ["pairs bm2"]="pair_none pair_one "
    ["pairs_fill mem"]="pair_none pair_one "
    ["pairs_fill bm2"]="pair_none pair_one "
    ["bridge aln"]="brdg "
    ["bridge bt2"]="brdg "
    ["pairs_se bt2"]="se_none "
    ["pairs_se mem"]="se_none "
)

# The rows of '.diag.tsv.gz', as 'file:category' pairs, for each layout.
rows_pe="file:category:main:all:main:unjudged:main:duplicate_all:"
rows_pe+="unal:all:unal:both_unaligned:unal:one_unaligned:unal:unjudged:"
rows_pe+="unal:duplicate_all:improp:all:improp:diff_ref:improp:orient_ff:"
rows_pe+="improp:orient_rr:improp:orient_rf:improp:dovetailed:improp:too_far:"
rows_pe+="improp:too_close:improp:failed_est:improp:other:improp:unjudged:"
rows_pe+="improp:duplicate_all:mapq:all:mapq:unjudged:mapq:duplicate_all:"
rows_se="file:category:main:all:main:unjudged:main:duplicate_all:unal:all:"
rows_se+="unal:unjudged:unal:duplicate_all:mapq:all:mapq:unjudged:"
rows_se+="mapq:duplicate_all:"


# Run every call on every input, in both modes, with and without '--req_flg',
# and check the invariants of '--unal_mode'.
for ent in "${ins[@]}"; do
    IFS='|' read -r inp set typ stm_fq <<< "${ent}"
    if [[ "${typ}" == "pe" ]]; then
        fq_1="${dir_fx}/fastq/pe/${stm_fq}_pe_R1.atria.fastq.gz"
        fq_2="${dir_fx}/fastq/pe/${stm_fq}_pe_R2.atria.fastq.gz"
    else
        fq_1="${dir_fx}/fastq/se/${stm_fq}_se.atria.fastq.gz"
        fq_2=""
    fi

    for cal_ent in "${cals[@]}"; do
        IFS='|' read -r cal arg typs <<< "${cal_ent}"
        if [[ " ${typs} " != *" ${typ} "* ]]; then continue; fi

        arg="${arg//@SET@/${set}}"
        read -r -a args <<< "${arg}"
        if [[ "${typ}" == "pe" ]]; then
            inputs=( --fq_1 "${fq_1}" --fq_2 "${fq_2}" )
            reqs=( no yes )
        else
            inputs=( --fq_1 "${fq_1}" )
            reqs=( no )
        fi

        sam_raw="${tmp}/${inp}.${cal}.raw.sam"
        if \
            ! run_raw "${cal}" "${set}" "${fq_1}" "${fq_2}" "${sam_raw}"
        then
            record_fail "${inp} ${cal}: the raw aligner call failed"
            continue
        fi

        for req in "${reqs[@]}"; do
            lbl="${inp} ${cal}"
            args_req=()
            if [[ "${req}" == "yes" ]]; then
                lbl+=" with '--req_flg'"
                args_req=( --req_flg )
            fi

            ok=true
            for mod in remove separate; do
                stm="${tmp}/${inp}.${cal}.${req}.${mod}"
                rc=0
                out="$(
                    run_helper \
                        --threads 1 \
                        "${args[@]}" \
                        "${args_req[@]}" \
                        --unal_mode "${mod}" \
                        "${inputs[@]}" \
                        --fil_out "${stm}.bam"
                )" || rc=$?
                printf '%s\n' "${out}" > "${stm}.log"

                if [[ "${rc}" -ne 0 || ! -s "${stm}.bam" ]]; then
                    cp "${stm}.log" "${dir_log}/unal_mode_${inp}_${cal}.log"
                    record_fail \
                        "${lbl} '--unal_mode ${mod}' failed (exit ${rc})"
                    ok=false
                fi
            done

            if [[ "${ok}" != "true" ]]; then continue; fi

            stm="${tmp}/${inp}.${cal}.${req}"
            fil_side="${stm}.separate.unal.bam"
            fil_imp="${stm}.separate.improp.bam"
            fil_mpq="${stm}.separate.mapq.bam"
            arr_dia=( "${fil_side}" )
            if [[ "${typ}" == "pe" ]]; then arr_dia+=( "${fil_imp}" ); fi
            arr_dia+=( "${fil_mpq}" )

            # Only 'separate' writes the diagnostic files, each indexed, plus
            # Markdup 2's statistics and the counts. Single-end data has no
            # improper file.
            ok_fil=true
            for fil in "${arr_dia[@]}"; do
                if [[ ! -s "${fil}" || ! -s "${fil}.bai" ]]; then
                    ok_fil=false
                fi
            done
            for sfx in unal.bam improp.bam mapq.bam all.markdup.txt.gz \
                diag.tsv.gz
            do
                if [[ -e "${stm}.remove.${sfx}" ]]; then ok_fil=false; fi
            done
            if [[
                ! -s "${stm}.separate.all.markdup.txt.gz"
                    || ! -s "${stm}.separate.diag.tsv.gz"
            ]]; then
                ok_fil=false
            fi
            if [[ "${typ}" == "se" && -e "${fil_imp}" ]]; then ok_fil=false; fi

            if [[ "${ok_fil}" == "true" ]]; then
                record_pass \
                    "${lbl}: only 'separate' writes the diagnostic files"
            else
                record_fail \
                    "${lbl}: diagnostic files missing or written by 'remove'"
            fi

            # Invariant: 'remove' and 'separate' give the same alignment file,
            # duplicate flags included, and the same duplicate statistics.
            if \
                cmp -s \
                    <(print_recs "${stm}.remove.bam") \
                    <(print_recs "${stm}.separate.bam") \
                && cmp -s \
                    <(gzip -cd "${stm}.remove.markdup.txt.gz" | sed 1d) \
                    <(gzip -cd "${stm}.separate.markdup.txt.gz" | sed 1d)
            then
                record_pass \
                    "${lbl}: the alignment file and its duplicate statistics" \
                    "are the same in both modes"
            else
                record_fail \
                    "${lbl}: the alignment file or its duplicate statistics" \
                    "differ between modes"
            fi

            # Invariants: under 'separate', every raw record lands once in one
            # of the four files, with each template whole. Under 'remove', each
            # template lands whole in the alignment file or nowhere.
            bad="$(
                print_split_templates \
                    "${sam_raw}" yes "${stm}.separate.bam" "${arr_dia[@]}"
            )"
            bad_rm="$(
                print_split_templates "${sam_raw}" no "${stm}.remove.bam"
            )"
            if [[ -z "${bad}" && -z "${bad_rm}" ]]; then
                record_pass "${lbl}: every template lands whole, once"
            else
                record_fail \
                    "${lbl}: templates split, lost, or not from the aligner:" \
                    "'${bad}${bad_rm}'"
            fi

            # Invariant: the side file holds only unaligned reads and their
            # mates.
            n_bad="$(
                run_samtools view -c --expr 'flag & 12 == 0' "${fil_side}"
            )"
            n_side="$(run_samtools view -c "${fil_side}")"
            if [[ "${n_bad}" -eq 0 ]]; then
                record_pass \
                    "${lbl}: the side file's ${n_side} records are unaligned" \
                    "reads and their mates"
            else
                record_fail \
                    "${lbl}: the side file holds ${n_bad} records that are" \
                    "neither unaligned nor mates of unaligned reads"
            fi

            # Invariants: the improper file holds only pairs with both mates
            # aligned, and none has two primary records still flagged proper.
            # Without '--req_flg', it is empty. Every template in the low-MAPQ
            # file has a primary record below '--mapq' (the default, 1).
            if [[ "${typ}" == "pe" ]]; then
                n_imp="$(run_samtools view -c "${fil_imp}")"
                n_bad="$(
                    run_samtools view -c --expr 'flag & 12' "${fil_imp}"
                )"
                n_pp="$(
                    run_samtools view -F 2308 -f 2 "${fil_imp}" \
                        | cut -f 1 \
                        | sort \
                        | uniq -d \
                        | wc -l \
                        | tr -d ' '
                )"
                if [[
                    "${n_bad}" -eq 0 && "${n_pp}" -eq 0
                        && ( "${req}" == "yes" || "${n_imp}" -eq 0 )
                ]]; then
                    record_pass \
                        "${lbl}: the improper file's ${n_imp} records are" \
                        "improper pairs with both mates aligned"
                else
                    record_fail \
                        "${lbl}: the improper file holds ${n_bad} records" \
                        "with an unaligned read or mate and ${n_pp} proper" \
                        "pairs (${n_imp} records)"
                fi
            fi

            n_mpq="$(run_samtools view -c "${fil_mpq}")"
            n_bad="$(
                LC_ALL=C comm -23 \
                    <(
                        run_samtools view "${fil_mpq}" \
                            | cut -f 1 \
                            | LC_ALL=C sort -u
                    ) \
                    <(
                        run_samtools view -F 2308 --expr 'mapq < 1' \
                            "${fil_mpq}" \
                            | cut -f 1 \
                            | LC_ALL=C sort -u
                    ) \
                    | wc -l \
                    | tr -d ' '
            )"

            if [[ "${n_bad}" -eq 0 ]]; then
                record_pass \
                    "${lbl}: the low-MAPQ file's ${n_mpq} records belong to" \
                    "templates with a primary record below MAPQ 1"
            else
                record_fail \
                    "${lbl}: the low-MAPQ file holds ${n_bad} templates" \
                    "whose primary records all reach MAPQ 1"
            fi

            # The table of counts: its rows, in order, and its 'all' rows equal
            # to each file's own counts, summing to the raw records.
            tsv="${stm}.separate.diag.tsv"
            gzip -cd "${tsv}.gz" > "${tsv}" 2> /dev/null || true
            got="$(cut -f 1,2 "${tsv}" | tr '\t\n' '::')"
            if [[ "${typ}" == "pe" ]]; then
                exp="${rows_pe}"
            else
                exp="${rows_se}"
            fi
            n_raw="$(run_samtools view -c "${sam_raw}")"
            n_sum="$(
                awk -F '\t' '$2 == "all" { n += $4 } END { print n + 0 }' \
                    "${tsv}"
            )"
            ok_tsv=true
            for fil in "${stm}.separate.bam" "${arr_dia[@]}"; do
                rol="${fil##*.separate.}"
                rol="${rol%.bam}"
                if [[ "${rol}" == "bam" ]]; then rol="main"; fi
                n_rec="$(run_samtools view -c "${fil}")"
                n_tpl="$(
                    run_samtools view "${fil}" \
                        | cut -f 1 \
                        | sort -u \
                        | wc -l \
                        | tr -d ' '
                )"
                row="$(
                    awk -F '\t' -v r="${rol}" \
                        '$1 == r && $2 == "all" { print $3 " " $4 }' "${tsv}"
                )"
                if [[ "${row}" != "${n_tpl} ${n_rec}" ]]; then
                    ok_tsv=false
                fi
            done

            if [[
                "${got}" == "${exp}" && "${n_sum}" == "${n_raw}"
                    && "${ok_tsv}" == "true"
            ]]; then
                record_pass \
                    "${lbl}: the counts' rows are in order, each file's" \
                    "'all' row matches it, and they sum to the ${n_raw} raw" \
                    "records"
            else
                record_fail \
                    "${lbl}: the counts' rows are '${got}', their 'all'" \
                    "records sum to ${n_sum} of ${n_raw}, and each file's" \
                    "'all' row matches it: ${ok_tsv}"
            fi

            # Under 'separate', Step 2 notes how many paired-end templates it
            # moved: those with raw records on both sides of Step 1's split,
            # some with an unaligned read or mate and some without. In
            # 'bridge', BWA-backtrack splits one. 'remove', and single-end
            # data, move nothing and print no note.
            got="$(
                grep -oE 'moved [0-9]+ template' "${stm}.separate.log" \
                    | grep -oE '[0-9]+' \
                    || true
            )"
            if [[ "${typ}" == "pe" ]]; then
                n_exp="$(
                    LC_ALL=C comm -12 \
                        <(
                            run_samtools view -F 12 "${sam_raw}" \
                                | cut -f 1 \
                                | LC_ALL=C sort -u
                        ) \
                        <(
                            run_samtools view --expr 'flag & 12' "${sam_raw}" \
                                | cut -f 1 \
                                | LC_ALL=C sort -u
                        ) \
                        | wc -l \
                        | tr -d ' '
                )"
                if [[
                    "${inp} ${cal}" == "bridge aln" && "${n_exp}" -ne 1
                ]]; then
                    record_fail \
                        "${lbl}: the raw output splits ${n_exp} templates," \
                        "expected 1"
                fi
            else
                n_exp=""
            fi

            n_rm="$(
                grep -c 'moved [0-9]* template' "${stm}.remove.log"
            )" || true
            if [[ "${got}" == "${n_exp}" && "${n_rm}" -eq 0 ]]; then
                record_pass \
                    "${lbl}: Step 2 notes the templates it moved:" \
                    "'${got:-none}'"
            else
                record_fail \
                    "${lbl}: Step 2 notes '${got:-none}' moved templates," \
                    "expected '${n_exp:-none}', and none under 'remove'"
            fi

            # The side file's read names, where pinned above.
            if [[
                "${req}" == "no" && -n "${side_exp[${inp} ${cal}]:-}"
            ]]; then
                got="$(print_names "${fil_side}")"
                if [[ "${got}" == "${side_exp[${inp} ${cal}]}" ]]; then
                    record_pass "${lbl}: the side file holds '${got% }'"
                else
                    record_fail \
                        "${lbl}: the side file holds '${got}', expected" \
                        "'${side_exp[${inp} ${cal}]}'"
                fi
            fi

            # With the default Bowtie 2 PE call, the side file holds exactly
            # the pairs Bowtie 2 reports as not aligned concordantly.
            if [[
                "${cal}" == "bt2" && "${typ}" == "pe" && "${req}" == "no"
            ]]; then
                n_cnc="$(
                    awk \
                        '/aligned concordantly 0 times/ { print $1; exit }' \
                        "${stm}.separate.log"
                )"
                n_tpl="$(
                    run_samtools view "${fil_side}" \
                        | cut -f 1 \
                        | sort -u \
                        | wc -l \
                        | tr -d ' '
                )"
                if [[ -n "${n_cnc}" && "${n_cnc}" -eq "${n_tpl}" ]]; then
                    record_pass \
                        "${lbl}: the side file's ${n_tpl} pairs match" \
                        "Bowtie 2's summary"
                else
                    record_fail \
                        "${lbl}: the side file has ${n_tpl} pairs; Bowtie 2" \
                        "reports '${n_cnc}' not aligned concordantly"
                fi
            fi
        done

        # Invariants: '--req_flg' never changes the side file, duplicate flags
        # included, since Markdup 2 judges the whole library either way. It
        # only moves records between the other three files, which together hold
        # the same records with or without it.
        if [[ "${typ}" == "pe" ]]; then
            stm="${tmp}/${inp}.${cal}"
            fil_no="${stm}.no.separate.unal.bam"
            fil_yes="${stm}.yes.separate.unal.bam"
            if [[ -s "${fil_no}" && -s "${fil_yes}" ]]; then
                if \
                    cmp -s \
                        <(print_recs "${fil_no}") \
                        <(print_recs "${fil_yes}") \
                    && cmp -s \
                        <(
                            for sfx in bam improp.bam mapq.bam; do
                                print_recs_nodup "${stm}.no.separate.${sfx}"
                            done | sort
                        ) \
                        <(
                            for sfx in bam improp.bam mapq.bam; do
                                print_recs_nodup "${stm}.yes.separate.${sfx}"
                            done | sort
                        )
                then
                    record_pass \
                        "${inp} ${cal}: '--req_flg' leaves the side file as" \
                        "it is and only moves records between the others"
                else
                    record_fail \
                        "${inp} ${cal}: '--req_flg' changes the side file or" \
                        "the records the other files hold"
                fi
            fi
        fi
    done
done


# BWA flags no pair proper on 'pairs', too few to estimate insert sizes, so
# '--req_flg' would leave its alignment file empty. On 'pairs_fill', BWA-MEM
# and BWA-MEM2 flag every filler pair and the concordant classes proper, and
# not 'far', 'dove', 'rf', 'ff', or 'chrs'. With '--req_flg', 'separate' keeps
# exactly the proper pairs in the alignment file, as 'remove' does.
cls_all="pair_chrs pair_conc pair_conc2 pair_cont pair_dove pair_far"
cls_all+=" pair_ff pair_ovlp pair_rf pair_rpt "
cls_pp="pair_conc pair_conc2 pair_cont pair_ovlp pair_rpt "
for cal in mem bm2; do
    n_pp="$(
        run_samtools view -c -f 2 -F 2304 \
            "${tmp}/pairs_fill.${cal}.raw.sam" 2> /dev/null
    )" || n_pp=0
    if [[ "${n_pp}" -ge 600 ]]; then
        record_pass \
            "pairs_fill ${cal}: the aligner flags ${n_pp} records proper"
    else
        record_fail \
            "pairs_fill ${cal}: the aligner flags '${n_pp}' records proper," \
            "expected at least 600"
    fi

    for req in no yes; do
        exp="${cls_all}"
        if [[ "${req}" == "yes" ]]; then exp="${cls_pp}"; fi
        fil_main="${tmp}/pairs_fill.${cal}.${req}.separate.bam"
        got="$(
            run_samtools view "${fil_main}" 2> /dev/null \
                | cut -f 1 \
                | sort -u \
                | grep -v '^fill_' \
                | tr '\n' ' '
        )" || true
        n_fil="$(
            run_samtools view "${fil_main}" 2> /dev/null \
                | cut -f 1 \
                | sort -u \
                | grep -c '^fill_'
        )" || true

        if [[ "${got}" == "${exp}" && "${n_fil}" == 300 ]]; then
            record_pass \
                "pairs_fill ${cal} (req_flg ${req}): the alignment file" \
                "holds '${exp% }' and the 300 filler pairs"
        else
            record_fail \
                "pairs_fill ${cal} (req_flg ${req}): the alignment file" \
                "holds '${got}' and ${n_fil} filler pairs, expected '${exp}'" \
                "and 300"
        fi
    done
done


# With '--bt2_relaxed', every pair with both mates aligned stays in the
# alignment file with Bowtie 2's own flags, whatever its class, including the
# repeat-mate pair that mate rescue places; with '--req_flg', only the proper
# pairs stay. 'fixmate' changes TLEN for some classes, not flags. Each name's
# flags are compared sorted, since the files' record orders differ.
fil_raw="${tmp}/pairs.bt2_relaxed.raw.sam"
for req in no yes; do
    fil_main="${tmp}/pairs.bt2_relaxed.${req}.separate.bam"
    got="$(
        {
            run_samtools view "${fil_raw}" \
                | awk -F '\t' '{ print "R\t" $1 "\t" $2 }'
            run_samtools view "${fil_main}" \
                | awk -F '\t' '
                    { print "M\t" $1 "\t" $2 - int($2 / 1024) % 2 * 1024 }
                '
        } \
            | sort -t $'\t' -k 1,1 -k 2,2 -k 3,3n \
            | awk -F '\t' '
                $1 == "R" { raw[$2] = raw[$2] " " $3; next }
                { got[$2] = got[$2] " " $3 }
                END {
                    for (nam in got) {
                        out = out (got[nam] == raw[nam] ? "" : "!") nam " "
                    }
                    print out
                }
            ' \
            | tr ' ' '\n' \
            | sed '/^$/d' \
            | sort \
            | tr '\n' ' '
    )"

    if [[ "${req}" == "no" ]]; then
        exp="pair_chrs pair_conc pair_conc2 pair_cont pair_dove pair_far"
        exp+=" pair_ff pair_ovlp pair_rf pair_rpt "
    else
        exp="pair_conc pair_conc2 pair_cont pair_ovlp pair_rpt "
    fi

    if [[ "${got}" == "${exp}" ]]; then
        record_pass \
            "bt2_relaxed (req_flg ${req}): the alignment file keeps" \
            "'${exp% }' with Bowtie 2's flags"
    else
        record_fail \
            "bt2_relaxed (req_flg ${req}): the alignment file holds" \
            "'${got}', expected '${exp}' ('!' marks changed flags)"
    fi
done


# Duplicates in the side file, from Markdup 2 with BWA-MEM:
# - 'dupA1': its aligned mate matches the first mate of 'dupB' and is marked.
# - 'dupA2' and 'dupA3' share a position but sit far apart on the tile; one is
#   marked.
# - 'dupA4' and 'dupA5' sit within 2,500 pixels; one is marked as optical.
#
# Markdup 1 reads the alignment file alone, so its statistics count none of
# them, as under 'remove'. Markdup 2's statistics and the table of counts count
# all 3.
stm="${tmp}/dup.mem.no"
got="$(
    run_samtools view -f 1024 "${stm}.separate.unal.bam" \
        | awk -F '\t' '
            {
                dt = ""
                for (i = 12; i <= NF; i++) {
                    if ($i ~ /^dt:Z:/) { dt = substr($i, 6) }
                }
                print $4 ":" dt
            }
        ' \
        | sort \
        | tr '\n' ' '
)"
if [[ "${got}" == "1501:LB 2201:LB 2601:SQ " ]]; then
    record_pass "dup mem: the side file's duplicates are marked as expected"
else
    record_fail \
        "dup mem: side file duplicates '${got}', expected" \
        "'1501:LB 2201:LB 2601:SQ '"
fi

n_rm="$(
    gzip -cd "${stm}.remove.markdup.txt.gz" \
        | awk -F ': ' '$1 == "DUPLICATE SINGLE" { print $2 }'
)"
n_sp="$(
    gzip -cd "${stm}.separate.markdup.txt.gz" \
        | awk -F ': ' '$1 == "DUPLICATE SINGLE" { print $2 }'
)"
n_all="$(
    gzip -cd "${stm}.separate.all.markdup.txt.gz" \
        | awk -F ': ' '$1 == "DUPLICATE SINGLE" { print $2 }'
)"
row="$(
    awk -F '\t' '$1 == "unal" && $2 == "duplicate_all" { print $3 " " $4 }' \
        <(gzip -cd "${stm}.separate.diag.tsv.gz")
)"
if [[
    "${n_rm}" == "0"
        && "${n_sp}" == "0"
        && "${n_all}" == "3"
        && "${row}" == "3 3"
]]; then
    record_pass \
        "dup mem: Markdup 1 counts no duplicates, and Markdup 2 and the" \
        "table count the side file's 3"
else
    record_fail \
        "dup mem: duplicates counted ${n_rm} under 'remove', ${n_sp} by" \
        "Markdup 1, ${n_all} by Markdup 2, and '${row}' in the table;" \
        "expected 0, 0, 3, and '3 3'"
fi


# Each diagnostic file carries the alignment file's record, then what it holds,
# its warning, and that Markdup 2 marked it. Markdup 2's optical line replaces
# Markdup 1's. BWA-MEM flags no 'pairs' pair proper, so with '--req_flg' every
# aligned pair is improper.
stm="${tmp}/pairs.mem.yes.separate"
declare -A arr_cnt_exp=(
    [unal]="unaligned reads and their mates, aligned or not"
    [improp]="pairs with both mates aligned that --req_flg rejects"
    [mapq]="templates that --mapq rejects and --req_flg does not"
)
re_opt=$'^@CO\tpipeline: (optical_distance=|skipped=optical duplicates)'
lin_main="$(print_rec_lines "${stm}.bam" | sed '$d')"
for kind in unal improp mapq; do
    exp="${lin_main}"
    exp+=$'\n@CO\tpipeline: content='"${arr_cnt_exp[${kind}]}"
    exp+=$'\n@CO\tpipeline: warning=not for signal or scaling'
    exp+=$'\n@CO\tpipeline: markdup=all files (Markdup 2)'
    got="$(print_rec_lines "${stm}.${kind}.bam")"
    lst="$(tail -n 1 <<< "${got}")"

    if [[
        "$(sed '$d' <<< "${got}")" == "${exp}"
            && "${lst}" =~ ${re_opt}
            && "$(print_rec_lines "${stm}.bam" | tail -n 1)" =~ ${re_opt}
    ]]; then
        record_pass "the '${kind}' file carries the record and what it holds"
    else
        record_fail "the '${kind}' file's record differs: '${got}'"
    fi
done

if \
    run_samtools view -H "${stm}.bam" \
        | grep -q -F -x $'@CO\tpipeline: --unal_mode=separate'
then
    record_pass "the alignment file records '--unal_mode=separate'"
else
    record_fail "the alignment file does not record '--unal_mode=separate'"
fi


# CRAM output with '--qname': each diagnostic file follows the alignment file
# into CRAM, with an index and a queryname-sorted copy. It holds the BAM run's
# records, and no BAM file is left. The relaxed call with '--req_flg' fills the
# improper file.
stm="${tmp}/cram"
rc=0
out="$(
    run_helper \
        --threads 1 \
        --aligner bowtie2 \
        --bt2_relaxed \
        --req_flg \
        --index "${dir_fx}/bowtie2/pairs" \
        --unal_mode separate \
        --ref_fa "${ref_par}" \
        --fq_1 "${dir_fx}/fastq/pe/pairs_pe_R1.atria.fastq.gz" \
        --fq_2 "${dir_fx}/fastq/pe/pairs_pe_R2.atria.fastq.gz" \
        --qname \
        --fil_out "${stm}.cram"
)" || rc=$?

ok=true
if [[ "${rc}" -ne 0 || ! -s "${stm}.diag.tsv.gz" ]]; then ok=false; fi
for kind in unal improp mapq; do
    if [[
        ! -s "${stm}.${kind}.cram"
            || ! -s "${stm}.${kind}.cram.crai"
            || ! -s "${stm}.${kind}.qnam.cram"
            || -e "${stm}.${kind}.bam"
            || -e "${stm}.${kind}.qnam.bam"
    ]]; then
        ok=false
    fi
done

if [[ "${ok}" == "true" ]]; then
    record_pass \
        "CRAM with '--qname' writes each diagnostic file, its index, and" \
        "copy, and no BAM file"
else
    printf '%s\n' "${out}" > "${dir_log}/unal_mode_cram.log"
    record_fail \
        "CRAM with '--qname' did not write the diagnostic files as CRAM" \
        "(exit ${rc})"
fi

stm_bam="${tmp}/pairs.bt2_relaxed.yes.separate"
ok=true
for kind in unal improp mapq; do
    if ! \
        {
            cmp -s \
                <(print_recs_nomd -T "${ref_par}" "${stm}.${kind}.cram") \
                <(print_recs_nomd "${stm_bam}.${kind}.bam") \
            && cmp -s \
                <(print_recs_nomd -T "${ref_par}" "${stm}.${kind}.qnam.cram") \
                <(print_recs_nomd "${stm_bam}.${kind}.bam") \
            && run_samtools view -H "${stm}.${kind}.qnam.cram" \
                | grep -q -F -x \
                    $'@CO\tpipeline: content='"${arr_cnt_exp[${kind}]}"
        }
    then
        ok=false
    fi
done

if \
    [[ "${ok}" == "true" ]] \
        && cmp -s "${stm}.diag.tsv.gz" "${stm_bam}.diag.tsv.gz"
then
    record_pass \
        "the CRAM diagnostic files and copies hold the BAM run's records," \
        "the copies carry their records, and the counts match"
else
    record_fail \
        "a CRAM diagnostic file, its copy, or the counts differ from the" \
        "BAM run's"
fi


# Determinism (#61): each diagnostic file is the same at two thread counts, in
# records and content lines. So are the '@PG' IDs that 'samtools merge' makes
# unique in the diagnostic files, since their suffixes are seeded. A suffix
# drawn from the clock would change with the pause between the runs. The
# alignment file is not merged, and its IDs match as well.
for cal in bt2 mem; do
    if [[ "${cal}" == "bt2" ]]; then
        args=( --aligner bowtie2 --index "${dir_fx}/bowtie2/pairs" )
    else
        args=( --aligner bwa --index "${dir_fx}/bwa/pairs.fa" )
    fi

    for thr in 1 3; do
        run_helper \
            --threads "${thr}" \
            "${args[@]}" \
            --unal_mode separate \
            --fq_1 "${dir_fx}/fastq/pe/pairs_pe_R1.atria.fastq.gz" \
            --fq_2 "${dir_fx}/fastq/pe/pairs_pe_R2.atria.fastq.gz" \
            --fil_out "${tmp}/det.${cal}.${thr}.bam" > /dev/null || true
        if [[ "${thr}" == "1" ]]; then sleep 1.1; fi
    done

    ok=true
    for sfx in bam unal.bam improp.bam mapq.bam; do
        ids_1="$(print_pg_ids "${tmp}/det.${cal}.1.${sfx}")"
        ids_3="$(print_pg_ids "${tmp}/det.${cal}.3.${sfx}")"

        if [[ -z "${ids_1}" || "${ids_1}" != "${ids_3}" ]]; then ok=false; fi
        if [[ "${sfx}" != "bam" && "${ids_1}" != *"-"* ]]; then ok=false; fi
    done

    if [[ "${ok}" == "true" ]]; then
        record_pass "${cal}: every file's '@PG' IDs match between runs"
    else
        record_fail "${cal}: the '@PG' IDs differ between runs"
    fi

    ok=true
    for sfx in unal.bam improp.bam mapq.bam; do
        if ! \
            {
                cmp -s \
                    <(run_samtools view "${tmp}/det.${cal}.1.${sfx}") \
                    <(run_samtools view "${tmp}/det.${cal}.3.${sfx}") \
                && cmp -s \
                    <(
                        run_samtools view -H --no-PG \
                            "${tmp}/det.${cal}.1.${sfx}" \
                            | grep $'^@CO\tpipeline: '
                    ) \
                    <(
                        run_samtools view -H --no-PG \
                            "${tmp}/det.${cal}.3.${sfx}" \
                            | grep $'^@CO\tpipeline: '
                    )
            }
        then
            ok=false
        fi
    done

    if [[ "${ok}" == "true" ]]; then
        record_pass \
            "${cal}: the diagnostic files are the same at 1 and 3 threads"
    else
        record_fail "${cal}: a diagnostic file differs between 1 and 3 threads"
    fi
done


# Determinism (#61) on the thread fixture. Its BWA-MEM output depends on the
# thread count unless '-K' fixes the batch. With '-K', one batch holds all its
# pairs, so the 20,000 long pairs are improper at every thread count. If the
# batch followed the thread count, they would be proper at one thread and move
# into the alignment file. So each file's records, in file order, and the
# counts must match at 1 and 3 threads.
fq_thr_1="${dir_fx}/fastq/pe/thr_pe_R1.atria.fastq.gz"
fq_thr_2="${dir_fx}/fastq/pe/thr_pe_R2.atria.fastq.gz"
for thr in 1 3; do
    rc=0
    out="$(
        run_helper \
            --threads "${thr}" \
            --aligner bwa \
            --index "${dir_fx}/bwa/thr.fa" \
            --req_flg \
            --unal_mode separate \
            --fq_1 "${fq_thr_1}" \
            --fq_2 "${fq_thr_2}" \
            --fil_out "${tmp}/thr.${thr}.bam"
    )" || rc=$?
    if [[ "${rc}" -ne 0 ]]; then
        printf '%s\n' "${out}" > "${dir_log}/unal_mode_thr_${thr}.log"
        record_fail "thr: the run at ${thr} thread(s) failed (exit ${rc})"
    fi
done

n_imp="$(
    awk -F '\t' '$1 == "improp" && $2 == "all" { print $3 }' \
        <(gzip -cd "${tmp}/thr.1.diag.tsv.gz" 2> /dev/null)
)" || true
if [[ "${n_imp}" == "20000" ]]; then
    record_pass "thr: the improper file holds the 20,000 long pairs"
else
    record_fail \
        "thr: the improper file holds '${n_imp}' pairs, expected the" \
        "20,000 long pairs"
fi

ok=true
for sfx in bam unal.bam improp.bam mapq.bam; do
    if ! \
        cmp -s \
            <(run_samtools view "${tmp}/thr.1.${sfx}") \
            <(run_samtools view "${tmp}/thr.3.${sfx}")
    then
        ok=false
    fi
done

if \
    [[ "${ok}" == "true" ]] \
        && cmp -s "${tmp}/thr.1.diag.tsv.gz" "${tmp}/thr.3.diag.tsv.gz"
then
    record_pass "thr: every file and the counts are the same at 1 and 3 threads"
else
    record_fail "thr: a file or the counts differ between 1 and 3 threads"
fi

# The alignment file's 'duplicate_all' row, which comes from Markdup 2's
# statistics, checked against a recount. The recount clears the four files'
# flags and marks them together again, in the order Step 4 merges them. It must
# mark the same records of the alignment file.
stm="${tmp}/thr.1"
run_samtools merge -s 0 -o "${tmp}/thr.dia.bam" \
    "${stm}.unal.bam" "${stm}.improp.bam" "${stm}.mapq.bam"
run_samtools view "${tmp}/thr.dia.bam" \
    | cut -f 1 \
    | sort -u \
    > "${tmp}/thr.dia.names.txt"
n_rec="$(
    run_samtools merge -s 0 -u -o - "${stm}.bam" "${tmp}/thr.dia.bam" \
        | run_samtools markdup -c -u - - 2> /dev/null \
        | run_samtools view -N "^${tmp}/thr.dia.names.txt" -f 1024 - \
        | tee "${tmp}/thr.dup.sam" \
        | wc -l \
        | tr -d ' '
)" || true
n_tpl="$(
    awk -F '\t' 'int($2 / 256) % 2 + int($2 / 2048) % 2 == 0 { print $1 }' \
        "${tmp}/thr.dup.sam" \
        | sort -u \
        | wc -l \
        | tr -d ' '
)"
row="$(
    awk -F '\t' '$1 == "main" && $2 == "duplicate_all" { print $3 " " $4 }' \
        <(gzip -cd "${stm}.diag.tsv.gz")
)"
if [[ "${row}" == "${n_tpl} ${n_rec}" && "${n_rec}" -gt 0 ]]; then
    record_pass \
        "thr: the alignment file's Markdup 2 duplicates ('${row}') match a" \
        "recount"
else
    record_fail \
        "thr: the alignment file's Markdup 2 duplicates are '${row}', a" \
        "recount gives '${n_tpl} ${n_rec}'"
fi


# The table's categories, pinned on the relaxed call. 'pair_none' has both
# mates unaligned and 'pair_one' one mate, so the side file holds 3 unaligned
# records. With '--req_flg', 'pair_chrs' is improper by reference, 'pair_ff' by
# orientation, 'pair_rf' by facing outward, 'pair_dove' by dovetailing, and
# 'pair_far' by spanning more than '--bt2_X' (1000).
#
# Bowtie 2's highest MAPQ in end-to-end mode is 42, which every aligned pair
# here gets. So '--mapq 43' fails every primary record. The improper pairs then
# fail both filters and must stay improper (precedence), each counted in the
# low-MAPQ columns. The proper pairs go to the low-MAPQ file, and the alignment
# file is empty.
stm="${tmp}/rows"
run_helper \
    --threads 1 \
    --aligner bowtie2 \
    --bt2_relaxed \
    --req_flg \
    --mapq 43 \
    --index "${dir_fx}/bowtie2/pairs" \
    --unal_mode separate \
    --fq_1 "${dir_fx}/fastq/pe/pairs_pe_R1.atria.fastq.gz" \
    --fq_2 "${dir_fx}/fastq/pe/pairs_pe_R2.atria.fastq.gz" \
    --fil_out "${stm}.bam" > "${stm}.log" 2>&1 || true

exp="main:all:0:0:NA:NA main:unjudged:0:0:NA:NA"
exp+=" main:duplicate_all:0:0:NA:NA unal:all:2:4:NA:NA"
exp+=" unal:both_unaligned:1:2:NA:NA unal:one_unaligned:1:2:NA:NA"
exp+=" unal:unjudged:2:3:NA:NA unal:duplicate_all:0:0:NA:NA"
exp+=" improp:all:5:10:5:10 improp:diff_ref:1:2:1:2"
exp+=" improp:orient_ff:1:2:1:2 improp:orient_rr:0:0:0:0"
exp+=" improp:orient_rf:1:2:1:2 improp:dovetailed:1:2:1:2"
exp+=" improp:too_far:1:2:1:2 improp:too_close:0:0:0:0"
exp+=" improp:failed_est:NA:NA:NA:NA improp:other:0:0:0:0"
exp+=" improp:unjudged:0:0:0:0 improp:duplicate_all:0:0:0:0"
exp+=" mapq:all:5:10:NA:NA mapq:unjudged:0:0:NA:NA"
exp+=" mapq:duplicate_all:0:0:NA:NA "
got="$(
    gzip -cd "${stm}.diag.tsv.gz" 2> /dev/null | sed 1d | tr '\t\n' ': '
)" || true
nam_dr="$(
    run_samtools view -F 2308 --expr 'refid != mrefid' "${stm}.improp.bam" \
        | cut -f 1 \
        | sort -u \
        | tr '\n' ' '
)"
nam_ss="$(
    run_samtools view -F 2308 --expr 'flag & 48 == 0 || flag & 48 == 48' \
        "${stm}.improp.bam" \
        | cut -f 1 \
        | sort -u \
        | tr '\n' ' '
)"

if [[
    "${got}" == "${exp}"
        && "${nam_dr}" == "pair_chrs "
        && "${nam_ss}" == "pair_ff "
]]; then
    record_pass "rows: the table counts each category as expected"
else
    record_fail \
        "rows: the table holds '${got}', expected '${exp}'; different" \
        "references '${nam_dr}', same orientation '${nam_ss}'"
fi

# The table on its own, from the hand-made files in 'sam/diag/', at MAPQ 30 and
# a Bowtie 2 '-X' of 500. Their recipe says what each record is. Each of these
# cases would change the table if the count were wrong:
# - 'i_xr' has both mates forward, but different references come first.
# - A pair at MAPQ 30 exactly is not low, and its supplementary record at MAPQ
#   0 neither makes it low nor gives it a reason.
# - 'i_q10' spans 510 bases, more than '-X', and is low.
# - 'i_dov' is dovetailed only because its reverse mate starts first, and
#   'i_dov_fe' only because its forward mate ends last.
# - Unaligned, secondary, and supplementary records are unjudged in every file.
# - The main file's Markdup 2 duplicates come from the statistics, less the
#   side file's one.
dir_dia="${dir_fx}/sam/diag"
arr_dia=(
    "${dir_dia}/main.sam" "${dir_dia}/unal.sam" "${dir_dia}/improp.sam"
    "${dir_dia}/mapq.sam"
)
exp_all="main:all:2:5:NA:NA main:unjudged:1:1:NA:NA"
exp_all+=" main:duplicate_all:1:2:NA:NA unal:all:3:6:NA:NA"
exp_all+=" unal:both_unaligned:1:2:NA:NA unal:one_unaligned:2:4:NA:NA"
exp_all+=" unal:unjudged:3:4:NA:NA unal:duplicate_all:1:1:NA:NA"
exp_all+=" improp:all:8:17:1:2 improp:diff_ref:1:2:0:0"
exp_all+=" improp:orient_ff:1:2:0:0 improp:orient_rr:1:2:0:0"
exp_all+=" improp:orient_rf:1:2:0:0 improp:dovetailed:2:4:0:0"
exp_bt2="${exp_all} improp:too_far:1:2:1:2 improp:too_close:0:0:0:0"
exp_bt2+=" improp:failed_est:NA:NA:NA:NA improp:other:1:3:0:0"
exp_bwa="${exp_all} improp:too_far:NA:NA:NA:NA"
exp_bwa+=" improp:too_close:NA:NA:NA:NA improp:failed_est:NA:NA:NA:NA"
exp_bwa+=" improp:other:2:5:1:2"
exp_end=" improp:unjudged:1:1:0:0 improp:duplicate_all:0:0:0:0"
exp_end+=" mapq:all:1:3:NA:NA mapq:unjudged:1:1:NA:NA"
exp_end+=" mapq:duplicate_all:0:0:NA:NA "

for al in bowtie2 bwa; do
    rc=0
    run_func _write_diag_tsv \
        align_fastqs false 1 "${al}" 500 pe 30 "${arr_dia[@]}" \
        "${dir_dia}/all.markdup.txt.gz" \
        "${tmp}/fx.${al}.diag.tsv.gz" > "${tmp}/fx.${al}.diag.log" 2>&1 \
        || rc=$?

    if [[ "${al}" == "bowtie2" ]]; then
        exp="${exp_bt2}${exp_end}"
    else
        exp="${exp_bwa}${exp_end}"
    fi
    got="$(
        gzip -cd "${tmp}/fx.${al}.diag.tsv.gz" 2> /dev/null \
            | sed 1d \
            | tr '\t\n' ': '
    )" || true
    hdr="$(
        gzip -cd "${tmp}/fx.${al}.diag.tsv.gz" 2> /dev/null \
            | head -n 1 \
            | tr '\t' ' '
    )" || true

    if [[
        "${rc}" -eq 0 && "${got}" == "${exp}"
            && "${hdr}" == "file category templates records"*" records_mapq_low"
    ]]; then
        record_pass "the table from the hand-made files is as expected (${al})"
    else
        record_fail \
            "the table from the hand-made files is '${got}' with header" \
            "'${hdr}' (${al}, exit ${rc}), expected '${exp}'"
    fi
done

# Statistics the other files cannot follow from: a duplicate total below the
# side file's duplicate, then an odd count of primary duplicates left for the
# main file. Only that row turns NA, with a warning, and the run goes on.
for bad in neg odd; do
    rc=0
    run_func _write_diag_tsv \
        align_fastqs false 1 bowtie2 500 pe 30 "${arr_dia[@]}" \
        "${dir_dia}/all_${bad}.markdup.txt.gz" \
        "${tmp}/fx.${bad}.diag.tsv.gz" > "${tmp}/fx.${bad}.diag.log" 2>&1 \
        || rc=$?

    got="$(
        gzip -cd "${tmp}/fx.${bad}.diag.tsv.gz" 2> /dev/null \
            | sed 1d \
            | tr '\t\n' ': '
    )" || true
    exp="${exp_bt2}${exp_end}"
    exp="${exp/main:duplicate_all:1:2:NA:NA/main:duplicate_all:NA:NA:NA:NA}"

    if [[
        "${rc}" -eq 0 && "${got}" == "${exp}"
    ]] && grep -q "'duplicate_all' row is NA" "${tmp}/fx.${bad}.diag.log"; then
        record_pass \
            "statistics that do not add up ('${bad}') give an NA row and a" \
            "warning"
    else
        record_fail \
            "statistics that do not add up ('${bad}') give '${got}' (exit" \
            "${rc}), expected the main file's 'duplicate_all' row NA and a" \
            "warning"
    fi
done

# The 'gawk' program on its own: with a minimum fragment of 300 bases, 'i_q30'
# (210 bases) is too close. Without 'MC' on a primary record, or without the
# first mate's primary record of 'i_rf', it stops.
fil_awk="${ROOT_REPO}/lib/awk/classify_improp_pairs.awk"
got="$(
    run_samtools view "${dir_dia}/improp.sam" \
        | gawk -f "${fil_awk}" \
            -v aligner=bowtie2 -v max_frag=500 -v min_frag=300 -v mapq=30 \
        | awk -F '\t' '$2 ~ /^(too_far|too_close|other)$/' \
        | tr '\t\n' ': '
)" || true
exp="improp:too_far:1:2:1:2 improp:too_close:1:3:0:0 improp:other:0:0:0:0 "
rc=0
run_samtools view "${dir_dia}/improp.sam" \
    | sed 's/\tMC:Z:10M$//' \
    | gawk -f "${fil_awk}" -v aligner=bowtie2 -v max_frag=500 -v mapq=30 \
        > /dev/null 2>&1 \
    || rc=$?
rc_r1=0
run_samtools view "${dir_dia}/improp.sam" \
    | awk -F '\t' '!($1 == "i_rf" && $2 == 81)' \
    | gawk -f "${fil_awk}" -v aligner=bowtie2 -v max_frag=500 -v mapq=30 \
        > /dev/null 2>&1 \
    || rc_r1=$?

if [[ "${got}" == "${exp}" && "${rc}" -ne 0 && "${rc_r1}" -ne 0 ]]; then
    record_pass \
        "the improper-pair program counts a short fragment as too close and" \
        "stops without 'MC' or a first mate"
else
    record_fail \
        "the improper-pair program gives '${got}', expected '${exp}', and" \
        "exits ${rc} without 'MC' and ${rc_r1} without a first mate"
fi


# Step 2 on its own, with '--req_flg', on a SAM file of these templates:
# - 't75': a #75 template.
# - 'supp': the same, with a supplementary record.
# - 'pair': a complete pair.
# - 'lone': an aligned R1 with no R2 record.
# - 'ffpp' and 'xchr': a same-orientation pair and a cross-reference pair, both
#   flagged proper.
# - 'lowq': a proper pair with a mate at MAPQ 0.
# - 'ffmq': a same-orientation pair flagged proper at MAPQ 0.
# - 'disc': a pair the aligner did not flag proper.
#
# The test makes Step 1's split itself, as Step 1 would. 'fixmate' unflags
# 'ffpp', 'ffmq', and 'xchr', so '--req_flg' rejects them after it in both
# modes. Under 'separate', both split templates move whole into Set U, and
# every rejected template is kept. 'lone', 'disc', and the unflagged pairs are
# improper. That includes 'ffmq', though it also fails MAPQ 30. 'lowq' is
# low-MAPQ at MAPQ 30.
declare -A arr_s2_exp=(
    ["remove 0"]="lowq:2 pair:2 "
    ["remove 30"]="pair:2 "
    ["separate 0"]="lowq:2 pair:2 "
    ["separate 30"]="pair:2 "
)
declare -A arr_s2_set=(
    ["U"]="supp:3 t75:2 "
    ["I"]="disc:2 ffmq:2 ffpp:2 lone:1 xchr:2 "
    ["Q 0"]=""
    ["Q 30"]="lowq:2 "
)
for mod in remove separate; do
    for mpq in 0 30; do
        stm="${tmp}/step2.${mod}.${mpq}"
        if [[ "${mod}" == "separate" ]]; then
            run_samtools view -b -F 12 -U "${stm}.setU.bam" -o "${stm}.bam" \
                "${sam_spl}"
        else
            run_samtools view -b -F 12 -f 2 -o "${stm}.bam" "${sam_spl}"
        fi

        rc=0
        out="$(
            run_func _sort_qname_bam \
                align_fastqs false 1 bwa "${mpq}" true "${mod}" pe \
                "${stm}.bam" "${stm}.qnam.bam"
        )" || rc=$?

        got="$(print_counts "${stm}.qnam.bam")"
        exp="${arr_s2_exp[${mod} ${mpq}]}"
        ok=true
        if [[ "${rc}" -ne 0 || "${got}" != "${exp}" ]]; then ok=false; fi

        if [[ "${mod}" == "separate" ]]; then
            got_u="$(print_counts "${stm}.setU.bam")"
            got_i="$(print_counts "${stm}.setI.bam")"
            got_q="$(print_counts "${stm}.setQ.bam")"
            nam="$(tr '\n' ' ' < "${stm}.unal.names.txt" 2> /dev/null)" \
                || true
            if [[
                "${got_u}" != "${arr_s2_set[U]}"
                    || "${got_i}" != "${arr_s2_set[I]}"
                    || "${got_q}" != "${arr_s2_set[Q ${mpq}]}"
                    || "${nam}" != "supp t75 "
            ]]; then
                ok=false
            fi
            got+="; U '${got_u}'; I '${got_i}'; Q '${got_q}'"
        elif [[ -e "${stm}.setI.bam" || -e "${stm}.setQ.bam" ]]; then
            ok=false
        fi

        if [[ "${ok}" == "true" ]]; then
            record_pass \
                "Step 2 alone under '${mod}' at MAPQ ${mpq} keeps" \
                "'${exp% }'"
        else
            printf '%s\n' "${out}" \
                > "${dir_log}/unal_mode_step2_${mod}_${mpq}.log"
            record_fail \
                "Step 2 alone under '${mod}' at MAPQ ${mpq} kept '${got}'," \
                "expected '${exp}' (exit ${rc})"
        fi

        # The routing moves 2 templates, 'supp' and 't75', with 3 of their
        # records in Set M ('supp' primary and supplementary, 't75' R1). So its
        # note tells templates from records. 'remove' moves nothing.
        got="$(grep -oE 'moved [0-9]+ template.*' <<< "${out}" || true)"
        if [[ "${mod}" == "separate" ]]; then
            exp="moved 2 template(s) split between the work files into Set U"
            exp+=" (3 record(s))."
        else
            exp=""
        fi

        if [[ "${got}" == "${exp}" ]]; then
            record_pass \
                "Step 2 alone under '${mod}' at MAPQ ${mpq} notes" \
                "'${got:-no moved templates}'"
        else
            record_fail \
                "Step 2 alone under '${mod}' at MAPQ ${mpq} notes" \
                "'${got:-nothing}', expected '${exp:-no note}'"
        fi
    done
done


# The AWK BED converter drops unaligned reads and their mates, so a
# queryname-sorted file with the side file's records merged in gives the same
# fragments as the alignment file alone.
stm="${tmp}/pairs.mem.no.separate"
run_samtools sort -n -o "${tmp}/bed.main.bam" "${stm}.bam"
run_samtools sort -n -o "${tmp}/bed.side.bam" "${stm}.unal.bam"
run_samtools merge -n -o "${tmp}/bed.both.bam" \
    "${tmp}/bed.main.bam" "${tmp}/bed.side.bam"

rc=0
for nam in both main; do
    run_func convert_alignments_bed_awk \
        1 "${tmp}/bed.${nam}.bam" "${tmp}/bed.${nam}.bed.gz" "" \
        "${tmp}/bed.${nam}.out" "${tmp}/bed.${nam}.err" > /dev/null \
        || rc=$?
done

n_bed="$(gzip -cd "${tmp}/bed.main.bed.gz" 2> /dev/null | wc -l | tr -d ' ')"
if [[
    "${rc}" -eq 0
        && "${n_bed}" -gt 0
        && -n "$(run_samtools view -c -f 4 "${tmp}/bed.both.bam")"
]] \
    && cmp -s \
        <(gzip -cd "${tmp}/bed.both.bed.gz") \
        <(gzip -cd "${tmp}/bed.main.bed.gz")
then
    record_pass "the awk BED converter gives no fragment for unaligned reads"
else
    record_fail \
        "the awk BED converter makes fragments from unaligned reads" \
        "(exit ${rc}; ${n_bed} fragments from the alignment file)"
fi


finish
