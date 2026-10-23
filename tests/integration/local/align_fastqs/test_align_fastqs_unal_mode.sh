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
# '--unal_mode remove' drops unaligned reads and their mates; 'separate' also
# writes them to a side file and leaves the alignment file as it is. Every
# aligner call runs on the 'pairs' (one pair per pair class), 'pairs_fill'
# (the pair classes among 300 concordant pairs, so BWA can flag pairs proper),
# 'bridge' (the BWA-backtrack lone mate of #75), and 'dup' (duplicates among
# the side file's aligned mates) fixtures, with and without '--req_flg', and
# its records are checked against the aligner's own output.
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
    "${dir_fx}/fastq/pe/bridge_pe_R2.atria.fastq.gz" || {
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
# helper's outputs can be checked against it. Arguments: the call's label, the
# fixture set, the second FASTQ (or empty), and the output SAM.
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
    esac 2> "${sam%.sam}.raw.log" > "${sam}"
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


# Check that each template of the raw output lands whole in one place: the
# alignment file, the side file, or neither (dropped by '--mapq' or
# '--req_flg'). Print nothing when it holds, or the offending names.
function print_split_templates() {
    local sam_raw="${1}"
    local fil_main="${2}"
    local fil_side="${3}"

    {
        print_keys "${sam_raw}" | sed 's/^/R\t/'
        print_keys "${fil_main}" | sed 's/^/M\t/'
        if [[ -n "${fil_side}" ]]; then
            print_keys "${fil_side}" | sed 's/^/S\t/'
        fi
    } | awk -F '\t' '
        {
            key = $2 SUBSEP $3 SUBSEP $4 SUBSEP $5
            if ($1 == "R") { n_raw[$2]++; raw[key] = 1 }
            else {
                n[$2, $1]++
                if (!(key in raw)) { bad[$2] = 1 }
                if (key in seen) { bad[$2] = 1 }
                seen[key] = 1
            }
        }
        END {
            for (nam in n_raw) {
                m = n[nam, "M"] + 0
                s = n[nam, "S"] + 0
                if ((m > 0 && s > 0) || (m > 0 && m != n_raw[nam]) \
                    || (s > 0 && s != n_raw[nam])) {
                    bad[nam] = 1
                }
            }
            for (nam in bad) { print nam }
        }
    ' | sort | tr '\n' ' '
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
        if ! \
            run_raw "${cal}" "${set}" "${fq_1}" "${fq_2}" "${sam_raw}"
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

            # Only 'separate' writes a side file, indexed.
            if [[
                ! -e "${stm}.remove.unal.bam"
                    && -s "${fil_side}"
                    && -s "${fil_side}.bai"
            ]]; then
                record_pass "${lbl}: only 'separate' writes a side file"
            else
                record_fail "${lbl}: side file missing or written by 'remove'"
            fi

            # Invariant: 'remove' and 'separate' give the same alignment file,
            # duplicate flags included.
            if cmp -s \
                <(print_recs "${stm}.remove.bam") \
                <(print_recs "${stm}.separate.bam")
            then
                record_pass \
                    "${lbl}: the alignment file is the same in both modes"
            else
                record_fail "${lbl}: the alignment file differs between modes"
            fi

            # Invariants: every raw record lands once, and each template lands
            # whole, in the alignment file, the side file, or neither.
            bad="$(
                print_split_templates \
                    "${sam_raw}" "${stm}.separate.bam" "${fil_side}"
            )"
            bad_rm="$(
                print_split_templates \
                    "${sam_raw}" "${stm}.remove.bam" ""
            )"
            if [[ -z "${bad}" && -z "${bad_rm}" ]]; then
                record_pass "${lbl}: every template lands whole, once"
            else
                record_fail \
                    "${lbl}: templates split or records not from the" \
                    "aligner: '${bad}${bad_rm}'"
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

            # The side file's read names, where pinned above.
            if [[ "${req}" == "no" && -n "${side_exp[${inp} ${cal}]:-}" ]]
            then
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
            if [[ "${cal}" == "bt2" && "${typ}" == "pe" && "${req}" == "no" ]]
            then
                n_cnc="$(
                    awk '/aligned concordantly 0 times/ { print $1; exit }' \
                        "${stm}.separate.log"
                )"
                n_tpl="$(
                    run_samtools view "${fil_side}" | cut -f 1 | sort -u \
                        | wc -l | tr -d ' '
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

        # Invariant: '--req_flg' never changes the side file's records. Their
        # duplicate flags can change, since they are marked against the whole
        # library, which '--req_flg' changes.
        if [[ "${typ}" == "pe" ]]; then
            stm="${tmp}/${inp}.${cal}"
            fil_no="${stm}.no.separate.unal.bam"
            fil_yes="${stm}.yes.separate.unal.bam"
            if [[ -s "${fil_no}" && -s "${fil_yes}" ]]; then
                if cmp -s \
                    <(print_recs_nodup "${fil_no}") \
                    <(print_recs_nodup "${fil_yes}")
                then
                    record_pass \
                        "${inp} ${cal}: '--req_flg' leaves the side file's" \
                        "records as they are"
                else
                    record_fail \
                        "${inp} ${cal}: '--req_flg' changes the side file's" \
                        "records"
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
                | cut -f 1 | sort -u | grep -v '^fill_' | tr '\n' ' '
        )" || true
        n_fil="$(
            run_samtools view "${fil_main}" 2> /dev/null \
                | cut -f 1 | sort -u | grep -c '^fill_'
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
        } | sort -t $'\t' -k 1,1 -k 2,2 -k 3,3n | awk -F '\t' '
            $1 == "R" { raw[$2] = raw[$2] " " $3; next }
            { got[$2] = got[$2] " " $3 }
            END {
                for (nam in got) {
                    out = out (got[nam] == raw[nam] ? "" : "!") nam " "
                }
                print out
            }
        ' | tr ' ' '\n' | sed '/^$/d' | sort | tr '\n' ' '
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


# Duplicates in the side file: with BWA-MEM, the aligned mate of 'dupA1'
# matches the first mate of 'dupB' and is marked; of 'dupA2' and 'dupA3', which
# share a position far apart on the tile, one is marked; of 'dupA4' and
# 'dupA5', within 2,500 pixels, one is marked as optical. Marking them leaves
# the alignment file's flags as they are (checked above), and the statistics
# count them.
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
        | sort | tr '\n' ' '
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
if [[ "${n_rm}" == "0" && "${n_sp}" == "3" ]]; then
    record_pass "dup mem: the duplicate statistics count the side file's 3"
else
    record_fail \
        "dup mem: duplicate statistics count ${n_rm} with 'remove' and" \
        "${n_sp} with 'separate', expected 0 and 3"
fi


# The side file carries the alignment file's record, plus what it holds.
stm="${tmp}/pairs.mem.no.separate"
got="$(
    diff \
        <(run_samtools view -H --no-PG "${stm}.bam" | grep '^@CO') \
        <(run_samtools view -H --no-PG "${stm}.unal.bam" | grep '^@CO') \
        | grep '^[<>]' || true
)"
exp="> @CO	pipeline: content=unaligned reads and their mates, aligned or not
> @CO	pipeline: warning=not for signal or scaling"
if [[ "${got}" == "${exp}" ]]; then
    record_pass "the side file carries the record and what it holds"
else
    record_fail "the side file's record differs: '${got}'"
fi

if \
    run_samtools view -H "${stm}.bam" \
        | grep -q -F -x $'@CO\tpipeline: --unal_mode=separate'
then
    record_pass "the alignment file records '--unal_mode=separate'"
else
    record_fail "the alignment file does not record '--unal_mode=separate'"
fi


# CRAM output with '--qname': the side file follows the alignment file into
# CRAM, with an index and a queryname-sorted copy, holding the BAM run's
# records.
stm="${tmp}/cram"
rc=0
out="$(
    run_helper \
        --threads 1 \
        --aligner bowtie2 \
        --index "${dir_fx}/bowtie2/pairs" \
        --unal_mode separate \
        --ref_fa "${ref_par}" \
        --fq_1 "${dir_fx}/fastq/pe/pairs_pe_R1.atria.fastq.gz" \
        --fq_2 "${dir_fx}/fastq/pe/pairs_pe_R2.atria.fastq.gz" \
        --qname \
        --fil_out "${stm}.cram"
)" || rc=$?

if [[
    "${rc}" -eq 0
    && -s "${stm}.unal.cram"
    && -s "${stm}.unal.cram.crai"
    && -s "${stm}.unal.qnam.cram"
    && ! -e "${stm}.unal.bam"
    && ! -e "${stm}.unal.qnam.bam"
]]; then
    record_pass "CRAM with '--qname' writes the side file, its index, and copy"
else
    printf '%s\n' "${out}" > "${dir_log}/unal_mode_cram.log"
    record_fail "CRAM with '--qname' did not write the side files (exit ${rc})"
fi

fil_bam="${tmp}/pairs.bt2.no.separate.unal.bam"
if \
    cmp -s \
        <(run_samtools view -T "${ref_par}" "${stm}.unal.cram" | sort) \
        <(print_recs "${fil_bam}") \
    && cmp -s \
        <(run_samtools view -T "${ref_par}" "${stm}.unal.qnam.cram" | sort) \
        <(print_recs "${fil_bam}")
then
    record_pass "the CRAM side file and its copy hold the BAM run's records"
else
    record_fail "the CRAM side file or its copy differs from the BAM run's"
fi

lin_cnt=$'@CO\tpipeline: content=unaligned reads and their mates,'
lin_cnt+=' aligned or not'
if \
    run_samtools view -H "${stm}.unal.qnam.cram" \
        | grep -q -F -x "${lin_cnt}"
then
    record_pass "the side file's '--qname' copy carries its record"
else
    record_fail "the side file's '--qname' copy lacks its record"
fi


# Determinism (#61): the side file is the same at two thread counts, records
# and content lines alike. In both files, the '@PG' IDs that 'samtools merge'
# makes unique are the same too: their suffixes are seeded, not drawn from the
# clock, which the pause between the runs would change.
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
    for sfx in bam unal.bam; do
        ids_1="$(print_pg_ids "${tmp}/det.${cal}.1.${sfx}")"
        ids_3="$(print_pg_ids "${tmp}/det.${cal}.3.${sfx}")"

        if [[ "${ids_1}" != *"-"* || "${ids_1}" != "${ids_3}" ]]; then
            ok=false
        fi
    done

    if [[ "${ok}" == "true" ]]; then
        record_pass "${cal}: both files' merged '@PG' IDs match between runs"
    else
        record_fail "${cal}: the merged '@PG' IDs differ between runs"
    fi

    if \
        cmp -s \
            <(run_samtools view "${tmp}/det.${cal}.1.unal.bam") \
            <(run_samtools view "${tmp}/det.${cal}.3.unal.bam") \
        && cmp -s \
            <(
                run_samtools view -H --no-PG "${tmp}/det.${cal}.1.unal.bam" \
                    | grep $'^@CO\tpipeline: '
            ) \
            <(
                run_samtools view -H --no-PG "${tmp}/det.${cal}.3.unal.bam" \
                    | grep $'^@CO\tpipeline: '
            )
    then
        record_pass "${cal}: the side file is the same at 1 and 3 threads"
    else
        record_fail "${cal}: the side file differs between 1 and 3 threads"
    fi
done


# Step 2 on its own, with '--req_flg', on a SAM file with a #75 template
# ('t75'), the same with a supplementary record ('supp'), a complete pair
# ('pair'), an aligned R1 with no R2 record ('lone'), and a same-strand pair
# flagged proper ('ffpp'). Step 1's split is made as Step 1 makes it. Under
# 'separate', both split templates move whole into Set A; 'lone', which has no
# record in Set A, is dropped by the template filter; and 'ffpp' is kept, since
# '-f 2' runs before 'fixmate' unflags it, all as under 'remove'.
for mod in remove separate; do
    stm="${tmp}/step2.${mod}"
    if [[ "${mod}" == "separate" ]]; then
        run_samtools view -b -F 12 -U "${stm}.setA.bam" -o "${stm}.bam" \
            "${sam_spl}"
    else
        run_samtools view -b -F 12 -f 2 -o "${stm}.bam" "${sam_spl}"
    fi

    rc=0
    out="$(
        run_func _sort_qname_bam \
            align_fastqs false 1 bwa 0 true "${mod}" pe \
            "${stm}.bam" "${stm}.qnam.bam"
    )" || rc=$?

    got="$(
        run_samtools view "${stm}.qnam.bam" | cut -f 1 | sort | uniq -c \
            | awk '{ printf "%s:%s ", $2, $1 }'
    )"
    if [[ "${mod}" == "separate" ]]; then
        exp="ffpp:2 pair:2 supp:3 t75:2 "
        nam="$(tr '\n' ' ' < "${stm}.unal.names.txt" 2> /dev/null || true)"
    else
        exp="ffpp:2 pair:2 "
        nam="supp t75 "
    fi

    if [[ "${rc}" -eq 0 && "${got}" == "${exp}" && "${nam}" == "supp t75 " ]]
    then
        record_pass "Step 2 alone under '${mod}' keeps '${exp% }'"
    else
        printf '%s\n' "${out}" > "${dir_log}/unal_mode_step2_${mod}.log"
        record_fail \
            "Step 2 alone under '${mod}' kept '${got}', expected '${exp}'" \
            "(exit ${rc}; Set A names '${nam}')"
    fi
done


# The awk BED converter drops unaligned reads and their mates, so a
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
