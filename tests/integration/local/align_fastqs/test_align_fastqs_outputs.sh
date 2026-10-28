#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_align_fastqs_outputs.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="align-fastqs shipped files"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Define fixture and output paths.
#
# Work files that a later step rewrites are written at compression level 1 and
# are not indexed; the files that ship are written at the default level, and
# Step 5 (or Step 6, for CRAM) indexes each coordinate-sorted one.
dir_fx="${ROOT_REPO}/tests/fixtures/align_fastqs"
idx_bt2="${dir_fx}/bowtie2/pairs"
idx_bwa="${dir_fx}/bwa/bridge.fa"
fq_pe_1="${dir_fx}/fastq/pe/pairs_pe_R1.atria.fastq.gz"
fq_pe_2="${dir_fx}/fastq/pe/pairs_pe_R2.atria.fastq.gz"
fq_se="${dir_fx}/fastq/se/pairs_se.atria.fastq.gz"
fq_brg_1="${dir_fx}/fastq/pe/bridge_pe_R1.atria.fastq.gz"
fq_brg_2="${dir_fx}/fastq/pe/bridge_pe_R2.atria.fastq.gz"
idx_fil="${dir_fx}/bwa/pairs.fa"
fq_fil_1="${dir_fx}/fastq/pe/pairs_fill_pe_R1.atria.fastq.gz"
fq_fil_2="${dir_fx}/fastq/pe/pairs_fill_pe_R2.atria.fastq.gz"
ref_fa="${dir_fx}/reference/pairs.fa"

tmp="${TEST_DIR_TMP}/align_fastqs_outputs"


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${tmp}"

require_env_project env_nam || {
    finish
    exit $?
}

# CRAM writing reads the reference's '.fai', which 'make.sh' builds; without
# it, Samtools would write one beside the fixture.
require_files_nonempty \
    "${ref_fa}" \
    "${ref_fa}.fai" \
    "${idx_bt2}.rev.2.bt2" \
    "${idx_bwa}.bwt" \
    "${fq_pe_1}" \
    "${fq_pe_2}" \
    "${fq_se}" \
    "${fq_brg_1}" \
    "${fq_brg_2}" \
    "${idx_fil}.bwt" \
    "${fq_fil_1}" \
    "${fq_fil_2}" || {
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


# Run Step 5 alone on a file, with the record held in the arguments.
function run_record() {
    # shellcheck disable=SC2016  # Expand in the child shell, not this one.
    "${TEST_BASH}" -c '
        # shellcheck source=lib/bash/core/source_helpers.sh
        source "${1}/lib/bash/core/source_helpers.sh"

        source_helpers "${1}/lib/bash" workflows/align_fastqs
        arr_cnt=( "${4}" )
        arr_run=( "${5}" )
        _write_align_record align_fastqs false 1 "${2}" arr_cnt arr_run "${3}"
    ' _ "${ROOT_REPO}" "$@" \
        $'@CO\tpipeline: chain=align_fastqs' \
        $'@CO\tpipeline-run: --threads=1'
}


# Print 'same <n>' or 'differ <n>' for the BGZF blocks holding two BAM files'
# records; see 'tests/support/compare_bam_blocks.py'.
function compare_record_blocks() {
    "${TEST_PYTHON}" "${ROOT_REPO}/tests/support/compare_bam_blocks.py" "$@"
}


# Check every file one call wrote: the set of names, the compression of each
# BAM file, and each index against a fresh one.
function check_outputs() {
    local lbl="${1}"
    local dir="${2}"
    local exp="${3}"
    local fil got res src

    got="$(cd "${dir}" && find . -type f ! -name '*.log' | LC_ALL=C sort)"
    if [[ "${got}" == "${exp}" ]]; then
        record_pass "${lbl}: writes exactly the expected files"
    else
        record_fail "${lbl}: writes other files than expected"
        diff <(printf '%s\n' "${exp}") <(printf '%s\n' "${got}") >&2 || true
    fi

    # A file without records has no record block to judge, such as an improper
    # file without '--req_flg'.
    for fil in "${dir}"/*.bam; do
        [[ -e "${fil}" ]] || continue
        run_samtools view -b --no-PG -o "${tmp}/default.bam" "${fil}"
        res="$(compare_record_blocks "${fil}" "${tmp}/default.bam")"

        if [[ "${res}" == "same "* && "${res}" != "same 0" ]]; then
            record_pass "${lbl}: '${fil##*/}' is at the default level"
        elif [[
            "${res}" == "same 0" && "$(run_samtools view -c "${fil}")" == 0
        ]]; then
            record_pass "${lbl}: '${fil##*/}' holds no records"
        else
            record_fail "${lbl}: '${fil##*/}' is not at the default level" \
                "(${res})"
        fi
    done

    for fil in "${dir}"/*.bai "${dir}"/*.crai; do
        [[ -e "${fil}" ]] || continue
        src="${fil%.*}"
        run_samtools index -o "${tmp}/fresh.${fil##*.}" "${src}"

        if cmp -s "${fil}" "${tmp}/fresh.${fil##*.}"; then
            record_pass "${lbl}: '${fil##*/}' equals a fresh index"
        else
            record_fail "${lbl}: '${fil##*/}' differs from a fresh index"
        fi
    done
}


# label|output extension|helper arguments|files expected beside 'x'
#
# With '--unal_mode separate', the BWA-backtrack call on 'pairs_fill' fills
# every diagnostic file: '--req_flg' rejects the improper pairs and MAPQ 30 a
# repeat pair.
arr_calls=(
    "Bowtie 2 PE 'remove' BAM|bam|--aligner bowtie2 --index ${idx_bt2}\
 --fq_1 ${fq_pe_1} --fq_2 ${fq_pe_2} --qname\
|bam bam.bai markdup.txt.gz qnam.bam"
    "Bowtie 2 PE 'separate' BAM|bam|--aligner bowtie2 --index ${idx_bt2}\
 --fq_1 ${fq_pe_1} --fq_2 ${fq_pe_2} --qname --unal_mode separate\
|bam bam.bai markdup.txt.gz qnam.bam unal.bam unal.bam.bai unal.qnam.bam\
 improp.bam improp.bam.bai improp.qnam.bam mapq.bam mapq.bam.bai\
 mapq.qnam.bam all.markdup.txt.gz diag.tsv.gz"
    "Bowtie 2 PE 'remove' CRAM|cram|--aligner bowtie2 --index ${idx_bt2}\
 --fq_1 ${fq_pe_1} --fq_2 ${fq_pe_2} --qname --ref_fa ${ref_fa}\
|cram cram.crai markdup.txt.gz qnam.cram"
    "Bowtie 2 PE 'separate' CRAM|cram|--aligner bowtie2 --index ${idx_bt2}\
 --fq_1 ${fq_pe_1} --fq_2 ${fq_pe_2} --qname --ref_fa ${ref_fa}\
 --unal_mode separate\
|cram cram.crai markdup.txt.gz qnam.cram unal.cram unal.cram.crai\
 unal.qnam.cram improp.cram improp.cram.crai improp.qnam.cram mapq.cram\
 mapq.cram.crai mapq.qnam.cram all.markdup.txt.gz diag.tsv.gz"
    "BWA-backtrack PE 'separate' BAM|bam|--aligner bwa --bwa_alg aln\
 --mapq 0 --index ${idx_bwa} --fq_1 ${fq_brg_1} --fq_2 ${fq_brg_2}\
 --unal_mode separate\
|bam bam.bai markdup.txt.gz unal.bam unal.bam.bai improp.bam improp.bam.bai\
 mapq.bam mapq.bam.bai all.markdup.txt.gz diag.tsv.gz"
    "BWA-backtrack PE 'separate' '--req_flg' BAM|bam|--aligner bwa\
 --bwa_alg aln --req_flg --mapq 30 --index ${idx_fil} --fq_1 ${fq_fil_1}\
 --fq_2 ${fq_fil_2} --qname --unal_mode separate\
|bam bam.bai markdup.txt.gz qnam.bam unal.bam unal.bam.bai unal.qnam.bam\
 improp.bam improp.bam.bai improp.qnam.bam mapq.bam mapq.bam.bai\
 mapq.qnam.bam all.markdup.txt.gz diag.tsv.gz"
    "Bowtie 2 SE 'separate' BAM|bam|--aligner bowtie2 --index ${idx_bt2}\
 --fq_1 ${fq_se} --unal_mode separate\
|bam bam.bai markdup.txt.gz unal.bam unal.bam.bai mapq.bam mapq.bam.bai\
 all.markdup.txt.gz diag.tsv.gz"
)

for idx in "${!arr_calls[@]}"; do
    IFS='|' read -r lbl ext arg fil <<< "${arr_calls[idx]}"
    read -r -a arr_args <<< "${arg}"
    read -r -a arr_fil <<< "${fil}"
    dir="${tmp}/call_${idx}"
    mkdir -p "${dir}"

    if ! \
        run_helper \
            --threads 2 \
            "${arr_args[@]}" \
            --fil_out "${dir}/x.${ext}" \
            > "${dir}/run.log" 2>&1
    then
        record_fail "${lbl}: the run failed; see '${dir}/run.log'"
        continue
    fi

    exp="$(printf './x.%s\n' "${arr_fil[@]}" | LC_ALL=C sort)"
    check_outputs "${lbl}" "${dir}" "${exp}"

    # In paired-end data, 'fixmate -c -m' writes each record's mate score
    # ('ms') and mate CIGAR ('MC'), and one template CIGAR ('ct') per pair; the
    # main file holds aligned pairs only, so every record carries 'ms' and
    # 'MC', and half of them 'ct'.
    if [[ "${lbl}" == *" PE "* ]]; then
        fil="${dir}/x.${ext}"
        n_rec="$(run_samtools view -c --reference "${ref_fa}" "${fil}")"
        n_bad="$(
            run_samtools view --reference "${ref_fa}" "${fil}" \
                | awk -F '\t' '
                    {
                        tag = ""
                        for (i = 12; i <= NF; i++) {
                            tag = tag " " substr($i, 1, 5)
                        }
                        if (tag !~ / ms:i:/ || tag !~ / MC:Z:/) n++
                        if (tag ~ / ct:Z:/) n_ct++
                    }
                    END { print n + (2 * n_ct != NR) }
                '
        )"

        if (( n_rec > 0 && n_bad == 0 )); then
            record_pass "${lbl}: ${n_rec} records carry 'fixmate' tags"
        else
            record_fail "${lbl}: of ${n_rec} records, ${n_bad} lack 'ms' or" \
                "'MC' (1 more if 'ct' is not on one record per pair)"
        fi
    fi
done


# Step 5 indexes a file it skips because it already carries a record, as no
# earlier step leaves an index; for CRAM output ('idx' false), it indexes none.
dir="${tmp}/skip"
mkdir -p "${dir}"

for idx in true false; do
    cp "${tmp}/call_0/x.bam" "${dir}/skip_${idx}.bam"
    rc=0
    out="$(run_record "${dir}/skip_${idx}.bam" "${idx}" 2>&1)" || rc=$?
    lbl="Step 5 on a file with a record, 'idx' ${idx}"

    if [[
        "${rc}" -ne 0
        || "${out}" != *"already carries a pipeline record"*
    ]]; then
        record_fail "${lbl}: exit ${rc}, or no skip message: ${out}"
    elif [[ "${idx}" == "false" ]]; then
        if [[ ! -e "${dir}/skip_${idx}.bam.bai" ]]; then
            record_pass "${lbl}: skipped and not indexed"
        else
            record_fail "${lbl}: indexed the file"
        fi
    else
        run_samtools index -o "${tmp}/fresh.bai" "${dir}/skip_${idx}.bam"

        if cmp -s "${dir}/skip_${idx}.bam.bai" "${tmp}/fresh.bai"; then
            record_pass "${lbl}: skipped, then indexed"
        else
            record_fail "${lbl}: skipped but left no correct index"
        fi
    fi
done


finish
