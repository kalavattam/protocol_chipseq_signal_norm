#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_align_fastqs_spaced_paths.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="align-fastqs spaced paths"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Define fixture and output paths. Each arm runs on the fixture paths and on
# spaced links to them; both runs must write the same files and records.
dir_fx="${ROOT_REPO}/tests/fixtures/align_fastqs"

tmp="${TEST_DIR_TMP}/align_fastqs_spaced_paths"
dir_in="${tmp}/in dir"
dir_pln="${tmp}/plain"
dir_spc="${tmp}/out dir"


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${dir_in}/bt2 idx" "${dir_in}/bwa idx" "${dir_in}/bm2 idx" \
    "${dir_pln}" "${dir_spc}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty \
    "${dir_fx}/fastq/se/tiny_se.atria.fastq.gz" \
    "${dir_fx}/fastq/pe/tiny_pe_R1.atria.fastq.gz" \
    "${dir_fx}/fastq/pe/tiny_pe_R2.atria.fastq.gz" \
    "${dir_fx}/reference/tiny.fa" \
    "${dir_fx}/bowtie2/tiny.1.bt2" \
    "${dir_fx}/bwa/tiny.fa.bwt" \
    "${dir_fx}/bwa-mem2/tiny.fa.0123" || {
    finish
    exit $?
}

# Link the inputs and indexes under directories whose names hold spaces.
for fil in \
    "${dir_fx}"/fastq/se/tiny_se.atria.fastq.gz \
    "${dir_fx}"/fastq/pe/tiny_pe_R[12].atria.fastq.gz \
    "${dir_fx}"/reference/tiny.fa*
do
    ln -s "${fil}" "${dir_in}/$(basename "${fil}")"
done

for fil in "${dir_fx}"/bowtie2/tiny*; do
    ln -s "${fil}" "${dir_in}/bt2 idx/$(basename "${fil}")"
done

for fil in "${dir_fx}"/bwa/tiny*; do
    ln -s "${fil}" "${dir_in}/bwa idx/$(basename "${fil}")"
done

for fil in "${dir_fx}"/bwa-mem2/tiny*; do
    ln -s "${fil}" "${dir_in}/bm2 idx/$(basename "${fil}")"
done


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
    ' _ "${ROOT_REPO}" --threads 1 "$@"
}


# Whether two alignment files hold the same records ('@PG' lines name the
# paths, so they are left out).
function same_records() {
    cmp -s \
        <(run_samtools view --no-PG -T "${dir_fx}/reference/tiny.fa" "${1}") \
        <(run_samtools view --no-PG -T "${dir_fx}/reference/tiny.fa" "${2}")
}


# Compare the files under two output directories: the same names, and, for
# alignment files, the same records.
function check_same_outputs() {
    local lbl="${1}"
    local dir_a="${2}"
    local dir_b="${3}"
    local fil nam n_dif=0

    if [[
        "$(cd "${dir_a}" && find . -type f | LC_ALL=C sort)" \
            != "$(cd "${dir_b}" && find . -type f | LC_ALL=C sort)"
    ]]; then
        record_fail "${lbl}: the two runs wrote different files"
        return
    fi

    while IFS= read -r fil; do
        nam="${fil#"${dir_a}"/}"
        case "${nam}" in
            *.bam|*.cram)
                if ! same_records "${fil}" "${dir_b}/${nam}"; then
                    n_dif=$(( n_dif + 1 ))
                fi
                ;;
        esac
    done < <(find "${dir_a}" -type f | LC_ALL=C sort)

    if (( n_dif == 0 )); then
        record_pass "${lbl}: same files and records"
    else
        record_fail "${lbl}: ${n_dif} alignment file(s) differ"
    fi
}


# label|arguments on fixture paths|arguments on spaced paths|output stem
fq_se_p="${dir_fx}/fastq/se/tiny_se.atria.fastq.gz"
fq_1_p="${dir_fx}/fastq/pe/tiny_pe_R1.atria.fastq.gz"
fq_2_p="${dir_fx}/fastq/pe/tiny_pe_R2.atria.fastq.gz"
fq_se_s="${dir_in}/tiny_se.atria.fastq.gz"
fq_1_s="${dir_in}/tiny_pe_R1.atria.fastq.gz"
fq_2_s="${dir_in}/tiny_pe_R2.atria.fastq.gz"

arms=(
    "bowtie2 PE, '--req_flg', '--qname'|bt2|pe|--aligner bowtie2 --req_flg --qname|bam"
    "BWA-backtrack PE, '--mapq 20'|bwa|pe|--aligner bwa --bwa_alg aln --mapq 20|bam"
    "BWA-MEM2 SE|bm2|se|--aligner bwa-mem2|bam"
    "BWA-MEM PE, CRAM, '--qname'|bwa|pe|--aligner bwa --qname|cram"
)

for idx in "${!arms[@]}"; do
    IFS='|' read -r lbl ind typ arg ext <<< "${arms[idx]}"
    read -r -a args <<< "${arg}"

    case "${ind}" in
        bt2)
            ind_p="${dir_fx}/bowtie2/tiny"
            ind_s="${dir_in}/bt2 idx/tiny"
            ;;

        bwa)
            ind_p="${dir_fx}/bwa/tiny.fa"
            ind_s="${dir_in}/bwa idx/tiny.fa"
            ;;

        bm2)
            ind_p="${dir_fx}/bwa-mem2/tiny.fa"
            ind_s="${dir_in}/bm2 idx/tiny.fa"
            ;;
    esac

    if [[ "${typ}" == "pe" ]]; then
        fq_p=( --fq_1 "${fq_1_p}" --fq_2 "${fq_2_p}" )
        fq_s=( --fq_1 "${fq_1_s}" --fq_2 "${fq_2_s}" )
    else
        fq_p=( --fq_1 "${fq_se_p}" )
        fq_s=( --fq_1 "${fq_se_s}" )
    fi

    ref_p=()
    ref_s=()
    if [[ "${ext}" == "cram" ]]; then
        ref_p=( --ref_fa "${dir_fx}/reference/tiny.fa" )
        ref_s=( --ref_fa "${dir_in}/tiny.fa" )
    fi

    mkdir -p "${dir_pln}/arm_${idx}" "${dir_spc}/arm ${idx}"
    rc_p=0
    rc_s=0

    run_helper "${args[@]}" --index "${ind_p}" "${ref_p[@]}" "${fq_p[@]}" \
        --fil_out "${dir_pln}/arm_${idx}/tiny.${ext}" > /dev/null 2>&1 \
        || rc_p=$?
    run_helper "${args[@]}" --index "${ind_s}" "${ref_s[@]}" "${fq_s[@]}" \
        --fil_out "${dir_spc}/arm ${idx}/tiny.${ext}" > /dev/null 2>&1 \
        || rc_s=$?

    if [[ "${rc_p}" -ne 0 || "${rc_s}" -ne 0 ]]; then
        record_fail "${lbl}: exit ${rc_p} on plain paths, ${rc_s} on spaced"
        continue
    fi

    check_same_outputs \
        "${lbl}" "${dir_pln}/arm_${idx}" "${dir_spc}/arm ${idx}"
done


# The drivers carry spaced paths too: 'execute' on a spaced input and output
# directory writes what the helper wrote above for the same arm.
dir_exe="${tmp}/exe dir"
mkdir -p "${dir_exe}/logs"

rc=0
"${TEST_BASH}" "${ROOT_REPO}/bin/execute_align_fastqs.sh" \
    --threads 1 \
    --csv_fil_in "${fq_1_s},${fq_2_s}" \
    --dir_out "${dir_exe}" \
    --dir_eo "${dir_exe}/logs" \
    --max_job 1 \
    --aligner bowtie2 \
    --index "${dir_in}/bt2 idx/tiny" \
    --req_flg \
    > /dev/null 2>&1 || rc=$?

if [[
    "${rc}" -eq 0
    && -s "${dir_exe}/tiny_pe.bam"
    && -s "${dir_exe}/tiny_pe.bam.bai"
]] && same_records "${dir_pln}/arm_0/tiny.bam" "${dir_exe}/tiny_pe.bam"; then
    record_pass "execute on spaced paths: same records as the helper"
else
    record_fail "execute on spaced paths: exit ${rc} or records differ"
fi


finish
