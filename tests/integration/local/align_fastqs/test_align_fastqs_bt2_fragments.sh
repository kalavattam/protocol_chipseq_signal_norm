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


# Define fixture and output paths. The fragment-length pairs on the repeat
# reference show which fragments each Bowtie 2 paired-end call keeps.
dir_fx="${ROOT_REPO}/tests/fixtures/align_fastqs"
idx_rpt="${dir_fx}/bowtie2/repeat"
fq_1="${dir_fx}/fastq/pe/frg_pe_R1.atria.fastq.gz"
fq_2="${dir_fx}/fastq/pe/frg_pe_R2.atria.fastq.gz"

tmp="${TEST_DIR_TMP}/align_fastqs_bt2_fragments"


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${tmp}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty \
    "${fq_1}" \
    "${fq_2}" \
    "${idx_rpt}.1.bt2" || {
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
