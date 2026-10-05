#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_filter_alignments_reference_order.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="filter-alignments reference order"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Define fixture and output paths. The input header lists 'SP_I' before the S.
# cerevisiae chromosomes, so a filter that drops references without renumbering
# them would move reads to the wrong chromosome.
dir_fx="${ROOT_REPO}/tests/fixtures/filter_alignments"
in_sam="${dir_fx}/sam/filter_sc_sp.sam"
ref_fa="${dir_fx}/reference/filter_sc_sp.fa"

tmp="${TEST_DIR_TMP}/filter_alignments_reference_order"
dir_log="${TEST_DIR_LOG}/filter_alignments"
sam_ord="${tmp}/reordered.sam"
in_bam="${tmp}/reordered.bam"


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${tmp}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty "${in_sam}" "${ref_fa}" "${ref_fa}.fai" || {
    finish
    exit $?
}

# Move the 'SP_I' '@SQ' line ahead of every other '@SQ' line.
awk '
    /^@SQ/ && /\tSN:SP_I\t/ { sq_first = $0; next }
    /^@SQ/ { sq[++n_sq] = $0; next }
    /^@/ { hd[++n_hd] = $0; next }
    { body[++n_bd] = $0 }
    END {
        for (i = 1; i <= n_hd; i++) if (hd[i] ~ /^@HD/) print hd[i]
        print sq_first
        for (i = 1; i <= n_sq; i++) print sq[i]
        for (i = 1; i <= n_hd; i++) if (hd[i] !~ /^@HD/) print hd[i]
        for (i = 1; i <= n_bd; i++) print body[i]
    }
' "${in_sam}" > "${sam_ord}"

if [[ "$(grep -m1 '^@SQ' "${sam_ord}")" != *$'\tSN:SP_I\t'* ]]; then
    record_fail "reordered fixture does not list 'SP_I' first"
    finish
    exit $?
fi

if ! \
    build_filter_alignments_fixture_bam \
        "${sam_ord}" \
        "${in_bam}" \
        "${dir_log}/reference_order_prepare.log" \
        "filter-alignments reordered BAM fixture"
then
    finish
    exit $?
fi


# Run one helper on the reordered input.
function run_filter() {
    # Expand in the child shell, not this one.
    # shellcheck disable=SC2016
    "${TEST_BASH}" -c '
        # shellcheck disable=SC1090
        source "${1}/lib/bash/core/source_helpers.sh"

        source_helpers "${1}/lib/bash" \
            core/check_args \
            core/check_inputs \
            core/format_outputs \
            workflows/filter_alignment
        shift 1
        "$@"
    ' _ "${ROOT_REPO}" "$@" 2>&1
}


# Print each record's read name and chromosome.
function print_qname_rname() {
    run_samtools view -T "${ref_fa}" "${1}" | cut -f 1,3 | sort
}


# Every read kept must sit on the chromosome it had in the input, and the
# expected reads must be kept.
pairs_in="$(print_qname_rname "${in_bam}")"

sp_all="r_sp_i r_sp_mito r_sp_mtr r_sp_tg"
rows=(
    "filter_alignment_sc|bam|--mito|r_sc_I r_sc_mito"
    "filter_alignment_sc|cram|--mito|r_sc_I r_sc_mito"
    "filter_alignment_sc|bam||r_sc_I"
    "filter_alignment_sp|bam|--mito --tg --mtr|${sp_all}"
)

for idx in "${!rows[@]}"; do
    IFS='|' read -r fnc ext arg exp <<< "${rows[idx]}"
    read -r -a args <<< "${arg}"
    read -r -a arr_exp <<< "${exp}"
    fil_out="${tmp}/out_${idx}.${ext}"
    lbl="${fnc} to ${ext^^} ${arg:-without options}"

    rc=0
    out="$(
        run_filter "${fnc}" \
            --threads 1 \
            --fil_in "${in_bam}" \
            --fil_out "${fil_out}" \
            --ref_fa "${ref_fa}" \
            "${args[@]}"
    )" || rc=$?

    if [[ "${rc}" -ne 0 || ! -s "${fil_out}" ]]; then
        printf '%s\n' "${out}" > "${dir_log}/reference_order_${idx}.log"
        record_fail "${lbl} failed with a reordered reference (exit ${rc})"
        continue
    fi

    pairs_out="$(print_qname_rname "${fil_out}")"
    ok=true

    # Each kept read keeps its chromosome.
    while IFS= read -r pair; do
        if ! grep -Fxq -- "${pair}" <<< "${pairs_in}"; then ok=false; fi
    done <<< "${pairs_out}"

    # Exactly the expected reads are kept.
    if [[
        "$(cut -f 1 <<< "${pairs_out}" | tr '\n' ' ')" \
            != "$(printf '%s\n' "${arr_exp[@]}" | sort | tr '\n' ' ')"
    ]]; then
        ok=false
    fi

    if [[ "${ok}" == "true" ]]; then
        record_pass "${lbl} keeps each read on its chromosome"
    else
        record_fail \
            "${lbl} moved or lost reads with a reordered reference:" \
            "$(tr '\n' ' ' <<< "${pairs_out}")"
    fi
done


finish
