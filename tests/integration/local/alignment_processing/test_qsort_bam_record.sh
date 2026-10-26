#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_qsort_bam_record.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="qsort BAM/CRAM pipeline record"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Define fixture and output paths. The input carries an earlier record, as an
# 'align_fastqs' output would, so the sort's record must follow it
# ('OUTPUT.ALIGNMENT.RECORD').
dir_fx="${ROOT_REPO}/tests/fixtures/calculate_scaling_factor"
fx_bam="${dir_fx}/bam/pe/IP_WT_G1_Hho1_6336.sc.bam"

ref_fa="${dir_fx}/reference/tiny.fa"

tmp="${TEST_DIR_TMP}/qsort_bam_record"
in_bam="${tmp}/in.bam"
in_cram="${tmp}/in.cram"
lin_aln=$'@CO\tpipeline: chain=align_fastqs'


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${tmp}"

require_env_project env_nam || {
    finish
    exit $?
}

# CRAM writing reads the reference's '.fai', which 'make.sh' builds; without
# it, Samtools would write one beside the fixture.
require_files_nonempty "${fx_bam}" "${ref_fa}" "${ref_fa}.fai" || {
    finish
    exit $?
}

# Give the input an earlier record, then make a CRAM copy of it.
{
    run_samtools view -H --no-PG "${fx_bam}"
    printf '%s\n' "${lin_aln}"
} > "${tmp}/in.header.sam"
run_samtools reheader "${tmp}/in.header.sam" "${fx_bam}" > "${in_bam}"
run_samtools view -C -T "${ref_fa}" -o "${in_cram}" "${in_bam}"


# Run the helper.
function run_qsort() {
    # shellcheck disable=SC2016  # Expand in the child shell, not this one.
    "${TEST_BASH}" -c '
        # shellcheck source=lib/bash/core/source_helpers.sh
        source "${1}/lib/bash/core/source_helpers.sh"

        source_helpers "${1}/lib/bash" workflows/process_alignments
        shift 1
        qsort_file_alignments "$@"
    ' _ "${ROOT_REPO}" "$@"
}


# Print an alignment file's '@CO' lines; reading adds no '@PG' line.
function print_co() {
    run_samtools view -H --no-PG -T "${ref_fa}" "${1}" | grep '^@CO' || true
}


# Print the expected path lines from 'key path' pairs.
function print_pth_exp() {
    while (( $# > 1 )); do
        printf '@CO\tpipeline-run: --%s=%s\n' "${1}" "${2}"
        shift 2
    done
}


# label|input|output|threads
rows=(
    "BAM|${in_bam}|${tmp}/out_bam/x.qnam.bam|1"
    "BAM, two threads|${in_bam}|${tmp}/out_bam_2/x.qnam.bam|2"
    "CRAM|${in_cram}|${tmp}/out_cram/x.qnam.cram|1"
)

for row in "${rows[@]}"; do
    IFS='|' read -r lbl fil_in fil_out thr <<< "${row}"
    dir_out="$(dirname "${fil_out}")"
    mkdir -p "${dir_out}"

    ref=""
    exp_cnt=$'@CO\tpipeline: chain=qsort_bam\n@CO\tpipeline: --fil_in='
    exp_cnt+="${fil_in##*/}"
    pth=( fil_in "${fil_in}" )
    if [[ "${fil_in}" == *.cram ]]; then
        ref="${ref_fa}"
        exp_cnt+=$'\n@CO\tpipeline: --ref_fa=tiny.fa'
        pth+=( ref_fa "${ref_fa}" )
    fi
    pth+=( fil_out "${fil_out}" )

    rc=0
    run_qsort \
        "${thr}" "${fil_in}" "${fil_out}" "${ref}" \
        "${dir_out}/qsort.out" "${dir_out}/qsort.err" \
        > "${dir_out}/run.log" 2>&1 || rc=$?
    if (( rc != 0 )) || [[ ! -s "${fil_out}" ]]; then
        record_fail "${lbl}: sort failed (exit ${rc})"
        continue
    fi

    assert_alignment_record \
        "${fil_out}" "${ref_fa}" qsort_bam "${thr}" \
        "${exp_cnt}" "$(print_pth_exp "${pth[@]}")" "" "${lbl}" \
        || true

    # The earlier record is kept, ahead of the sort's.
    co="$(print_co "${fil_out}")"
    if [[
        "${co}" == "${lin_aln}"$'\n@CO\tpipeline: chain=qsort_bam'*
    ]]; then
        record_pass "${lbl}: the input's earlier record is kept first"
    else
        record_fail "${lbl}: the input's earlier record is lost or misplaced"
    fi

    # A direct call to the library gives the library function alone.
    if [[
        "${co}" == *$'calls=process_alignments.sh::qsort_file_alignments\n'*
    ]]; then
        record_pass "${lbl}: the call path is the library alone"
    else
        record_fail "${lbl}: the call path is wrong for a direct call"
    fi

    # The sort still holds every record, sorted by queryname.
    if [[
        "$(run_samtools view -c -T "${ref_fa}" "${fil_out}")" \
            == "$(run_samtools view -c -T "${ref_fa}" "${fil_in}")"
        && "$(run_samtools view -H -T "${ref_fa}" "${fil_out}")" \
            == *$'\tSO:queryname'*
    ]]; then
        record_pass "${lbl}: every record is kept, sorted by queryname"
    else
        record_fail "${lbl}: records lost, or the file is not name-sorted"
    fi

    if [[
        -n "$(
            find "${dir_out}" \
                -name '*.header.sam' \
                -o -name '*.record.bam' \
                -o -name '*.tmp.*'
        )"
    ]]; then
        record_fail "${lbl}: the record's temporary files are left behind"
    else
        record_pass "${lbl}: the record's temporary files are removed"
    fi
done


# Two thread counts give the same content lines and different run lines.
co_1="$(print_co "${tmp}/out_bam/x.qnam.bam")"
co_2="$(print_co "${tmp}/out_bam_2/x.qnam.bam")"
if [[
    -n "${co_1}"
        && "$(grep 'pipeline: ' <<< "${co_1}")" \
            == "$(grep 'pipeline: ' <<< "${co_2}")"
        && "$(grep 'pipeline-run: ' <<< "${co_1}")" \
            != "$(grep 'pipeline-run: ' <<< "${co_2}")"
]]; then
    record_pass "two thread counts: same content lines, different run lines"
else
    record_fail "two thread counts: content lines differ or run lines match"
fi


# Through the driver: 'submit' calls the library, so the call path names both,
# and the record names the requested environment.
dir_drv="${tmp}/driver"
mkdir -p "${dir_drv}/out" "${dir_drv}/logs"
rc=0
"${TEST_BASH}" "${ROOT_REPO}/bin/submit_qsort_bam.sh" \
    --env_nam "${env_nam}" \
    --dir_scr "${ROOT_REPO}/bin" \
    --threads 1 \
    --csv_fil_in "${in_bam}" \
    --dir_out "${dir_drv}/out" \
    --dir_eo "${dir_drv}/logs" \
    --nam_job test_qsort_record \
    > "${dir_drv}/run.log" 2>&1 || rc=$?

fil_drv="${dir_drv}/out/in.qnam.bam"
if (( rc != 0 )) || [[ ! -s "${fil_drv}" ]]; then
    record_fail "driver run failed (exit ${rc}); see '${dir_drv}/run.log'"
else
    assert_alignment_record \
        "${fil_drv}" "" qsort_bam 1 \
        $'@CO\tpipeline: chain=qsort_bam\n@CO\tpipeline: --fil_in=in.bam' \
        "$(print_pth_exp fil_in "${in_bam}" fil_out "${fil_drv}")" \
        "" "driver" \
        || true

    co="$(print_co "${fil_drv}")"
    if [[
        "${co}" == *$'pipeline-run: calls=submit_qsort_bam.sh > '*
            && "${co}" == *$'process_alignments.sh::qsort_file_alignments\n'*
            && "${co}" == *$'pipeline-run: env='"${env_nam}"$'\n'*
    ]]; then
        record_pass "driver: the call path and environment are recorded"
    else
        record_fail "driver: the call path or environment line is wrong"
        grep 'pipeline-run: ' <<< "${co}" >&2 || true
    fi
fi


finish
