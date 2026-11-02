#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_filter_alignments_record.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="filter-alignments pipeline record"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Define fixture and output paths. The input already has the record line that
# 'align_fastqs' writes, and the filter's record must come after it
# ('OUTPUT.ALIGNMENT.RECORD').
dir_fx="${ROOT_REPO}/tests/fixtures/filter_alignments"
ref_fa="${dir_fx}/reference/filter_sc_sp.fa"

in_bam="${dir_fx}/bam/prior_record.bam"
in_cram="${dir_fx}/cram/prior_record.cram"

tmp="${TEST_DIR_TMP}/filter_alignments_record"
dir_log="${TEST_DIR_LOG}/filter_alignments"
chr_sc="I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI"
lin_aln=$'@CO\tpipeline: chain=align_fastqs'


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${tmp}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

# CRAM writing reads the reference's '.fai', which 'make.sh' builds; without
# it, Samtools would write one beside the fixture.
require_files_nonempty \
    "${in_bam}" \
    "${in_cram}" \
    "${ref_fa}" \
    "${ref_fa}.fai" || {
    finish
    exit $?
}


# Run one helper.
function run_filter() {
    # shellcheck disable=SC2016  # Expand in the child shell, not this one.
    "${TEST_BASH}" -c '
        # shellcheck source=lib/bash/core/source_helpers.sh
        source "${1}/lib/bash/core/source_helpers.sh"

        source_helpers "${1}/lib/bash" \
            core/check_args \
            core/check_inputs \
            core/format_outputs \
            workflows/filter_alignment
        shift 1
        "$@"
    ' _ "${ROOT_REPO}" "$@"
}


# Print an alignment file's '@CO' lines; reading adds no '@PG' line.
function print_co() {
    run_samtools view -H --no-PG -T "${ref_fa}" "${1}" | grep '^@CO' || true
}


# Print the content lines of a record.
function print_cnt() {
    grep $'^@CO\tpipeline: ' <<< "${1}" || true
}


# Print the run lines of a record.
function print_run() {
    grep $'^@CO\tpipeline-run: ' <<< "${1}" || true
}


# Print '@CO' lines from 'key=value' arguments under one prefix.
function print_lines() {
    local pfx="${1}"
    shift 1
    printf "@CO\t${pfx}: %s\n" "$@"
}


# Print the expected path lines from 'key path' pairs.
function print_pth_exp() {
    while (( $# > 1 )); do
        printf '@CO\tpipeline-run: --%s=%s\n' "${1}" "${2}"
        shift 2
    done
}


# Add one row, as label|helper|input|output|options|expected content lines; the
# content lines are the arguments after the fifth, joined with ';'.
function add_row() {
    local lbl="${1}"
    local fnc="${2}"
    local fil="${3}"
    local out="${4}"
    local opt="${5}"
    local IFS=';'

    shift 5
    rows+=( "${lbl}|${fnc}|${fil}|${out}|${opt}|$*" )
}


rows=()

add_row \
    "S. cerevisiae, BAM to BAM" \
    filter_alignment_sc \
    "${in_bam}" \
    sc.bam \
    "" \
    --retain=sc \
    --mito=false \
    "chromosomes=${chr_sc}" \
    --fil_in=prior_record.bam

add_row \
    "S. cerevisiae, BAM to CRAM, '--mito'" \
    filter_alignment_sc \
    "${in_bam}" \
    sc_mito.cram \
    --mito \
    --retain=sc \
    --mito=true \
    "chromosomes=${chr_sc} Mito" \
    --fil_in=prior_record.bam \
    --ref_fa=filter_sc_sp.fa

add_row \
    "S. pombe, CRAM to BAM, '--tg --mtr'" \
    filter_alignment_sp \
    "${in_cram}" \
    sp.bam \
    "--tg --mtr" \
    --retain=sp \
    --mito=false \
    --tg=true \
    --mtr=true \
    "chromosomes=SP_II_TG SP_I SP_II SP_III SP_MTR" \
    --fil_in=prior_record.cram \
    --ref_fa=filter_sc_sp.fa

for idx in "${!rows[@]}"; do
    IFS='|' read -r lbl fnc fil_in nam opt lin <<< "${rows[idx]}"
    read -r -a args <<< "${opt}"
    IFS=';' read -r -a arr_lin <<< "chain=filter_alignments;${lin}"
    fil_out="${tmp}/out_${idx}/${nam}"
    ref=()
    pth=( fil_in "${fil_in}" )
    if [[ "${fil_in}" == *.cram || "${nam}" == *.cram ]]; then
        ref=( --ref_fa "${ref_fa}" )
        pth+=( ref_fa "${ref_fa}" )
    fi
    pth+=( fil_out "${fil_out}" )
    mkdir -p "${tmp}/out_${idx}"

    rc=0
    run_filter "${fnc}" \
        --threads 1 \
        --fil_in "${fil_in}" \
        --fil_out "${fil_out}" \
        "${ref[@]}" \
        "${args[@]}" \
        > "${tmp}/out_${idx}/run.log" 2>&1 || rc=$?
    if (( rc != 0 )) || [[ ! -s "${fil_out}" ]]; then
        record_fail "${lbl}: run failed (exit ${rc})"
        continue
    fi

    assert_alignment_record \
        "${fil_out}" "${ref_fa}" filter_alignments 1 \
        "$(print_lines pipeline "${arr_lin[@]}")" \
        "$(print_pth_exp "${pth[@]}")" \
        "" "${lbl}" \
        || true

    # The earlier record is kept, ahead of the filter's.
    co="$(print_co "${fil_out}")"
    if [[
        "${co}" == "${lin_aln}"$'\n'*$'@CO\tpipeline: chain=filter_alignments'*
    ]]; then
        record_pass "${lbl}: the input's earlier record is kept first"
    else
        record_fail "${lbl}: the input's earlier record is lost or misplaced"
    fi

    # A direct call to the library gives the library function alone.
    if [[
        "${co}" == *$'pipeline-run: calls=filter_alignment.sh::'"${fnc}"$'\n'*
    ]]; then
        record_pass "${lbl}: the call path is the library alone"
    else
        record_fail "${lbl}: the call path is wrong for a direct call"
    fi
done


# Filtering a filtered file appends a second record; the shared check reads the
# last one.
fil_two="${tmp}/twice.sc.bam"
run_filter filter_alignment_sc \
    --threads 1 \
    --fil_in "${tmp}/out_0/sc.bam" \
    --fil_out "${fil_two}" \
    > "${tmp}/twice.log" 2>&1 || true

if [[
    "$(
        print_co "${fil_two}" \
            | grep -c $'^@CO\tpipeline: chain=filter_alignments$'
    )" == 2
]]; then
    record_pass "filtering twice keeps both filter records"
else
    record_fail "filtering twice lost or merged a filter record"
fi


# Two thread counts give the same content lines and different run lines.
for thr in 1 2; do
    mkdir -p "${tmp}/thr_${thr}"
    run_filter filter_alignment_sc \
        --threads "${thr}" \
        --fil_in "${in_bam}" \
        --fil_out "${tmp}/thr_${thr}/x.sc.bam" \
        > "${tmp}/thr_${thr}/run.log" 2>&1 || true
done

co_1="$(print_co "${tmp}/thr_1/x.sc.bam")"
co_2="$(print_co "${tmp}/thr_2/x.sc.bam")"
if [[
    -n "${co_1}"
        && "$(print_cnt "${co_1}")" == "$(print_cnt "${co_2}")"
        && "$(print_run "${co_1}")" != "$(print_run "${co_2}")"
        && "${co_2}" == *'pipeline-run: --threads=2'*
]]; then
    record_pass "two thread counts: same content lines, different run lines"
else
    record_fail "two thread counts: content lines differ or run lines match"
fi


# Through the drivers: 'execute' runs 'submit', which calls the library, so the
# call path names all three, and the record names the requested environment.
dir_drv="${tmp}/driver"
mkdir -p "${dir_drv}/out" "${dir_drv}/logs"
rc=0
"${TEST_BASH}" "${ROOT_REPO}/bin/execute_filter_alignments.sh" \
    --env_nam "${env_nam}" \
    --threads 1 \
    --csv_fil_in "${in_bam}" \
    --dir_out "${dir_drv}/out" \
    --retain sc \
    --dir_eo "${dir_drv}/logs" \
    --nam_job test_filter_record \
    --max_job 1 \
    > "${dir_drv}/run.log" 2>&1 || rc=$?

fil_drv="${dir_drv}/out/prior_record.sc.bam"
if (( rc != 0 )) || [[ ! -s "${fil_drv}" ]]; then
    record_fail "driver run failed (exit ${rc}); see '${dir_drv}/run.log'"
else
    assert_alignment_record \
        "${fil_drv}" "${ref_fa}" filter_alignments 1 \
        "$(
            print_lines pipeline \
                chain=filter_alignments --retain=sc --mito=false \
                "chromosomes=${chr_sc}" --fil_in=prior_record.bam
        )" \
        "$(print_pth_exp fil_in "${in_bam}" fil_out "${fil_drv}")" \
        "" "driver" \
        || true

    co="$(print_co "${fil_drv}")"
    if [[
        "${co}" == *$'pipeline-run: calls=execute_filter_alignments.sh > '*
            && "${co}" == *$'submit_filter_alignments.sh > '*
            && "${co}" == *$'filter_alignment.sh::filter_alignment_sc\n'*
            && "${co}" == *$'pipeline-run: env='"${env_nam}"$'\n'*
    ]]; then
        record_pass "driver: the call path and environment are recorded"
    else
        record_fail "driver: the call path or environment line is wrong"
        grep 'pipeline-run: ' <<< "${co}" >&2 || true
    fi
fi


# Through an alias: 'execute_filter_bams.sh' hands over to
# 'execute_filter_alignments.sh', and the call path names it first.
dir_als="${tmp}/alias"
mkdir -p "${dir_als}/out" "${dir_als}/logs"
rc=0
"${TEST_BASH}" "${ROOT_REPO}/bin/execute_filter_bams.sh" \
    --env_nam "${env_nam}" \
    --threads 1 \
    --csv_fil_in "${in_bam}" \
    --dir_out "${dir_als}/out" \
    --retain sc \
    --dir_eo "${dir_als}/logs" \
    --nam_job test_filter_record_alias \
    --max_job 1 \
    > "${dir_als}/run.log" 2>&1 || rc=$?

exp_als="execute_filter_bams.sh > execute_filter_alignments.sh >"
exp_als+=" submit_filter_alignments.sh >"
exp_als+=" filter_alignment.sh::filter_alignment_sc"
if [[
    "${rc}" -eq 0
        && "$(print_co "${dir_als}/out/prior_record.sc.bam")" \
            == *$'pipeline-run: calls='"${exp_als}"$'\n'*
]]; then
    record_pass "alias: the call path names the alias first"
else
    record_fail "alias: the call path is wrong (exit ${rc})"
fi


finish
