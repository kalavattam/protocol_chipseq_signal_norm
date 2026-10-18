#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_align_fastqs_record.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="align-fastqs pipeline record"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Define fixture and output paths. Every output, the '--qname' copy and CRAM
# included, must carry the pipeline record ('OUTPUT.ALIGNMENT.RECORD').
dir_fx="${ROOT_REPO}/tests/fixtures/align_fastqs"
idx_bt2="${dir_fx}/bowtie2/tiny"
idx_bwa="${dir_fx}/bwa/tiny.fa"
idx_bm2="${dir_fx}/bwa-mem2/tiny.fa"
in_pe_1="${dir_fx}/fastq/pe/tiny_pe_R1.atria.fastq.gz"
in_pe_2="${dir_fx}/fastq/pe/tiny_pe_R2.atria.fastq.gz"
in_se="${dir_fx}/fastq/se/tiny_se.atria.fastq.gz"

tmp="${TEST_DIR_TMP}/align_fastqs_record"
ref_fa="${tmp}/tiny.fa"
re_opt='^@CO'$'\t''pipeline: (optical_distance=2500|'
re_opt+='skipped=optical duplicates [(]not used; )'


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${tmp}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty \
    "${dir_fx}/reference/tiny.fa" \
    "${idx_bt2}.1.bt2" \
    "${idx_bwa}.bwt" \
    "${idx_bm2}.0123" \
    "${in_pe_1}" \
    "${in_pe_2}" \
    "${in_se}" || {
    finish
    exit $?
}

# CRAM writing indexes the reference, so use a copy outside the fixtures.
cp "${dir_fx}/reference/tiny.fa" "${ref_fa}"


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


# Print an alignment file's '@CO' lines; reading adds no '@PG' line.
function print_co() {
    run_samtools view -H --no-PG -T "${ref_fa}" "${1}" | grep '^@CO' || true
}


# Print the content lines, or the run lines, of a record.
function print_cnt() {
    grep $'^@CO\tpipeline: ' <<< "${1}" || true
}

function print_run() {
    grep $'^@CO\tpipeline-run: ' <<< "${1}" || true
}


# Print the content lines a row expects, leaving out the commit and optical
# lines, which are checked on their own.
function print_cnt_exp() {
    local lin

    for lin in "$@"; do
        printf '@CO\tpipeline: %s\n' "${lin}"
    done
}


# Check one output's record: the shared checks, with the optical line, whose
# reason depends on the reads, checked here.
function check_record() {
    local lbl="${1}"
    local fil="${2}"
    local exp_cnt="${3}"
    local exp_pth="${4}"
    local thr="${5}"

    assert_alignment_record \
        "${fil}" "${ref_fa}" align_fastqs "${thr}" \
        "${exp_cnt}" "${exp_pth}" "${re_opt}" "${lbl}" \
        || true

    # The optical line ends the content lines, right before the run lines.
    if
        print_co "${fil}" \
            | awk -v re="${re_opt}" '
                /^@CO\tpipeline-run: / { exit !(prv ~ re) }
                { prv = $0 }
            '
    then
        record_pass "${lbl}: optical line ends the content lines"
    else
        record_fail "${lbl}: optical line missing or out of place"
    fi
}


# Print the value that follows an option in an argument list.
function print_arg_val() {
    local opt="${1}"
    shift 1

    while (( $# > 1 )); do
        if [[ "${1}" == "${opt}" ]]; then
            printf '%s\n' "${2}"
            return 0
        fi
        shift 1
    done
}


# Print the expected path lines.
function print_pth_exp() {
    local key

    while (( $# > 1 )); do
        key="${1}"
        printf '@CO\tpipeline-run: --%s=%s\n' "${key}" "${2}"
        shift 2
    done
}


# Add one row, as label|output|arguments|expected content lines; the content
# lines are the arguments after the third, joined with ';'.
function add_row() {
    local lbl="${1}"
    local out="${2}"
    local arg="${3}"
    local IFS=';'

    shift 3
    rows+=( "${lbl}|${out}|${arg}|chain=align_fastqs;$*" )
}


# One row per aligner and layout.
a_pe="--fq_1 ${in_pe_1} --fq_2 ${in_pe_2}"
a_se="--fq_1 ${in_se}"
a_ref="--ref_fa ${ref_fa}"
k_pe=( --fq_1=tiny_pe_R1.atria.fastq.gz --fq_2=tiny_pe_R2.atria.fastq.gz )
k_se=( --fq_1=tiny_se.atria.fastq.gz "skipped=fixmate (single-end)" )
rows=()

add_row \
    "Bowtie 2 PE, '--qname'" \
    bt2_pe.bam \
    "--aligner bowtie2 --index ${idx_bt2} ${a_pe} --qname" \
    --aligner=bowtie2 \
    --bt2_mode=end-to-end \
    bt2_seed=0 \
    --bt2_X=1000 \
    --mapq=1 \
    --req_flg=false \
    --qname=true \
    --index=tiny \
    "${k_pe[@]}"

add_row \
    "Bowtie 2 PE, '--bt2_inherited'" \
    bt2_inh.bam \
    "--aligner bowtie2 --index ${idx_bt2} ${a_pe} --bt2_inherited --req_flg" \
    --aligner=bowtie2 \
    --bt2_mode=end-to-end \
    bt2_seed=0 \
    --bt2_inherited=true \
    --mapq=1 \
    --req_flg=true \
    --qname=false \
    --index=tiny \
    "${k_pe[@]}"

add_row \
    "Bowtie 2 SE, local" \
    bt2_se.bam \
    "--aligner bowtie2 --bt2_mode local --index ${idx_bt2} ${a_se}" \
    --aligner=bowtie2 \
    --bt2_mode=local \
    bt2_seed=0 \
    --mapq=1 \
    --qname=false \
    --index=tiny \
    "${k_se[@]}"

add_row \
    "BWA-MEM PE, '--mapq 0'" \
    mem_pe.bam \
    "--aligner bwa --index ${idx_bwa} ${a_pe} --mapq 0" \
    --aligner=bwa \
    --bwa_alg=mem \
    bwa_K=100000000 \
    --mapq=0 \
    --req_flg=false \
    --qname=false \
    --index=tiny.fa \
    "${k_pe[@]}" \
    "skipped=MAPQ filter (--mapq=0)"

add_row \
    "BWA-backtrack SE" \
    aln_se.bam \
    "--aligner bwa --bwa_alg aln --index ${idx_bwa} ${a_se} --mapq 20" \
    --aligner=bwa \
    --bwa_alg=aln \
    --mapq=20 \
    --qname=false \
    --index=tiny.fa \
    "${k_se[@]}"

add_row \
    "BWA-MEM2 PE, CRAM, '--qname'" \
    bm2_pe.cram \
    "--aligner bwa-mem2 --index ${idx_bm2} ${a_ref} ${a_pe} --qname" \
    --aligner=bwa-mem2 \
    bwa_K=100000000 \
    --mapq=1 \
    --req_flg=false \
    --qname=true \
    --index=tiny.fa \
    --ref_fa=tiny.fa \
    "${k_pe[@]}"

for idx in "${!rows[@]}"; do
    IFS='|' read -r lbl nam arg lin <<< "${rows[idx]}"
    read -r -a args <<< "${arg}"
    IFS=';' read -r -a arr_lin <<< "${lin}"
    dir_out="${tmp}/out_${idx}"
    fil_out="${dir_out}/${nam}"
    mkdir -p "${dir_out}"

    rc=0
    run_helper --threads 1 "${args[@]}" --fil_out "${fil_out}" \
        > "${dir_out}/run.log" 2>&1 || rc=$?
    if (( rc != 0 )) || [[ ! -s "${fil_out}" ]]; then
        record_fail "${lbl}: run failed (exit ${rc}); see '${dir_out}/run.log'"
        continue
    fi

    # The full paths, in the order the record gives them.
    pth=()
    for key in index fq_1 fq_2 ref_fa; do
        val="$(print_arg_val "--${key}" "${args[@]}")"
        if [[ -n "${val}" ]]; then pth+=( "${key}" "${val}" ); fi
    done
    pth+=( fil_out "${fil_out}" )

    exp_cnt="$(print_cnt_exp "${arr_lin[@]}")"
    exp_pth="$(print_pth_exp "${pth[@]}")"
    check_record "${lbl}" "${fil_out}" "${exp_cnt}" "${exp_pth}" 1

    # A direct call to the library gives the library function alone.
    if [[
        "$(print_co "${fil_out}")" \
            == *$'pipeline-run: calls=align_fastqs.sh::align_fastqs\n'*
    ]]; then
        record_pass "${lbl}: the call path is the library alone"
    else
        record_fail "${lbl}: the call path is wrong for a direct call"
    fi

    # The tools' own '@PG' lines are kept: the aligner first, the record's
    # 'samtools reheader' among them.
    pg="$(
        run_samtools view -H --no-PG -T "${ref_fa}" "${fil_out}" \
            | grep '^@PG' || true
    )"
    if [[
        "$(head -n 1 <<< "${pg}" | cut -f 3)" =~ ^PN:(bowtie2|bwa|bwa-mem2)$
        && "${pg}" == *"CL:samtools reheader "*
    ]]; then
        record_pass "${lbl}: tools' '@PG' lines are kept"
    else
        record_fail "${lbl}: an '@PG' line is missing"
    fi

    # The index matches the rewritten file: fetching each reference by region
    # finds every record.
    if [[ "${fil_out}" == *.bam ]]; then
        n_all="$(run_samtools view -c "${fil_out}")"
        n_reg=0
        while read -r chr; do
            n="$(run_samtools view -c "${fil_out}" "${chr}" 2> /dev/null)" \
                || n=-1
            n_reg=$(( n_reg + n ))
        done < <(
            run_samtools view -H --no-PG "${fil_out}" \
                | sed -n $'s/^@SQ\tSN:\\([^\t]*\\).*/\\1/p'
        )
        if (( n_all > 0 && n_reg == n_all )); then
            record_pass "${lbl}: the index finds every record by region"
        else
            record_fail "${lbl}: the index found ${n_reg} of ${n_all} records"
        fi
    fi

    # The '--qname' copy carries the same record.
    if [[ "${arg}" == *"--qname"* ]]; then
        fil_qnm="${fil_out%.*}.qnam.${fil_out##*.}"
        if [[
            -s "${fil_qnm}"
            && "$(print_co "${fil_qnm}")" == "$(print_co "${fil_out}")"
        ]]; then
            record_pass "${lbl}: the '--qname' copy carries the same record"
        else
            record_fail "${lbl}: the '--qname' copy's record differs"
        fi
    fi

    if [[
        -n "$(find "${dir_out}" -name '*.header.sam' -o -name '*.record.bam')"
    ]]; then
        record_fail "${lbl}: the record's temporary files are left behind"
    else
        record_pass "${lbl}: the record's temporary files are removed"
    fi
done


# A dry run prints the content lines a real run writes.
dir_dry="${tmp}/dry"
mkdir -p "${dir_dry}"
rc=0
out="$(
    run_helper \
        --dry_run \
        --threads 1 \
        --aligner bwa \
        --bwa_alg aln \
        --index "${idx_bwa}" \
        --fq_1 "${in_se}" \
        --mapq 20 \
        --fil_out "${dir_dry}/aln_se.bam" 2> /dev/null
)" || rc=$?

# Decode the printed lines, "    $'@CO\tpipeline: ...' \", into header lines.
# shellcheck disable=SC2016  # The program is for awk, not the shell.
dec="$(
    awk -v q="'" '
        substr($0, 1, 6) == "    $" q {
            lin = substr($0, 7)
            sub(q " \\\\$", "", lin)
            if (index(lin, "@CO\\tpipeline") == 1) {
                sub(/\\t/, "\t", lin)
                print lin
            }
        }
    ' <<< "${out}"
)"
got="$(print_cnt "${dec}" | grep -v -E 'commit=|optical_distance=' || true)"
exp="$(
    print_cnt "$(print_co "${tmp}/out_4/aln_se.bam")" \
        | grep -v -E 'commit=|optical_distance=|skipped=optical ' \
        || true
)"
if [[ "${rc}" -eq 0 && -n "${exp}" && "${got}" == "${exp}" ]]; then
    record_pass "dry run prints the content lines a real run writes"
else
    record_fail "dry run content lines differ from a real run's (exit ${rc})"
    diff <(printf '%s\n' "${exp}") <(printf '%s\n' "${got}") >&2 || true
fi

# It prints the same run lines too, by key: their values, such as the start
# time and the paths, differ between the two runs.
function print_run_keys() {
    print_run "${1}" | sed -E $'s/^@CO\tpipeline-run: ([^=]+)=.*$/\\1/'
}

got="$(print_run_keys "${dec}")"
exp="$(print_run_keys "$(print_co "${tmp}/out_4/aln_se.bam")")"
if [[ -n "${exp}" && "${got}" == "${exp}" ]]; then
    record_pass "dry run prints the run lines a real run writes, by key"
else
    record_fail "dry run run lines differ from a real run's, by key"
    diff <(printf '%s\n' "${exp}") <(printf '%s\n' "${got}") >&2 || true
fi


# Step 5 skips a file that already carries a record, so writing twice leaves
# one record.
fil_two="${tmp}/twice.bam"
cp "${tmp}/out_4/aln_se.bam" "${fil_two}"
rc=0
# shellcheck disable=SC2016  # Expand in the child shell, not this one.
"${TEST_BASH}" -c '
    # shellcheck source=lib/bash/core/source_helpers.sh
    source "${1}/lib/bash/core/source_helpers.sh"

    source_helpers "${1}/lib/bash" workflows/align_fastqs
    arr_cnt=( "${3}" )
    arr_run=( "${4}" )
    _write_align_record align_fastqs false 1 "${2}" arr_cnt arr_run
' _ "${ROOT_REPO}" "${fil_two}" \
    $'@CO\tpipeline: chain=align_fastqs' \
    $'@CO\tpipeline-run: --threads=1' \
    > "${tmp}/twice.log" 2>&1 \
    || rc=$?

if [[
    "${rc}" -eq 0
    && "$(print_co "${fil_two}" \
        | grep -c $'^@CO\tpipeline: chain=align_fastqs$')" == 1
    && "$(cat "${tmp}/twice.log")" == *"already carries a pipeline record"*
]]; then
    record_pass "writing the record twice leaves one record"
else
    record_fail "writing the record twice added a second record (exit ${rc})"
fi


# Two thread counts give the same content lines and different run lines.
for thr in 1 2; do
    mkdir -p "${tmp}/thr_${thr}"
    run_helper \
        --threads "${thr}" \
        --aligner bowtie2 \
        --index "${idx_bt2}" \
        --fq_1 "${in_pe_1}" \
        --fq_2 "${in_pe_2}" \
        --fil_out "${tmp}/thr_${thr}/x.bam" \
        > "${tmp}/thr_${thr}/run.log" 2>&1 || true
done

co_1="$(print_co "${tmp}/thr_1/x.bam")"
co_2="$(print_co "${tmp}/thr_2/x.bam")"
if [[
    -n "${co_1}"
    && "$(print_cnt "${co_1}")" == "$(print_cnt "${co_2}")"
    && "$(print_run "${co_1}")" != "$(print_run "${co_2}")"
    && "$(print_run "${co_2}")" == *$'pipeline-run: --threads=2'*
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
"${TEST_BASH}" "${ROOT_REPO}/bin/execute_align_fastqs.sh" \
    --env_nam "${env_nam}" \
    --threads 1 \
    --aligner bowtie2 \
    --bt2_mode global \
    --mapq 0 \
    --index "${idx_bt2}" \
    --csv_fil_in "${in_se}" \
    --dir_out "${dir_drv}/out" \
    --out_ext bam \
    --sfx_se ".atria.fastq.gz" \
    --sfx_pe "_R1.atria.fastq.gz" \
    --dir_eo "${dir_drv}/logs" \
    --nam_job test_align_record \
    --max_job 1 \
    > "${dir_drv}/run.log" 2>&1 || rc=$?

fil_drv="${dir_drv}/out/tiny_se.bam"
if (( rc != 0 )) || [[ ! -s "${fil_drv}" ]]; then
    record_fail "driver run failed (exit ${rc}); see '${dir_drv}/run.log'"
else
    add_row \
        "driver" \
        "" \
        "" \
        --aligner=bowtie2 \
        --bt2_mode=global \
        bt2_seed=0 \
        --mapq=0 \
        --qname=false \
        --index=tiny \
        "${k_se[@]}" \
        "skipped=MAPQ filter (--mapq=0)"
    IFS='|' read -r _ _ _ lin <<< "${rows[-1]}"
    IFS=';' read -r -a arr_lin <<< "${lin}"

    check_record "driver" "${fil_drv}" \
        "$(print_cnt_exp "${arr_lin[@]}")" \
        "$(print_pth_exp index "${idx_bt2}" fq_1 "${in_se}" \
            fil_out "${fil_drv}")" \
        1

    co="$(print_co "${fil_drv}")"
    if [[
        "${co}" == *$'pipeline-run: calls=execute_align_fastqs.sh > '*
        && "${co}" == *$'submit_align_fastqs.sh > '*
        && "${co}" == *$'align_fastqs.sh::align_fastqs\n'*
        && "${co}" == *$'pipeline-run: env='"${env_nam}"$'\n'*
    ]]; then
        record_pass "driver: the call path and environment are recorded"
    else
        record_fail "driver: the call path or environment line is wrong"
        print_run "${co}" >&2
    fi
fi


finish
