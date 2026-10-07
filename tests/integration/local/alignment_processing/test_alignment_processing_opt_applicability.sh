#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_alignment_processing_opt_applicability.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="alignment-processing option applicability"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Define fixture and output paths. Per 'HELP.PARAMETER.APPLICABILITY',
# '--ref_fa' without CRAM input is ignored with a warning and not checked, so
# a missing path is not refused; with a mixed list, only CRAM entries get it.
dir_fx="${ROOT_REPO}/tests/fixtures/calculate_scaling_factor"
in_bam="${dir_fx}/bam/pe/IP_WT_G1_Hho1_6336.sc.bam"
in_cram="${dir_fx}/cram/pe/IP_WT_G1_Hho1_6337.sc.cram"
ref_fa="${dir_fx}/reference/tiny.fa"

tmp="${TEST_DIR_TMP}/alignment_processing_opt_applicability"
q_bam="${tmp}/qsort_0/IP_WT_G1_Hho1_6336.sc.qnam.bam"


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${tmp}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty "${in_bam}" "${in_cram}" "${ref_fa}" || {
    finish
    exit $?
}


# Run one 'submit' driver ('qsort_bam' or 'convert_bam_bed'), writing to its
# own output directory.
function run_submit() {
    local drv="${1}" dir_out="${2}"
    shift 2

    mkdir -p "${dir_out}/logs"
    "${TEST_BASH}" "${ROOT_REPO}/bin/submit_${drv}.sh" \
        --env_nam "${env_nam}" \
        --dir_scr "${ROOT_REPO}/bin" \
        --dir_out "${dir_out}" \
        --dir_eo "${dir_out}/logs" \
        "$@" 2>&1
}


# label|driver|inputs|arguments|expected outputs. Rows given '--ref_fa' with
# BAM only must warn once, as the driver drops the reference before the helper;
# mixed rows must not warn.
mix="${in_bam},${in_cram}"
ref="--ref_fa ${ref_fa}"
no_ref="--ref_fa ${tmp}/no.fa"

rows=(
    "--ref_fa with BAM only|qsort_bam|${in_bam}|${ref}|1"
    "missing --ref_fa, BAM only|qsort_bam|${in_bam}|${no_ref}|1"
    "--ref_fa with a CRAM entry|qsort_bam|${mix}|${ref}|2"
    "--ref_fa with BAM only|convert_bam_bed|${in_bam}|${ref}|1"
    "missing --ref_fa, BAM only|convert_bam_bed|${in_bam}|${no_ref}|1"
    "--ref_fa with a CRAM entry|convert_bam_bed|${mix}|${ref}|2"
    "--ref_fa with BAM only, AWK|convert_bam_bed|${q_bam}|${ref} --use_awk|1"
)

for idx in "${!rows[@]}"; do
    IFS='|' read -r lbl drv csv arg n_exp <<< "${rows[idx]}"
    read -r -a args <<< "${arg}"
    dir_out="${tmp}/${drv%%_*}_${idx}"
    lbl="submit_${drv} ${lbl}"

    rc=0
    out="$(
        run_submit "${drv}" "${dir_out}" \
            --threads 1 --csv_fil_in "${csv}" "${args[@]}"
    )" || rc=$?

    n_out="$(find "${dir_out}" -maxdepth 1 -type f | wc -l | tr -d ' ')"
    n_msg="$(grep -c 'has no effect' <<< "${out}" || true)"

    if [[ "${csv}" == *.cram ]]; then
        n_msg_exp=0
    else
        n_msg_exp=1
    fi

    if [[
        "${rc}" -eq 0
        && "${n_out}" -eq "${n_exp}"
        && "${n_msg}" -eq "${n_msg_exp}"
    ]]; then
        if (( n_msg_exp )); then
            record_pass "${lbl} warns once and still runs"
        else
            record_pass "${lbl} runs with no message"
        fi
    else
        record_fail \
            "${lbl}: exit ${rc}, ${n_out} of ${n_exp} outputs," \
            "${n_msg} of ${n_msg_exp} messages"
    fi
done


# 'convert_bam_bed' passes '--threads' to the Python converter. A converter
# that writes its arguments shows what it got.
fil_py="${tmp}/print_args.py"
cat > "${fil_py}" << 'PY'
import sys

with open(sys.argv[sys.argv.index("--fil_out") + 1], "w") as fil:
    fil.write(" ".join(sys.argv[1:]) + "\n")
PY

dir_out="${tmp}/convert_threads"
rc=0
out="$(
    run_submit convert_bam_bed "${dir_out}" \
        --threads 3 --csv_fil_in "${in_bam}" --pth_scr_py "${fil_py}"
)" || rc=$?
args_py="$(cat "${dir_out}/IP_WT_G1_Hho1_6336.sc.bed.gz" 2> /dev/null || true)"

if [[ "${rc}" -eq 0 && " ${args_py} " == *" --threads 3 "* ]]; then
    record_pass "submit_convert_bam_bed passes '--threads' to the converter"
else
    record_fail \
        "submit_convert_bam_bed did not pass '--threads' (exit ${rc};" \
        "converter got '${args_py}')"
fi


# Call a helper directly from a child shell.
function run_helper() {
    # shellcheck disable=SC2016  # Expand in the child shell, not this one.
    "${TEST_BASH}" -c '
        # shellcheck disable=SC1090
        source "${1}/lib/bash/core/source_helpers.sh"

        source_helpers "${1}/lib/bash" \
            core/check_args \
            core/check_inputs \
            core/format_outputs \
            workflows/process_alignments
        shift 1
        "$@"
    ' _ "${ROOT_REPO}" "$@" 2>&1
}

# Each helper warns about a reference given with BAM input and still runs,
# without checking it, so a missing path is not refused. The Python helper
# does not pass it to the converter: helper and leading arguments|output.
dir_hlp="${tmp}/helper"
mkdir -p "${dir_hlp}"

rows_hlp=(
    "qsort_file_alignments 1|${dir_hlp}/q.bam"
    "convert_alignments_bed_awk 1|${dir_hlp}/a.bed.gz"
    "convert_alignments_bed_python 1 ${fil_py}|${dir_hlp}/p.bed.gz"
)

for row in "${rows_hlp[@]}"; do
    IFS='|' read -r arg fil_out <<< "${row}"
    read -r -a args <<< "${arg}"
    lbl="${args[0]}"
    fil_in="${in_bam}"
    if [[ "${lbl}" == *awk ]]; then fil_in="${q_bam}"; fi

    rc=0
    out="$(
        run_helper "${args[@]}" "${fil_in}" "${fil_out}" "${tmp}/no.fa" \
            "${dir_hlp}/out.txt" "${dir_hlp}/err.txt"
    )" || rc=$?

    if [[
        "${rc}" -eq 0
        && -s "${fil_out}"
        && "${out}" == *"'ref_fa' has no effect without CRAM input"*
    ]] && ! grep -qa -- '--ref_fa' "${fil_out}"; then
        record_pass "${lbl} warns about 'ref_fa' with BAM and still runs"
    else
        record_fail \
            "${lbl} did not warn, failed, or passed 'ref_fa' on (exit ${rc})"
    fi
done

finish
