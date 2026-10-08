#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_filter_alignments_opt_applicability.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="filter-alignments option applicability"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Define fixture and output paths. Per 'HELP.PARAMETER.APPLICABILITY', '--tg'
# and '--mtr' are refused with '--retain sc', and '--ref_fa' is ignored with a
# warning when no input or output is CRAM.
dir_fx="${ROOT_REPO}/tests/fixtures/filter_alignments"
in_sam="${dir_fx}/sam/filter_sc_sp.sam"
ref_fa="${dir_fx}/reference/filter_sc_sp.fa"

tmp="${TEST_DIR_TMP}/filter_alignments_opt_applicability"
dir_in="${tmp}/input"
dir_log="${TEST_DIR_LOG}/filter_alignments"
in_bam="${dir_in}/bam_in.bam"
in_cram="${dir_in}/cram_in.cram"
msg_ref="'--ref_fa' has no effect without CRAM input or output"


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${dir_in}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty "${in_sam}" "${ref_fa}" "${ref_fa}.fai" || {
    finish
    exit $?
}

if ! \
    build_filter_alignments_fixture_bam \
        "${in_sam}" \
        "${in_bam}" \
        "${dir_log}/opt_applicability_prepare_bam.log" \
        "filter-alignments applicability BAM fixture" \
|| ! \
    build_filter_alignments_fixture_cram \
        "${in_sam}" \
        "${ref_fa}" \
        "${in_cram}" \
        "${dir_log}/opt_applicability_prepare_cram.log" \
        "filter-alignments applicability CRAM fixture"
then
    finish
    exit $?
fi


# 'execute' rows run as dry runs: label|input list|arguments|act. A refused row
# names the option; a warned row warns and does not pass '--ref_fa' on; a
# control row neither refuses nor warns, and passes on what it was given.
rows_exe=(
    "--tg with sc|${in_bam}|--retain sc --tg|refuse --tg"
    "--mtr with sc|${in_bam}|--retain sc --mtr|refuse --mtr"
    "--tg --mtr with sp|${in_bam}|--retain sp --tg --mtr|pass --tg --mtr"
    "--ref_fa with BAM only|${in_bam}|--ref_fa ${ref_fa}|warn"
    "--ref_fa to CRAM output|${in_bam}|--out_ext cram --ref_fa ${ref_fa}|pass"
    "--ref_fa with a CRAM entry|${in_bam},${in_cram}|--ref_fa ${ref_fa}|pass"
)

for idx in "${!rows_exe[@]}"; do
    IFS='|' read -r lbl csv arg act <<< "${rows_exe[idx]}"
    read -r -a args <<< "${arg}"
    read -r -a acts <<< "${act}"
    dir_out="${tmp}/exe_${idx}"
    mkdir -p "${dir_out}/logs"

    rc=0
    out="$(
        "${TEST_BASH}" "${ROOT_REPO}/bin/execute_filter_alignments.sh" \
            --dry_run \
            --threads 1 \
            --csv_fil_in "${csv}" \
            --dir_out "${dir_out}" \
            --max_job 1 \
            "${args[@]}" 2>&1
    )" || rc=$?

    cmd="$(grep -m1 'submit_filter_alignments.sh --env' <<< "${out}" || true)"
    lbl="execute ${lbl}"

    case "${acts[0]}" in
        refuse)
            if [[
                "${rc}" -ne 0
                && "${out}" == *"error("*"'${acts[1]}' is for '--retain sp'"*
            ]]; then
                record_pass "${lbl} is refused by name"
            else
                record_fail "${lbl} was not refused by name (exit ${rc})"
            fi
            ;;

        warn)
            if [[
                "${rc}" -eq 0
                && "${out}" == *"warning("*"${msg_ref}"*
                && -n "${cmd}"
                && "${cmd}" != *" --ref_fa "*
            ]]; then
                record_pass "${lbl} warns and is not passed on"
            else
                record_fail \
                    "${lbl} did not warn, or was passed on (exit ${rc})"
            fi
            ;;

        pass)
            # Every option the row gives is passed on.
            ok=true
            for opt in "${args[@]}"; do
                if [[ "${opt}" == --* && "${cmd}" != *" ${opt} "* ]]; then
                    ok=false
                fi
            done

            if [[
                "${rc}" -eq 0
                && "${ok}" == "true"
                && "${out}" != *"has no effect"*
            ]]; then
                record_pass "${lbl} is passed on with no message"
            else
                record_fail "${lbl} was dropped or flagged (exit ${rc})"
            fi
            ;;
    esac
done

# Run in parallel, a mixed list gets one command per entry; with BAM output,
# only the CRAM entry gets the reference, so the BAM entry's log stays free of
# a warning. With CRAM output, both entries get it.
for out_ext in bam cram; do
    dir_out="${tmp}/exe_mix_parallel_${out_ext}"
    mkdir -p "${dir_out}/logs"

    rc=0
    out="$(
        "${TEST_BASH}" "${ROOT_REPO}/bin/execute_filter_alignments.sh" \
            --dry_run \
            --threads 4 \
            --csv_fil_in "${in_bam},${in_cram}" \
            --dir_out "${dir_out}" \
            --out_ext "${out_ext}" \
            --ref_fa "${ref_fa}" \
            --max_job 2 2>&1
    )" || rc=$?

    cmds="$(grep 'submit_filter_alignments.sh --env' <<< "${out}" || true)"
    cmd_bam="$(grep -F -- "--csv_fil_in ${in_bam} " <<< "${cmds}" || true)"
    cmd_cram="$(grep -F -- "--csv_fil_in ${in_cram} " <<< "${cmds}" || true)"
    want_bam="no"
    if [[ "${out_ext}" == "cram" ]]; then want_bam="yes"; fi
    has_bam="no"
    if [[ "${cmd_bam}" == *" --ref_fa "* ]]; then has_bam="yes"; fi

    if [[
        "${rc}" -eq 0
        && -n "${cmd_bam}"
        && "${cmd_cram}" == *" --ref_fa "*
        && "${has_bam}" == "${want_bam}"
    ]]; then
        record_pass \
            "execute passes '--ref_fa' per parallel entry with" \
            "'--out_ext ${out_ext}'"
    else
        record_fail \
            "execute passed '--ref_fa' wrongly per parallel entry with" \
            "'--out_ext ${out_ext}' (exit ${rc})"
    fi
done

# A time limit is for Slurm jobs only, so a local run warns about it.
dir_out="${tmp}/exe_time"
mkdir -p "${dir_out}/logs"

rc=0
out="$(
    "${TEST_BASH}" "${ROOT_REPO}/bin/execute_filter_alignments.sh" \
        --dry_run \
        --threads 1 \
        --csv_fil_in "${in_bam}" \
        --dir_out "${dir_out}" \
        --max_job 1 \
        --time 1:00:00 2>&1
)" || rc=$?

if [[
    "${rc}" -eq 0
    && "${out}" == *"warning("*"'--time' has no effect without '--slurm'"*
]]; then
    record_pass "execute warns about '--time' without '--slurm'"
else
    record_fail "execute did not warn about '--time' without '--slurm'"
fi


# Run 'submit' on a list, writing to its own output directory.
function run_submit() {
    local dir_out="${1}"
    shift 1

    mkdir -p "${dir_out}/logs"
    "${TEST_BASH}" "${ROOT_REPO}/bin/submit_filter_alignments.sh" \
        --env_nam "${env_nam}" \
        --dir_scr "${ROOT_REPO}/bin" \
        --threads 1 \
        --dir_out "${dir_out}" \
        --dir_eo "${dir_out}/logs" \
        "$@" 2>&1
}


# Print the '@PG' lines of a filtered output.
function print_pg() {
    run_samtools view -H "${1}" | grep '^@PG' || true
}


# A refusal in 'submit' stops before any file is filtered.
for opt in --tg --mtr; do
    dir_out="${tmp}/sub_refuse_${opt#--}"
    rc=0
    out="$(
        run_submit "${dir_out}" --csv_fil_in "${in_bam}" --retain sc "${opt}"
    )" || rc=$?

    if [[
        "${rc}" -ne 0
        && "${out}" == *"error("*"'${opt}' is for '--retain sp'"*
        && ! -e "${dir_out}/bam_in.sc.bam"
    ]]; then
        record_pass "submit '${opt}' with sc is refused before filtering"
    else
        record_fail "submit '${opt}' with sc was not refused (exit ${rc})"
    fi
done

# With BAM only, 'submit' warns, filters, and records no reference; an ignored
# file is not checked, so a missing one does not stop the run.
dir_out="${tmp}/sub_bam_ref"
rc=0
out="$(
    run_submit "${dir_out}" \
        --csv_fil_in "${in_bam}" \
        --ref_fa "${tmp}/missing.fa"
)" || rc=$?
pg="$(print_pg "${dir_out}/bam_in.sc.bam")"

if [[
    "${rc}" -eq 0
    && "${out}" == *"warning("*"${msg_ref}"*
    && "${pg}" == *"retain=sc"*
    && "${pg}" != *"ref_fa="*
]]; then
    record_pass "submit '--ref_fa' with BAM only warns and is not recorded"
else
    record_fail \
        "submit '--ref_fa' with BAM only did not warn, or was recorded" \
        "(exit ${rc})"
fi

# A mixed list keeps the reference for the CRAM entry only, with no warning
# from the driver or from the helper's per-entry logs.
dir_out="${tmp}/sub_mixed_ref"
rc=0
out="$(
    run_submit "${dir_out}" \
        --csv_fil_in "${in_bam},${in_cram}" \
        --ref_fa "${ref_fa}"
)" || rc=$?
pg_bam="$(print_pg "${dir_out}/bam_in.sc.bam")"
pg_cram="$(print_pg "${dir_out}/cram_in.sc.bam")"

if [[
    "${rc}" -eq 0
    && "${out}" != *"has no effect"*
    && "${pg_bam}" == *"retain=sc"*
    && "${pg_bam}" != *"ref_fa="*
    && "${pg_cram}" == *"ref_fa="*
]] && ! grep -rq "has no effect" "${dir_out}/logs"; then
    record_pass "submit passes '--ref_fa' only to a mixed list's CRAM entry"
else
    record_fail \
        "submit '--ref_fa' with a mixed list warned, or reached the BAM" \
        "entry (exit ${rc})"
fi


# The helpers check on their own: validation warns when called with no CRAM and
# stays quiet with CRAM input, and '@PG' records the reference only with CRAM.
function run_helper() {
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
    ' _ "${ROOT_REPO}" "$@" 2>&1
}

rc=0
out="$(
    run_helper \
        _validate_args_filter_alignment \
            filter_alignment_sc 1 "${in_bam}" "${tmp}/helper.bam" "${ref_fa}"
)" || rc=$?
if [[
    "${rc}" -eq 0
    && "${out}" == *"warning("*"::filter_alignment_sc): ${msg_ref}"*
]]; then
    record_pass "helper warns about '--ref_fa' with no CRAM"
else
    record_fail \
        "helper did not warn about '--ref_fa' with no CRAM (exit ${rc})"
fi

rc=0
out="$(
    run_helper \
        _validate_args_filter_alignment \
            filter_alignment_sc 1 "${in_cram}" "${tmp}/helper.bam" "${ref_fa}"
)" || rc=$?
if [[ "${rc}" -eq 0 && "${out}" != *"has no effect"* ]]; then
    record_pass "helper accepts '--ref_fa' with CRAM input silently"
else
    record_fail "helper flagged '--ref_fa' with CRAM input (exit ${rc})"
fi

pg_bam="$(
    run_helper \
        _build_filter_pg_cl \
            filter_alignment_sc sc \
            "${in_bam}" "${tmp}/helper.bam" \
            false false false false "${ref_fa}"
)"
pg_cram="$(
    run_helper \
        _build_filter_pg_cl \
            filter_alignment_sc sc \
            "${in_bam}" "${tmp}/helper.cram" \
            false false false false "${ref_fa}"
)"
if [[
    "${pg_bam}" == *"out_ext=bam"*
    && "${pg_bam}" != *"ref_fa="*
    && "${pg_cram}" == *"ref_fa=${ref_fa}"*
]]; then
    record_pass "helper '@PG' records '--ref_fa' only with CRAM"
else
    record_fail "helper '@PG' recorded '--ref_fa' wrongly"
fi


finish
