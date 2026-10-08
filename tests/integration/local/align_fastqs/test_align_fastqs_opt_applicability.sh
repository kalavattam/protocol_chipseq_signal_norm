#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_align_fastqs_opt_applicability.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="align-fastqs option applicability"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Define fixture and output paths. Per 'HELP.PARAMETER.APPLICABILITY', each
# aligner refuses the other's option; '--ref_fa' without CRAM and an all-SE
# '--req_flg' are ignored with a warning.
dir_fx="${ROOT_REPO}/tests/fixtures/align_fastqs"
in_se="${dir_fx}/fastq/se/tiny_se.atria.fastq.gz"
in_pe_1="${dir_fx}/fastq/pe/tiny_pe_R1.atria.fastq.gz"
in_pe_2="${dir_fx}/fastq/pe/tiny_pe_R2.atria.fastq.gz"
in_pe="${in_pe_1},${in_pe_2}"
ref_fa="${dir_fx}/reference/tiny.fa"
idx_bt2="${dir_fx}/bowtie2/tiny"
idx_bwa="${dir_fx}/bwa/tiny.fa"
idx_bm2="${dir_fx}/bwa-mem2/tiny.fa"

tmp="${TEST_DIR_TMP}/align_fastqs_opt_applicability"
msg_ref="'--ref_fa' has no effect"
msg_req="'--req_flg' has no effect with single-end input"
msg_bm2="'--bwa_alg' has no effect with '--aligner bwa-mem2'"

# Each aligner with its index.
a_bt2="--aligner bowtie2 --index ${idx_bt2}"
a_bwa="--aligner bwa --index ${idx_bwa}"
a_bm2="--aligner bwa-mem2 --index ${idx_bm2}"


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${tmp}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty \
    "${in_se}" \
    "${in_pe_1}" \
    "${in_pe_2}" \
    "${ref_fa}" \
    "${idx_bt2}.1.bt2" \
    "${idx_bwa}" \
    "${idx_bm2}" || {
    finish
    exit $?
}


# 'execute' rows are dry runs: label|input ('se' or 'mix')|arguments|act.
# Refused rows name the option; warned rows warn and drop it; control rows do
# neither and pass on the option they name.
ref_cram="--ref_fa ${ref_fa} --out_ext cram"
rows_exe=(
    "--bt2_mode with bwa|se|${a_bwa} --bt2_mode local|refuse --bt2_mode"
    "--bt2_mode with bwa-mem2|se|${a_bm2} --bt2_mode local|refuse --bt2_mode"
    "--bwa_alg with bowtie2|se|${a_bt2} --bwa_alg aln|refuse --bwa_alg"
    "default --bwa_alg with bowtie2|se|${a_bt2} --bwa_alg mem|refuse --bwa_alg"
    "--bwa_alg mem with bwa-mem2|se|${a_bm2} --bwa_alg mem|warn --bwa_alg"
    "--bwa_alg aln with bwa-mem2|se|${a_bm2} --bwa_alg aln|refuse --bwa_alg"
    "--bt2_mode with bowtie2|se|${a_bt2} --bt2_mode local|pass --bt2_mode"
    "--ref_fa with BAM output|se|${a_bt2} --ref_fa ${ref_fa}|warn --ref_fa"
    "--ref_fa with CRAM output|se|${a_bt2} ${ref_cram}|pass --ref_fa"
    "--req_flg with SE only|se|${a_bt2} --req_flg|warn --req_flg"
    "--req_flg with a PE entry|mix|${a_bt2} --req_flg|pass --req_flg"
    "--time without --slurm|se|${a_bt2} --time 1:00:00|warn --time"
)

for idx in "${!rows_exe[@]}"; do
    IFS='|' read -r lbl inp arg act <<< "${rows_exe[idx]}"
    csv="${in_se}"
    if [[ "${inp}" == "mix" ]]; then csv="${in_se};${in_pe}"; fi
    read -r -a args <<< "${arg}"
    read -r -a acts <<< "${act}"
    dir_out="${tmp}/exe_${idx}"

    mkdir -p "${dir_out}/logs"

    rc=0
    out="$(
        "${TEST_BASH}" "${ROOT_REPO}/bin/execute_align_fastqs.sh" \
            --dry_run \
            --threads 1 \
            --csv_fil_in "${csv}" \
            --dir_out "${dir_out}" \
            --max_job 1 \
            "${args[@]}" 2>&1
    )" || rc=$?

    cmd="$(grep -m1 'submit_align_fastqs.sh --env' <<< "${out}" || true)"
    lbl="execute ${lbl}"
    opt="${acts[1]:-}"

    case "${acts[0]}" in
        refuse)
            if [[
                "${rc}" -ne 0
                && "${out}" == *"error("*"'${opt}' is for '--aligner "*
            ]]; then
                record_pass "${lbl} is refused by name"
            else
                record_fail "${lbl} was not refused by name (exit ${rc})"
            fi
            ;;

        warn)
            if [[
                "${rc}" -eq 0
                && "${out}" == *"warning("*"'${opt}' has no effect"*
                && -n "${cmd}"
                && "${cmd}" != *" ${opt}"*
            ]]; then
                record_pass "${lbl} warns and is not passed on"
            else
                record_fail \
                    "${lbl} did not warn, or was passed on (exit ${rc})"
            fi
            ;;

        pass)
            if [[
                "${rc}" -eq 0
                && -n "${cmd}"
                && ( -z "${opt}" || "${cmd}" == *" ${opt}"* )
                && "${out}" != *"has no effect"*
            ]]; then
                record_pass "${lbl} is accepted with no message"
            else
                record_fail "${lbl} was dropped or flagged (exit ${rc})"
            fi
            ;;
    esac
done

# Run in parallel, a mixed list gets one command per entry, and only the
# paired-end entry gets '--req_flg', so the single-end entry's log stays free
# of a warning.
dir_out="${tmp}/exe_mix_parallel"
read -r -a args <<< "${a_bt2}"
mkdir -p "${dir_out}/logs"

rc=0
out="$(
    "${TEST_BASH}" "${ROOT_REPO}/bin/execute_align_fastqs.sh" \
        --dry_run \
        --threads 4 \
        --csv_fil_in "${in_se};${in_pe}" \
        --dir_out "${dir_out}" \
        --max_job 2 \
        "${args[@]}" \
        --req_flg 2>&1
)" || rc=$?

cmd_se="$(
    grep 'submit_align_fastqs.sh --env' <<< "${out}" \
        | grep -F -- "--csv_fil_in ${in_se} " \
        || true
)"
cmd_pe="$(
    grep 'submit_align_fastqs.sh --env' <<< "${out}" \
        | grep -F -- "${in_pe_2}" \
        || true
)"

if [[
    "${rc}" -eq 0
    && -n "${cmd_se}"
    && -n "${cmd_pe}"
    && "${cmd_se}" != *" --req_flg"*
    && "${cmd_pe}" == *" --req_flg"*
]]; then
    record_pass "execute passes '--req_flg' only to a parallel PE entry"
else
    record_fail "execute passed '--req_flg' to a parallel SE entry (exit ${rc})"
fi


# Run 'submit' on a list, writing to its own output directory.
function run_submit() {
    local dir_out="${1}"
    shift 1

    mkdir -p "${dir_out}/logs"
    "${TEST_BASH}" "${ROOT_REPO}/bin/submit_align_fastqs.sh" \
        --env_nam "${env_nam}" \
        --dir_scr "${ROOT_REPO}/bin" \
        --threads 1 \
        --dir_out "${dir_out}" \
        --dir_eo "${dir_out}/logs" \
        --sfx_se ".atria.fastq.gz" \
        --sfx_pe "_R1.atria.fastq.gz" \
        "$@" 2>&1
}


# A refusal in 'submit' stops before any alignment: option|arguments|message.
rows_ref=(
    "--bt2_mode|${a_bwa} --bt2_mode local|'--bt2_mode' is for '--aligner "
    "--bwa_alg|${a_bt2} --bwa_alg mem|'--bwa_alg' is for '--aligner "
    "--bwa_alg aln|${a_bm2} --bwa_alg aln|'--bwa_alg' is for '--aligner bwa';"
)

for idx in "${!rows_ref[@]}"; do
    IFS='|' read -r opt arg msg <<< "${rows_ref[idx]}"
    read -r -a args <<< "${arg}"
    dir_out="${tmp}/sub_refuse_${idx}"
    rc=0
    out="$(
        run_submit "${dir_out}" --csv_fil_in "${in_se}" "${args[@]}"
    )" || rc=$?

    if [[
        "${rc}" -ne 0
        && "${out}" == *"error("*"${msg}"*
        && ! -e "${dir_out}/tiny_se.bam"
    ]]; then
        record_pass "submit '${opt}' is refused before aligning"
    else
        record_fail "submit '${opt}' was not refused (exit ${rc})"
    fi
done

# With BAM output, only single-end entries, or '--bwa_alg mem' for 'bwa-mem2',
# 'submit' warns and aligns.
rows_wrn=(
    "--ref_fa|${msg_ref}|${a_bt2} --ref_fa ${ref_fa}"
    "--req_flg|${msg_req}|${a_bt2} --req_flg"
    "--bwa_alg mem|${msg_bm2}|${a_bm2} --bwa_alg mem"
)

for idx in "${!rows_wrn[@]}"; do
    IFS='|' read -r opt msg arg <<< "${rows_wrn[idx]}"
    read -r -a args <<< "${arg}"
    dir_out="${tmp}/sub_warn_${idx}"
    rc=0
    out="$(
        run_submit "${dir_out}" --csv_fil_in "${in_se}" "${args[@]}"
    )" || rc=$?

    if [[
        "${rc}" -eq 0
        && "${out}" == *"warning("*"${msg}"*
        && -s "${dir_out}/tiny_se.bam"
    ]] && ! grep -rq "has no effect" "${dir_out}/logs"; then
        record_pass "submit '${opt}' warns once and aligns"
    else
        record_fail \
            "submit '${opt}' did not warn, or did not align (exit ${rc})"
    fi
done

# A mixed list passes '--req_flg' only to its paired-end entry, so no layer
# warns; the trace shows how often it was passed.
dir_out="${tmp}/sub_mixed"
mkdir -p "${dir_out}/logs"

rc=0
out="$(
    "${TEST_BASH}" -x "${ROOT_REPO}/bin/submit_align_fastqs.sh" \
        --env_nam "${env_nam}" \
        --dir_scr "${ROOT_REPO}/bin" \
        --threads 1 \
        --aligner bowtie2 \
        --index "${idx_bt2}" \
        --csv_fil_in "${in_se};${in_pe}" \
        --dir_out "${dir_out}" \
        --dir_eo "${dir_out}/logs" \
        --sfx_se ".atria.fastq.gz" \
        --sfx_pe "_R1.atria.fastq.gz" \
        --req_flg 2>&1
)" || rc=$?
n_req="$(grep -c "cmd_aln+=(--req_flg)" <<< "${out}" || true)"

if [[
    "${rc}" -eq 0
    && "${n_req}" -eq 1
    && "${out}" != *"has no effect"*
    && -s "${dir_out}/tiny_se.bam"
    && -s "${dir_out}/tiny_pe.bam"
]] && ! grep -rq "has no effect" "${dir_out}/logs"; then
    record_pass "submit passes '--req_flg' only to a mixed list's PE entry"
else
    record_fail \
        "submit '--req_flg' with a mixed list warned, or reached the SE" \
        "entry (exit ${rc}; passed ${n_req} times)"
fi

# 'submit' passes '--mapq' on every time, so its default of 1 and an explicit 0
# both reach the helper; the mixed run above gave no '--mapq'. The trace shows
# each helper call with its arguments expanded.
patn_call='^\++ align_fastqs '
n_mq1="$(
    grep -E "${patn_call}" <<< "${out}" | grep -cE -- '--mapq 1( |$)' || true
)"
dir_out="${tmp}/sub_mapq0"
mkdir -p "${dir_out}/logs"
out="$(
    "${TEST_BASH}" -x "${ROOT_REPO}/bin/submit_align_fastqs.sh" \
        --env_nam "${env_nam}" \
        --dir_scr "${ROOT_REPO}/bin" \
        --threads 1 \
        --aligner bowtie2 \
        --index "${idx_bt2}" \
        --csv_fil_in "${in_se}" \
        --dir_out "${dir_out}" \
        --dir_eo "${dir_out}/logs" \
        --sfx_se ".atria.fastq.gz" \
        --sfx_pe "_R1.atria.fastq.gz" \
        --mapq 0 2>&1
)" || true
n_mq0="$(
    grep -E "${patn_call}" <<< "${out}" | grep -cE -- '--mapq 0( |$)' || true
)"

if [[ "${n_mq1}" -eq 2 && "${n_mq0}" -eq 1 ]]; then
    record_pass "submit passes its '--mapq' default and an explicit 0 on"
else
    record_fail \
        "submit did not pass '--mapq' on (default ${n_mq1} of 2, explicit 0" \
        "${n_mq0} of 1)"
fi


# The helper checks on its own.
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
    ' _ "${ROOT_REPO}" --threads 1 --fq_1 "${in_se}" "$@" 2>&1
}


# label|arguments|act|message
ref_bt2="'--bt2_mode' is for '--aligner bowtie2'"
ref_bwa="'--bwa_alg' is for '--aligner bwa'"
rows_hlp=(
    "--bt2_mode with bwa|${a_bwa} --bt2_mode local|refuse|${ref_bt2}"
    "default --bwa_alg with bowtie2|${a_bt2} --bwa_alg mem|refuse|${ref_bwa}"
    "--bwa_alg aln with bwa-mem2|${a_bm2} --bwa_alg aln|refuse|${ref_bwa}"
    "--bwa_alg mem with bwa-mem2|${a_bm2} --bwa_alg mem|warn|${msg_bm2}"
    "--req_flg with SE input|${a_bt2} --req_flg|warn|${msg_req}"
    "--ref_fa with BAM output|${a_bt2} --ref_fa ${ref_fa}|warn|${msg_ref}"
    "--bt2_mode in capitals|${a_bt2} --bt2_mode Local|pass|"
)

for idx in "${!rows_hlp[@]}"; do
    IFS='|' read -r lbl arg act msg <<< "${rows_hlp[idx]}"
    read -r -a args <<< "${arg}"
    fil_out="${tmp}/helper_${idx}.bam"
    lbl="helper ${lbl}"

    rc=0
    out="$(run_helper "${args[@]}" --fil_out "${fil_out}")" || rc=$?

    case "${act}" in
        refuse)
            if [[
                "${rc}" -ne 0
                && "${out}" == *"error("*"::align_fastqs): ${msg}"*
                && ! -e "${fil_out}"
            ]]; then
                record_pass "${lbl} is refused by name"
            else
                record_fail "${lbl} was not refused by name (exit ${rc})"
            fi
            ;;

        warn)
            if [[
                "${rc}" -eq 0
                && "${out}" == *"warning("*"::align_fastqs): ${msg}"*
                && -s "${fil_out}"
            ]]; then
                record_pass "${lbl} warns and aligns"
            else
                record_fail \
                    "${lbl} did not warn, or did not align (exit ${rc})"
            fi
            ;;

        pass)
            if [[ "${rc}" -eq 0 && -s "${fil_out}" ]]; then
                record_pass "${lbl} is accepted"
            else
                record_fail "${lbl} was not accepted (exit ${rc})"
            fi
            ;;
    esac
done


# Called directly with no '--mapq', the helper filters at its default of 1; a
# dry run prints the 'samtools view' call it builds.
out="$(
    run_helper \
        --dry_run \
        --aligner bowtie2 \
        --index "${idx_bt2}" \
        --fil_out "${tmp}/helper_mapq.bam"
)" || true

if [[ "${out}" == *"| samtools view \\"*"-F 4 \\"$'\n'"        -q 1 \\"* ]]; then
    record_pass "helper defaults '--mapq' to 1"
else
    record_fail "helper did not default '--mapq' to 1"
fi

# The 'bt2' alias is gone, and every layer words an invalid aligner alike.
msg_aln="'--aligner' must be 'bowtie2', 'bwa', or 'bwa-mem2': 'bt2'."
arg_bt2=( --aligner bt2 --index "${idx_bt2}" )

declare -A out_lyr
out_lyr[execute]="$(
    "${TEST_BASH}" "${ROOT_REPO}/bin/execute_align_fastqs.sh" \
        --dry_run \
        --threads 1 \
        --csv_fil_in "${in_se}" \
        --dir_out "${tmp}/exe_0" \
        --max_job 1 \
        "${arg_bt2[@]}" 2>&1
)" || true
out_lyr[submit]="$(
    run_submit "${tmp}/sub_bt2" --csv_fil_in "${in_se}" "${arg_bt2[@]}"
)" || true
out_lyr[helper]="$(
    run_helper "${arg_bt2[@]}" --fil_out "${tmp}/helper_bt2.bam"
)" || true

for lyr in execute submit helper; do
    if [[ "${out_lyr[${lyr}]}" == *"error("*"${msg_aln}"* ]]; then
        record_pass "${lyr} refuses '--aligner bt2' with the shared message"
    else
        record_fail "${lyr} did not refuse '--aligner bt2' as expected"
    fi
done

# Validation lets through a Bowtie 2 index stem with no index files; Step 1
# then fails, which must stop the helper with one error in the helper's name.
fil_out="${tmp}/helper_step1.bam"
rc=0
out="$(
    run_helper \
        --aligner bowtie2 \
        --index "${tmp}/no_index" \
        --fil_out "${fil_out}"
)" || rc=$?
n_err="$(grep -c '^error(' <<< "${out}" || true)"

if [[
    "${rc}" -ne 0
    && "${n_err}" -eq 1
    && "${out}" == *"::align_fastqs): Step #1: failed to align single-end"*
    && ! -e "${fil_out}"
]]; then
    record_pass "helper stops at a Step 1 failure with one error by name"
else
    record_fail \
        "helper did not stop at a Step 1 failure with one error by name" \
        "(exit ${rc}; ${n_err} errors)"
fi


finish
