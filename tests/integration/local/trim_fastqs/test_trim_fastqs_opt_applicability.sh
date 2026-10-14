#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_trim_fastqs_opt_applicability.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="trim-fastqs option applicability"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Define fixture and output paths. Per 'HELP.PARAMETER.APPLICABILITY', '--time'
# without '--slurm' and a suffix for a layout no entry has are ignored with a
# warning.
dir_fx="${ROOT_REPO}/tests/fixtures/trim_fastqs"
in_se="${dir_fx}/fastq/se/tiny_se.fastq.gz"
in_pe_1="${dir_fx}/fastq/pe/tiny_pe_R1.fastq.gz"
in_pe_2="${dir_fx}/fastq/pe/tiny_pe_R2.fastq.gz"
in_pe="${in_pe_1},${in_pe_2}"

msg_sse="'--sfx_se' has no effect with paired-end input"
msg_spe="'--sfx_pe' has no effect with single-end input"
sfx_se=".fastq.gz"
sfx_pe="_R1.fastq.gz"

tmp="${TEST_DIR_TMP}/trim_fastqs_opt_applicability"


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${tmp}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty "${in_se}" "${in_pe_1}" "${in_pe_2}" || {
    finish
    exit $?
}


# Dry-run 'execute' locally, writing to its own output directory.
function run_execute() {
    local dir_out="${1}"
    shift 1

    mkdir -p "${dir_out}/logs"
    "${TEST_BASH}" "${ROOT_REPO}/bin/execute_trim_fastqs.sh" \
        --dry_run \
        --env_nam "${env_nam}" \
        --threads 1 \
        --csv_fil_in "${in_se}" \
        --dir_out "${dir_out}" \
        --sfx_se "${sfx_se}" \
        --dir_eo "${dir_out}/logs" \
        --max_job 1 \
        "$@" 2>&1
}


# A time limit is for Slurm jobs only, so a local run warns about it.
rc=0
out="$(run_execute "${tmp}/exe_time" --time 1:00:00)" || rc=$?

if [[
    "${rc}" -eq 0
    && "${out}" == *"warning("*"'--time' has no effect without '--slurm'"*
]]; then
    record_pass "execute warns about '--time' without '--slurm'"
else
    record_fail "execute did not warn about '--time' without '--slurm'"
fi

# Without '--time', the default limit draws no warning.
rc=0
out="$(run_execute "${tmp}/exe_default")" || rc=$?

if [[ "${rc}" -eq 0 && "${out}" != *"has no effect"* ]]; then
    record_pass "execute does not warn about the default time limit"
else
    record_fail "execute warned about the default time limit (exit ${rc})"
fi


# 'execute' rows are dry runs: label|input ('se', 'pe', or 'mix')|arguments|
# act. Warned rows warn and drop the option; control rows do neither and pass
# it on; dropped rows do not warn but pass on no default for the option.
rows_exe=(
    "--sfx_pe with SE only|se|--sfx_pe ${sfx_pe}|warn --sfx_pe"
    "--sfx_se with PE only|pe|--sfx_se ${sfx_se}|warn --sfx_se"
    "--sfx_pe with a PE entry|mix|--sfx_pe ${sfx_pe}|pass --sfx_pe"
    "--sfx_se with an SE entry|mix|--sfx_se ${sfx_se}|pass --sfx_se"
    "default --sfx_pe with SE only|se||drop --sfx_pe"
    "default --sfx_se with PE only|pe||drop --sfx_se"
)

for idx in "${!rows_exe[@]}"; do
    IFS='|' read -r lbl inp arg act <<< "${rows_exe[idx]}"
    csv="${in_se}"
    if [[ "${inp}" == "pe" ]]; then csv="${in_pe}"; fi
    if [[ "${inp}" == "mix" ]]; then csv="${in_se};${in_pe}"; fi
    read -r -a args <<< "${arg}"
    read -r -a acts <<< "${act}"
    dir_out="${tmp}/exe_sfx_${idx}"

    mkdir -p "${dir_out}/logs"

    rc=0
    out="$(
        "${TEST_BASH}" "${ROOT_REPO}/bin/execute_trim_fastqs.sh" \
            --dry_run \
            --env_nam "${env_nam}" \
            --threads 1 \
            --csv_fil_in "${csv}" \
            --dir_out "${dir_out}" \
            --dir_eo "${dir_out}/logs" \
            --max_job 1 \
            "${args[@]}" 2>&1
    )" || rc=$?

    cmd="$(grep -m1 'submit_trim_fastqs.sh --env' <<< "${out}" || true)"
    lbl="execute ${lbl}"
    opt="${acts[1]}"

    case "${acts[0]}" in
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
                && "${cmd}" == *" ${opt}"*
                && "${out}" != *"has no effect"*
            ]]; then
                record_pass "${lbl} is accepted with no message"
            else
                record_fail "${lbl} was dropped or flagged (exit ${rc})"
            fi
            ;;

        drop)
            if [[
                "${rc}" -eq 0
                && -n "${cmd}"
                && "${cmd}" != *" ${opt}"*
                && "${out}" != *"has no effect"*
            ]]; then
                record_pass "${lbl} is not passed on, with no message"
            else
                record_fail "${lbl} was passed on or flagged (exit ${rc})"
            fi
            ;;
    esac
done


# Run per entry, each command passes on only the suffix for its own layout.
dir_out="${tmp}/exe_sfx_entry"
mkdir -p "${dir_out}/logs"

rc=0
out="$(
    "${TEST_BASH}" "${ROOT_REPO}/bin/execute_trim_fastqs.sh" \
        --dry_run \
        --env_nam "${env_nam}" \
        --threads 1 \
        --csv_fil_in "${in_se};${in_pe}" \
        --dir_out "${dir_out}" \
        --dir_eo "${dir_out}/logs" \
        --max_job 2 \
        --sfx_se "${sfx_se}" \
        --sfx_pe "${sfx_pe}" 2>&1
)" || rc=$?

pat="submit_trim_fastqs.sh --env"
cmd_se="$(grep -m1 "${pat}.*tiny_se\." <<< "${out}" || true)"
cmd_pe="$(grep -m1 "${pat}.*tiny_pe_R1" <<< "${out}" || true)"

if [[
    "${rc}" -eq 0
    && "${cmd_se}" == *" --sfx_se "* && "${cmd_se}" != *" --sfx_pe "*
    && "${cmd_pe}" == *" --sfx_pe "* && "${cmd_pe}" != *" --sfx_se "*
    && "${out}" != *"has no effect"*
]]; then
    record_pass "execute passes each entry only its own suffix"
else
    record_fail \
        "execute passed an entry the other layout's suffix (exit ${rc})"
fi


# 'submit' rows trim for real, and each is one of three kinds. Refused rows
# leave out the suffix for a layout an entry has: they must fail before
# anything is trimmed. Warned rows add a suffix for a layout no entry has: they
# must warn once and still trim. Control rows give only the needed suffix: they
# must trim with no message. Refused rows need no Atria, so they always run;
# the others run only with 'RUN_ATRIA=1'.
#
# Row format: label|input ('se' or 'pe')|arguments|act|message
o_se="--sfx_se ${sfx_se}"
o_pe="--sfx_pe ${sfx_pe}"
rows_sfx=(
    "no --sfx_se with an SE entry|se|${o_pe}|refuse|'sfx_se' is empty"
    "no --sfx_pe with a PE entry|pe|${o_se}|refuse|'sfx_pe' is empty"
)

if is_atria_enabled; then
    rows_sfx+=(
        "--sfx_pe with SE only|se|${o_se} ${o_pe}|warn|${msg_spe}"
        "--sfx_se with PE only|pe|${o_se} ${o_pe}|warn|${msg_sse}"
        "only --sfx_se with SE only|se|${o_se}|pass|"
        "only --sfx_pe with PE only|pe|${o_pe}|pass|"
    )
else
    record_skip \
        "Atria submit suffix checks disabled; set RUN_ATRIA=1 to enable"
fi

for idx in "${!rows_sfx[@]}"; do
    IFS='|' read -r lbl inp arg act msg <<< "${rows_sfx[idx]}"
    csv="${in_se}"
    stm="tiny_se"
    if [[ "${inp}" == "pe" ]]; then csv="${in_pe}"; stm="tiny_pe"; fi
    read -r -a args <<< "${arg}"
    dir_out="${tmp}/sub_sfx_${idx}"
    lbl="submit ${lbl}"

    mkdir -p "${dir_out}/logs"

    rc=0
    out="$(
        "${TEST_BASH}" "${ROOT_REPO}/bin/submit_trim_fastqs.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/bin" \
            --threads 1 \
            --csv_fil_in "${csv}" \
            --dir_out "${dir_out}" \
            --dir_eo "${dir_out}/logs" \
            "${args[@]}" 2>&1
    )" || rc=$?

    n_fq="$(
        find "${dir_out}" -maxdepth 1 -type f -name "${stm}*.fastq.gz" | wc -l
    )"

    case "${act}" in
        refuse)
            if [[
                "${rc}" -ne 0
                && "${out}" == *"error("*"${msg}"*
                && "${n_fq}" -eq 0
            ]]; then
                record_pass "${lbl} is refused before trimming"
            else
                record_fail "${lbl} was not refused (exit ${rc})"
            fi
            ;;

        warn)
            if [[
                "${rc}" -eq 0
                && "${out}" == *"warning(submit_trim_fastqs.sh"*"${msg}"*
                && "${out}" != *"::check_string_fastqs): "*
                && "${n_fq}" -gt 0
            ]]; then
                record_pass "${lbl} warns once and trims"
            else
                record_fail \
                    "${lbl} did not warn once, or did not trim (exit ${rc})"
            fi
            ;;

        pass)
            if [[
                "${rc}" -eq 0
                && "${out}" != *"has no effect"*
                && "${n_fq}" -gt 0
            ]]; then
                record_pass "${lbl} trims with no message"
            else
                record_fail "${lbl} was flagged, or did not trim (exit ${rc})"
            fi
            ;;
    esac
done

finish
