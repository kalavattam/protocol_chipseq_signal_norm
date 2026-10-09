#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_compute_signal_opt_applicability.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="compute-signal option applicability"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Define fixture and output paths. Each row passes one option to a mode it does
# not apply to, per 'HELP.PARAMETER.APPLICABILITY'; a warned option is not
# checked, so rows with a missing file or a malformed value still warn.
dir_fx="${ROOT_REPO}/tests/fixtures/compute_signal"
in_se="${dir_fx}/bam/se/tiny_se.bam"
fil_A="${dir_fx}/bedgraph/ratio_A.bdg.gz"
fil_B="${dir_fx}/bedgraph/ratio_B.bdg.gz"
ref_fa="${dir_fx}/reference/tiny.fa"

tmp="${TEST_DIR_TMP}/compute_signal_opt_applicability"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_err}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty \
    "${in_se}" \
    "${fil_A}" \
    "${fil_B}" \
    "${ref_fa}" || {
    finish
    exit $?
}


# Print the arguments that select a mode and its inputs.
function args_mode() {
    case "${1}" in
        signal) printf '%s\n' --mode signal --csv_fil_in "${in_se}" ;;
        coord)
            printf '%s\n' --mode coord --typ_out bed.gz --csv_fil_in "${in_se}"
            ;;
        ratio)
            printf '%s\n' \
                --mode ratio --csv_fil_A "${fil_A}" --csv_fil_B "${fil_B}"
            ;;
    esac
}


# 'execute' rows run as dry runs: mode|option|value|refuse or warn.
rows_exe=(
    "ratio|--csv_fil_in|${in_se}|refuse"
    "signal|--csv_fil_A|${fil_A}|refuse"
    "coord|--csv_fil_B|${fil_B}|refuse"
    "coord|--method|frag|refuse"
    "ratio|--siz_bin|20|refuse"
    "coord|--siz_bin|20|refuse"
    "coord|--csv_scl_fct|2|refuse"
    "ratio|--csv_usr_frg|20|refuse"
    "signal|--csv_dep_min|5|refuse"
    "coord|--csv_pseudo|1:1|refuse"
    "signal|--eps|0.5|refuse"
    "coord|--skip_00|pre_scale|refuse"
    "signal|--strict_bins||refuse"
    "coord|--drp_nan||refuse"
    "signal|--skp_pfx|I|refuse"
    "coord|--track||refuse"
    "coord|--typ_sig|count|refuse"
    "ratio|--ref_fa|${ref_fa}|warn"
    "signal|--chr_siz|${ref_fa}.fai|warn"
    "ratio|--engine|window|warn"
    "coord|--siz_win|20|warn"
    "signal|--siz_win|20|warn"
    "coord|--dp|3|warn"
    "ratio|--no_report||warn"
    "ratio|--ref_fa|${tmp}/missing.fa|warn"
    "signal|--chr_siz|${tmp}/missing.sizes|warn"
    "signal|--siz_win|abc|warn"
    "coord|--dp|abc|warn"
)

for row in "${rows_exe[@]}"; do
    IFS='|' read -r mode opt val act <<< "${row}"
    mapfile -t args < <(args_mode "${mode}")
    args+=( "${opt}" )
    if [[ -n "${val}" ]]; then args+=( "${val}" ); fi

    rc=0
    out="$(
        "${TEST_BASH}" "${ROOT_REPO}/bin/execute_compute_signal.sh" \
            --dry_run \
            "${args[@]}" \
            --dir_out "${dir_out}" \
            --dir_eo "${dir_err}" \
            --max_job 1 \
            --threads 1 2>&1
    )" || rc=$?

    cmd="$(grep -m1 'submit_compute_signal.sh --env' <<< "${out}" || true)"
    lbl="execute '${opt}' with '--mode ${mode}'"

    if [[ "${act}" == "refuse" ]]; then
        if [[
            "${rc}" -ne 0
            && "${out}" == *"error("*"'${opt}' is for '--mode "*
        ]]; then
            record_pass "${lbl} is refused by name"
        else
            record_fail "${lbl} was not refused by name (exit ${rc})"
        fi
    elif [[
        "${rc}" -eq 0
        && "${out}" == *"warning("*"'${opt}' has no effect with '--"*
        && -n "${cmd}"
        && "${cmd}" != *" ${opt} "*
    ]]; then
        record_pass "${lbl} warns and is not passed on"
    else
        record_fail "${lbl} did not warn, or was passed on (exit ${rc})"
    fi
done


# 'submit' rows: a refusal stops before any tool runs; a warning lets the run
# finish without passing the option to the tool.
function args_submit() {
    local mode="${1}"
    local nam="${2}"

    printf '%s\n' \
        --env_nam "${env_nam}" \
        --dir_scr "${ROOT_REPO}/bin" \
        --threads 1 \
        --dir_eo "${dir_err}" \
        --nam_job "opt_appl_${nam}"

    case "${mode}" in
        signal)
            printf '%s\n' \
                --mode signal --csv_fil_in "${in_se}" \
                --csv_fil_out "${dir_out}/${nam}.bedGraph.gz" --no_report
            ;;
        coord)
            printf '%s\n' \
                --mode coord --csv_fil_in "${in_se}" \
                --csv_fil_out "${dir_out}/${nam}.bed.gz"
            ;;
        ratio)
            printf '%s\n' \
                --mode ratio --csv_fil_A "${fil_A}" --csv_fil_B "${fil_B}" \
                --csv_fil_out "${dir_out}/${nam}.bedGraph.gz"
            ;;
    esac
}

rows_sub=(
    "ratio|--csv_fil_in|${in_se}|refuse"
    "signal|--csv_fil_A|${fil_A}|refuse"
    "coord|--csv_fil_B|${fil_B}|refuse"
    "coord|--method|frag|refuse"
    "ratio|--siz_bin|20|refuse"
    "coord|--siz_bin|20|refuse"
    "coord|--csv_scl_fct|2|refuse"
    "ratio|--csv_usr_frg|20|refuse"
    "signal|--csv_dep_min|5|refuse"
    "coord|--csv_pseudo|1:1|refuse"
    "signal|--eps|0.5|refuse"
    "coord|--skip_00|pre_scale|refuse"
    "signal|--strict_bins||refuse"
    "coord|--drp_nan||refuse"
    "signal|--skp_pfx|I|refuse"
    "signal|--track||refuse"
    "coord|--typ_sig|count|refuse"
    "ratio|--csv_report_n_frg|${dir_out}/x.n_frg.txt|refuse"
    "coord|--csv_report_n_ovlp|${dir_out}/x.n_ovlp.txt|refuse"
    "ratio|--ref_fa|${ref_fa}|warn"
    "signal|--chr_siz|${ref_fa}.fai|warn"
    "ratio|--engine|window|warn"
    "coord|--siz_win|20|warn"
    "signal|--siz_win|20|warn"
    "coord|--dp|3|warn"
    "ratio|--no_report||warn"
    "ratio|--ref_fa|${tmp}/missing.fa|warn"
    "signal|--chr_siz|${tmp}/missing.sizes|warn"
    "signal|--siz_win|abc|warn"
    "coord|--dp|abc|warn"
)

for row in "${rows_sub[@]}"; do
    IFS='|' read -r mode opt val act <<< "${row}"
    nam="sub_${mode}_${opt#--}"
    mapfile -t args < <(args_submit "${mode}" "${nam}")
    args+=( "${opt}" )
    if [[ -n "${val}" ]]; then args+=( "${val}" ); fi

    rc=0
    out="$(
        "${TEST_BASH}" "${ROOT_REPO}/bin/submit_compute_signal.sh" \
            "${args[@]}" 2>&1
    )" || rc=$?

    cmd="$(grep -m1 'run_py compute_signal' <<< "${out}" || true)"
    lbl="submit '${opt}' with '--mode ${mode}'"

    if [[ "${act}" == "refuse" ]]; then
        if [[
            "${rc}" -ne 0
            && "${out}" == *"error("*"'${opt}' is for '--mode "*
            && -z "${cmd}"
        ]]; then
            record_pass "${lbl} is refused by name before any tool runs"
        else
            record_fail "${lbl} was not refused by name (exit ${rc})"
        fi
    elif [[
        "${rc}" -eq 0
        && "${out}" == *"warning("*"'${opt}' has no effect with '--"*
        && -n "${cmd}"
        && "${cmd}" != *" ${opt} "*
    ]]; then
        record_pass "${lbl} warns and is not passed to the tool"
    else
        record_fail "${lbl} did not warn, or reached the tool (exit ${rc})"
    fi
done


# An explicit report list contradicts '--no_report' in 'submit'.
rc=0
out="$(
    "${TEST_BASH}" "${ROOT_REPO}/bin/submit_compute_signal.sh" \
        --env_nam "${env_nam}" \
        --dir_scr "${ROOT_REPO}/bin" \
        --threads 1 \
        --dir_eo "${dir_err}" \
        --mode signal \
        --csv_fil_in "${in_se}" \
        --csv_fil_out "${dir_out}/contra.bedGraph.gz" \
        --csv_report_n_frg "${dir_out}/contra.n_frg.txt" \
        --no_report 2>&1
)" || rc=$?

if [[
    "${rc}" -ne 0
    && "${out}" == *"'--no_report' and an explicit report list contradict"*
]]; then
    record_pass "submit refuses '--no_report' with an explicit report list"
else
    record_fail \
        "submit did not refuse '--no_report' with an explicit report list" \
        "(exit ${rc})"
fi

# A report-only run builds no track, so the track's settings are refused or
# warned about and never passed on, and a warned one's value is not checked;
# the empty row is the control.
rows_rpt=(
    "--method|frag|refuse"
    "--csv_scl_fct|2|refuse"
    "--engine|window|warn"
    "--siz_win|20|warn"
    "--dp|3|warn"
    "--engine|bogus|warn"
    "--siz_win|abc|warn"
    "--dp|abc|warn"
    "||none"
)

for row in "${rows_rpt[@]}"; do
    IFS='|' read -r opt val act <<< "${row}"
    args=()
    if [[ -n "${opt}" ]]; then args+=( "${opt}" "${val}" ); fi

    rc=0
    out="$(
        "${TEST_BASH}" "${ROOT_REPO}/bin/execute_compute_signal.sh" \
            --dry_run \
            --mode signal \
            --csv_fil_in "${in_se}" \
            --report_only \
            "${args[@]}" \
            --dir_out "${dir_out}" \
            --dir_eo "${dir_err}" \
            --max_job 1 \
            --threads 1 2>&1
    )" || rc=$?

    cmd="$(grep -m1 'submit_compute_signal.sh --env' <<< "${out}" || true)"
    lbl="execute '${opt:-nothing extra}${val:+ ${val}}' with '--report_only'"
    msg="'${opt}' has no effect with '--report_only'"

    case "${act}" in
        refuse)
            if [[
                "${rc}" -ne 0
                && "${out}" == *"error("*"'${opt}' is for track output"*
            ]]; then
                record_pass "${lbl} is refused by name"
            else
                record_fail "${lbl} was not refused by name (exit ${rc})"
            fi
            ;;
        *)
            if [[
                "${rc}" -eq 0
                && -n "${cmd}"
                && "${cmd}" != *" --method "*
                && "${cmd}" != *" --csv_scl_fct "*
                && "${cmd}" != *" --engine "*
                && "${cmd}" != *" --siz_win "*
                && "${cmd}" != *" --dp "*
                && (
                    "${act}" == "none"
                    && "${out}" != *"has no effect"*
                    || "${out}" == *"warning("*"${msg}"*
                )
            ]]; then
                record_pass "${lbl} passes on no track setting"
            else
                record_fail "${lbl} warned wrongly, or passed a setting on"
            fi
            ;;
    esac

    rc=0
    out="$(
        "${TEST_BASH}" "${ROOT_REPO}/bin/submit_compute_signal.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/bin" \
            --threads 1 \
            --dir_eo "${dir_err}" \
            --nam_job "opt_appl_rpt_${opt#--}" \
            --mode signal \
            --csv_fil_in "${in_se}" \
            --csv_report_n_frg "${dir_out}/rpt_${opt#--}.n_frg.txt" \
            --csv_report_n_ovlp "${dir_out}/rpt_${opt#--}.n_ovlp.txt" \
            "${args[@]}" 2>&1
    )" || rc=$?

    cmd="$(grep -m1 'run_py compute_signal' <<< "${out}" || true)"
    lbl="submit '${opt:-nothing extra}' without '--csv_fil_out'"
    msg="'${opt}' has no effect without '--csv_fil_out'"

    case "${act}" in
        refuse)
            if [[
                "${rc}" -ne 0
                && "${out}" == *"error("*"'${opt}' is for track output"*
                && -z "${cmd}"
            ]]; then
                record_pass "${lbl} is refused by name before any tool runs"
            else
                record_fail "${lbl} was not refused by name (exit ${rc})"
            fi
            ;;
        *)
            if [[
                "${rc}" -eq 0
                && -n "${cmd}"
                && "${cmd}" != *" --method "*
                && "${cmd}" != *" --scl_fct "*
                && "${cmd}" != *" --engine "*
                && "${cmd}" != *" --siz_win "*
                && "${cmd}" != *" --dp "*
                && (
                    "${act}" == "none"
                    && "${out}" != *"has no effect"*
                    || "${out}" == *"warning("*"${msg}"*
                )
            ]]; then
                record_pass "${lbl} passes the tool no track setting"
            else
                record_fail "${lbl} warned wrongly, or passed a setting on"
            fi
            ;;
    esac
done

# With no track to derive from, 'submit' writes only the report list given and
# refuses '--siz_bin', which then sizes nothing.
dir_one="${dir_out}/one_list"
mkdir -p "${dir_one}"

for siz in "" 20; do
    args=()
    if [[ -n "${siz}" ]]; then args+=( --siz_bin "${siz}" ); fi

    rc=0
    out="$(
        "${TEST_BASH}" "${ROOT_REPO}/bin/submit_compute_signal.sh" \
            --env_nam "${env_nam}" \
            --dir_scr "${ROOT_REPO}/bin" \
            --threads 1 \
            --dir_eo "${dir_err}" \
            --nam_job "opt_appl_one_list" \
            --mode signal \
            --csv_fil_in "${in_se}" \
            --csv_report_n_frg "${dir_one}/se.n_frg.txt" \
            "${args[@]}" 2>&1
    )" || rc=$?

    cmd="$(grep -m1 'run_py compute_signal' <<< "${out}" || true)"

    if [[ -z "${siz}" ]]; then
        if [[
            "${rc}" -eq 0
            && -s "${dir_one}/se.n_frg.txt"
            && ! -e "${dir_one}/se.n_ovlp.txt"
            && "${cmd}" != *" --report_n_ovlp"*
            && "${cmd}" != *" --siz_bin "*
        ]]; then
            record_pass "submit writes only the one report list given"
        else
            record_fail "submit with one report list failed (exit ${rc})"
        fi
    elif [[
        "${rc}" -ne 0
        && "${out}" == *"'--siz_bin' is for track output or"*
        && -z "${cmd}"
    ]]; then
        record_pass "submit refuses '--siz_bin' without an overlap report"
    else
        record_fail "submit took '--siz_bin' without an overlap report"
    fi
done

# Options in their own modes stay silent: no applicability message at all.
rc=0
out="$(
    "${TEST_BASH}" "${ROOT_REPO}/bin/execute_compute_signal.sh" \
        --dry_run \
        --mode signal \
        --csv_fil_in "${in_se}" \
        --engine window \
        --siz_win 20 \
        --dp 3 \
        --dir_out "${dir_out}" \
        --dir_eo "${dir_err}" \
        --max_job 1 \
        --threads 1 2>&1
)" || rc=$?

if [[
    "${rc}" -eq 0
    && "${out}" != *"has no effect"*
    && "${out}" == *"--siz_win 20"*
]]; then
    record_pass "execute passes '--siz_win' to the 'window' engine silently"
else
    record_fail "execute warned about, or dropped, an applicable option"
fi

# A time limit is for Slurm jobs only, so a local run warns about it.
rc=0
out="$(
    "${TEST_BASH}" "${ROOT_REPO}/bin/execute_compute_signal.sh" \
        --dry_run \
        --mode signal \
        --csv_fil_in "${in_se}" \
        --time 1:00:00 \
        --dir_out "${dir_out}" \
        --dir_eo "${dir_err}" \
        --max_job 1 \
        --threads 1 2>&1
)" || rc=$?

if [[
    "${rc}" -eq 0
    && "${out}" == *"'--time' has no effect without '--slurm'"*
]]; then
    record_pass "execute warns about '--time' without '--slurm'"
else
    record_fail "execute did not warn about '--time' without '--slurm'"
fi

# A whole-count 'edger' run must pass its tools nothing they ignore, such as
# the fragment counts only the fractional types read.
dir_edg="${tmp}/edger_count"
mkdir -p "${dir_edg}/logs"

rc=0
out="$(
    "${TEST_BASH}" "${ROOT_REPO}/bin/execute_compute_signal.sh" \
        --mode ratio \
        --csv_fil_A "${dir_fx}/bedgraph/count/tiny_se.bedGraph" \
        --csv_fil_B "${dir_fx}/bedgraph/count/tiny_pe.bedGraph" \
        --csv_pseudo edger \
        --typ_sig count \
        --dir_out "${dir_edg}" \
        --dir_eo "${dir_edg}/logs" \
        --max_job 1 \
        --threads 1 2>&1
)" || rc=$?

if [[
    "${rc}" -eq 0
    && "${out}" != *"has no effect"*
]] && ! grep -rq "has no effect" "${dir_edg}/logs"; then
    record_pass "a whole-count 'edger' ratio run passes nothing ignored"
else
    record_fail \
        "a whole-count 'edger' ratio run failed or passed an ignored option" \
        "(exit ${rc}); see $(print_relpath "${dir_edg}/logs")"
fi


unset row mode opt val act args rc out cmd lbl nam dir_edg
finish
