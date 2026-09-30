#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_compute_signal_parallel.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6);
# - Anthropic Claude Code (Opus 5, Opus 5.5).
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="execute compute-signal GNU Parallel"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


if ! \
    is_parallel_enabled
then
    record_skip \
        "GNU Parallel compute-signal check disabled; set RUN_PARALLEL=1 to" \
        "enable"
    finish
    exit $?
fi

# Define fixture and output paths for lightweight execute-layer config checks.
dir_fx="${ROOT_REPO}/tests/fixtures/compute_signal/bedgraph"
fil_A="${dir_fx}/ratio_A.bdg"
fil_B="${dir_fx}/ratio_B.bdg"

tmp="${TEST_DIR_TMP}/execute_compute_signal_parallel"
dir_in="${tmp}/in"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/compute_signal"
fil_A_1="${dir_in}/ratio_A_1.bdg"
fil_A_2="${dir_in}/ratio_A_2.bdg"

cfg_parallel="${dir_err}/test_execute_compute_parallel.config_parallel.txt"
fil_out_wet_1="${dir_out}/exec_parallel_wet_ratio_A_1.bdg"
fil_out_wet_2="${dir_out}/exec_parallel_wet_ratio_A_2.bdg"

log_env_parallel="${dir_log}/execute_compute_signal_parallel_env.log"
log_dry_run="${dir_log}/execute_compute_signal_parallel_dry_run.log"
log_wet_run="${dir_log}/execute_compute_signal_parallel_wet_run.log"


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${dir_in}" "${dir_out}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty \
    "${fil_A}" \
    "${fil_B}" || {
    finish
    exit $?
}

cp "${fil_A}" "${fil_A_1}"
cp "${fil_A}" "${fil_A_2}"

require_files_nonempty \
    "${fil_A_1}" \
    "${fil_A_2}" || {
    finish
    exit $?
}

# GNU Parallel dry-run config should invoke non-executable submit scripts via
# Bash.
# shellcheck disable=SC2154
if ! \
    require_env_parallel \
        "${env_nam}" \
        "${log_env_parallel}"
then
    finish
    exit $?
fi

if \
    run_capture \
        "execute compute-signal GNU Parallel dry-run" \
        "${log_dry_run}" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/execute_compute_signal.sh" \
            --dry_run \
            --env_nam "${env_nam}" \
            --threads 2 \
            --mode ratio \
            --method linear \
            --csv_fil_A "${fil_A}" \
            --csv_fil_B "${fil_B}" \
            --dir_out "${dir_out}" \
            --typ_out bdg \
            --prefix exec_parallel \
            --eps 0 \
            --dp 3 \
            --dir_eo "${dir_err}" \
            --nam_job "test_execute_compute_parallel" \
            --max_job 2
then
    record_pass "execute_compute_signal.sh GNU Parallel dry-run exits 0"
else
    record_fail \
        "execute_compute_signal.sh GNU Parallel dry-run failed; see" \
        "$(print_relpath "${log_dry_run}")"
fi

assert_file_nonempty \
    "${cfg_parallel}" \
    "execute GNU Parallel config"

if [[ -s "${cfg_parallel}" ]]; then
    assert_pattern_found \
        "${cfg_parallel}" \
        "${TEST_BASH} ${ROOT_REPO}/bin/submit_compute_signal.sh" \
        "execute GNU Parallel config uses Bash-prefixed submit command"
fi


# Real GNU Parallel execution should produce one output per ratio job.
if \
    run_capture \
        "execute compute-signal GNU Parallel wet run" \
        "${log_wet_run}" \
        "${TEST_BASH}" "${ROOT_REPO}/bin/execute_compute_signal.sh" \
            --env_nam "${env_nam}" \
            --threads 2 \
            --mode ratio \
            --method linear \
            --csv_fil_A "${fil_A_1},${fil_A_2}" \
            --csv_fil_B "${fil_B},${fil_B}" \
            --dir_out "${dir_out}" \
            --typ_out bdg \
            --prefix exec_parallel_wet \
            --eps 0 \
            --dp 3 \
            --dir_eo "${dir_err}" \
            --nam_job "test_execute_compute_parallel_wet" \
            --max_job 2
then
    record_pass "execute_compute_signal.sh GNU Parallel wet run exits 0"
else
    record_fail \
        "execute_compute_signal.sh GNU Parallel wet run failed; see" \
        "$(print_relpath "${log_wet_run}")"
fi

assert_file_nonempty \
    "${fil_out_wet_1}" \
    "execute GNU Parallel wet output 1"

assert_file_nonempty \
    "${fil_out_wet_2}" \
    "execute GNU Parallel wet output 2"

for out in "${fil_out_wet_1}" "${fil_out_wet_2}"; do
    if [[ -s "${out}" ]]; then
        assert_pattern_found \
            "${out}" \
            $'^I\t0\t10\t2$' \
            "$(basename "${out}") has I:0-10 = 2"

        assert_pattern_found \
            "${out}" \
            $'^I\t40\t50\t4$' \
            "$(basename "${out}") has I:40-50 = 4"

        assert_pattern_found \
            "${out}" \
            $'^I\t60\t70\t0.333$' \
            "$(basename "${out}") has I:60-70 = 0.333"
    fi
done


# A mixed pseudocount list under GNU Parallel: each job must get its own
# element, and the literal pair must run untouched.
# shellcheck source=lib/bash/core/format_outputs.sh
source "${ROOT_REPO}/lib/bash/core/format_outputs.sh"

dir_cnt="${ROOT_REPO}/tests/fixtures/compute_signal/bedgraph/count"
trk_se="${dir_cnt}/tiny_se.bedGraph"
trk_pe="${dir_cnt}/tiny_pe.bedGraph"
dir_mix="${tmp}/mixed"
mkdir -p "${dir_mix}"

exp_mix="$(
    "${TEST_MANAGED_PYTHON}" \
        -m protocol_chipseq_signal_norm.cli.compute_pseudo \
        --method edger \
        --typ_sig count \
        --prior_count 1 \
        --fil_A "${trk_pe}" \
        --fil_B "${trk_se}" \
        --n_frg_A "$(< "$(derive_report_path "${trk_pe}" n_frg)")" \
        --n_frg_B "$(< "$(derive_report_path "${trk_se}" n_frg)")" \
        --n_ovlp_A "$(< "$(derive_report_path "${trk_pe}" n_ovlp)")" \
        --n_ovlp_B "$(< "$(derive_report_path "${trk_se}" n_ovlp)")" \
        2>/dev/null
)"

rc_mix=0
"${TEST_BASH}" "${ROOT_REPO}/bin/execute_compute_signal.sh" \
    --env_nam "${env_nam}" \
    --threads 2 \
    --mode ratio \
    --method log2 \
    --typ_sig count \
    --prior_count 1 \
    --csv_fil_A "${trk_se},${trk_pe}" \
    --csv_fil_B "${trk_pe},${trk_se}" \
    --dir_out "${dir_mix}" \
    --typ_out bedGraph \
    --prefix mix \
    --dir_eo "${dir_mix}" \
    --nam_job "test_execute_compute_parallel_mix" \
    --max_job 2 \
    --csv_scl_fct NA,NA \
    --csv_dep_min NA,NA \
    --csv_pseudo "1:1,edger" > /dev/null 2>&1 || rc_mix=$?

pfx_mix_log="${dir_mix}/test_execute_compute_parallel_mix"
log_mix_lit="${pfx_mix_log}.mix_tiny_se.stderr.txt"
log_mix_edg="${pfx_mix_log}.mix_tiny_pe.stderr.txt"

if [[
    "${rc_mix}" -eq 0
    && -n "${exp_mix}"
    && -s "${dir_mix}/mix_tiny_se.bedGraph"
    && -s "${dir_mix}/mix_tiny_pe.bedGraph"
]] \
    && grep -qF -- "--pseudo 1:1 " "${log_mix_lit}" 2>/dev/null \
    && grep -qF -- "--pseudo ${exp_mix} " "${log_mix_edg}" 2>/dev/null
then
    record_pass \
        "GNU Parallel with '1:1,edger' gives each pair its own element"
else
    record_fail \
        "GNU Parallel with '1:1,edger' did not give each pair its own" \
        "element (exit ${rc_mix}); see $(print_relpath "${dir_mix}")"
fi
unset dir_cnt trk_se trk_pe dir_mix exp_mix rc_mix pfx_mix_log log_mix_lit
unset log_mix_edg


finish
