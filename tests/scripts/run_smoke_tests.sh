#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: run_smoke_tests.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


#  Require Bash >= 4.4 before doing any work
if [[ -z "${BASH_VERSION:-}" ]]; then
    echo "error(run_smoke_tests.sh):" \
        "this test suite requires Bash >= 4.4." >&2
    exit 1
elif ((
    BASH_VERSINFO[0] < 4 || ( BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] < 4 )
)); then
    echo "error(run_smoke_tests.sh):" \
        "this test suite requires Bash >= 4.4; current version is" \
        "'${BASH_VERSION}'." >&2
    exit 1
fi

#  Run with unset-variable checking while still collecting failed groups
set -u

#  Set paths for the smoke-test runner
dir_run="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
dir_rep="$(cd "${dir_run}/../.." > /dev/null 2>&1 && pwd)"
dir_out="${dir_rep}/tests/output"
dir_log="${dir_out}/logs"

mkdir -p "${dir_log}"


#  Define smoke-test groups ===================================================
tests=(
    "${dir_run}/smoke/test_shell_syntax.sh"
    "${dir_run}/smoke/test_python_startup.sh"
    "${dir_run}/smoke/test_help_output.sh"
    "${dir_run}/smoke/test_startup_sources.sh"
    "${dir_run}/smoke/test_dry_run_commands.sh"
    "${dir_run}/smoke/test_submit_align_fastqs_bowtie2.sh"
    "${dir_run}/smoke/test_execute_align_fastqs_bowtie2.sh"
    "${dir_run}/smoke/test_submit_align_fastqs_bwa_mem.sh"
    "${dir_run}/smoke/test_execute_align_fastqs_bwa_mem.sh"
    "${dir_run}/smoke/test_submit_align_fastqs_bwa_mem2.sh"
    "${dir_run}/smoke/test_execute_align_fastqs_bwa_mem2.sh"
    "${dir_run}/smoke/test_submit_align_fastqs_bwa_aln.sh"
    "${dir_run}/smoke/test_execute_align_fastqs_bwa_aln.sh"
    "${dir_run}/smoke/test_execute_align_fastqs_parallel.sh"
    "${dir_run}/smoke/test_execute_align_fastqs_slurm.sh"
    "${dir_run}/smoke/test_submit_trim_fastqs_se.sh"
    "${dir_run}/smoke/test_submit_trim_fastqs_pe.sh"
    "${dir_run}/smoke/test_execute_trim_fastqs_se.sh"
    "${dir_run}/smoke/test_execute_trim_fastqs_pe.sh"
    "${dir_run}/smoke/test_execute_trim_fastqs_parallel.sh"
    "${dir_run}/smoke/test_submit_compute_signal_ratio.sh"
    "${dir_run}/smoke/test_execute_compute_signal_ratio.sh"
    "${dir_run}/smoke/test_execute_compute_signal_parallel.sh"
    "${dir_run}/smoke/test_execute_compute_signal_bam.sh"
    "${dir_run}/smoke/test_execute_compute_signal_cram.sh"
    "${dir_run}/smoke/test_submit_compute_signal_bam.sh"
    "${dir_run}/smoke/test_submit_compute_signal_cram.sh"
    "${dir_run}/smoke/test_execute_download_fastqs_se_local.sh"
    "${dir_run}/smoke/test_execute_download_fastqs_pe_local.sh"
    "${dir_run}/smoke/test_execute_download_fastqs_mixed_local.sh"
    "${dir_run}/smoke/test_submit_filter_bams_bam.sh"
    "${dir_run}/smoke/test_execute_filter_bams_bam.sh"
    "${dir_run}/smoke/test_help_style.sh"
)

n_pass=0
n_fail=0
n_total=0


#  Run smoke-test groups ======================================================
printf 'Running smoke tests from %s\n' "${dir_rep}"
printf 'Logs: %s\n\n\n' "${dir_log}"

for scr_tst in "${tests[@]}"; do
    n_total=$(( n_total + 1 ))
    name="$(basename "${scr_tst}" .sh)"
    log="${dir_log}/${name}.log"

    printf '%s\n' "==== ${name} ===="
    if "${BASH}" "${scr_tst}" > "${log}" 2>&1; then
        cat "${log}"
        printf 'GROUP PASS: %s\n\n\n' "${name}"
        n_pass=$(( n_pass + 1 ))
    else
        cat "${log}"
        printf 'GROUP FAIL: %s; see %s\n\n\n' "${name}" "${log}" >&2
        n_fail=$(( n_fail + 1 ))
    fi
done


#  Summarize results ==========================================================
printf '==== Smoke Test Summary ====\n'
printf 'Groups run:    %d\n' "${n_total}"
printf 'Groups passed: %d\n' "${n_pass}"
printf 'Groups failed: %d\n' "${n_fail}"

if (( n_fail > 0 )); then
    exit 1
fi

exit 0
