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


#  Define functions ==========================================================
#  Print a smoke-test runner note
function runner_note() {
    echo "note(run_smoke_tests.sh):" "$@" >&2
}


#  Print a smoke-test runner error and stop
function runner_die() {
    echo "error(run_smoke_tests.sh):" "$@" >&2
    exit 1
}


#  Check whether each required fixture file exists and is non-empty
function has_fixture_set() {
    local arr_ref="${1:-}"
    local file=""

    if [[ -z "${arr_ref}" ]]; then
        runner_die "'has_fixture_set' requires an array reference."
    fi

    local -n arr_file="${arr_ref}"

    for file in "${arr_file[@]}"; do
        if [[ ! -s "${file}" ]]; then
            return 1
        fi
    done
}


#  Run one fixture generator
function run_fixture_generator() {
    local label="${1:-fixture}"
    local script="${2:-}"
    local env_req="${3:-}"
    local -a arr_cmd=()

    if [[ ! -r "${script}" ]]; then
        runner_die "${label} generator is not readable: '${script}'."
    fi

    runner_note "generating missing ${label} fixtures."

    if [[ -z "${env_req}" || "${CONDA_DEFAULT_ENV:-}" == "${env_req}" ]]; then
        arr_cmd=( "${BASH}" "${script}" )
    elif \
        command -v conda > /dev/null 2>&1
    then
        arr_cmd=( conda run -n "${env_req}" bash "${script}" )
    else
        runner_die \
            "${label} fixture generation requires active environment" \
            "'${env_req}' or Conda in PATH."
    fi

    if ! \
        "${arr_cmd[@]}"
    then
        runner_die "${label} fixture generation failed: '${script}'."
    fi
}


#  Generate one fixture set only when required files are missing
function ensure_fixture_set() {
    local label="${1:-fixture}"
    local script="${2:-}"
    local arr_ref="${3:-}"
    local env_req="${4:-}"
    local file=""

    if \
        has_fixture_set "${arr_ref}"
    then
        runner_note "using existing ${label} fixtures."
        return 0
    fi

    local -n arr_file="${arr_ref}"

    runner_note "missing required ${label} fixture files:"
    for file in "${arr_file[@]}"; do
        if [[ ! -s "${file}" ]]; then
            printf '  %s\n' "${file}" >&2
        fi
    done

    run_fixture_generator "${label}" "${script}" "${env_req}"

    if ! \
        has_fixture_set "${arr_ref}"
    then
        runner_die "${label} fixture generation did not create required files."
    fi

    runner_note "generated ${label} fixtures."
}


#  Set paths for the smoke-test runner ========================================
dir_run="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
dir_rep="$(cd "${dir_run}/../.." > /dev/null 2>&1 && pwd)"
dir_out="${dir_rep}/tests/outputs"
dir_log="${dir_out}/logs"

mkdir -p "${dir_log}"


#  Define fixture sets generated automatically when missing ==================
# shellcheck disable=SC2034
{
    dir_fix_aln="${dir_rep}/tests/align_fastqs/fixtures"
    arr_fix_aln=(
        "${dir_fix_aln}/reference/tiny.fa"
        "${dir_fix_aln}/fastq/se/tiny_se.atria.fastq.gz"
        "${dir_fix_aln}/fastq/pe/tiny_pe_R1.atria.fastq.gz"
        "${dir_fix_aln}/fastq/pe/tiny_pe_R2.atria.fastq.gz"
        "${dir_fix_aln}/bowtie2/tiny.1.bt2"
        "${dir_fix_aln}/bowtie2/tiny.2.bt2"
        "${dir_fix_aln}/bowtie2/tiny.3.bt2"
        "${dir_fix_aln}/bowtie2/tiny.4.bt2"
        "${dir_fix_aln}/bowtie2/tiny.rev.1.bt2"
        "${dir_fix_aln}/bowtie2/tiny.rev.2.bt2"
        "${dir_fix_aln}/bwa/tiny.fa"
        "${dir_fix_aln}/bwa/tiny.fa.amb"
        "${dir_fix_aln}/bwa/tiny.fa.ann"
        "${dir_fix_aln}/bwa/tiny.fa.bwt"
        "${dir_fix_aln}/bwa/tiny.fa.pac"
        "${dir_fix_aln}/bwa/tiny.fa.sa"
        "${dir_fix_aln}/bwa-mem2/tiny.fa"
        "${dir_fix_aln}/bwa-mem2/tiny.fa.0123"
        "${dir_fix_aln}/bwa-mem2/tiny.fa.amb"
        "${dir_fix_aln}/bwa-mem2/tiny.fa.ann"
        "${dir_fix_aln}/bwa-mem2/tiny.fa.bwt.2bit.64"
        "${dir_fix_aln}/bwa-mem2/tiny.fa.pac"
    )

    dir_fix_scl="${dir_rep}/tests/calculate_scaling_factor/fixtures"
    arr_fix_scl=(
        "${dir_fix_scl}/bam/se/IP_A.sc.bam"
        "${dir_fix_scl}/bam/se/IP_A.sc.bam.bai"
        "${dir_fix_scl}/bam/se/IP_B.sc.bam"
        "${dir_fix_scl}/bam/se/IP_B.sc.bam.bai"
        "${dir_fix_scl}/bam/se/in_A.sc.bam"
        "${dir_fix_scl}/bam/se/in_A.sc.bam.bai"
        "${dir_fix_scl}/bam/se/in_B.sc.bam"
        "${dir_fix_scl}/bam/se/in_B.sc.bam.bai"
        "${dir_fix_scl}/bam/se/IP_A.sp.bam"
        "${dir_fix_scl}/bam/se/IP_A.sp.bam.bai"
        "${dir_fix_scl}/bam/se/IP_B.sp.bam"
        "${dir_fix_scl}/bam/se/IP_B.sp.bam.bai"
        "${dir_fix_scl}/bam/se/in_A.sp.bam"
        "${dir_fix_scl}/bam/se/in_A.sp.bam.bai"
        "${dir_fix_scl}/bam/se/in_B.sp.bam"
        "${dir_fix_scl}/bam/se/in_B.sp.bam.bai"
        "${dir_fix_scl}/reference/tiny.fa"
        "${dir_fix_scl}/reference/tiny.fa.fai"
        "${dir_fix_scl}/bam/pe/IP_A.sc.bam"
        "${dir_fix_scl}/bam/pe/IP_A.sc.bam.bai"
        "${dir_fix_scl}/bam/pe/IP_B.sc.bam"
        "${dir_fix_scl}/bam/pe/IP_B.sc.bam.bai"
        "${dir_fix_scl}/bam/pe/in_A.sc.bam"
        "${dir_fix_scl}/bam/pe/in_A.sc.bam.bai"
        "${dir_fix_scl}/bam/pe/in_B.sc.bam"
        "${dir_fix_scl}/bam/pe/in_B.sc.bam.bai"
        "${dir_fix_scl}/bam/pe/IP_A.sp.bam"
        "${dir_fix_scl}/bam/pe/IP_A.sp.bam.bai"
        "${dir_fix_scl}/bam/pe/IP_B.sp.bam"
        "${dir_fix_scl}/bam/pe/IP_B.sp.bam.bai"
        "${dir_fix_scl}/bam/pe/in_A.sp.bam"
        "${dir_fix_scl}/bam/pe/in_A.sp.bam.bai"
        "${dir_fix_scl}/bam/pe/in_B.sp.bam"
        "${dir_fix_scl}/bam/pe/in_B.sp.bam.bai"
        "${dir_fix_scl}/cram/se/IP_A.sc.cram"
        "${dir_fix_scl}/cram/se/IP_B.sc.cram"
        "${dir_fix_scl}/cram/se/in_A.sc.cram"
        "${dir_fix_scl}/cram/se/in_B.sc.cram"
        "${dir_fix_scl}/cram/se/IP_A.sp.cram"
        "${dir_fix_scl}/cram/se/IP_B.sp.cram"
        "${dir_fix_scl}/cram/se/in_A.sp.cram"
        "${dir_fix_scl}/cram/se/in_B.sp.cram"
        "${dir_fix_scl}/cram/pe/IP_A.sc.cram"
        "${dir_fix_scl}/cram/pe/in_A.sc.cram"
        "${dir_fix_scl}/cram/pe/IP_A.sp.cram"
        "${dir_fix_scl}/cram/pe/in_A.sp.cram"
        "${dir_fix_scl}/cram/pe/IP_B.sc.cram"
        "${dir_fix_scl}/cram/pe/in_B.sc.cram"
        "${dir_fix_scl}/cram/pe/IP_B.sp.cram"
        "${dir_fix_scl}/cram/pe/in_B.sp.cram"
        "${dir_fix_scl}/parts/example_scaling_factors.spike.tsv.part.000000"
        "${dir_fix_scl}/parts/example_scaling_factors.spike.tsv.part.000002"
        "${dir_fix_scl}/parts/example_scaling_factors.siq.tsv.part.000000"
        "${dir_fix_scl}/parts/example_scaling_factors.siq.tsv.part.000002"
        "${dir_fix_scl}/parts/malformed_scaling_factors.spike.tsv.part.000003"
        "${dir_fix_scl}/parts/header_scaling_factors.spike.tsv.part.000004"
        "${dir_fix_scl}/parts/duplicate_index_A.spike.tsv.part.000005"
        "${dir_fix_scl}/parts/duplicate_index_B.spike.tsv.part.000005"
    )

    dir_fix_sig="${dir_rep}/tests/compute_signal/fixtures"
    arr_fix_sig=(
        "${dir_fix_sig}/bedgraph/ratio_A.bdg"
        "${dir_fix_sig}/bedgraph/ratio_B.bdg"
        "${dir_fix_sig}/bedgraph/ratio_headers_A.bdg"
        "${dir_fix_sig}/bedgraph/ratio_headers_B.bdg"
        "${dir_fix_sig}/reference/tiny.fa"
        "${dir_fix_sig}/reference/tiny.fa.fai"
        "${dir_fix_sig}/bam/se/tiny_se.bam"
        "${dir_fix_sig}/bam/se/tiny_se.bam.bai"
        "${dir_fix_sig}/bam/pe/tiny_pe.bam"
        "${dir_fix_sig}/bam/pe/tiny_pe.bam.bai"
        "${dir_fix_sig}/cram/se/tiny_se.cram"
        "${dir_fix_sig}/cram/se/tiny_se.cram.crai"
        "${dir_fix_sig}/cram/pe/tiny_pe.cram"
        "${dir_fix_sig}/cram/pe/tiny_pe.cram.crai"
    )

    dir_fix_dwn="${dir_rep}/tests/download_fastqs/fixtures"
    arr_fix_dwn=(
        "${dir_fix_dwn}/source/se/tiny_download_se.fastq.gz"
        "${dir_fix_dwn}/source/pe/tiny_download_pe_R1.fastq.gz"
        "${dir_fix_dwn}/source/pe/tiny_download_pe_R2.fastq.gz"
        "${dir_fix_dwn}/metadata/local_se.template.tsv"
        "${dir_fix_dwn}/metadata/local_pe.template.tsv"
        "${dir_fix_dwn}/metadata/local_mixed.template.tsv"
    )

    dir_fix_flt="${dir_rep}/tests/filter_alignments/fixtures"
    arr_fix_flt=(
        "${dir_fix_flt}/sam/filter_sc_sp.sam"
        "${dir_fix_flt}/reference/filter_sc_sp.fa"
        "${dir_fix_flt}/reference/filter_sc_sp.fa.fai"
    )

    dir_fix_trm="${dir_rep}/tests/trim_fastqs/fixtures"
    arr_fix_trm=(
        "${dir_fix_trm}/fastq/se/tiny_se.fastq.gz"
        "${dir_fix_trm}/fastq/pe/tiny_pe_R1.fastq.gz"
        "${dir_fix_trm}/fastq/pe/tiny_pe_R2.fastq.gz"
    )
}

#  Define smoke-test groups ===================================================
arr_tst=(
    "${dir_run}/smoke/test_shell_syntax.sh"
    "${dir_run}/smoke/test_python_startup.sh"
    "${dir_run}/smoke/test_help_output.sh"
    "${dir_run}/smoke/test_startup_sources.sh"
    "${dir_run}/smoke/test_dry_run_commands.sh"
    "${dir_run}/smoke/test_combine_parts_scaling_factor.sh"
    "${dir_run}/smoke/test_write_header.sh"
    "${dir_run}/smoke/test_submit_calculate_scaling_factor_spike.sh"
    "${dir_run}/smoke/test_execute_calculate_scaling_factor_spike.sh"
    "${dir_run}/smoke/test_clean_test_outputs.sh"
    "${dir_run}/smoke/test_install_envs_layout.sh"
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
    "${dir_run}/smoke/test_execute_download_fastqs_parallel.sh"
    "${dir_run}/smoke/test_submit_filter_alignments_bam.sh"
    "${dir_run}/smoke/test_execute_filter_alignments_bam.sh"
    "${dir_run}/smoke/test_submit_filter_alignments_cram_input_bam.sh"
    "${dir_run}/smoke/test_execute_filter_alignments_cram_input_bam.sh"
    "${dir_run}/smoke/test_submit_filter_alignments_bam_to_cram.sh"
    "${dir_run}/smoke/test_execute_filter_alignments_bam_to_cram.sh"
    "${dir_run}/smoke/test_submit_filter_alignments_cram_to_cram.sh"
    "${dir_run}/smoke/test_execute_filter_alignments_cram_to_cram.sh"
    "${dir_run}/smoke/test_execute_filter_alignments_cram_parallel.sh"
    "${dir_run}/smoke/test_help_style.sh"
)

n_pass=0
n_fail=0
n_totl=0


#  Generate missing fixture sets =============================================
ensure_fixture_set \
    "align-fastqs" \
    "${dir_rep}/tests/align_fastqs/scripts/make_fixtures.sh" \
    arr_fix_aln \
    env_protocol

ensure_fixture_set \
    "calculate-scaling-factor" \
    "${dir_rep}/tests/calculate_scaling_factor/scripts/make_fixtures.sh" \
    arr_fix_scl \
    env_protocol

ensure_fixture_set \
    "compute-signal" \
    "${dir_rep}/tests/compute_signal/scripts/make_fixtures.sh" \
    arr_fix_sig \
    env_protocol

ensure_fixture_set \
    "download-fastqs" \
    "${dir_rep}/tests/download_fastqs/scripts/make_fixtures.sh" \
    arr_fix_dwn

ensure_fixture_set \
    "filter-alignments" \
    "${dir_rep}/tests/filter_alignments/scripts/make_fixtures.sh" \
    arr_fix_flt

ensure_fixture_set \
    "trim-fastqs" \
    "${dir_rep}/tests/trim_fastqs/scripts/make_fixtures.sh" \
    arr_fix_trm


#  Run smoke-test groups ======================================================
printf 'Running smoke tests from %s\n' "${dir_rep}"
printf 'Logs: %s\n\n\n' "${dir_log}"

for scr_tst in "${arr_tst[@]}"; do
    n_totl=$(( n_totl + 1 ))
    name="$(basename "${scr_tst}" .sh)"
    log="${dir_log}/${name}.log"

    printf '%s\n' "==== ${name} ===="
    if \
        "${BASH}" "${scr_tst}" > "${log}" 2>&1
    then
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
printf 'Groups run:    %d\n' "${n_totl}"
printf 'Groups passed: %d\n' "${n_pass}"
printf 'Groups failed: %d\n' "${n_fail}"

if (( n_fail > 0 )); then
    exit 1
fi

exit 0
