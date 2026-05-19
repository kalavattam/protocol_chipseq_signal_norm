#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_filter_bams_cram_parallel.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="execute filter-bams CRAM GNU Parallel"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

rec_section "${TEST_NAME}"

if ! is_parallel_enabled; then
    rec_skip \
        "GNU Parallel filter-bams CRAM check disabled;" \
        "set RUN_PARALLEL=1 to enable"
    finish
    exit $?
fi

dir_fx="${ROOT_REPO}/tests/filter_bams/fixtures"
in_sam="${dir_fx}/sam/filter_sc_sp.sam"
ref_fa="${dir_fx}/reference/filter_sc_sp.fa"
ref_fai="${ref_fa}.fai"

tmp="${TEST_DIR_TMP}/execute_filter_bams_cram_parallel"
dir_in="${tmp}/input"
dir_out_sc="${tmp}/out/sc"
dir_out_sp="${tmp}/out/sp"
dir_err_sc="${tmp}/logs/sc"
dir_err_sp="${tmp}/logs/sp"
dir_log="${TEST_DIR_LOG}/filter_bams"
in_cram_1="${dir_in}/filter_sc_sp_parallel_1.cram"
in_cram_2="${dir_in}/filter_sc_sp_parallel_2.cram"
csv_infile="${in_cram_1},${in_cram_2}"
nam_job_sc="test_execute_filter_bams_cram_parallel_sc"
nam_job_sp="test_execute_filter_bams_cram_parallel_sp"
config_sc="${dir_err_sc}/${nam_job_sc}.config_parallel.txt"
config_sp="${dir_err_sp}/${nam_job_sp}.config_parallel.txt"

rm -rf "${tmp}"
mkdir -p "${dir_in}" "${dir_out_sc}" "${dir_out_sp}" \
    "${dir_err_sc}" "${dir_err_sp}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_exist "${in_sam}" "${ref_fa}" "${ref_fai}" || {
    finish
    exit $?
}

# shellcheck disable=SC2154
if ! \
    require_parallel_env \
        "${env_nam}" \
        "${dir_log}/execute_filter_bams_cram_parallel_env.log"
then
    finish
    exit $?
fi


#  Build two deterministic CRAM inputs with distinct basenames
log="${dir_log}/execute_filter_bams_cram_parallel_prepare_1.log"
if ! \
    prepare_filter_bams_cram_fixture \
        "${in_sam}" "${ref_fa}" "${in_cram_1}" "${log}" \
        "execute filter-bams GNU Parallel CRAM fixture 1"
then
    finish
    exit $?
fi

log="${dir_log}/execute_filter_bams_cram_parallel_prepare_2.log"
if ! \
    prepare_filter_bams_cram_fixture \
        "${in_sam}" "${ref_fa}" "${in_cram_2}" "${log}" \
        "execute filter-bams GNU Parallel CRAM fixture 2"
then
    finish
    exit $?
fi


#  Run execute_filter_bams.sh through local GNU Parallel for two CRAM inputs
function run_parallel_case() {
    local retain="${1:-}"
    local nam_job="${2:-}"
    local dir_out_lcl="${3:-}"
    local dir_err_lcl="${4:-}"
    local log_lcl="${5:-}"

    shift 5

    if \
        run_capture \
            "execute filter-bams CRAM GNU Parallel retain=${retain}" \
            "${log_lcl}" \
            "${TEST_BASH}" "${ROOT_REPO}/scripts/execute_filter_bams.sh" \
                --threads 2 \
                --csv_infile "${csv_infile}" \
                --dir_out "${dir_out_lcl}" \
                --out_ext cram \
                --retain "${retain}" \
                --ref_fa "${ref_fa}" \
                --err_out "${dir_err_lcl}" \
                --nam_job "${nam_job}" \
                --max_job 2 \
                "$@"
    then
        rec_pass \
            "execute_filter_bams.sh CRAM GNU Parallel retain=${retain} exits 0"
    else
        rec_fail \
            "execute_filter_bams.sh CRAM GNU Parallel retain=${retain} failed;" \
            "see $(rec_relpath "${log_lcl}")"
    fi
}


function assert_parallel_config() {
    local config="${1:-}"
    local label="${2:-execute filter-bams GNU Parallel config}"

    assert_file_nonempty "${config}" "${label}"

    if [[ -s "${config}" ]]; then
        assert_grep_pattern \
            "${config}" \
            "${TEST_BASH} ${ROOT_REPO}/scripts/submit_filter_bams.sh" \
            "${label} uses Bash-prefixed submit command"
    fi
}


function assert_parallel_sp_output() {
    local sample="${1:-}"
    local outfile="${dir_out_sp}/${sample}.sp.cram"

    assert_file_nonempty "${outfile}" \
        "execute GNU Parallel retain=sp ${sample} CRAM output"
    assert_cram_index "${outfile}" \
        "execute GNU Parallel retain=sp ${sample} CRAI index"
    assert_file_exists \
        "${dir_err_sp}/${nam_job_sp}.${sample}.stdout.txt" \
        "execute GNU Parallel retain=sp ${sample} submit stdout log"
    assert_file_exists \
        "${dir_err_sp}/${nam_job_sp}.${sample}.stderr.txt" \
        "execute GNU Parallel retain=sp ${sample} submit stderr log"

    if [[ -s "${outfile}" ]]; then
        if \
            run_capture \
                "quickcheck execute filter-bams GNU Parallel retain=sp ${sample}" \
                "${dir_log}/execute_filter_bams_cram_parallel_sp_${sample}_quickcheck.log" \
                run_samtools quickcheck "${outfile}"
        then
            rec_pass \
                "execute GNU Parallel retain=sp ${sample} CRAM passes samtools quickcheck"
        else
            rec_fail \
                "execute GNU Parallel retain=sp ${sample} CRAM fails samtools quickcheck"
        fi

        assert_cram_count "${outfile}" "${ref_fa}" SP_I 1 \
            "execute GNU Parallel retain=sp ${sample} keeps SP_I" \
            "${dir_out_sp}/${sample}.SP_I.count.txt"
        assert_cram_count "${outfile}" "${ref_fa}" SP_II_TG 1 \
            "execute GNU Parallel retain=sp ${sample} keeps SP_II_TG" \
            "${dir_out_sp}/${sample}.SP_II_TG.count.txt"
        assert_cram_count "${outfile}" "${ref_fa}" SP_MTR 1 \
            "execute GNU Parallel retain=sp ${sample} keeps SP_MTR" \
            "${dir_out_sp}/${sample}.SP_MTR.count.txt"
        assert_cram_count "${outfile}" "${ref_fa}" SP_Mito 1 \
            "execute GNU Parallel retain=sp ${sample} keeps SP_Mito" \
            "${dir_out_sp}/${sample}.SP_Mito.count.txt"
        assert_cram_count "${outfile}" "${ref_fa}" I 0 \
            "execute GNU Parallel retain=sp ${sample} omits S. cerevisiae chromosomes" \
            "${dir_out_sp}/${sample}.I.count.txt"
        assert_cram_count "${outfile}" "${ref_fa}" chrUn 0 \
            "execute GNU Parallel retain=sp ${sample} omits unrelated contigs" \
            "${dir_out_sp}/${sample}.chrUn.count.txt"
    fi
}


function assert_parallel_sc_output() {
    local sample="${1:-}"
    local outfile="${dir_out_sc}/${sample}.sc.cram"

    assert_file_nonempty "${outfile}" \
        "execute GNU Parallel retain=sc ${sample} CRAM output"
    assert_cram_index "${outfile}" \
        "execute GNU Parallel retain=sc ${sample} CRAI index"
    assert_file_exists \
        "${dir_err_sc}/${nam_job_sc}.${sample}.stdout.txt" \
        "execute GNU Parallel retain=sc ${sample} submit stdout log"
    assert_file_exists \
        "${dir_err_sc}/${nam_job_sc}.${sample}.stderr.txt" \
        "execute GNU Parallel retain=sc ${sample} submit stderr log"

    if [[ -s "${outfile}" ]]; then
        if \
            run_capture \
                "quickcheck execute filter-bams GNU Parallel retain=sc ${sample}" \
                "${dir_log}/execute_filter_bams_cram_parallel_sc_${sample}_quickcheck.log" \
                run_samtools quickcheck "${outfile}"
        then
            rec_pass \
                "execute GNU Parallel retain=sc ${sample} CRAM passes samtools quickcheck"
        else
            rec_fail \
                "execute GNU Parallel retain=sc ${sample} CRAM fails samtools quickcheck"
        fi

        assert_cram_count "${outfile}" "${ref_fa}" I 1 \
            "execute GNU Parallel retain=sc ${sample} keeps chromosome I" \
            "${dir_out_sc}/${sample}.I.count.txt"
        assert_cram_count "${outfile}" "${ref_fa}" Mito 0 \
            "execute GNU Parallel retain=sc ${sample} omits Mito without --mito" \
            "${dir_out_sc}/${sample}.Mito.count.txt"
        assert_cram_count "${outfile}" "${ref_fa}" SP_I 0 \
            "execute GNU Parallel retain=sc ${sample} omits S. pombe chromosomes" \
            "${dir_out_sc}/${sample}.SP_I.count.txt"
        assert_cram_count "${outfile}" "${ref_fa}" chrUn 0 \
            "execute GNU Parallel retain=sc ${sample} omits unrelated contigs" \
            "${dir_out_sc}/${sample}.chrUn.count.txt"
    fi
}

run_parallel_case \
    sp \
    "${nam_job_sp}" \
    "${dir_out_sp}" \
    "${dir_err_sp}" \
    "${dir_log}/execute_filter_bams_cram_parallel_sp.log" \
    --tg \
    --mtr \
    --mito

assert_parallel_config "${config_sp}" \
    "execute filter-bams GNU Parallel retain=sp config"
assert_parallel_sp_output "filter_sc_sp_parallel_1"
assert_parallel_sp_output "filter_sc_sp_parallel_2"

run_parallel_case \
    sc \
    "${nam_job_sc}" \
    "${dir_out_sc}" \
    "${dir_err_sc}" \
    "${dir_log}/execute_filter_bams_cram_parallel_sc.log"

assert_parallel_config "${config_sc}" \
    "execute filter-bams GNU Parallel retain=sc config"
assert_parallel_sc_output "filter_sc_sp_parallel_1"
assert_parallel_sc_output "filter_sc_sp_parallel_2"

finish
