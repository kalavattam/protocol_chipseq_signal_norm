#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_calculate_scaling_factor_spike.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="execute calculate-scaling-factor spike"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


function run_case_spike() {
    local cas="${1:-}"
    local arr_cmd_nam="${2:-}"
    local row_0="${3:-}"
    local row_1="${4:-}"
    local method="${5:-}"
    local tail_0="${6:-}"
    local tail_1="${7:-}"
    shift 7 || true

    if [[ -n "${method}" ]]; then
        set -- --method "${method}" "$@"
    fi

    run_case_scaling_factor_execute \
        "${cas}" \
        spike \
        "${arr_cmd_nam}" \
        "${dir_out}" \
        "${dir_log}" \
        spike \
        test_execute_calculate_scaling_factor_spike \
        $'^main_ip\tspike_ip\tmain_in\tspike_in\tspike\tcoef\tnum_mp\tnum_sp\tnum_mn\tnum_sn$' \
        "${row_0}" \
        "${row_1}" \
        "${tail_0}" \
        "${tail_1}" \
        "$@"
}



#  Define fixture, output, and worker-input paths
scr_exe="${ROOT_REPO}/bin/execute_calculate_scaling_factor.sh"
dir_fix="${ROOT_REPO}/tests/fixtures/calculate_scaling_factor"
dir_bam="${dir_fix}/bam"
dir_bam_se="${dir_bam}/se"
dir_bam_pe="${dir_bam}/pe"
dir_cram_se="${dir_fix}/cram/se"
dir_cram_pe="${dir_fix}/cram/pe"
ref_fa="${dir_fix}/reference/tiny.fa"

tmp="${TEST_DIR_TMP}/execute_calculate_scaling_factor_spike"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/calculate_scaling_factor"

bam_se_mip_0="${dir_bam_se}/IP_WT_log_Brn1_rep1.sc.bam"
bam_se_mip_1="${dir_bam_se}/IP_WT_log_Brn1_rep2.sc.bam"
bam_se_min_0="${dir_bam_se}/in_WT_log_Brn1_rep1.sc.bam"
bam_se_min_1="${dir_bam_se}/in_WT_log_Brn1_rep2.sc.bam"
bam_se_sip_0="${dir_bam_se}/IP_WT_log_Brn1_rep1.sp.bam"
bam_se_sip_1="${dir_bam_se}/IP_WT_log_Brn1_rep2.sp.bam"
bam_se_sin_0="${dir_bam_se}/in_WT_log_Brn1_rep1.sp.bam"
bam_se_sin_1="${dir_bam_se}/in_WT_log_Brn1_rep2.sp.bam"

bai_se_mip_0="${bam_se_mip_0}.bai"
bai_se_mip_1="${bam_se_mip_1}.bai"
bai_se_min_0="${bam_se_min_0}.bai"
bai_se_min_1="${bam_se_min_1}.bai"
bai_se_sip_0="${bam_se_sip_0}.bai"
bai_se_sip_1="${bam_se_sip_1}.bai"
bai_se_sin_0="${bam_se_sin_0}.bai"
bai_se_sin_1="${bam_se_sin_1}.bai"

bam_pe_mip_0="${dir_bam_pe}/IP_WT_G1_Hho1_6336.sc.bam"
bam_pe_mip_1="${dir_bam_pe}/IP_WT_G1_Hho1_6337.sc.bam"
bam_pe_min_0="${dir_bam_pe}/in_WT_G1_Hho1_6336.sc.bam"
bam_pe_min_1="${dir_bam_pe}/in_WT_G1_Hho1_6337.sc.bam"
bam_pe_sip_0="${dir_bam_pe}/IP_WT_G1_Hho1_6336.sp.bam"
bam_pe_sip_1="${dir_bam_pe}/IP_WT_G1_Hho1_6337.sp.bam"
bam_pe_sin_0="${dir_bam_pe}/in_WT_G1_Hho1_6336.sp.bam"
bam_pe_sin_1="${dir_bam_pe}/in_WT_G1_Hho1_6337.sp.bam"

bai_pe_mip_0="${bam_pe_mip_0}.bai"
bai_pe_mip_1="${bam_pe_mip_1}.bai"
bai_pe_min_0="${bam_pe_min_0}.bai"
bai_pe_min_1="${bam_pe_min_1}.bai"
bai_pe_sip_0="${bam_pe_sip_0}.bai"
bai_pe_sip_1="${bam_pe_sip_1}.bai"
bai_pe_sin_0="${bam_pe_sin_0}.bai"
bai_pe_sin_1="${bam_pe_sin_1}.bai"

cram_se_mip_0="${dir_cram_se}/IP_WT_log_Brn1_rep1.sc.cram"
cram_se_mip_1="${dir_cram_se}/IP_WT_log_Brn1_rep2.sc.cram"
cram_se_min_0="${dir_cram_se}/in_WT_log_Brn1_rep1.sc.cram"
cram_se_min_1="${dir_cram_se}/in_WT_log_Brn1_rep2.sc.cram"
cram_se_sip_0="${dir_cram_se}/IP_WT_log_Brn1_rep1.sp.cram"
cram_se_sip_1="${dir_cram_se}/IP_WT_log_Brn1_rep2.sp.cram"
cram_se_sin_0="${dir_cram_se}/in_WT_log_Brn1_rep1.sp.cram"
cram_se_sin_1="${dir_cram_se}/in_WT_log_Brn1_rep2.sp.cram"

crai_se_mip_0="${cram_se_mip_0}.crai"
crai_se_mip_1="${cram_se_mip_1}.crai"
crai_se_min_0="${cram_se_min_0}.crai"
crai_se_min_1="${cram_se_min_1}.crai"
crai_se_sip_0="${cram_se_sip_0}.crai"
crai_se_sip_1="${cram_se_sip_1}.crai"
crai_se_sin_0="${cram_se_sin_0}.crai"
crai_se_sin_1="${cram_se_sin_1}.crai"

cram_pe_mip_0="${dir_cram_pe}/IP_WT_G1_Hho1_6336.sc.cram"
cram_pe_mip_1="${dir_cram_pe}/IP_WT_G1_Hho1_6337.sc.cram"
cram_pe_min_0="${dir_cram_pe}/in_WT_G1_Hho1_6336.sc.cram"
cram_pe_min_1="${dir_cram_pe}/in_WT_G1_Hho1_6337.sc.cram"
cram_pe_sip_0="${dir_cram_pe}/IP_WT_G1_Hho1_6336.sp.cram"
cram_pe_sip_1="${dir_cram_pe}/IP_WT_G1_Hho1_6337.sp.cram"
cram_pe_sin_0="${dir_cram_pe}/in_WT_G1_Hho1_6336.sp.cram"
cram_pe_sin_1="${dir_cram_pe}/in_WT_G1_Hho1_6337.sp.cram"

crai_pe_mip_0="${cram_pe_mip_0}.crai"
crai_pe_mip_1="${cram_pe_mip_1}.crai"
crai_pe_min_0="${cram_pe_min_0}.crai"
crai_pe_min_1="${cram_pe_min_1}.crai"
crai_pe_sip_0="${cram_pe_sip_0}.crai"
crai_pe_sip_1="${cram_pe_sip_1}.crai"
crai_pe_sin_0="${cram_pe_sin_0}.crai"
crai_pe_sin_1="${cram_pe_sin_1}.crai"

fil_out_dry="${dir_out}/scaling.dry.spike.tsv"
fil_out_default="${dir_out}/scaling.default.spike.tsv"
fil_out_cram_missing_ref="${dir_out}/scaling.cram_missing_ref.spike.tsv"
fil_out_invalid_method="${dir_out}/scaling.invalid_method.spike.tsv"

prt_default_0="${fil_out_default}.part.000000"
prt_default_1="${fil_out_default}.part.000001"

nam_job_dry="test_execute_calculate_scaling_factor_spike_dry"
nam_job_default="test_execute_calculate_scaling_factor_spike_default"
nam_job_cram_missing_ref="test_execute_calculate_scaling_factor_spike_cram_missing_ref"
nam_job_invalid_method="test_execute_calculate_scaling_factor_spike_invalid_method"

log_dry="${dir_log}/execute_spike_dry_run.log"
log_cram_missing_ref="${dir_log}/execute_spike_cram_missing_ref.log"
log_invalid_method="${dir_log}/execute_spike_invalid_method.log"
log_existing="${dir_log}/execute_spike_existing.log"
log_force_no_parts="${dir_log}/execute_spike_force_no_parts.log"

log_out_default_0="${dir_err}/${nam_job_default}.IP_WT_log_Brn1_rep1.sc.stdout.txt"
log_err_default_0="${dir_err}/${nam_job_default}.IP_WT_log_Brn1_rep1.sc.stderr.txt"
log_out_default_1="${dir_err}/${nam_job_default}.IP_WT_log_Brn1_rep2.sc.stdout.txt"
log_err_default_1="${dir_err}/${nam_job_default}.IP_WT_log_Brn1_rep2.sc.stderr.txt"

log_err_pe_bam_0="${dir_err}/test_execute_calculate_scaling_factor_spike_pe_bam.IP_WT_G1_Hho1_6336.sc.stderr.txt"
log_err_pe_bam_1="${dir_err}/test_execute_calculate_scaling_factor_spike_pe_bam.IP_WT_G1_Hho1_6337.sc.stderr.txt"

log_err_se_cram_0="${dir_err}/test_execute_calculate_scaling_factor_spike_se_cram.IP_WT_log_Brn1_rep1.sc.stderr.txt"

log_err_pe_cram_0="${dir_err}/test_execute_calculate_scaling_factor_spike_pe_cram.IP_WT_G1_Hho1_6336.sc.stderr.txt"
log_err_pe_cram_1="${dir_err}/test_execute_calculate_scaling_factor_spike_pe_cram.IP_WT_G1_Hho1_6337.sc.stderr.txt"

log_err_mixed_layout_0="${dir_err}/test_execute_calculate_scaling_factor_spike_mixed_layout.IP_WT_log_Brn1_rep1.sc.stderr.txt"
log_err_mixed_layout_1="${dir_err}/test_execute_calculate_scaling_factor_spike_mixed_layout.IP_WT_G1_Hho1_6337.sc.stderr.txt"
log_err_mixed_format_1="${dir_err}/test_execute_calculate_scaling_factor_spike_mixed_format.IP_WT_G1_Hho1_6337.sc.stderr.txt"

row_bam_se_0="${bam_se_mip_0}"$'\t'"${bam_se_sip_0}"$'\t'"${bam_se_min_0}"$'\t'"${bam_se_sin_0}"
row_bam_se_1="${bam_se_mip_1}"$'\t'"${bam_se_sip_1}"$'\t'"${bam_se_min_1}"$'\t'"${bam_se_sin_1}"

row_bam_pe_0="${bam_pe_mip_0}"$'\t'"${bam_pe_sip_0}"$'\t'"${bam_pe_min_0}"$'\t'"${bam_pe_sin_0}"
row_bam_pe_1="${bam_pe_mip_1}"$'\t'"${bam_pe_sip_1}"$'\t'"${bam_pe_min_1}"$'\t'"${bam_pe_sin_1}"

row_cram_se_0="${cram_se_mip_0}"$'\t'"${cram_se_sip_0}"$'\t'"${cram_se_min_0}"$'\t'"${cram_se_sin_0}"
row_cram_se_1="${cram_se_mip_1}"$'\t'"${cram_se_sip_1}"$'\t'"${cram_se_min_1}"$'\t'"${cram_se_sin_1}"

row_cram_pe_0="${cram_pe_mip_0}"$'\t'"${cram_pe_sip_0}"$'\t'"${cram_pe_min_0}"$'\t'"${cram_pe_sin_0}"
row_cram_pe_1="${cram_pe_mip_1}"$'\t'"${cram_pe_sip_1}"$'\t'"${cram_pe_min_1}"$'\t'"${cram_pe_sin_1}"

row_mix_lyt_bam_0="${row_bam_se_0}"
row_mix_lyt_bam_1="${row_bam_pe_1}"

row_mix_fmt_bam_se_cram_pe_0="${row_bam_se_0}"
row_mix_fmt_bam_se_cram_pe_1="${cram_pe_mip_1}"$'\t'"${cram_pe_sip_1}"$'\t'"${cram_pe_min_1}"$'\t'"${cram_pe_sin_1}"

# shellcheck disable=2054
arr_cmd_bam_se=(
    "${TEST_BASH}" "${scr_exe}"
        --threads 1
        --mode spike
        --aln_typ auto
        --csv_mip "${bam_se_mip_0},${bam_se_mip_1}"
        --csv_min "${bam_se_min_0},${bam_se_min_1}"
        --csv_sip "${bam_se_sip_0},${bam_se_sip_1}"
        --csv_sin "${bam_se_sin_0},${bam_se_sin_1}"
        --dir_eo "${dir_err}"
        --max_job 1
)

# shellcheck disable=2034
{
    arr_cmd_bam_pe=(
        "${TEST_BASH}" "${scr_exe}"
            --threads 1
            --mode spike
            --aln_typ auto
            --csv_mip "${bam_pe_mip_0},${bam_pe_mip_1}"
            --csv_min "${bam_pe_min_0},${bam_pe_min_1}"
            --csv_sip "${bam_pe_sip_0},${bam_pe_sip_1}"
            --csv_sin "${bam_pe_sin_0},${bam_pe_sin_1}"
            --dir_eo "${dir_err}"
            --max_job 1
    )

    arr_cmd_cram_se=(
        "${TEST_BASH}" "${scr_exe}"
            --threads 1
            --mode spike
            --aln_typ auto
            --ref_fa "${ref_fa}"
            --csv_mip "${cram_se_mip_0},${cram_se_mip_1}"
            --csv_min "${cram_se_min_0},${cram_se_min_1}"
            --csv_sip "${cram_se_sip_0},${cram_se_sip_1}"
            --csv_sin "${cram_se_sin_0},${cram_se_sin_1}"
            --dir_eo "${dir_err}"
            --max_job 1
    )

    arr_cmd_cram_pe=(
        "${TEST_BASH}" "${scr_exe}"
            --threads 1
            --mode spike
            --aln_typ auto
            --ref_fa "${ref_fa}"
            --csv_mip "${cram_pe_mip_0},${cram_pe_mip_1}"
            --csv_min "${cram_pe_min_0},${cram_pe_min_1}"
            --csv_sip "${cram_pe_sip_0},${cram_pe_sip_1}"
            --csv_sin "${cram_pe_sin_0},${cram_pe_sin_1}"
            --dir_eo "${dir_err}"
            --max_job 1
    )

    arr_cmd_mix_lyt_bam=(
        "${TEST_BASH}" "${scr_exe}"
            --threads 1
            --mode spike
            --aln_typ auto
            --csv_mip "${bam_se_mip_0},${bam_pe_mip_1}"
            --csv_min "${bam_se_min_0},${bam_pe_min_1}"
            --csv_sip "${bam_se_sip_0},${bam_pe_sip_1}"
            --csv_sin "${bam_se_sin_0},${bam_pe_sin_1}"
            --dir_eo "${dir_err}"
            --max_job 1
    )

    arr_cmd_mix_fmt_bam_se_cram_pe=(
        "${TEST_BASH}" "${scr_exe}"
            --threads 1
            --mode spike
            --aln_typ auto
            --ref_fa "${ref_fa}"
            --csv_mip "${bam_se_mip_0},${cram_pe_mip_1}"
            --csv_min "${bam_se_min_0},${cram_pe_min_1}"
            --csv_sip "${bam_se_sip_0},${cram_pe_sip_1}"
            --csv_sin "${bam_se_sin_0},${cram_pe_sin_1}"
            --dir_eo "${dir_err}"
            --max_job 1
    )
}


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

if ! \
    require_files_nonempty \
        "${bam_se_mip_0}" "${bam_se_mip_1}" \
        "${bai_se_mip_0}" "${bai_se_mip_1}" \
        "${bam_se_min_0}" "${bam_se_min_1}" \
        "${bai_se_min_0}" "${bai_se_min_1}" \
        "${bam_se_sip_0}" "${bam_se_sip_1}" \
        "${bai_se_sip_0}" "${bai_se_sip_1}" \
        "${bam_se_sin_0}" "${bam_se_sin_1}" \
        "${bai_se_sin_0}" "${bai_se_sin_1}" \
        "${bam_pe_mip_0}" "${bam_pe_mip_1}" \
        "${bai_pe_mip_0}" "${bai_pe_mip_1}" \
        "${bam_pe_min_0}" "${bam_pe_min_1}" \
        "${bai_pe_min_0}" "${bai_pe_min_1}" \
        "${bam_pe_sip_0}" "${bam_pe_sip_1}" \
        "${bai_pe_sip_0}" "${bai_pe_sip_1}" \
        "${bam_pe_sin_0}" "${bam_pe_sin_1}" \
        "${bai_pe_sin_0}" "${bai_pe_sin_1}" \
        "${cram_se_mip_0}" "${cram_se_mip_1}" \
        "${crai_se_mip_0}" "${crai_se_mip_1}" \
        "${cram_se_min_0}" "${cram_se_min_1}" \
        "${crai_se_min_0}" "${crai_se_min_1}" \
        "${cram_se_sip_0}" "${cram_se_sip_1}" \
        "${crai_se_sip_0}" "${crai_se_sip_1}" \
        "${cram_se_sin_0}" "${cram_se_sin_1}" \
        "${crai_se_sin_0}" "${crai_se_sin_1}" \
        "${cram_pe_mip_0}" "${cram_pe_mip_1}" \
        "${crai_pe_mip_0}" "${crai_pe_mip_1}" \
        "${cram_pe_min_0}" "${cram_pe_min_1}" \
        "${crai_pe_min_0}" "${crai_pe_min_1}" \
        "${cram_pe_sip_0}" "${cram_pe_sip_1}" \
        "${crai_pe_sip_0}" "${crai_pe_sip_1}" \
        "${cram_pe_sin_0}" "${cram_pe_sin_1}" \
        "${crai_pe_sin_0}" "${crai_pe_sin_1}" \
        "${ref_fa}"
then
    finish
    exit $?
fi


#  Run one serial spike-in calculation case and assert its assembled table
#  Dry-run execution should report assembly without writing a final table
if \
    run_capture \
        "execute calculate-scaling-factor spike dry-run" \
        "${log_dry}" \
        "${arr_cmd_bam_se[@]}" \
        --env_nam "${env_nam}" \
        --fil_out "${fil_out_dry}" \
        --nam_job "${nam_job_dry}" \
        --dry_run
then
    record_pass "execute_calculate_scaling_factor.sh --dry_run exits 0"
else
    record_fail \
        "execute_calculate_scaling_factor.sh --dry_run failed; see" \
        "$(print_relpath "${log_dry}")"
fi

if [[ ! -e "${fil_out_dry}" ]]; then
    record_pass "execute_calculate_scaling_factor.sh --dry_run writes no TSV"
else
    record_fail "execute_calculate_scaling_factor.sh --dry_run wrote TSV"
fi

assert_pattern_found \
    "${log_dry}" \
    'combine_parts_scaling_factor.sh' \
    "execute_calculate_scaling_factor.sh --dry_run reports combiner command"

assert_pattern_found \
    "${log_dry}" \
    'write_header.sh' \
    "execute_calculate_scaling_factor.sh --dry_run reports header command"

assert_pattern_found \
    "${log_dry}" \
    'write_header.sh.*--fil_in.*--in_place' \
    "execute_calculate_scaling_factor.sh --dry_run reports in-place headering"

assert_pattern_found \
    "${log_dry}" \
    "${TEST_BASH} ${ROOT_REPO}/bin/submit_calculate_scaling_factor.sh" \
    "execute_calculate_scaling_factor.sh --dry_run reports Bash-prefixed" \
    "submit command"

assert_pattern_found \
    "${log_dry}" \
    '--method chiprx_alpha_ratio' \
    "execute_calculate_scaling_factor.sh --dry_run resolves default method"

assert_pattern_found \
    "${log_dry}" \
    '--idx_out 0' \
    "execute_calculate_scaling_factor.sh --dry_run reports first part index"

assert_pattern_found \
    "${log_dry}" \
    '--idx_out 1' \
    "execute_calculate_scaling_factor.sh --dry_run reports second part index"


#  Default method should use alignment-derived counts and automatic SE detection
run_case_spike \
    default \
    arr_cmd_bam_se \
    "${row_bam_se_0}" \
    "${row_bam_se_1}" \
    "" \
    $'2\tchiprx_alpha_ratio\t3\t1\t2\t2' \
    $'0.5\tchiprx_alpha_ratio\t2\t2\t3\t1'

run_case_spike \
    no_header \
    arr_cmd_bam_se \
    "${row_bam_se_0}" \
    "${row_bam_se_1}" \
    "" \
    $'2\tchiprx_alpha_ratio\t3\t1\t2\t2' \
    $'0.5\tchiprx_alpha_ratio\t2\t2\t3\t1' \
    --no_header

assert_file_exists \
    "${log_out_default_0}" \
    "execute scaling-factor spike default first submit stdout log exists"

assert_file_exists \
    "${log_err_default_0}" \
    "execute scaling-factor spike default first submit stderr log exists"

assert_file_exists \
    "${log_out_default_1}" \
    "execute scaling-factor spike default second submit stdout log exists"

assert_file_exists \
    "${log_err_default_1}" \
    "execute scaling-factor spike default second submit stderr log exists"

assert_pattern_found \
    "${log_err_default_0}" \
    'idx_out=0' \
    "execute scaling-factor spike default first submit uses idx_out=0"

assert_pattern_found \
    "${log_err_default_1}" \
    'idx_out=1' \
    "execute scaling-factor spike default second submit uses idx_out=1"

assert_pattern_found \
    "${log_err_default_0}" \
    'typ_mp=se' \
    "execute scaling-factor spike default auto-detects first main IP as SE"

assert_pattern_found \
    "${log_err_default_0}" \
    'num_mp=3' \
    "execute scaling-factor spike default counts first main IP alignments"

assert_pattern_found \
    "${log_err_default_0}" \
    'num_sp=1' \
    "execute scaling-factor spike default counts first spike IP alignments"

assert_pattern_found \
    "${log_err_default_0}" \
    'num_mn=2' \
    "execute scaling-factor spike default counts first main input alignments"

assert_pattern_found \
    "${log_err_default_0}" \
    'num_sn=2' \
    "execute scaling-factor spike default counts first spike input alignments"


#  Canonical spike-in methods and one accepted alias should propagate to TSVs
run_case_spike \
    fractional \
    arr_cmd_bam_se \
    "${row_bam_se_0}" \
    "${row_bam_se_1}" \
    fractional \
    $'2\tfractional\t3\t1\t2\t2' \
    $'0.5\tfractional\t2\t2\t3\t1'

run_case_spike \
    ratio_alias \
    arr_cmd_bam_se \
    "${row_bam_se_0}" \
    "${row_bam_se_1}" \
    chiprx_ratio \
    $'2\tchiprx_alpha_ratio\t3\t1\t2\t2' \
    $'0.5\tchiprx_alpha_ratio\t2\t2\t3\t1'

run_case_spike \
    alpha_ip \
    arr_cmd_bam_se \
    "${row_bam_se_0}" \
    "${row_bam_se_1}" \
    chiprx_alpha_ip \
    $'1000000\tchiprx_alpha_ip\t3\t1\t2\t2' \
    $'500000\tchiprx_alpha_ip\t2\t2\t3\t1'

run_case_spike \
    alpha_in \
    arr_cmd_bam_se \
    "${row_bam_se_0}" \
    "${row_bam_se_1}" \
    chiprx_alpha_in \
    $'500000\tchiprx_alpha_in\t3\t1\t2\t2' \
    $'1000000\tchiprx_alpha_in\t2\t2\t3\t1'

run_case_spike \
    rxinput \
    arr_cmd_bam_se \
    "${row_bam_se_0}" \
    "${row_bam_se_1}" \
    rxinput_alpha \
    $'500000\trxinput_alpha\t3\t1\t2\t2' \
    $'125000\trxinput_alpha\t2\t2\t3\t1'


#  Alignment variants should preserve counts under automatic type detection
run_case_spike \
    pe_bam \
    arr_cmd_bam_pe \
    "${row_bam_pe_0}" \
    "${row_bam_pe_1}" \
    "" \
    $'2\tchiprx_alpha_ratio\t3\t1\t2\t2' \
    $'0.5\tchiprx_alpha_ratio\t2\t2\t3\t1'

assert_pattern_found \
    "${log_err_pe_bam_0}" \
    'typ_mp=pe' \
    "execute scaling-factor spike PE BAM auto-detects first main IP as PE"

assert_pattern_found \
    "${log_err_pe_bam_1}" \
    'typ_mp=pe' \
    "execute scaling-factor spike PE BAM auto-detects second main IP as PE"

run_case_spike \
    se_cram \
    arr_cmd_cram_se \
    "${row_cram_se_0}" \
    "${row_cram_se_1}" \
    "" \
    $'2\tchiprx_alpha_ratio\t3\t1\t2\t2' \
    $'0.5\tchiprx_alpha_ratio\t2\t2\t3\t1'

assert_pattern_found \
    "${log_err_se_cram_0}" \
    'typ_mp=se' \
    "execute scaling-factor spike SE CRAM auto-detects first main IP as SE"

run_case_spike \
    pe_cram \
    arr_cmd_cram_pe \
    "${row_cram_pe_0}" \
    "${row_cram_pe_1}" \
    "" \
    $'2\tchiprx_alpha_ratio\t3\t1\t2\t2' \
    $'0.5\tchiprx_alpha_ratio\t2\t2\t3\t1'

assert_pattern_found \
    "${log_err_pe_cram_0}" \
    'typ_mp=pe' \
    "execute scaling-factor spike PE CRAM auto-detects first main IP as PE"

assert_pattern_found \
    "${log_err_pe_cram_1}" \
    'typ_mp=pe' \
    "execute scaling-factor spike PE CRAM auto-detects second main IP as PE"

run_case_spike \
    mixed_layout \
    arr_cmd_mix_lyt_bam \
    "${row_mix_lyt_bam_0}" \
    "${row_mix_lyt_bam_1}" \
    "" \
    $'2\tchiprx_alpha_ratio\t3\t1\t2\t2' \
    $'0.5\tchiprx_alpha_ratio\t2\t2\t3\t1'

assert_pattern_found \
    "${log_err_mixed_layout_0}" \
    'typ_mp=se' \
    "execute scaling-factor spike mixed layout keeps first main IP as SE"

assert_pattern_found \
    "${log_err_mixed_layout_1}" \
    'typ_mp=pe' \
    "execute scaling-factor spike mixed layout keeps second main IP as PE"

run_case_spike \
    mixed_format \
    arr_cmd_mix_fmt_bam_se_cram_pe \
    "${row_mix_fmt_bam_se_cram_pe_0}" \
    "${row_mix_fmt_bam_se_cram_pe_1}" \
    "" \
    $'2\tchiprx_alpha_ratio\t3\t1\t2\t2' \
    $'0.5\tchiprx_alpha_ratio\t2\t2\t3\t1'

assert_pattern_found \
    "${log_err_mixed_format_1}" \
    'typ_mp=pe' \
    "execute scaling-factor spike mixed BAM/CRAM list reads PE CRAM input"


#  CRAM input without a reference FASTA should fail before worker execution
if \
    run_capture \
        "execute calculate-scaling-factor spike CRAM missing reference" \
        "${log_cram_missing_ref}" \
        "${arr_cmd_bam_se[@]}" \
        --env_nam "${env_nam}" \
        --csv_mip "${cram_se_mip_0},${cram_se_mip_1}" \
        --csv_min "${cram_se_min_0},${cram_se_min_1}" \
        --csv_sip "${cram_se_sip_0},${cram_se_sip_1}" \
        --csv_sin "${cram_se_sin_0},${cram_se_sin_1}" \
        --fil_out "${fil_out_cram_missing_ref}" \
        --nam_job "${nam_job_cram_missing_ref}"
then
    record_fail \
        "execute_calculate_scaling_factor.sh CRAM without ref unexpectedly" \
        "pass"
else
    assert_pattern_found \
        "${log_cram_missing_ref}" \
        "'--ref_fa' is required" \
        "execute_calculate_scaling_factor.sh rejects CRAM without ref_fa"
fi


#  Broadcast depth overrides should replace per-file Samtools counts
run_case_spike \
    dep_broadcast \
    arr_cmd_bam_se \
    "${row_bam_se_0}" \
    "${row_bam_se_1}" \
    fractional \
    $'1\tfractional\t10\t2\t20\t4' \
    $'1\tfractional\t10\t2\t20\t4' \
    --csv_dep_mip 10 \
    --csv_dep_min 20 \
    --csv_dep_sip 2 \
    --csv_dep_sin 4


#  Unknown spike-in method should fail before running workers
if \
    run_capture \
        "execute calculate-scaling-factor spike invalid method" \
        "${log_invalid_method}" \
        "${arr_cmd_bam_se[@]}" \
        --env_nam "${env_nam}" \
        --method not_a_method \
        --fil_out "${fil_out_invalid_method}" \
        --nam_job "${nam_job_invalid_method}"
then
    record_fail \
        "execute_calculate_scaling_factor.sh invalid method unexpectedly pass"
else
    assert_pattern_found \
        "${log_invalid_method}" \
        'is not recognized' \
        "execute_calculate_scaling_factor.sh rejects invalid spike method"
fi


#  Existing final output should fail unless '--force' is supplied
if \
    run_capture \
        "execute calculate-scaling-factor spike existing output" \
        "${log_existing}" \
        "${arr_cmd_bam_se[@]}" \
        --env_nam "${env_nam}" \
        --fil_out "${fil_out_default}" \
        --nam_job "${nam_job_default}"
then
    record_fail \
        "execute_calculate_scaling_factor.sh existing output unexpectedly pass"
else
    assert_pattern_found \
        "${log_existing}" \
        'output file already exists' \
        "execute_calculate_scaling_factor.sh rejects existing output"
fi


#  '--force' should replace output and '--no_parts' should remove worker part
#+ files
if \
    run_capture \
        "execute calculate-scaling-factor spike force no-parts" \
        "${log_force_no_parts}" \
        "${arr_cmd_bam_se[@]}" \
        --env_nam "${env_nam}" \
        --fil_out "${fil_out_default}" \
        --nam_job "${nam_job_default}" \
        --force \
        --no_parts
then
    record_pass \
        "execute_calculate_scaling_factor.sh --force --no_parts exits 0"
else
    record_fail \
        "execute_calculate_scaling_factor.sh --force --no_parts failed; see" \
        "$(print_relpath "${log_force_no_parts}")"
fi

assert_file_nonempty \
    "${fil_out_default}" \
    "execute scaling-factor spike replaced final TSV"

if [[ ! -e "${prt_default_0}" && ! -e "${prt_default_1}" ]]; then
    record_pass "execute scaling-factor spike --no_parts removes worker parts"
else
    record_fail "execute scaling-factor spike --no_parts retained worker parts"
fi

finish
