#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_calculate_scaling_factor_spike.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="execute calculate-scaling-factor spike"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

print_section "${TEST_NAME}"

#  Define fixture, output, and worker-input paths
scr_exe="${ROOT_REPO}/scripts/execute_calculate_scaling_factor.sh"
dir_fix="${ROOT_REPO}/tests/calculate_scaling_factor/fixtures"
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

bam_pe_mip_0="${dir_bam_pe}/IP_WT_G1_Hho1_6336.sc.bam"
bam_pe_mip_1="${dir_bam_pe}/IP_WT_G1_Hho1_6337.sc.bam"
bam_pe_min_0="${dir_bam_pe}/in_WT_G1_Hho1_6336.sc.bam"
bam_pe_min_1="${dir_bam_pe}/in_WT_G1_Hho1_6337.sc.bam"
bam_pe_sip_0="${dir_bam_pe}/IP_WT_G1_Hho1_6336.sp.bam"
bam_pe_sip_1="${dir_bam_pe}/IP_WT_G1_Hho1_6337.sp.bam"
bam_pe_sin_0="${dir_bam_pe}/in_WT_G1_Hho1_6336.sp.bam"
bam_pe_sin_1="${dir_bam_pe}/in_WT_G1_Hho1_6337.sp.bam"

cram_se_mip_0="${dir_cram_se}/IP_WT_log_Brn1_rep1.sc.cram"
cram_se_mip_1="${dir_cram_se}/IP_WT_log_Brn1_rep2.sc.cram"
cram_se_min_0="${dir_cram_se}/in_WT_log_Brn1_rep1.sc.cram"
cram_se_min_1="${dir_cram_se}/in_WT_log_Brn1_rep2.sc.cram"
cram_se_sip_0="${dir_cram_se}/IP_WT_log_Brn1_rep1.sp.cram"
cram_se_sip_1="${dir_cram_se}/IP_WT_log_Brn1_rep2.sp.cram"
cram_se_sin_0="${dir_cram_se}/in_WT_log_Brn1_rep1.sp.cram"
cram_se_sin_1="${dir_cram_se}/in_WT_log_Brn1_rep2.sp.cram"

cram_pe_mip_0="${dir_cram_pe}/IP_WT_G1_Hho1_6336.sc.cram"
cram_pe_mip_1="${dir_cram_pe}/IP_WT_G1_Hho1_6337.sc.cram"
cram_pe_min_0="${dir_cram_pe}/in_WT_G1_Hho1_6336.sc.cram"
cram_pe_min_1="${dir_cram_pe}/in_WT_G1_Hho1_6337.sc.cram"
cram_pe_sip_0="${dir_cram_pe}/IP_WT_G1_Hho1_6336.sp.cram"
cram_pe_sip_1="${dir_cram_pe}/IP_WT_G1_Hho1_6337.sp.cram"
cram_pe_sin_0="${dir_cram_pe}/in_WT_G1_Hho1_6336.sp.cram"
cram_pe_sin_1="${dir_cram_pe}/in_WT_G1_Hho1_6337.sp.cram"

rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

if ! \
    require_files_nonempty \
        "${bam_se_mip_0}" "${bam_se_mip_1}" \
        "${bam_se_mip_0}.bai" "${bam_se_mip_1}.bai" \
        "${bam_se_min_0}" "${bam_se_min_1}" \
        "${bam_se_min_0}.bai" "${bam_se_min_1}.bai" \
        "${bam_se_sip_0}" "${bam_se_sip_1}" \
        "${bam_se_sip_0}.bai" "${bam_se_sip_1}.bai" \
        "${bam_se_sin_0}" "${bam_se_sin_1}" \
        "${bam_se_sin_0}.bai" "${bam_se_sin_1}.bai" \
        "${bam_pe_mip_0}" "${bam_pe_mip_1}" \
        "${bam_pe_mip_0}.bai" "${bam_pe_mip_1}.bai" \
        "${bam_pe_min_0}" "${bam_pe_min_1}" \
        "${bam_pe_min_0}.bai" "${bam_pe_min_1}.bai" \
        "${bam_pe_sip_0}" "${bam_pe_sip_1}" \
        "${bam_pe_sip_0}.bai" "${bam_pe_sip_1}.bai" \
        "${bam_pe_sin_0}" "${bam_pe_sin_1}" \
        "${bam_pe_sin_0}.bai" "${bam_pe_sin_1}.bai" \
        "${cram_se_mip_0}" "${cram_se_mip_1}" \
        "${cram_se_min_0}" "${cram_se_min_1}" \
        "${cram_se_sip_0}" "${cram_se_sip_1}" \
        "${cram_se_sin_0}" "${cram_se_sin_1}" \
        "${cram_pe_mip_0}" "${cram_pe_mip_1}" \
        "${cram_pe_min_0}" "${cram_pe_min_1}" \
        "${cram_pe_sip_0}" "${cram_pe_sip_1}" \
        "${cram_pe_sin_0}" "${cram_pe_sin_1}" \
        "${ref_fa}"
then
    finish
    exit $?
fi

# shellcheck disable=2054
arr_cmd_bam_se=(
    "${TEST_BASH}" "${scr_exe}"
        --threads 1
        --mode spike
        --csv_mip "${bam_se_mip_0},${bam_se_mip_1}"
        --csv_min "${bam_se_min_0},${bam_se_min_1}"
        --csv_sip "${bam_se_sip_0},${bam_se_sip_1}"
        --csv_sin "${bam_se_sin_0},${bam_se_sin_1}"
        --aln_typ auto
        --err_out "${dir_err}"
        --max_job 1
)

# shellcheck disable=2034
{
    arr_cmd_bam_pe=(
        "${TEST_BASH}" "${scr_exe}"
            --threads 1
            --mode spike
            --csv_mip "${bam_pe_mip_0},${bam_pe_mip_1}"
            --csv_min "${bam_pe_min_0},${bam_pe_min_1}"
            --csv_sip "${bam_pe_sip_0},${bam_pe_sip_1}"
            --csv_sin "${bam_pe_sin_0},${bam_pe_sin_1}"
            --aln_typ auto
            --err_out "${dir_err}"
            --max_job 1
    )

    arr_cmd_cram_se=(
        "${TEST_BASH}" "${scr_exe}"
            --threads 1
            --mode spike
            --csv_mip "${cram_se_mip_0},${cram_se_mip_1}"
            --csv_min "${cram_se_min_0},${cram_se_min_1}"
            --csv_sip "${cram_se_sip_0},${cram_se_sip_1}"
            --csv_sin "${cram_se_sin_0},${cram_se_sin_1}"
            --aln_typ auto
            --ref_fa "${ref_fa}"
            --err_out "${dir_err}"
            --max_job 1
    )

    arr_cmd_cram_pe=(
        "${TEST_BASH}" "${scr_exe}"
            --threads 1
            --mode spike
            --csv_mip "${cram_pe_mip_0},${cram_pe_mip_1}"
            --csv_min "${cram_pe_min_0},${cram_pe_min_1}"
            --csv_sip "${cram_pe_sip_0},${cram_pe_sip_1}"
            --csv_sin "${cram_pe_sin_0},${cram_pe_sin_1}"
            --aln_typ auto
            --ref_fa "${ref_fa}"
            --err_out "${dir_err}"
            --max_job 1
    )

    arr_cmd_mix_lyt_bam=(
        "${TEST_BASH}" "${scr_exe}"
            --threads 1
            --mode spike
            --csv_mip "${bam_se_mip_0},${bam_pe_mip_1}"
            --csv_min "${bam_se_min_0},${bam_pe_min_1}"
            --csv_sip "${bam_se_sip_0},${bam_pe_sip_1}"
            --csv_sin "${bam_se_sin_0},${bam_pe_sin_1}"
            --aln_typ auto
            --err_out "${dir_err}"
            --max_job 1
    )

    arr_cmd_mix_fmt_bam_se_cram_pe=(
        "${TEST_BASH}" "${scr_exe}"
            --threads 1
            --mode spike
            --csv_mip "${bam_se_mip_0},${cram_pe_mip_1}"
            --csv_min "${bam_se_min_0},${cram_pe_min_1}"
            --csv_sip "${bam_se_sip_0},${cram_pe_sip_1}"
            --csv_sin "${bam_se_sin_0},${cram_pe_sin_1}"
            --aln_typ auto
            --ref_fa "${ref_fa}"
            --err_out "${dir_err}"
            --max_job 1
    )
}

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


#  Run one serial spike-in calculation case and assert its assembled table
function run_case_spike() {
    local cas="${1:-}"
    local arr_cmd_nam="${2:-}"
    local row_0="${3:-}"
    local row_1="${4:-}"
    local method="${5:-}"
    local tail_0="${6:-}"
    local tail_1="${7:-}"
    shift 7 || true

    local -n arr_cmd_ref="${arr_cmd_nam}"

    local expt_hdr=true
    local fil_out="${dir_out}/scaling.${cas}.spike.tsv"
    local prt_0="${fil_out}.part.000000"
    local prt_1="${fil_out}.part.000001"
    local nam_job="test_execute_calculate_scaling_factor_spike_${cas}"
    local log="${dir_log}/execute_spike_${cas}.log"
    local -a arr_case=(
        "${arr_cmd_ref[@]}"
        --fil_out "${fil_out}"
        --nam_job "${nam_job}"
    )

    if [[ -n "${method}" ]]; then
        arr_case+=( --method "${method}" )
    fi

    for arg in "$@"; do
        if [[ "${arg}" =~ ^--no[_-]header$ ]]; then
            expt_hdr=false
        fi
    done

    arr_case+=( "$@" )

    if \
        run_capture \
            "execute calculate-scaling-factor spike ${cas}" \
            "${log}" \
            "${arr_case[@]}"
    then
        record_pass \
            "execute_calculate_scaling_factor.sh spike ${cas} exits 0"
    else
        record_fail \
            "execute_calculate_scaling_factor.sh spike ${cas} failed; see" \
            "$(print_relpath "${log}")"
    fi

    assert_file_nonempty \
        "${fil_out}" \
        "execute scaling-factor spike ${cas} final TSV"

    assert_file_nonempty \
        "${prt_0}" \
        "execute scaling-factor spike ${cas} first retained part"

    assert_file_nonempty \
        "${prt_1}" \
        "execute scaling-factor spike ${cas} second retained part"

    if [[ "${expt_hdr}" == "true" ]]; then
        assert_pattern_found \
            "${fil_out}" \
            $'^main_ip\tspike_ip\tmain_in\tspike_in\tspike\tcoef\tnum_mp\tnum_sp\tnum_mn\tnum_sn$' \
            "execute scaling-factor spike ${cas} final TSV has core header"
    else
        assert_pattern_absent \
            "${fil_out}" \
            $'^main_ip\tspike_ip\tmain_in\tspike_in\tspike\tcoef\tnum_mp\tnum_sp\tnum_mn\tnum_sn$' \
            "execute scaling-factor spike ${cas} final TSV omits core header"
    fi

    assert_pattern_found \
        "${fil_out}" \
        "^${row_0}"$'\t'"${tail_0}"'$' \
        "execute scaling-factor spike ${cas} final TSV has first row"

    assert_pattern_found \
        "${fil_out}" \
        "^${row_1}"$'\t'"${tail_1}"'$' \
        "execute scaling-factor spike ${cas} final TSV has second row"
}


#  Dry-run execution should report assembly without writing a final table
fil_out_dry="${dir_out}/scaling.dry.spike.tsv"
log="${dir_log}/execute_spike_dry_run.log"
if \
    run_capture \
        "execute calculate-scaling-factor spike dry-run" \
        "${log}" \
        "${arr_cmd_bam_se[@]}" \
        --fil_out "${fil_out_dry}" \
        --nam_job test_execute_calculate_scaling_factor_spike_dry \
        --dry_run
then
    record_pass "execute_calculate_scaling_factor.sh --dry_run exits 0"
else
    record_fail \
        "execute_calculate_scaling_factor.sh --dry_run failed; see" \
        "$(print_relpath "${log}")"
fi

if [[ ! -e "${fil_out_dry}" ]]; then
    record_pass "execute_calculate_scaling_factor.sh --dry_run writes no TSV"
else
    record_fail "execute_calculate_scaling_factor.sh --dry_run wrote TSV"
fi

assert_pattern_found \
    "${log}" \
    'combine_parts_scaling_factor.sh' \
    "execute_calculate_scaling_factor.sh --dry_run reports combiner command"

assert_pattern_found \
    "${log}" \
    'write_header.sh' \
    "execute_calculate_scaling_factor.sh --dry_run reports header command"

assert_pattern_found \
    "${log}" \
    "${TEST_BASH} ${ROOT_REPO}/scripts/submit_calculate_scaling_factor.sh" \
    "execute_calculate_scaling_factor.sh --dry_run reports Bash-prefixed" \
    "submit command"

assert_pattern_found \
    "${log}" \
    '--method chiprx_alpha_ratio' \
    "execute_calculate_scaling_factor.sh --dry_run resolves default method"

assert_pattern_found \
    "${log}" \
    '--idx_out 0' \
    "execute_calculate_scaling_factor.sh --dry_run reports first part index"

assert_pattern_found \
    "${log}" \
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

fil_out="${dir_out}/scaling.default.spike.tsv"
prt_0="${fil_out}.part.000000"
prt_1="${fil_out}.part.000001"
nam_job="test_execute_calculate_scaling_factor_spike_default"
log_out_0="${dir_err}/${nam_job}.IP_WT_log_Brn1_rep1.sc.stdout.txt"
log_err_0="${dir_err}/${nam_job}.IP_WT_log_Brn1_rep1.sc.stderr.txt"
log_out_1="${dir_err}/${nam_job}.IP_WT_log_Brn1_rep2.sc.stdout.txt"
log_err_1="${dir_err}/${nam_job}.IP_WT_log_Brn1_rep2.sc.stderr.txt"

assert_file_exists \
    "${log_out_0}" \
    "execute scaling-factor spike default first submit stdout log exists"

assert_file_exists \
    "${log_err_0}" \
    "execute scaling-factor spike default first submit stderr log exists"

assert_file_exists \
    "${log_out_1}" \
    "execute scaling-factor spike default second submit stdout log exists"

assert_file_exists \
    "${log_err_1}" \
    "execute scaling-factor spike default second submit stderr log exists"

assert_pattern_found \
    "${log_err_0}" \
    'idx_out=0' \
    "execute scaling-factor spike default first submit uses idx_out=0"

assert_pattern_found \
    "${log_err_1}" \
    'idx_out=1' \
    "execute scaling-factor spike default second submit uses idx_out=1"

assert_pattern_found \
    "${log_err_0}" \
    'typ_mp=se' \
    "execute scaling-factor spike default auto-detects first main IP as SE"

assert_pattern_found \
    "${log_err_0}" \
    'num_mp=3' \
    "execute scaling-factor spike default counts first main IP alignments"

assert_pattern_found \
    "${log_err_0}" \
    'num_sp=1' \
    "execute scaling-factor spike default counts first spike IP alignments"

assert_pattern_found \
    "${log_err_0}" \
    'num_mn=2' \
    "execute scaling-factor spike default counts first main input alignments"

assert_pattern_found \
    "${log_err_0}" \
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

log_err_0="${dir_err}/test_execute_calculate_scaling_factor_spike_pe_bam.IP_WT_G1_Hho1_6336.sc.stderr.txt"
log_err_1="${dir_err}/test_execute_calculate_scaling_factor_spike_pe_bam.IP_WT_G1_Hho1_6337.sc.stderr.txt"

assert_pattern_found \
    "${log_err_0}" \
    'typ_mp=pe' \
    "execute scaling-factor spike PE BAM auto-detects first main IP as PE"

assert_pattern_found \
    "${log_err_1}" \
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

log_err_0="${dir_err}/test_execute_calculate_scaling_factor_spike_se_cram.IP_WT_log_Brn1_rep1.sc.stderr.txt"

assert_pattern_found \
    "${log_err_0}" \
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

log_err_0="${dir_err}/test_execute_calculate_scaling_factor_spike_pe_cram.IP_WT_G1_Hho1_6336.sc.stderr.txt"
log_err_1="${dir_err}/test_execute_calculate_scaling_factor_spike_pe_cram.IP_WT_G1_Hho1_6337.sc.stderr.txt"

assert_pattern_found \
    "${log_err_0}" \
    'typ_mp=pe' \
    "execute scaling-factor spike PE CRAM auto-detects first main IP as PE"

assert_pattern_found \
    "${log_err_1}" \
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

log_err_0="${dir_err}/test_execute_calculate_scaling_factor_spike_mixed_layout.IP_WT_log_Brn1_rep1.sc.stderr.txt"
log_err_1="${dir_err}/test_execute_calculate_scaling_factor_spike_mixed_layout.IP_WT_G1_Hho1_6337.sc.stderr.txt"

assert_pattern_found \
    "${log_err_0}" \
    'typ_mp=se' \
    "execute scaling-factor spike mixed layout keeps first main IP as SE"

assert_pattern_found \
    "${log_err_1}" \
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

log_err_1="${dir_err}/test_execute_calculate_scaling_factor_spike_mixed_format.IP_WT_G1_Hho1_6337.sc.stderr.txt"

assert_pattern_found \
    "${log_err_1}" \
    'typ_mp=pe' \
    "execute scaling-factor spike mixed BAM/CRAM list reads PE CRAM input"


#  CRAM input without a reference FASTA should fail before worker execution
log="${dir_log}/execute_spike_cram_missing_ref.log"
if \
    run_capture \
        "execute calculate-scaling-factor spike CRAM missing reference" \
        "${log}" \
        "${arr_cmd_bam_se[@]}" \
        --csv_mip "${cram_se_mip_0},${cram_se_mip_1}" \
        --csv_min "${cram_se_min_0},${cram_se_min_1}" \
        --csv_sip "${cram_se_sip_0},${cram_se_sip_1}" \
        --csv_sin "${cram_se_sin_0},${cram_se_sin_1}" \
        --fil_out "${dir_out}/scaling.cram_missing_ref.spike.tsv" \
        --nam_job test_execute_calculate_scaling_factor_spike_cram_missing_ref
then
    record_fail \
        "execute_calculate_scaling_factor.sh CRAM without ref unexpectedly" \
        "pass"
else
    assert_pattern_found \
        "${log}" \
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
    --dep_mip 10 \
    --dep_min 20 \
    --dep_sip 2 \
    --dep_sin 4


#  Unknown spike-in method should fail before running workers
log="${dir_log}/execute_spike_invalid_method.log"
if \
    run_capture \
        "execute calculate-scaling-factor spike invalid method" \
        "${log}" \
        "${arr_cmd_bam_se[@]}" \
        --method not_a_method \
        --fil_out "${dir_out}/scaling.invalid_method.spike.tsv" \
        --nam_job test_execute_calculate_scaling_factor_spike_invalid_method
then
    record_fail \
        "execute_calculate_scaling_factor.sh invalid method unexpectedly pass"
else
    assert_pattern_found \
        "${log}" \
        'is not recognized' \
        "execute_calculate_scaling_factor.sh rejects invalid spike method"
fi


#  Existing final output should fail unless '--force' is supplied
log="${dir_log}/execute_spike_existing.log"
if \
    run_capture \
        "execute calculate-scaling-factor spike existing output" \
        "${log}" \
        "${arr_cmd_bam_se[@]}" \
        --fil_out "${fil_out}" \
        --nam_job "${nam_job}"
then
    record_fail \
        "execute_calculate_scaling_factor.sh existing output unexpectedly pass"
else
    assert_pattern_found \
        "${log}" \
        'output file already exists' \
        "execute_calculate_scaling_factor.sh rejects existing output"
fi


#  '--force' should replace output and '--no_parts' should remove worker part
#+ files
log="${dir_log}/execute_spike_force_no_parts.log"
if \
    run_capture \
        "execute calculate-scaling-factor spike force no-parts" \
        "${log}" \
        "${arr_cmd_bam_se[@]}" \
        --fil_out "${fil_out}" \
        --nam_job "${nam_job}" \
        --force \
        --no_parts
then
    record_pass \
        "execute_calculate_scaling_factor.sh --force --no_parts exits 0"
else
    record_fail \
        "execute_calculate_scaling_factor.sh --force --no_parts failed; see" \
        "$(print_relpath "${log}")"
fi

assert_file_nonempty \
    "${fil_out}" \
    "execute scaling-factor spike replaced final TSV"

if [[ ! -e "${prt_0}" && ! -e "${prt_1}" ]]; then
    record_pass "execute scaling-factor spike --no_parts removes worker parts"
else
    record_fail "execute scaling-factor spike --no_parts retained worker parts"
fi

finish
