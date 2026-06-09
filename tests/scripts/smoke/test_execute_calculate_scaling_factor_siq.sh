#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_calculate_scaling_factor_siq.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="execute calculate-scaling-factor siQ"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

print_section "${TEST_NAME}"

#  Define fixture, output, and worker-input paths
scr_exe="${ROOT_REPO}/scripts/execute_calculate_scaling_factor.sh"
cfg_met="${ROOT_REPO}/data/raw/docs/parse_metadata_siqchip.yml"
dir_fix="${ROOT_REPO}/tests/calculate_scaling_factor/fixtures"
dir_bam="${dir_fix}/bam"
dir_bam_se="${dir_bam}/se"
dir_bam_pe="${dir_bam}/pe"
dir_met="${dir_fix}/metadata"
tbl_met="${dir_met}/measurements_siqchip.tsv"

tmp="${TEST_DIR_TMP}/execute_calculate_scaling_factor_siq"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/calculate_scaling_factor"

bam_pe_mip_0="${dir_bam_pe}/IP_WT_G1_Hho1_6336.sc.bam"
bam_pe_mip_1="${dir_bam_pe}/IP_WT_G1_Hho1_6337.sc.bam"
bam_pe_min_0="${dir_bam_pe}/in_WT_G1_Hho1_6336.sc.bam"
bam_pe_min_1="${dir_bam_pe}/in_WT_G1_Hho1_6337.sc.bam"

bam_se_mip_0="${dir_bam_se}/IP_WT_log_Brn1_rep1.sc.bam"
bam_se_mip_1="${dir_bam_se}/IP_WT_log_Brn1_rep2.sc.bam"
bam_se_min_0="${dir_bam_se}/in_WT_log_Brn1_rep1.sc.bam"
bam_se_min_1="${dir_bam_se}/in_WT_log_Brn1_rep2.sc.bam"

rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

if ! \
    require_files_nonempty \
        "${scr_exe}" \
        "${cfg_met}" \
        "${tbl_met}" \
        "${bam_pe_mip_0}" \
        "${bam_pe_mip_0}.bai" \
        "${bam_pe_mip_1}" \
        "${bam_pe_mip_1}.bai" \
        "${bam_pe_min_0}" \
        "${bam_pe_min_0}.bai" \
        "${bam_pe_min_1}" \
        "${bam_pe_min_1}.bai" \
        "${bam_se_mip_0}" \
        "${bam_se_mip_0}.bai" \
        "${bam_se_mip_1}" \
        "${bam_se_mip_1}.bai" \
        "${bam_se_min_0}" \
        "${bam_se_min_0}.bai" \
        "${bam_se_min_1}" \
        "${bam_se_min_1}.bai"
then
    finish
    exit $?
fi

# shellcheck disable=2054
arr_cmd_bam_pe=(
    "${TEST_BASH}" "${scr_exe}"
        --threads 1
        --mode siq
        --csv_mip "${bam_pe_mip_0},${bam_pe_mip_1}"
        --csv_min "${bam_pe_min_0},${bam_pe_min_1}"
        --aln_typ auto
        --tbl_met "${tbl_met}"
        --cfg_met "${cfg_met}"
        --eqn 6nd
        --err_out "${dir_err}"
        --max_job 1
)

arr_cmd_bam_se=(
    "${TEST_BASH}" "${scr_exe}"
        --threads 1
        --mode siq
        --csv_mip "${bam_se_mip_0},${bam_se_mip_1}"
        --csv_min "${bam_se_min_0},${bam_se_min_1}"
        --aln_typ auto
        --len_def 150
        --tbl_met "${tbl_met}"
        --cfg_met "${cfg_met}"
        --eqn 6nd
        --err_out "${dir_err}"
        --max_job 1
)

row_bam_pe_0="${bam_pe_mip_0}"$'\t'"${bam_pe_min_0}"
row_bam_pe_1="${bam_pe_mip_1}"$'\t'"${bam_pe_min_1}"
row_bam_se_0="${bam_se_mip_0}"$'\t'"${bam_se_min_0}"
row_bam_se_1="${bam_se_mip_1}"$'\t'"${bam_se_min_1}"


#  Run one serial siQ calculation case and assert its assembled table
function run_case_siq() {
    local cas="${1:-}"
    local arr_cmd_nam="${2:-}"
    local row_0="${3:-}"
    local row_1="${4:-}"
    local tail_0="${5:-}"
    local tail_1="${6:-}"
    shift 6 || true

    local -n arr_cmd_ref="${arr_cmd_nam}"

    local expt_hdr=true
    local fil_out="${dir_out}/scaling.${cas}.siq.tsv"
    local prt_0="${fil_out}.part.000000"
    local prt_1="${fil_out}.part.000001"
    local nam_job="test_execute_calculate_scaling_factor_siq_${cas}"
    local log="${dir_log}/execute_siq_${cas}.log"
    local -a arr_case=(
        "${arr_cmd_ref[@]}"
        --fil_out "${fil_out}"
        --nam_job "${nam_job}"
    )

    for arg in "$@"; do
        if [[ "${arg}" =~ ^--no[_-]header$ ]]; then
            expt_hdr=false
        fi
    done

    arr_case+=( "$@" )

    if \
        run_capture \
            "execute calculate-scaling-factor siQ ${cas}" \
            "${log}" \
            "${arr_case[@]}"
    then
        record_pass "execute_calculate_scaling_factor.sh siQ ${cas} exits 0"
    else
        record_fail \
            "execute_calculate_scaling_factor.sh siQ ${cas} failed; see" \
            "$(print_relpath "${log}")"
    fi

    assert_file_nonempty \
        "${fil_out}" \
        "execute scaling-factor siQ ${cas} final TSV"

    assert_file_nonempty \
        "${prt_0}" \
        "execute scaling-factor siQ ${cas} first retained part"

    assert_file_nonempty \
        "${prt_1}" \
        "execute scaling-factor siQ ${cas} second retained part"

    if [[ "${expt_hdr}" == "true" ]]; then
        assert_pattern_found \
            "${fil_out}" \
            $'^fil_ip\tfil_in\tsiq\teqn\tmass_ip\tmass_in\tvol_all\tvol_in\tdep_ip\tdep_in\tlen_ip\tlen_in$' \
            "execute scaling-factor siQ ${cas} final TSV has core header"
    else
        assert_pattern_absent \
            "${fil_out}" \
            $'^fil_ip\tfil_in\tsiq\teqn\tmass_ip\tmass_in\tvol_all\tvol_in\tdep_ip\tdep_in\tlen_ip\tlen_in$' \
            "execute scaling-factor siQ ${cas} final TSV omits core header"
    fi

    assert_pattern_found \
        "${fil_out}" \
        "^${row_0}"$'\t'"${tail_0}"'$' \
        "execute scaling-factor siQ ${cas} final TSV has first row"

    assert_pattern_found \
        "${fil_out}" \
        "^${row_1}"$'\t'"${tail_1}"'$' \
        "execute scaling-factor siQ ${cas} final TSV has second row"
}


#  Dry-run execution should report worker, combiner, and header commands
fil_out_dry="${dir_out}/scaling.dry.siq.tsv"
log="${dir_log}/execute_siq_dry_run.log"
if \
    run_capture \
        "execute calculate-scaling-factor siQ dry-run" \
        "${log}" \
        "${arr_cmd_bam_pe[@]}" \
        --fil_out "${fil_out_dry}" \
        --nam_job test_execute_calculate_scaling_factor_siq_dry \
        --dry_run
then
    record_pass "execute_calculate_scaling_factor.sh siQ --dry_run exits 0"
else
    record_fail \
        "execute_calculate_scaling_factor.sh siQ --dry_run failed; see" \
        "$(print_relpath "${log}")"
fi

if [[ ! -e "${fil_out_dry}" ]]; then
    record_pass "execute_calculate_scaling_factor.sh siQ --dry_run writes no TSV"
else
    record_fail "execute_calculate_scaling_factor.sh siQ --dry_run wrote TSV"
fi

assert_pattern_found \
    "${log}" \
    'combine_parts_scaling_factor.sh' \
    "execute_calculate_scaling_factor.sh siQ --dry_run reports combiner"

assert_pattern_found \
    "${log}" \
    'write_header.sh' \
    "execute_calculate_scaling_factor.sh siQ --dry_run reports header"

assert_pattern_found \
    "${log}" \
    "${TEST_BASH} ${ROOT_REPO}/scripts/submit_calculate_scaling_factor.sh" \
    "execute_calculate_scaling_factor.sh siQ --dry_run reports submit"

assert_pattern_found \
    "${log}" \
    '--idx_out 0' \
    "execute_calculate_scaling_factor.sh siQ --dry_run reports first index"

assert_pattern_found \
    "${log}" \
    '--idx_out 1' \
    "execute_calculate_scaling_factor.sh siQ --dry_run reports second index"


#  PE BAM input should compute lengths from paired-end fragment sizes
run_case_siq \
    pe_bam \
    arr_cmd_bam_pe \
    "${row_bam_pe_0}" \
    "${row_bam_pe_1}" \
    $'0.002660098522167487697376\t6nd\t2.7\t72.5\t300\t20\t3\t2\t20\t20' \
    $'0.00440373436674299824356\t6nd\t5\t81.1\t300\t20\t2\t3\t20\t20'

fil_out="${dir_out}/scaling.pe_bam.siq.tsv"
prt_0="${fil_out}.part.000000"
prt_1="${fil_out}.part.000001"
nam_job="test_execute_calculate_scaling_factor_siq_pe_bam"
log_err_0="${dir_err}/${nam_job}.IP_WT_G1_Hho1_6336.sc.stderr.txt"
log_err_1="${dir_err}/${nam_job}.IP_WT_G1_Hho1_6337.sc.stderr.txt"

assert_pattern_found \
    "${log_err_0}" \
    'typ_ip=pe' \
    "execute scaling-factor siQ PE BAM auto-detects first IP as PE"

assert_pattern_found \
    "${log_err_1}" \
    'typ_ip=pe' \
    "execute scaling-factor siQ PE BAM auto-detects second IP as PE"


#  SE BAM input should use explicit default fragment lengths
run_case_siq \
    se_bam \
    arr_cmd_bam_se \
    "${row_bam_se_0}" \
    "${row_bam_se_1}" \
    $'0.004761904761904762334312\t6nd\t4\t60\t300\t20\t3\t2\t150\t150' \
    $'0.004761904761904761466951\t6nd\t5\t75\t300\t20\t2\t3\t150\t150'

nam_job="test_execute_calculate_scaling_factor_siq_se_bam"
log_err_0="${dir_err}/${nam_job}.IP_WT_log_Brn1_rep1.sc.stderr.txt"
log_err_1="${dir_err}/${nam_job}.IP_WT_log_Brn1_rep2.sc.stderr.txt"

assert_pattern_found \
    "${log_err_0}" \
    'typ_ip=se' \
    "execute scaling-factor siQ SE BAM auto-detects first IP as SE"

assert_pattern_found \
    "${log_err_1}" \
    'typ_ip=se' \
    "execute scaling-factor siQ SE BAM auto-detects second IP as SE"

run_case_siq \
    no_header \
    arr_cmd_bam_pe \
    "${row_bam_pe_0}" \
    "${row_bam_pe_1}" \
    $'0.002660098522167487697376\t6nd\t2.7\t72.5\t300\t20\t3\t2\t20\t20' \
    $'0.00440373436674299824356\t6nd\t5\t81.1\t300\t20\t2\t3\t20\t20' \
    --no_header


#  Existing final output should fail unless '--force' is supplied
fil_out="${dir_out}/scaling.pe_bam.siq.tsv"
prt_0="${fil_out}.part.000000"
prt_1="${fil_out}.part.000001"
nam_job="test_execute_calculate_scaling_factor_siq_pe_bam"
log="${dir_log}/execute_siq_existing.log"
if \
    run_capture \
        "execute calculate-scaling-factor siQ existing output" \
        "${log}" \
        "${arr_cmd_bam_pe[@]}" \
        --fil_out "${fil_out}" \
        --nam_job "${nam_job}"
then
    record_fail \
        "execute_calculate_scaling_factor.sh siQ existing output" \
        "unexpectedly passed"
else
    assert_pattern_found \
        "${log}" \
        'output file already exists' \
        "execute_calculate_scaling_factor.sh siQ rejects existing output"
fi


#  '--force' should replace output and '--no_parts' should remove worker parts
log="${dir_log}/execute_siq_force_no_parts.log"
if \
    run_capture \
        "execute calculate-scaling-factor siQ force no-parts" \
        "${log}" \
        "${arr_cmd_bam_pe[@]}" \
        --fil_out "${fil_out}" \
        --nam_job "${nam_job}" \
        --force \
        --no_parts
then
    record_pass \
        "execute_calculate_scaling_factor.sh siQ --force --no_parts exits 0"
else
    record_fail \
        "execute_calculate_scaling_factor.sh siQ --force --no_parts failed;" \
        "see $(print_relpath "${log}")"
fi

assert_file_nonempty \
    "${fil_out}" \
    "execute scaling-factor siQ replaced final TSV"

if [[ ! -e "${prt_0}" && ! -e "${prt_1}" ]]; then
    record_pass "execute scaling-factor siQ --no_parts removes worker parts"
else
    record_fail "execute scaling-factor siQ --no_parts retained worker parts"
fi

finish
