#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_calculate_scaling_factor_siq.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5) were used in development.
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
dir_cram="${dir_fix}/cram"
dir_bam_se="${dir_bam}/se"
dir_bam_pe="${dir_bam}/pe"
dir_cram_pe="${dir_cram}/pe"
dir_met="${dir_fix}/metadata"

tbl_met="${dir_met}/measurements_siqchip.tsv"
tbl_gz="${dir_met}/measurements_siqchip.tsv.gz"
tbl_als="${dir_met}/measurements_siqchip_aliases.tsv"
tbl_pfx="${dir_met}/measurements_siqchip_skip_prefixes.tsv"
tbl_col="${dir_met}/measurements_siqchip_alias_collision.tsv"
tbl_dup="${dir_met}/measurements_siqchip_duplicate_match.tsv"
tbl_trt="${dir_met}/measurements_siqchip_treatment.tsv"
cfg_rgx="${dir_fix}/config/parse_metadata_siqchip_regex.yml"
cfg_trt="${dir_fix}/config/parse_metadata_siqchip_match_treatment.yml"
ref_fa="${dir_fix}/reference/tiny.fa"

tmp="${TEST_DIR_TMP}/execute_calculate_scaling_factor_siq"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/calculate_scaling_factor"

bam_pe_mip_0="${dir_bam_pe}/IP_WT_G1_Hho1_6336.sc.bam"
bam_pe_mip_1="${dir_bam_pe}/IP_WT_G1_Hho1_6337.sc.bam"
bam_pe_min_0="${dir_bam_pe}/in_WT_G1_Hho1_6336.sc.bam"
bam_pe_min_1="${dir_bam_pe}/in_WT_G1_Hho1_6337.sc.bam"
hu_pe_mip_0="${dir_bam_pe}/IP_WT_G1_HU_Hho1_6336.sc.bam"
hu_pe_mip_1="${dir_bam_pe}/IP_WT_G1_HU_Hho1_6337.sc.bam"
hu_pe_min_0="${dir_bam_pe}/in_WT_G1_HU_Hho1_6336.sc.bam"
hu_pe_min_1="${dir_bam_pe}/in_WT_G1_HU_Hho1_6337.sc.bam"

cram_pe_mip_0="${dir_cram_pe}/IP_WT_G1_Hho1_6336.sc.cram"
cram_pe_mip_1="${dir_cram_pe}/IP_WT_G1_Hho1_6337.sc.cram"
cram_pe_min_0="${dir_cram_pe}/in_WT_G1_Hho1_6336.sc.cram"
cram_pe_min_1="${dir_cram_pe}/in_WT_G1_Hho1_6337.sc.cram"

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
        "${cfg_rgx}" \
        "${cfg_trt}" \
        "${tbl_met}" \
        "${tbl_gz}" \
        "${tbl_als}" \
        "${tbl_pfx}" \
        "${tbl_col}" \
        "${tbl_dup}" \
        "${tbl_trt}" \
        "${bam_pe_mip_0}" \
        "${bam_pe_mip_0}.bai" \
        "${bam_pe_mip_1}" \
        "${bam_pe_mip_1}.bai" \
        "${bam_pe_min_0}" \
        "${bam_pe_min_0}.bai" \
        "${bam_pe_min_1}" \
        "${bam_pe_min_1}.bai" \
        "${hu_pe_mip_0}" \
        "${hu_pe_mip_0}.bai" \
        "${hu_pe_mip_1}" \
        "${hu_pe_mip_1}.bai" \
        "${hu_pe_min_0}" \
        "${hu_pe_min_0}.bai" \
        "${hu_pe_min_1}" \
        "${hu_pe_min_1}.bai" \
        "${cram_pe_mip_0}" \
        "${cram_pe_mip_1}" \
        "${cram_pe_min_0}" \
        "${cram_pe_min_1}" \
        "${ref_fa}" \
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

# shellcheck disable=SC2034
{
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

    arr_cmd_cram_pe=(
        "${TEST_BASH}" "${scr_exe}"
            --threads 1
            --mode siq
            --csv_mip "${cram_pe_mip_0},${cram_pe_mip_1}"
            --csv_min "${cram_pe_min_0},${cram_pe_min_1}"
            --aln_typ auto
            --ref_fa "${ref_fa}"
            --tbl_met "${tbl_met}"
            --cfg_met "${cfg_met}"
            --eqn 6nd
            --err_out "${dir_err}"
            --max_job 1
    )

    arr_cmd_mxd_lyt=(
        "${TEST_BASH}" "${scr_exe}"
            --threads 1
            --mode siq
            --csv_mip "${bam_se_mip_0},${bam_pe_mip_1}"
            --csv_min "${bam_se_min_0},${bam_pe_min_1}"
            --aln_typ auto
            --len_def 150
            --tbl_met "${tbl_met}"
            --cfg_met "${cfg_met}"
            --eqn 6nd
            --err_out "${dir_err}"
            --max_job 1
    )

    arr_cmd_mxd_fmt=(
        "${TEST_BASH}" "${scr_exe}"
            --threads 1
            --mode siq
            --csv_mip "${bam_pe_mip_0},${cram_pe_mip_1}"
            --csv_min "${bam_pe_min_0},${cram_pe_min_1}"
            --aln_typ auto
            --ref_fa "${ref_fa}"
            --tbl_met "${tbl_met}"
            --cfg_met "${cfg_met}"
            --eqn 6nd
            --err_out "${dir_err}"
            --max_job 1
    )

    arr_cmd_trt=(
        "${TEST_BASH}" "${scr_exe}"
            --threads 1
            --mode siq
            --csv_mip "${hu_pe_mip_0},${hu_pe_mip_1}"
            --csv_min "${hu_pe_min_0},${hu_pe_min_1}"
            --aln_typ auto
            --tbl_met "${tbl_trt}"
            --cfg_met "${cfg_trt}"
            --eqn 6nd
            --err_out "${dir_err}"
            --max_job 1
    )
}

row_bam_pe_0="${bam_pe_mip_0}"$'\t'"${bam_pe_min_0}"
row_bam_pe_1="${bam_pe_mip_1}"$'\t'"${bam_pe_min_1}"
row_hu_pe_0="${hu_pe_mip_0}"$'\t'"${hu_pe_min_0}"
row_hu_pe_1="${hu_pe_mip_1}"$'\t'"${hu_pe_min_1}"
row_cram_pe_0="${cram_pe_mip_0}"$'\t'"${cram_pe_min_0}"
row_cram_pe_1="${cram_pe_mip_1}"$'\t'"${cram_pe_min_1}"
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

    run_case_scaling_factor_execute \
        "${cas}" \
        siQ \
        "${arr_cmd_nam}" \
        "${dir_out}" \
        "${dir_log}" \
        siq \
        test_execute_calculate_scaling_factor_siq \
        $'^fil_ip\tfil_in\tsiq\teqn\tmass_ip\tmass_in\tvol_all\tvol_in\tdep_ip\tdep_in\tlen_ip\tlen_in$' \
        "${row_0}" \
        "${row_1}" \
        "${tail_0}" \
        "${tail_1}" \
        "$@"
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


#  PE CRAM input should work when an explicit reference is supplied
run_case_siq \
    pe_cram \
    arr_cmd_cram_pe \
    "${row_cram_pe_0}" \
    "${row_cram_pe_1}" \
    $'0.002660098522167487697376\t6nd\t2.7\t72.5\t300\t20\t3\t2\t20\t20' \
    $'0.00440373436674299824356\t6nd\t5\t81.1\t300\t20\t2\t3\t20\t20'

nam_job="test_execute_calculate_scaling_factor_siq_pe_cram"
log_err_0="${dir_err}/${nam_job}.IP_WT_G1_Hho1_6336.sc.stderr.txt"
log_err_1="${dir_err}/${nam_job}.IP_WT_G1_Hho1_6337.sc.stderr.txt"

assert_pattern_found \
    "${log_err_0}" \
    'typ_ip=pe' \
    "execute scaling-factor siQ PE CRAM auto-detects first IP as PE"

assert_pattern_found \
    "${log_err_1}" \
    "ref_fa=${ref_fa}" \
    "execute scaling-factor siQ PE CRAM forwards ref_fa"


#  CRAM input should fail clearly when no reference is supplied
log="${dir_log}/execute_siq_cram_missing_ref.log"
if \
    run_capture \
        "execute calculate-scaling-factor siQ CRAM missing reference" \
        "${log}" \
        "${TEST_BASH}" "${scr_exe}" \
            --threads 1 \
            --mode siq \
            --csv_mip "${cram_pe_mip_0},${cram_pe_mip_1}" \
            --csv_min "${cram_pe_min_0},${cram_pe_min_1}" \
            --aln_typ auto \
            --tbl_met "${tbl_met}" \
            --cfg_met "${cfg_met}" \
            --eqn 6nd \
            --fil_out "${dir_out}/scaling.cram_missing_ref.siq.tsv" \
            --err_out "${dir_err}" \
            --nam_job test_execute_calculate_scaling_factor_siq_cram_missing_ref \
            --max_job 1
then
    record_fail \
        "execute_calculate_scaling_factor.sh siQ CRAM without reference" \
        "unexpectedly passed"
else
    assert_pattern_found \
        "${log}" \
        "'--ref_fa' is required" \
        "execute_calculate_scaling_factor.sh rejects CRAM without ref_fa"
fi


#  Mixed SE/PE BAM lists should keep per-file auto-detection
run_case_siq \
    mixed_layout \
    arr_cmd_mxd_lyt \
    "${row_bam_se_0}" \
    "${row_bam_pe_1}" \
    $'0.004761904761904762334312\t6nd\t4\t60\t300\t20\t3\t2\t150\t150' \
    $'0.00440373436674299824356\t6nd\t5\t81.1\t300\t20\t2\t3\t20\t20'

nam_job="test_execute_calculate_scaling_factor_siq_mixed_layout"
log_err_0="${dir_err}/${nam_job}.IP_WT_log_Brn1_rep1.sc.stderr.txt"
log_err_1="${dir_err}/${nam_job}.IP_WT_G1_Hho1_6337.sc.stderr.txt"

assert_pattern_found \
    "${log_err_0}" \
    'typ_ip=se' \
    "execute scaling-factor siQ mixed layout keeps first IP as SE"

assert_pattern_found \
    "${log_err_1}" \
    'typ_ip=pe' \
    "execute scaling-factor siQ mixed layout keeps second IP as PE"


#  Mixed BAM/CRAM lists should route CRAM samples with the shared reference
run_case_siq \
    mixed_format \
    arr_cmd_mxd_fmt \
    "${row_bam_pe_0}" \
    "${row_cram_pe_1}" \
    $'0.002660098522167487697376\t6nd\t2.7\t72.5\t300\t20\t3\t2\t20\t20' \
    $'0.00440373436674299824356\t6nd\t5\t81.1\t300\t20\t2\t3\t20\t20'

nam_job="test_execute_calculate_scaling_factor_siq_mixed_format"
log_err_1="${dir_err}/${nam_job}.IP_WT_G1_Hho1_6337.sc.stderr.txt"

assert_pattern_found \
    "${log_err_1}" \
    "ref_fa=${ref_fa}" \
    "execute scaling-factor siQ mixed BAM/CRAM list forwards ref_fa"


#  Metadata parser variants should propagate through the execute wrapper
run_case_siq \
    gzip_metadata \
    arr_cmd_bam_pe \
    "${row_bam_pe_0}" \
    "${row_bam_pe_1}" \
    $'0.002660098522167487697376\t6nd\t2.7\t72.5\t300\t20\t3\t2\t20\t20' \
    $'0.00440373436674299824356\t6nd\t5\t81.1\t300\t20\t2\t3\t20\t20' \
    --tbl_met "${tbl_gz}"

run_case_siq \
    alias_metadata \
    arr_cmd_bam_pe \
    "${row_bam_pe_0}" \
    "${row_bam_pe_1}" \
    $'0.002660098522167487697376\t6nd\t2.7\t72.5\t300\t20\t3\t2\t20\t20' \
    $'0.00440373436674299824356\t6nd\t5\t81.1\t300\t20\t2\t3\t20\t20' \
    --tbl_met "${tbl_als}"

run_case_siq \
    skip_prefixes \
    arr_cmd_bam_pe \
    "${row_bam_pe_0}" \
    "${row_bam_pe_1}" \
    $'0.002660098522167487697376\t6nd\t2.7\t72.5\t300\t20\t3\t2\t20\t20' \
    $'0.00440373436674299824356\t6nd\t5\t81.1\t300\t20\t2\t3\t20\t20' \
    --tbl_met "${tbl_pfx}"

run_case_siq \
    regex_config \
    arr_cmd_bam_pe \
    "${row_bam_pe_0}" \
    "${row_bam_pe_1}" \
    $'0.002660098522167487697376\t6nd\t2.7\t72.5\t300\t20\t3\t2\t20\t20' \
    $'0.00440373436674299824356\t6nd\t5\t81.1\t300\t20\t2\t3\t20\t20' \
    --cfg_met "${cfg_rgx}"

run_case_siq \
    treatment \
    arr_cmd_trt \
    "${row_hu_pe_0}" \
    "${row_hu_pe_1}" \
    $'0.007963320463320463019063\t6nd\t9.9\t88.8\t300\t20\t3\t2\t20\t20' \
    $'0.007936507936507936067372\t6nd\t11.1\t99.9\t300\t20\t2\t3\t20\t20'


#  Metadata parser failures should propagate through execute-layer logs
log="${dir_log}/execute_siq_alias_collision.log"
fil_err="${dir_err}/test_execute_calculate_scaling_factor_siq_alias_collision.IP_WT_G1_Hho1_6336.sc.stderr.txt"
if \
    run_capture \
        "execute calculate-scaling-factor siQ alias collision" \
        "${log}" \
        "${TEST_BASH}" "${scr_exe}" \
            --threads 1 \
            --mode siq \
            --csv_mip "${bam_pe_mip_0}" \
            --csv_min "${bam_pe_min_0}" \
            --aln_typ auto \
            --tbl_met "${tbl_col}" \
            --cfg_met "${cfg_met}" \
            --eqn 6nd \
            --fil_out "${dir_out}/scaling.alias_collision.siq.tsv" \
            --err_out "${dir_err}" \
            --nam_job test_execute_calculate_scaling_factor_siq_alias_collision \
            --max_job 1
then
    record_fail \
        "execute_calculate_scaling_factor.sh siQ alias collision" \
        "unexpectedly passed"
else
    assert_pattern_found \
        "${fil_err}" \
        "multiple columns that normalize to the same canonical name" \
        "execute_calculate_scaling_factor.sh siQ alias collision fails clearly"
fi

log="${dir_log}/execute_siq_duplicate_match.log"
fil_err="${dir_err}/test_execute_calculate_scaling_factor_siq_duplicate_match.IP_WT_G1_Hho1_6336.sc.stderr.txt"
if \
    run_capture \
        "execute calculate-scaling-factor siQ duplicate match" \
        "${log}" \
        "${TEST_BASH}" "${scr_exe}" \
            --threads 1 \
            --mode siq \
            --csv_mip "${bam_pe_mip_0}" \
            --csv_min "${bam_pe_min_0}" \
            --aln_typ auto \
            --tbl_met "${tbl_dup}" \
            --cfg_met "${cfg_met}" \
            --eqn 6nd \
            --fil_out "${dir_out}/scaling.duplicate_match.siq.tsv" \
            --err_out "${dir_err}" \
            --nam_job test_execute_calculate_scaling_factor_siq_duplicate_match \
            --max_job 1
then
    record_fail \
        "execute_calculate_scaling_factor.sh siQ duplicate match" \
        "unexpectedly passed"
else
    assert_pattern_found \
        "${fil_err}" \
        "Multiple matching rows found" \
        "execute_calculate_scaling_factor.sh siQ duplicate match fails clearly"
fi


#  Equation selection should propagate through the execute wrapper
run_case_siq \
    eqn5 \
    arr_cmd_bam_pe \
    "${row_bam_pe_0}" \
    "${row_bam_pe_1}" \
    $'0.001655172413793103407958\t5\t2.7\t72.5\t300\t20\t3\t2\t20\t20' \
    $'0.006165228113440198234874\t5\t5\t81.1\t300\t20\t2\t3\t20\t20' \
    --eqn 5

run_case_siq \
    eqn5nd \
    arr_cmd_bam_pe \
    "${row_bam_pe_0}" \
    "${row_bam_pe_1}" \
    $'0.002482758620689655328778\t5nd\t2.7\t72.5\t300\t20\t3\t2\t20\t20' \
    $'0.004110152075626798823249\t5nd\t5\t81.1\t300\t20\t2\t3\t20\t20' \
    --eqn 5nd

run_case_siq \
    eqn6 \
    arr_cmd_bam_pe \
    "${row_bam_pe_0}" \
    "${row_bam_pe_1}" \
    $'0.001773399014778325203864\t6\t2.7\t72.5\t300\t20\t3\t2\t20\t20' \
    $'0.00660560155011449736534\t6\t5\t81.1\t300\t20\t2\t3\t20\t20' \
    --eqn 6


#  Broadcast length and depth overrides should reach each submit worker
run_case_siq \
    len_override \
    arr_cmd_bam_pe \
    "${row_bam_pe_0}" \
    "${row_bam_pe_1}" \
    $'0.005320197044334976262114\t6nd\t2.7\t72.5\t300\t20\t3\t2\t100\t200' \
    $'0.00880746873348599648712\t6nd\t5\t81.1\t300\t20\t2\t3\t100\t200' \
    --len_mip 100 \
    --len_min 200

run_case_siq \
    dep_override \
    arr_cmd_bam_pe \
    "${row_bam_pe_0}" \
    "${row_bam_pe_1}" \
    $'0.005320197044334976262114\t6\t2.7\t72.5\t300\t20\t10\t20\t20\t20' \
    $'0.00880746873348599648712\t6\t5\t81.1\t300\t20\t10\t20\t20\t20' \
    --eqn 6 \
    --dep_mip 10 \
    --dep_min 20


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
