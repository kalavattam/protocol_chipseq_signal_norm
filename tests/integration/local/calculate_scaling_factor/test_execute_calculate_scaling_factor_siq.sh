#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_calculate_scaling_factor_siq.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="execute calculate-scaling-factor siQ"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


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


#  Define fixture, output, and worker-input paths
scr_exe="${ROOT_REPO}/bin/execute_calculate_scaling_factor.sh"
cfg_met="${ROOT_REPO}/data/raw/docs/parse_metadata_siqchip.yml"

dir_fix="${ROOT_REPO}/tests/fixtures/calculate_scaling_factor"
dir_bam="${dir_fix}/bam"
dir_cram="${dir_fix}/cram"
dir_bam_se="${dir_bam}/se"
dir_bam_pe="${dir_bam}/pe"
dir_cram_pe="${dir_cram}/pe"
dir_met="${dir_fix}/metadata"

tbl_met="${dir_met}/measurements_siqchip.tsv"
tbl_gz="${dir_met}/measurements_siqchip.tsv.gz"
tbl_lib="${dir_met}/measurements_siqchip_lib_volume.tsv"
tbl_dup="${dir_met}/measurements_siqchip_duplicate_match.tsv"
tbl_pre="${dir_met}/measurements_siqchip_precomputed.tsv"
ref_fa="${dir_fix}/reference/tiny.fa"

tmp="${TEST_DIR_TMP}/execute_calculate_scaling_factor_siq"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/calculate_scaling_factor"

bam_pe_mip_0="${dir_bam_pe}/IP_WT_G1_Hho1_6336.sc.bam"
bam_pe_mip_1="${dir_bam_pe}/IP_WT_G1_Hho1_6337.sc.bam"
bam_pe_min_0="${dir_bam_pe}/in_WT_G1_Hho1_6336.sc.bam"
bam_pe_min_1="${dir_bam_pe}/in_WT_G1_Hho1_6337.sc.bam"

bai_pe_mip_0="${bam_pe_mip_0}.bai"
bai_pe_mip_1="${bam_pe_mip_1}.bai"
bai_pe_min_0="${bam_pe_min_0}.bai"
bai_pe_min_1="${bam_pe_min_1}.bai"

bam_pe_mip_0_hu="${dir_bam_pe}/IP_WT_G1_HU_Hho1_6336.sc.bam"
bam_pe_mip_1_hu="${dir_bam_pe}/IP_WT_G1_HU_Hho1_6337.sc.bam"
bam_pe_min_0_hu="${dir_bam_pe}/in_WT_G1_HU_Hho1_6336.sc.bam"
bam_pe_min_1_hu="${dir_bam_pe}/in_WT_G1_HU_Hho1_6337.sc.bam"

bai_pe_mip_0_hu="${bam_pe_mip_0_hu}.bai"
bai_pe_mip_1_hu="${bam_pe_mip_1_hu}.bai"
bai_pe_min_0_hu="${bam_pe_min_0_hu}.bai"
bai_pe_min_1_hu="${bam_pe_min_1_hu}.bai"

bam_se_mip_0="${dir_bam_se}/IP_WT_log_Brn1_rep1.sc.bam"
bam_se_mip_1="${dir_bam_se}/IP_WT_log_Brn1_rep2.sc.bam"
bam_se_min_0="${dir_bam_se}/in_WT_log_Brn1_rep1.sc.bam"
bam_se_min_1="${dir_bam_se}/in_WT_log_Brn1_rep2.sc.bam"

bai_se_mip_0="${bam_se_mip_0}.bai"
bai_se_mip_1="${bam_se_mip_1}.bai"
bai_se_min_0="${bam_se_min_0}.bai"
bai_se_min_1="${bam_se_min_1}.bai"

cram_pe_mip_0="${dir_cram_pe}/IP_WT_G1_Hho1_6336.sc.cram"
cram_pe_mip_1="${dir_cram_pe}/IP_WT_G1_Hho1_6337.sc.cram"
cram_pe_min_0="${dir_cram_pe}/in_WT_G1_Hho1_6336.sc.cram"
cram_pe_min_1="${dir_cram_pe}/in_WT_G1_Hho1_6337.sc.cram"

crai_pe_mip_0="${cram_pe_mip_0}.crai"
crai_pe_mip_1="${cram_pe_mip_1}.crai"
crai_pe_min_0="${cram_pe_min_0}.crai"
crai_pe_min_1="${cram_pe_min_1}.crai"

fil_out_dry="${dir_out}/scaling.dry.siq.tsv"
fil_out_pe_bam="${dir_out}/scaling.pe_bam.siq.tsv"
fil_out_cram_missing_ref="${dir_out}/scaling.cram_missing_ref.siq.tsv"
fil_out_invalid_eqn="${dir_out}/scaling.invalid_eqn.siq.tsv"
fil_out_method_not_applicable="${dir_out}/scaling.method_not_applicable.siq.tsv"
fil_out_duplicate_match="${dir_out}/scaling.duplicate_match.siq.tsv"

prt_pe_bam_0="${fil_out_pe_bam}.part.000000"
prt_pe_bam_1="${fil_out_pe_bam}.part.000001"

nam_job_dry="test_execute_calculate_scaling_factor_siq_dry"
nam_job_pe_bam="test_execute_calculate_scaling_factor_siq_pe_bam"
nam_job_pe_cram="test_execute_calculate_scaling_factor_siq_pe_cram"
nam_job_cram_missing_ref="test_execute_calculate_scaling_factor_siq_cram_missing_ref"
nam_job_invalid_eqn="test_execute_calculate_scaling_factor_siq_invalid_eqn"
nam_job_method_not_applicable="test_execute_calculate_scaling_factor_siq_method_not_applicable"
nam_job_mixed_layout="test_execute_calculate_scaling_factor_siq_mixed_layout"
nam_job_mixed_format="test_execute_calculate_scaling_factor_siq_mixed_format"
nam_job_duplicate_match="test_execute_calculate_scaling_factor_siq_duplicate_match"
nam_job_se_bam="test_execute_calculate_scaling_factor_siq_se_bam"

log_dry="${dir_log}/execute_siq_dry_run.log"
log_cram_missing_ref="${dir_log}/execute_siq_cram_missing_ref.log"
log_invalid_eqn="${dir_log}/execute_siq_invalid_eqn.log"
log_method_not_applicable="${dir_log}/execute_siq_method_not_applicable.log"
log_duplicate_match="${dir_log}/execute_siq_duplicate_match.log"
log_existing="${dir_log}/execute_siq_existing.log"
log_force_no_parts="${dir_log}/execute_siq_force_no_parts.log"

log_err_pe_bam_0="${dir_err}/${nam_job_pe_bam}.IP_WT_G1_Hho1_6336.sc.stderr.txt"
log_err_pe_bam_1="${dir_err}/${nam_job_pe_bam}.IP_WT_G1_Hho1_6337.sc.stderr.txt"
log_err_pe_cram_0="${dir_err}/${nam_job_pe_cram}.IP_WT_G1_Hho1_6336.sc.stderr.txt"
log_err_pe_cram_1="${dir_err}/${nam_job_pe_cram}.IP_WT_G1_Hho1_6337.sc.stderr.txt"
log_err_mixed_layout_0="${dir_err}/${nam_job_mixed_layout}.IP_WT_log_Brn1_rep1.sc.stderr.txt"
log_err_mixed_layout_1="${dir_err}/${nam_job_mixed_layout}.IP_WT_G1_Hho1_6337.sc.stderr.txt"
log_err_mixed_format_1="${dir_err}/${nam_job_mixed_format}.IP_WT_G1_Hho1_6337.sc.stderr.txt"
log_err_duplicate_match="${dir_err}/${nam_job_duplicate_match}.IP_WT_G1_Hho1_6336.sc.stderr.txt"
log_err_se_bam_0="${dir_err}/${nam_job_se_bam}.IP_WT_log_Brn1_rep1.sc.stderr.txt"
log_err_se_bam_1="${dir_err}/${nam_job_se_bam}.IP_WT_log_Brn1_rep2.sc.stderr.txt"

row_bam_pe_0="${bam_pe_mip_0}"$'\t'"${bam_pe_min_0}"
row_bam_pe_1="${bam_pe_mip_1}"$'\t'"${bam_pe_min_1}"
row_cram_pe_0="${cram_pe_mip_0}"$'\t'"${cram_pe_min_0}"
row_cram_pe_1="${cram_pe_mip_1}"$'\t'"${cram_pe_min_1}"
row_bam_se_0="${bam_se_mip_0}"$'\t'"${bam_se_min_0}"
row_bam_se_1="${bam_se_mip_1}"$'\t'"${bam_se_min_1}"

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
        --dir_eo "${dir_err}"
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
            --dir_eo "${dir_err}"
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
            --dir_eo "${dir_err}"
            --max_job 1
    )

    arr_cmd_lib=(
        "${TEST_BASH}" "${scr_exe}"
            --threads 1
            --mode siq
            --csv_mip "${bam_pe_mip_0},${bam_pe_mip_1}"
            --csv_min "${bam_pe_min_0},${bam_pe_min_1}"
            --aln_typ auto
            --tbl_met "${tbl_lib}"
            --cfg_met "${cfg_met}"
            --eqn 6nd
            --dir_eo "${dir_err}"
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
            --dir_eo "${dir_err}"
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
        "${scr_exe}" \
        "${cfg_met}" \
        "${tbl_met}" \
        "${tbl_gz}" \
        "${tbl_lib}" \
        "${tbl_dup}" \
        "${tbl_pre}" \
        "${bam_pe_mip_0}" \
        "${bai_pe_mip_0}" \
        "${bam_pe_mip_1}" \
        "${bai_pe_mip_1}" \
        "${bam_pe_min_0}" \
        "${bai_pe_min_0}" \
        "${bam_pe_min_1}" \
        "${bai_pe_min_1}" \
        "${bam_pe_mip_0_hu}" \
        "${bai_pe_mip_0_hu}" \
        "${bam_pe_mip_1_hu}" \
        "${bai_pe_mip_1_hu}" \
        "${bam_pe_min_0_hu}" \
        "${bai_pe_min_0_hu}" \
        "${bam_pe_min_1_hu}" \
        "${bai_pe_min_1_hu}" \
        "${cram_pe_mip_0}" \
        "${crai_pe_mip_0}" \
        "${cram_pe_mip_1}" \
        "${crai_pe_mip_1}" \
        "${cram_pe_min_0}" \
        "${crai_pe_min_0}" \
        "${cram_pe_min_1}" \
        "${crai_pe_min_1}" \
        "${ref_fa}" \
        "${bam_se_mip_0}" \
        "${bai_se_mip_0}" \
        "${bam_se_mip_1}" \
        "${bai_se_mip_1}" \
        "${bam_se_min_0}" \
        "${bai_se_min_0}" \
        "${bam_se_min_1}" \
        "${bai_se_min_1}"
then
    finish
    exit $?
fi


#  Run one serial siQ calculation case and assert its assembled table
#  Dry-run execution should report worker, combiner, and header commands
if \
    run_capture \
        "execute calculate-scaling-factor siQ dry-run" \
        "${log_dry}" \
        "${arr_cmd_bam_pe[@]}" \
        --env_nam "${env_nam}" \
        --fil_out "${fil_out_dry}" \
        --nam_job "${nam_job_dry}" \
        --dry_run
then
    record_pass "execute_calculate_scaling_factor.sh siQ --dry_run exits 0"
else
    record_fail \
        "execute_calculate_scaling_factor.sh siQ --dry_run failed; see" \
        "$(print_relpath "${log_dry}")"
fi

if [[ ! -e "${fil_out_dry}" ]]; then
    record_pass "execute_calculate_scaling_factor.sh siQ --dry_run writes no TSV"
else
    record_fail "execute_calculate_scaling_factor.sh siQ --dry_run wrote TSV"
fi

assert_pattern_found \
    "${log_dry}" \
    'combine_parts_scaling_factor.sh' \
    "execute_calculate_scaling_factor.sh siQ --dry_run reports combiner"

assert_pattern_found \
    "${log_dry}" \
    'write_header.sh' \
    "execute_calculate_scaling_factor.sh siQ --dry_run reports header"

assert_pattern_found \
    "${log_dry}" \
    'write_header.sh.*--fil_in.*--in_place' \
    "execute_calculate_scaling_factor.sh siQ --dry_run reports in-place headering"

assert_pattern_found \
    "${log_dry}" \
    "${TEST_BASH} ${ROOT_REPO}/bin/submit_calculate_scaling_factor.sh" \
    "execute_calculate_scaling_factor.sh siQ --dry_run reports submit"

assert_pattern_found \
    "${log_dry}" \
    '--idx_out 0' \
    "execute_calculate_scaling_factor.sh siQ --dry_run reports first index"

assert_pattern_found \
    "${log_dry}" \
    '--idx_out 1' \
    "execute_calculate_scaling_factor.sh siQ --dry_run reports second index"


#  PE BAM input should compute lengths from paired-end fragment sizes
run_case_siq \
    pe_bam \
    arr_cmd_bam_pe \
    "${row_bam_pe_0}" \
    "${row_bam_pe_1}" \
    $'0.001912211397724232486706\t6nd\t2.7\t72.5\t300\t20\t3\t2\t626\t450' \
    $'0.002902612244746139314594\t6nd\t5\t81.1\t300\t20\t2\t3\t663\t437'

assert_pattern_found \
    "${log_err_pe_bam_0}" \
    'typ_ip=pe' \
    "execute scaling-factor siQ PE BAM auto-detects first IP as PE"

assert_pattern_found \
    "${log_err_pe_bam_1}" \
    'typ_ip=pe' \
    "execute scaling-factor siQ PE BAM auto-detects second IP as PE"


#  Metadata library-loading volumes should correct alpha through execute
run_case_scaling_factor_execute \
    pe_lib_volume \
    siQ \
    arr_cmd_lib \
    "${dir_out}" \
    "${dir_log}" \
    siq \
    test_execute_calculate_scaling_factor_siq \
    $'^fil_ip\tfil_in\tsiq\teqn\tmass_ip\tmass_in\tvol_all\tvol_in\tdep_ip\tdep_in\tlen_ip\tlen_in\tlib_vol_ip\tlib_vol_in$' \
    "${row_bam_pe_0}" \
    "${row_bam_pe_1}" \
    $'0.003824422795448464973411\t6nd\t2.7\t72.5\t300\t20\t3\t2\t626\t450\t2\t4' \
    $'0.005805224489492278629188\t6nd\t5\t81.1\t300\t20\t2\t3\t663\t437\t2\t4'


#  PE CRAM input should work when an explicit reference is supplied
run_case_siq \
    pe_cram \
    arr_cmd_cram_pe \
    "${row_cram_pe_0}" \
    "${row_cram_pe_1}" \
    $'0.001912211397724232486706\t6nd\t2.7\t72.5\t300\t20\t3\t2\t626\t450' \
    $'0.002902612244746139314594\t6nd\t5\t81.1\t300\t20\t2\t3\t663\t437'

assert_pattern_found \
    "${log_err_pe_cram_0}" \
    'typ_ip=pe' \
    "execute scaling-factor siQ PE CRAM auto-detects first IP as PE"

assert_pattern_found \
    "${log_err_pe_cram_1}" \
    "ref_fa=${ref_fa}" \
    "execute scaling-factor siQ PE CRAM forwards ref_fa"


#  CRAM input should fail clearly when no reference is supplied
if \
    run_capture \
        "execute calculate-scaling-factor siQ CRAM missing reference" \
        "${log_cram_missing_ref}" \
        "${TEST_BASH}" "${scr_exe}" \
            --env_nam "${env_nam}" \
            --threads 1 \
            --mode siq \
            --csv_mip "${cram_pe_mip_0},${cram_pe_mip_1}" \
            --csv_min "${cram_pe_min_0},${cram_pe_min_1}" \
            --aln_typ auto \
            --tbl_met "${tbl_met}" \
            --cfg_met "${cfg_met}" \
            --eqn 6nd \
            --fil_out "${fil_out_cram_missing_ref}" \
            --dir_eo "${dir_err}" \
            --nam_job "${nam_job_cram_missing_ref}" \
            --max_job 1
then
    record_fail \
        "execute_calculate_scaling_factor.sh siQ CRAM without reference" \
        "unexpectedly passed"
else
    assert_pattern_found \
        "${log_cram_missing_ref}" \
        "'--ref_fa' is required" \
        "execute_calculate_scaling_factor.sh rejects CRAM without ref_fa"
fi


#  Invalid equation identifiers should fail before processing
if \
    run_capture \
        "execute calculate-scaling-factor siQ invalid equation" \
        "${log_invalid_eqn}" \
        "${TEST_BASH}" "${scr_exe}" \
            --env_nam "${env_nam}" \
            --threads 1 \
            --mode siq \
            --csv_mip "${bam_pe_mip_0},${bam_pe_mip_1}" \
            --csv_min "${bam_pe_min_0},${bam_pe_min_1}" \
            --aln_typ auto \
            --tbl_met "${tbl_met}" \
            --cfg_met "${cfg_met}" \
            --eqn 7 \
            --fil_out "${fil_out_invalid_eqn}" \
            --dir_eo "${dir_err}" \
            --nam_job "${nam_job_invalid_eqn}" \
            --max_job 1
then
    record_fail \
        "execute_calculate_scaling_factor.sh siQ invalid equation" \
        "unexpectedly passed"
else
    assert_pattern_found \
        "${log_invalid_eqn}" \
        "equation ('--eqn') was assigned '7'" \
        "execute_calculate_scaling_factor.sh rejects invalid siQ equation"
fi


#  Spike-in method arguments are invalid in siQ mode
if \
    run_capture \
        "execute calculate-scaling-factor siQ method not applicable" \
        "${log_method_not_applicable}" \
        "${arr_cmd_bam_pe[@]}" \
        --env_nam "${env_nam}" \
        --method fractional \
        --fil_out "${fil_out_method_not_applicable}" \
        --nam_job "${nam_job_method_not_applicable}"
then
    record_fail \
        "execute_calculate_scaling_factor.sh siQ method argument" \
        "unexpectedly passed"
else
    assert_pattern_found \
        "${log_method_not_applicable}" \
        "'--method' may be used only when '--mode spike' is active" \
        "execute_calculate_scaling_factor.sh rejects method under siQ mode"
fi


#  Mixed SE/PE BAM lists should keep per-file auto-detection
run_case_siq \
    mixed_layout \
    arr_cmd_mxd_lyt \
    "${row_bam_se_0}" \
    "${row_bam_pe_1}" \
    $'0.004761904761904762334312\t6nd\t4\t60\t300\t20\t3\t2\t150\t150' \
    $'0.002902612244746139314594\t6nd\t5\t81.1\t300\t20\t2\t3\t663\t437'

assert_pattern_found \
    "${log_err_mixed_layout_0}" \
    'typ_ip=se' \
    "execute scaling-factor siQ mixed layout keeps first IP as SE"

assert_pattern_found \
    "${log_err_mixed_layout_1}" \
    'typ_ip=pe' \
    "execute scaling-factor siQ mixed layout keeps second IP as PE"


#  Mixed BAM/CRAM lists should route CRAM samples with the shared reference
run_case_siq \
    mixed_format \
    arr_cmd_mxd_fmt \
    "${row_bam_pe_0}" \
    "${row_cram_pe_1}" \
    $'0.001912211397724232486706\t6nd\t2.7\t72.5\t300\t20\t3\t2\t626\t450' \
    $'0.002902612244746139314594\t6nd\t5\t81.1\t300\t20\t2\t3\t663\t437'

assert_pattern_found \
    "${log_err_mixed_format_1}" \
    "ref_fa=${ref_fa}" \
    "execute scaling-factor siQ mixed BAM/CRAM list forwards ref_fa"


#  Metadata parser variants should propagate through the execute wrapper
run_case_siq \
    gzip_metadata \
    arr_cmd_bam_pe \
    "${row_bam_pe_0}" \
    "${row_bam_pe_1}" \
    $'0.001912211397724232486706\t6nd\t2.7\t72.5\t300\t20\t3\t2\t626\t450' \
    $'0.002902612244746139314594\t6nd\t5\t81.1\t300\t20\t2\t3\t663\t437' \
    --tbl_met "${tbl_gz}"

if \
    run_capture \
        "execute calculate-scaling-factor siQ duplicate match" \
        "${log_duplicate_match}" \
        "${TEST_BASH}" "${scr_exe}" \
            --env_nam "${env_nam}" \
            --threads 1 \
            --mode siq \
            --csv_mip "${bam_pe_mip_0}" \
            --csv_min "${bam_pe_min_0}" \
            --aln_typ auto \
            --tbl_met "${tbl_dup}" \
            --cfg_met "${cfg_met}" \
            --eqn 6nd \
            --fil_out "${fil_out_duplicate_match}" \
            --dir_eo "${dir_err}" \
            --nam_job "${nam_job_duplicate_match}" \
            --max_job 1
then
    record_fail \
        "execute_calculate_scaling_factor.sh siQ duplicate match" \
        "unexpectedly passed"
else
    assert_pattern_found \
        "${log_err_duplicate_match}" \
        "Multiple metadata rows matched" \
        "execute_calculate_scaling_factor.sh siQ duplicate match fails clearly"
fi


#  Equation selection should propagate through the execute wrapper
run_case_siq \
    eqn5 \
    arr_cmd_bam_pe \
    "${row_bam_pe_0}" \
    "${row_bam_pe_1}" \
    $'0.001189820425250633388267\t5\t2.7\t72.5\t300\t20\t3\t2\t626\t450' \
    $'0.004063657142644594953695\t5\t5\t81.1\t300\t20\t2\t3\t663\t437' \
    --eqn 5

run_case_siq \
    eqn5nd \
    arr_cmd_bam_pe \
    "${row_bam_pe_0}" \
    "${row_bam_pe_1}" \
    $'0.001784730637875950407661\t5nd\t2.7\t72.5\t300\t20\t3\t2\t626\t450' \
    $'0.002709104761763063591584\t5nd\t5\t81.1\t300\t20\t2\t3\t663\t437' \
    --eqn 5nd

run_case_siq \
    eqn6 \
    arr_cmd_bam_pe \
    "${row_bam_pe_0}" \
    "${row_bam_pe_1}" \
    $'0.001274807598482821657804\t6\t2.7\t72.5\t300\t20\t3\t2\t626\t450' \
    $'0.004353918367119209188731\t6\t5\t81.1\t300\t20\t2\t3\t663\t437' \
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
    $'0.00382442279544846453973\t6\t2.7\t72.5\t300\t20\t10\t20\t626\t450' \
    $'0.005805224489492278629188\t6\t5\t81.1\t300\t20\t10\t20\t663\t437' \
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

assert_pattern_found \
    "${log_err_se_bam_0}" \
    'typ_ip=se' \
    "execute scaling-factor siQ SE BAM auto-detects first IP as SE"

assert_pattern_found \
    "${log_err_se_bam_1}" \
    'typ_ip=se' \
    "execute scaling-factor siQ SE BAM auto-detects second IP as SE"

run_case_siq \
    no_header \
    arr_cmd_bam_pe \
    "${row_bam_pe_0}" \
    "${row_bam_pe_1}" \
    $'0.001912211397724232486706\t6nd\t2.7\t72.5\t300\t20\t3\t2\t626\t450' \
    $'0.002902612244746139314594\t6nd\t5\t81.1\t300\t20\t2\t3\t663\t437' \
    --no_header


#  Existing final output should fail unless '--force' is supplied
if \
    run_capture \
        "execute calculate-scaling-factor siQ existing output" \
        "${log_existing}" \
        "${arr_cmd_bam_pe[@]}" \
        --env_nam "${env_nam}" \
        --fil_out "${fil_out_pe_bam}" \
        --nam_job "${nam_job_pe_bam}"
then
    record_fail \
        "execute_calculate_scaling_factor.sh siQ existing output" \
        "unexpectedly passed"
else
    assert_pattern_found \
        "${log_existing}" \
        'output file already exists' \
        "execute_calculate_scaling_factor.sh siQ rejects existing output"
fi


#  '--force' should replace output and '--no_parts' should remove worker parts
if \
    run_capture \
        "execute calculate-scaling-factor siQ force no-parts" \
        "${log_force_no_parts}" \
        "${arr_cmd_bam_pe[@]}" \
        --env_nam "${env_nam}" \
        --fil_out "${fil_out_pe_bam}" \
        --nam_job "${nam_job_pe_bam}" \
        --force \
        --no_parts
then
    record_pass \
        "execute_calculate_scaling_factor.sh siQ --force --no_parts exits 0"
else
    record_fail \
        "execute_calculate_scaling_factor.sh siQ --force --no_parts failed;" \
        "see $(print_relpath "${log_force_no_parts}")"
fi

assert_file_nonempty \
    "${fil_out_pe_bam}" \
    "execute scaling-factor siQ replaced final TSV"

if [[ ! -e "${prt_pe_bam_0}" && ! -e "${prt_pe_bam_1}" ]]; then
    record_pass "execute scaling-factor siQ --no_parts removes worker parts"
else
    record_fail "execute scaling-factor siQ --no_parts retained worker parts"
fi

finish
