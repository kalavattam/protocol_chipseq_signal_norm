#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_execute_compute_signal_bam.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6);
# - Anthropic Claude Code (Opus 5).
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="execute compute-signal BAM"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Define fixture and output paths for the execute-to-submit BAM path.
dir_fx="${ROOT_REPO}/tests/fixtures/compute_signal"
in_se="${dir_fx}/bam/se/tiny_se.bam"
in_pe="${dir_fx}/bam/pe/tiny_pe.bam"

tmp="${TEST_DIR_TMP}/execute_compute_signal_bam"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"
dir_log="${TEST_DIR_LOG}/compute_signal"

fil_out_se_signal="${dir_out}/tiny_se.bdg"
log_se_signal="${dir_log}/execute_compute_signal_bam_se_signal.log"

fil_out_se_coord="${dir_out}/tiny_se.bed"
log_se_coord="${dir_log}/execute_compute_signal_bam_se_coord.log"

fil_out_pe_signal="${dir_out}/tiny_pe.bdg"
log_pe_signal="${dir_log}/execute_compute_signal_bam_pe_signal.log"

fil_out_pe_coord="${dir_out}/tiny_pe.bed"
log_pe_coord="${dir_log}/execute_compute_signal_bam_pe_coord.log"

fil_out_se_signal_scaled="${dir_out}/scaled.tiny_se.bdg"
log_se_signal_scaled="${dir_log}/execute_compute_signal_bam_se_signal_scaled.log"

fil_out_se_signal_usr_frg="${dir_out}/usr_frg_signal.tiny_se.bdg"
log_se_signal_usr_frg="${dir_log}/execute_compute_signal_bam_se_signal_usr_frg.log"

fil_out_se_coord_usr_frg="${dir_out}/usr_frg_coord.tiny_se.bed"
log_se_coord_usr_frg="${dir_log}/execute_compute_signal_bam_se_coord_usr_frg.log"


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_err}" "${dir_log}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty \
    "${in_se}" \
    "${in_pe}" || {
    finish
    exit $?
}


# Signal mode: two 10-bp SE alignments produce chromosome-I bedGraph bins.
run_case_compute_signal \
    execute \
    bam \
    "se_signal" \
    "signal" \
    "${in_se}" \
    "bdg" \
    "${log_se_signal}" \
    "${dir_out}" \
    "${dir_err}" \
    "" \
    --method unadj \
    --siz_bin 10 \
    --engine window \
    --siz_win 20 \
    --csv_scl_fct NA

assert_file_nonempty \
    "${fil_out_se_signal}" \
    "execute SE signal bedGraph output"

if [[ -s "${fil_out_se_signal}" ]]; then
    assert_pattern_found \
        "${fil_out_se_signal}" \
        $'^I\t0\t10\t10$' \
        "execute SE signal output has chromosome-I bin I:0-10"

    assert_pattern_found \
        "${fil_out_se_signal}" \
        $'^I\t20\t30\t10$' \
        "execute SE signal output has chromosome-I bin I:20-30"
fi


# Coord mode: the same SE alignments emit BED-like processed fragments.
run_case_compute_signal \
    execute \
    bam \
    "se_coord" \
    "coord" \
    "${in_se}" \
    "bed" \
    "${log_se_coord}" \
    "${dir_out}" \
    "${dir_err}" \
    ""

assert_file_nonempty \
    "${fil_out_se_coord}" \
    "execute SE coord BED output"

if [[ -s "${fil_out_se_coord}" ]]; then
    assert_pattern_found \
        "${fil_out_se_coord}" \
        $'^I\t0\t10\t10$' \
        "execute SE coord output has chromosome-I fragment I:0-10"

    assert_pattern_found \
        "${fil_out_se_coord}" \
        $'^I\t20\t30\t10$' \
        "execute SE coord output has chromosome-I fragment I:20-30"
fi


# Signal mode: two PE fragments cover bins from I:10-60.
run_case_compute_signal \
    execute \
    bam \
    "pe_signal" \
    "signal" \
    "${in_pe}" \
    "bdg" \
    "${log_pe_signal}" \
    "${dir_out}" \
    "${dir_err}" \
    "" \
    --method unadj \
    --siz_bin 10 \
    --csv_scl_fct NA

assert_file_nonempty \
    "${fil_out_pe_signal}" \
    "execute PE signal bedGraph output"

if [[ -s "${fil_out_pe_signal}" ]]; then
    assert_pattern_found \
        "${fil_out_pe_signal}" \
        $'^I\t10\t20\t10$' \
        "execute PE signal output has chromosome-I bin I:10-20"

    assert_pattern_found \
        "${fil_out_pe_signal}" \
        $'^I\t40\t50\t10$' \
        "execute PE signal output has chromosome-I bin I:40-50"
fi


# Coord mode: PE output emits one BED-like row per leftmost proper pair.
run_case_compute_signal \
    execute \
    bam \
    "pe_coord" \
    "coord" \
    "${in_pe}" \
    "bed" \
    "${log_pe_coord}" \
    "${dir_out}" \
    "${dir_err}" \
    ""

assert_file_nonempty \
    "${fil_out_pe_coord}" \
    "execute PE coord BED output"

if [[ -s "${fil_out_pe_coord}" ]]; then
    assert_pattern_found \
        "${fil_out_pe_coord}" \
        $'^I\t10\t40\t30$' \
        "execute PE coord output has chromosome-I fragment I:10-40"

    assert_pattern_found \
        "${fil_out_pe_coord}" \
        $'^I\t40\t60\t20$' \
        "execute PE coord output has chromosome-I fragment I:40-60"
fi


# Scaling factor and prefix propagation: raw 10-bp SE bins scaled by 2.
run_case_compute_signal \
    execute \
    bam \
    "se_signal_scaled" \
    "signal" \
    "${in_se}" \
    "bdg" \
    "${log_se_signal_scaled}" \
    "${dir_out}" \
    "${dir_err}" \
    "" \
    --method unadj \
    --siz_bin 10 \
    --csv_scl_fct 2 \
    --prefix scaled

assert_file_nonempty \
    "${fil_out_se_signal_scaled}" \
    "execute scaled SE signal output"

if [[ -s "${fil_out_se_signal_scaled}" ]]; then
    assert_pattern_found \
        "${fil_out_se_signal_scaled}" \
        $'^I\t0\t10\t20$' \
        "execute scaled SE signal output has I:0-10 = 20"

    assert_pattern_found \
        "${fil_out_se_signal_scaled}" \
        $'^I\t20\t30\t20$' \
        "execute scaled SE signal output has I:20-30 = 20"
fi


# Fixed SE fragment length extends reads to 20 bp in signal mode.
run_case_compute_signal \
    execute \
    bam \
    "se_signal_usr_frg" \
    "signal" \
    "${in_se}" \
    "bdg" \
    "${log_se_signal_usr_frg}" \
    "${dir_out}" \
    "${dir_err}" \
    "" \
    --method unadj \
    --siz_bin 10 \
    --csv_scl_fct NA \
    --csv_usr_frg 20 \
    --prefix usr_frg_signal

assert_file_nonempty \
    "${fil_out_se_signal_usr_frg}" \
    "execute usr_frg SE signal output"

if [[ -s "${fil_out_se_signal_usr_frg}" ]]; then
    assert_pattern_found \
        "${fil_out_se_signal_usr_frg}" \
        $'^I\t0\t10\t10$' \
        "execute usr_frg SE signal output has I:0-10 = 10"

    assert_pattern_found \
        "${fil_out_se_signal_usr_frg}" \
        $'^I\t10\t20\t20$' \
        "execute usr_frg SE signal output has I:10-20 = 20"

    assert_pattern_found \
        "${fil_out_se_signal_usr_frg}" \
        $'^I\t20\t30\t10$' \
        "execute usr_frg SE signal output has I:20-30 = 10"
fi


# Fixed SE fragment length is reflected in coord-mode BED intervals.
run_case_compute_signal \
    execute \
    bam \
    "se_coord_usr_frg" \
    "coord" \
    "${in_se}" \
    "bed" \
    "${log_se_coord_usr_frg}" \
    "${dir_out}" \
    "${dir_err}" \
    "" \
    --csv_usr_frg 20 \
    --prefix usr_frg_coord

assert_file_nonempty \
    "${fil_out_se_coord_usr_frg}" \
    "execute usr_frg SE coord output"

if [[ -s "${fil_out_se_coord_usr_frg}" ]]; then
    assert_pattern_found \
        "${fil_out_se_coord_usr_frg}" \
        $'^I\t0\t20\t20$' \
        "execute usr_frg SE coord output has chromosome-I fragment I:0-20"

    assert_pattern_found \
        "${fil_out_se_coord_usr_frg}" \
        $'^I\t10\t30\t20$' \
        "execute usr_frg SE coord output has chromosome-I fragment I:10-30"
fi


# Forwarding contract for '--siz_win' and the retired '--chr_siz'. These assert
# against the command the wrapper emitted, not merely that a run succeeded: a
# silently dropped option would still produce output.
log_fwd="${tmp}/logs/test_execute_compute_bam_se_signal.tiny_se.stderr.txt"

if [[ -s "${log_fwd}" ]]; then
    assert_pattern_found \
        "${log_fwd}" \
        "--siz_win *20" \
        "execute signal forwards the supplied '--siz_win'"

    assert_pattern_absent \
        "${log_fwd}" \
        "--chr_siz" \
        "execute signal no longer forwards '--chr_siz'"
fi

log_coord="${tmp}/logs/test_execute_compute_bam_se_coord.tiny_se.stderr.txt"

if [[ -s "${log_coord}" ]]; then
    assert_pattern_absent \
        "${log_coord}" \
        "--siz_win" \
        "execute coord omits '--siz_win'"

    assert_pattern_absent \
        "${log_coord}" \
        "--engine" \
        "execute coord omits '--engine'"
fi


# Validation contract for '--siz_win': a nonpositive or nonnumeric value is
# rejected before any job is built.
for val_bad in 0 abc; do
    out_bad="$(
        bash "${ROOT_REPO}/bin/execute_compute_signal.sh" \
            --dry_run \
            --mode signal \
            --csv_fil_in "${in_se}" \
            --dir_out "${dir_out}" \
            --dir_eo "${dir_err}" \
            --siz_bin 10 \
            --siz_win "${val_bad}" 2>&1
    )" || true

    if [[ "${out_bad}" == *"'--siz_win' was assigned"* ]]; then
        record_pass "execute rejects an invalid '--siz_win ${val_bad}'"
    else
        record_fail "execute did not reject '--siz_win ${val_bad}' itself"
    fi
done


# Default contract: a signal case that names neither option still forwards both
# wrapper defaults, so a dropped default cannot pass unnoticed.
log_dflt="${tmp}/logs/test_execute_compute_bam_pe_signal.tiny_pe.stderr.txt"
if [[ -s "${log_dflt}" ]]; then
    assert_pattern_found \
        "${log_dflt}" \
        "--engine *chrom" \
        "execute forwards the default '--engine chrom'"

    assert_pattern_found \
        "${log_dflt}" \
        "--siz_win *100000" \
        "execute forwards the default '--siz_win 100000'"
fi


# Validation contract for '--engine': only 'chrom' and 'window' are accepted.
for eng_bad in chrm windowed bogus; do
    out_bad="$(
        bash "${ROOT_REPO}/bin/execute_compute_signal.sh" \
            --dry_run \
            --mode signal \
            --csv_fil_in "${in_se}" \
            --dir_out "${dir_out}" \
            --dir_eo "${dir_err}" \
            --siz_bin 10 \
            --engine "${eng_bad}" 2>&1
    )" || true

    if [[ "${out_bad}" == *"'--engine' must be 'chrom' or 'window'"* ]]; then
        record_pass "execute rejects an invalid '--engine ${eng_bad}'"
    else
        record_fail "execute did not reject '--engine ${eng_bad}' itself"
    fi
done


# Report flags: the wrapper derives '<track>.N.txt' and '<track>.L.txt' from
# each output name, so assert those paths and their contents, not just success.
dir_rep="${tmp}/reports"
mkdir -p "${dir_rep}"

bash "${ROOT_REPO}/bin/execute_compute_signal.sh" \
    --mode signal \
    --csv_fil_in "${in_se},${in_pe}" \
    --dir_out "${dir_rep}" \
    --dir_eo "${dir_err}" \
    --typ_out bedGraph \
    --siz_bin 10 \
    --method unadj \
    --report_N \
    --report_L \
    > /dev/null 2>&1 || true

for samp in tiny_se tiny_pe; do
    assert_file_nonempty \
        "${dir_rep}/${samp}.bedGraph" \
        "execute report run still writes the ${samp} track"

    assert_file_nonempty \
        "${dir_rep}/${samp}.N.txt" \
        "execute derives ${samp}.N.txt beside the track"

    assert_file_nonempty \
        "${dir_rep}/${samp}.L.txt" \
        "execute derives ${samp}.L.txt beside the track"
done

# SE holds two 10-bp fragments, one bin each; PE spans five bins. Distinct
# values prove each sample got its own report path rather than the first.
assert_file_exact_line "${dir_rep}/tiny_se.N.txt" "2" \
    "execute SE fragment count is 2"
assert_file_exact_line "${dir_rep}/tiny_se.L.txt" "2" \
    "execute SE spanned-bin count is 2"
assert_file_exact_line "${dir_rep}/tiny_pe.L.txt" "5" \
    "execute PE spanned-bin count is 5, distinct from SE"

dir_only="${tmp}/report_only"
mkdir -p "${dir_only}"

bash "${ROOT_REPO}/bin/execute_compute_signal.sh" \
    --mode signal \
    --csv_fil_in "${in_se}" \
    --dir_out "${dir_only}" \
    --dir_eo "${dir_err}" \
    --typ_out bedGraph \
    --siz_bin 10 \
    --method unadj \
    --report_N \
    --report_L \
    --report_only \
    > /dev/null 2>&1 || true

assert_file_nonempty \
    "${dir_only}/tiny_se.N.txt" \
    "execute report-only writes the fragment count"

if [[ -e "${dir_only}/tiny_se.bedGraph" ]]; then
    record_fail "execute report-only unexpectedly wrote a track"
else
    record_pass "execute report-only writes no track"
fi

out_bad="$(
    bash "${ROOT_REPO}/bin/execute_compute_signal.sh" \
        --mode signal \
        --csv_fil_in "${in_se}" \
        --dir_out "${dir_only}" \
        --dir_eo "${dir_err}" \
        --siz_bin 10 \
        --report_only 2>&1
)" || true

if [[ "${out_bad}" == *"'--report_only' requires"* ]]; then
    record_pass "execute rejects '--report_only' with nothing to report"
else
    record_fail "execute accepted '--report_only' with no report flag"
fi


finish
