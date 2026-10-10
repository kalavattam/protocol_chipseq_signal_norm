#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_check_formats.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="file-format recognition helpers and their callers"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Define fixture and output paths. Format names and file suffixes match in any
# letter case; the short bedGraph spellings 'bdg' and 'bg' are refused.
dir_fx="${ROOT_REPO}/tests/fixtures/compute_signal"
in_se="${dir_fx}/bam/se/tiny_se.bam"
fil_A="${dir_fx}/bedgraph/ratio_A.bedGraph.gz"
fil_B="${dir_fx}/bedgraph/ratio_B.bedGraph.gz"

tmp="${TEST_DIR_TMP}/check_formats"
dir_out="${tmp}/out"
dir_err="${tmp}/logs"


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${dir_out}" "${dir_err}"

require_env_project env_nam || {
    finish
    exit $?
}

require_files_nonempty "${in_se}" "${fil_A}" "${fil_B}" || {
    finish
    exit $?
}


# Call one helper in a fresh shell and print its stdout, stderr, and status.
function run_fmt() {
    # shellcheck disable=SC2016  # Expand in the child shell, not this one.
    "${TEST_BASH}" -c '
        # shellcheck source=lib/bash/core/source_helpers.sh
        source "${1}/lib/bash/core/source_helpers.sh"

        source_helpers "${1}/lib/bash" core/check_formats
        shift 1

        rc=0
        "$@" || rc=$?
        printf "rc=%s\n" "${rc}"
    ' _ "${ROOT_REPO}" "$@" 2>&1
}


# Helper rows: function|argument|status|output. A refusal row expects its error
# text; an unrecognized row expects no output.
rows_fn=(
    "canonicalize_fmt|bedGraph|0|bedGraph"
    "canonicalize_fmt|BEDGRAPH.gz|0|bedGraph.gz"
    "canonicalize_fmt|bEdGrApH|0|bedGraph"
    "canonicalize_fmt|CRAM|0|cram"
    "canonicalize_fmt|bed.GZ|0|bed.gz"
    "canonicalize_fmt|bdg|1|'bdg' is not accepted; use 'bedGraph'."
    "canonicalize_fmt|bg.gz|1|'bg.gz' is not accepted; use 'bedGraph.gz'."
    "canonicalize_fmt|bam.gz|2|"
    "canonicalize_fmt|txt|2|"
    "fmt_of_path|x.CrAm|0|cram"
    "fmt_of_path|dir/x.BEDGRAPH.gz|0|bedGraph.gz"
    "fmt_of_path|x.bed.gz|0|bed.gz"
    "fmt_of_path|x.bam.gz|2|"
    "fmt_of_path|d.d/noext|2|"
    "fmt_of_path|x.BDG|1|'.BDG' is not accepted; the file name must end in"
    "check_fmt_path|x.bedGraph|0|"
    "check_fmt_path|/dev/fd/63|0|"
    "check_fmt_path|x.bg.gz|1|end in '.bedGraph.gz'."
    "strip_fmt_suffix|IP_x.BEDGRAPH.gz|0|IP_x"
    "strip_fmt_suffix|x.CRAM|0|x"
    "strip_fmt_suffix|x.SAM|0|x"
    "strip_fmt_suffix|x.txt|0|x.txt"
)

for row in "${rows_fn[@]}"; do
    IFS='|' read -r fnc arg rc_exp out_exp <<< "${row}"
    out="$(run_fmt "${fnc}" "${arg}")"
    lbl="${fnc} '${arg}'"

    # A success prints exactly the expected output; a refusal contains its
    # message.
    want="rc=${rc_exp}"
    if [[ -n "${out_exp}" && "${rc_exp}" == "0" ]]; then
        want="${out_exp}"$'\n'"${want}"
    fi

    if [[
        "${out}" == "${want}"
        || (
            "${rc_exp}" != "0"
            && -n "${out_exp}"
            && "${out}" == *"${out_exp}"*
            && "${out}" == *"rc=${rc_exp}"
        )
    ]]; then
        record_pass "${lbl} gives status ${rc_exp}"
    else
        record_fail "${lbl} gave '${out//$'\n'/ | }'; expected ${rc_exp}"
    fi
done


# Dry-run 'execute' with one set of arguments, printing all output.
function run_execute() {
    "${TEST_BASH}" "${ROOT_REPO}/bin/execute_compute_signal.sh" \
        --dry_run \
        --dir_out "${dir_out}" \
        --dir_eo "${dir_err}" \
        --max_job 1 \
        --threads 1 \
        "$@" 2>&1
}

# A short '--typ_out' is refused before any work.
rc=0
out="$(run_execute --mode signal --csv_fil_in "${in_se}" --typ_out bdg.gz)" \
    || rc=$?

if [[
    "${rc}" -ne 0
    && "${out}" == *"'bdg.gz' is not accepted; use 'bedGraph.gz'."*
]]; then
    record_pass "execute refuses '--typ_out bdg.gz'"
else
    record_fail "execute did not refuse '--typ_out bdg.gz' (exit ${rc})"
fi

# A '--typ_out' in another letter case is used in its canonical spelling.
rc=0
out="$(run_execute --mode signal --csv_fil_in "${in_se}" --typ_out BEDGRAPH)" \
    || rc=$?
cmd="$(grep -m1 'submit_compute_signal.sh --env' <<< "${out}" || true)"

if [[
    "${rc}" -eq 0
    && "${cmd}" == *"tiny_se.bedGraph "*
    && "${out}" != *"Coercing"*
]]; then
    record_pass "execute writes '--typ_out BEDGRAPH' as '.bedGraph'"
else
    record_fail "execute did not normalize '--typ_out BEDGRAPH' (exit ${rc})"
fi

# A ratio input with a short suffix is refused by name, before any check of the
# file itself.
rc=0
out="$(
    run_execute --mode ratio --csv_fil_A "${tmp}/IP.bdg" --csv_fil_B "${fil_B}"
)" || rc=$?

if [[
    "${rc}" -ne 0
    && "${out}" == *"'.bdg' is not accepted; the file name must end in"*
]]; then
    record_pass "execute refuses a '.bdg' ratio input"
else
    record_fail "execute did not refuse a '.bdg' ratio input (exit ${rc})"
fi


# Run 'submit' on a ratio pair, writing to its own output directory.
function run_submit() {
    "${TEST_BASH}" "${ROOT_REPO}/bin/submit_compute_signal.sh" \
        --env_nam "${env_nam}" \
        --dir_scr "${ROOT_REPO}/bin" \
        --threads 1 \
        --dir_eo "${dir_err}" \
        --mode ratio \
        "$@" 2>&1
}

# 'submit' refuses a short suffix on an input or an output, before any work:
# label|input A|output|message.
rows_sub=(
    "input|${tmp}/IP.bg.gz|${dir_out}/r.bedGraph|'.bg.gz' is not accepted"
    "output|${fil_A}|${dir_out}/r.bdg|'.bdg' is not accepted"
)

for row in "${rows_sub[@]}"; do
    IFS='|' read -r lbl in_A fil_out msg <<< "${row}"

    rc=0
    out="$(
        run_submit \
            --csv_fil_A "${in_A}" \
            --csv_fil_B "${fil_B}" \
            --csv_fil_out "${fil_out}"
    )" || rc=$?

    if [[
        "${rc}" -ne 0
        && "${out}" == *"${msg}"*
        && "${out}" != *"run_py compute_signal_ratio"*
    ]]; then
        record_pass "submit refuses a short-suffix ${lbl} before any work"
    else
        record_fail "submit did not refuse a short-suffix ${lbl} (exit ${rc})"
    fi
done


# Callers recognize formats in any letter case. Links with upper-case suffixes
# stand in for files a user named that way.
dir_up="${tmp}/upper"
dir_cs="${ROOT_REPO}/tests/fixtures/calculate_scaling_factor"
in_cram="${dir_cs}/cram/pe/IP_WT_G1_Hho1_6337.sc.cram"
ref_cs="${dir_cs}/reference/tiny.fa"
ref_se="${dir_fx}/reference/tiny.fa"

mkdir -p "${dir_up}/logs" "${dir_up}/q"
ln -s "${in_se}" "${dir_up}/tiny_se.BAM"
ln -s "${dir_fx}/cram/se/tiny_se.cram" "${dir_up}/tiny_se.CRAM"
ln -s "${dir_fx}/cram/se/tiny_se.cram.crai" "${dir_up}/tiny_se.CRAM.crai"
ln -s "${in_cram}" "${dir_up}/IP_up.sc.CRAM"
printf 'chrI\t0\t10\t0.5\nchrI\t10\t20\t0.5\n' \
    | gzip -c > "${dir_up}/unity.bedGraph.GZ"


# Call a sourced helper library function in a fresh shell, printing its stdout
# and status.
function run_lib() {
    local lib="${1}"
    shift 1

    # shellcheck disable=SC2016  # Expand in the child shell, not this one.
    "${TEST_BASH}" -c '
        # shellcheck source=lib/bash/core/source_helpers.sh
        source "${1}/lib/bash/core/source_helpers.sh"

        source_helpers "${1}/lib/bash" "${2}"
        shift 2

        rc=0
        "$@" || rc=$?
        printf "rc=%s\n" "${rc}"
    ' _ "${ROOT_REPO}" "${lib}" "$@" 2> /dev/null
}

# Check one helper call's stdout and status: library, expected output, then the
# function and its arguments.
function check_lib() {
    local lib="${1}" out_exp="${2}"
    local out want="rc=0"
    shift 2

    if [[ -n "${out_exp}" ]]; then want="${out_exp}"$'\n'"${want}"; fi
    out="$(dir_eo="${dir_up}" nam_job=j run_lib "${lib}" "$@")"

    if [[ "${out}" == "${want}" ]]; then
        record_pass "${1} reads an upper-case suffix"
    else
        record_fail "${1} gave '${out//$'\n'/ | }'"
    fi
}

check_lib core/wrap_cmd \
    "${dir_up}/j.IP_x.stdout.txt;${dir_up}/j.IP_x.stderr.txt" \
    get_submit_logs d/IP_x.BEDGRAPH.GZ

check_lib core/format_outputs \
    d/x.n_frg.txt \
    derive_report_path d/x.bedGraph.GZ n_frg

check_lib core/check_unity \
    "" \
    check_unity "${dir_up}/unity.bedGraph.GZ" 0.99 1.01 true

check_lib workflows/process_region \
    "" \
    check_region_bdg "${dir_up}/unity.bedGraph.GZ" chrI:1-20

# A signal output takes the input's stem without its upper-case suffix.
rc=0
out="$(run_execute --mode signal --csv_fil_in "${dir_up}/tiny_se.BAM")" \
    || rc=$?
cmd="$(grep -m1 'submit_compute_signal.sh --env' <<< "${out}" || true)"

if [[ "${rc}" -eq 0 && "${cmd}" == *"/tiny_se.bedGraph.gz "* ]]; then
    record_pass "execute derives 'tiny_se.bedGraph.gz' from 'tiny_se.BAM'"
else
    record_fail "execute did not strip '.BAM' (exit ${rc})"
fi

# 'compute_signal' reads a reference for an upper-case CRAM name, without a
# note that it has no effect.
rc=0
out="$(
    PYTHONDONTWRITEBYTECODE=1 \
    "${TEST_MANAGED_PYTHON}" \
        -m protocol_chipseq_signal_norm.cli.compute_signal \
        --fil_in "${dir_up}/tiny_se.CRAM" \
        --ref_fa "${ref_se}" \
        --fil_out "${dir_up}/tiny_se.bedGraph" \
        --siz_bin 10 \
        --method unadj 2>&1
)" || rc=$?

if [[
    "${rc}" -eq 0
    && -s "${dir_up}/tiny_se.bedGraph"
    && "${out}" != *"'--ref_fa'"*
]]; then
    record_pass "compute_signal takes '--ref_fa' for 'tiny_se.CRAM'"
else
    record_fail "compute_signal did not read 'tiny_se.CRAM' (exit ${rc})"
fi

# 'execute' filter accepts '--out_ext' in any letter case and passes on the
# canonical spelling.
rc=0
out="$(
    "${TEST_BASH}" "${ROOT_REPO}/bin/execute_filter_alignments.sh" \
        --dry_run \
        --threads 1 \
        --csv_fil_in "${in_se}" \
        --dir_out "${dir_out}" \
        --dir_eo "${dir_err}" \
        --max_job 1 \
        --out_ext CRAM \
        --ref_fa "${ref_se}" 2>&1
)" || rc=$?
cmd="$(grep -m1 'submit_filter_alignments.sh --env' <<< "${out}" || true)"

if [[ "${rc}" -eq 0 && "${cmd}" == *" --out_ext cram "* ]]; then
    record_pass "execute filter passes '--out_ext CRAM' as 'cram'"
else
    record_fail "execute filter kept '--out_ext CRAM' (exit ${rc})"
fi

# Sort and convert an upper-case CRAM name: each step gets the reference, and
# outputs take the canonical suffixes.
q_up="${dir_up}/IP_up.sc.qnam.cram"
rc=0
"${TEST_BASH}" "${ROOT_REPO}/bin/submit_qsort_bam.sh" \
    --env_nam "${env_nam}" \
    --dir_scr "${ROOT_REPO}/bin" \
    --threads 1 \
    --csv_fil_in "${dir_up}/IP_up.sc.CRAM" \
    --ref_fa "${ref_cs}" \
    --dir_out "${dir_up}" \
    --dir_eo "${dir_up}/logs" \
    --nam_job test_check_formats_qsort \
    > /dev/null 2>&1 || rc=$?

if [[ "${rc}" -eq 0 && -s "${q_up}" ]]; then
    record_pass "qsort writes 'IP_up.sc.qnam.cram' from 'IP_up.sc.CRAM'"
else
    record_fail "qsort did not write '${q_up##*/}' (exit ${rc})"
fi

# A separate directory keeps the link apart from the sorted file on a file
# system that ignores letter case.
ln -s "${q_up}" "${dir_up}/q/IP_up.sc.QNAM.CRAM"
rc=0
"${TEST_BASH}" "${ROOT_REPO}/bin/submit_convert_bam_bed.sh" \
    --env_nam "${env_nam}" \
    --dir_scr "${ROOT_REPO}/bin" \
    --threads 1 \
    --csv_fil_in "${dir_up}/q/IP_up.sc.QNAM.CRAM" \
    --ref_fa "${ref_cs}" \
    --dir_out "${dir_out}" \
    --dir_eo "${dir_up}/logs" \
    --nam_job test_check_formats_convert \
    > /dev/null 2>&1 || rc=$?

if [[ "${rc}" -eq 0 && -s "${dir_out}/IP_up.sc.bed.gz" ]]; then
    record_pass "convert writes 'IP_up.sc.bed.gz' from 'IP_up.sc.QNAM.CRAM'"
else
    record_fail "convert did not write 'IP_up.sc.bed.gz' (exit ${rc})"
fi

unset row fnc arg rc_exp out_exp out want lbl rc cmd msg in_A fil_out
unset dir_up dir_cs in_cram ref_cs ref_se q_up
finish
