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

unset row fnc arg rc_exp out_exp out want lbl rc cmd msg in_A fil_out
finish
