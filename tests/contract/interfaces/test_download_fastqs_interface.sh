#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_download_fastqs_interface.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.6);
# - Anthropic Claude Code (Opus 5.5).
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="download-fastqs interface"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Require one command to succeed.
function expect_success() {
    local label="${1:-command}"
    local log="${2:-}"
    shift 2

    if run_capture "${label}" "${log}" "$@"; then
        record_pass "${label} exits 0"
    else
        record_fail "${label} failed; see $(print_relpath "${log}")"
    fi
}


# Require one command to fail.
function expect_failure() {
    local label="${1:-command}"
    local log="${2:-}"
    shift 2

    if run_capture "${label}" "${log}" "$@"; then
        record_fail "${label} unexpectedly exits 0"
    else
        record_pass "${label} exits nonzero"
    fi
}


# Assert that one worker directory remains empty.
function assert_dir_empty() {
    local dir="${1:-}"
    local label="${2:-directory}"

    if [[ -n "$(find "${dir}" -mindepth 1 -print -quit)" ]]; then
        record_fail "${label} contains a worker-side effect"
    else
        record_pass "${label} has no worker-side effects"
    fi
}


# Run a command from one explicit working directory.
function run_from_dir() {
    local dir="${1:-}"
    shift

    (
        cd "${dir}" || exit 1
        "$@"
    )
}


# Assert one fixed-string occurrence count.
function assert_count_fixed() {
    local file="${1:-}"
    local pattern="${2:-}"
    local expected="${3:-0}"
    local label="${4:-pattern count}"
    local observed=0

    observed="$(grep -F -c -- "${pattern}" "${file}" || true)"
    if [[ "${observed}" == "${expected}" ]]; then
        record_pass "${label}"
    else
        record_fail \
            "${label}; expected ${expected}, observed ${observed}"
    fi
}


print_section "${TEST_NAME}"

scr_exe="${ROOT_REPO}/bin/execute_download_fastqs.sh"
scr_sub="${ROOT_REPO}/bin/submit_download_fastqs.sh"
dir_scr="${ROOT_REPO}/bin"

tmp="${TEST_DIR_TMP}/download_fastqs_interface"
dir_cap="${tmp}/capture"
dir_out="${tmp}/worker/out"
dir_sym="${tmp}/worker/links"
dir_eo="${tmp}/worker/logs"
fil_missing="${tmp}/missing.tsv"
dir_fx="${ROOT_REPO}/tests/fixtures/download_fastqs"
dir_mtd="${dir_fx}/metadata"
env_missing="phase1b_interface_missing_env"

rm -rf "${tmp}"
mkdir -p "${dir_cap}" "${dir_out}" "${dir_sym}" "${dir_eo}"

require_files_nonempty \
    "${dir_mtd}/dry_se.tsv" \
    "${dir_mtd}/worker_failure.tsv" \
    "${dir_fx}/tool/wget_ok/wget" \
    "${dir_fx}/tool/wget_fail/wget" || {
    finish
    exit $?
}


# Canonical execute spellings parse before deterministic validation failure.
log="${dir_cap}/execute_canonical.log"

expect_failure "canonical execute aliases" "${log}" \
    "${TEST_BASH}" "${scr_exe}" \
        --env_nam "${env_missing}" \
        -fi "${fil_missing}" \
        --fil_in "${fil_missing}" \
        --dir_out "${dir_out}" \
        --dir_sym "${dir_sym}" \
        -deo "${dir_eo}" \
        --dir_eo "${dir_eo}"

assert_pattern_absent \
    "${log}" \
    'unknown option/parameter' \
    "canonical execute aliases are accepted"


# Systematic hidden hyphen spellings remain accepted but undocumented.
log="${dir_cap}/execute_hidden_hyphens.log"

expect_failure "hidden execute hyphen aliases" "${log}" \
    "${TEST_BASH}" "${scr_exe}" \
        --dry-run \
        --env-nam "${env_missing}" \
        --fil-in "${fil_missing}" \
        --dir-out "${dir_out}" \
        --dir-sym "${dir_sym}" \
        --dir-symlink "${dir_sym}" \
        --nam-job phase1b_interface \
        --dir-eo "${dir_eo}"

assert_pattern_absent \
    "${log}" \
    'unknown option/parameter' \
    "hidden execute hyphen aliases are accepted"


# Canonical and hidden submit bootstrap options reach required-value checks.
for opt in -ds --dir_scr --dir-scr; do
    nam="${opt#-}"
    log="${dir_cap}/submit_${nam//-/_}.log"

    expect_failure "submit ${opt} acceptance" "${log}" \
        "${TEST_BASH}" "${scr_sub}" "${opt}" "${dir_scr}"

    assert_pattern_found \
        "${log}" \
        "'srr' is empty or unset" \
        "submit ${opt} is accepted"

    assert_pattern_absent \
        "${log}" \
        'unknown parameter\|required option' \
        "submit ${opt} avoids bootstrap option errors"
done


# The retired positional worker payload is rejected.
log="${dir_cap}/submit_positional.log"

expect_failure "submit positional payload" "${log}" \
    "${TEST_BASH}" "${scr_sub}" \
        --dir_scr "${dir_scr}" \
        SRR_PHASE1B "file://${tmp}/never_exists.fastq.gz" NA \
        "${dir_out}" "${dir_sym}" sample_phase1b "${dir_eo}" phase1b

assert_pattern_found \
    "${log}" \
    "unknown option/parameter passed: 'SRR_PHASE1B'" \
    "submit rejects the retired positional payload"


# Malformed MD5 and size values are rejected before any download.
arr_args_sub=(
    --dir_scr "${dir_scr}"
    --srr SRR_PHASE1B
    --url_1 "file://${tmp}/never_exists.fastq.gz"
    --dir_out "${dir_out}"
    --dir_sym "${dir_sym}"
    --nam_cus sample_phase1b
    --dir_eo "${dir_eo}"
    --nam_job phase1b
)

arr_rows_bad=(
    "md5_1|--md5_1|xyz|'--md5_1' must be 32 hexadecimal characters"
    "md5_2|-m2|0123|'--md5_2' must be 32 hexadecimal characters"
    "siz_1|--siz-1|abc|'--siz_1' must be a positive integer"
    "siz_2|-s2|0|'--siz_2' must be a positive integer"
)

for row in "${arr_rows_bad[@]}"; do
    IFS='|' read -r nam opt val pat <<< "${row}"
    log="${dir_cap}/submit_bad_${nam}.log"

    expect_failure \
        "submit malformed ${nam}" \
        "${log}" \
        "${TEST_BASH}" "${scr_sub}" \
            "${arr_args_sub[@]}" \
            --url_2 "file://${tmp}/never_exists_R2.fastq.gz" \
            "${opt}" \
            "${val}"

    assert_pattern_found \
        "${log}" \
        "${pat}" \
        "submit rejects a malformed ${nam} passed as '${opt}'"

    assert_pattern_absent \
        "${log}" \
        'Downloading' \
        "submit checks ${nam} before downloading"
done


# Second-mate checks are refused for a single-end entry, before any download;
# an explicit 'NA' counts as not given.
dir_se2="${tmp}/single_end_r2"
mkdir -p "${dir_se2}/out" "${dir_se2}/links" "${dir_se2}/logs"

arr_args_se2=(
    --dir_scr "${dir_scr}"
    --srr SRR_SE_R2
    --url_1 "file://${tmp}/never_exists.fastq.gz"
    --dir_out "${dir_se2}/out"
    --dir_sym "${dir_se2}/links"
    --nam_cus sample_se_r2
    --dir_eo "${dir_se2}/logs"
    --nam_job se_r2
)

arr_rows_se=(
    "md5_2|--md5_2|0123456789abcdef0123456789abcdef|refused"
    "siz_2|-s2|42|refused"
    "md5_2_na|--md5_2|NA|accepted"
)

for row in "${arr_rows_se[@]}"; do
    IFS='|' read -r nam opt val exp <<< "${row}"
    log="${dir_cap}/submit_se_${nam}.log"

    "${TEST_BASH}" "${scr_sub}" "${arr_args_se2[@]}" "${opt}" "${val}" \
        > "${log}" 2>&1 || true

    if [[ "${exp}" == "refused" ]]; then
        assert_pattern_found \
            "${log}" \
            "'--${nam}' is for paired-end entries ('--url_2')" \
            "submit refuses '${opt}' for a single-end entry"

        assert_pattern_absent \
            "${log}" \
            'Downloading' \
            "submit refuses '${opt}' before downloading"
    else
        assert_pattern_absent \
            "${log}" \
            'is for paired-end entries' \
            "submit accepts '${opt} NA' for a single-end entry"
    fi
done


# Retired aliases remain rejected.
arr_retired_aliases=(-i "--in""file" -eo)

for opt in "${arr_retired_aliases[@]}"; do
    nam="${opt#-}"
    log="${dir_cap}/retired_${nam//-/_}.log"

    expect_failure "retired execute alias ${opt}" "${log}" \
        "${TEST_BASH}" "${scr_exe}" "${opt}" "${fil_missing}"

    assert_pattern_found \
        "${log}" \
        "unknown option/parameter passed: '${opt}'" \
        "retired execute alias ${opt} is rejected"
done


# Rendered help exposes canonical public aliases and exact callables only.
log_exe_hlp="${dir_cap}/execute_help.log"
log_sub_hlp="${dir_cap}/submit_help.log"

expect_success \
    "execute canonical help" \
    "${log_exe_hlp}" \
    "${TEST_BASH}" "${scr_exe}" --help

expect_success \
    "submit canonical help" \
    "${log_sub_hlp}" \
    "${TEST_BASH}" "${scr_sub}" --help

for file in "${log_exe_hlp}" "${log_sub_hlp}"; do
    assert_pattern_found \
        "${file}" \
        '-h, --help' \
        "rendered help documents -h, --help"

    assert_pattern_absent \
        "${file}" \
        '--hlp' \
        "rendered help omits hidden --hlp"
done

for opt in \
    --dry-run --env-nam --fil-in --dir-out --dir-sym --dir-symlink \
    --nam-job --dir-eo --dir-scr
do
    assert_pattern_absent \
        "${log_exe_hlp}" \
        "${opt}" \
        "execute help omits hidden ${opt}"
done

for opt in --hlp --dir-scr; do
    assert_pattern_absent \
        "${log_sub_hlp}" \
        "${opt}" \
        "submit help omits hidden ${opt}"
done

for dep in \
    awk bash basename cat conda cut dirname grep ln parallel rm sbatch wget
do
    assert_pattern_found \
        "${log_exe_hlp}" \
        "[+-] ${dep}" \
        "execute help documents callable ${dep}"
done

for dep in bash basename cat dirname ln wget; do
    assert_pattern_found \
        "${log_sub_hlp}" \
        "[+-] ${dep}" \
        "submit help documents callable ${dep}"
done

assert_pattern_absent \
    "${log_sub_hlp}" \
    '--dry_run' \
    "submit help does not imply a nonexistent dry-run mode"

assert_pattern_found \
    "${log_sub_hlp}" \
    "Maintained entrypoint directory" \
    "submit help identifies the maintained bootstrap directory"


# Hidden help remains accepted.
expect_success \
    "execute hidden --hlp" \
    "${dir_cap}/execute_hlp.log" \
    "${TEST_BASH}" "${scr_exe}" --hlp


# Submit help is terminal from every required argument position.
expect_success \
    "submit help first" \
    "${dir_cap}/submit_help_first.log" \
    "${TEST_BASH}" "${scr_sub}" --help

expect_success \
    "submit short help first" \
    "${dir_cap}/submit_h_first.log" \
    "${TEST_BASH}" "${scr_sub}" -h

expect_success \
    "submit hidden help first" \
    "${dir_cap}/submit_hlp_first.log" \
    "${TEST_BASH}" "${scr_sub}" --hlp

expect_success \
    "submit help after junk" \
    "${dir_cap}/submit_help_junk.log" \
    "${TEST_BASH}" "${scr_sub}" junk --help

expect_success \
    "submit help after dir_scr" \
    "${dir_cap}/submit_help_after_dir_scr.log" \
    "${TEST_BASH}" "${scr_sub}" --dir_scr "${dir_scr}" --help

expect_success \
    "submit help before dir_scr" \
    "${dir_cap}/submit_help_before_dir_scr.log" \
    "${TEST_BASH}" "${scr_sub}" --help --dir_scr "${dir_scr}"


# A complete worker payload remains inert when followed by help.
log_terminal="${dir_cap}/submit_help_full_payload.log"

expect_success \
    "submit terminal help with worker payload" \
    "${log_terminal}" \
    "${TEST_BASH}" "${scr_sub}" \
        --dir_scr "${dir_scr}" \
        --srr SRR_PHASE1B \
        --url_1 "file://${tmp}/never_exists.fastq.gz" \
        --dir_out "${dir_out}" \
        --dir_sym "${dir_sym}" \
        --nam_cus sample_phase1b \
        --dir_eo "${dir_eo}" \
        --nam_job phase1b \
        --help

for file in "${dir_cap}"/submit_*help*.log "${dir_cap}/submit_h_first.log"; do
    assert_pattern_found \
        "${file}" \
        '^Usage$' \
        "terminal submit help renders Usage"

    assert_pattern_absent \
        "${file}" \
        'error(\|is empty or unset\|Downloading\|Symlinking' \
        "terminal submit help bypasses validation and downloads"
done

assert_dir_empty "${dir_out}" "submit help output directory"
assert_dir_empty "${dir_sym}" "submit help symlink directory"
assert_dir_empty "${dir_eo}" "submit help worker-log directory"


# Dry runs cover SE, PE, mixed, and duplicate-accession planning.
dir_dry="${tmp}/dry_run"
mkdir -p "${dir_dry}/out" "${dir_dry}/links" "${dir_dry}/logs"

for dry_case in se pe mixed duplicate; do
    log="${dir_cap}/dry_${dry_case}.log"

    expect_success \
        "execute ${dry_case} dry run" \
        "${log}" \
        "${TEST_BASH}" "${scr_exe}" \
            --dry_run \
            --fil_in "${dir_mtd}/dry_${dry_case}.tsv" \
            --dir_out "${dir_dry}/out" \
            --dir_sym "${dir_dry}/links" \
            --dir_eo "${dir_dry}/logs"
done

assert_count_fixed \
    "${dir_cap}/dry_se.log" \
    "/submit_download_fastqs.sh" \
    1 \
    "SE dry run plans one worker"

assert_count_fixed \
    "${dir_cap}/dry_pe.log" \
    "/submit_download_fastqs.sh" \
    1 \
    "PE dry run plans one worker"

assert_pattern_absent \
    "${dir_cap}/dry_se.log" \
    '--url_2\|--md5_\|--siz_' \
    "SE dry run passes no second-URL, MD5, or size options"

assert_pattern_found \
    "${dir_cap}/dry_pe.log" \
    '--url_2' \
    "PE dry run passes '--url_2'"

assert_count_fixed \
    "${dir_cap}/dry_mixed.log" \
    "/submit_download_fastqs.sh" \
    2 \
    "mixed dry run plans two workers"

assert_count_fixed \
    "${dir_cap}/dry_duplicate.log" \
    "/submit_download_fastqs.sh" \
    1 \
    "duplicate-accession dry run plans one unique worker"

# A Slurm dry run prints one 'sbatch' call that wraps GNU Parallel; the
# 'sbatch' stub passes the tool check here and fails if the dry run submits.
dir_tol="${ROOT_REPO}/tests/fixtures/slurm_dry_run/tool"
log="${dir_cap}/dry_slurm.log"

expect_success \
    "execute Slurm dry run" \
    "${log}" \
    env PATH="${dir_tol}:${PATH}" \
    "${TEST_BASH}" "${scr_exe}" \
        --dry_run \
        --slurm \
        --threads 2 \
        --fil_in "${dir_mtd}/dry_mixed.tsv" \
        --dir_out "${dir_dry}/out" \
        --dir_sym "${dir_dry}/links" \
        --dir_eo "${dir_dry}/logs"

assert_count_fixed \
    "${log}" \
    "sbatch --job-name=" \
    1 \
    "Slurm dry run prints one 'sbatch' call"

assert_count_fixed \
    "${log}" \
    "--wrap=parallel\\ --jobs\\ 2\\ \\<\\ " \
    1 \
    "Slurm dry run wraps 'parallel --jobs 2' around the config"

expect_failure \
    "execute duplicate custom-name dry run" \
    "${dir_cap}/dry_duplicate_custom.log" \
    "${TEST_BASH}" "${scr_exe}" \
        --dry_run \
        --fil_in "${dir_mtd}/dry_duplicate_custom.tsv" \
        --dir_out "${dir_dry}/out" \
        --dir_sym "${dir_dry}/links" \
        --dir_eo "${dir_dry}/logs"

expect_failure \
    "execute conflicting-accession dry run" \
    "${dir_cap}/dry_conflicting.log" \
    "${TEST_BASH}" "${scr_exe}" \
        --dry_run \
        --fil_in "${dir_mtd}/dry_conflicting.tsv" \
        --dir_out "${dir_dry}/out" \
        --dir_sym "${dir_dry}/links" \
        --dir_eo "${dir_dry}/logs"

assert_pattern_found \
    "${dir_cap}/dry_duplicate_custom.log" \
    "duplicate custom_name" \
    "dry run rejects a duplicate custom name"

assert_pattern_found \
    "${dir_cap}/dry_conflicting.log" \
    "conflicting.*FASTQ URL" \
    "dry run rejects conflicting accession URLs"


# Fake wget proves relative paths and worker-failure propagation without I/O.
# The metadata directory is linked under the relative root, so that inputs are
# also named by relative paths and still read in place.
dir_fake="${dir_fx}/tool/wget_ok"
dir_rel="${tmp}/relative"
mkdir -p \
    "${dir_rel}/submit/out" \
    "${dir_rel}/submit/links" \
    "${dir_rel}/submit/logs" \
    "${dir_rel}/execute/out" \
    "${dir_rel}/execute/links" \
    "${dir_rel}/execute/logs" \
    "${dir_rel}/failure/out" \
    "${dir_rel}/failure/links" \
    "${dir_rel}/failure/logs"

ln -s "${dir_mtd}" "${dir_rel}/metadata"


expect_success \
    "submit relative-directory fake download" \
    "${dir_cap}/submit_relative.log" \
    run_from_dir "${tmp}" \
        env PATH="${dir_fake}:${PATH}" \
        "${TEST_BASH}" "${scr_sub}" \
            --dir_scr "${dir_scr}" \
            --srr SRR_REL_SUB \
            --url_1 https://example.invalid/relative.fastq.gz \
            --dir_out relative/submit/out \
            --dir_sym relative/submit/links \
            --nam_cus relative_submit \
            --dir_eo relative/submit/logs \
            --nam_job relative_submit

assert_custom_symlink \
    "${dir_rel}/submit/links/relative_submit.fastq.gz" \
    "direct submit relative-directory symlink"

if [[
    "${dir_rel}/submit/links/relative_submit.fastq.gz" \
    -ef "${dir_rel}/submit/out/SRR_REL_SUB.fastq.gz"
]]; then
    record_pass "direct submit relative symlink resolves to its output"
else
    record_fail "direct submit relative symlink target is broken"
fi

# A URL without a scheme, as in ENA's 'fastq_ftp' values, is fetched over
# 'https://'.
mkdir -p \
    "${dir_rel}/scheme/out" \
    "${dir_rel}/scheme/links" \
    "${dir_rel}/scheme/logs"

expect_success \
    "submit fake download of a URL without a scheme" \
    "${dir_cap}/submit_no_scheme.log" \
    env PATH="${dir_fake}:${PATH}" \
    "${TEST_BASH}" "${scr_sub}" \
        --dir_scr "${dir_scr}" \
        --srr SRR_SCHEME \
        --url_1 ftp.example.invalid/vol1/fastq/scheme.fastq.gz \
        --dir_out "${dir_rel}/scheme/out" \
        --dir_sym "${dir_rel}/scheme/links" \
        --nam_cus scheme \
        --dir_eo "${dir_rel}/scheme/logs" \
        --nam_job scheme

assert_file_exact_line \
    "${dir_rel}/scheme/out/wget_url.txt" \
    "https://ftp.example.invalid/vol1/fastq/scheme.fastq.gz" \
    "submit adds 'https://' to a URL without a scheme"

expect_success \
    "execute relative-directory fake download" \
    "${dir_cap}/execute_relative.log" \
    run_from_dir "${tmp}" \
        env PATH="${dir_fake}:${PATH}" \
        "${TEST_BASH}" "${scr_exe}" \
            --fil_in relative/metadata/relative_execute.tsv \
            --dir_out relative/execute/out \
            --dir_sym relative/execute/links \
            --dir_eo relative/execute/logs

assert_custom_symlink \
    "${dir_rel}/execute/links/relative_execute.fastq.gz" \
    "execute relative-directory symlink"

if [[
    "${dir_rel}/execute/links/relative_execute.fastq.gz" \
    -ef "${dir_rel}/execute/out/SRR_REL_EXEC.fastq.gz"
]]; then
    record_pass "execute relative symlink resolves to its output"
else
    record_fail "execute relative symlink target is broken"
fi

dir_fake="${dir_fx}/tool/wget_fail"

expect_failure \
    "execute propagates fake worker failure" \
    "${dir_cap}/execute_worker_failure.log" \
    run_from_dir "${tmp}" \
        env PATH="${dir_fake}:${PATH}" \
        "${TEST_BASH}" "${scr_exe}" \
            --fil_in relative/metadata/worker_failure.tsv \
            --dir_out relative/failure/out \
            --dir_sym relative/failure/links \
            --dir_eo relative/failure/logs

if [[
    ! -e "${dir_rel}/failure/out/SRR_FAIL.fastq.gz"
        && ! -e "${dir_rel}/failure/out/SRR_FAIL.fastq.gz.part"
]]; then
    record_pass "failed worker leaves no FASTQ or partial file"
else
    record_fail "failed worker leaves a FASTQ or partial file"
fi

if [[ ! -e "${dir_rel}/failure/links/failed_worker.fastq.gz" ]]; then
    record_pass "failed worker does not create a custom-name symlink"
else
    record_fail "failed worker unexpectedly creates a custom-name symlink"
fi


# 'link_fastqs' refuses a missing accession file and creates no link.
mkdir -p "${dir_rel}/link/out" "${dir_rel}/link/links"
log="${dir_cap}/link_fastqs_missing.log"

# shellcheck disable=SC2016  # Expand in the child shell, not this one.
expect_failure \
    "link_fastqs with a missing accession file" \
    "${log}" \
    "${TEST_BASH}" -c '
        source "${1}/lib/bash/workflows/download_fastqs.sh"
        link_fastqs "${2}/out" SRR_LINK "${2}/links" missing_link se
    ' _ "${ROOT_REPO}" "${dir_rel}/link"

assert_pattern_found \
    "${log}" \
    "downloaded FASTQ is missing or empty" \
    "link_fastqs reports a missing accession file"

if [[ ! -L "${dir_rel}/link/links/missing_link.fastq.gz" ]]; then
    record_pass "link_fastqs creates no link to a missing file"
else
    record_fail "link_fastqs created a link to a missing file"
fi


finish
