#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_download_fastqs_interface.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="download-fastqs interface"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


#  Require one command to succeed
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


#  Require one command to fail
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


#  Assert that one worker directory remains empty
function assert_dir_empty() {
    local dir="${1:-}"
    local label="${2:-directory}"

    if [[ -n "$(find "${dir}" -mindepth 1 -print -quit)" ]]; then
        record_fail "${label} contains a worker-side effect"
    else
        record_pass "${label} has no worker-side effects"
    fi
}


#  Run a command from one explicit working directory
function run_from_dir() {
    local dir="${1:-}"
    shift

    (
        cd "${dir}" || exit 1
        "$@"
    )
}


#  Assert one fixed-string occurrence count
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
env_missing="phase1b_interface_missing_env"

rm -rf "${tmp}"
mkdir -p "${dir_cap}" "${dir_out}" "${dir_sym}" "${dir_eo}"


#  Canonical execute spellings parse before deterministic validation failure
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


#  Systematic hidden hyphen spellings remain accepted but undocumented
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


#  Canonical and hidden submit bootstrap options reach positional validation
for opt in -ds --dir_scr --dir-scr; do
    log="${dir_cap}/submit_${opt#-}.log"
    log="${log//-/_}"
    expect_failure "submit ${opt} acceptance" "${log}" \
        "${TEST_BASH}" "${scr_sub}" "${opt}" "${dir_scr}"
    assert_pattern_found \
        "${log}" \
        'requires 8 positional arguments' \
        "submit ${opt} is accepted"
    assert_pattern_absent \
        "${log}" \
        'unknown parameter\|required option' \
        "submit ${opt} avoids bootstrap option errors"
done


#  Retired aliases remain rejected
retired_aliases=(-i "--in""file" -eo)
for opt in "${retired_aliases[@]}"; do
    log="${dir_cap}/retired_${opt#-}.log"
    log="${log//-/_}"
    expect_failure "retired execute alias ${opt}" "${log}" \
        "${TEST_BASH}" "${scr_exe}" "${opt}" "${fil_missing}"
    assert_pattern_found \
        "${log}" \
        "unknown option/parameter passed: '${opt}'" \
        "retired execute alias ${opt} is rejected"
done


#  Rendered help exposes canonical public aliases and exact callables only
log_exe_hlp="${dir_cap}/execute_help.log"
log_sub_hlp="${dir_cap}/submit_help.log"
expect_success "execute canonical help" "${log_exe_hlp}" \
    "${TEST_BASH}" "${scr_exe}" --help
expect_success "submit canonical help" "${log_sub_hlp}" \
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


#  Hidden help remains accepted
expect_success "execute hidden --hlp" "${dir_cap}/execute_hlp.log" \
    "${TEST_BASH}" "${scr_exe}" --hlp


#  Submit help is terminal from every required argument position
expect_success "submit help first" "${dir_cap}/submit_help_first.log" \
    "${TEST_BASH}" "${scr_sub}" --help
expect_success "submit short help first" "${dir_cap}/submit_h_first.log" \
    "${TEST_BASH}" "${scr_sub}" -h
expect_success "submit hidden help first" "${dir_cap}/submit_hlp_first.log" \
    "${TEST_BASH}" "${scr_sub}" --hlp
expect_success "submit help after junk" "${dir_cap}/submit_help_junk.log" \
    "${TEST_BASH}" "${scr_sub}" junk --help
expect_success \
    "submit help after dir_scr" \
    "${dir_cap}/submit_help_after_dir_scr.log" \
    "${TEST_BASH}" "${scr_sub}" --dir_scr "${dir_scr}" --help
expect_success \
    "submit help before dir_scr" \
    "${dir_cap}/submit_help_before_dir_scr.log" \
    "${TEST_BASH}" "${scr_sub}" --help --dir_scr "${dir_scr}"


#  A complete worker payload remains inert when followed by help
log_terminal="${dir_cap}/submit_help_full_payload.log"
expect_success "submit terminal help with worker payload" "${log_terminal}" \
    "${TEST_BASH}" "${scr_sub}" \
        --dir_scr "${dir_scr}" \
        SRR_PHASE1B "file://${tmp}/never_exists.fastq.gz" NA \
        "${dir_out}" "${dir_sym}" sample_phase1b "${dir_eo}" phase1b \
        --help

for file in "${dir_cap}"/submit_*help*.log "${dir_cap}/submit_h_first.log"; do
    assert_pattern_found \
        "${file}" \
        '^Usage$' \
        "terminal submit help renders Usage"
    assert_pattern_absent \
        "${file}" \
        'error(\|requires 8 positional\|Downloading\|Symlinking' \
        "terminal submit help bypasses validation and downloads"
done

assert_dir_empty "${dir_out}" "submit help output directory"
assert_dir_empty "${dir_sym}" "submit help symlink directory"
assert_dir_empty "${dir_eo}" "submit help worker-log directory"


#  Dry runs cover SE, PE, mixed, and duplicate-accession planning
dir_dry="${tmp}/dry_run"
mkdir -p "${dir_dry}/out" "${dir_dry}/links" "${dir_dry}/logs"

cat > "${dir_dry}/se.tsv" << 'EOM'
run_accession	custom_name	fastq_https
SRR_DRY_SE	dry_se	https://example.invalid/dry_se.fastq.gz
EOM

cat > "${dir_dry}/pe.tsv" << 'EOM'
run_accession	custom_name	fastq_https
SRR_DRY_PE	dry_pe	https://example.invalid/dry_pe_R1.fastq.gz;https://example.invalid/dry_pe_R2.fastq.gz
EOM

cat > "${dir_dry}/mixed.tsv" << 'EOM'
run_accession	custom_name	fastq_https
SRR_DRY_SE	dry_se	https://example.invalid/dry_se.fastq.gz
SRR_DRY_PE	dry_pe	https://example.invalid/dry_pe_R1.fastq.gz;https://example.invalid/dry_pe_R2.fastq.gz
EOM

cat > "${dir_dry}/duplicate.tsv" << 'EOM'
run_accession	custom_name	fastq_https
SRR_DRY_DUP	dry_dup_a	https://example.invalid/dry_dup.fastq.gz
SRR_DRY_DUP	dry_dup_b	https://example.invalid/dry_dup.fastq.gz
EOM

cat > "${dir_dry}/duplicate_custom.tsv" << 'EOM'
run_accession	custom_name	fastq_https
SRR_DRY_CUS_A	dry_duplicate	https://example.invalid/a.fastq.gz
SRR_DRY_CUS_B	dry_duplicate	https://example.invalid/b.fastq.gz
EOM

cat > "${dir_dry}/conflicting.tsv" << 'EOM'
run_accession	custom_name	fastq_https
SRR_DRY_CONFLICT	dry_conflict_a	https://example.invalid/a.fastq.gz
SRR_DRY_CONFLICT	dry_conflict_b	https://example.invalid/b.fastq.gz
EOM

for dry_case in se pe mixed duplicate; do
    log="${dir_cap}/dry_${dry_case}.log"
    expect_success "execute ${dry_case} dry run" "${log}" \
        "${TEST_BASH}" "${scr_exe}" \
            --dry_run \
            --fil_in "${dir_dry}/${dry_case}.tsv" \
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

expect_failure \
    "execute duplicate custom-name dry run" \
    "${dir_cap}/dry_duplicate_custom.log" \
    "${TEST_BASH}" "${scr_exe}" \
        --dry_run \
        --fil_in "${dir_dry}/duplicate_custom.tsv" \
        --dir_out "${dir_dry}/out" \
        --dir_sym "${dir_dry}/links" \
        --dir_eo "${dir_dry}/logs"
expect_failure \
    "execute conflicting-accession dry run" \
    "${dir_cap}/dry_conflicting.log" \
    "${TEST_BASH}" "${scr_exe}" \
        --dry_run \
        --fil_in "${dir_dry}/conflicting.tsv" \
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


#  Fake wget proves relative paths and worker-failure propagation without I/O
dir_fake="${tmp}/fake_bin"
dir_rel="${tmp}/relative"
mkdir -p \
    "${dir_fake}" \
    "${dir_rel}/submit/out" \
    "${dir_rel}/submit/links" \
    "${dir_rel}/submit/logs" \
    "${dir_rel}/execute/out" \
    "${dir_rel}/execute/links" \
    "${dir_rel}/execute/logs" \
    "${dir_rel}/failure/out" \
    "${dir_rel}/failure/links" \
    "${dir_rel}/failure/logs"

cat > "${dir_fake}/wget" << 'EOM'
#!/usr/bin/env bash
set -euo pipefail
output=""
while (( $# > 0 )); do
    case "${1}" in
        -O)
            output="${2}"
            shift 2
            ;;

        *)
            shift
            ;;
    esac
done
printf '%s\n' '@synthetic-fastq' 'ACGT' '+' 'IIII' > "${output}"
EOM
chmod +x "${dir_fake}/wget"

cat > "${dir_rel}/execute/input.tsv" << 'EOM'
run_accession	custom_name	fastq_https
SRR_REL_EXEC	relative_execute	https://example.invalid/relative.fastq.gz
EOM

expect_success \
    "submit relative-directory fake download" \
    "${dir_cap}/submit_relative.log" \
    run_from_dir "${tmp}" \
        env PATH="${dir_fake}:${PATH}" \
        "${TEST_BASH}" "${scr_sub}" \
            --dir_scr "${dir_scr}" \
            SRR_REL_SUB \
            https://example.invalid/relative.fastq.gz \
            NA \
            relative/submit/out \
            relative/submit/links \
            relative_submit \
            relative/submit/logs \
            relative_submit

assert_custom_symlink \
    "${dir_rel}/submit/links/relative_submit.fastq.gz" \
    "direct submit relative-directory symlink"
if [[ \
    "${dir_rel}/submit/links/relative_submit.fastq.gz" \
    -ef "${dir_rel}/submit/out/SRR_REL_SUB.fastq.gz" \
]]; then
    record_pass "direct submit relative symlink resolves to its output"
else
    record_fail "direct submit relative symlink target is broken"
fi

expect_success \
    "execute relative-directory fake download" \
    "${dir_cap}/execute_relative.log" \
    run_from_dir "${tmp}" \
        env PATH="${dir_fake}:${PATH}" \
        "${TEST_BASH}" "${scr_exe}" \
            --fil_in relative/execute/input.tsv \
            --dir_out relative/execute/out \
            --dir_sym relative/execute/links \
            --dir_eo relative/execute/logs

assert_custom_symlink \
    "${dir_rel}/execute/links/relative_execute.fastq.gz" \
    "execute relative-directory symlink"
if [[ \
    "${dir_rel}/execute/links/relative_execute.fastq.gz" \
    -ef "${dir_rel}/execute/out/SRR_REL_EXEC.fastq.gz" \
]]; then
    record_pass "execute relative symlink resolves to its output"
else
    record_fail "execute relative symlink target is broken"
fi

cat > "${dir_fake}/wget" << 'EOM'
#!/usr/bin/env bash
set -euo pipefail
output=""
while (( $# > 0 )); do
    case "${1}" in
        -O)
            output="${2}"
            shift 2
            ;;

        *)
            shift
            ;;
    esac
done
printf '%s\n' '@partial-fastq' 'TGCA' '+' 'IIII' > "${output}"
exit 42
EOM

cat > "${dir_rel}/failure/input.tsv" << 'EOM'
run_accession	custom_name	fastq_https
SRR_FAIL	failed_worker	https://example.invalid/failure.fastq.gz
EOM

expect_failure \
    "execute propagates fake worker failure" \
    "${dir_cap}/execute_worker_failure.log" \
    run_from_dir "${tmp}" \
        env PATH="${dir_fake}:${PATH}" \
        "${TEST_BASH}" "${scr_exe}" \
            --fil_in relative/failure/input.tsv \
            --dir_out relative/failure/out \
            --dir_sym relative/failure/links \
            --dir_eo relative/failure/logs

assert_file_nonempty \
    "${dir_rel}/failure/out/SRR_FAIL.fastq.gz" \
    "failed fake worker leaves diagnostic partial output"
if [[ ! -e "${dir_rel}/failure/links/failed_worker.fastq.gz" ]]; then
    record_pass "failed worker does not create a custom-name symlink"
else
    record_fail "failed worker unexpectedly creates a custom-name symlink"
fi

finish
