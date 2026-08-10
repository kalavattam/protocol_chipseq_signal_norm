#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: run_tests.sh
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


if [[ -z "${BASH_VERSION:-}" ]]; then
    echo "error(run_tests.sh): Bash >= 4.4 is required." >&2
    exit 1
elif ((
    BASH_VERSINFO[0] < 4 ||
    (BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] < 4)
)); then
    echo \
        "error(run_tests.sh): Bash >= 4.4 is required; found" \
        "'${BASH_VERSION}'." >&2
    exit 1
fi

set -u


function die() {
    echo "error(run_tests.sh):" "$@" >&2
    exit 1
}


function usage() {
    cat << EOM
Usage: bash tests/run_tests.sh [--list] [group ...]

Groups:
  unit                  Python unit tests
  contract              Repository and interface contracts
  integration-local     Fixture-backed local integrations
  integration-parallel  GNU Parallel integrations (RUN_PARALLEL=1)
  integration-slurm     Non-wet Slurm integration (RUN_SLURM=1)
  all-safe              unit, contract, local, and parallel groups

With no group, the runner uses 'all-safe'. The separately confirmed two-job wet validation is never discovered by this runner.

--list prints the deduplicated selection without creating artifacts, generating fixtures, or executing tests.
EOM
}


function gate_enabled() {
    local name="${1}"
    local value="${!name:-}"

    value="$(normalize_bool "${value:-false}" "${name}")" ||
        die "invalid Boolean gate ${name}='${!name}'"
    [[ "${value}" == "true" ]]
}


function ensure_fixture() {
    local feature="${1}"
    local sentinel="${2}"
    local environment="${3:-}"
    local recipe="${repo_root}/tests/fixtures/${feature}/make.sh"
    local -a arr_command=( "${managed_bash}" "${recipe}" )

    [[ -s "${sentinel}" ]] && return 0

    [[ -r "${recipe}" ]] || die "missing fixture recipe: ${recipe}"

    if [[ -n "${environment}" && "${environment}" != "env_protocol" ]]; then
        die "unsupported fixture environment: ${environment}"
    fi

    echo "note(run_tests.sh): generating ${feature} fixtures." >&2

    "${arr_command[@]}" || die "fixture generation failed: ${feature}"

    [[ -s "${sentinel}" ]] || die "fixture sentinel is missing: ${sentinel}"
}


function ensure_integration_fixtures() {
    ensure_fixture align_fastqs \
        "${repo_root}/tests/fixtures/align_fastqs/reference/tiny.fa" \
        env_protocol

    ensure_fixture calculate_scaling_factor \
        "${repo_root}/tests/fixtures/calculate_scaling_factor/reference/tiny.fa" \
        env_protocol

    ensure_fixture compute_signal \
        "${repo_root}/tests/fixtures/compute_signal/bam/se/tiny_se.bam" \
        env_protocol

    if gate_enabled RUN_DOWNLOAD; then
        ensure_fixture download_fastqs \
            "${repo_root}/tests/fixtures/download_fastqs/metadata/local_se.template.tsv"
    fi

    ensure_fixture filter_alignments \
        "${repo_root}/tests/fixtures/filter_alignments/reference/filter_sc_sp.fa" \
        env_protocol

    ensure_fixture trim_fastqs \
        "${repo_root}/tests/fixtures/trim_fastqs/fastq/se/tiny_se.fastq.gz"
}


# Unit-scoped fixtures are prepared separately from the integration set,
# because unit tests run without the local, Parallel, and Slurm gates that
# decide whether workflow fixtures are needed at all.
function ensure_unit_fixtures() {
    # Header fixtures are literal text, so this recipe needs no environment.
    ensure_fixture ai_attribution \
        "${repo_root}/tests/fixtures/ai_attribution/source/multi_vendor.sh"

    # JSON form fixtures are literal text for the same reason.
    ensure_fixture json_source_form \
        "${repo_root}/tests/fixtures/json_source_form/source/canonical.json"
}


function add_shell_tests() {
    local root="${1}"
    local file=""

    while IFS= read -r file; do
        arr_shell_tests+=( "${file}" )
    done < <(find "${root}" -type f -name 'test_*.sh' -print | LC_ALL=C sort)
}


function add_group() {
    local group="${1}"

    case "${group}" in
        -h|--hlp|--help)
            usage
            exit 0
            ;;

        unit)
            run_unit=1
            ;;

        contract)
            add_shell_tests "${repo_root}/tests/contract"
            ;;

        integration-local)
            selected_local=1
            add_shell_tests "${repo_root}/tests/integration/local"
            ;;

        integration-parallel)
            selected_parallel=1
            add_shell_tests "${repo_root}/tests/integration/parallel"
            ;;

        integration-slurm)
            selected_slurm=1
            arr_shell_tests+=(
                "${repo_root}/tests/integration/slurm/align_fastqs/test_execute_align_fastqs.sh"
            )
            ;;

        all-safe)
            add_group unit
            add_group contract
            add_group integration-local
            add_group integration-parallel
            ;;

        *)
            die "unknown group: ${group}"
            ;;
    esac
}


function validate_group() {
    local group="${1}"

    case "${group}" in
        -h|--hlp|--help|unit|contract|integration-local|integration-parallel|\
        integration-slurm|all-safe) return 0 ;;
        *) die "unknown group: ${group}" ;;
    esac
}


repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=lib/bash/core/check_args.sh
source "${repo_root}/lib/bash/core/check_args.sh" ||
    die "failed to source the canonical Boolean normalizer"

list_only=0
if [[ "${1:-}" == "--list" ]]; then
    list_only=1
    shift
fi
run_unit=0
need_fixtures=0
selected_local=0
selected_parallel=0
selected_slurm=0

declare -a arr_shell_tests=()
declare -a arr_groups=( "${@:-all-safe}" )

for group in "${arr_groups[@]}"; do
    validate_group "${group}"
done

for group in "${arr_groups[@]}"; do
    add_group "${group}"
done

declare -A arr_seen_tests=()
declare -a arr_unique_tests=()
for test_script in "${arr_shell_tests[@]}"; do
    test_script="$(
        cd "$(dirname "${test_script}")" > /dev/null 2>&1 &&
            printf '%s/%s\n' "$(pwd -P)" "$(basename "${test_script}")"
    )"
    if [[ -z "${arr_seen_tests[${test_script}]:-}" ]]; then
        arr_seen_tests["${test_script}"]=1
        arr_unique_tests+=( "${test_script}" )
    fi
done
arr_shell_tests=( "${arr_unique_tests[@]}" )

if (( list_only == 1 )); then
    (( run_unit == 1 )) && printf '%s\n' unit
    for test_script in "${arr_shell_tests[@]}"; do
        printf '%s\n' "${test_script#"${repo_root}/"}"
    done
    exit 0
fi

if [[ -v TEST_ARTIFACT_ROOT && -z "${TEST_ARTIFACT_ROOT}" ]]; then
    die "TEST_ARTIFACT_ROOT must not be empty"
fi

artifact_root="${TEST_ARTIFACT_ROOT-${repo_root}/artifacts/tests}"
if [[ "${artifact_root}" != /* || "${artifact_root}" == "/" ]]; then
    die "TEST_ARTIFACT_ROOT must be an absolute non-root path"
fi

artifact_parent="$(
    cd "$(dirname "${artifact_root}")" > /dev/null 2>&1 && pwd -P
)" || die "TEST_ARTIFACT_ROOT parent does not exist"

artifact_root="${artifact_parent}/$(basename "${artifact_root}")"
if [[ -L "${artifact_root}" ]]; then
    die "TEST_ARTIFACT_ROOT must not be a symlink"
fi

case "${artifact_root}/" in
    "${repo_root}/"*)
        if [[ "${artifact_root}" != "${repo_root}/artifacts/tests" ]]; then
            die "TEST_ARTIFACT_ROOT must be outside the repository"
        fi
        ;;
esac
export TEST_ARTIFACT_ROOT="${artifact_root}"

log_root="${artifact_root}/logs"
mkdir -p "${log_root}"

# shellcheck disable=SC2016
if [[ "${CONDA_DEFAULT_ENV:-}" == "env_protocol" && \
    -n "${CONDA_PREFIX:-}" ]]
then
    target_prefix="${CONDA_PREFIX}"
else
    conda_executable="${CONDA_EXE:-}"
    if [[ -z "${conda_executable}" || ! -x "${conda_executable}" ]]; then
        conda_executable="$(command -v conda || true)"
    fi

    [[ -n "${conda_executable}" ]] \
        || die "Conda is required to resolve env_protocol"
    target_prefix="$(
        ${conda_executable} run -n env_protocol /bin/sh -c \
            'printf "%s\n" "$CONDA_PREFIX"'
        )" || die "failed to resolve env_protocol prefix"
fi

managed_bash="${target_prefix}/bin/bash"
managed_python="${target_prefix}/bin/python"

[[ -x "${managed_bash}" ]] \
    || die "managed Bash is unavailable: ${managed_bash}"
[[ -x "${managed_python}" ]] \
    || die "managed Python is unavailable: ${managed_python}"

export CONDA_DEFAULT_ENV=env_protocol
export CONDA_PREFIX="${target_prefix}"
export PATH="${target_prefix}/bin:${PATH}"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONPYCACHEPREFIX="${artifact_root}/pycache"
export PYTEST_ADDOPTS="${PYTEST_ADDOPTS:+${PYTEST_ADDOPTS} }-o cache_dir=${artifact_root}/pytest_cache"
export TEST_MANAGED_PREFIX="${target_prefix}"
export TEST_MANAGED_BASH="${managed_bash}"
export TEST_MANAGED_PYTHON="${managed_python}"

if (( selected_local == 1 )); then
    need_fixtures=1
fi
if (( selected_parallel == 1 )) && gate_enabled RUN_PARALLEL; then
    need_fixtures=1
fi
if (( selected_slurm == 1 )) && gate_enabled RUN_SLURM; then
    need_fixtures=1
fi

if (( need_fixtures == 1 )); then
    ensure_integration_fixtures
fi

passed=0
failed=0

if (( run_unit == 1 )); then
    ensure_unit_fixtures

    arr_unit_python=( "${managed_python}" )

    "${managed_python}" -c 'import pytest' > /dev/null 2>&1 || \
        die "pytest is unavailable in env_protocol"

    echo "==== unit ===="

    if "${arr_unit_python[@]}" -m pytest -q \
        "${repo_root}/tests/unit" \
        > "${log_root}/unit.log" 2>&1
    then
        cat "${log_root}/unit.log"
        passed=$(( passed + 1 ))
    else
        cat "${log_root}/unit.log"
        failed=$(( failed + 1 ))
    fi
fi

for test_script in "${arr_shell_tests[@]}"; do
    relative="${test_script#"${repo_root}/tests/"}"
    log="${log_root}/${relative%.sh}.log"

    mkdir -p "$(dirname "${log}")"

    echo "==== ${relative} ===="

    if "${managed_bash}" "${test_script}" > "${log}" 2>&1; then
        cat "${log}"
        passed=$(( passed + 1 ))
    else
        cat "${log}"
        failed=$(( failed + 1 ))
    fi
done

printf 'Groups passed: %d\n' "${passed}"
printf 'Groups failed: %d\n' "${failed}"
(( failed == 0 ))
