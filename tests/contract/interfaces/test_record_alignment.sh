#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_record_alignment.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="alignment-record helpers"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Define output paths. 'print_record_commit' names the checkout holding the
# libraries, so a copy of the core libraries goes into a throwaway Git
# repository, and another into a directory outside any repository.
tmp="${TEST_DIR_TMP}/record_alignment"
dir_git="${tmp}/repo"
dir_nog="${tmp}/nogit"


print_section "${TEST_NAME}"

rm -rf "${tmp}"
mkdir -p "${dir_git}/lib/bash" "${dir_nog}/lib/bash"

cp -R "${ROOT_REPO}/lib/bash/core" "${dir_git}/lib/bash/"
cp -R "${ROOT_REPO}/lib/bash/core" "${dir_nog}/lib/bash/"

git -C "${dir_git}" init -q
git -C "${dir_git}" add lib
git -C "${dir_git}" \
    -c user.name=test -c user.email=test@example.com \
    commit -q -m "fixture"
sha="$(git -C "${dir_git}" rev-parse HEAD)"


# Call one helper from a copy of the core libraries in a fresh shell: 'commit'
# prints the commit, and 'run' prints the run lines for the arguments that
# follow. Git stops looking for a repository at 'tmp', so the copy outside the
# throwaway repository is outside every repository.
function run_rec() {
    local dir="${1}"
    shift 1

    # shellcheck disable=SC2016  # Expand in the child shell, not this one.
    GIT_CEILING_DIRECTORIES="${tmp}" "${TEST_BASH}" -c '
        # shellcheck source=lib/bash/core/source_helpers.sh
        source "${1}/lib/bash/core/source_helpers.sh"

        source_helpers "${1}/lib/bash" core/record_alignment
        if [[ "${2}" == "commit" ]]; then
            print_record_commit
        else
            shift 2
            arr=()
            add_record_run_lines arr "$@" || exit $?
            printf "%s\n" "${arr[@]}"
        fi
    ' _ "${dir}" "$@"
}


# Print the run lines for the given arguments, one per line.
function print_run() {
    run_rec "${dir_git}" run "$@"
}


# The commit: the full hash when the checkout is clean, with '-dirty' only when
# a tracked file changes, and 'unknown' outside a repository.
got="$(run_rec "${dir_git}" commit)"
if [[ "${got}" == "${sha}" ]]; then
    record_pass "a clean checkout gives its full commit hash"
else
    record_fail "a clean checkout gave '${got}', not '${sha}'"
fi

touch "${dir_git}/untracked.txt"
got="$(run_rec "${dir_git}" commit)"
if [[ "${got}" == "${sha}" ]]; then
    record_pass "an untracked file does not mark the commit dirty"
else
    record_fail "an untracked file gave '${got}'"
fi

printf '\n' >> "${dir_git}/lib/bash/core/record_alignment.sh"
got="$(run_rec "${dir_git}" commit)"
if [[ "${got}" == "${sha}-dirty" ]]; then
    record_pass "a changed tracked file marks the commit dirty"
else
    record_fail "a changed tracked file gave '${got}'"
fi

got="$(run_rec "${dir_nog}" commit)"
if [[ "${got}" == "unknown" ]]; then
    record_pass "outside a Git checkout, the commit is 'unknown'"
else
    record_fail "outside a Git checkout, the commit was '${got}'"
fi


# The run lines: threads, start time, host, call path, and environment, in that
# order; then the paths in the order given, relative ones made absolute and
# empty ones skipped.
got="$(
    cd "${tmp}" \
        && unset SLURM_JOB_ID SLURM_ARRAY_TASK_ID \
        && unset PIPELINE_CALLS \
        && CONDA_DEFAULT_ENV=env_x print_run fnc_x 3 \
            fil_in in.bam ref_fa "" fil_out /abs/out.bam
)"
keys="$(sed -E $'s/^@CO\tpipeline-run: ([^=]+)=.*$/\\1/' <<< "${got}")"
if [[
    "${keys}" == $'--threads\nstart\nhost\ncalls\nenv\n--fil_in\n--fil_out'
    && "${got}" == *$'pipeline-run: calls=fnc_x\n'*
    && "${got}" == *$'pipeline-run: --threads=3\n'*
    && "${got}" == *$'pipeline-run: env=env_x\n'*
    && "${got}" == *$'pipeline-run: --fil_in='"${tmp}/in.bam"$'\n'*
    && "${got}" == *$'pipeline-run: --fil_out=/abs/out.bam'
]]; then
    record_pass "run lines come in order, with relative paths made absolute"
else
    record_fail "run lines are wrong:"
    printf '%s\n' "${got}" >&2
fi

got="$(
    unset CONDA_DEFAULT_ENV \
        && PIPELINE_CALLS='a.sh > b.sh' \
            SLURM_JOB_ID=123 SLURM_ARRAY_TASK_ID=4 print_run fnc_x 1
)"
if [[
    "${got}" == *$'pipeline-run: env=none\n'*
    && "${got}" == *$'pipeline-run: slurm_job=123\n'*
    && "${got}" == *$'pipeline-run: slurm_task=4'
    && "${got}" == *$'pipeline-run: calls=a.sh > b.sh > fnc_x\n'*
]]; then
    record_pass "run lines give the call path, Slurm job and task, and no env"
else
    record_fail "run lines miss the call path, Slurm job, task, or no env:"
    printf '%s\n' "${got}" >&2
fi

# In an array task, the job is the array's ID, not the task's own job ID.
got="$(
    SLURM_JOB_ID=124 SLURM_ARRAY_JOB_ID=123 SLURM_ARRAY_TASK_ID=1 \
        print_run fnc_x 1
)"
if [[
    "${got}" == *$'pipeline-run: slurm_job=123\n'*
    && "${got}" == *$'pipeline-run: slurm_task=1'
    && "${got}" != *'slurm_job=124'*
]]; then
    record_pass "an array task gives the array's job ID and its task"
else
    record_fail "an array task's job line is not the array's ID:"
    grep 'slurm_' <<< "${got}" >&2 || true
fi

rc=0
print_run fnc_x 1 fil_in > /dev/null 2>&1 || rc=$?
if (( rc == 1 )); then
    record_pass "a key without a path is refused"
else
    record_fail "a key without a path returned ${rc}, not 1"
fi


finish
