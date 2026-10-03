#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_install_preseq.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="preseq installer"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


# Build a stand-in environment prefix: the compiler, 'make', and 'htslib' the
# installer looks for, none of them functional. The installer must stop at its
# checks before it would use any of them.
function make_fake_prefix() {
    local pth="${1:?}"

    rm -rf "${pth}"
    mkdir -p "${pth}/bin" "${pth}/include/htslib" "${pth}/lib"
    printf '#!/usr/bin/env bash\nexit 1\n' > "${pth}/bin/clang++"
    printf '#!/usr/bin/env bash\nexit 1\n' > "${pth}/bin/make"
    chmod +x "${pth}/bin/clang++" "${pth}/bin/make"
    : > "${pth}/include/htslib/sam.h"
    : > "${pth}/lib/libhts.so.3"
}


function run_fake_env() {
    local label="${1:?}"
    local log="${2:?}"
    local prefix="${3:?}"
    shift 3

    run_capture \
        "${label}" \
        "${log}" \
        env \
            CONDA_DEFAULT_ENV=env_qc \
            "CONDA_PREFIX=${prefix}" \
            CXX= \
            "${TEST_BASH}" "${script}" "$@"
}


script="${ROOT_REPO}/install/scripts/install_preseq.sh"
dir_log="${TEST_DIR_LOG}/install_preseq"
mkdir -p "${dir_log}"

print_section "${TEST_NAME}"

# Help.
log_help="${dir_log}/help.log"
if \
    run_capture \
        "preseq installer help" \
        "${log_help}" \
        "${TEST_BASH}" "${script}" --help
then
    record_pass "preseq installer help exits 0"
    assert_pattern_found \
        "${log_help}" "^Usage" "preseq installer help prints Usage"
else
    record_fail "preseq installer help failed"
fi

# Dry run: the plan names the release, its checksum, the patch, and the
# configure flag, and nothing is activated or fetched.
log_dry="${dir_log}/dry_run.log"
if \
    run_capture \
        "preseq installer dry run" \
        "${log_dry}" \
        "${TEST_BASH}" "${script}" --dry_run
then
    record_pass "preseq installer dry run exits 0"
    assert_pattern_found \
        "${log_dry}" \
        "releases/download/v3\.2\.0/preseq-3\.2\.0\.tar\.gz" \
        "dry run names the pinned release"
    assert_pattern_found \
        "${log_dry}" \
        "95b81c9054e0d651de398585c7e96b807ad98f0bdc541b3e46665febbe2134d9" \
        "dry run reports the pinned checksum"
    assert_pattern_found \
        "${log_dry}" \
        "#include <cstdint>.*src/smithlab_cpp/sam_record\.hpp" \
        "dry run reports the header patch"
    assert_pattern_found \
        "${log_dry}" "--enable-hts" "dry run configures with htslib"
    assert_pattern_found \
        "${log_dry}" "env_nam=env_qc" "dry run defaults to env_qc"
else
    record_fail "preseq installer dry run failed"
fi

# Dry run with a local tarball.
fil_any="${TEST_DIR_TMP}/install_preseq_any.tar.gz"
printf 'not a tarball\n' > "${fil_any}"
log_dry_tar="${dir_log}/dry_run_tarball.log"
if \
    run_capture \
        "preseq installer dry run with tarball" \
        "${log_dry_tar}" \
        "${TEST_BASH}" "${script}" --dry_run --fil_tar "${fil_any}"
then
    assert_pattern_found \
        "${log_dry_tar}" \
        "would copy '${fil_any}' and verify its SHA-256" \
        "dry run with a local tarball copies instead of downloading"
else
    record_fail "preseq installer dry run with a tarball failed"
fi

# Argument errors.
for case_args in \
    "unknown|--not_a_real_option|unknown option/parameter passed" \
    "if_exists|--if_exists maybe|invalid '--if_exists' value" \
    "threads|--threads 0|must be an integer greater than or equal to 1" \
    "tarball|--fil_tar ${TEST_DIR_TMP}/absent.tar.gz|tarball does not exist"
do
    IFS='|' read -r nam args expect <<< "${case_args}"
    log_arg="${dir_log}/arg_${nam}.log"
    # shellcheck disable=SC2086
    if \
        run_capture \
            "preseq installer rejects ${nam}" \
            "${log_arg}" \
            "${TEST_BASH}" "${script}" --dry_run ${args}
    then
        record_fail "preseq installer accepted invalid ${nam} input"
    else
        assert_pattern_found \
            "${log_arg}" "${expect}" "preseq installer rejects invalid ${nam}"
    fi
done

# A tarball whose checksum does not match is refused before any build.
prefix="${TEST_DIR_TMP}/install_preseq_prefix"
make_fake_prefix "${prefix}"
log_sha="${dir_log}/checksum_mismatch.log"
if \
    run_fake_env \
        "preseq installer checksum mismatch" \
        "${log_sha}" \
        "${prefix}" \
        --fil_tar "${fil_any}" \
        --dir_tmp "${TEST_DIR_TMP}/sha"
then
    record_fail "preseq installer built from a tarball with a bad checksum"
else
    assert_pattern_found \
        "${log_sha}" \
        "SHA-256 mismatch" \
        "preseq installer refuses a tarball with a bad checksum"
    if [[ ! -e "${prefix}/bin/preseq" ]]; then
        record_pass "checksum mismatch installs nothing"
    else
        record_fail "checksum mismatch left a preseq in the prefix"
    fi
fi

# An environment without a compiler is reported, not worked around.
rm -f "${prefix}/bin/clang++"
log_cxx="${dir_log}/no_compiler.log"
if \
    run_fake_env \
        "preseq installer without a compiler" \
        "${log_cxx}" \
        "${prefix}" \
        --fil_tar "${fil_any}"
then
    record_fail "preseq installer proceeded without a compiler"
else
    assert_pattern_found \
        "${log_cxx}" \
        "no C++ compiler was found" \
        "preseq installer requires the environment's compiler"
fi

# An existing installation: 'fail' leaves it alone, 'reuse' accepts the pinned
# version, and 'reuse' refuses another version.
make_fake_prefix "${prefix}"
printf '#!/usr/bin/env bash\necho "Version: 3.2.0" >&2\n' \
    > "${prefix}/bin/preseq"
chmod +x "${prefix}/bin/preseq"

log_fail="${dir_log}/existing_fail.log"
if \
    run_fake_env \
        "preseq installer existing fail" \
        "${log_fail}" \
        "${prefix}" \
        --fil_tar "${fil_any}"
then
    record_fail "preseq installer replaced an existing install under 'fail'"
else
    assert_pattern_found \
        "${log_fail}" \
        "'preseq' is already installed" \
        "'--if_exists fail' stops at an existing install"
fi

log_reuse="${dir_log}/existing_reuse.log"
if \
    run_fake_env \
        "preseq installer existing reuse" \
        "${log_reuse}" \
        "${prefix}" \
        --fil_tar "${fil_any}" \
        --if_exists reuse
then
    assert_pattern_found \
        "${log_reuse}" \
        "Reusing the installed 'preseq' 3\.2\.0" \
        "'--if_exists reuse' keeps a matching install"
else
    record_fail "'--if_exists reuse' rejected a matching install"
fi

printf '#!/usr/bin/env bash\necho "Version: 3.1.2" >&2\n' \
    > "${prefix}/bin/preseq"
log_reuse_old="${dir_log}/existing_reuse_old.log"
if \
    run_fake_env \
        "preseq installer reuse of another version" \
        "${log_reuse_old}" \
        "${prefix}" \
        --fil_tar "${fil_any}" \
        --if_exists reuse
then
    record_fail "'--if_exists reuse' kept a mismatched version"
else
    assert_pattern_found \
        "${log_reuse_old}" \
        "is not version 3\.2\.0" \
        "'--if_exists reuse' refuses a mismatched version"
fi

finish "$@"
