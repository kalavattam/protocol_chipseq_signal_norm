#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_install_atria.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="Atria installer"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"

script="${ROOT_REPO}/install/scripts/install_atria.sh"
tool="${ROOT_REPO}/tests/fixtures/install_atria/tool"
log="${TEST_DIR_LOG}/install_atria.log"

print_section "${TEST_NAME}"

for mode in fail reuse update; do
    mode_log="${TEST_DIR_LOG}/install_atria_${mode}.log"
    if \
        run_capture \
            "Atria ${mode} dry run" \
            "${mode_log}" \
            "${TEST_BASH}" "${script}" \
                --dry_run \
                --if_exists "${mode}" \
                --dir_install "${TEST_DIR_TMP}/${mode}"
    then
        record_pass "Atria ${mode} dry run exits 0"
    else
        record_fail "Atria ${mode} dry run failed"
    fi
done

latest_log="${TEST_DIR_LOG}/install_atria_latest.log"
if run_capture "Atria latest dry run" "${latest_log}" \
    "${TEST_BASH}" "${script}" --dry_run --if_exists update \
    --v_atria latest --dir_install "${TEST_DIR_TMP}/latest"
then
    assert_pattern_found \
        "${latest_log}" \
        "would resolve the latest stable" \
        "latest dry run does not query upstream"
else
    record_fail "Atria latest dry run failed"
fi

function assert_active_atria() {
    local dir_install="${1:?}"
    local nam_build="${2:?}"
    local version="${3:?}"
    local pth_current="${dir_install}/Atria/current"

    if [[ -L "${pth_current}" ]] \
        && [[ "$(readlink "${pth_current}")" == "${nam_build}" ]]
    then
        record_pass "Atria current points to ${nam_build}"
    else
        record_fail "Atria current does not point to ${nam_build}"
    fi

    if [[ -x "${pth_current}/bin/atria" ]] \
        && [[ "$("${pth_current}/bin/atria" --version)" == "${version}" ]]
    then
        record_pass "Atria current reports ${version}"
    else
        record_fail "Atria current does not report ${version}"
    fi
}


install="${TEST_DIR_TMP}/install"
fake_log="${TEST_DIR_LOG}/install_atria_fake.log"
snippet="${TEST_DIR_TMP}/install_atria_path.sh"

rm -rf "${install}"
: > "${fake_log}"

if \
    PATH="${tool}:${PATH}" \
    ATRIA_FAKE_LOG="${fake_log}" \
    ATRIA_FAKE_VERSION=v4.1.4 \
    CONDA_DEFAULT_ENV=env_protocol \
    "${TEST_BASH}" "${script}" \
        --if_exists update \
        --v_atria 4.1.4 \
        --dir_install "${install}" \
        > "${log}" 2>&1
then
    assert_active_atria "${install}" "atria-4.1.4" "v4.1.4"
else
    record_fail \
        "fake v4.1.4 installation failed; see $(print_relpath "${log}")"
fi

: > "${fake_log}"
: > "${install}/Atria/Manifest.toml"
cat > "${snippet}" << 'EOM'
export PATH="${PATH}:/legacy/julia-1.8.5/bin"
export PATH="/legacy/Atria/current/bin:${PATH}"
EOM

if \
    PATH="${tool}:${PATH}" \
    ATRIA_FAKE_LOG="${fake_log}" \
    CONDA_DEFAULT_ENV=env_protocol \
    "${TEST_BASH}" "${script}" \
        --if_exists update \
        --v_atria latest \
        --v_julia 1.8.5 \
        --path_snippet "${snippet}" \
        --dir_install "${install}" \
        > "${log}" 2>&1
then
    assert_pattern_found \
        "${fake_log}" \
        'git ls-remote --tags --refs' \
        "latest resolves tags through Git"
    assert_pattern_found \
        "${fake_log}" \
        'git checkout --detach v4.1.5' \
        "latest checks out the highest stable tag"
    assert_pattern_found \
        "${log}" \
        'Atria tag v4.1.5 (reported v4.1.5)' \
        "summary separates Atria tag and reported version"
    assert_active_atria "${install}" "atria-4.1.5" "v4.1.5"

    if [[ ! -e "${install}/Atria/Manifest.toml" ]]; then
        record_pass "Atria tag update clears the generated dependency state"
    else
        record_fail "Atria tag update retained generated dependency state"
    fi

    if [[ -d "${install}/julia-1.8.5" ]] \
        && [[ -d "${install}/Atria/atria-4.1.4" ]]
    then
        record_pass "Atria update retains inactive Julia and Atria builds"
    else
        record_fail "Atria update did not retain inactive builds"
    fi

    if \
        grep -Fq \
            '# >>> protocol_chipseq_signal_norm install_atria.sh >>>' \
            "${snippet}" \
        && grep -Fq \
            "${install}/julia-1.8.5/bin" \
            "${snippet}"
    then
        record_pass "Atria update writes the managed PATH block"
    else
        record_fail "Atria update did not write the managed PATH block"
    fi

    resolved_paths="$(
        # shellcheck disable=SC2030
        PATH="/usr/bin:/bin"

        # shellcheck disable=SC1090
        source "${snippet}"
        printf '%s\n%s\n' "$(command -v julia)" "$(command -v atria)"
    )"
    expected_paths="${install}/julia-1.8.5/bin/julia
${install}/Atria/current/bin/atria"

    if [[ "${resolved_paths}" == "${expected_paths}" ]]
    then
        record_pass "managed PATH block supersedes legacy entries"
    else
        record_fail "managed PATH block did not supersede legacy entries"
    fi
else
    record_fail \
        "fake latest installation failed; see $(print_relpath "${log}")"
fi

: > "${fake_log}"
# shellcheck disable=SC2031
if PATH="${tool}:${PATH}" ATRIA_FAKE_LOG="${fake_log}" \
    CONDA_DEFAULT_ENV=env_protocol \
    "${TEST_BASH}" "${script}" --if_exists update --v_atria latest \
        --v_julia 1.9.4 --path_snippet "${snippet}" \
        --dir_install "${install}" > "${log}" 2>&1
then
    assert_active_atria "${install}" "atria-4.1.5_julia-1.9.4" "v4.1.5"

    if [[ "$(< "${install}/Atria/atria-4.1.5_julia-1.9.4/.protocol_chipseq_julia_version")" \
        == "1.9.4" ]]
    then
        record_pass "Julia upgrade creates a matching Atria build"
    else
        record_fail "Julia upgrade did not create a matching Atria build"
    fi
else
    record_fail "fake Julia upgrade failed; see $(print_relpath "${log}")"
fi

: > "${fake_log}"
# shellcheck disable=SC2031
if PATH="${tool}:${PATH}" ATRIA_FAKE_LOG="${fake_log}" \
    CONDA_DEFAULT_ENV=env_protocol \
    "${TEST_BASH}" "${script}" --if_exists update --v_atria latest \
        --v_julia 1.8.5 --path_snippet "${snippet}" \
        --dir_install "${install}" > "${log}" 2>&1
then
    assert_pattern_absent \
        "${fake_log}" \
        '^tar ' \
        "Julia downgrade reuses the retained Julia build"
    assert_pattern_absent \
        "${fake_log}" \
        '^git checkout' \
        "Julia downgrade does not check out Atria"
    assert_active_atria "${install}" "atria-4.1.5" "v4.1.5"

    if [[ "$(< "${install}/Atria/atria-4.1.5/.protocol_chipseq_julia_version")" \
        == "1.8.5" ]]
    then
        record_pass "Julia downgrade reactivates a matching Atria build"
    else
        record_fail "Julia downgrade did not reactivate matching Atria build"
    fi
else
    record_fail "fake Julia downgrade failed; see $(print_relpath "${log}")"
fi

fail_install="${TEST_DIR_TMP}/fail_existing_install"
fail_snippet="${TEST_DIR_TMP}/fail_existing_path.sh"
rm -rf "${fail_install}"
mkdir -p "${fail_install}/Atria/.git"
printf '%s\n' "Atria sentinel" > "${fail_install}/Atria/sentinel"
printf '%s\n' "PATH sentinel" > "${fail_snippet}"
: > "${fake_log}"

# shellcheck disable=SC2031
if \
    PATH="${tool}:${PATH}" \
    ATRIA_FAKE_LOG="${fake_log}" \
    CONDA_DEFAULT_ENV=env_protocol \
    "${TEST_BASH}" "${script}" \
        --if_exists fail \
        --path_snippet "${fail_snippet}" \
        --dir_install "${fail_install}" \
        > "${log}" 2>&1
then
    record_fail "Atria fail mode unexpectedly changed an existing install"
else
    assert_pattern_found \
        "${log}" \
        "nothing was changed" \
        "Atria fail mode reports its non-mutation guarantee"

    if [[ "$(< "${fail_install}/Atria/sentinel")" == "Atria sentinel" ]] \
        && [[ ! -e "${fail_install}/julia-1.8.5" ]] \
        && [[ "$(< "${fail_snippet}")" == "PATH sentinel" ]]
    then
        record_pass "Atria fail mode leaves existing installation state unchanged"
    else
        record_fail "Atria fail mode changed existing installation state"
    fi

    if [[ ! -s "${fake_log}" ]]; then
        record_pass "Atria fail mode invokes no fake installation command"
    else
        record_fail "Atria fail mode invoked an installation command"
    fi
fi

: > "${fake_log}"
# shellcheck disable=SC2031
if \
    PATH="${tool}:${PATH}" \
    ATRIA_FAKE_LOG="${fake_log}" \
    CONDA_DEFAULT_ENV=env_protocol \
    "${TEST_BASH}" "${script}" \
        --if_exists reuse \
        --v_atria 4.1.5 \
        --v_julia 1.8.5 \
        --dir_install "${install}" \
        > "${log}" 2>&1
then
    assert_pattern_absent \
        "${fake_log}" \
        '^git checkout|^tar |^curl ' \
        "Atria reuse does not check out or rebuild matching components"
else
    record_fail "matching Atria reuse failed; see $(print_relpath "${log}")"
fi

: > "${fake_log}"
# shellcheck disable=SC2031
if \
    PATH="${tool}:${PATH}" \
    ATRIA_FAKE_LOG="${fake_log}" \
    CONDA_DEFAULT_ENV=env_protocol \
    "${TEST_BASH}" "${script}" \
        --if_exists reuse \
        --v_atria 4.1.4 \
        --v_julia 1.8.5 \
        --dir_install "${install}" \
        > "${log}" 2>&1
then
    record_fail "mismatching Atria reuse unexpectedly succeeded"
else
    assert_pattern_absent \
        "${fake_log}" \
        '^git checkout|^tar |^curl ' \
        "mismatching Atria reuse leaves components unchanged"
    assert_active_atria "${install}" "atria-4.1.5" "v4.1.5"
fi

cat > "${install}/julia-1.8.5/bin/julia" << 'EOM'
#!/usr/bin/env bash

printf 'julia version 0.0.0\n'
EOM
chmod +x "${install}/julia-1.8.5/bin/julia"
: > "${fake_log}"

# shellcheck disable=SC2031
if \
    PATH="${tool}:${PATH}" \
    ATRIA_FAKE_LOG="${fake_log}" \
    ATRIA_FAKE_TAR_FAIL=true \
    CONDA_DEFAULT_ENV=env_protocol \
    "${TEST_BASH}" "${script}" \
        --if_exists update \
        --v_atria 4.1.5 \
        --v_julia 1.8.5 \
        --dir_install "${install}" \
        > "${log}" 2>&1
then
    record_fail "failed staged Julia update unexpectedly succeeded"
else
    if \
        [[ "$("${install}/julia-1.8.5/bin/julia")" == "julia version 0.0.0" ]] \
        && ! compgen -G "${install}/julia-1.8.5.invalid.*" > /dev/null
    then
        record_pass "failed staged Julia update preserves the selected Julia"
    else
        record_fail "failed staged Julia update changed the selected Julia"
    fi

    assert_active_atria "${install}" "atria-4.1.5" "v4.1.5"
fi

: > "${fake_log}"
# shellcheck disable=SC2031
if \
    PATH="${tool}:${PATH}" \
    ATRIA_FAKE_LOG="${fake_log}" \
    CONDA_DEFAULT_ENV=env_protocol \
    "${TEST_BASH}" "${script}" \
        --if_exists update \
        --v_atria 4.1.5 \
        --v_julia 1.8.5 \
        --dir_install "${install}" \
        > "${log}" 2>&1
then
    if \
        [[ "$("${install}/julia-1.8.5/bin/julia" --version)" == "julia version 1.8.5" ]] \
        && compgen -G "${install}/julia-1.8.5.invalid.*" > /dev/null
    then
        record_pass "Atria update quarantines and replaces invalid Julia"
    else
        record_fail "Atria update did not retain and replace invalid Julia"
    fi
else
    record_fail "staged Julia update failed; see $(print_relpath "${log}")"
fi

malformed_snippet="${TEST_DIR_TMP}/install_atria_malformed_path.sh"
cat > "${malformed_snippet}" << 'EOM'
# <<< protocol_chipseq_signal_norm install_atria.sh <<<
# >>> protocol_chipseq_signal_norm install_atria.sh >>>
EOM
# shellcheck disable=SC2031
if PATH="${tool}:${PATH}" ATRIA_FAKE_LOG="${fake_log}" \
    CONDA_DEFAULT_ENV=env_protocol \
    "${TEST_BASH}" "${script}" --if_exists update --v_atria latest \
        --v_julia 1.8.5 --path_snippet "${malformed_snippet}" \
        --dir_install "${install}" > "${log}" 2>&1
then
    record_fail "malformed managed PATH block unexpectedly succeeded"
else
    assert_pattern_found \
        "${log}" \
        "has malformed managed block delimiters" \
        "malformed managed PATH block is rejected"

    if [[ "$(< "${malformed_snippet}")" == $'# <<< protocol_chipseq_signal_norm install_atria.sh <<<\n# >>> protocol_chipseq_signal_norm install_atria.sh >>>' ]]; then
        record_pass "malformed managed PATH block is unchanged"
    else
        record_fail "malformed managed PATH block was changed"
    fi
fi

bad_install="${TEST_DIR_TMP}/bad_install"
rm -rf "${bad_install}"
: > "${fake_log}"
# shellcheck disable=SC2031
if \
    PATH="${tool}:${PATH}" \
    ATRIA_FAKE_LOG="${fake_log}" \
    ATRIA_FAKE_DIR_VERSION=v4.1.5 \
    ATRIA_FAKE_VERSION=v4.1.4 \
    CONDA_DEFAULT_ENV=env_protocol \
    "${TEST_BASH}" "${script}" \
        --if_exists update \
        --v_atria 4.1.5 \
        --dir_install "${bad_install}" \
        > "${log}" 2>&1
then
    record_fail "wrong-version fake update unexpectedly succeeded"
else
    assert_pattern_found \
        "${log}" \
        "does not appear to match requested tag 'v4.1.5'" \
        "wrong-version fake update reports the mismatch"

    if [[ ! -e "${bad_install}/Atria/current" ]]; then
        record_pass "wrong-version fake update does not activate a build"
    else
        record_fail "wrong-version fake update activated a build"
    fi
fi

prune_log="${TEST_DIR_LOG}/install_atria_prune.log"
if \
    run_capture "Atria rejects retired prune option" "${prune_log}" \
    "${TEST_BASH}" "${script}" \
        --dry_run \
        --prune \
        --dir_install "${TEST_DIR_TMP}/prune"
then
    record_fail "retired Atria prune option unexpectedly succeeded"
else
    assert_pattern_found \
        "${prune_log}" \
        "unknown option/parameter passed: '--prune'" \
        "Atria rejects retired prune option"
fi

finish "$@"
