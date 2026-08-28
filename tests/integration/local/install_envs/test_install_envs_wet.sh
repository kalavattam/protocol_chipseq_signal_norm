#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_install_envs_wet.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="install envs wet"

# Source shared test helpers.
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"

script="${ROOT_REPO}/install/scripts/install_envs.sh"
env_nam="env_siqchip"
log_create="${TEST_DIR_LOG}/install_envs_wet_create.log"
log_update="${TEST_DIR_LOG}/install_envs_wet_update.log"
dir_wet=""
dir_home=""
dir_envs=""
dir_pkgs=""
pkg_mgr=""


function cleanup_wet_root() {
    if [[ -n "${dir_wet}" && -d "${dir_wet}" ]]; then
        rm -rf "${dir_wet}"
    fi
}


function run_isolated_install() {
    local log="${1:?}"
    shift

    run_capture \
        "install_envs isolated ${env_nam}" \
        "${log}" \
        env \
            CONDA_DEFAULT_ENV= \
            CONDA_PREFIX= \
            "HOME=${dir_home}" \
            "XDG_CACHE_HOME=${dir_home}/.cache" \
            "CONDA_ENVS_PATH=${dir_envs}" \
            "CONDA_PKGS_DIRS=${dir_pkgs}" \
            "${TEST_BASH}" "${script}" "$@"
}


function run_isolated_tool() {
    env \
        CONDA_DEFAULT_ENV= \
        CONDA_PREFIX= \
        "HOME=${dir_home}" \
        "XDG_CACHE_HOME=${dir_home}/.cache" \
        "CONDA_ENVS_PATH=${dir_envs}" \
        "CONDA_PKGS_DIRS=${dir_pkgs}" \
        "${pkg_mgr}" run -n "${env_nam}" "$@"
}


print_section "${TEST_NAME}"

if ! gate_value="$(normalize_test_gate RUN_INSTALL_ENVS)"; then
    record_fail "RUN_INSTALL_ENVS has an invalid Boolean value"
    finish "$@"
    exit $?
fi

if [[ "${gate_value}" != "true" ]]; then
    record_skip \
        "install_envs wet check disabled; set RUN_INSTALL_ENVS=1 to enable"
    finish "$@"
    exit $?
fi

if check_cmd_exists mamba; then
    pkg_mgr="mamba"
elif check_cmd_exists conda; then
    pkg_mgr="conda"
else
    record_skip "install_envs wet check requires mamba or conda on PATH"
    finish "$@"
    exit $?
fi

dir_wet="$(mktemp -d "${TEST_DIR_TMP}/install_envs_wet.XXXXXX")"
dir_home="${dir_wet}/home"
dir_envs="${dir_wet}/envs"
dir_pkgs="${dir_wet}/pkgs"
mkdir -p "${dir_home}"
trap cleanup_wet_root EXIT

if \
    run_isolated_install \
        "${log_create}" \
        --env_nam "${env_nam}" \
        --if_exists update \
        --yes
then
    if [[ -d "${dir_envs}/${env_nam}" ]]; then
        record_pass "install_envs creates env_siqchip in the disposable root"
    else
        record_fail "install_envs did not create env_siqchip in the temp root"
    fi

    if run_isolated_tool samtools --version > /dev/null 2>&1; then
        record_pass "created disposable env_siqchip runs samtools"
    else
        record_fail "created disposable env_siqchip cannot run samtools"
    fi
else
    record_fail \
        "disposable env_siqchip creation failed; see" \
        "$(print_relpath "${log_create}")"
    finish "$@"
    exit $?
fi

if \
    run_isolated_install \
        "${log_update}" \
        --env_nam "${env_nam}" \
        --if_exists update \
        --update_package samtools=1.24 \
        --yes
then
    if run_isolated_tool samtools --version > /dev/null 2>&1; then
        record_pass "install_envs updates disposable env_siqchip"
    else
        record_fail "updated disposable env_siqchip cannot run samtools"
    fi
else
    record_fail \
        "disposable env_siqchip update failed; see" \
        "$(print_relpath "${log_update}")"
fi

finish "$@"
