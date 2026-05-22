#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_install_envs_layout.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5.5) was used in development.
#
# Distributed under the MIT license.


TEST_NAME="install envs layout"

#  Source shared smoke-test helpers
# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." > /dev/null 2>&1 && pwd
)/lib/test_helpers.sh"

rec_section "${TEST_NAME}"

dir_log="${TEST_DIR_LOG}/install_envs"
mkdir -p "${dir_log}"

scr_install="${ROOT_REPO}/install/scripts/install_envs.sh"
scr_entry="${ROOT_REPO}/install/scripts/install_envs_entrypoint.sh"
scr_help="${ROOT_REPO}/scripts/functions/help/help_install_envs.sh"
yml_analyze="${ROOT_REPO}/install/envs/env_analyze.yml"
yml_protocol="${ROOT_REPO}/install/envs/env_protocol.yml"
yml_siqchip="${ROOT_REPO}/install/envs/env_siqchip.yml"


function assert_readable_yaml() {
    local file="${1:-}"
    local label="${2:-YAML}"

    if [[ -f "${file}" && -r "${file}" ]]; then
        rec_pass "${label} exists and is readable"
    else
        rec_fail "${label} missing or unreadable: $(rec_relpath "${file}")"
    fi
}


function assert_install_dry_run() {
    local env_lcl="${1:-}"
    local yml_lcl="${2:-}"
    local log_lcl="${3:-}"
    shift 3

    if ! \
        run_capture \
            "install_envs ${env_lcl} dry-run" "${log_lcl}" \
            "${TEST_BASH}" "${scr_install}" "$@"
    then
        rec_fail \
            "install_envs.sh ${env_lcl} dry-run failed; see" \
            "$(rec_relpath "${log_lcl}")"
        return
    fi

    rec_pass "install_envs.sh ${env_lcl} dry-run exits 0"
    assert_grep_pattern \
        "${log_lcl}" \
        "YAML: ${yml_lcl}" \
        "install_envs.sh ${env_lcl} dry-run reports YAML path"
    assert_grep_pattern \
        "${log_lcl}" \
        "env create -f" \
        "install_envs.sh ${env_lcl} dry-run uses env create -f"
}


for spec in \
    "bash -n install_envs:${TEST_BASH}:${scr_install}" \
    "sh -n install_envs_entrypoint:sh:${scr_entry}" \
    "bash -n help_install_envs:${TEST_BASH}:${scr_help}"
do
    label="${spec%%:*}"
    rest="${spec#*:}"
    shell_cmd="${rest%%:*}"
    file="${rest#*:}"
    log="${dir_log}/${label// /_}.log"

    if run_capture "${label}" "${log}" "${shell_cmd}" -n "${file}"; then
        rec_pass "${label}"
    else
        rec_fail "${label}; see $(rec_relpath "${log}")"
    fi
done

for spec in \
    "install_envs help:${TEST_BASH}:${scr_install}" \
    "install_envs_entrypoint help:sh:${scr_entry}"
do
    label="${spec%%:*}"
    rest="${spec#*:}"
    shell_cmd="${rest%%:*}"
    file="${rest#*:}"
    log="${dir_log}/${label// /_}.log"

    if run_capture "${label}" "${log}" "${shell_cmd}" "${file}" --help; then
        rec_pass "${label} exits 0"
        assert_grep_pattern "${log}" '^Usage:' "${label} prints Usage"
    else
        rec_fail "${label} failed; see $(rec_relpath "${log}")"
    fi
done

assert_readable_yaml "${yml_analyze}" "env_analyze YAML"
assert_readable_yaml "${yml_protocol}" "env_protocol YAML"
assert_readable_yaml "${yml_siqchip}" "env_siqchip YAML"

if check_cmd_exists mamba || check_cmd_exists conda; then
    log="${dir_log}/env_siqchip_dry_run.log"
    assert_install_dry_run \
        "env_siqchip" "${yml_siqchip}" "${log}" \
        --dry_run \
        --env_nam env_siqchip \
        --yes
    assert_grep_pattern \
        "${log}" \
        "--yes" \
        "install_envs.sh env_siqchip dry-run includes --yes"

    log="${dir_log}/env_analyze_dry_run.log"
    assert_install_dry_run \
        "env_analyze" "${yml_analyze}" "${log}" \
        --dry_run \
        --env_nam env_analyze \
        --yes

    log="${dir_log}/env_protocol_channels_dry_run.log"
    assert_install_dry_run \
        "env_protocol channels" "${yml_protocol}" "${log}" \
        --dry_run \
        --env_nam env_protocol \
        --channels fhcc-main,fhcc-bioconda \
        --override_channels \
        --yes
    assert_grep_pattern \
        "${log}" \
        "--override-channels" \
        "install_envs.sh channel dry-run includes --override-channels"
    assert_grep_pattern \
        "${log}" \
        "-c fhcc-main" \
        "install_envs.sh channel dry-run includes fhcc-main"
    assert_grep_pattern \
        "${log}" \
        "-c fhcc-bioconda" \
        "install_envs.sh channel dry-run includes fhcc-bioconda"
    assert_grep_pattern \
        "${log}" \
        "--yes" \
        "install_envs.sh channel dry-run includes --yes"
else
    rec_skip \
        "install_envs.sh command-construction dry-runs require mamba or conda" \
        "on PATH"
fi

log="${dir_log}/install_override_channels_no_channels.log"
if \
    run_capture \
        "install_envs override_channels without channels" "${log}" \
        "${TEST_BASH}" "${scr_install}" \
            --dry_run \
            --env_nam env_protocol \
            --override_channels
then
    rec_fail "install_envs.sh --override_channels without --channels unexpectedly succeeded"
else
    assert_grep_pattern \
        "${log}" \
        "--override_channels.*requires.*--channels" \
        "install_envs.sh rejects --override_channels without --channels"
fi

log="${dir_log}/entrypoint_override_channels_no_channels.log"
if check_cmd_exists mamba || check_cmd_exists conda; then
    if \
        run_capture \
            "entrypoint override_channels without channels" "${log}" \
            sh "${scr_entry}" \
                --dry_run \
                --env_nam env_protocol \
                --override_channels
    then
        rec_fail \
            "install_envs_entrypoint.sh --override_channels without --channels" \
            "unexpectedly succeeded"
    else
        assert_grep_pattern \
            "${log}" \
            "--override_channels.*requires.*--channels" \
            "install_envs_entrypoint.sh rejects --override_channels without --channels"
    fi
else
    rec_skip \
        "entrypoint negative override-channel check requires mamba or conda" \
        "for handoff"
fi

log="${dir_log}/install_if_exis_bogus.log"
if \
    run_capture \
        "install_envs invalid if_exis" "${log}" \
        "${TEST_BASH}" "${scr_install}" \
            --dry_run \
            --env_nam env_protocol \
            --if_exis bogus
then
    rec_fail "install_envs.sh invalid --if_exis unexpectedly succeeded"
else
    assert_grep_pattern \
        "${log}" \
        "invalid '--if_exis' value" \
        "install_envs.sh rejects invalid --if_exis"
fi

finish
