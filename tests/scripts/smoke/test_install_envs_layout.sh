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

print_section "${TEST_NAME}"

dir_log="${TEST_DIR_LOG}/install_envs"
mkdir -p "${dir_log}"

scr_inl="${ROOT_REPO}/install/scripts/install_envs.sh"
scr_ent="${ROOT_REPO}/install/scripts/install_envs_entrypoint.sh"
scr_hlp="${ROOT_REPO}/scripts/functions/help/help_install_envs.sh"
yml_anl="${ROOT_REPO}/install/envs/env_analyze.yml"
yml_prt="${ROOT_REPO}/install/envs/env_protocol.yml"
yml_siq="${ROOT_REPO}/install/envs/env_siqchip.yml"


function assert_readable_yaml() {
    local file="${1:-}"
    local label="${2:-YAML}"

    if [[ -f "${file}" && -r "${file}" ]]; then
        record_pass "${label} exists and is readable"
    else
        record_fail \
            "${label} missing or unreadable: $(print_relpath "${file}")"
    fi
}


function assert_install_dry_run() {
    local env_lcl="${1:-}"
    local yml_lcl="${2:-}"
    local log_lcl="${3:-}"
    shift 3

    if ! \
        run_capture \
            "install_envs ${env_lcl} dry-run" \
            "${log_lcl}" \
            "${TEST_BASH}" "${scr_inl}" "$@"
    then
        record_fail \
            "install_envs.sh ${env_lcl} dry-run failed; see" \
            "$(print_relpath "${log_lcl}")"
        return
    fi

    record_pass "install_envs.sh ${env_lcl} dry-run exits 0"
    assert_pattern_found \
        "${log_lcl}" \
        "YAML: ${yml_lcl}" \
        "install_envs.sh ${env_lcl} dry-run reports YAML path"
    assert_pattern_found \
        "${log_lcl}" \
        "env create -f" \
        "install_envs.sh ${env_lcl} dry-run uses env create -f"
}


for spec in \
    "bash -n install_envs:${TEST_BASH}:${scr_inl}" \
    "sh -n install_envs_entrypoint:sh:${scr_ent}" \
    "bash -n help_install_envs:${TEST_BASH}:${scr_hlp}"
do
    label="${spec%%:*}"
    rest="${spec#*:}"
    shell_cmd="${rest%%:*}"
    file="${rest#*:}"
    log="${dir_log}/${label// /_}.log"

    if \
        run_capture \
            "${label}" \
            "${log}" \
            "${shell_cmd}" -n "${file}"
    then
        record_pass "${label}"
    else
        record_fail "${label}; see $(print_relpath "${log}")"
    fi
done

for spec in \
    "install_envs help:${TEST_BASH}:${scr_inl}" \
    "install_envs_entrypoint help:sh:${scr_ent}"
do
    label="${spec%%:*}"
    rest="${spec#*:}"
    shell_cmd="${rest%%:*}"
    file="${rest#*:}"
    log="${dir_log}/${label// /_}.log"

    if \
        run_capture \
            "${label}" \
            "${log}" \
            "${shell_cmd}" "${file}" --help
    then
        record_pass "${label} exits 0"
        assert_pattern_found \
            "${log}" \
            '^Usage:' \
            "${label} prints Usage"
    else
        record_fail "${label} failed; see $(print_relpath "${log}")"
    fi
done

assert_readable_yaml "${yml_anl}" "env_analyze YAML"
assert_readable_yaml "${yml_prt}" "env_protocol YAML"
assert_readable_yaml "${yml_siq}" "env_siqchip YAML"

if \
    check_cmd_exists mamba || check_cmd_exists conda
then
    log="${dir_log}/env_siqchip_dry_run.log"
    assert_install_dry_run \
        "env_siqchip" "${yml_siq}" "${log}" \
        --dry_run \
        --env_nam env_siqchip \
        --yes
    assert_pattern_found \
        "${log}" \
        "--yes" \
        "install_envs.sh env_siqchip dry-run includes --yes"

    log="${dir_log}/env_analyze_dry_run.log"
    assert_install_dry_run \
        "env_analyze" "${yml_anl}" "${log}" \
        --dry_run \
        --env_nam env_analyze \
        --yes

    log="${dir_log}/env_protocol_channels_dry_run.log"
    assert_install_dry_run \
        "env_protocol channels" "${yml_prt}" "${log}" \
        --dry_run \
        --env_nam env_protocol \
        --channels fhcc-main,fhcc-bioconda \
        --override_channels \
        --yes
    assert_pattern_found \
        "${log}" \
        "--override-channels" \
        "install_envs.sh channel dry-run includes --override-channels"
    assert_pattern_found \
        "${log}" \
        "-c fhcc-main" \
        "install_envs.sh channel dry-run includes fhcc-main"
    assert_pattern_found \
        "${log}" \
        "-c fhcc-bioconda" \
        "install_envs.sh channel dry-run includes fhcc-bioconda"
    assert_pattern_found \
        "${log}" \
        "--yes" \
        "install_envs.sh channel dry-run includes --yes"
else
    record_skip \
        "install_envs.sh command-construction dry-runs require mamba or" \
        "conda on PATH"
fi

log="${dir_log}/install_override_channels_no_channels.log"
if \
    run_capture \
        "install_envs override_channels without channels" \
        "${log}" \
        "${TEST_BASH}" "${scr_inl}" \
            --dry_run \
            --env_nam env_protocol \
            --override_channels
then
    record_fail \
        "install_envs.sh --override_channels without --channels unexpectedly" \
        "succeeded"
else
    assert_pattern_found \
        "${log}" \
        "--override_channels.*requires.*--channels" \
        "install_envs.sh rejects --override_channels without --channels"
fi

log="${dir_log}/entrypoint_override_channels_no_channels.log"
if \
    check_cmd_exists mamba || check_cmd_exists conda
then
    if \
        run_capture \
            "entrypoint override_channels without channels" \
            "${log}" \
            sh "${scr_ent}" \
                --dry_run \
                --env_nam env_protocol \
                --override_channels
    then
        record_fail \
            "install_envs_entrypoint.sh --override_channels without" \
            "--channels" \
            "unexpectedly succeeded"
    else
        assert_pattern_found \
            "${log}" \
            "--override_channels.*requires.*--channels" \
            "install_envs_entrypoint.sh rejects --override_channels without --channels"
    fi
else
    record_skip \
        "entrypoint negative override-channel check requires mamba or conda" \
        "for handoff"
fi

log="${dir_log}/install_if_exis_bogus.log"
if \
    run_capture \
        "install_envs invalid if_exis" \
        "${log}" \
        "${TEST_BASH}" "${scr_inl}" \
            --dry_run \
            --env_nam env_protocol \
            --if_exis bogus
then
    record_fail "install_envs.sh invalid --if_exis unexpectedly succeeded"
else
    assert_pattern_found \
        "${log}" \
        "invalid '--if_exis' value" \
        "install_envs.sh rejects invalid --if_exis"
fi

finish
