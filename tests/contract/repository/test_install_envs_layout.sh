#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_install_envs_layout.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in design, development,
# and documentation, with all output reviewed, edited, and approved by the
# author.
#
# Distributed under the MIT license.


set -euo pipefail

TEST_NAME="install envs layout"

#  Source shared test helpers
# shellcheck source=tests/support/test_helpers.sh
source "$(
    git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel
)/tests/support/test_helpers.sh"


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

    if [[ " $* " == *" --if_exists reuse "* ]]; then
        assert_pattern_found \
            "${log_lcl}" \
            "would reuse environment '${env_lcl}'" \
            "install_envs.sh ${env_lcl} dry-run reports reuse"
    elif [[ " $* " == *" --if_exists update "* ]]; then
        assert_pattern_found \
            "${log_lcl}" \
            "would update environment '${env_lcl}' from its YAML" \
            "install_envs.sh ${env_lcl} dry-run reports update"
        assert_pattern_found \
            "${log_lcl}" \
            "install --freeze-installed -n ${env_lcl}" \
            "install_envs.sh ${env_lcl} update freezes installed packages"
        assert_pattern_found \
            "${log_lcl}" \
            "--override-channels -c conda-forge -c bioconda" \
            "install_envs.sh ${env_lcl} update uses YAML channels"
    else
        assert_pattern_found \
            "${log_lcl}" \
            "env create -f" \
            "install_envs.sh ${env_lcl} dry-run uses env create -f"
    fi
}


function assert_install_alias() {
    local interface="${1:?}"
    local shell_cmd="${2:?}"
    local script_path="${3:?}"
    local alias="${4:?}"
    local value="${5:-}"
    local accepted="${6:?}"
    local safe_alias="${alias//[^A-Za-z0-9]/_}"
    local log="${dir_log}/${interface}_${safe_alias}.log"
    local -a arr_command=( "${shell_cmd}" "${script_path}" "${alias}" )

    if [[ -n "${value}" ]]; then
        arr_command+=( "${value}" )
    fi
    arr_command+=( --not_a_real_option )

    run_capture \
        "${interface} parser ${alias}" \
        "${log}" \
        "${arr_command[@]}" || true

    if [[ "${accepted}" == "true" ]]; then
        assert_pattern_found \
            "${log}" \
            "not_a_real_option" \
            "${interface} accepts ${alias} before rejecting the sentinel"
        assert_pattern_absent \
            "${log}" \
            "Unknown.*${alias}|unknown option/parameter passed: '${alias}'" \
            "${interface} does not reject accepted alias ${alias}"
    else
        assert_pattern_found \
            "${log}" \
            "${alias}" \
            "${interface} rejects retired alias ${alias}"
        assert_pattern_absent \
            "${log}" \
            "not_a_real_option" \
            "${interface} rejects ${alias} before the sentinel"
    fi
}


function assert_help_alias_absent() {
    local file="${1:?}"
    local alias="${2:?}"

    if grep -Eq -- "(^|[[:space:],])${alias}([,[:space:]]|$)" "${file}"
    then
        record_fail \
            "$(print_relpath "${file}") unexpectedly exposes ${alias}"
    else
        record_pass "$(print_relpath "${file}") hides ${alias}"
    fi
}


dir_log="${TEST_DIR_LOG}/install_envs"

scr_inl="${ROOT_REPO}/install/scripts/install_envs.sh"
scr_ent="${ROOT_REPO}/install/scripts/install_envs_entrypoint.sh"
scr_atr="${ROOT_REPO}/install/scripts/install_atria.sh"
scr_hlp="${ROOT_REPO}/lib/bash/help/help_install_envs.sh"

yml_anl="${ROOT_REPO}/install/envs/env_analyze.yml"
yml_prt="${ROOT_REPO}/install/envs/env_protocol.yml"
yml_siq="${ROOT_REPO}/install/envs/env_siqchip.yml"

log_atria_dry_run="${dir_log}/install_atria_dry_run.log"
log_siqchip_dry_run="${dir_log}/env_siqchip_dry_run.log"
log_analyze_dry_run="${dir_log}/env_analyze_dry_run.log"
log_protocol_channels_dry_run="${dir_log}/env_protocol_channels_dry_run.log"
log_protocol_update_dry_run="${dir_log}/env_protocol_update_dry_run.log"
log_override_no_channels="${dir_log}/install_override_channels_no_channels.log"
log_entrypoint_override_no_channels="${dir_log}/entrypoint_override_channels_no_channels.log"
log_if_exists_bogus="${dir_log}/install_if_exists_bogus.log"
log_install_help="${dir_log}/install_envs_help_aliases.log"
log_entrypoint_help="${dir_log}/install_envs_entrypoint_help_aliases.log"


print_section "${TEST_NAME}"

mkdir -p "${dir_log}"

for spec in \
    "bash -n install_envs:${TEST_BASH}:${scr_inl}" \
    "sh -n install_envs_entrypoint:sh:${scr_ent}" \
    "bash -n install_atria:${TEST_BASH}:${scr_atr}" \
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
    "install_envs_entrypoint help:sh:${scr_ent}" \
    "install_atria help:${TEST_BASH}:${scr_atr}"
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
            '^Usage$' \
            "${label} prints Usage"
    else
        record_fail "${label} failed; see $(print_relpath "${log}")"
    fi
done

run_capture \
    "install_envs canonical help aliases" \
    "${log_install_help}" \
    "${TEST_BASH}" "${scr_inl}" --help || true
run_capture \
    "install_envs_entrypoint canonical help aliases" \
    "${log_entrypoint_help}" \
    sh "${scr_ent}" --help || true

for log in "${log_install_help}" "${log_entrypoint_help}"; do
    assert_pattern_found \
        "${log}" \
        '^  -ch, --channels : list of str$' \
        "$(print_relpath "${log}") exposes canonical channel aliases"
    assert_pattern_found \
        "${log}" \
        '^  -oc, --override_channels : flag$' \
        "$(print_relpath "${log}") exposes canonical override aliases"
    for alias in \
        --channel \
        --channel_list \
        --channel-list \
        --override_channel \
        --override-channel \
        --override-channels
    do
        assert_help_alias_absent "${log}" "${alias}"
    done
done

for interface in bash entrypoint; do
    if [[ "${interface}" == "bash" ]]; then
        shell_cmd="${TEST_BASH}"
        script_path="${scr_inl}"
    else
        shell_cmd="sh"
        script_path="${scr_ent}"
    fi

    assert_install_alias \
        "${interface}" "${shell_cmd}" "${script_path}" -ch conda-forge true
    assert_install_alias \
        "${interface}" "${shell_cmd}" "${script_path}" \
        --channels conda-forge true
    assert_install_alias \
        "${interface}" "${shell_cmd}" "${script_path}" -oc "" true
    assert_install_alias \
        "${interface}" "${shell_cmd}" "${script_path}" \
        --override_channels "" true
    assert_install_alias \
        "${interface}" "${shell_cmd}" "${script_path}" \
        --override-channels "" true

    for alias in \
        --channel \
        --channel_list \
        --channel-list \
        --override_channel \
        --override-channel
    do
        assert_install_alias \
            "${interface}" "${shell_cmd}" "${script_path}" \
            "${alias}" "" false
    done
done

if \
    run_capture \
        "install_atria dry-run" \
        "${log_atria_dry_run}" \
        "${TEST_BASH}" "${scr_atr}" \
            --dry_run \
            --if_exists reuse \
            --dir_install "${TEST_DIR_TMP}/install_atria_layout"
then
    record_pass "install_atria.sh dry-run exits 0"
    assert_pattern_found \
        "${log_atria_dry_run}" \
        "v_atria=4.1.5" \
        "install_atria.sh dry-run reports Atria 4.1.5"

    assert_pattern_found \
        "${log_atria_dry_run}" \
        "tag_atr=v4.1.5" \
        "install_atria.sh dry-run maps Atria tag v4.1.5"

    assert_pattern_found \
        "${log_atria_dry_run}" \
        "atria-4.1.5/bin" \
        "install_atria.sh dry-run reports provisional Atria bin path"
else
    record_fail \
        "install_atria.sh dry-run failed; see" \
        "$(print_relpath "${log_atria_dry_run}")"
fi

assert_readable_yaml "${yml_anl}" "env_analyze YAML"
assert_readable_yaml "${yml_prt}" "env_protocol YAML"
assert_readable_yaml "${yml_siq}" "env_siqchip YAML"
assert_pattern_found \
    "${yml_prt}" \
    '^  - shellcheck=0\.10\.0' \
    "env_protocol pins ShellCheck 0.10.0"

if \
    check_cmd_exists mamba || check_cmd_exists conda
then
    assert_install_dry_run \
        "env_siqchip" "${yml_siq}" "${log_siqchip_dry_run}" \
        --dry_run \
        --env_nam env_siqchip \
        --yes

    assert_pattern_found \
        "${log_siqchip_dry_run}" \
        "--yes" \
        "install_envs.sh env_siqchip dry-run includes --yes"

    assert_install_dry_run \
        "env_analyze" "${yml_anl}" "${log_analyze_dry_run}" \
        --dry_run \
        --env_nam env_analyze \
        --yes

    assert_install_dry_run \
        "env_protocol" "${yml_prt}" \
        "${log_protocol_channels_dry_run}" \
        --dry_run \
        --env_nam env_protocol \
        --if_exists reuse \
        --channels fhcc-main,fhcc-bioconda \
        --override_channels \
        --yes

    assert_pattern_found \
        "${log_protocol_channels_dry_run}" \
        "Package command:.* run -n env_protocol python -m pip install" \
        "install_envs.sh reuse dry-run reports package refresh"

    assert_pattern_found \
        "${log_protocol_channels_dry_run}" \
        "--editable ${ROOT_REPO}" \
        "install_envs.sh package refresh targets repository root"

    assert_pattern_found \
        "${log_protocol_channels_dry_run}" \
        "--no-deps.*--no-build-isolation" \
        "install_envs.sh package refresh disables dependency and build isolation"

    assert_install_dry_run \
        "env_protocol" "${yml_prt}" \
        "${log_protocol_update_dry_run}" \
        --dry_run \
        --env_nam env_protocol \
        --if_exists update \
        --update_package shellcheck=0.10.0 \
        --yes

    assert_pattern_found \
        "${log_protocol_update_dry_run}" \
        "shellcheck=0.10.0" \
        "install_envs.sh targeted update includes declared ShellCheck pin"

    assert_pattern_absent \
        "${log_protocol_update_dry_run}" \
        " asciigenome " \
        "install_envs.sh targeted update excludes unrelated YAML packages"
else
    record_skip \
        "install_envs.sh command-construction dry-runs require mamba or" \
        "conda on PATH"
fi

if \
    run_capture \
        "install_envs override_channels without channels" \
        "${log_override_no_channels}" \
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
        "${log_override_no_channels}" \
        "--override_channels.*requires.*--channels" \
        "install_envs.sh rejects --override_channels without --channels"
fi

if \
    check_cmd_exists mamba || check_cmd_exists conda
then
    if \
        run_capture \
            "entrypoint override_channels without channels" \
            "${log_entrypoint_override_no_channels}" \
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
            "${log_entrypoint_override_no_channels}" \
            "--override_channels.*requires.*--channels" \
            "install_envs_entrypoint.sh rejects --override_channels without --channels"
    fi
else
    record_skip \
        "entrypoint negative override-channel check requires mamba or conda" \
        "for handoff"
fi

if \
    run_capture \
        "install_envs invalid if_exists" \
        "${log_if_exists_bogus}" \
        "${TEST_BASH}" "${scr_inl}" \
            --dry_run \
            --env_nam env_protocol \
            --if_exists bogus
then
    record_fail "install_envs.sh invalid --if_exists unexpectedly succeeded"
else
    assert_pattern_found \
        "${log_if_exists_bogus}" \
        "invalid '--if_exists' value" \
        "install_envs.sh rejects invalid --if_exists"
fi

finish
