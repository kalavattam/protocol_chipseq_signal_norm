#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: test_install_envs_layout.sh
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


set -euo pipefail

TEST_NAME="install envs layout"

# Source shared test helpers.
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
            "install -n ${env_lcl}" \
            "install_envs.sh ${env_lcl} update reconciles to its YAML"
        assert_pattern_found \
            "${log_lcl}" \
            "install -n ${env_lcl} -c conda-forge -c bioconda" \
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

log_siqchip_dry_run="${dir_log}/env_siqchip_dry_run.log"
log_analyze_dry_run="${dir_log}/env_analyze_dry_run.log"
log_protocol_channels_dry_run="${dir_log}/env_protocol_channels_dry_run.log"
log_protocol_update_dry_run="${dir_log}/env_protocol_update_dry_run.log"
log_override_no_channels="${dir_log}/install_override_channels_no_channels.log"
log_entrypoint_override_no_channels="${dir_log}/entrypoint_override_channels_no_channels.log"
log_if_exists_bogus="${dir_log}/install_if_exists_bogus.log"
log_render_create="${dir_log}/install_render_create.log"
log_render_update="${dir_log}/install_render_update.log"
log_implied_update="${dir_log}/install_implied_update.log"
log_update_conflict="${dir_log}/install_update_conflict.log"
log_channels_additive="${dir_log}/install_channels_additive.log"
log_tmp_lifecycle="${dir_log}/install_tmp_lifecycle.log"
log_stop_render="${dir_log}/install_stop_render.log"
log_condarc="${dir_log}/install_condarc.log"
log_condarc_none="${dir_log}/install_condarc_none.log"
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

# Channel overrides must reach both the create and the update path. The create
# path cannot express them as flags, because 'env create' rejects
# '--override-channels' and '-c', so they are rendered into a temporary copy of
# the environment YAML instead. The update path uses 'install', which accepts
# the flags directly. Both were broken before 2026-08-12: creation emitted
# flags the subcommand refuses, and the update path rejected '--channels'
# outright.
chn_mirror="https://example.invalid/conda-forge/,https://example.invalid/bioconda/"

if \
    check_cmd_exists mamba || check_cmd_exists conda
then
    if \
        run_capture \
            "install_envs render create channels" \
            "${log_render_create}" \
            "${TEST_BASH}" "${scr_inl}" \
                --dry_run \
                --env_nam env_siqchip \
                --channels "${chn_mirror}" \
                --override_channels
    then
        assert_pattern_found \
            "${log_render_create}" \
            "Rendered YAML:" \
            "install_envs.sh create renders a YAML when channels are given"
        assert_pattern_found \
            "${log_render_create}" \
            "  - https://example.invalid/conda-forge/" \
            "install_envs.sh create renders the supplied channels"
        assert_pattern_found \
            "${log_render_create}" \
            "  - nodefaults" \
            "install_envs.sh create renders nodefaults for override"
        assert_pattern_absent \
            "${log_render_create}" \
            "env create .*--override-channels" \
            "install_envs.sh create passes no channel flags to env create"
        assert_pattern_absent \
            "${log_render_create}" \
            "env create -f ${ROOT_REPO}/install/envs" \
            "install_envs.sh create installs from the rendered copy"
        assert_pattern_found \
            "${log_render_create}" \
            "Rendered condarc:" \
            "install_envs.sh create renders a condarc when channels are given"
        assert_pattern_found \
            "${log_render_create}" \
            "env CONDARC=.* env create" \
            "install_envs.sh create runs env create under the rendered condarc"
    else
        record_fail "install_envs.sh create with channels unexpectedly failed"
    fi

    if \
        run_capture \
            "install_envs render update channels" \
            "${log_render_update}" \
            "${TEST_BASH}" "${scr_inl}" \
                --dry_run \
                --env_nam env_protocol \
                --if_exists update \
                --channels "${chn_mirror}" \
                --override_channels
    then
        assert_pattern_found \
            "${log_render_update}" \
            "-c https://example.invalid/conda-forge/" \
            "install_envs.sh update accepts --channels"
        assert_pattern_absent \
            "${log_render_update}" \
            "-c conda-forge" \
            "install_envs.sh update drops YAML channels when overridden"
        assert_pattern_found \
            "${log_render_update}" \
            "env CONDARC=.* install -n" \
            "install_envs.sh update runs install under the rendered condarc"
    else
        record_fail "install_envs.sh update with channels unexpectedly failed"
    fi
else
    record_skip \
        "install_envs.sh channel-rendering dry-runs require mamba or conda" \
        "on PATH"
fi

# The spec below names a version because 'env_protocol.yml' pins that package,
# and '--update_package' matches a declared specification exactly rather than
# by name. An unpinned package is named bare; a pinned one must carry its pin.
# '--update_package' implies '--if_exists update', because it is meaningful in
# no other mode. An explicitly different '--if_exists' is still a conflict.
# 'update' also no longer freezes installed packages: freezing cannot reconcile
# an environment against a YAML that declares a different version, because the
# solver reports a conflict rather than changing the package.
if \
    check_cmd_exists mamba || check_cmd_exists conda
then
    if \
        run_capture \
            "install_envs update_package implies update" \
            "${log_implied_update}" \
            "${TEST_BASH}" "${scr_inl}" \
                --dry_run \
                --env_nam env_protocol \
                --update_package samtools=1.24
    then
        assert_pattern_found \
            "${log_implied_update}" \
            "would update environment 'env_protocol'" \
            "install_envs.sh --update_package implies --if_exists update"
        assert_pattern_absent \
            "${log_implied_update}" \
            "--freeze-installed" \
            "install_envs.sh update does not freeze installed packages"
    else
        record_fail \
            "install_envs.sh --update_package without --if_exists" \
            "unexpectedly failed"
    fi
else
    record_skip \
        "install_envs.sh implied-update dry-run requires mamba or conda on" \
        "PATH"
fi

if \
    run_capture \
        "install_envs update_package conflicting if_exists" \
        "${log_update_conflict}" \
        "${TEST_BASH}" "${scr_inl}" \
            --dry_run \
            --env_nam env_protocol \
            --if_exists reuse \
            --update_package samtools
then
    record_fail \
        "install_envs.sh --update_package with --if_exists reuse" \
        "unexpectedly succeeded"
else
    assert_pattern_found \
        "${log_update_conflict}" \
        "conflicts with '--if_exists reuse'" \
        "install_envs.sh rejects --update_package with a conflicting mode"
fi

# '--channels' adds channels ahead of the declared ones, which is what '-c'
# means to conda; '--override_channels' is what makes the supplied list
# exclusive. Assert both directions, because collapsing them into "replace"
# silently drops 'conda-forge' and 'bioconda' from every install.
if \
    check_cmd_exists mamba || check_cmd_exists conda
then
    if \
        run_capture \
            "install_envs additive channels" \
            "${log_channels_additive}" \
            "${TEST_BASH}" "${scr_inl}" \
                --dry_run \
                --env_nam env_siqchip \
                --channels https://example.invalid/cf/
    then
        assert_pattern_found \
            "${log_channels_additive}" \
            "  - https://example.invalid/cf/" \
            "install_envs.sh --channels adds the supplied channel"
        assert_pattern_found \
            "${log_channels_additive}" \
            "  - conda-forge" \
            "install_envs.sh --channels retains declared channels"
        assert_pattern_absent \
            "${log_channels_additive}" \
            "nodefaults" \
            "install_envs.sh --channels alone adds no nodefaults"
    else
        record_fail "install_envs.sh additive-channel dry-run failed"
    fi

    # The rendered YAML is a temporary artifact. A dry run keeps it so the
    # caller can read it; nothing may be written beside the tracked YAML, and
    # the tracked YAML itself must never be modified.
    pth_declared="${ROOT_REPO}/install/envs/env_siqchip.yml"
    hsh_before="$(shasum "${pth_declared}" | cut -d' ' -f1)"

    if \
        run_capture \
            "install_envs rendered yaml lifecycle" \
            "${log_tmp_lifecycle}" \
            "${TEST_BASH}" "${scr_inl}" \
                --dry_run \
                --env_nam env_siqchip \
                --channels https://example.invalid/cf/ \
                --override_channels
    then
        pth_rendered="$(
            grep '^Rendered YAML: ' "${log_tmp_lifecycle}" \
                | sed 's|^Rendered YAML: ||'
        )"

        if [[ -n "${pth_rendered}" && -f "${pth_rendered}" ]]; then
            record_pass \
                "install_envs.sh dry run retains the rendered YAML"
        else
            record_fail \
                "install_envs.sh dry run did not retain a readable rendered" \
                "YAML"
        fi

        case "${pth_rendered}" in
            "${ROOT_REPO}"/*)
                record_fail \
                    "install_envs.sh rendered YAML was written inside the" \
                    "repository: '${pth_rendered}'"
                ;;
            *)
                record_pass \
                    "install_envs.sh renders outside the repository tree"
                ;;
        esac

        hsh_after="$(shasum "${pth_declared}" | cut -d' ' -f1)"

        if [[ "${hsh_before}" == "${hsh_after}" ]]; then
            record_pass "install_envs.sh leaves the tracked YAML unmodified"
        else
            record_fail "install_envs.sh modified the tracked YAML"
        fi

        rm -f "${pth_rendered}"
    else
        record_fail "install_envs.sh rendered-YAML dry-run failed"
    fi
else
    record_skip \
        "install_envs.sh channel and rendering dry-runs require mamba or" \
        "conda on PATH"
fi

# 'mirrored_channels' is not a channel source, so '--override_channels' cannot
# reach it: it redirects a supplied channel URL whose final path segment
# matches a mirrored name, and the packages are then fetched from the mirror
# list instead. Miniforge ships one claiming 'conda-forge' for 'anaconda.org',
# which makes every install at a channel-proxying site read the requested
# mirror and download from a host that may be unreachable. The rendered condarc
# empties it. Assert the content, not merely that a file was named, and assert
# that nothing is rendered when no channels are supplied — a configuration file
# imposed on an install nobody asked to redirect is its own defect.
if \
    check_cmd_exists mamba || check_cmd_exists conda
then
    if \
        run_capture \
            "install_envs rendered condarc" \
            "${log_condarc}" \
            "${TEST_BASH}" "${scr_inl}" \
                --dry_run \
                --env_nam env_siqchip \
                --channels https://example.invalid/cf/ \
                --override_channels
    then
        # Tolerate no match. Under 'set -e' a failing 'grep' in a command
        # substitution ends the run, which would turn a missing condarc into an
        # aborted suite rather than a reported failure — and every assertion
        # after this point would silently never execute.
        pth_condarc_rendered="$(
            grep '^Rendered condarc: ' "${log_condarc}" \
                | sed 's|^Rendered condarc: ||' \
                || true
        )"

        if [[ -n "${pth_condarc_rendered}" && -f "${pth_condarc_rendered}" ]]
        then
            record_pass \
                "install_envs.sh dry run retains the rendered condarc"

            if \
                grep -q 'mirrored_channels: {}' "${pth_condarc_rendered}"
            then
                record_pass \
                    "install_envs.sh rendered condarc empties" \
                    "mirrored_channels"
            else
                record_fail \
                    "install_envs.sh rendered condarc does not empty" \
                    "mirrored_channels"
            fi

            case "${pth_condarc_rendered}" in
                "${ROOT_REPO}"/*)
                    record_fail \
                        "install_envs.sh rendered condarc was written inside" \
                        "the repository: '${pth_condarc_rendered}'"
                    ;;
                *)
                    record_pass \
                        "install_envs.sh renders the condarc outside the" \
                        "repository tree"
                    ;;
            esac

            rm -f "${pth_condarc_rendered}"
        else
            record_fail \
                "install_envs.sh dry run did not retain a readable rendered" \
                "condarc"
        fi
    else
        record_fail "install_envs.sh rendered-condarc dry run failed"
    fi

    if \
        run_capture \
            "install_envs no condarc without channels" \
            "${log_condarc_none}" \
            "${TEST_BASH}" "${scr_inl}" \
                --dry_run \
                --env_nam env_siqchip
    then
        assert_pattern_absent \
            "${log_condarc_none}" \
            "Rendered condarc:" \
            "install_envs.sh renders no condarc without --channels"
        assert_pattern_absent \
            "${log_condarc_none}" \
            "CONDARC=" \
            "install_envs.sh sets no CONDARC without --channels"
    else
        record_fail "install_envs.sh no-channel dry run failed"
    fi
else
    record_skip \
        "install_envs.sh condarc-rendering dry-runs require mamba or conda" \
        "on PATH"
fi

# The stop path installs nothing and removes the rendered file on the way out,
# so it must not name a path the caller cannot then read.
if \
    run_capture \
        "install_envs stop path rendered report" \
        "${log_stop_render}" \
        "${TEST_BASH}" "${scr_inl}" \
            --dry_run \
            --env_nam env_protocol \
            --channels https://example.invalid/cf/ \
            --override_channels
then
    assert_pattern_found \
        "${log_stop_render}" \
        "would stop because environment" \
        "install_envs.sh reports the stop path"
    assert_pattern_absent \
        "${log_stop_render}" \
        "Rendered YAML:" \
        "install_envs.sh names no rendered YAML on the stop path"
    assert_pattern_absent \
        "${log_stop_render}" \
        "Rendered condarc:" \
        "install_envs.sh names no rendered condarc on the stop path"
else
    record_fail "install_envs.sh stop-path dry run exited non-zero"
fi

finish
