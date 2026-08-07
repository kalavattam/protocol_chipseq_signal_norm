#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: install_envs.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
# GPT-5.6) were used in design, development, and documentation, with all output
# reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


#  Require Bash >= 4.4 before doing any work
if [[ -z "${BASH_VERSION:-}" ]]; then
    echo "error(shell):" \
        "this script must be run under Bash >= 4.4." >&2
    exit 1
elif (( BASH_VERSINFO[0] < 4 || (BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] < 4) )); then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "this script requires Bash >= 4.4; current version is" \
        "'${BASH_VERSION}'." >&2
    exit 1
fi

#  Run in safe mode, exiting on errors, unset variables, and pipe failures
set -euo pipefail

#  Set paths to installation support and repository directories
dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
dir_ins="$(cd "${dir_scr}/.." > /dev/null 2>&1 && pwd)"
dir_rep="$(cd "${dir_ins}/.." > /dev/null 2>&1 && pwd)"


#  Source shared helpers
function source_helpers_script() {
    local fnc_src

    dir_fnc="${dir_rep}/lib/bash"
    fnc_src="${dir_fnc}/core/source_helpers.sh"

    if [[ ! -f "${fnc_src}" ]]; then
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "script not found: '${fnc_src}'." >&2
        return 1
    fi

    # shellcheck disable=SC1090
    source "${fnc_src}" || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source '${fnc_src}'." >&2
        return 1
    }

    source_helpers "${dir_fnc}" \
        check_args \
        check_env \
        check_inputs \
        format_outputs \
        handle_env \
        help/help_install_envs \
        || {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to source required helper scripts." >&2
            return 1
        }
}


# shellcheck disable=SC2120
function check_pkg_mgr() {
    local arg="${1:-}"
    local show_help

    show_help=$(cat << EOM
Usage
-----
  check_pkg_mgr
    [--help]

  Check that either Mamba or Conda is installed and available in PATH.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

Returns
-------
  0 if Mamba or Conda is available; otherwise, 1 and an error message.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - mamba or conda

Examples
--------
  1. Verify that Mamba or Conda is available on the current PATH.
    '''bash
    check_pkg_mgr
    '''

  2. Confirm that positional arguments are rejected.
    '''bash
    if ! check_pkg_mgr conda; then
        printf '%s\n' 'unexpected argument rejected as expected'
    fi
    '''
EOM
    )

    if [[ "${arg}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -n "${arg}" ]]; then
        echo_err "unexpected argument to 'check_pkg_mgr()': '${arg}'."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    if command -v mamba >/dev/null 2>&1; then
        return 0
    elif command -v conda >/dev/null 2>&1; then
        return 0
    else
        echo_err "neither Mamba nor Conda is installed on the system."
        echo >&2
        echo \
            "Mamba is a package manager that makes package installations" \
            "faster and more reliable in comparison to Conda." >&2
        echo >&2
        echo \
            "For Mamba installation instructions, please check the following" \
            "link: https://github.com/mamba-org/mamba#installation" >&2

        return 1
    fi
}


function init_arg_defs() {
    dry_run=false
    env_nam=""
    if_exists="fail"
    channels=""
    override_channels=false
    yes=false
    pth_yml=""
    env_action=create
    pkg_mgr=""

    arr_channels=()
    cmd=()
    packages=()
    requested_packages=()
    cmd_pkg=()
    declare -ga arr_channels cmd packages requested_packages cmd_pkg
}


function parse_args() {
    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        help_install_envs >&2
        return 2
    fi

    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
            -dr|--dry|--dry[_-]run)
                dry_run=true
                shift 1
                ;;

            -en|--env|--env[_-]nam)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_install_envs >&2
                    return 1
                }
                env_nam="${2}"
                shift 2
                ;;

            -ie|--if[_-]exists)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_install_envs >&2
                    return 1
                }
                if_exists="${2}"
                shift 2
                ;;

            -up|--update[_-]package)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_install_envs >&2
                    return 1
                }
                requested_packages+=( "${2}" )
                shift 2
                ;;

            -ch|--channels)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_install_envs >&2
                    return 1
                }
                channels="${2}"
                shift 2
                ;;

            -oc|--override[_-]channels)
                override_channels=true
                shift 1
                ;;

            -y|--yes)
                yes=true
                shift 1
                ;;

            *)
                echo_err "unknown option/parameter passed: '${1}'."
                echo >&2
                help_install_envs >&2
                return 1
                ;;
        esac
    done
}


function validate_args() {
    local channel=""

    validate_var "env_nam" "${env_nam}" || return 1
    validate_var "if_exists" "${if_exists}" || return 1

    if [[ -n "${channels}" ]]; then
        IFS=',' read -r -a arr_channels <<< "${channels}"

        for channel in "${arr_channels[@]}"; do
            if [[ -z "${channel}" ]]; then
                echo_err \
                    "invalid '--channels' value: '${channels}'. Channel" \
                    "names must be comma-delimited and non-empty."
                return 1
            fi
        done
    else
        arr_channels=()
    fi

    if [[
        "${override_channels}" == "true" && "${#arr_channels[@]}" -eq 0
    ]]; then
        echo_err "'--override_channels' requires '--channels'."
        return 1
    fi

    case "${if_exists}" in
        fail|reuse|update) : ;;
        *)
            echo_err \
                "invalid '--if_exists' value: '${if_exists}'. Must be 'fail'" \
                "'reuse', or 'update'."
            return 1
            ;;
    esac

    if [[ "${if_exists}" == "update" && ${#arr_channels[@]} -gt 0 ]]; then
        echo_err "'--if_exists update' uses channels from the environment YAML."
        return 1
    fi

    if [[
        "${if_exists}" != "update" && ${#requested_packages[@]} -gt 0
    ]]; then
        echo_err "'--update_package' requires '--if_exists update'."
        return 1
    fi

    case "${env_nam}" in
        env_align|env_analyze|env_protocol|env_repro|env_siqchip) : ;;
        *)
            ## NOTE: 'env_align' and 'env_repro' are not exposed to users ##
            echo_err \
                "invalid environment name specified. Must be 'env_analyze'," \
                "'env_protocol', or 'env_siqchip'."
            return 1
            ;;
    esac
}


function resolve_env_definition() {
    case "${env_nam}" in
        env_analyze|env_protocol|env_siqchip)
            pth_yml="${dir_rep}/install/envs/${env_nam}.yml"

            if [[ ! -f "${pth_yml}" ]]; then
                echo_err "environment YAML does not exist: '${pth_yml}'."
                return 1
            fi

            if [[ ! -r "${pth_yml}" ]]; then
                echo_err "environment YAML is not readable: '${pth_yml}'."
                return 1
            fi
            ;;
    esac
}


function load_yaml_update_specs() {
    local line=""
    local section=""
    local value=""
    local -a yaml_channels=()
    local -a yaml_packages=()

    while IFS= read -r line || [[ -n "${line}" ]]; do
        case "${line}" in
            channels:)
                section=channels
                continue
                ;;

            dependencies:)
                section=dependencies
                continue
                ;;

            [![:space:]]*)
                section=""
                ;;
        esac

        if [[ "${line}" =~ ^[[:space:]]{2}-[[:space:]]+(.+)$ ]]; then
            value="${BASH_REMATCH[1]%%#*}"
            value="${value%"${value##*[![:space:]]}"}"
            value="${value#"${value%%[![:space:]]*}"}"

            if [[ -z "${value}" || "${value}" == *:* ]]; then
                echo_err \
                    "unsupported nested or empty YAML entry in '${pth_yml}':" \
                    "'${line}'."
                return 1
            fi

            case "${section}" in
                channels) yaml_channels+=( "${value}" ) ;;
                dependencies) yaml_packages+=( "${value}" ) ;;
                *)
                    echo_err \
                        "unexpected YAML list entry in '${pth_yml}':" \
                        "'${line}'."
                    return 1
                    ;;
            esac
        elif [[ "${line}" =~ ^[[:space:]]+-[[:space:]]+ ]]; then
            echo_err \
                "unsupported YAML indentation in '${pth_yml}': '${line}'."
            return 1
        fi
    done < "${pth_yml}"

    if (( ${#yaml_channels[@]} == 0 || ${#yaml_packages[@]} == 0 )); then
        echo_err \
            "environment YAML must declare non-empty channels and" \
            "dependencies: '${pth_yml}'."
        return 1
    fi

    arr_channels=( "${yaml_channels[@]}" )
    packages=( "${yaml_packages[@]}" )
}


function build_install_cmd() {
    local channel=""

    # shellcheck disable=SC2119
    check_pkg_mgr || return 1

    if command -v mamba >/dev/null 2>&1; then
        pkg_mgr=mamba
        if [[ -n "${pth_yml}" ]]; then
            cmd=( mamba env create -f "${pth_yml}" )
        else
            cmd=( mamba create -n "${env_nam}" )
        fi
    else
        pkg_mgr=conda
        if [[ -n "${pth_yml}" ]]; then
            cmd=( conda env create -f "${pth_yml}" )
        else
            cmd=( conda create -n "${env_nam}" )
        fi
    fi

    if [[ "${override_channels}" == "true" ]]; then
        cmd+=( --override-channels )
    fi

    for channel in "${arr_channels[@]}"; do
        cmd+=( -c "${channel}" )
    done

    if [[ "${yes}" == "true" ]]; then
        cmd+=( --yes )
    fi

    if [[ -z "${pth_yml}" && "${env_nam}" == "env_align" ]]; then
        packages=(  ## NOTE: Retained for old work; not exposed in the docs ##
            bamtools
            bbmap
            bedtools
            bowtie2
            bwa
            datamash
            fastqc
            gawk
            gnuplot
            macs3
            minimap
            mosdepth
            parallel
            picard
            preseq
            rename
            samtools
            subread
            tree
            ucsc-bedgraphtobigwig
            ucsc-bedsort
            ucsc-facount
            wget
        )
    elif [[ -z "${pth_yml}" && "${env_nam}" == "env_repro" ]]; then
        packages=(  ## NOTE: Not exposing this to users in the docs ##
            bc
            bowtie2=2.3.4.2  ## NOTE: Explicitly pinning old version ##
            deeptools=3.3.1  ## NOTE: Explicitly pinning old version ##
            gawk
            ipython
            parallel
            pbzip2
            pigz
            python=3.6       ## NOTE: Explicitly pinning old version ##
            rename
            samtools=1.9     ## NOTE: Explicitly pinning old version ##
            tree
            wget
        )
    fi

    if [[ "${env_nam}" == "env_protocol" ]]; then
        cmd_pkg=(
            "${pkg_mgr}" run -n "${env_nam}"
                python -m pip install
                    --editable "${dir_rep}"
                    --no-deps
                    --no-build-isolation
        )
    fi
}


function print_dry_run() {
    if [[ "${env_action}" == "reuse" ]]; then
        echo "dryrun($(basename "${BASH_SOURCE[0]}")):" \
            "would reuse environment '${env_nam}'."
    elif [[ "${env_action}" == "update" ]]; then
        echo "dryrun($(basename "${BASH_SOURCE[0]}")):" \
            "would update environment '${env_nam}' from its YAML."
    elif [[ "${env_action}" == "stop" ]]; then
        echo "dryrun($(basename "${BASH_SOURCE[0]}")):" \
            "would stop because environment '${env_nam}' exists."
    else
        echo "dryrun($(basename "${BASH_SOURCE[0]}")):" \
            "would create environment '${env_nam}'."
    fi

    if [[ -n "${pth_yml}" ]]; then
        echo "YAML: ${pth_yml}"
    fi

    if [[ "${env_action}" == "create" || "${env_action}" == "update" ]]; then
        printf 'Environment command:'
        for tok in "${cmd[@]}" "${packages[@]}"; do
            printf ' %q' "${tok}"
        done
        printf '\n'
    fi

    if [[ "${env_action}" != "stop" ]] && (( ${#cmd_pkg[@]} > 0 )); then
        printf 'Package command:'
        for tok in "${cmd_pkg[@]}"; do
            printf ' %q' "${tok}"
        done
        printf '\n'
    fi
}


function handle_existing_env() {
    local found=false
    local requested=""
    local yaml_package=""
    local -a yaml_packages=()

    if ! \
        check_env_installed "${env_nam}" "true"
    then
        return 0
    fi

    case "${if_exists}" in
        fail)
            env_action=stop
            echo_err \
                "an environment with the name '${env_nam}' is already" \
                "installed."
            echo >&2

            if [[ "${dry_run}" == "true" ]]; then
                print_dry_run
                echo >&2
                echo "dryrun($(basename "${BASH_SOURCE[0]}")):" \
                    "a non-dry run would stop here unless '--if_exists reuse'" \
                    "was specified." >&2
                echo >&2
                return 2
            fi

            echo \
                "Nothing was changed. To reuse the existing environment," \
                "rerun with '--if_exists reuse'." >&2
            echo >&2

            echo \
                "To rebuild '${env_nam}', remove the existing environment" \
                "first, then rerun this script." >&2
            echo >&2

            echo \
                "For example, run 'mamba env remove -n \"${env_nam}\"'" \
                "prior to rerunning the script." >&2
            return 1
            ;;

        reuse)
            env_action=reuse
            echo \
                "Environment '${env_nam}' already exists; reusing it because" \
                "'--if_exists reuse' was specified."
            return 0
            ;;

        update)
            if [[ -z "${pth_yml}" ]]; then
                echo_err \
                    "'--if_exists update' requires a YAML-backed environment."
                return 1
            fi

            env_action=update

            load_yaml_update_specs || return 1

            if (( ${#requested_packages[@]} > 0 )); then
                yaml_packages=( "${packages[@]}" )
                packages=()

                for requested in "${requested_packages[@]}"; do
                    found=false

                    for yaml_package in "${yaml_packages[@]}"; do
                        if [[ "${requested}" == "${yaml_package}" ]]; then
                            found=true
                            break
                        fi
                    done

                    if [[ "${found}" != "true" ]]; then
                        echo_err \
                            "requested update package is not declared exactly" \
                            "in '${pth_yml}': '${requested}'."
                        return 1
                    fi

                    packages+=( "${requested}" )
                done
            fi

            cmd=(
                "${pkg_mgr}" install
                    --freeze-installed
                    -n "${env_nam}"
                    --override-channels
            )

            for channel in "${arr_channels[@]}"; do
                cmd+=( -c "${channel}" )
            done

            if [[ "${yes}" == "true" ]]; then
                cmd+=( --yes )
            fi

            echo \
                "Environment '${env_nam}' already exists; updating it because" \
                "'--if_exists update' was specified."

            return 0
            ;;
    esac
}


function warn_install_duration() {
    case "${env_nam}" in
        env_align|env_analyze)
            echo_warn \
                "creating '${env_nam}' may take some time given the large" \
                "number of packages to install. Do not worry if little or no" \
                "apparent progress is shown after 10, 20, or even 30" \
                "minutes. The installation should eventually complete."
        ;;
    esac
}


function run_install() {
    if [[ "${env_action}" == "update" ]]; then
        echo "Updating environment '${env_nam}'."
    else
        echo "Creating environment '${env_nam}'."
        warn_install_duration
    fi

    #MAYBE: change function from "private" to "public"
    _handle_env_deactivate

    if [[ "${env_action}" == "update" && "${yes}" == "true" ]]; then
        CONDA_ALWAYS_YES=true "${cmd[@]}" "${packages[@]}" || {
            echo_err \
                "failed to update environment '${env_nam}'. Please check the" \
                "error message(s) above."
            return 1
        }
    elif ! "${cmd[@]}" "${packages[@]}"; then
        echo_err \
            "failed to ${env_action} environment '${env_nam}'. Please check" \
            "the error message(s) above."
        return 1
    fi
}


function install_pkg_editable() {
    (( ${#cmd_pkg[@]} > 0 )) || return 0

    echo "Refreshing editable package in '${env_nam}'."
    mkdir -p "${dir_rep}/artifacts/tests/package"
    if ! "${cmd_pkg[@]}"; then
        echo_err \
            "failed to install the editable package in '${env_nam}'."
        return 1
    fi
}


function main() {
    local rc=0

    source_helpers_script
    init_arg_defs

    parse_args "$@" || rc=$?
    if (( rc == 2 )); then
        return 0
    elif (( rc != 0 )); then
        return "${rc}"
    fi

    validate_args
    resolve_env_definition
    build_install_cmd

    handle_existing_env || rc=$?
    if (( rc == 2 )); then
        return 0
    elif (( rc != 0 )); then
        return "${rc}"
    fi

    if [[ "${dry_run}" == "true" ]]; then
        print_dry_run
        return 0
    fi

    if [[ "${env_action}" == "create" || "${env_action}" == "update" ]]; then
        run_install || return 1
    fi

    install_pkg_editable
}


main "$@"
