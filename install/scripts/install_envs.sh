#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: install_envs.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
#   GPT-5.6);
# - Anthropic Claude Code (Opus 5).
#
# Distributed under the MIT license.


# Require Bash >= 4.4 before doing any work.
if [[ -z "${BASH_VERSION:-}" ]]; then
    echo "error(shell):" \
        "this script must be run under Bash >= 4.4." >&2
    exit 1
elif ((
    BASH_VERSINFO[0] < 4 || (BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] < 4)
)); then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "this script requires Bash >= 4.4; current version is" \
        "'${BASH_VERSION}'." >&2
    exit 1
fi

# Run in safe mode, exiting on errors, unset variables, and pipe failures.
set -euo pipefail

# Set paths to installation support and repository directories.
dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
dir_ins="$(cd "${dir_scr}/.." > /dev/null 2>&1 && pwd)"
dir_rep="$(cd "${dir_ins}/.." > /dev/null 2>&1 && pwd)"


# Source shared helpers.
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
    if_exists_explicit=false
    channels=""
    override_channels=false
    yes=false
    pth_yml=""
    pth_yml_eff=""
    pth_condarc=""
    dir_tmp=""
    env_action=create
    pkg_mgr=""
    pkg_mgr_v=""

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
                if_exists_explicit=true
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

    # '--update_package' is meaningful only when reconciling an existing
    # environment, so it implies '--if_exists update' rather than requiring the
    # caller to name both. Refuse only a genuine conflict, where the caller
    # asked for some other behaviour explicitly.
    if (( ${#requested_packages[@]} > 0 )); then
        if [[ "${if_exists_explicit}" != "true" ]]; then
            if_exists=update
        elif [[ "${if_exists}" != "update" ]]; then
            echo_err \
                "'--update_package' conflicts with '--if_exists" \
                "${if_exists}'. Package selections apply only to" \
                "'--if_exists update'."
            return 1
        fi
    fi

    case "${env_nam}" in
        env_align|env_analyze|env_protocol|env_repro|env_siqchip) : ;;
        *)
            # Note: 'env_align' and 'env_repro' are not exposed to users.
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

    # Channels supplied with '--channels' are added ahead of the declared
    # ones, matching what '-c' means to conda itself: higher priority, not
    # replacement. '--override_channels' is what makes the supplied list
    # exclusive, so only then are the declared channels dropped.
    if [[ "${override_channels}" != "true" ]]; then
        arr_channels+=( "${yaml_channels[@]}" )
    fi

    packages=( "${yaml_packages[@]}" )
}


function ensure_tmp_dir() {
    if [[ -n "${dir_tmp}" && -d "${dir_tmp}" ]]; then
        return 0
    fi

    # Create a directory and name files inside it, rather than using a 'mktemp'
    # template. BSD and GNU 'mktemp' disagree about templates whose 'XXXXXX' is
    # not final, and conda infers a file's format from its extension, so the
    # names have to survive intact.
    dir_tmp="$(mktemp -d)" || {
        echo_err \
            "failed to create a temporary directory for rendered" \
            "installation files."
        return 1
    }
}


function render_condarc() {
    if [[ -n "${pth_condarc}" ]]; then
        return 0
    fi

    if (( ${#arr_channels[@]} == 0 )); then
        return 0
    fi

    ensure_tmp_dir || return 1

    pth_condarc="${dir_tmp}/condarc"

    # 'mirrored_channels' maps a channel name onto a mirror list, so a supplied
    # URL whose final path segment matches that name is fetched from the list
    # instead. Miniforge ships one sending 'conda-forge' to 'anaconda.org',
    # unreachable at a proxying site. No channel flag reaches it.
    #
    # Written unconditionally: conda 24.7.1 and mamba 1.5.9 do not know the
    # setting and accept a file carrying it without complaint.
    #
    # 'CONDARC' replaces the caller's configuration rather than merging, so
    # this file also drops the channels that configuration declares. Measured
    # 2026-08-13, and it is what makes creation exclusive; adding a 'channels:'
    # key here would silently change that.
    printf 'mirrored_channels: {}\n' > "${pth_condarc}" || {
        echo_err "failed to render a package-manager configuration file."
        return 1
    }
}


function render_env_yaml() {
    local channel=""
    local in_channels=false
    local line=""

    pth_yml_eff="${pth_yml}"

    if [[ -z "${pth_yml}" ]] || (( ${#arr_channels[@]} == 0 )); then
        return 0
    fi

    # 'conda env create' and 'mamba env create' reject '--override-channels'
    # and '-c' on the supported baseline — confirmed on conda 24.7.1 and mamba
    # 1.5.9 — so channel selection cannot be expressed as flags on that
    # subcommand. Mamba 2.x does accept them, but relying on that would make
    # the supported versions diverge for no gain. Rewrite the declared
    # 'channels:' block instead and install from the rewritten copy, which
    # every version reads identically.
    ensure_tmp_dir || return 1

    pth_yml_eff="${dir_tmp}/${env_nam}.yml"

    while IFS= read -r line || [[ -n "${line}" ]]; do
        if [[ "${line}" =~ ^channels:[[:space:]]*$ ]]; then
            printf 'channels:\n'

            for channel in "${arr_channels[@]}"; do
                printf '  - %s\n' "${channel}"
            done

            # 'nodefaults' is the environment-file spelling of
            # '--override-channels'; without it conda still consults the
            # 'defaults' channel named in a user's '.condarc'.
            if [[ "${override_channels}" == "true" ]]; then
                printf '  - nodefaults\n'
            fi

            in_channels=true
            continue
        fi

        if [[ "${in_channels}" == "true" ]]; then
            if [[ "${line}" =~ ^[[:space:]]+-[[:space:]] ]]; then
                # Declared channels survive below the supplied ones unless the
                # caller asked for the supplied list to be exclusive.
                if [[ "${override_channels}" == "true" ]]; then
                    continue
                fi
            elif [[ "${line}" =~ ^[^[:space:]] ]]; then
                in_channels=false
            fi
        fi

        printf '%s\n' "${line}"
    done < "${pth_yml}" > "${pth_yml_eff}" || {
        echo_err "failed to render '${pth_yml}' to '${pth_yml_eff}'."
        return 1
    }
}


function cleanup_rendered() {
    if [[ -n "${dir_tmp}" && -d "${dir_tmp}" ]]; then
        rm -rf "${dir_tmp}"
    fi
}


function check_pkg_mgr_v() {
    local floor_major=""
    local floor_minor=""
    local version=""
    local major=""
    local minor=""

    # Take the first line. Mamba 1.5.x prints its own version and then the
    # bundled conda version on a second line, so reading the whole output
    # yields two numbers and comparing the wrong one is easy.
    version="$(
        "${pkg_mgr}" --version 2>/dev/null \
            | head -n 1 \
            | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?' \
            | head -n 1
    )"

    if [[ -z "${version}" ]]; then
        echo_warn \
            "could not determine the '${pkg_mgr}' version; proceeding" \
            "without a version check."
        return 0
    fi

    pkg_mgr_v="${version}"

    # The floors are the versions published with the Bio-protocol manuscript.
    # Anything at or above them is supported; nothing below is tested.
    if [[ "${pkg_mgr}" == "mamba" ]]; then
        floor_major=1
        floor_minor=5
    else
        floor_major=24
        floor_minor=7
    fi

    major="${version%%.*}"
    minor="${version#*.}"
    minor="${minor%%.*}"

    if \
           (( major < floor_major )) \
        || { (( major == floor_major )) && (( minor < floor_minor )); }
    then
        echo_err \
            "'${pkg_mgr}' ${version} is older than the supported minimum," \
            "${floor_major}.${floor_minor}."
        echo >&2
        echo \
            "This repository is tested against 'conda' >= 24.7 and 'mamba'" \
            ">= 1.5, the versions published with the Bio-protocol" \
            "manuscript. Please update before installing." >&2
        return 1
    fi
}


function resolve_pkg_mgr() {
    # shellcheck disable=SC2119
    check_pkg_mgr || return 1

    if command -v mamba >/dev/null 2>&1; then
        pkg_mgr=mamba
    else
        pkg_mgr=conda
    fi

    check_pkg_mgr_v || return 1
}


function init_cmd_with_condarc() {
    cmd=()

    if [[ -n "${pth_condarc}" ]]; then
        # Set for this command alone. The caller's own configuration is left
        # untouched, and the assignment is part of the command the dry run
        # prints, so what is reported is what runs.
        cmd=( env "CONDARC=${pth_condarc}" )
    fi
}


function build_create_cmd() {
    local channel=""

    validate_var "pkg_mgr" "${pkg_mgr}" || return 1

    init_cmd_with_condarc

    if [[ -n "${pth_yml}" ]]; then
        cmd+=( "${pkg_mgr}" env create -f "${pth_yml_eff}" )
    else
        cmd+=( "${pkg_mgr}" create -n "${env_nam}" )
    fi

    # Channel flags belong only to the package-list form. On the YAML form the
    # channels are already rendered into the file that '-f' points at, and the
    # 'env create' subcommand rejects these flags outright.
    if [[ -z "${pth_yml}" ]]; then
        if [[ "${override_channels}" == "true" ]]; then
            cmd+=( --override-channels )
        fi

        for channel in "${arr_channels[@]}"; do
            cmd+=( -c "${channel}" )
        done
    fi

    if [[ "${yes}" == "true" ]]; then
        cmd+=( --yes )
    fi

    if [[ -z "${pth_yml}" && "${env_nam}" == "env_align" ]]; then
        packages=(  # Note: retained for old work; not exposed in the docs.
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
        packages=(  # Note: Not exposing this to users in the docs.
            bc
            bowtie2=2.3.4.2  # Note: explicitly pinning old version.
            deeptools=3.3.1  # Note: explicitly pinning old version.
            gawk
            ipython
            parallel
            pbzip2
            pigz
            python=3.6       # Note: explicitly pinning old version.
            rename
            samtools=1.9     # Note: explicitly pinning old version.
            tree
            wget
        )
    fi

}


function build_pkg_cmd() {
    validate_var "pkg_mgr" "${pkg_mgr}" || return 1

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

    if [[ -n "${pkg_mgr}" ]]; then
        echo "Package manager: ${pkg_mgr} ${pkg_mgr_v:-unknown}"
    fi

    if [[ -n "${pth_yml}" ]]; then
        echo "YAML: ${pth_yml}"
    fi

    # Report the rendered file only when a command would install from it. On
    # the stop path no install happens and the file is removed on the way out,
    # so naming it here would print a path that no longer exists.
    if [[
        "${env_action}" != "stop" \
        && -n "${pth_yml_eff}" \
        && "${pth_yml_eff}" != "${pth_yml}"
    ]]; then
        echo "Rendered YAML: ${pth_yml_eff}"
        echo "Rendered channels:"

        # Read these back from the rendered file rather than reprinting the
        # supplied list. Without '--override_channels' the declared channels
        # are retained below the supplied ones, so the supplied list alone
        # would understate what the solver is actually given.
        sed -n '/^channels:$/,/^[^[:space:]]/p' "${pth_yml_eff}" \
            | grep -E '^[[:space:]]+-[[:space:]]' \
            || true
    fi

    # Report the rendered configuration wherever a command would use it. It is
    # what stops a mirrored channel name from redirecting the supplied
    # channels, so a reader checking why an install reached a given host needs
    # to see that it was in force.
    if [[ "${env_action}" != "stop" && -n "${pth_condarc}" ]]; then
        echo "Rendered condarc: ${pth_condarc}"
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

            # No '--freeze-installed' here. Freezing means "add what is missing
            # and change nothing else", which cannot reconcile an environment
            # against a YAML whose declared version differs from the installed
            # one — the solver reports a conflict instead. Scope is narrowed
            # with '--update_package', which is explicit about what may change,
            # rather than by refusing changes.
            render_condarc || return 1

            init_cmd_with_condarc
            cmd+=( "${pkg_mgr}" install -n "${env_nam}" )

            # '--override-channels' is conda's spelling of this script's
            # '--override_channels', so it is passed only when the caller asked
            # for it. Passing it unconditionally silently applied the flag to
            # every update, which both contradicted the interface and made the
            # update path disagree with creation.
            if [[ "${override_channels}" == "true" ]]; then
                cmd+=( --override-channels )
            fi

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

    # Maybe: change function from "private" to "public".
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

    validate_args          || return 1
    resolve_env_definition || return 1
    resolve_pkg_mgr        || return 1
    build_pkg_cmd          || return 1

    # Decide what will happen before building anything for it. The update
    # branch constructs its own command, and the reuse and stop branches
    # construct none, so rendering or building ahead of this point would do
    # work for actions that never run.
    handle_existing_env || rc=$?
    if (( rc == 2 )); then
        return 0
    elif (( rc != 0 )); then
        return "${rc}"
    fi

    if [[ "${env_action}" == "create" ]]; then
        render_condarc   || return 1
        render_env_yaml  || return 1
        build_create_cmd || return 1
    fi

    # A dry run retains the rendered file so the caller can read exactly what
    # would have been installed from.
    if [[ "${dry_run}" == "true" ]]; then
        print_dry_run
        return 0
    fi

    if [[ "${env_action}" == "create" || "${env_action}" == "update" ]]; then
        if ! run_install; then
            # Retain the rendered files on failure; they are the record of
            # which channels were actually offered to the solver, and of
            # whether channel redirection was suppressed.
            if [[ -n "${dir_tmp}" && -d "${dir_tmp}" ]]; then
                if [[
                    -n "${pth_yml_eff}" \
                    && "${pth_yml_eff}" != "${pth_yml}"
                ]]; then
                    echo_err \
                        "the rendered environment YAML was kept for" \
                        "inspection: '${pth_yml_eff}'."
                fi

                if [[ -n "${pth_condarc}" ]]; then
                    echo_err \
                        "the rendered package-manager configuration was kept" \
                        "for inspection: '${pth_condarc}'."
                fi
            fi

            return 1
        fi
    fi

    cleanup_rendered
    install_pkg_editable
}


main "$@"
