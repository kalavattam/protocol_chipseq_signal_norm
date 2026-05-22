#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: install_envs.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
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


#  Source and define functions ================================================
dir_fnc="${dir_rep}/scripts/functions"
fnc_src="${dir_fnc}/source_helpers.sh"

if [[ ! -f "${fnc_src}" ]]; then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "script not found: '${fnc_src}'." >&2
    exit 1
fi

# shellcheck disable=SC1090
source "${fnc_src}" || {
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "failed to source '${fnc_src}'." >&2
    exit 1
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
        exit 1
    }

unset fnc_src


# shellcheck disable=SC2120
function check_pkg_mgr() {
    local arg="${1:-}"
    local show_help

    show_help=$(cat << EOM
Usage:
  check_pkg_mgr [-h|--hlp|--help]

Description:
  Check that either Mamba or Conda is installed and available in PATH.

Returns:
  0 if Mamba or Conda is available; otherwise, 1 and an error message.

Dependency:
  Bash >= 4.4
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


#  Parse arguments ============================================================
#  Initialize variables along with default assignments
dry_run=false
env_nam=""
if_exis="fail"
channels=""
override_channels=false
yes=false

#  Parse arguments
if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
    help_install_envs >&2
    exit 0
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
                exit 1
            }
            env_nam="${2}"
            shift 2
            ;;

        -ie|--if[_-]ex|--if[_-]exis|--if[_-]exists)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_install_envs >&2
                exit 1
            }
            if_exis="${2}"
            shift 2
            ;;

        -ch|--channel|--channels|--channel[_-]list)
            require_optarg "${1}" "${2:-}" "main" || {
                echo >&2
                help_install_envs >&2
                exit 1
            }
            channels="${2}"
            shift 2
            ;;

        -oc|--override[_-]channel|--override[_-]channels)
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
            exit 1
            ;;
    esac
done


#  Check that required arguments are provided, appropriate, formatted, etc.
validate_var "env_nam" "${env_nam}"
validate_var "if_exis" "${if_exis}"

declare -a arr_channels
if [[ -n "${channels}" ]]; then
    IFS=',' read -r -a arr_channels <<< "${channels}"

    for channel in "${arr_channels[@]}"; do
        if [[ -z "${channel}" ]]; then
            echo_err \
                "invalid '--channels' value: '${channels}'. Channel names" \
                "must be comma-delimited and non-empty."
            exit 1
        fi
    done
else
    arr_channels=()
fi

if [[ "${override_channels}" == "true" && "${#arr_channels[@]}" -eq 0 ]]; then
    echo_err "'--override_channels' requires '--channels'."
    exit 1
fi

case "${if_exis}" in
    fail|reuse) : ;;
    *)
        echo_err \
            "invalid '--if_exis' value: '${if_exis}'. Must be 'fail' or" \
            "'reuse'."
        exit 1
        ;;
esac

case "${env_nam}" in
    env_align|env_analyze|env_protocol|env_repro|env_siqchip) : ;;
    *)
        ## NOTE: 'env_align' and 'env_repro' are not exposed to users ##
        echo_err \
            "invalid environment name specified. Must be 'env_analyze'," \
            "'env_protocol', or 'env_siqchip'."
        exit 1
        ;;
esac

#  Resolve environment definition and package list ============================
pth_yml=""

case "${env_nam}" in
    env_analyze|env_protocol|env_siqchip)
        pth_yml="${dir_rep}/install/envs/${env_nam}.yml"

        if [[ ! -f "${pth_yml}" ]]; then
            echo_err "environment YAML does not exist: '${pth_yml}'."
            exit 1
        fi

        if [[ ! -r "${pth_yml}" ]]; then
            echo_err "environment YAML is not readable: '${pth_yml}'."
            exit 1
        fi
        ;;
esac

#  Construct the package manager command =====================================
declare -a cmd packages

#  Check that supported package manager is in PATH
# shellcheck disable=SC2119
check_pkg_mgr || exit 1

if command -v mamba >/dev/null 2>&1; then
    if [[ -n "${pth_yml}" ]]; then
        cmd=( mamba env create -f "${pth_yml}" )
    else
        cmd=( mamba create -n "${env_nam}" )
    fi
else
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

function print_dry_run() {
    echo "dryrun($(basename "${BASH_SOURCE[0]}")):" \
        "would create environment '${env_nam}'."

    if [[ -n "${pth_yml}" ]]; then
        echo "YAML: ${pth_yml}"
    fi

    printf 'Command:'
    for tok in "${cmd[@]}" "${packages[@]}"; do
        printf ' %q' "${tok}"
    done
    printf '\n'
}

if \
    check_env_installed "${env_nam}" "true"
then
    case "${if_exis}" in
        fail)
            echo_err \
                "an environment with the name '${env_nam}' is already" \
                "installed."
            echo >&2

            if [[ "${dry_run}" == "true" ]]; then
                print_dry_run
                echo >&2
                echo "dryrun($(basename "${BASH_SOURCE[0]}")):" \
                    "a non-dry run would stop here unless '--if_exis reuse'" \
                    "was specified." >&2
                echo >&2
                exit 0
            fi

            echo \
                "Nothing was changed. To reuse the existing environment," \
                "rerun with '--if_exis reuse'." >&2
            echo >&2

            echo \
                "To rebuild '${env_nam}', remove the existing environment" \
                "first, then rerun this script." >&2
            echo >&2

            echo \
                "For example, run 'mamba env remove -n \"${env_nam}\"'" \
                "prior to rerunning the script." >&2
            exit 1
            ;;

        reuse)
            echo \
                "Environment '${env_nam}' already exists; reusing it because" \
                "'--if_exis reuse' was specified."
            exit 0
            ;;
    esac
fi

#  In dry-run mode, print resolved command and exit without installing
if [[ "${dry_run}" == "true" ]]; then
    print_dry_run
    exit 0
fi

#  Do the main work ===========================================================
echo "Creating environment '${env_nam}'."

#  Warn about potentially time-consuming installations
case "${env_nam}" in
    env_align|env_analyze)
        echo_warn \
            "creating '${env_nam}' may take some time given the large number" \
            "of packages to install. Do not worry if little or no apparent" \
            "progress is shown after 10, 20, or even 30 minutes. The" \
            "installation should eventually complete."
    ;;
esac

#  If not in base environment, deactivate current environment
_handle_env_deactivate  #MAYBE: change function from "private" to "public"

#  Run the environment installation
if ! "${cmd[@]}" "${packages[@]}"; then
    echo_err \
        "failed to create environment '${env_nam}'. Please check the error" \
        "message(s) above."
    exit 1
fi
