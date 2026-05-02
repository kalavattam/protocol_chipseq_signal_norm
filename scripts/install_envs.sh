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


#  Require Bash >= 5 before doing any work
if [[ -z "${BASH_VERSION:-}" ]]; then
    echo "error(shell):" \
        "this script must be run under Bash >= 5." >&2
    exit 1
elif (( BASH_VERSINFO[0] < 5 )); then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "this script requires Bash >= 5; current version is" \
        "'${BASH_VERSION}'." >&2
    exit 1
fi

#  Run in safe mode, exiting on errors, unset variables, and pipe failures
set -euo pipefail

#  Set path to the 'scripts' directory
dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"


#  Source and define functions ================================================
# arr_fnc=(  #TODO: record in 'help_install_envs' before deleting
#     check_env
#     format_outputs
#     handle_env
#     help/help_install_envs
#     # check_env_installed  ## check_env ##
#     # echo_err           ## format_outputs ##
#     # echo_warn         ## format_outputs ##
#     # handle_env           ## handle_env ##
#     # help/help_install_envs
# )

dir_fnc="${dir_scr}/functions"
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
  Bash >= 5
EOM
    )

    if [[ "${arg}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -n "${arg}" ]]; then
        echo \
            "Error: Unexpected argument to 'check_pkg_mgr()':" \
            "'${arg}'." >&2
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    if command -v mamba >/dev/null 2>&1; then
        return 0
    elif command -v conda >/dev/null 2>&1; then
        return 0
    else
        echo "Error: Neither Mamba nor Conda is installed on the system." >&2
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

if \
    check_env_installed "${env_nam}" "true"
then
    echo_err \
        "an environment with the name '${env_nam}' is already installed."
    exit 1
fi

#  Check that supported package manager is in PATH
# shellcheck disable=SC2119
check_pkg_mgr || exit 1


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

#  Construct the package manager command
declare -a cmd packages

if command -v mamba >/dev/null 2>&1; then
    cmd=( mamba create -n "${env_nam}" )
else
    cmd=( conda create -n "${env_nam}" )
fi

if [[ "${yes}" == "true" ]]; then
    cmd+=( --yes )
fi

#  Assign an array of packages to install
#TODO: switch to external YAML files for env package lists
if [[ "${env_nam}" == "env_align" ]]; then
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
elif [[ "${env_nam}" == "env_analyze" ]]; then
    packages=(
        bioconductor-annotationdbi
        bioconductor-chipqc
        bioconductor-chipseeker
        bioconductor-clusterprofiler
        bioconductor-deseq2
        bioconductor-diffbind
        bioconductor-edger
        bioconductor-enhancedvolcano
        bioconductor-genomicfeatures
        bioconductor-genomicranges
        bioconductor-ihw
        bioconductor-iranges
        bioconductor-pcatools
        bioconductor-sva
        datamash
        deeptools
        gawk
        gnuplot
        ipython
        pandas
        parallel
        pbzip2
        phantompeakqualtools
        pigz
        r-argparse
        r-dendextend
        r-devtools
        r-ggalt
        r-ggpubr
        r-ggrepel
        r-ggsci
        r-pheatmap
        r-plotly
        r-readxl
        r-rjson
        r-tidyverse
        r-upsetr
        r-venneuler
        r-writexl
        r-xml2
        rename
        tree
    )
elif [[ "${env_nam}" == "env_protocol" ]]; then
    packages=(
        asciigenome  ## NOTE: Added since publication in Bio-protocol ##
        bash         ## NOTE: Added since publication in Bio-protocol ##
        bc
        bowtie2
        bwa          ## NOTE: Added since publication in Bio-protocol ##
        bwa-mem2     ## NOTE: Added since publication in Bio-protocol ##
        datamash     ## NOTE: Added since publication in Bio-protocol ##
        fastqc
        gawk
        gnuplot      ## NOTE: Added since publication in Bio-protocol ## ## TODO: remove from 'env_protocol' ##
        ipython
        matplotlib
        multiqc
        parallel
        pbzip2
        pigz
        pysam
        python=3.11  ## NOTE: Restrict to v3.11 for 'sequali' installation ##
        pyyaml       ## NOTE: Made explicit since Bio-protocol publication ##
        r-argparse
        r-ggsci
        r-plotly
        rename
        samtools
        sequali
        tree         ## NOTE: Added since publication in Bio-protocol ##
        wget
    )
elif [[ "${env_nam}" == "env_repro" ]]; then
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
elif [[ "${env_nam}" == "env_siqchip" ]]; then
    packages=(
        bc
        bedtools
        datamash
        gfortran
        gnuplot
        parallel
        samtools
        tree
        ucsc-bedclip
        ucsc-bedgraphtobigwig
    )
fi

#  In dry-run mode, print resolved command and exit without installing
if [[ "${dry_run}" == "true" ]]; then
    echo "Dry run: would create environment '${env_nam}'."

    printf 'Command:'
    for tok in "${cmd[@]}" "${packages[@]}"; do
        printf ' %q' "${tok}"
    done
    printf '\n'

    exit 0
fi

#  Run the environment installation
if ! "${cmd[@]}" "${packages[@]}"; then
    echo_err \
        "failed to create environment '${env_nam}'. Please check the error" \
        "message(s) above."
    exit 1
fi
