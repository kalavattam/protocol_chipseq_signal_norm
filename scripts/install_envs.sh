#!/bin/bash
# -*- coding: utf-8 -*-
# 
# Copyright 2024-2025 by Kris Alavattam
# Email: kalavattam@gmail.com
# 
# Distributed under the MIT license.
#
# Script: install_envs.sh
# Description: Sets up Mamba environments used in the workflow.


#  Run script in interactive mode (true) or command-line mode (false)
interactive=false

#  Exit on errors, unset variables, or pipe failures if not in "interactive
#+ mode"
if ! ${interactive:-false}; then set -euo pipefail; fi

#  Set the path to the "scripts" directory
if ${interactive:-false}; then
    ## WARNING: If 'interactive=true', change path as needed ##
    dir_scr="${HOME}/repos/protocol_chipseq_signal_norm/scripts"
else
    dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi


#  Source and define functions ================================================
#  Set the path to the "functions" directory
dir_fnc="${dir_scr}/functions"

# shellcheck disable=SC1090
for script in \
    check_installed_env.sh \
    check_installed_mamba.sh \
    echo_error.sh \
    echo_warning.sh \
    handle_env.sh \
    help/help_install_envs.sh
do
    source "${dir_fnc}/${script}"
done


#  Set values for interactive mode
function set_interactive() {
    #  Hardcoded argument assignments
    env_nam="env_align"
    yes=false
}


#  Parse arguments ============================================================
#  Initialize variables along with default assignments
env_nam=""
yes=false

#  Parse arguments
if [[ -z "${1:-}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    help_install_envs
    if ! ${interactive}; then exit 0; fi
fi

if ${interactive:-false}; then
    set_interactive
else
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
           -en|--env_nam) env_nam="${2}"; shift 2 ;;
            -y|--yes)     yes=true;       shift 1 ;;
            *)
                echo "## Unknown parameter passed: '${1}' ##" >&2
                echo "" >&2
                help_install_envs >&2
                exit 1
                ;;
        esac
    done
fi


#  Check that required arguments are provided, appropriate, formatted, etc.
if [[ -z "${env_nam}" ]]; then
    echo_error "Environment name was not specified."
fi

case "${env_nam}" in
    env_align|env_analyze|env_protocol|env_repro|env_siqchip) : ;;
    *)
        ## NOTE: 'env_align' and 'env_repro' are not exposed to users ##
        echo_error \
            "Invalid environment name specified. Must be 'env_analyze'," \
            "'env_protocol', or 'env_siqchip'."
        ;;
esac

if check_installed_env "${env_nam}"; then
    echo_error \
        "An environment with the name '${env_nam}' is already installed."
fi

#  Check that dependencies are in PATH
check_installed_mamba


#  Do the main work ===========================================================
echo "Creating environment '${env_nam}'."

case "${env_nam}" in
    env_align|env_analyze)
        echo_warning \
            "Creating '${env_nam}' will take some time given the >100" \
            "packages to install. Don't worry if this script is still" \
            "running without any apparent progress after 10, 20, or even 30" \
            "minutes. It will eventually complete."
    ;;
esac

#  If not in base environment, then deactivate current environment
handle_env_deactivate

#  Construct the mamba command
cmd="mamba create -n ${env_nam}"

if ${yes}; then cmd+=" --yes"; fi

#  Assign an array of packages to install
if [[ "${env_nam}" == "env_align" ]]; then
    packages=(  ## NOTE: For old work; not exposing this in the docs ##
        bamtools
        bbmap
        bedtools
        bowtie2
        bwa
        datamash
        fastqc
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
        bc
        bowtie2
        datamash  ## NOTE: Added since publication in Bio-protocol ##
        fastqc
        gawk
        ipython
        matplotlib
        multiqc
        parallel
        pbzip2
        pigz
        pysam
        python=3.11  # Restrict to version 3.11 for 'sequali' installation
        r-argparse
        r-ggsci
        r-plotly
        rename
        samtools
        sequali
        tree  ## NOTE: Added since publication in Bio-protocol ##
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

#  Run the mamba environment installation
if ! eval "${cmd} ${packages[*]}"; then
    echo_error \
        "Failed to create environment '${env_nam}'. Please check the error" \
        "message(s) above."
    exit 1
fi
