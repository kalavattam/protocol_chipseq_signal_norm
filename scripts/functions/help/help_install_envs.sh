#!/bin/bash
# -*- coding: utf-8 -*-
# 
# Copyright 2024-2025 by Kris Alavattam
# Email: kalavattam@gmail.com
# 
# Distributed under the MIT license.
#
# Script: help_install_envs.sh
# Description: Defines the help text printer for install_envs.sh.


function help_install_envs() {
cat << EOM
install_envs.sh
    --env_nam <str> [--yes]

Description:
    Sets up one of multiple Mamba environments with necessary programs and
    dependencies used in this project. For more details, see "Notes" below.

Arguments:
     -h, --help     Display this help message and exit.
    -en, --env_nam  Mamba environment to create: 'env_analyze', 'env_protocol',
                    or 'env_siqchip'.
     -y, --yes      Automatically answer yes to mamba prompts.

Dependencies:
    - Bash
    - Mamba

Notes:
    - The call to 'mamba' will be adapted to allow Rosetta translation if the
      script detects that the system has "arm64" architecture.
    - The following packages are installed via a call to 'mamba create':
        + env_analyze
            - bioconductor-annotationdbi
            - bioconductor-chipqc
            - bioconductor-chipseeker
            - bioconductor-clusterprofiler
            - bioconductor-deseq2
            - bioconductor-diffbind
            - bioconductor-edger
            - bioconductor-enhancedvolcano
            - bioconductor-genomicfeatures
            - bioconductor-genomicranges
            - bioconductor-ihw
            - bioconductor-iranges
            - bioconductor-pcatools
            - bioconductor-sva
            - datamash
            - deeptools
            - ipython
            - pandas
            - parallel
            - pbzip2
            - phantompeakqualtools
            - pigz
            - r-argparse
            - r-dendextend
            - r-devtools
            - r-ggalt
            - r-ggpubr
            - r-ggrepel
            - r-ggsci
            - r-pheatmap
            - r-plotly
            - r-readxl
            - r-rjson
            - r-tidyverse
            - r-upsetr
            - r-venneuler
            - r-writexl
            - r-xml2
            - rename
            - tree
        + env_protocol
            - bc
            - bowtie2
            - datamash  ## NOTE: Added since publication in Bio-protocol ##
            - fastqc
            - gawk
            - ipython
            - matplotlib
            - multiqc
            - parallel
            - pbzip2
            - pigz
            - pysam
            - python=3.11  # Restrict to version 3.11 for 'sequali'
            - r-argparse
            - r-ggsci
            - r-plotly
            - rename
            - samtools
            - sequali
            - tree  ## NOTE: Added since publication in Bio-protocol ##
            - wget
        + env_siqchip
            - bc
            - bedtools
            - datamash
            - gfortran
            - gnuplot
            - parallel
            - samtools
            - tree
            - ucsc-bedclip
            - ucsc-bedgraphtobigwig
    - Given that >100 packages and dependencies may be installed depending on
      the environment specified, the execution and completion of the 'mamba
      create' operation will take >10 minutes if making use of cached packages
      from previous installations/separate environments. If this is a fresh
      installation that does not make use of cached packages, then the creation
      of the environment may take even longer: e.g., more than 20 or 30
      minutes.

Example:
    \`\`\`bash
    bash install_envs.sh
        --env_nam "env_protocol"
        --yes
    \`\`\`
EOM
}