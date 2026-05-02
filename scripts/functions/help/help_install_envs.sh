#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: help_install_envs.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
#
# Distributed under the MIT license.


function help_install_envs() {
    cat << EOM
Usage:
    install_envs.sh
        [--help] [--dry_run] --env_nam {env_analyze,env_protocol,env_siqchip} [--yes]

Description:
    Use Mamba (preferred) or Conda to set up one of a few environments containing the programs and dependencies used in this project. See "Notes" for more details.

Arguments:
     -h, --help     Display this help message and exit.
    -dr, --dry_run  Print resolved installation command and exit without installing.
    -en, --env_nam  Package manager environment to create: 'env_analyze', 'env_protocol', or 'env_siqchip'.
     -y, --yes      Automatically answer yes to package-manager prompts.

Dependencies:
    - Bash >= 5
    - Conda
    - Mamba (optional, preferred)

Notes:
    - The following packages are installed via a call to 'mamba create' (preferred) or 'conda create':
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
            - gawk
            - gnuplot
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
            - asciigenome  ## NOTE: Added since publication in Bio-protocol ##
            - bash         ## NOTE: Added since publication in Bio-protocol ##
            - bc
            - bowtie2
            - bwa          ## NOTE: Added since publication in Bio-protocol ##
            - bwa-mem2     ## NOTE: Added since publication in Bio-protocol ##
            - datamash     ## NOTE: Added since publication in Bio-protocol ##
            - fastqc
            - gawk
            - gnuplot      ## NOTE: Added since publication in Bio-protocol ##  ## TODO: remove from 'env_protocol' ##
            - ipython
            - matplotlib
            - multiqc
            - parallel
            - pbzip2
            - pigz
            - pysam
            - python=3.11  ## NOTE: Restrict to v3.11 for 'sequali' installation ##
            - pyyaml       ## NOTE: Made explicit since Bio-protocol publication ##
            - r-argparse
            - r-ggsci
            - r-plotly
            - rename
            - samtools
            - sequali
            - tree         ## NOTE: Added since publication in Bio-protocol ##
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
    - Depending on the specified environment, a large number of packages and dependencies may need to be installed. As a result, the 'mamba create' or 'conda create' operation can take more than 10 minutes, especially for a fresh installation that does not make use of cached packages. In such cases, environment creation may take even longer (e.g., more than 20 or 30 minutes).

Examples:
    1. Dry-run environment installation for the main workflow, including '--yes' in the resolved package-manager command
    '''bash
    bash install_envs.sh
        --dry_run
        --env_nam "env_protocol"
        --yes
    '''

    2. Run environment installation for the main workflow, using '--yes' to automatically confirm package-manager prompts
    '''bash
    bash install_envs.sh
        --env_nam "env_protocol"
        --yes
    '''

    3. Install environment needed to run siQ-ChIP (see github.com/BradleyDickson/siQ-ChIP or github.com/kalavattam/siQ-ChIP/tree/protocol)
    '''bash
    bash install_envs.sh --env_nam "env_siqchip"
    '''
EOM
}


#  Hidden environments:
#+
#+ + env_align  ## NOTE: Retained for old work; not exposed in the docs ##
#+     - bamtools
#+     - bbmap
#+     - bedtools
#+     - bowtie2
#+     - bwa
#+     - datamash
#+     - fastqc
#+     - gawk
#+     - gnuplot
#+     - macs3
#+     - minimap
#+     - mosdepth
#+     - parallel
#+     - picard
#+     - preseq
#+     - rename
#+     - samtools
#+     - subread
#+     - tree
#+     - ucsc-bedgraphtobigwig
#+     - ucsc-bedsort
#+     - ucsc-facount
#+     - wget
#+ + env_repro  ## NOTE: Not exposing this to users in the docs ##
#+     - bc
#+     - bowtie2=2.3.4.2  ## NOTE: Explicitly pinning old version ##
#+     - deeptools=3.3.1  ## NOTE: Explicitly pinning old version ##
#+     - gawk
#+     - ipython
#+     - parallel
#+     - pbzip2
#+     - pigz
#+     - python=3.6       ## NOTE: Explicitly pinning old version ##
#+     - rename
#+     - samtools=1.9     ## NOTE: Explicitly pinning old version ##
#+     - tree
#+     - wget
