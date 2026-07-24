#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: help_install_envs.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
# GPT-5.6) were used in development and documentation.
#
# Distributed under the MIT license.


# TODO: decide whether 'gnuplot' should remain in the 'env_protocol' definition
function help_install_envs() {
    cat << EOM
Usage
-----
  install/scripts/install_envs.sh
    [--help] [--dry_run] --env_nam <choice> [--if_exists <choice>] [--update_package <spec>] [--channels <csv>] [--override_channels] [--yes]

    Use Mamba (preferred) or Conda to create repository environments from YAML files under 'install/envs/'. See "Notes" for more details.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  -dr, --dry, --dry_run : flag
    Run script in dry-run mode. Print resolved installation command and exit without installing.

  -en, --env, --env_nam : {'env_analyze', 'env_protocol', 'env_siqchip'}
    Conda environment to activate. Environment to create: 'env_analyze', 'env_protocol', or 'env_siqchip'.

  -ie, --if_exists : {'fail', 'reuse', 'update'}
    What to do if the requested environment already exists: fail, reuse without dependency reconciliation, or update from its YAML (default: 'fail').

  -up, --update_package : str
    With '--if_exists update', install only this exact YAML-declared package specification. Repeat to select more than one package; omit to reconcile every declared dependency.

  -ch, --channels : list of str
    Comma-delimited package-manager channels to append to the environment creation command.

  -oc, --override_channels : flag
    Add '--override-channels' to the package-manager command.

  -y, --yes : flag
    Automatically answer yes to package-manager prompts.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - mamba or conda

  - User-facing environments are created from YAML files under 'install/envs/' via 'mamba env create -f <yaml>' (preferred) or 'conda env create -f <yaml>'.
  - See the following YAML files for environment package lists:

      | environment  | file path                     |
      | :----        | :----                         |
      | env_analyze  | install/envs/env_analyze.yml  |
      | env_protocol | install/envs/env_protocol.yml |
      | env_siqchip  | install/envs/env_siqchip.yml  |

  - Depending on the specified environment, a large number of packages and dependencies may need to be installed. As a result, the 'mamba env create' or 'conda env create' operation can take more than 10 minutes, especially for a fresh installation that does not make use of cached packages. In such cases, environment creation may take even longer (e.g., more than 20 or 30 minutes).
  - With '--if_exists reuse', this script skips environment creation. For 'env_protocol', it refreshes the repository's managed editable Python-package installation; it does not otherwise reconcile environment dependencies.
  - With '--if_exists update', this script reconciles a YAML-backed existing environment with installed packages frozen and without pruning, then refreshes the managed editable package when applicable. YAML channels remain authoritative. Use repeatable '--update_package' selections for a reviewed incremental transaction.
  - Use '--channels' and '--override_channels' to apply site-specific channels, such as institutional mirrors, in place of or in addition to channels listed in the YAML files.

Examples
--------
  1. Preview creation of the main workflow environment with automatic confirmation.
    '''bash
    bash install/scripts/install_envs.sh \\
        --dry_run \\
        --env_nam env_protocol \\
        --yes
    '''

  2. Preview creation through site-specific channels that replace YAML channels.
    '''bash
    bash install/scripts/install_envs.sh \\
        --dry_run \\
        --env_nam env_protocol \\
        --channels fhcc-main,fhcc-bioconda \\
        --override_channels \\
        --yes
    '''
EOM
}


#  Hidden environments:
#+
#+ + env_align  ## NOTE: Retained for old work; not exposed in the docs ##
#+   - bamtools
#+   - bbmap
#+   - bedtools
#+   - bowtie2
#+   - bwa
#+   - datamash
#+   - fastqc
#+   - gawk
#+   - gnuplot
#+   - macs3
#+   - minimap
#+   - mosdepth
#+   - parallel
#+   - picard
#+   - preseq
#+   - rename
#+   - samtools
#+   - subread
#+   - tree
#+   - ucsc-bedgraphtobigwig
#+   - ucsc-bedsort
#+   - ucsc-facount
#+   - wget
#+ + env_repro  ## NOTE: Not exposing this to users in the docs ##
#+   - bc
#+   - bowtie2=2.3.4.2  ## NOTE: Explicitly pinning old version ##
#+   - deeptools=3.3.1  ## NOTE: Explicitly pinning old version ##
#+   - gawk
#+   - ipython
#+   - parallel
#+   - pbzip2
#+   - pigz
#+   - python=3.6       ## NOTE: Explicitly pinning old version ##
#+   - rename
#+   - samtools=1.9     ## NOTE: Explicitly pinning old version ##
#+   - tree
#+   - wget
