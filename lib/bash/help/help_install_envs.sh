#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: help_install_envs.sh
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


# TODO: decide whether 'gnuplot' should remain in 'env_protocol'. Probably not.
# TODO FIXME: under Usage, do a sensical breakdown of option types per line.
function help_install_envs() {
    cat >&2 << EOM
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
    What to do if the requested environment already exists (default: 'fail'). 'fail' stops without changing anything, and reports how to reuse or rebuild the environment instead. 'reuse' leaves the environment as it is, reconciling no dependencies. 'update' creates a missing YAML-backed environment or reconciles an existing one to its YAML: declared packages are installed, and an installed version is changed where the YAML declares a different one, which may mean a downgrade; packages the YAML no longer lists are left in place rather than pruned.

  -up, --update_package : str
    Reconcile only this exact YAML-declared package specification, rather than every declared dependency. Repeat to select more than one. Implies '--if_exists update', so naming both is optional; supplying it alongside a different '--if_exists' value is an error. Use it to bound what a transaction may change.

  -ch, --channels : list of str
    Comma-delimited package-manager channels, searched ahead of the channels declared by the selected environment YAML. Declared channels are retained unless '--override_channels' is also given. On environment creation the resulting list is rendered into a temporary copy of the YAML, which is what the package manager installs from; the tracked YAML is never modified. Applies to creation and to '--if_exists update' alike. A temporary package-manager configuration is also rendered and applied to that command alone, disabling any configured channel redirection, so that a supplied channel is fetched from the host it names rather than from a mirror declared elsewhere; the caller's own configuration is left untouched.

  -oc, --override_channels : flag
    Search only the channels given by '--channels', dropping those declared by the environment YAML and the package manager's 'defaults'. Requires '--channels'. On environment creation the declared channels are omitted from the rendered copy and a 'nodefaults' entry is added; on '--if_exists update' the declared channels are omitted from the command.

  -y, --yes : flag
    Automatically answer yes to package-manager prompts.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - mamba >= 1.5 or conda >= 24.7

  - User-facing environments are created from YAML files under 'install/envs/' via 'mamba env create -f <yaml>' (preferred) or 'conda env create -f <yaml>'.
  - See the following YAML files for environment package lists:

      | environment  | file path                     |
      | :---         | :---                          |
      | env_analyze  | install/envs/env_analyze.yml  |
      | env_protocol | install/envs/env_protocol.yml |
      | env_siqchip  | install/envs/env_siqchip.yml  |

  - Depending on the specified environment, a large number of packages and dependencies may need to be installed. As a result, the 'mamba env create' or 'conda env create' operation can take more than 10 minutes, especially for a fresh installation that does not make use of cached packages. In such cases, environment creation may take even longer (e.g., more than 20 or 30 minutes).
  - With '--if_exists reuse', this script skips environment creation. For 'env_protocol', it refreshes the repository's managed editable Python-package installation; it does not otherwise reconcile environment dependencies.
  - With '--if_exists update', this script reconciles a YAML-backed existing environment to its YAML, then refreshes the managed editable package when applicable. Declared versions win: a package whose YAML specification differs from the installed one is changed to match, which may mean a downgrade. Packages the YAML no longer declares are left in place rather than pruned. Channels come from the YAML, with any '--channels' values searched ahead of them, or replacing them entirely under '--override_channels'. Use repeatable '--update_package' selections to bound the transaction to specific packages.
  - Use '--channels' to install from site-specific channels, such as an institutional mirror, searched ahead of those listed in the YAML files, and add '--override_channels' to use the supplied channels alone. This is the supported path at sites that proxy package channels, and it is preferred over configuring channels globally: it is explicit, it is recorded in the command that produced the environment, and it leaves the rest of a user's package-manager configuration untouched.

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


# TODO: need to pin some packages in 'env_align'.
show_hidden=false
if [[ "${show_hidden}" == true ]]; then
    cat >&2 << EOM
  Hidden environments:
    - env_align          # Note: reproduce my 2023-4 processing and analysis.
      + bamtools
      + bbmap
      + bedtools
      + bowtie2
      + bwa
      + datamash
      + fastqc
      + gawk
      + gnuplot
      + macs3
      + minimap
      + mosdepth
      + parallel
      + picard
      + preseq
      + rename
      + samtools
      + subread
      + tree
      + ucsc-bedgraphtobigwig
      + ucsc-bedsort
      + ucsc-facount
      + wget
    - env_repro          # Note: reproduce work from previous lab members.
      + bc
      + bowtie2=2.3.4.2  # Note: explicitly pinning old version.
      + deeptools=3.3.1  # Note: explicitly pinning old version.
      + gawk
      + ipython
      + parallel
      + pbzip2
      + pigz
      + python=3.6       # Note: explicitly pinning old version.
      + rename
      + samtools=1.9     # Note: explicitly pinning old version.
      + tree
      + wget
EOM
fi
