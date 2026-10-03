#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: help_install_preseq.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


# shellcheck disable=SC2154
function help_install_preseq() {
    cat >&2 << EOM
Usage
-----
  install_preseq.sh
    [--help] [--dry_run]
    [--env_nam <str>]
    [--dir_tmp <dir>] [--fil_tar <file>]
    [--threads <int>] [--if_exists <spec>]

  Build 'preseq' ${VERSION_PRESEQ} from its release tarball and install it into a Conda environment, by default 'env_qc'.

  No Bioconda build of 'preseq' currently coexists with the packages in 'env_qc' on Linux, so this script builds it against that environment's own compiler, 'make', and 'htslib'.

Parameters
----------
  -h, --help : flag
    Show this help message and exit.

  -dr, --dry, --dry_run : flag
    Run script in dry-run mode. Print the resolved plan without activating the environment, downloading, building, or installing.

  -en, --env, --env_nam : str
    Conda environment to activate. 'preseq' is built against this environment and installed into its prefix (default: '${env_nam}').

  -dt, --tmp, --dir_tmp : dir
    Working directory for the tarball and the build. If unset, a temporary directory is created under the system temp location and removed on exit.

  -ft, --fil_tar : file
    Local copy of the 'preseq' ${VERSION_PRESEQ} release tarball, used instead of downloading it. Its SHA-256 is verified in the same way.

  -th, --thr, --threads : int
    Number of threads to use. Passed to 'make' as its job count (default: ${threads}).

  -ie, --if_exists : {'fail', 'reuse', 'update'}
    What to do if 'preseq' is already installed in the environment (default: '${if_exists}'). 'fail' stops without changing anything; 'reuse' keeps an installation that reports version ${VERSION_PRESEQ} and fails otherwise; 'update' rebuilds and reinstalls it.

Notes
-----
  Runtime requirements:
    - awk
    - bash >= 4.4
    - conda or mamba (when the requested environment is not active)
    - curl (when '--fil_tar' is not specified)
    - grep
    - mktemp
    - mv
    - Network access (when neither '--dry_run' nor '--fil_tar' is specified)
    - sha256sum or shasum (when '--dry_run' is not specified)
    - tar
    - The environment's C++ compiler, 'make', and 'htslib', which 'install/envs/env_qc.yml' declares

  - The tarball's SHA-256 is fixed in this script and matches the Bioconda recipe for 'preseq' ${VERSION_PRESEQ}. A mismatch stops the script before anything is built.
  - The script adds '#include <cstdint>' to 'src/smithlab_cpp/sam_record.hpp'. GCC 13 and later no longer provide that header indirectly, and 'preseq' ${VERSION_PRESEQ} does not compile on Linux without it.
  - 'preseq' is configured with '--enable-hts', so it reads BAM input as well as duplicate-count histograms.
  - The installed program is checked for its version and by running 'lc_extrap' on a small synthetic histogram.
  - Create the environment first, for example with 'bash install/scripts/install_envs.sh --env_nam env_qc'.

Examples
--------
  1. Preview the build into 'env_qc'.
    '''bash
    bash install/scripts/install_preseq.sh \\
        --dry_run
    '''

  2. Build into 'env_qc' from a local copy of the release tarball, rebuilding any existing installation.
    '''bash
    bash install/scripts/install_preseq.sh \\
        --env_nam env_qc \\
        --fil_tar "\${HOME}/downloads/preseq-${VERSION_PRESEQ}.tar.gz" \\
        --threads 4 \\
        --if_exists update
    '''
EOM
}
