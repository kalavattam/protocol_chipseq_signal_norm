#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: help_install_atria.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in design, development,
# and documentation, with all output reviewed, edited, and approved by the
# author.
#
# Distributed under the MIT license.


# shellcheck disable=SC2154
function help_install_atria() {
    cat << EOM >&2
Usage
-----
  install_atria.sh
    [--help] [--dry_run]
    [--env_nam <str>]
    [--dir_install <dir>] [--dir_tmp <dir>]
    [--v_julia <spec>] [--v_atria <str>]
    [--path_snippet <file>] [--if_exists <spec>]

  Install Julia and Atria into a user-controlled directory without 'sudo'.

  The protocol-tested/default Julia version is 1.8.5. Upstream Atria docs describe Atria as tested with Julia 1.8 and 1.9, and recommend Julia 1.8.5 for speed relative to Julia 1.9.

Parameters
----------
  -h, --help : flag
    Show this help message and exit.

  -dr, --dry, --dry_run : flag
    Run script in dry-run mode. Print planned actions, derived paths, and URLs without downloading, extracting, cloning, building, or writing a PATH snippet.

  -en, --env, --env_nam : str
    Conda environment to activate before checking Atria runtime dependencies and building Atria (default: '${env_nam}').

  -di, --dir_install : dir
    User-controlled installation directory for Julia and Atria (default: '${dir_install}').

  -dt, --tmp, --dir_tmp : dir
    Working directory for downloaded files. If unset, a temporary directory is created under the system temp location and cleaned up on exit.

  -vj, --v_julia, --julia_version : str
    Julia version to install (default: '${v_julia}'). Verified SHA-256 mappings are bundled for Julia 1.8.0-1.8.5 and 1.9.0-1.9.4 on supported Linux and macOS x86_64 and aarch64 targets.

  -va, --v_atria, --atria_version : str
    Atria version to install (default: '${v_atria}'; the corresponding Git tag is expected to be 'v${v_atria}').

  -ps, --path_snippet : file
    Append PATH export lines to the requested file. This may be a shell configuration file such as '${HOME}/.bashrc' or '${HOME}/.zshrc', or a temporary snippet file for, e.g., testing. Without this option, the PATH lines are printed at the end.

  -ie, --if_exists : {'fail', 'reuse'}
    What to do if Julia and/or Atria already exist in the requested installation directory: 'fail' or 'reuse' (default: '${if_exists}').

Notes
-----
  Runtime requirements:
    - A compatible Conda environment providing the listed tools
    - bash >= 4.4
    - conda (when the requested environment is not active)
    - cp
    - curl
    - find
    - git
    - grep
    - head
    - mkdir
    - mktemp
    - Network access (when '--dry_run' is not specified)
    - pbzip2
    - pigz
    - rm
    - Rscript
    - sha256sum or shasum (when '--dry_run' is not specified)
    - sort
    - tail
    - tar

  - Supported Julia archive targets:
    + Linux glibc x86_64
    + Linux glibc aarch64 / arm64
    + macOS x86_64
    + macOS arm64 / aarch64
  - A non-dry run requires network access to download Julia and clone/fetch Atria.
  - BSD tar, which is commonly installed by default on macOS systems, is acceptable; extraction is verified by checking the expected Julia executable afterward.
  - Rscript, pigz, and pbzip2 are expected to come from the active project environment.

Examples
--------
  1. Preview the default Julia and Atria installation under a user-controlled directory.
    '''bash
    bash install/scripts/install_atria.sh \\
        --dry_run \\
        --dir_install "\${HOME}/opt/atria"
    '''

  2. Preview pinned versions, an explicit work directory, and reusable existing installs.
    '''bash
    bash install/scripts/install_atria.sh \\
        --dry_run \\
        --env_nam env_protocol \\
        --dir_install "\${HOME}/opt/atria" \\
        --dir_tmp "\${TMPDIR:-/tmp}/atria-build" \\
        --v_julia 1.8.5 \\
        --v_atria 3.2.2 \\
        --if_exists reuse
    '''
EOM
}
