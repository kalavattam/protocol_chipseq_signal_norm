#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: make.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-5.6);
# - Anthropic Claude Code (Opus 5.5).
#
# Distributed under the MIT license.


# Require Bash >= 4.4 before doing any work.
if [[ -z "${BASH_VERSION:-}" ]] || (( BASH_VERSINFO[0] < 4 )); then
    echo "error(make.sh):" \
        "Bash >= 4.4 is required." >&2
    exit 1
fi

# Run in safe mode, exiting on errors, unset variables, and pipe failures.
set -euo pipefail

# Resolve paths relative to 'tests/fixtures'.
dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source shared fixture-generation helpers.
# shellcheck source=tests/support/fixture_helpers.sh
source "${dir_scr}/../../support/fixture_helpers.sh"

# Declare generated package-manager tool paths.
dir_tool="${dir_scr}/tool"
fil_tool="${dir_tool}/fake_package_manager"
arr_tool_nam=( "conda" "mamba" )
arr_gen=( "${fil_tool}" )

for name in "${arr_tool_nam[@]}"; do
    arr_gen+=( "${dir_tool}/${name}" )
done

# Declare the starting states of the fake environment list, which the fake
# package manager appends to as it creates environments.
dir_sta="${dir_scr}/state"
arr_sta_nam=( "none" "env_protocol" "env_siqchip" )

for name in "${arr_sta_nam[@]}"; do
    arr_gen+=( "${dir_sta}/${name}.txt" )
done

# Remove stale outputs so regeneration is idempotent.
rm_files "${dir_scr}" "${arr_gen[@]}"
mkdirs "${dir_tool}" "${dir_sta}"

# Install the stub package manager under the name the tests call.
cp "${dir_scr}/stub_package_manager.sh" "${fil_tool}"
chmod +x "${fil_tool}"
for name in "${arr_tool_nam[@]}"; do
    ln -sf fake_package_manager "${dir_tool}/${name}"
done

# Write the starting states: no environment, or one named environment.
: > "${dir_sta}/none.txt"
printf '%s\n' env_protocol > "${dir_sta}/env_protocol.txt"
printf '%s\n' env_siqchip > "${dir_sta}/env_siqchip.txt"

succeed "generated install_envs fixtures under ${dir_scr}"
