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
    echo "error(make.sh): Bash >= 4.4 is required." >&2
    exit 1
fi

# Run in safe mode, exiting on errors, unset variables, and pipe failures.
set -euo pipefail

# Resolve paths relative to 'tests/fixtures'.
dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source shared fixture-generation helpers.
# shellcheck source=tests/support/fixture_helpers.sh
source "${dir_scr}/../../support/fixture_helpers.sh"

# Declare generated tool paths.
dir_tool="${dir_scr}/tool"
fil_tool="${dir_tool}/fake_tool"
arr_tool_nam=(git curl sha256sum tar julia pigz pbzip2 Rscript)
arr_gen=( "${fil_tool}" )

for name in "${arr_tool_nam[@]}"; do
    arr_gen+=( "${dir_tool}/${name}" )
done

# Declare the PATH snippets, sentinel contents, and stale state that the Atria
# test plants for the installer.
dir_snp="${dir_scr}/snippet"
fil_snp_lgc="${dir_snp}/legacy_path.sh"
fil_snp_bad="${dir_snp}/malformed_path.sh"
fil_snp_snt="${dir_snp}/sentinel_path.sh"
fil_snt_atr="${dir_snp}/sentinel_atria.txt"
fil_stl_mft="${dir_snp}/stale_manifest.toml"
fil_stl_jul="${dir_snp}/julia_0.0.0"
arr_gen+=(
    "${fil_snp_lgc}"
    "${fil_snp_bad}"
    "${fil_snp_snt}"
    "${fil_snt_atr}"
    "${fil_stl_mft}"
    "${fil_stl_jul}"
)

# Remove stale outputs so regeneration is idempotent.
rm_files "${dir_scr}" "${arr_gen[@]}"
mkdirs "${dir_tool}" "${dir_snp}"

# Install the stub tool under the name the tests call.
cp "${dir_scr}/stub_tool.sh" "${fil_tool}"
chmod +x "${fil_tool}"
for name in "${arr_tool_nam[@]}"; do
    ln -sf fake_tool "${dir_tool}/${name}"
done

# Write a PATH snippet holding the exports of an earlier installation that
# carried no markers.
{
    printf '%s\n' 'export PATH="${PATH}:/legacy/julia-1.8.5/bin"'
    printf '%s\n' 'export PATH="/legacy/Atria/current/bin:${PATH}"'
} > "${fil_snp_lgc}"

# Write a PATH snippet whose installer markers are in the wrong order.
{
    printf '%s\n' '# <<< protocol_chipseq_signal_norm install_atria.sh <<<'
    printf '%s\n' '# >>> protocol_chipseq_signal_norm install_atria.sh >>>'
} > "${fil_snp_bad}"

# Write sentinel contents for a snippet and an Atria checkout that a failed
# installation must leave unchanged.
printf '%s\n' "PATH sentinel" > "${fil_snp_snt}"
printf '%s\n' "Atria sentinel" > "${fil_snt_atr}"

# Write the generated dependency state that an Atria tag update must clear.
printf '%s\n' "# stale generated dependency state" > "${fil_stl_mft}"

# Install a selected Julia that reports a version no release has. A failed
# Julia update must keep it.
cp "${dir_scr}/stub_julia.sh" "${fil_stl_jul}"
chmod +x "${fil_stl_jul}"

succeed "generated install_atria fixtures under ${dir_scr}"
