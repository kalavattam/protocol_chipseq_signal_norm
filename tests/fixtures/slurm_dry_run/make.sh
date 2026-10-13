#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: make.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


# Require Bash >= 4.4 before doing any work.
if [[ -z "${BASH_VERSION:-}" ]]; then
    echo "error(make.sh):" \
        "this script must be run under Bash >= 4.4." >&2
    exit 1
elif ((
    BASH_VERSINFO[0] < 4 || ( BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] < 4 )
)); then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "this script requires Bash >= 4.4; current version is" \
        "'${BASH_VERSION}'." >&2
    exit 1
fi

# Run in safe mode, exiting on errors, unset variables, and pipe failures.
set -euo pipefail


# Resolve paths relative to 'tests/fixtures'.
dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
dir_fix="${dir_scr}"

# Source shared fixture-generation helpers.
# shellcheck source=tests/support/fixture_helpers.sh
source "${dir_scr}/../../support/fixture_helpers.sh"

# Declare every generated path up front. 'tool/' holds the 'sbatch' stub that
# lets a driver's Slurm dry run pass its tool check on a machine without Slurm.
dir_tol="${dir_fix}/tool"

fil_sbt="${dir_tol}/sbatch"


# Remove stale outputs so regeneration is idempotent.
rm_files "${dir_fix}" "${fil_sbt}"

mkdirs "${dir_tol}"


# Write the 'sbatch' stub. A dry run prints the 'sbatch' call and must never
# submit it, so the stub fails if it is ever called.
cat << 'EOM' > "${fil_sbt}"
#!/usr/bin/env bash

echo "error(sbatch): test stub called; a dry run must not submit." >&2
exit 1
EOM

chmod +x "${fil_sbt}"


succeed "generated Slurm dry-run fixtures under ${dir_fix}"
