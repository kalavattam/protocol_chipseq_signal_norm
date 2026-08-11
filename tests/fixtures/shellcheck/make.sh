#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: make.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5) was used in design, development, and
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

# Declare every generated path up front
dir_script="${dir_fix}/script"
dir_tool="${dir_fix}/tool"
fil_bash="${dir_script}/bash.sh"
fil_posix="${dir_script}/posix.sh"
fil_stub="${dir_tool}/fake_shellcheck.sh"


# Remove stale outputs so regeneration is idempotent.
rm_files "${dir_fix}" "${fil_bash}" "${fil_posix}" "${fil_stub}"

mkdirs "${dir_script}" "${dir_tool}"


# Two scripts whose shebang is the only thing under test. The runner copies
# them into a synthetic repository and must split them by declared language, so
# each body is the shortest program that a shebang can head.
cat << 'EOM' > "${fil_bash}"
#!/usr/bin/env bash
printf 'bash fixture\n'
EOM

cat << 'EOM' > "${fil_posix}"
#!/bin/sh
printf 'POSIX fixture\n'
EOM

# One ShellCheck stub, standing in for the real tool so the runner's executable
# selection, structured-finding parsing, and exit-status handling are exercised
# without applying repository findings. The stub is written in two heredocs
# around one assembled line: its findings payload is a single line of JSON
# output, and building that line from fragments keeps the recipe itself inside
# the shell line-length budget.
cat << 'EOM' > "${fil_stub}"
#!/usr/bin/env bash
set -u

if [[ "${1:-}" == "--version" ]]; then
    cat << 'VERSION'
ShellCheck - shell script analysis tool
version: 0.10.0
license: GNU General Public License, version 3
website: https://www.shellcheck.net
VERSION
    exit 0
fi

if [[ -n "${FAKE_SHELLCHECK_LOG:-}" ]]; then
    printf '%s\n' "$*" >> "${FAKE_SHELLCHECK_LOG}"
fi

comment='{"file":"fixture.sh","line":1,"column":1,'
comment+='"level":"warning","code":9999,'
comment+='"message":"fixture finding","fix":null}'

status="${FAKE_SHELLCHECK_STATUS:-0}"
case "${status}" in
    0)
        printf '{"comments":[]}\n'
        ;;
    1)
        printf '{"comments":[%s]}\n' "${comment}"
        ;;
    *)
        echo "fixture infrastructure failure" >&2
        ;;
esac
exit "${status}"
EOM

chmod +x "${fil_stub}"


succeed "generated ShellCheck runner fixtures under ${dir_fix}"
