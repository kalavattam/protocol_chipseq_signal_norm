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

# Declare every generated path up front. Each 'prefix/' state stands in for a
# managed environment prefix: the compiler, 'make', and 'htslib' the installer
# looks for, none of them functional.
dir_pfx="${dir_fix}/prefix"
dir_tar="${dir_fix}/tarball"
arr_sta=( "bare" "no_compiler" "installed" )

fil_tar="${dir_tar}/not_a_tarball.tar.gz"
arr_gen=( "${fil_tar}" )

for sta in "${arr_sta[@]}"; do
    arr_gen+=(
        "${dir_pfx}/${sta}/bin/make"
        "${dir_pfx}/${sta}/bin/clang++"
        "${dir_pfx}/${sta}/bin/preseq"
        "${dir_pfx}/${sta}/include/htslib/sam.h"
        "${dir_pfx}/${sta}/lib/libhts.so.3"
    )
done


# Remove stale outputs so regeneration is idempotent.
rm_files "${dir_fix}" "${arr_gen[@]}"

mkdirs "${dir_tar}"


# Write a file named like a tarball that is not one, so its checksum never
# matches.
printf 'not a tarball\n' > "${fil_tar}"

# Write each prefix state. 'no_compiler' lacks 'clang++', and 'installed'
# already holds a 'preseq' that reports version 3.2.0. The stubs exit 1, since
# the installer must stop at its checks before it would use them.
for sta in "${arr_sta[@]}"; do
    pth="${dir_pfx}/${sta}"

    mkdirs "${pth}/bin" "${pth}/include/htslib" "${pth}/lib"
    printf '#!/usr/bin/env bash\nexit 1\n' > "${pth}/bin/make"
    : > "${pth}/include/htslib/sam.h"
    : > "${pth}/lib/libhts.so.3"

    if [[ "${sta}" != "no_compiler" ]]; then
        printf '#!/usr/bin/env bash\nexit 1\n' > "${pth}/bin/clang++"
    fi

    if [[ "${sta}" == "installed" ]]; then
        printf '#!/usr/bin/env bash\necho "Version: 3.2.0" >&2\n' \
            > "${pth}/bin/preseq"
    fi

    chmod +x "${pth}/bin/"*
done


succeed "generated install_preseq fixtures under ${dir_fix}"
