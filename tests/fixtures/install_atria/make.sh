#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: make.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
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

# Remove stale outputs so regeneration is idempotent.
rm_files "${dir_scr}" "${arr_gen[@]}"
mkdirs "${dir_tool}"

cat << 'EOM' > "${fil_tool}"
#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: fake_tool
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.6) were used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


# Require Bash >= 4.4 before doing any work.
if [[ -z "${BASH_VERSION:-}" ]] || (( BASH_VERSINFO[0] < 4 )); then
    echo "error(fake_tool): Bash >= 4.4 is required." >&2
    exit 1
fi

# Run in safe mode, exiting on errors, unset variables, and pipe failures.
set -euo pipefail

# Identify the fake command selected through its symlink.
cmd_nam="${0##*/}"

# Record every invocation for the contract test.
function log_invocation() {
    {
        printf '%s' "${cmd_nam}"
        printf ' %q' "$@"
        printf '\n'
    } >> "${ATRIA_FAKE_LOG}"
}


# Return the stable Atria tags available to the resolver.
function print_tags() {
    cat << 'TAGS'
111 refs/tags/v4.1.4
222 refs/tags/v4.1.5
333 refs/tags/v4.1.5-dev
TAGS
}


# Emulate the Git operations exercised by the installer contract.
function handle_git() {
    local dir_dest=""

    case "${1:-}" in
        ls-remote) print_tags ;;

        clone)
            dir_dest="${!#}"
            mkdir -p "${dir_dest}/.git"
            ;;

        fetch) : ;;

        checkout) touch .git/current ;;

        rev-parse)
            if [[ "${2:-}" == "HEAD" && -f .git/current ]]; then
                printf '222\n'
            elif [[ "${2:-}" == "HEAD" ]]; then
                printf '111\n'
            else
                printf '222\n'
            fi
            ;;

        rev-list) printf '222\n' ;;
    esac
}


# Create the Julia executable that produces the fixture Atria executable.
function write_julia_binary() {
    local dir_destination="${1}"

    mkdir -p "${dir_destination}/julia-1.8.5/bin"
    cat > "${dir_destination}/julia-1.8.5/bin/julia" << 'JULIA'
#!/usr/bin/env bash

if [[ "${1:-}" == "--version" ]]; then
    printf 'julia version 1.8.5\n'
    exit 0
fi

mkdir -p "$(pwd)/atria-4.1.5/bin"

cat > "$(pwd)/atria-4.1.5/bin/atria" << 'ATRIA'
#!/usr/bin/env bash

if [[ "${1:-}" == "--version" ]]; then
    printf 'v4.1.5\n'
else
    printf 'Atria help\n'
fi
ATRIA

chmod +x "$(pwd)/atria-4.1.5/bin/atria"
JULIA

    chmod +x "${dir_destination}/julia-1.8.5/bin/julia"
}


# Emulate extraction of the Julia archive into the final command argument.
function handle_tar() {
    local dir_destination="${!#}"

    write_julia_binary "${dir_destination}"
}


# Emulate curl's output-file behavior used by the installer.
function handle_curl() {
    local fil_out=""

    while (( $# > 0 )); do
        case "${1}" in
            -o)
                fil_out="${2:-}"
                shift 2
                ;;

            *) shift ;;
        esac
    done

    : > "${fil_out}"
}


log_invocation "$@"

case "${cmd_nam}" in
    git) handle_git "$@" ;;

    curl) handle_curl "$@" ;;

    sha256sum)
        cat > /dev/null
        printf 'fixture: OK\n'
        ;;

    tar) handle_tar "$@" ;;

    julia) printf 'julia version 0.0.0\n' ;;

    pigz|pbzip2|Rscript) : ;;
esac
EOM

chmod +x "${fil_tool}"
for name in "${arr_tool_nam[@]}"; do
    ln -sf fake_tool "${dir_tool}/${name}"
done

succeed "generated install_atria fixtures under ${dir_scr}"
