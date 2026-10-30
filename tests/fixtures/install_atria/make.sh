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

        checkout) printf '%s\n' "${3:-}" > .git/current ;;

        rev-parse)
            if [[ "${2:-}" == "HEAD" && -f .git/current ]] \
                && [[ "$(< .git/current)" == "v4.1.5" ]]
            then
                printf '222\n'
            elif [[ "${2:-}" == "HEAD" ]]; then
                printf '111\n'
            else
                printf '222\n'
            fi
            ;;

        rev-list)
            if [[ "${!#}" == "v4.1.4" ]]; then
                printf '111\n'
            else
                printf '222\n'
            fi
            ;;
    esac
}


# Create the Julia executable that produces the fixture Atria executable.
function write_julia_binary() {
    local dir_destination="${1}"
    local julia_version="${2}"

    mkdir -p "${dir_destination}/julia-${julia_version}/bin"
    cat > "${dir_destination}/julia-${julia_version}/bin/julia" << JULIA
#!/usr/bin/env bash

if [[ "\${1:-}" == "--version" ]]; then
    printf 'julia version ${julia_version}\n'
    exit 0
fi

fake_version="\${ATRIA_FAKE_VERSION:-v4.1.5}"
dir_version="\${ATRIA_FAKE_DIR_VERSION:-\${fake_version}}"
dir_app="\$(pwd)/atria-\${dir_version#v}"

if [[ -d "\${dir_app}" ]]; then
    dir_app+="_julia-${julia_version}"
fi

mkdir -p "\${dir_app}/bin"

cat > "\${dir_app}/bin/julia" << APP_JULIA
#!/usr/bin/env bash

if [[ "\\\${1:-}" == "--version" ]]; then
    printf 'julia version ${julia_version}\n'
fi
APP_JULIA

cat > "\${dir_app}/bin/atria" << ATRIA
#!/usr/bin/env bash

if [[ "\\\${1:-}" == "--version" ]]; then
    printf '%s\n' "\${fake_version}"
else
    printf 'Atria help\n'
fi
ATRIA

chmod +x "\${dir_app}/bin/atria"
chmod +x "\${dir_app}/bin/julia"
JULIA

    chmod +x "${dir_destination}/julia-${julia_version}/bin/julia"
}


# Emulate extraction of the Julia archive into the final command argument.
function handle_tar() {
    local arg=""
    local dir_destination="${!#}"
    local tarball=""
    local julia_version=""

    if [[ "${ATRIA_FAKE_TAR_FAIL:-false}" == "true" ]]; then
        echo "error(fake_tool): forced Julia extraction failure." >&2
        return 1
    fi

    for arg in "$@"; do
        if [[ "${arg##*/}" == julia-*.tar.gz ]]; then
            tarball="${arg##*/}"
        fi
    done

    if [[ "${tarball}" =~ ^julia-([0-9]+\.[0-9]+\.[0-9]+)- ]]; then
        julia_version="${BASH_REMATCH[1]}"
    else
        echo "error(fake_tool): could not derive Julia version." >&2
        return 1
    fi

    write_julia_binary "${dir_destination}" "${julia_version}"
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

# Write a selected Julia that reports a version no release has. A failed Julia
# update must keep it.
cat << 'EOM' > "${fil_stl_jul}"
#!/usr/bin/env bash

printf 'julia version 0.0.0\n'
EOM
chmod +x "${fil_stl_jul}"

succeed "generated install_atria fixtures under ${dir_scr}"
