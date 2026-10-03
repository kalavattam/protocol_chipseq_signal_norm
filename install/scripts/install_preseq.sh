#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: install_preseq.sh
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
    echo "error(shell):" \
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


# Resolve paths to installation support and repository directories.
dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
dir_ins="$(cd "${dir_scr}/.." > /dev/null 2>&1 && pwd)"
dir_rep="$(cd "${dir_ins}/.." > /dev/null 2>&1 && pwd)"

# The release built here, with the checksum the Bioconda recipe also uses.
readonly VERSION_PRESEQ="3.2.0"
readonly URL_PRESEQ="https://github.com/smithlabcode/preseq/releases/download/v${VERSION_PRESEQ}/preseq-${VERSION_PRESEQ}.tar.gz"
readonly SHA256_PRESEQ="95b81c9054e0d651de398585c7e96b807ad98f0bdc541b3e46665febbe2134d9"

# The header that relies on '<cstdint>' arriving indirectly, which GCC 13 and
# later no longer provide.
readonly PTH_PATCH="src/smithlab_cpp/sam_record.hpp"


# Source shared helpers.
function source_helpers_script() {
    local fnc_src

    dir_fnc="${dir_rep}/lib/bash"
    fnc_src="${dir_fnc}/core/source_helpers.sh"

    if [[ ! -f "${fnc_src}" ]]; then
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "script not found: '${fnc_src}'." >&2
        return 1
    fi

    # shellcheck disable=SC1090
    source "${fnc_src}" || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source '${fnc_src}'." >&2
        return 1
    }

    source_helpers "${dir_fnc}" \
        check_args \
        check_env \
        check_inputs \
        check_numbers \
        format_outputs \
        handle_env \
        help/help_install_preseq \
        || {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to source required helper scripts." >&2
            return 1
        }
}


function init_arg_defs() {
    dry_run=false
    env_nam="env_qc"
    dir_tmp=""
    tmp_auto=false
    fil_tar=""
    threads=2
    if_exists="fail"
    sha_cmd=""
    pth_prefix=""
    pth_cxx=""
    pth_make=""
}


# Define utility functions.
function echo_dry() {
    printf 'dryrun(%s): %s\n' "$(basename "${BASH_SOURCE[0]}")" "$*"
}


function check_pkg_mgr() {
    if \
        command -v mamba > /dev/null 2>&1
    then
        return 0
    elif \
        command -v conda > /dev/null 2>&1
    then
        return 0
    fi

    echo_err "neither 'mamba' nor 'conda' is available in PATH."
    return 1
}


function select_sha_cmd() {
    if \
        command -v sha256sum > /dev/null 2>&1
    then
        sha_cmd="sha256sum"
    elif \
        command -v shasum > /dev/null 2>&1
    then
        sha_cmd="shasum"
    else
        echo_err "neither 'sha256sum' nor 'shasum' is available in PATH."
        return 1
    fi
}


function verify_sha256() {
    local file="${1:-}"
    local expect="${2:-}"

    validate_var "file" "${file}" || return 1
    validate_var "expect" "${expect}" || return 1

    case "${sha_cmd}" in
        sha256sum)
            (
                cd "$(dirname "${file}")"
                printf '%s  %s\n' "${expect}" "$(basename "${file}")" \
                    | sha256sum -c -
            )
            ;;

        shasum)
            (
                cd "$(dirname "${file}")"
                printf '%s  %s\n' "${expect}" "$(basename "${file}")" \
                    | shasum -a 256 -c -
            )
            ;;

        *)
            echo_err "unsupported SHA-256 command '${sha_cmd}'."
            return 1
            ;;
    esac
}


function prepare_tmp_dir() {
    if [[ -n "${dir_tmp}" ]]; then
        mkdir -p "${dir_tmp}"
        tmp_auto=false
    else
        dir_tmp="$(mktemp -d "${TMPDIR:-/tmp}/install_preseq.XXXXXX")"
        tmp_auto=true
    fi
}


function cleanup_tmp_dir() {
    if [[
        "${tmp_auto}" == "true" && -n "${dir_tmp}" && -d "${dir_tmp}"
    ]]; then
        rm -rf "${dir_tmp}"
    fi
}


# Prefer the compiler the environment's activation exports, then the
# environment's own compiler wrappers; never a system compiler, whose headers
# and runtime don't need to match the environment's 'htslib'.
function find_env_cxx() {
    local cand=""

    if [[ -n "${CXX:-}" ]]; then
        cand="$(command -v "${CXX}" 2> /dev/null || true)"

        if [[ -n "${cand}" && "${cand}" == "${pth_prefix}/"* ]]; then
            pth_cxx="${cand}"
            return 0
        fi
    fi

    for cand in \
        "${pth_prefix}"/bin/*-conda-linux-gnu-c++ \
        "${pth_prefix}"/bin/*-apple-darwin*-clang++ \
        "${pth_prefix}/bin/clang++" \
        "${pth_prefix}/bin/g++"
    do
        if [[ -x "${cand}" ]]; then
            pth_cxx="${cand}"
            return 0
        fi
    done

    echo_err \
        "no C++ compiler was found in '${pth_prefix}/bin'. Create or update" \
        "'${env_nam}' from its YAML, which declares 'cxx-compiler'."
    return 1
}


function find_env_make() {
    if [[ -x "${pth_prefix}/bin/make" ]]; then
        pth_make="${pth_prefix}/bin/make"
        return 0
    fi

    echo_err \
        "'make' was not found in '${pth_prefix}/bin'. Create or update" \
        "'${env_nam}' from its YAML, which declares 'make'."
    return 1
}


function check_env_htslib() {
    local lib=""

    if [[ ! -f "${pth_prefix}/include/htslib/sam.h" ]]; then
        echo_err \
            "'htslib' headers were not found under '${pth_prefix}/include'." \
            "Create or update '${env_nam}' from its YAML, which declares" \
            "'htslib'."
        return 1
    fi

    for lib in "${pth_prefix}"/lib/libhts.so* "${pth_prefix}"/lib/libhts.*dylib
    do
        if [[ -e "${lib}" ]]; then
            return 0
        fi
    done

    echo_err \
        "the 'htslib' library was not found under '${pth_prefix}/lib'." \
        "Create or update '${env_nam}' from its YAML, which declares 'htslib'."
    return 1
}


function parse_args() {
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
            -h|--hlp|--help)
                help_install_preseq
                return 2
                ;;

            -dr|--dry|--dry[_-]run)
                dry_run=true
                shift 1
                ;;

            -en|--env|--env[_-]nam)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_install_preseq
                    return 1
                }
                env_nam="${2}"
                shift 2
                ;;

            -dt|--tmp|--dir[_-]tmp)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_install_preseq
                    return 1
                }
                dir_tmp="${2}"
                shift 2
                ;;

            -ft|--fil[_-]tar)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_install_preseq
                    return 1
                }
                fil_tar="${2}"
                shift 2
                ;;

            -th|--thr|--threads)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_install_preseq
                    return 1
                }
                threads="${2}"
                shift 2
                ;;

            -ie|--if[_-]exists)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_install_preseq
                    return 1
                }
                if_exists="${2}"
                shift 2
                ;;

            *)
                echo_err "unknown option/parameter passed: '${1}'."
                echo >&2
                help_install_preseq
                return 1
                ;;
        esac
    done
}


function validate_args() {
    validate_var "env_nam" "${env_nam}"     || return 1
    validate_var "if_exists" "${if_exists}" || return 1
    check_int_pos "${threads}" "threads"    || return 1

    case "${if_exists}" in
        fail|reuse|update) : ;;
        *)
            echo_err \
                "invalid '--if_exists' value: '${if_exists}'. Must be" \
                "'fail', 'reuse', or 'update'."
            return 1
            ;;
    esac

    if [[ -n "${fil_tar}" && ! -f "${fil_tar}" ]]; then
        echo_err "tarball does not exist: '${fil_tar}'."
        return 1
    fi
}


function report_plan() {
    echo
    echo "Resolved installation plan:"
    echo "  - env_nam=${env_nam}"
    echo "  - prefix=${pth_prefix:-<prefix of ${env_nam}>}"
    echo "  - version=${VERSION_PRESEQ}"
    echo "  - source=${fil_tar:-${URL_PRESEQ}}"
    echo "  - sha256=${SHA256_PRESEQ}"
    echo "  - patch=add '#include <cstdint>' to ${PTH_PATCH}"
    echo "  - configure=--prefix=<prefix> --enable-hts"
    echo "  - dir_tmp=${dir_tmp:-AUTO}"
    echo "  - threads=${threads}"
    echo "  - if_exists=${if_exists}"
    echo
}


function check_runtime_req() {
    local cmd=""

    if [[ "${dry_run}" == "true" ]]; then
        echo_dry "would confirm 'mamba' or 'conda' is available."
        echo_dry \
            "would confirm environment '${env_nam}' exists and activate it."
        echo_dry \
            "would confirm '${env_nam}' provides a C++ compiler, 'make', and" \
            "'htslib' headers and library."
        echo_dry \
            "would select 'sha256sum' or 'shasum' for SHA-256 verification."
        return 0
    fi

    if [[ "${CONDA_DEFAULT_ENV:-}" == "${env_nam}" ]]; then
        echo "Environment '${env_nam}' is already active; reusing it."
    else
        check_pkg_mgr || return 1

        if ! \
            check_env_installed "${env_nam}"
        then
            echo_err \
                "environment '${env_nam}' appears not to be installed." \
                "Create it with 'install/scripts/install_envs.sh' first."
            return 1
        fi

        handle_env "${env_nam}" || return 1
    fi

    validate_var_dir "CONDA_PREFIX" "${CONDA_PREFIX:-}" || return 1
    pth_prefix="${CONDA_PREFIX}"

    find_env_cxx     || return 1
    find_env_make    || return 1
    check_env_htslib || return 1

    for cmd in awk grep mktemp mv tar; do
        if ! \
            command -v "${cmd}" > /dev/null 2>&1
        then
            echo_err "required command '${cmd}' is not available in PATH."
            return 1
        fi
    done

    if [[ -z "${fil_tar}" ]]; then
        check_pgrm_path curl || return 1
    fi

    select_sha_cmd || return 1
}


# Return 0 to proceed with a build, or 3 when a verified installation is reused
# and nothing remains to do.
function check_existing_install() {
    local pth_bin="${pth_prefix}/bin/preseq"

    if [[ ! -e "${pth_bin}" ]]; then
        return 0
    fi

    case "${if_exists}" in
        fail)
            echo_err "'preseq' is already installed: '${pth_bin}'."
            echo_err \
                "nothing was changed. To keep it, rerun with" \
                "'--if_exists reuse'; to rebuild it, use '--if_exists update'."
            return 1
            ;;

        reuse)
            verify_preseq_version "${pth_bin}" || {
                echo_err \
                    "the installed 'preseq' is not version" \
                    "${VERSION_PRESEQ}. Rerun with '--if_exists update' to" \
                    "rebuild it."
                return 1
            }
            echo "Reusing the installed 'preseq' ${VERSION_PRESEQ}:" \
                "'${pth_bin}'."
            return 3
            ;;

        update)
            echo "Rebuilding 'preseq' ${VERSION_PRESEQ} over '${pth_bin}'."
            return 0
            ;;
    esac
}


function fetch_tarball() {
    local pth_tar="${dir_tmp}/preseq-${VERSION_PRESEQ}.tar.gz"

    if [[ -n "${fil_tar}" ]]; then
        cp "${fil_tar}" "${pth_tar}"
    else
        curl -fsSL -o "${pth_tar}" "${URL_PRESEQ}" || {
            echo_err "failed to download '${URL_PRESEQ}'."
            return 1
        }
    fi

    verify_sha256 "${pth_tar}" "${SHA256_PRESEQ}" || {
        echo_err \
            "SHA-256 mismatch for '${pth_tar}'; refusing to build from it."
        return 1
    }

    tar -xzf "${pth_tar}" -C "${dir_tmp}"
}


# Insert '#include <cstdint>' before the header's first include, unless the
# header already has it.
function patch_cstdint() {
    local pth_hdr="${dir_tmp}/preseq-${VERSION_PRESEQ}/${PTH_PATCH}"
    local pth_new="${pth_hdr}.patched"

    if [[ ! -f "${pth_hdr}" ]]; then
        echo_err "header to patch not found: '${pth_hdr}'."
        return 1
    fi

    if \
        grep -q '^#include <cstdint>' "${pth_hdr}"
    then
        return 0
    fi

    if ! \
        grep -q '^#include' "${pth_hdr}"
    then
        echo_err "no '#include' line found in '${pth_hdr}'."
        return 1
    fi

    awk '
        !done && /^#include/ { print "#include <cstdint>"; done = 1 }
        { print }
    ' "${pth_hdr}" > "${pth_new}"
    mv "${pth_new}" "${pth_hdr}"
}


function build_preseq() {
    local dir_src="${dir_tmp}/preseq-${VERSION_PRESEQ}"

    (
        cd "${dir_src}"
        CXX="${pth_cxx}" \
        CXXFLAGS="-O3 -isystem ${pth_prefix}/include" \
        LDFLAGS="-L${pth_prefix}/lib -Wl,-rpath,${pth_prefix}/lib" \
            ./configure --prefix="${pth_prefix}" --enable-hts
        "${pth_make}" -j "${threads}"
        "${pth_make}" install
    ) || {
        echo_err "failed to build 'preseq' in '${dir_src}'."
        return 1
    }
}


function verify_preseq_version() {
    local pth_bin="${1:?}"
    local out=""

    out="$("${pth_bin}" 2>&1 || true)"
    [[ "${out}" == *"Version: ${VERSION_PRESEQ}"* ]]
}


# Run 'lc_extrap' on a small synthetic histogram: a successful extrapolation
# exercises the arithmetic core that the gates rely on.
function verify_preseq_run() {
    local pth_bin="${pth_prefix}/bin/preseq"
    local fil_hist="${dir_tmp}/smoke_hist.txt"
    local fil_out="${dir_tmp}/smoke_lc.txt"

    verify_preseq_version "${pth_bin}" || {
        echo_err \
            "'${pth_bin}' does not report version ${VERSION_PRESEQ}."
        return 1
    }

    awk 'BEGIN {
        for (j = 1; j <= 20; j++) printf "%d %d\n", j, int(1e6 * exp(-0.6 * j))
    }' > "${fil_hist}"

    "${pth_bin}" lc_extrap -H "${fil_hist}" -o "${fil_out}" -s 1e6 -e 1e7 \
        > /dev/null 2>&1 || {
        echo_err "'preseq lc_extrap' failed on the smoke-test histogram."
        return 1
    }

    if ! \
        head -n 1 "${fil_out}" | grep -q '^TOTAL_READS'
    then
        echo_err "'preseq lc_extrap' wrote unexpected output: '${fil_out}'."
        return 1
    fi
}


function report_dry_run_install() {
    echo_dry "would use working directory '${dir_tmp:-AUTO}'."

    if [[ -n "${fil_tar}" ]]; then
        echo_dry "would copy '${fil_tar}' and verify its SHA-256."
    else
        echo_dry "would download and verify '${URL_PRESEQ}'."
    fi

    echo_dry "would add '#include <cstdint>' to '${PTH_PATCH}'."
    echo_dry \
        "would configure with '--enable-hts' against the 'htslib' in" \
        "'${env_nam}', then build with ${threads} job(s) and install into" \
        "its prefix."
    echo_dry \
        "would confirm the installed 'preseq' reports version" \
        "${VERSION_PRESEQ} and runs 'lc_extrap' on a synthetic histogram."
}


function run_install() {
    local rc=0

    check_existing_install || rc=$?
    if (( rc == 3 )); then
        return 0
    elif (( rc != 0 )); then
        return "${rc}"
    fi

    prepare_tmp_dir
    trap cleanup_tmp_dir EXIT

    fetch_tarball
    patch_cstdint
    build_preseq
    verify_preseq_run

    echo
    echo "success($(basename "${BASH_SOURCE[0]}")):" \
        "installed and verified 'preseq' ${VERSION_PRESEQ} in" \
        "'${pth_prefix}'."
}


function main() {
    local rc=0

    source_helpers_script
    init_arg_defs

    parse_args "$@" || rc=$?
    if (( rc == 2 )); then
        return 0
    elif (( rc != 0 )); then
        return "${rc}"
    fi

    validate_args
    check_runtime_req
    report_plan

    if [[ "${dry_run}" == "true" ]]; then
        report_dry_run_install
        return 0
    fi

    run_install
}


main "$@"
