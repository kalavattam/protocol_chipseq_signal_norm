#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: install_atria.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in development and
# documentation.
#
# Distributed under the MIT license.


#  Require Bash >= 4.4 before doing any work
if [[ -z "${BASH_VERSION:-}" ]]; then
    echo "error(shell): this script must be run under Bash >= 4.4." >&2
    exit 1
elif ((
    BASH_VERSINFO[0] < 4 || ( BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] < 4 )
)); then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "this script requires Bash >= 4.4; current version is" \
        "'${BASH_VERSION}'." >&2
    exit 1
fi

#  Run in safe mode, exiting on errors, unset variables, and pipe failures
set -euo pipefail


#  Resolve paths to installation support and repository directories
dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
dir_ins="$(cd "${dir_scr}/.." > /dev/null 2>&1 && pwd)"
dir_rep="$(cd "${dir_ins}/.." > /dev/null 2>&1 && pwd)"


#  Source shared helpers
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
        format_outputs \
        handle_env \
        help/help_install_atria \
        || {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to source required helper scripts." >&2
            return 1
        }
}


function init_arg_defs() {
    dry_run=false
    env_nam="env_protocol"
    dir_install="${HOME}"
    dir_tmp=""
    tmp_auto=false
    v_julia="1.8.5"
    v_atria="4.1.5"
    path_snippet=""
    if_exists="fail"
    sha_cmd=""
}


#  Define utility functions
function echo_dry() {
    echo "dryrun($(basename "${BASH_SOURCE[0]}")):" "$@"
}


function run_or_print() {
    if [[ "${dry_run}" == "true" ]]; then
        printf 'dryrun(%s):' "$(basename "${BASH_SOURCE[0]}")"
        printf ' %q' "$@"
        printf '\n'
    else
        "$@"
    fi
}


function check_pkg_mgr() {
    if command -v mamba > /dev/null 2>&1; then
        return 0
    elif command -v conda > /dev/null 2>&1; then
        return 0
    fi

    echo_err "neither 'mamba' nor 'conda' is available in PATH."
    return 1
}


function select_sha_cmd() {
    if command -v sha256sum > /dev/null 2>&1; then
        sha_cmd="sha256sum"
    elif command -v shasum > /dev/null 2>&1; then
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
                printf \
                    '%s  %s\n' \
                    "${expect}" \
                    "$(basename "${file}")" \
                    | sha256sum -c -
            )
            ;;

        shasum)
            (
                cd "$(dirname "${file}")"
                printf \
                    '%s  %s\n' \
                    "${expect}" \
                    "$(basename "${file}")" \
                    | shasum -a 256 -c -
            )
            ;;

        *)
            echo_err "unsupported SHA-256 command '${sha_cmd}'."
            return 1
            ;;
    esac
}


function extract_julia_tar() {
    local tarball="${1:-}"
    local dir_dst="${2:-}"
    local log_tar=""

    validate_var_file "tarball" "${tarball}" || return 1
    validate_var_dir  "dir_dst" "${dir_dst}" || return 1

    log_tar="$(mktemp "${TMPDIR:-/tmp}/install_atria_tar.XXXXXX")"

    if \
        tar -xzf "${tarball}" -C "${dir_dst}" > "${log_tar}" 2>&1
    then
        rm -f "${log_tar}"
        return 0
    fi

    if [[ -x "${jul_bin}" ]]; then
        echo_warn \
            "tar returned a non-zero exit status while extracting Julia, but" \
            "the expected Julia executable was found: '${jul_bin}'."
        echo_warn \
            "Continuing because older macOS bsdtar versions may warn about" \
            "malformed pax extended attributes after otherwise successful" \
            "extraction."
        cat "${log_tar}" >&2
        rm -f "${log_tar}"
        return 0
    fi

    cat "${log_tar}" >&2
    rm -f "${log_tar}"
    echo_err "failed to extract Julia archive '${tarball}'."
    return 1
}


function verify_julia_exec() {
    local pth_jul="${1:-}"
    local jul_ver=""

    validate_var_file "pth_jul" "${pth_jul}" || return 1

    jul_ver="$("${pth_jul}" --version)"
    echo "${jul_ver}"

    if [[ "${jul_ver}" != "julia version ${v_julia}" ]]; then
        echo_err \
            "Julia executable exists at '${pth_jul}', but it does not match" \
            "the requested version '${v_julia}'. Reported version:" \
            "'${jul_ver}'."
        return 1
    fi
}


#  Intentionally set global Julia target variables used by later installer
#+ steps
function map_julia_target() {
    sys_os="$(uname -s)"
    sys_ar="$(uname -m)"

    jul_os=""
    jul_ar_dir=""
    jul_ar_fil=""
    jul_256=""

    case "${sys_os}:${sys_ar}" in
        Linux:x86_64|Linux:amd64)
            jul_os="linux"
            jul_ar_dir="x64"
            jul_ar_fil="x86_64"
            ;;

        Linux:aarch64|Linux:arm64)
            jul_os="linux"
            jul_ar_dir="aarch64"
            jul_ar_fil="aarch64"
            ;;

        Darwin:x86_64)
            jul_os="mac"
            jul_ar_dir="x64"
            jul_ar_fil="mac64"
            ;;

        Darwin:arm64|Darwin:aarch64)
            jul_os="mac"
            jul_ar_dir="aarch64"
            jul_ar_fil="macaarch64"
            ;;

        *)
            echo_err \
                "unsupported OS/architecture for Julia archive mapping:" \
                "'${sys_os}:${sys_ar}'. Supported targets are Linux x86_64," \
                "Linux aarch64/arm64, macOS x86_64, and macOS arm64/aarch64."
            return 1
            ;;
    esac

    v_min="${v_julia%.*}"

    case "${jul_os}" in
        linux) jul_tar="julia-${v_julia}-${jul_os}-${jul_ar_fil}.tar.gz" ;;
        mac)   jul_tar="julia-${v_julia}-${jul_ar_fil}.tar.gz"           ;;
        *)
            echo_err "unsupported Julia OS token '${jul_os}'."
            return 1
            ;;
    esac

    jul_url="https://julialang-s3.julialang.org/bin/${jul_os}/${jul_ar_dir}/${v_min}/${jul_tar}"

    case "${v_julia}:${jul_os}:${jul_ar_fil}" in
        1.9.4:linux:x86_64)   jul_256="07d20c4c2518833e2265ca0acee15b355463361aa4efdab858dad826cf94325c" ;;
        1.9.4:linux:aarch64)  jul_256="541d0c5a9378f8d2fc384bb8595fc6ffe20d61054629a6e314fb2f8dfe2f2ade" ;;
        1.9.4:mac:mac64)      jul_256="67eec264f6afc9e9bf72c0f62c84d91c2ebdfaed6a0aa11606e3c983d278b441" ;;
        1.9.4:mac:macaarch64) jul_256="67542975e86102eec95bc4bb7c30c5d8c7ea9f9a0b388f0e10f546945363b01a" ;;

        1.9.3:linux:x86_64)   jul_256="d76670cc9ba3e0fd4c1545dd3d00269c0694976a1176312795ebce1692d323d1" ;;
        1.9.3:linux:aarch64)  jul_256="55437879f6b98470d96c4048b922501b643dfffb8865abeb90c7333a83df7524" ;;
        1.9.3:mac:mac64)      jul_256="6eea87748424488226090d1e7d553e72ab106a873d63c732fc710a3d080abb97" ;;
        1.9.3:mac:macaarch64) jul_256="f518e38d7bd5b37766fb051916bd295993aa4b52a47018f4c98b5fde721ced87" ;;

        1.9.2:linux:x86_64)   jul_256="4c2d799f442d7fe718827b19da2bacb72ea041b9ce55f24eee7b1313f57c4383" ;;
        1.9.2:linux:aarch64)  jul_256="682397f8895149f0e283f0b27bffc6694033bdfb19f9366c80f6efdf3685f27c" ;;
        1.9.2:mac:mac64)      jul_256="a2e8eb31a89b26e4a99349303aeff8e8ee780144bbdb1f7eda6f41024d42cadb" ;;
        1.9.2:mac:macaarch64) jul_256="77c71ff8cb1fcdb84097e86a9fb579a8b34d8e7fd8e24d43107042e0fb988b76" ;;

        1.9.1:linux:x86_64)   jul_256="cde14a58f899251f30cfced87055626f44845780659ebe8d50cbc4c67b31997c" ;;
        1.9.1:linux:aarch64)  jul_256="b643ccd3e2a5960f7ce7055243743d0a39badda3974bce3d77861dd363badd10" ;;
        1.9.1:mac:mac64)      jul_256="49368ddaef4e37ed606808a6ca58ba0f1a4451a27b8201aed4d9d7b24c276817" ;;
        1.9.1:mac:macaarch64) jul_256="9e3e02ca6546513dce265379abe957cb2b5b0ccf4066219486da0eb872ddcebc" ;;

        1.9.0:linux:x86_64)   jul_256="00c614466ef9809c2eb23480e38d196a2c577fff2730c4f83d135b913d473359" ;;
        1.9.0:linux:aarch64)  jul_256="0a14315b53acd97f22d26d4a8fd2c237e524e95c3bec98d2d78b54d80c2bc364" ;;
        1.9.0:mac:mac64)      jul_256="00bc4c27eeb1bebd350597312ff0919176315fd3199c63ec963fb41e3b04bfaf" ;;
        1.9.0:mac:macaarch64) jul_256="53e62770a6990d5a89e7a001ef68aa25de25126a3be838200c4c9a705daea37c" ;;

        1.8.5:linux:x86_64)   jul_256="e71a24816e8fe9d5f4807664cbbb42738f5aa9fe05397d35c81d4c5d649b9d05" ;;
        1.8.5:linux:aarch64)  jul_256="a1f637b44c71ea9bc96d7c3ef347724c054a1e5227b980adebfc33599e5153a4" ;;
        1.8.5:mac:mac64)      jul_256="a1a859eda7fb41a0b55467339a11c3c1c0df78b27d1e160e80bc6758b3d8dae0" ;;
        1.8.5:mac:macaarch64) jul_256="ea85e0489c36324c4da62163aa1b82fcf2f52f72d173ee7dd213a3a92992cab7" ;;

        1.8.4:linux:x86_64)   jul_256="f0427a4d7910c47dc7c31f65ba7ecaafedbbc0eceb39c320a37fa33598004fd5" ;;
        1.8.4:linux:aarch64)  jul_256="dc4798c1ce8768fa35972e8b149ca3a85fc69e1074b609a72b2cfed5c4aa7050" ;;
        1.8.4:mac:mac64)      jul_256="597d4ec4f12241e78c75e220cebd455fe2af935a36f276d222a419616553663a" ;;
        1.8.4:mac:macaarch64) jul_256="06e81151d76ccd5ec0bd59cd51dde49cc1a15b1386624b4a61557cdbc5ae6d09" ;;

        1.8.3:linux:x86_64)   jul_256="33c3b09356ffaa25d3331c3646b1f2d4b09944e8f93fcb994957801b8bbf58a9" ;;
        1.8.3:linux:aarch64)  jul_256="dbffb134a413b712d4a8e1ee8e665ea55edb0865719a1bad9979123d6433acc9" ;;
        1.8.3:mac:mac64)      jul_256="f3367de05f681b884219854c304f1a85420000c8c137a98cd358b7fe5070dc84" ;;
        1.8.3:mac:macaarch64) jul_256="f57acd021e7e7c0a40d29967a0f680ca77c5988d856e1cc220982c6bfa3964ff" ;;

        1.8.2:linux:x86_64)   jul_256="671cf3a450b63a717e1eedd7f69087e3856f015b2e146cb54928f19a3c05e796" ;;
        1.8.2:linux:aarch64)  jul_256="f91c276428ffb30acc209e0eb3e70b1c91260e887e11d4b66f5545084b530547" ;;
        1.8.2:mac:mac64)      jul_256="7395f1c49e3c4dbc3714aef2c6cf310addd4c1072a12a8c4ca7f568c67acb15d" ;;
        1.8.2:mac:macaarch64) jul_256="a4fc5caa5bf2ac353f779c177a9a72ee91497a092959d58a16263041a68b92d3" ;;

        1.8.1:linux:x86_64)   jul_256="33054ee647ee8a4fb54fc05110e07e0b53e04591fe53d0a4cb4c7ed7a05e91f1" ;;
        1.8.1:linux:aarch64)  jul_256="ba06837ac2899547bbb799989f11464fecd6782226871c3b7a48619481042679" ;;
        1.8.1:mac:mac64)      jul_256="a1cc8dbd2a02b7ef436ca1450bab831a68b74ce2431c2e1043dc3354780e0bb2" ;;
        1.8.1:mac:macaarch64) jul_256="5a0f49416fd40760cd188d90db9742087fd3d2f86e725b5d31ed8ce91a2331ce" ;;

        1.8.0:linux:x86_64)   jul_256="e80d732ccb7f79e000d798cb8b656dc3641ab59516d6e4e52e16765017892a00" ;;
        1.8.0:linux:aarch64)  jul_256="e003cfb8680af1a65c3be55b53a48cc5186300adaaba8926209800b4d1f4ca7a" ;;
        1.8.0:mac:mac64)      jul_256="a77055e5005d05d43fce4dc51e78f664b4802138a432ab84c0adecc453dc88a5" ;;
        1.8.0:mac:macaarch64) jul_256="9c911f93405445ee71f6ebfef2b00f9ed9d4880b4bfa6c36fe865b75c2a46fbd" ;;
    esac
}


function check_julia_version() {
    case "${v_julia}" in
        1.8.5)
            return 0
            ;;

        1.8.[0-4])
            echo_warn \
                "Julia ${v_julia} was requested. The protocol-tested/default" \
                "version is 1.8.5."
            return 0
            ;;

        1.9.[0-4])
            echo_warn \
                "Julia ${v_julia} was requested. Upstream Atria docs report" \
                "Julia 1.8.5 is faster than Julia 1.9.x; this protocol" \
                "defaults to Julia 1.8.5."
            return 0
            ;;

        *)
            echo_err \
                "unsupported Julia version ${v_julia}. This installer" \
                "currently supports Julia 1.8.0-1.8.5 and 1.9.0-1.9.4," \
                "with Julia 1.8.5 as the protocol-tested/default version."
            return 1
            ;;
    esac
}


#  Intentionally sets global Atria tag variables used by later installer steps
function map_atria_tag() {
    if [[ "${v_atria}" == v* ]]; then
        tag_atr="${v_atria}"
        v_atria="${v_atria#v}"
    else
        tag_atr="v${v_atria}"
    fi
}


function check_req_cmd() {
    local cmd=""
    local rc=0

    #  Use 'curl' here because it is available by default on macOS and common
    #+ on Linux
    for cmd in cp curl find git grep head mkdir mktemp rm sort tail; do
        if ! command -v "${cmd}" > /dev/null 2>&1; then
            echo_err "required command '${cmd}' is not available in PATH."
            rc=1
        fi
    done

    if ! command -v tar > /dev/null 2>&1; then
        echo_err \
            "required command 'tar' is not available in PATH. This script" \
            "uses 'tar' to extract the Julia archive. BSD tar, which is" \
            "commonly installed by default on macOS systems, is acceptable;" \
            "extraction is verified by checking the expected Julia" \
            "executable afterward."
        rc=1
    fi

    return "${rc}"
}


function prepare_tmp_dir() {
    if [[ -n "${dir_tmp}" ]]; then
        mkdir -p "${dir_tmp}"
        tmp_auto=false
    else
        dir_tmp="$(mktemp -d "${TMPDIR:-/tmp}/install_atria.XXXXXX")"
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


function install_julia() {
    local pth_path_jul=""

    #  Allow global variable initialization here
    jul_dir="${dir_install}/julia-${v_julia}"
    jul_bin="${jul_dir}/bin/julia"
    jul_tar_pth="${dir_tmp}/${jul_tar}"

    if [[ -d "${jul_dir}" ]]; then
        case "${if_exists}" in
            fail)
                echo_err \
                    "Julia install directory already exists: '${jul_dir}'."
                echo_err \
                    "Nothing was changed. To reuse it, rerun with '--if_exists" \
                    "reuse'."
                return 1
                ;;

            reuse)
                echo \
                    "Julia install directory already exists; verifying and" \
                    "reusing '${jul_dir}'." >&2
                verify_julia_exec "${jul_bin}" || return 1
                return 0
                ;;
        esac
    fi

    if [[ "${if_exists}" == "reuse" ]]; then
        if pth_path_jul="$(command -v julia 2>/dev/null)"; then
            if \
                verify_julia_exec "${pth_path_jul}"
            then
                jul_bin="${pth_path_jul}"
                jul_dir="$(cd "$(dirname "${jul_bin}")/.." && pwd)"
                echo \
                    "Matching Julia executable found on PATH; reusing" \
                    "'${jul_bin}'." >&2
                return 0
            fi

            echo_warn \
                "Julia was found on PATH, but it does not match the" \
                "requested version '${v_julia}'; installing Julia under" \
                "'${dir_install}'."
        fi
    fi

    if [[ -z "${jul_256}" ]]; then
        echo_err \
            "no bundled SHA-256 mapping for Julia ${v_julia} on" \
            "'${sys_os}:${sys_ar}'. Refusing to download without a" \
            "known checksum."
        return 1
    fi

    if [[ ! -s "${jul_tar_pth}" ]]; then
        run_or_print curl -L --fail -o "${jul_tar_pth}" "${jul_url}"
    else
        echo "Julia tarball already exists; verifying '${jul_tar_pth}'."
    fi

    verify_sha256 "${jul_tar_pth}" "${jul_256}"
    extract_julia_tar "${jul_tar_pth}" "${dir_install}"

    verify_julia_exec "${jul_bin}" || return 1
}


function checkout_atria() {
    local pth_path_atr=""

    dir_prg="${dir_install}"
    dir_atr="${dir_prg}/Atria"

    if [[ "${if_exists}" == "reuse" && ! -e "${dir_atr}" ]]; then
        if pth_path_atr="$(command -v atria 2>/dev/null)"; then
            if verify_atria_exec "${pth_path_atr}"; then
                pth_atr="${pth_path_atr}"
                pth_bin="$(dirname "${pth_atr}")"
                echo \
                    "Matching Atria executable found on PATH; skipping" \
                    "repository checkout for '${dir_atr}'." >&2
                return 0
            fi
        fi
    fi

    mkdir -p "${dir_prg}"

    if [[ -d "${dir_atr}/.git" ]]; then
        case "${if_exists}" in
            fail)
                echo_err \
                    "Atria repository already exists: '${dir_atr}'."
                echo_err \
                    "Nothing was changed. To reuse it, rerun with '--if_exists" \
                    "reuse'."
                return 1
                ;;

            reuse)
                echo \
                    "Atria repository already exists; verifying and reusing" \
                    "'${dir_atr}'." >&2
                ;;
        esac
    elif [[ -e "${dir_atr}" ]]; then
        echo_err \
            "Atria install path exists but is not a Git repository:" \
            "'${dir_atr}'. Move it aside or choose another '--dir_install'."
        return 1
    else
        run_or_print \
            git clone "https://github.com/cihga39871/Atria.git" "${dir_atr}"
    fi

    if [[ "${dry_run}" == "true" ]]; then
        echo_dry "would check out Atria tag '${tag_atr}'."
        return 0
    fi

    (
        cd "${dir_atr}"
        git fetch --tags

        if \
            git rev-parse --verify --quiet "refs/tags/${tag_atr}" > /dev/null
        then
            git checkout --detach "${tag_atr}"
        else
            echo_err "could not find Atria tag '${tag_atr}'."
            return 1
        fi
    )
}


function copy_suitesparse() {
    local dir_src="${jul_dir}/lib/julia"
    local dir_dst="${dir_bld}/lib/julia"
    local lib=""
    local stem=""

    if [[ ! -d "${dir_dst}" ]]; then
        echo_err \
            "could not find Atria 'lib/julia' directory: '${dir_dst}'."
        return 1
    fi

    for stem in \
        libamd \
        libbtf \
        libcamd \
        libccolamd \
        libcholmod \
        libcolamd \
        libklu \
        libldl \
        librbio \
        libspqr \
        libsuitesparseconfig \
        libumfpack
    do
        for lib in "${dir_src}/${stem}".*; do
            [[ -e "${lib}" ]] || continue
            cp -p "${lib}" "${dir_dst}/"
        done
    done
}


function find_atria_dir() {
    dir_bld="$(
        find "${dir_atr}" \
            -maxdepth 1 \
            -type d \
            \( \
                -name "atria-${v_atria}" \
                -o -name "atria-${v_atria}-*" \
                -o -name "atria-${v_atria}_*" \
            \) \
            -print \
            | sort \
            | tail -n 1
    )"

    if [[ -z "${dir_bld}" ]]; then
        dir_bld="${dir_atr}/atria-${v_atria}"
    fi

    if [[ ! -d "${dir_bld}" ]]; then
        echo_err \
            "failed to locate Atria build directory under '${dir_atr}'."
        return 1
    fi
}


function has_atria_library_error() {
    local file="${1:-}"
    local patn="cholmod|suitesparse|libcholmod|libamd"

    patn="${patn}|could not load library|library not loaded"

    validate_var_file "file" "${file}" || return 1

    grep -qiE "${patn}" "${file}"
}


function verify_atria_exec() {
    local pth_atr="${1:-}"
    local log_chk="${2:-}"
    local tmp_hlp=""
    local tmp_ver=""
    local use_tmp_log=false

    validate_var_file "pth_atr" "${pth_atr}" || return 1

    if [[ -z "${log_chk}" ]]; then
        log_chk="$(mktemp "${TMPDIR:-/tmp}/atria_check.XXXXXX")"
        use_tmp_log=true
    else
        : > "${log_chk}"
    fi

    tmp_ver="$(mktemp "${TMPDIR:-/tmp}/atria_version.XXXXXX")"
    tmp_hlp="$(mktemp "${TMPDIR:-/tmp}/atria_help.XXXXXX")"

    if ! "${pth_atr}" --version > "${tmp_ver}" 2>&1; then
        cat "${tmp_ver}" > "${log_chk}"
        cat "${tmp_ver}" >&2
        rm -f "${tmp_ver}" "${tmp_hlp}"
        [[ "${use_tmp_log}" == "true" ]] && rm -f "${log_chk}"
        echo_err "Atria '--version' failed for '${pth_atr}'."
        return 1
    fi

    cat "${tmp_ver}" > "${log_chk}"

    if \
        grep -qi 'error' "${tmp_ver}" \
        || has_atria_library_error "${tmp_ver}"
    then
        cat "${tmp_ver}" >&2
        rm -f "${tmp_ver}" "${tmp_hlp}"
        [[ "${use_tmp_log}" == "true" ]] && rm -f "${log_chk}"
        echo_err "Atria emitted an error during executable verification."
        return 1
    fi

    if ! grep -Fq "v${v_atria}" "${tmp_ver}"; then
        cat "${tmp_ver}" >&2
        rm -f "${tmp_ver}" "${tmp_hlp}"
        [[ "${use_tmp_log}" == "true" ]] && rm -f "${log_chk}"
        echo_err \
            "Atria executable exists at '${pth_atr}', but it does not" \
            "appear to match requested version '${v_atria}'."
        return 1
    fi

    if ! "${pth_atr}" --help > "${tmp_hlp}" 2>&1; then
        {
            echo
            cat "${tmp_hlp}"
        } >> "${log_chk}"
        cat "${tmp_hlp}" >&2
        rm -f "${tmp_ver}" "${tmp_hlp}"
        [[ "${use_tmp_log}" == "true" ]] && rm -f "${log_chk}"
        echo_err "Atria '--help' failed for '${pth_atr}'."
        return 1
    fi

    if has_atria_library_error "${tmp_hlp}"; then
        {
            echo
            cat "${tmp_hlp}"
        } >> "${log_chk}"
        cat "${tmp_hlp}" >&2
        rm -f "${tmp_ver}" "${tmp_hlp}"
        [[ "${use_tmp_log}" == "true" ]] && rm -f "${log_chk}"
        echo_err "Atria emitted a library error during help verification."
        return 1
    fi

    cat "${tmp_ver}"
    rm -f "${tmp_ver}" "${tmp_hlp}"
    [[ "${use_tmp_log}" == "true" ]] && rm -f "${log_chk}"
    return 0
}


function retry_suitesparse_fallback() {
    local log_chk="${1:-}"
    local stem=""

    validate_var_file "log_chk" "${log_chk}" || return 1

    if [[ "${sys_os}" != "Darwin" ]]; then
        return 1
    fi

    if ! has_atria_library_error "${log_chk}"; then
        return 1
    fi

    echo_warn \
        "Atria verification reported a macOS SuiteSparse/CHOLMOD library" \
        "error. Copying Julia SuiteSparse libraries into the Atria build and" \
        "retrying verification."

    copy_suitesparse || return 1

    for stem in libamd libcholmod libsuitesparseconfig; do
        if ! \
            find "${dir_bld}/lib/julia" \
                -name "${stem}*" \
                -print \
                -quit \
                | grep -q .
        then
            echo_err \
                "failed to copy '${stem}' libraries into the Atria build."
            return 1
        fi
    done
}


function build_atria() {
    local log_chk=""

    if [[ "${dry_run}" == "true" ]]; then
        echo_dry \
            "would build Atria with '${jul_bin}' build_atria.jl."
        pth_bin="${dir_atr}/atria-${v_atria}/bin"
        return 0
    fi

    if \
        find_atria_dir > /dev/null 2>&1
    then
        pth_atr="${dir_bld}/bin/atria"

        if [[ -f "${pth_atr}" ]]; then
            pth_bin="$(dirname "${pth_atr}")"

            if \
                verify_atria_exec "${pth_atr}"
            then
                echo \
                    "Atria executable already exists and verifies; reusing" \
                    "'${pth_atr}'."
                return 0
            fi

            echo_warn \
                "existing Atria executable failed verification; rebuilding" \
                "Atria under '${dir_atr}'."
        fi
    fi

    if [[ "${if_exists}" == "reuse" ]]; then
        if \
            pth_atr="$(command -v atria 2>/dev/null)"
        then
            if \
                verify_atria_exec "${pth_atr}"
            then
                pth_bin="$(dirname "${pth_atr}")"
                echo \
                    "Matching Atria executable found on PATH; reusing" \
                    "'${pth_atr}'."
                return 0
            fi

            echo_warn \
                "Atria was found on PATH, but it did not verify as version" \
                "'${v_atria}'; rebuilding Atria under '${dir_atr}'."
        fi
    fi

    (
        cd "${dir_atr}"
        "${jul_bin}" build_atria.jl
    )

    find_atria_dir

    pth_atr="${dir_bld}/bin/atria"

    if [[ ! -f "${pth_atr}" ]]; then
        echo_err \
            "failed to locate built Atria executable under '${dir_atr}'."
        return 1
    fi

    pth_bin="$(dirname "${pth_atr}")"
    log_chk="$(mktemp "${TMPDIR:-/tmp}/atria_check.XXXXXX")"

    if verify_atria_exec "${pth_atr}" "${log_chk}"; then
        rm -f "${log_chk}"
        return 0
    fi

    if retry_suitesparse_fallback "${log_chk}"; then
        rm -f "${log_chk}"
        verify_atria_exec "${pth_atr}" || return 1
        return 0
    fi

    rm -f "${log_chk}"
    return 1
}


function print_write_path_snippet() {
    local snippet

    snippet="$(
        cat << EOM
#  Julia
export PATH="\${PATH}:${jul_dir}/bin"

#  Atria
export PATH="\${PATH}:${pth_bin}"
EOM
    )"

    echo
    echo "Add the following lines to your shell configuration or a sourced snippet:"
    echo
    printf '%s\n' "${snippet}"

    if [[ -n "${path_snippet}" ]]; then
        if [[ "${dry_run}" == "true" ]]; then
            echo
            echo_dry \
                "would append PATH snippet to '${path_snippet}'."
        else
            mkdir -p "$(dirname "${path_snippet}")"

            if \
                   [[ -f "${path_snippet}" ]] \
                && grep -Fq "${jul_dir}/bin" "${path_snippet}" \
                && grep -Fq "${pth_bin}" "${path_snippet}"
            then
                echo
                echo "PATH snippet already appears to exist in '${path_snippet}'."
            else
                {
                    echo
                    printf '%s\n' "${snippet}"
                } >> "${path_snippet}"
                echo
                echo "PATH snippet appended to '${path_snippet}'."
            fi
        fi
    fi
}


function parse_args() {
    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        help_install_atria >&2
        return 2
    fi

    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
            -h|--hlp|--help)
                help_install_atria
                return 2
                ;;

            -dr|--dry|--dry[_-]run)
                dry_run=true
                shift 1
                ;;

            -en|--env|--env[_-]nam)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_install_atria
                    return 1
                }
                env_nam="${2}"
                shift 2
                ;;

            -di|--dir[_-]install)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_install_atria
                    return 1
                }
                dir_install="${2}"
                shift 2
                ;;

            -dt|--tmp|--dir[_-]tmp)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_install_atria
                    return 1
                }
                dir_tmp="${2}"
                shift 2
                ;;

            -vj|--v[_-]julia|--julia[_-]version)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_install_atria
                    return 1
                }
                v_julia="${2}"
                shift 2
                ;;

            -va|--v[_-]atria|--atria[_-]version)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_install_atria
                    return 1
                }
                v_atria="${2}"
                shift 2
                ;;

            -ps|--path[_-]snippet)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_install_atria
                    return 1
                }
                path_snippet="${2}"
                shift 2
                ;;

            -ie|--if[_-]exists)
                require_optarg "${1}" "${2:-}" "main" || {
                    echo >&2
                    help_install_atria
                    return 1
                }
                if_exists="${2}"
                shift 2
                ;;

            *)
                echo_err "unknown option/parameter passed: '${1}'."
                echo >&2
                help_install_atria
                return 1
                ;;
        esac
    done
}


function validate_args() {
    validate_var "env_nam" "${env_nam}" || return 1
    validate_var "dir_install" "${dir_install}" || return 1
    validate_var "v_julia" "${v_julia}" || return 1
    validate_var "v_atria" "${v_atria}" || return 1
    validate_var "if_exists" "${if_exists}" || return 1

    case "${if_exists}" in
        fail|reuse) : ;;
        *)
            echo_err \
                "invalid '--if_exists' value: '${if_exists}'. Must be 'fail'" \
                "or 'reuse'."
            return 1
            ;;
    esac

    if [[ "${dry_run}" == "false" && ! -d "${dir_install}" ]]; then
        mkdir -p "${dir_install}"
    fi

    if [[ "${dry_run}" == "true" ]]; then
        dir_install="${dir_install/#\~/${HOME}}"
    else
        validate_var_dir "dir_install" "${dir_install}" || return 1
        dir_install="$(cd "${dir_install}" > /dev/null 2>&1 && pwd)"
    fi

    check_julia_version
    map_atria_tag
    map_julia_target
}


function report_plan() {
    if [[ "${dry_run}" == "true" ]]; then
        echo_dry \
            "no downloads, extraction, cloning, building, or PATH snippet" \
            "writes will be performed."
    fi

    echo
    echo "Resolved installation plan:"
    echo "  - env_nam=${env_nam}"
    echo "  - dir_install=${dir_install}"
    echo "  - dir_tmp=${dir_tmp:-AUTO}"
    echo "  - sys_os=${sys_os}"
    echo "  - sys_ar=${sys_ar}"
    echo "  - v_julia=${v_julia}"
    echo "  - v_min=${v_min}"
    echo "  - jul_os=${jul_os}"
    echo "  - jul_ar_dir=${jul_ar_dir}"
    echo "  - jul_ar_fil=${jul_ar_fil}"
    echo "  - jul_tar=${jul_tar}"
    echo "  - jul_url=${jul_url}"
    echo "  - jul_256=${jul_256:-UNSET}"
    echo "  - v_atria=${v_atria}"
    echo "  - tag_atr=${tag_atr}"
    echo "  - dir_atr=${dir_install}/Atria"
    echo "  - path_snippet=${path_snippet:-UNSET}"
    echo "  - if_exists=${if_exists}"
    echo
}


function check_checksum_mapping() {
    if [[ -n "${jul_256}" ]]; then
        return 0
    fi

    if [[ "${dry_run}" == "true" ]]; then
        echo_warn \
            "no bundled SHA-256 mapping was found for Julia ${v_julia} on" \
            "'${sys_os}:${sys_ar}'. A non-dry run will refuse to download" \
            "Julia until the checksum mapping is added."
        return 0
    fi

    echo_err \
        "no bundled SHA-256 mapping was found for Julia ${v_julia} on" \
        "'${sys_os}:${sys_ar}'. Refusing to download Julia without" \
        "a known checksum."
    return 1
}


function check_runtime_req() {
    if [[ "${dry_run}" == "true" ]]; then
        echo_dry "would confirm 'mamba' or 'conda' is available."
        echo_dry \
            "would confirm environment '${env_nam}' exists and activate it."
        echo_dry \
            "would confirm 'pigz', 'pbzip2', and 'Rscript' are available in" \
            "'${env_nam}'."
        echo_dry \
            "would confirm 'curl', 'git', 'grep', 'sort', and 'tar', among" \
            "other programs, are available."
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
                "Install or update the project environment before installing" \
                "Atria."
            return 1
        fi

        handle_env "${env_nam}" || return 1
    fi

    check_pgrm_path pigz || {
        echo_err \
            "'pigz' is missing from '${env_nam}'. Install/update the project" \
            "environment; this script does not install pigz with a system" \
            "package manager."
        return 1
    }
    check_pgrm_path pbzip2 || {
        echo_err \
            "'pbzip2' is missing from '${env_nam}'. Install/update the" \
            "project environment; this script does not install pbzip2 with" \
            "a system package manager."
        return 1
    }
    check_pgrm_path Rscript || {
        echo_err \
            "'Rscript' is missing from '${env_nam}'. Install/update the" \
            "project environment; this script expects R support from the" \
            "active project environment."
        return 1
    }
    check_req_cmd || return 1
    select_sha_cmd || return 1
}


function report_dry_run_install() {
    jul_dir="${dir_install}/julia-${v_julia}"
    jul_bin="${jul_dir}/bin/julia"
    dir_prg="${dir_install}"
    dir_atr="${dir_prg}/Atria"
    pth_bin="${dir_atr}/atria-${v_atria}/bin"

    echo_dry "would use working directory '${dir_tmp:-AUTO}'."

    if [[ -n "${jul_256}" ]]; then
        echo_dry "would download and verify '${jul_url}'."
        echo_dry \
            "would extract Julia under '${dir_install}', or reuse matching" \
            "Julia from '${jul_dir}' or PATH if '--if_exists reuse' was" \
            "specified."
    else
        echo_dry \
            "would refuse non-dry-run download of '${jul_url}' until a" \
            "SHA-256 mapping is added."
    fi

    echo_dry \
        "would discover the final Atria bin path after build; provisional" \
        "path is '${pth_bin}'."
    echo_dry \
        "would clone Atria under '${dir_atr}', or reuse a matching existing" \
        "Atria install under '${dir_atr}' or on PATH if '--if_exists reuse'" \
        "was specified."
    echo_dry "would check out Atria tag '${tag_atr}'."
    echo_dry "would build Atria using '${jul_bin}'."
    print_write_path_snippet
}


function run_install() {
    prepare_tmp_dir
    trap cleanup_tmp_dir EXIT

    install_julia
    checkout_atria
    build_atria
    print_write_path_snippet

    echo
    echo "success($(basename "${BASH_SOURCE[0]}")):" \
        "installed/verified Julia ${v_julia} and Atria ${v_atria}."
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
    report_plan
    check_checksum_mapping
    check_runtime_req

    if [[ "${dry_run}" == "true" ]]; then
        report_dry_run_install
        return 0
    fi

    run_install
}


main "$@"
