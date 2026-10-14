#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: download_fastqs.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# The following were used in design, development, and documentation, with all
# output reviewed, edited, and approved by the author:
# - OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
#   GPT-5.6);
# - Anthropic Claude Code (Opus 5.5).
#
# Distributed under the MIT license.


# derive_fastq_path
# verify_fastq
# download_fastq
# print_fastq_names_form
# restore_fastq_names
# link_fastqs
# download_fastqs


# Require Bash >= 4.4 before defining functions.
if [[ -z "${BASH_VERSION:-}" ]]; then
    echo "error(shell):" \
        "this script must be sourced or run under Bash >= 4.4." >&2

    if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
        return 1
    else
        exit 1
    fi
elif ((
    BASH_VERSINFO[0] < 4 || ( BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] < 4 )
)); then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "this script requires Bash >= 4.4; current version is" \
        "'${BASH_VERSION}'." >&2

    if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
        return 1
    else
        exit 1
    fi
fi

# Source required helper functions if needed.
{
    _dir_src_dl="$(
        cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
    )"

    # shellcheck source=lib/bash/core/source_helpers.sh
    source "${_dir_src_dl}/../core/source_helpers.sh" || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source '${_dir_src_dl}/../core/source_helpers.sh'." >&2

        if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
            return 1
        else
            exit 1
        fi
    }

    source_helpers "${_dir_src_dl}/.." \
        check_args check_env check_inputs format_outputs || {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to source required helper dependencies." >&2

            if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
                return 1
            else
                exit 1
            fi
        }

    unset _dir_src_dl
}


# Print the path of one accession- or custom-named FASTQ file.
function derive_fastq_path() {
    local dir="${1:-}"
    local stem="${2:-}"
    local mate="${3:-}"
    local show_help

    show_help=$(cat << EOM
Usage
-----
  derive_fastq_path
    [--help] dir stem mate

  Print the path of one FASTQ file named for a run accession or a custom name.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  dir : dir
    Directory that holds or will hold the file.

  2  stem : str
    Run accession or custom name.

  3  mate : str
    '1' or '2' for a paired-end mate, or 'NA' for single-end data.

Returns
-------
  Prints '<dir>/<stem>_R1.fastq.gz', '<dir>/<stem>_R2.fastq.gz', or '<dir>/<stem>.fastq.gz' and returns 0; returns 1 when an argument is missing or 'mate' is not '1', '2', or 'NA'.

Notes
-----
  Runtime requirements:
    bash >= 4.4

Examples
--------
  1. Print the path of a paired-end second mate.
    '''bash
    derive_fastq_path "\${dir_out}" SRR32203047 2
    '''

  2. Print the path of a single-end custom-named link.
    '''bash
    derive_fastq_path "\${dir_sym}" IP_WT_Q_Brn1_rep2 NA
    '''
EOM
    )

    if [[ "${dir}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${dir}" || -z "${stem}" || -z "${mate}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "expected three arguments: 'dir', 'stem', and 'mate'."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    case "${mate}" in
        1|2) printf '%s\n' "${dir}/${stem}_R${mate}.fastq.gz" ;;
        NA)  printf '%s\n' "${dir}/${stem}.fastq.gz"          ;;
        *)
            echo_err_func "${FUNCNAME[0]}" \
                "'mate' must be '1', '2', or 'NA': '${mate}'."
            return 1
            ;;
    esac
}


# Check a downloaded file against its expected size and MD5; without an MD5,
# test its gzip integrity. Print the reason for a failure.
function verify_fastq() {
    local fil="${1:-}"
    local md5="${2:-}"
    local siz="${3:-}"
    local val
    local show_help

    show_help=$(cat << EOM
Usage
-----
  verify_fastq
    [--help] fil md5 siz

  Check one downloaded FASTQ file against its expected size and MD5, or test its gzip integrity when no MD5 is given.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  fil : file
    Downloaded file to check.

  2  md5 : str
    Expected MD5, or 'NA' to test gzip integrity instead.

  3  siz : str
    Expected size in bytes, or 'NA' to skip the size check.

Returns
-------
  Returns 0 when every check passes; otherwise prints the reason to stdout and returns 1.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - gzip (when 'md5' is 'NA')
    - md5sum (when 'md5' is not 'NA')
    - wc

  - The size is checked before the MD5. MD5 values are compared without regard to letter case.

Examples
--------
  1. Check a download against the MD5 and size from an ENA filereport.
    '''bash
    verify_fastq \\
        "\${dir_out}/SRR32203047_R1.fastq.gz.part" \\
        6c0bc341f73b0ac0d345fa62cb80a044 \\
        270090667
    '''

  2. Test the gzip integrity of a download without an MD5.
    '''bash
    verify_fastq "\${dir_out}/SRR7175377.fastq.gz.part" NA NA
    '''
EOM
    )

    if [[ "${fil}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${fil}" || -z "${md5}" || -z "${siz}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "expected three arguments: 'fil', 'md5', and 'siz'."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    if [[ "${siz}" != "NA" ]]; then
        val="$(wc -c < "${fil}")"
        val="${val//[[:space:]]/}"

        if [[ "${val}" != "${siz}" ]]; then
            echo "size is ${val} bytes, expected ${siz}"
            return 1
        fi
    fi

    if [[ "${md5}" != "NA" ]]; then
        val="$(md5sum < "${fil}")"
        val="${val%% *}"

        if [[ "${val,,}" != "${md5,,}" ]]; then
            echo "MD5 is ${val}, expected ${md5,,}"
            return 1
        fi
    elif ! \
        gzip -t < "${fil}" 2> /dev/null
    then
        echo "file fails 'gzip -t'"
        return 1
    fi
}


# Download one FASTQ URL to a partial file, verify it, and move it into place;
# a file that fails verification is downloaded again, up to 'n_try' times.
function download_fastq() {
    local acc="${1:-}"
    local url="${2:-}"
    local fil_out="${3:-}"
    local md5="${4:-}"
    local siz="${5:-}"
    local log_out="${6:-}"
    local log_err="${7:-}"
    local fil_prt n_try i_try reason
    local show_help
    local -a arr_cmd_wget

    show_help=$(cat << EOM
Usage
-----
  download_fastq
    [--help] acc url fil_out md5 siz log_out log_err

  Download one FASTQ file to a partial file ('.part' appended to 'fil_out'), verify it, and move it to 'fil_out'.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  acc : str
    Run accession, used in progress messages.

  2  url : str
    FASTQ URL; a URL without a scheme is downloaded over 'https://'.

  3  fil_out : file
    Output FASTQ file.

  4  md5 : str
    Expected MD5, or 'NA'.

  5  siz : str
    Expected size in bytes, or 'NA'.

  6  log_out : file
    File that receives the stdout of 'wget'.

  7  log_err : file
    File that receives the stderr of 'wget'.

Returns
-------
  Returns 0 when the file is downloaded, verified, and moved into place; otherwise removes the partial file, prints an error, and returns 1.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - gzip (when 'md5' is 'NA')
    - md5sum (when 'md5' is not 'NA')
    - mv
    - Network access
    - rm
    - wc
    - wget

  - 'wget' makes up to 5 tries, including when the connection is refused. A failed download is not attempted again.
  - A file that fails 'verify_fastq()' is downloaded again, up to 3 attempts in all.
  - An existing 'fil_out' is replaced only by a file that passes verification.

Examples
--------
  1. Download and verify one paired-end mate from ENA.
    '''bash
    download_fastq \\
        SRR32203047 \\
        "ftp.sra.ebi.ac.uk/vol1/fastq/SRR322/047/SRR32203047/SRR32203047_1.fastq.gz" \\
        "\${dir_out}/SRR32203047_R1.fastq.gz" \\
        6c0bc341f73b0ac0d345fa62cb80a044 \\
        270090667 \\
        "\${dir_eo}/SRR32203047_R1.stdout.txt" \\
        "\${dir_eo}/SRR32203047_R1.stderr.txt"
    '''

  2. Download one single-end file without an MD5, checking only its gzip integrity.
    '''bash
    download_fastq \\
        SRR7175377 \\
        "\${url_se}" \\
        "\${dir_out}/SRR7175377.fastq.gz" \\
        NA \\
        NA \\
        "\${dir_eo}/SRR7175377.stdout.txt" \\
        "\${dir_eo}/SRR7175377.stderr.txt"
    '''
EOM
    )

    if [[ "${acc}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif (( $# != 7 )); then
        echo_err_func "${FUNCNAME[0]}" \
            "expected seven arguments: 'acc', 'url', 'fil_out', 'md5'," \
            "'siz', 'log_out', and 'log_err'."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    fil_prt="${fil_out}.part"
    n_try=3

    # ENA's 'fastq_ftp' values carry no scheme, and 'wget' would then use HTTP
    # on port 80.
    if [[ "${url}" != *://* ]]; then
        url="https://${url}"
    fi

    arr_cmd_wget=(
        wget
            --progress=dot:mega
            --tries=5
            --retry-connrefused
            --waitretry=10
            -O "${fil_prt}"
            "${url}"
    )

    for (( i_try = 1; i_try <= n_try; i_try++ )); do
        echo \
            "Downloading ${acc} from ${url} (attempt ${i_try} of" \
            "${n_try})." >&2
        printf '%q ' "${arr_cmd_wget[@]}" >&2
        echo >&2

        if ! \
            "${arr_cmd_wget[@]}" > "${log_out}" 2> "${log_err}"
        then
            rm -f -- "${fil_prt}"
            echo_err "failed to download '${url}'; see '${log_err}'."
            return 1
        fi

        if \
            reason="$(verify_fastq "${fil_prt}" "${md5}" "${siz}")"
        then
            mv -f -- "${fil_prt}" "${fil_out}" || {
                rm -f -- "${fil_prt}"
                echo_err "failed to move '${fil_prt}' to '${fil_out}'."
                return 1
            }
            return 0
        fi

        if (( i_try < n_try )); then
            echo_warn \
                "download of '${url}' failed verification: ${reason};" \
                "downloading it again."
        fi
    done

    rm -f -- "${fil_prt}"
    echo_err \
        "download of '${url}' failed verification ${n_try} times: ${reason}."
    return 1
}


# Print the read-name form of a FASTQ file from its first header: 'archive'
# when the name needs restoring, 'lost' when the archive did not store it,
# 'restored' when it was restored before, 'other' otherwise.
function print_fastq_names_form() {
    local fil="${1:-}"
    local acc="${2:-}"
    local hdr nam cmt
    local show_help

    show_help=$(cat << EOM
Usage
-----
  print_fastq_names_form
    [--help] fil acc

  Print the read-name form of one gzipped FASTQ file, judged from its first header.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  fil : file
    Gzipped FASTQ file.

  2  acc : str
    Run accession of the file.

Returns
-------
  Prints one word to stdout and returns 0:
    - 'archive': the name is '<acc>.<n>' and the comment holds another name, as in '@<acc>.<n> <Illumina name>/<m>'; the name needs restoring.
    - 'lost': the name is '<acc>.<n>' and there is no comment, or the comment is '<n>' or '<n>/<m>'; the archive did not store the Illumina name.
    - 'restored': the name is not '<acc>.<n>' and the comment is; the name was restored before.
    - 'other': the name is not '<acc>.<n>' and the comment is not either, as in a FASTQ file straight from the sequencer; the name needs nothing.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - gzip
    - head

Examples
--------
  1. Check whether a downloaded mate still has the archive's headers.
    '''bash
    print_fastq_names_form "\${dir_out}/SRR32203047_R1.fastq.gz" SRR32203047
    '''

  2. Display the helper's argument contract.
    '''bash
    print_fastq_names_form --help
    '''
EOM
    )

    if [[ "${fil}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${fil}" || -z "${acc}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "expected two arguments: 'fil' and 'acc'."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    # 'head' closes the pipe early, so 'gzip' may exit on SIGPIPE.
    hdr="$(gzip -cd < "${fil}" 2> /dev/null | head -n 1)" || true
    nam="${hdr%% *}"
    nam="${nam#@}"
    cmt=""

    if [[ "${hdr}" == *" "* ]]; then
        cmt="${hdr#* }"
    fi

    # Only an '<acc>.<n>' name needs restoring; any other name is left as is.
    if [[ ! "${nam}" =~ ^"${acc}"\.[0-9]+$ ]]; then
        if [[ "${cmt}" =~ ^"${acc}"\.[0-9]+$ ]]; then
            echo restored
        else
            echo other
        fi
    elif [[ -z "${cmt}" || "${cmt}" =~ ^[0-9]+(/[12])?$ ]]; then
        echo lost
    else
        echo archive
    fi
}


# Restore the Illumina read names that ENA moves into the header comment; the
# file(s) of one entry are rewritten in full or left as downloaded.
function restore_fastq_names() {
    local acc="${1:-}"
    local fq_1="${2:-}"
    local fq_2="${3:-}"
    local fil form tmp_1 tmp_2 res
    local show_help
    local -a arr_fil arr_awk

    show_help=$(cat << EOM
Usage
-----
  restore_fastq_names
    [--help] acc fq_1 fq_2

  Restore the Illumina read names of one downloaded single- or paired-end FASTQ entry from the header comments written by ENA.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  acc : str
    Run accession of the entry.

  2  fq_1 : file
    First or only gzipped FASTQ file.

  3  fq_2 : file
    Second gzipped FASTQ file, or 'NA' for single-end data.

Returns
-------
  Returns 0 when the names are restored, already restored, not stored by the archive, or left as downloaded after a failed check; returns 1 when a file is missing, empty, unreadable, or cannot be replaced.

Notes
-----
  Runtime requirements:
    - awk
    - bash >= 4.4
    - gzip
    - head
    - mv
    - rm

  - ENA writes '@<acc>.<n> <Illumina name>/<m>'; each header becomes '@<Illumina name> <acc>.<n>'. Sequence and quality lines are copied unchanged.
  - Every header of every file must have that form, with spot numbers '<n>' in order from 1, mate suffix '<m>' matching the file (optional for single-end data), and equal names in both mates. Otherwise both files stay as downloaded, and one warning gives the count of failed records and the first failed header.
  - Each file is rewritten to '<file stem>.names.tmp.gz' beside it and moved over the download only when every record passes.
  - A file whose names the archive did not store ('@<acc>.<n> <n>/1'), whose names are already restored, or whose names are not '<acc>.<n>' at all, such as a FASTQ file straight from the sequencer, is left as is with one note and is not read past its first header.

Examples
--------
  1. Restore the names of one downloaded paired-end entry.
    '''bash
    restore_fastq_names \\
        SRR32203047 \\
        "\${dir_out}/SRR32203047_R1.fastq.gz" \\
        "\${dir_out}/SRR32203047_R2.fastq.gz"
    '''

  2. Restore the names of one downloaded single-end entry.
    '''bash
    restore_fastq_names SRR7175377 "\${dir_out}/SRR7175377.fastq.gz" NA
    '''
EOM
    )

    if [[ "${acc}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif (( $# != 3 )); then
        echo_err_func "${FUNCNAME[0]}" \
            "expected three arguments: 'acc', 'fq_1', and 'fq_2'."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    arr_fil=( "${fq_1}" )

    if [[ "${fq_2}" != "NA" ]]; then
        arr_fil+=( "${fq_2}" )
    fi

    for fil in "${arr_fil[@]}"; do
        if [[ ! -s "${fil}" ]]; then
            echo_err "FASTQ file is missing or empty: '${fil}'."
            return 1
        fi
    done

    form="$(print_fastq_names_form "${fq_1}" "${acc}")" || return 1

    case "${form}" in
        lost)
            echo_note \
                "${acc}: read names were not stored by the archive; left" \
                "as downloaded."
            return 0
            ;;

        restored)
            echo_note "${acc}: read names are already restored."
            return 0
            ;;

        other)
            echo_note \
                "${acc}: read names are not in the archive's '<run>.<n>'" \
                "form, so they are left as is."
            return 0
            ;;
    esac

    echo "Restoring the read names of ${acc}." >&2

    tmp_1="${fq_1%.fastq.gz}.names.tmp.gz"
    tmp_2=""

    if [[ "${fq_2}" != "NA" ]]; then
        tmp_2="${fq_2%.fastq.gz}.names.tmp.gz"
    fi

    # Read both mates in step, writing each through its own 'gzip -c'; paths
    # come from the environment, so the shell never reparses them. On success,
    # print 'ok <records>'; otherwise print the counts and first failed header
    # as 'bad <failed> <records> <first header>' and exit 1.
    # shellcheck disable=SC2016  # The program is for awk, not the shell.
    arr_awk=(
        awk
            -v acc="${acc}"
            -v pe="${tmp_2:+1}"
            '
            function check_header(hdr, mate,    sp, pfx, nam) {
                if (substr(hdr, 1, 1) != "@") { return "" }
                sp = index(hdr, " ")
                if (sp == 0) { return "" }
                pfx = substr(hdr, 2, sp - 2)
                nam = substr(hdr, sp + 1)
                if (pfx != acc "." rec) { return "" }

                if (nam ~ /\/[12]$/) {
                    if (substr(nam, length(nam)) != mate) { return "" }
                    nam = substr(nam, 1, length(nam) - 2)
                } else if (pe) {
                    return ""
                }

                if (nam ~ /[ \t]/) { return "" }
                if (gsub(/:/, ":", nam) != 6) { return "" }
                return nam
            }

            function mark_bad(txt) {
                n_bad++
                if (n_bad == 1) { first = txt }
            }

            BEGIN {
                cmd_in = "gzip -cd < \"$FQ_2\""
                cmd_1 = "gzip -c > \"$TMP_1\""
                cmd_2 = "gzip -c > \"$TMP_2\""
            }

            {
                ln_1 = $0

                if (pe && (cmd_in | getline ln_2) <= 0) {
                    mark_bad("second mate ends first")
                    short = 1
                    exit
                }

                if ((NR - 1) % 4 != 0) {
                    if (!n_bad) {
                        print ln_1 | cmd_1
                        if (pe) { print ln_2 | cmd_2 }
                    }
                    next
                }

                rec++
                nam_1 = check_header(ln_1, 1)
                nam_2 = pe ? check_header(ln_2, 2) : nam_1

                if (nam_1 == "" || nam_2 == "" || nam_1 != nam_2) {
                    mark_bad(pe ? ln_1 " | " ln_2 : ln_1)
                    next
                }

                if (!n_bad) {
                    print "@" nam_1 " " acc "." rec | cmd_1
                    if (pe) { print "@" nam_2 " " acc "." rec | cmd_2 }
                }
            }

            END {
                if (pe && !short && (cmd_in | getline ln_2) > 0) {
                    mark_bad("first mate ends first")
                }
                if (NR % 4 != 0) { mark_bad("incomplete last record") }

                close(cmd_1)
                if (pe) { close(cmd_2) }

                if (n_bad) {
                    print "bad", n_bad, rec, first
                    exit 1
                }
                print "ok", rec
            }
            '
    )

    if \
        res="$(
            set -o pipefail
            gzip -cd < "${fq_1}" \
                | FQ_2="${fq_2}" TMP_1="${tmp_1}" TMP_2="${tmp_2}" \
                    "${arr_awk[@]}"
        )"
    then
        mv -f -- "${tmp_1}" "${fq_1}" || {
            rm -f -- "${tmp_1}" ${tmp_2:+"${tmp_2}"}
            echo_err "failed to move '${tmp_1}' to '${fq_1}'."
            return 1
        }

        if [[ -n "${tmp_2}" ]]; then
            mv -f -- "${tmp_2}" "${fq_2}" || {
                rm -f -- "${tmp_2}"
                echo_err "failed to move '${tmp_2}' to '${fq_2}'."
                return 1
            }
        fi

        return 0
    fi

    rm -f -- "${tmp_1}" ${tmp_2:+"${tmp_2}"}

    if [[ "${res}" =~ ^bad\ ([0-9]+)\ ([0-9]+)\ (.*)$ ]]; then
        echo_warn \
            "${acc}: ${BASH_REMATCH[1]} of ${BASH_REMATCH[2]} records do not" \
            "have the archive's read-name form, so the names were not" \
            "restored and the file(s) stay as downloaded; first failed" \
            "header: '${BASH_REMATCH[3]}'."
        return 0
    fi

    echo_err "failed to read or rewrite the FASTQ file(s) of '${acc}'."
    return 1
}


# Link the accession-named FASTQ file(s) of one entry to its custom name.
function link_fastqs() {
    local dir_out="${1:-}"
    local acc="${2:-}"
    local dir_sym="${3:-}"
    local cus="${4:-}"
    local layout="${5:-}"
    local mate src dst
    local show_help
    local -a arr_mate

    show_help=$(cat << EOM
Usage
-----
  link_fastqs
    [--help] dir_out acc dir_sym cus layout

  Create the custom-name symlink(s) for one downloaded single- or paired-end FASTQ entry.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  dir_out : dir
    Directory of the accession-named FASTQ file(s).

  2  acc : str
    Run accession.

  3  dir_sym : dir
    Directory for the custom-named symlink(s).

  4  cus : str
    Custom name.

  5  layout : str
    'pe' for paired-end data ('_R1' and '_R2' files) or 'se' for single-end data.

Returns
-------
  Returns 0 when every link is created; otherwise prints an error and returns 1.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - ln

  - Each accession-named file must exist and be non-empty; an existing link is replaced ('ln -sf').

Examples
--------
  1. Link one downloaded paired-end entry to its custom name.
    '''bash
    link_fastqs "\${dir_out}" SRR32203047 "\${dir_sym}" in_WT_Q_Hho1_6337 pe
    '''

  2. Display the helper's argument contract.
    '''bash
    link_fastqs --help
    '''
EOM
    )

    if [[ "${dir_out}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif (( $# != 5 )); then
        echo_err_func "${FUNCNAME[0]}" \
            "expected five arguments: 'dir_out', 'acc', 'dir_sym', 'cus'," \
            "and 'layout'."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    case "${layout}" in
        pe) arr_mate=( 1 2 ) ;;
        se) arr_mate=( NA )  ;;
        *)
            echo_err_func "${FUNCNAME[0]}" \
                "'layout' must be 'pe' or 'se': '${layout}'."
            return 1
            ;;
    esac

    for mate in "${arr_mate[@]}"; do
        src="$(derive_fastq_path "${dir_out}" "${acc}" "${mate}")" || return 1
        dst="$(derive_fastq_path "${dir_sym}" "${cus}" "${mate}")" || return 1

        if [[ ! -s "${src}" ]]; then
            echo_err "downloaded FASTQ is missing or empty: '${src}'."
            return 1
        fi

        ln -sf "${src}" "${dst}" || {
            echo_err "failed to create symlink '${dst}' -> '${src}'."
            return 1
        }
    done
}


# Download, verify, and link one single- or paired-end FASTQ entry.
function download_fastqs() {
    local srr url_1 url_2 md5_1 md5_2 siz_1 siz_2
    local dir_out dir_sym nam_cus dir_eo nam_job keep_names use_existing
    local dir nam val idx mate layout pfx_log fil form md5 siz reason
    local show_help
    local -a arr_mate arr_url arr_md5 arr_siz arr_fq

    show_help=$(cat << EOM
Usage
-----
  download_fastqs
    [--help]
    --srr <str> --url_1 <str> [--url_2 <str>]
    [--md5_1 <str>] [--md5_2 <str>] [--siz_1 <int>] [--siz_2 <int>]
    --dir_out <dir> --dir_sym <dir> [--use_existing] [--keep_names]
    --nam_cus <str> --dir_eo <dir> --nam_job <str>

  Download one single- or paired-end FASTQ entry, verify each file, restore its Illumina read names, and create its custom-name symlink(s).

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  -sr, --srr : str
    Run accession; downloaded files are named for it.

  -u1, --url_1 : str
    URL of the first FASTQ file.

  -u2, --url_2 : str
    URL of the second FASTQ file for paired-end data (default: 'NA', for single-end data).

  -m1, --md5_1 : str
    Expected MD5 of the first FASTQ file (default: 'NA', not checked).

  -m2, --md5_2 : str
    Expected MD5 of the second FASTQ file (default: 'NA', not checked); for paired-end entries only.

  -s1, --siz_1 : int
    Expected size of the first FASTQ file in bytes (default: 'NA', not checked).

  -s2, --siz_2 : int
    Expected size of the second FASTQ file in bytes (default: 'NA', not checked); for paired-end entries only.

  -do, --dir_out : dir
    Output directory for accession-named FASTQ files.

  -dy, --dir_sym, --dir_symlink : dir
    Directory for custom-named symlinks.

  -ue, --use_existing : flag
    Check and keep a FASTQ file already in '--dir_out' instead of downloading it again.

  -kn, --keep_names : flag
    Keep the read names as the archive wrote them; do not restore the Illumina names.

  -nc, --nam_cus : str
    Custom name for the symlink(s).

  -deo, --dir_eo : dir
    Directory for the 'wget' stdout and stderr log files.

  -nj, --nam_job : str
    Job name, used as the log file prefix.

Returns
-------
  Returns 0 when every file is downloaded, verified, and linked; otherwise prints an error and returns 1.

Notes
-----
  Runtime requirements:
    - awk (unless '--keep_names' is given)
    - bash >= 4.4
    - gzip
    - head
    - ln
    - md5sum (when '--md5_1' or '--md5_2' is given)
    - mv
    - Network access
    - rm
    - wc
    - wget

  - Directories must exist; they are resolved to absolute physical paths, so the symlinks do not depend on the working directory.
  - See 'download_fastq()' for verification and retries, and 'restore_fastq_names()' for the read names.
  - With '--use_existing', a file already in '--dir_out' is checked rather than downloaded again. The check uses the file's MD5 and size when they are given, and 'gzip -t' when they are not.
    + A file whose names are already restored no longer matches the archive's MD5, so it gets 'gzip -t' only.
    + A file that fails its check fails the entry and is not downloaded again.
    + A missing file is downloaded as usual.

Examples
--------
  1. Download, verify, and link one paired-end entry with its ENA MD5 values and sizes.
    '''bash
    download_fastqs \\
        --srr SRR32203047 \\
        --url_1 "\${url_r1}" \\
        --url_2 "\${url_r2}" \\
        --md5_1 6c0bc341f73b0ac0d345fa62cb80a044 \\
        --md5_2 cd9d7b3b45cb059f5f38d36f99c42655 \\
        --siz_1 270090667 \\
        --siz_2 271411747 \\
        --dir_out "\${dir_out}" \\
        --dir_sym "\${dir_sym}" \\
        --nam_cus in_WT_Q_Hho1_6337 \\
        --dir_eo "\${dir_eo}" \\
        --nam_job download_fastqs
    '''

  2. Download and link one single-end entry without an MD5.
    '''bash
    download_fastqs \\
        --srr SRR7175377 \\
        --url_1 "\${url_se}" \\
        --dir_out "\${dir_out}" \\
        --dir_sym "\${dir_sym}" \\
        --nam_cus IP_WT_Q_Brn1_rep2 \\
        --dir_eo "\${dir_eo}" \\
        --nam_job download_fastqs
    '''
EOM
    )

    srr=""
    url_1=""
    url_2="NA"
    md5_1="NA"
    md5_2="NA"
    siz_1="NA"
    siz_2="NA"
    dir_out=""
    dir_sym=""
    nam_cus=""
    dir_eo=""
    nam_job=""
    keep_names=false
    use_existing=false

    while (( $# > 0 )); do
        case "${1}" in
            -h|--hlp|--help)
                echo "${show_help}" >&2
                return 0
                ;;

            -sr|--srr)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                srr="${2}"
                shift 2
                ;;

            -u1|--url[_-]1)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                url_1="${2}"
                shift 2
                ;;

            -u2|--url[_-]2)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                url_2="${2}"
                shift 2
                ;;

            -m1|--md5[_-]1)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                md5_1="${2}"
                shift 2
                ;;

            -m2|--md5[_-]2)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                md5_2="${2}"
                shift 2
                ;;

            -s1|--siz[_-]1)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                siz_1="${2}"
                shift 2
                ;;

            -s2|--siz[_-]2)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                siz_2="${2}"
                shift 2
                ;;

            -do|--dir[_-]out)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                dir_out="${2}"
                shift 2
                ;;

            -dy|--dir[_-]sym|--dir[_-]symlink)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                dir_sym="${2}"
                shift 2
                ;;

            -ue|--use[_-]existing)
                use_existing=true
                shift 1
                ;;

            -kn|--keep[_-]names)
                keep_names=true
                shift 1
                ;;

            -nc|--nam[_-]cus)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                nam_cus="${2}"
                shift 2
                ;;

            -deo|--dir[_-]eo)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                dir_eo="${2}"
                shift 2
                ;;

            -nj|--nam[_-]job)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                nam_job="${2}"
                shift 2
                ;;

            *)
                echo_err_func "${FUNCNAME[0]}" \
                    "unknown option/parameter passed: '${1}'."
                echo >&2
                echo "${show_help}" >&2
                return 1
                ;;
        esac
    done

    validate_var "srr"     "${srr}"     || return 1
    validate_var "url_1"   "${url_1}"   || return 1
    validate_var "url_2"   "${url_2}"   || return 1
    validate_var "md5_1"   "${md5_1}"   || return 1
    validate_var "md5_2"   "${md5_2}"   || return 1
    validate_var "siz_1"   "${siz_1}"   || return 1
    validate_var "siz_2"   "${siz_2}"   || return 1
    validate_var "dir_out" "${dir_out}" || return 1
    validate_var "dir_sym" "${dir_sym}" || return 1
    validate_var "nam_cus" "${nam_cus}" || return 1
    validate_var "dir_eo"  "${dir_eo}"  || return 1
    validate_var "nam_job" "${nam_job}" || return 1

    # The second mate's checks apply only to paired-end entries; 'NA' counts as
    # not given.
    if [[ "${url_2}" == "NA" ]]; then
        for nam in md5_2 siz_2; do
            if [[ "${!nam}" != "NA" ]]; then
                echo_err \
                    "'--${nam}' is for paired-end entries ('--url_2'); it" \
                    "has no effect without '--url_2'."
                return 1
            fi
        done
    fi

    for dir in dir_out dir_sym dir_eo; do
        val="${!dir}"

        if [[ ! -d "${val}" ]]; then
            echo_err "directory does not exist: '${val}'"
            return 1
        fi

        val="$(cd "${val}" > /dev/null 2>&1 && pwd -P)" || {
            echo_err "failed to resolve directory: '${!dir}'"
            return 1
        }
        printf -v "${dir}" '%s' "${val}"
    done

    for nam in md5_1 md5_2; do
        val="${!nam}"

        if [[ "${val}" != "NA" && ! "${val}" =~ ^[0-9A-Fa-f]{32}$ ]]; then
            echo_err \
                "'--${nam}' must be 32 hexadecimal characters or 'NA':" \
                "'${val}'."
            return 1
        fi
    done

    for nam in siz_1 siz_2; do
        val="${!nam}"

        if [[ "${val}" != "NA" && ! "${val}" =~ ^[1-9][0-9]*$ ]]; then
            echo_err \
                "'--${nam}' must be a positive integer number of bytes or" \
                "'NA': '${val}'."
            return 1
        fi
    done

    check_pgrm_path gzip || return 1
    check_pgrm_path head || return 1
    check_pgrm_path ln   || return 1
    check_pgrm_path wget || return 1

    if [[ "${keep_names}" == "false" ]]; then
        check_pgrm_path awk || return 1
    fi

    if [[ "${md5_1}" != "NA" || "${md5_2}" != "NA" ]]; then
        check_pgrm_path md5sum || return 1
    fi

    if [[ "${url_2}" != "NA" ]]; then
        layout=pe
        arr_mate=( 1 2 )
    else
        layout=se
        arr_mate=( NA )
    fi

    arr_url=( "${url_1}" "${url_2}" )
    arr_md5=( "${md5_1}" "${md5_2}" )
    arr_siz=( "${siz_1}" "${siz_2}" )

    arr_fq=( NA NA )

    for idx in "${!arr_mate[@]}"; do
        mate="${arr_mate[idx]}"
        fil="$(derive_fastq_path "${dir_out}" "${srr}" "${mate}")" || return 1
        arr_fq[idx]="${fil}"

        if [[ "${mate}" == "NA" ]]; then
            pfx_log="${dir_eo}/${nam_job}.${srr}"
        else
            pfx_log="${dir_eo}/${nam_job}.${srr}_R${mate}"
        fi

        # Check a file already in place instead of downloading it. Restored
        # names no longer match the archive's MD5; do only 'gzip -t' then.
        if [[ "${use_existing}" == "true" && -e "${fil}" ]]; then
            echo "Checking existing ${fil} for ${srr}." >&2
            form="$(print_fastq_names_form "${fil}" "${srr}")" || return 1
            md5="${arr_md5[idx]}"
            siz="${arr_siz[idx]}"

            if [[ "${form}" == "restored" ]]; then
                md5=NA
                siz=NA
            fi

            if ! \
                reason="$(verify_fastq "${fil}" "${md5}" "${siz}")"
            then
                echo_err \
                    "existing file '${fil}' fails its check: ${reason};" \
                    "remove it or run without '--use_existing'."
                return 1
            fi

            continue
        fi

        download_fastq \
            "${srr}" \
            "${arr_url[idx]}" \
            "${fil}" \
            "${arr_md5[idx]}" \
            "${arr_siz[idx]}" \
            "${pfx_log}.stdout.txt" \
            "${pfx_log}.stderr.txt" \
            || return 1
    done

    if [[ "${keep_names}" == "false" ]]; then
        restore_fastq_names "${srr}" "${arr_fq[0]}" "${arr_fq[1]}" || return 1
    fi

    echo "Symlinking ${srr} to ${nam_cus}." >&2
    link_fastqs "${dir_out}" "${srr}" "${dir_sym}" "${nam_cus}" "${layout}"
}
