#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: align_fastqs_separate.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Anthropic Claude Code (Opus 5.5) was used in design, development, and
# documentation, with all output reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


# _route_unal_bam
# _route_filt_bam
# _sort_coord_diag_bam
# _mark_dup_all_bam
# _split_diag_bam
# _count_diag_cat
# _write_diag_tsv


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

# The AWK program that '_write_diag_tsv' runs on the improper file.
_fil_awk_improp="$(
    cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
)/../../awk/classify_improp_pairs.awk"


# Step 2 with '--unal_mode separate': sort Step 1's unaligned work file ("Set
# U") by queryname and list its read names. For paired-end data, also sort the
# aligned reads ("Set M") by queryname, move every template with records in
# both work files wholly into Set U, and fix mates in Set U.
function _route_unal_bam() {
    local func="${1:-}"
    local dry_run="${2:-}"
    local threads="${3:-}"
    local fq_2="${4:-}"
    local fil_wrk="${5:-}"
    local bam_qnam="${6:-}"
    local bam_set_u="${7:-}"
    local txt_nam="${fil_wrk%.bam}.unal.names.txt"
    local bam_srt bam_mov bam_mate n_rec n_tpl
    local -a arr_cmd arr_view arr_cut arr_uniq arr_srt arr_mrg arr_fix
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _route_unal_bam
    [--help] func dry_run threads fq_2 fil_wrk bam_qnam bam_set_u

  Sort Step 1's unaligned work file ("Set U") by queryname and list its read names (Step 2 of 'align_fastqs()' with '--unal_mode separate'). The list goes to the path of 'fil_wrk' with '.unal.names.txt' for '.bam', and Step 4 uses it to split Set U back out.

  For paired-end data, also sort the BAM work file ("Set M") by queryname into 'bam_qnam'. As it goes, move every record whose read name is in Set U into Set U, so that no template has records in both. Print a note with the number of templates and records moved. Then fix mates in the newly filtered Set U.

  For single-end data, only Set U is sorted and listed. Each read is its own template, so no records move and no mates are fixed, and the caller sorts Set M itself.

  A template is moved by its read name, not by a count of its records, so the rule holds for any number of records per read.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  func : str
    Name of the calling function for diagnostics.

  2  dry_run : bool
    Whether to print the commands to stdout instead of running them.

  3  threads : int
    Number of threads to use.

  4  fq_2 : file
    Second FASTQ input file, or empty for single-end data.

  5  fil_wrk : file
    BAM work file from Step 1 ("Set M"); names the read-name list.

  6  bam_qnam : file
    Queryname-sorted BAM file to write for paired-end data: Set M without the moved records.

  7  bam_set_u : file
    Set U from Step 1, rewritten sorted by queryname, with the moved records and fixed mates.

Returns
-------
  Returns 0 on success or after printing the commands for a dry run; prints an error and returns 1 otherwise.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - cut
    - mv
    - samtools
    - uniq
    - wc

Examples
--------
  1. Move the templates split between the two work files of a paired-end run into Set U.
    '''bash
    _route_unal_bam \\
        align_fastqs \\
        false \\
        2 \\
        tests/fixtures/align_fastqs/fastq/pe/bridge_pe_R2.atria.fastq.gz \\
        work/bridge.bam \\
        work/bridge.qnam.bam \\
        work/bridge.setU.bam
    '''

  2. Display the helper's argument contract.
    '''bash
    _route_unal_bam --help
    '''
EOM
    )

    if [[ "${func}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    # Sort Set U by queryname, as Set M is. Like every work file that a later
    # step rewrites, it is written at compression level 1.
    bam_srt="${bam_set_u%.bam}.qnam.bam"
    arr_cmd=(
        samtools sort
            -@ "${threads}"
            -n
            -O "bam,level=1"
            -o "${bam_srt}"
            "${bam_set_u}"
    )

    if ! \
        _run_or_print_cmd "${dry_run}" arr_cmd
    then
        echo_err_func "${func}" \
            "Step #2: failed to queryname-sort the unaligned-read work file."
        return 1
    fi

    arr_cmd=( mv -f "${bam_srt}" "${bam_set_u}" )

    if ! \
        _run_or_print_cmd "${dry_run}" arr_cmd
    then
        echo_err_func "${func}" \
            "Step #2: failed to rename the queryname-sorted unaligned-read" \
            "work file."
        return 1
    fi

    # List Set U's read names once each; its records are grouped by name.
    arr_view=( samtools view -@ "${threads}" "${bam_set_u}" )
    arr_cut=( cut -f 1 )
    arr_uniq=( uniq )

    if [[ "${dry_run}" == "true" ]]; then
        _print_dry_run_cmd --out "${txt_nam}" arr_view arr_cut arr_uniq
    elif ! (
        set -o pipefail
        "${arr_view[@]}" | "${arr_cut[@]}" | "${arr_uniq[@]}" > "${txt_nam}"
    ); then
        rm -f "${txt_nam}"
        echo_err_func "${func}" \
            "Step #2: failed to list the read names in the unaligned-read" \
            "work file."
        return 1
    fi

    # A single-end read is one template, so none is split between the files.
    if [[ -z "${fq_2}" ]]; then return 0; fi

    # Sort Set M by queryname and, in the same stream, move every record whose
    # read name is in Set U. For example, a BWA-backtrack partner of an
    # unaligned mate (#75) joins its mate.
    bam_mov="${fil_wrk%.bam}.moved.bam"
    arr_srt=(
        samtools sort
            -@ "${threads}"
            -n
            -u
            -o -
            "${fil_wrk}"
    )
    arr_cmd=(
        samtools view
            -@ "${threads}"
            -O "bam,level=1"
            -N "${txt_nam}"
            -o "${bam_mov}"
            -U "${bam_qnam}"
            -
    )

    if [[ "${dry_run}" == "true" ]]; then
        _print_dry_run_cmd arr_srt arr_cmd
    elif ! (
        set -o pipefail
        "${arr_srt[@]}" | "${arr_cmd[@]}"
    ); then
        echo_err_func "${func}" \
            "Step #2: failed to queryname-sort the BAM work file and move" \
            "templates split between the work files."
        return 1
    fi

    # Report how many templates and records were moved, so that each run shows
    # how many templates Step 1 split between the work files, such as
    # BWA-backtrack's partners of unaligned mates (#75). The moved file is
    # sorted by queryname, so 'uniq' counts its templates.
    if [[ "${dry_run}" != "true" ]]; then
        if ! \
            n_rec="$(samtools view -c -@ "${threads}" "${bam_mov}")"
        then
            echo_err_func "${func}" \
                "Step #2: failed to count the moved records."
            return 1
        fi

        if ! \
            n_tpl="$(
                set -o pipefail
                samtools view -@ "${threads}" "${bam_mov}" \
                    | cut -f 1 \
                    | uniq \
                    | wc -l
            )"
        then
            echo_err_func "${func}" \
                "Step #2: failed to count the moved templates."
            return 1
        fi

        echo_note_func "${func}" \
            "Step #2: moved ${n_tpl// /} template(s) split between the work" \
            "files into Set U (${n_rec} record(s))."
    fi

    # Merge the moved templates into Set U and fix mates in the same stream.
    # Every template in Set U is now whole, so fixing mates here gives the same
    # result as fixing them beside Set M. 'merge' gives clashing '@PG' IDs
    # random suffixes, and '-s 0' seeds them so the header is reproducible.
    bam_mate="${bam_set_u%.bam}.mate.bam"
    arr_mrg=(
        samtools merge
            -@ "${threads}"
            -n
            -s 0
            -u
            -o -
            "${bam_set_u}"
            "${bam_mov}"
    )
    arr_fix=(
        samtools fixmate
            -@ "${threads}"
            -c
            -m
            -O "bam,level=1"
            -
            "${bam_mate}"
    )

    if [[ "${dry_run}" == "true" ]]; then
        _print_dry_run_cmd arr_mrg arr_fix
    elif ! (
        set -o pipefail
        "${arr_mrg[@]}" | "${arr_fix[@]}"
    ); then
        echo_err_func "${func}" \
            "Step #2: failed to merge the moved templates into the" \
            "unaligned-read work file and fix its mate pairs."
        return 1
    fi

    arr_cmd=( mv -f "${bam_mate}" "${bam_set_u}" )

    if ! \
        _run_or_print_cmd "${dry_run}" arr_cmd
    then
        echo_err_func "${func}" \
            "Step #2: failed to rename the mate-fixed unaligned-read work" \
            "file."
        return 1
    fi

    arr_cmd=( rm -f "${bam_mov}" )

    if ! \
        _run_or_print_cmd "${dry_run}" arr_cmd
    then
        echo_err_func "${func}" \
            "Step #2: failed to delete the moved-template work file."
        return 1
    fi
}


# Step 2 with '--unal_mode separate': sort the templates that the template
# filter dropped ("Set R") into improper pairs ("Set I") and low-MAPQ templates
# ("Set Q"), then list the read names of every diagnostic work file.
function _route_filt_bam() {
    local func="${1:-}"
    local dry_run="${2:-}"
    local threads="${3:-}"
    local req_flg="${4:-}"
    local fq_2="${5:-}"
    local fil_wrk="${6:-}"
    local stm="${fil_wrk%.bam}"
    local bam_rej="${stm}.setR.bam"
    local bam_rej0="${stm}.setI0.bam"
    local bam_imp="${stm}.setI1.bam"
    local bam_set_i="${stm}.setI.bam"
    local bam_set_q="${stm}.setQ.bam"
    local txt_unal="${stm}.unal.names.txt"
    local txt_imp="${stm}.setI.names.txt"
    local txt_mapq="${stm}.setQ.names.txt"
    local txt_diag="${stm}.diag.names.txt"
    local bam_mate
    local -a arr_cmd arr_view arr_cut arr_uniq arr_awk
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _route_filt_bam
    [--help] func dry_run threads req_flg fq_2 fil_wrk

  Sort the templates that Step 2's template filter dropped into two queryname-sorted work files beside 'fil_wrk': the improper pairs and the low-MAPQ templates. Then list the read names of each diagnostic work file (Step 2 of 'align_fastqs()' with '--unal_mode separate').

  A dropped template whose two primary records still carry the proper-pair flag failed only '--mapq', so it goes to the low-MAPQ work file (the path of 'fil_wrk' with '.setQ.bam' for '.bam'). Every other dropped template goes to the improper work file ('.setI.bam'), as do the records that '--req_flg' kept out of mate fixing, so a template that fails both filters is improper.

  Without 'req_flg' or for single-end data, every dropped template goes to the low-MAPQ work file.

  Single-end data gets no improper work file, and paired-end data without 'req_flg' gets an empty one.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  func : str
    Name of the calling function for diagnostics.

  2  dry_run : bool
    Whether to print the commands to stdout instead of running them.

  3  threads : int
    Number of threads to use.

  4  req_flg : bool
    Whether '--req_flg' applied to the template filter, which it does only for paired-end data.

  5  fq_2 : file
    Second FASTQ input file, or empty for single-end data.

  6  fil_wrk : file
    BAM work file of Step 2; names every work file and name list.

Returns
-------
  Returns 0 on success or after printing the commands for a dry run; prints an error and returns 1 otherwise.

Notes
-----
  Runtime requirements:
    - awk
    - bash >= 4.4
    - cat
    - cut
    - mv
    - rm
    - samtools
    - uniq

  - Reads the dropped templates from the path of 'fil_wrk' with '.setR.bam' for '.bam' and, with 'req_flg', the records kept out of mate fixing from '.setI0.bam'; deletes both.
  - Writes the read names of the improper and low-MAPQ work files, once each, with '.setI.names.txt' and '.setQ.names.txt' for '.bam'.
  - It also writes the names of every diagnostic work file, including the routing's unaligned-read list ('.unal.names.txt'), with '.diag.names.txt'.

Examples
--------
  1. Sort the dropped templates of a paired-end run with '--req_flg'.
    '''bash
    _route_filt_bam \\
        align_fastqs \\
        false \\
        2 \\
        true \\
        tests/fixtures/align_fastqs/fastq/pe/pairs_pe_R2.atria.fastq.gz \\
        work/pairs.bam
    '''

  2. Display the helper's argument contract.
    '''bash
    _route_filt_bam --help
    '''
EOM
    )

    if [[ "${func}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    arr_cut=( cut -f 1 )
    arr_uniq=( uniq )

    if [[ -z "${fq_2}" ]]; then
        # A single-end read has no mate to be improper with, so every dropped
        # read failed '--mapq'.
        arr_cmd=( mv -f "${bam_rej}" "${bam_set_q}" )

        if ! \
            _run_or_print_cmd "${dry_run}" arr_cmd
        then
            echo_err_func "${func}" \
                "Step #2: failed to rename the low-MAPQ work file."
            return 1
        fi

        arr_view=( samtools view -@ "${threads}" "${bam_set_q}" )

        if [[ "${dry_run}" == "true" ]]; then
            _print_dry_run_cmd --out "${txt_mapq}" arr_view arr_cut arr_uniq
        elif ! (
            set -o pipefail
            "${arr_view[@]}" | "${arr_cut[@]}" | "${arr_uniq[@]}" \
                > "${txt_mapq}"
        ); then
            echo_err_func "${func}" \
                "Step #2: failed to list the read names in the low-MAPQ work" \
                "file."
            return 1
        fi
    else
        # Fix mates in the records that '--req_flg' kept out of mate fixing, as
        # Markdup 2 needs the tags that 'fixmate -m' writes.
        if [[ "${req_flg}" == "true" ]]; then
            bam_mate="${bam_rej0%.bam}.mate.bam"
            arr_cmd=(
                samtools fixmate
                    -@ "${threads}"
                    -c
                    -m
                    -O "bam,level=1"
                    "${bam_rej0}"
                    "${bam_mate}"
            )

            if ! \
                _run_or_print_cmd "${dry_run}" arr_cmd
            then
                echo_err_func "${func}" \
                    "Step #2: failed to fix mate pairs in the records that" \
                    "'--req_flg' kept out."
                return 1
            fi

            arr_cmd=( mv -f "${bam_mate}" "${bam_rej0}" )

            if ! \
                _run_or_print_cmd "${dry_run}" arr_cmd
            then
                echo_err_func "${func}" \
                    "Step #2: failed to rename the mate-fixed records that" \
                    "'--req_flg' kept out."
                return 1
            fi
        fi

        # List the dropped templates that failed only '--mapq'. With '--
        # req_flg', these are the ones whose two primary records still carry
        # 0x2; without it, they're all of them. The file is grouped by name, so
        # 'uniq' counts each template's records.
        if [[ "${req_flg}" == "true" ]]; then
            arr_view=(
                samtools view -@ "${threads}" -F 2308 -f 2 "${bam_rej}"
            )
            arr_uniq=( uniq -c )
            # shellcheck disable=SC2016  # The program is AWK.
            arr_awk=( awk '$1 == 2 { print $2 }' )

            if [[ "${dry_run}" == "true" ]]; then
                _print_dry_run_cmd \
                    --out "${txt_mapq}" arr_view arr_cut arr_uniq arr_awk
            elif ! (
                set -o pipefail
                "${arr_view[@]}" \
                    | "${arr_cut[@]}" \
                    | "${arr_uniq[@]}" \
                    | "${arr_awk[@]}" \
                    > "${txt_mapq}"
            ); then
                echo_err_func "${func}" \
                    "Step #2: failed to list the dropped templates that" \
                    "failed only '--mapq'."
                return 1
            fi
        else
            arr_view=( samtools view -@ "${threads}" "${bam_rej}" )

            if [[ "${dry_run}" == "true" ]]; then
                _print_dry_run_cmd \
                    --out "${txt_mapq}" arr_view arr_cut arr_uniq
            elif ! (
                set -o pipefail
                "${arr_view[@]}" | "${arr_cut[@]}" | "${arr_uniq[@]}" \
                    > "${txt_mapq}"
            ); then
                echo_err_func "${func}" \
                    "Step #2: failed to list the dropped templates."
                return 1
            fi
        fi

        arr_cmd=(
            samtools view
                -@ "${threads}"
                -O "bam,level=1"
                -N "${txt_mapq}"
                -o "${bam_set_q}"
                -U "${bam_imp}"
                "${bam_rej}"
        )

        if ! \
            _run_or_print_cmd "${dry_run}" arr_cmd
        then
            echo_err_func "${func}" \
                "Step #2: failed to split the dropped templates into the" \
                "improper and low-MAPQ work files."
            return 1
        fi

        # Join the records that '--req_flg' kept out of mate fixing to the
        # improper pairs. '-s 0' seeds the random suffixes that 'merge' gives
        # clashing '@PG' IDs, so the header is reproducible.
        if [[ "${req_flg}" == "true" ]]; then
            arr_cmd=(
                samtools merge
                    -@ "${threads}"
                    -n
                    -s 0
                    -O "bam,level=1"
                    -o "${bam_set_i}"
                    "${bam_rej0}"
                    "${bam_imp}"
            )
        else
            arr_cmd=( mv -f "${bam_imp}" "${bam_set_i}" )
        fi

        if ! \
            _run_or_print_cmd "${dry_run}" arr_cmd
        then
            echo_err_func "${func}" \
                "Step #2: failed to write the improper work file."
            return 1
        fi

        arr_cmd=( rm -f "${bam_rej}" "${bam_rej0}" "${bam_imp}" )

        if ! \
            _run_or_print_cmd "${dry_run}" arr_cmd
        then
            echo_err_func "${func}" \
                "Step #2: failed to delete the dropped-template work files."
            return 1
        fi

        # List the improper work file's read names once each; its records are
        # grouped by name.
        arr_view=( samtools view -@ "${threads}" "${bam_set_i}" )
        arr_uniq=( uniq )

        if [[ "${dry_run}" == "true" ]]; then
            _print_dry_run_cmd --out "${txt_imp}" arr_view arr_cut arr_uniq
        elif ! (
            set -o pipefail
            "${arr_view[@]}" | "${arr_cut[@]}" | "${arr_uniq[@]}" \
                > "${txt_imp}"
        ); then
            echo_err_func "${func}" \
                "Step #2: failed to list the read names in the improper work" \
                "file."
            return 1
        fi
    fi

    # Step 4 keeps only these names from Markdup 2, which also reads the main
    # file.
    arr_cmd=( cat "${txt_unal}" "${txt_mapq}" )
    if [[ -n "${fq_2}" ]]; then arr_cmd+=( "${txt_imp}" ); fi

    if [[ "${dry_run}" == "true" ]]; then
        _print_dry_run_cmd --out "${txt_diag}" arr_cmd
    elif ! \
        "${arr_cmd[@]}" > "${txt_diag}"
    then
        echo_err_func "${func}" \
            "Step #2: failed to list the read names of the diagnostic work" \
            "files."
        return 1
    fi
}


# Step 3 with '--unal_mode separate': sort the diagnostic work files of Step 2
# together by coordinates, for Markdup 2.
function _sort_coord_diag_bam() {
    local func="${1:-}"
    local dry_run="${2:-}"
    local threads="${3:-}"
    local fq_2="${4:-}"
    local fil_wrk="${5:-}"
    local stm="${fil_wrk%.bam}"
    local bam_diag="${stm}.diag.bam"
    local -a arr_in arr_mrg arr_srt arr_cmd
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _sort_coord_diag_bam
    [--help] func dry_run threads fq_2 fil_wrk

  Merge the queryname-sorted diagnostic work files of Step 2 beside 'fil_wrk', sort them by coordinates into one work file, and delete them (Step 3 of 'align_fastqs()' with '--unal_mode separate').

  The inputs are the unaligned-read, improper, and low-MAPQ work files ('.setU.bam', '.setI.bam', and '.setQ.bam' for '.bam').

  The output is the path of 'fil_wrk' with '.diag.bam' for '.bam'.

  Single-end data has no improper work file, so only the other two are merged.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  func : str
    Name of the calling function for diagnostics.

  2  dry_run : bool
    Whether to print the commands to stdout instead of running them.

  3  threads : int
    Number of threads to use.

  4  fq_2 : file
    Second FASTQ input file, or empty for single-end data.

  5  fil_wrk : file
    BAM work file of the main output; names the diagnostic work files.

Returns
-------
  Returns 0 on success or after printing the commands for a dry run; prints an error and returns 1 otherwise.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - rm
    - samtools

  - The work file is written at compression level 1, since Step 4 reads it whole and splits it again.
  - 'samtools merge' renames clashing '@PG' IDs with random suffixes; '-s 0' seeds them, making the header reproducible.

Examples
--------
  1. Sort the diagnostic work files of a paired-end run.
    '''bash
    _sort_coord_diag_bam \\
        align_fastqs \\
        false \\
        2 \\
        tests/fixtures/align_fastqs/fastq/pe/pairs_pe_R2.atria.fastq.gz \\
        work/pairs.bam
    '''

  2. Display the helper's argument contract.
    '''bash
    _sort_coord_diag_bam --help
    '''
EOM
    )

    if [[ "${func}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    arr_in=( "${stm}.setU.bam" )
    if [[ -n "${fq_2}" ]]; then arr_in+=( "${stm}.setI.bam" ); fi
    arr_in+=( "${stm}.setQ.bam" )

    arr_mrg=(
        samtools merge
            -@ "${threads}"
            -n
            -s 0
            -u
            -o -
            "${arr_in[@]}"
    )
    arr_srt=(
        samtools sort
            -@ "${threads}"
            -O "bam,level=1"
            -o "${bam_diag}"
            -
    )

    if [[ "${dry_run}" == "true" ]]; then
        _print_dry_run_cmd arr_mrg arr_srt
    elif ! (
        set -o pipefail
        "${arr_mrg[@]}" | "${arr_srt[@]}"
    ); then
        echo_err_func "${func}" \
            "Step #3: failed to merge and coordinate-sort the diagnostic" \
            "work files."
        return 1
    fi

    arr_cmd=( rm -f "${arr_in[@]}" )

    if ! \
        _run_or_print_cmd "${dry_run}" arr_cmd
    then
        echo_err_func "${func}" \
            "Step #3: failed to delete the queryname-sorted diagnostic work" \
            "files."
        return 1
    fi
}


# Step 4 with '--unal_mode separate': mark duplicates in all the files together
# (Markdup 2), keeping the marked records of the diagnostic work file.
function _mark_dup_all_bam() {
    local func="${1:-}"
    local dry_run="${2:-}"
    local threads="${3:-}"
    local fil_wrk="${4:-}"
    local txt_all="${5:-}"
    local dist_nam="${6:-}"
    local stm="${fil_wrk%.bam}"
    local bam_diag="${stm}.diag.bam"
    local bam_mrk="${stm}.diag.mark.bam"
    local txt_diag="${stm}.diag.names.txt"
    local dist_opt rc_opt
    local -a arr_mrg arr_mrk arr_view arr_cmd arr_opt
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _mark_dup_all_bam
    [--help] func dry_run threads fil_wrk txt_all dist_nam

  Mark duplicates in the main BAM work file and the diagnostic work file beside it ('.diag.bam' for '.bam') together. Keep only the diagnostic records, now marked, in the diagnostic work file (Step 4 of 'align_fastqs()' with '--unal_mode separate': Markdup 2).

  The main work file is read, not changed. Step 4 has already marked it alone (Markdup 1), as under '--unal_mode remove'. Those flags need no clearing: 'samtools markdup' only sets the duplicate flag, the main file's records are dropped from this pass's output, and its statistics count only the duplicates it marks itself. The statistics describe all the files together.

  Optical duplicates are told apart ('-d 2500') only when '_check_optical_names()' allows it for both files' read names.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  func : str
    Name of the calling function for diagnostics.

  2  dry_run : bool
    Whether to print the commands to stdout instead of running them.

  3  threads : int
    Number of threads to use.

  4  fil_wrk : file
    Coordinate-sorted BAM work file of the main output, marked alone by Step 4.

  5  txt_all : file
    Path for the duplicate statistics, without '.gz'; the file is written gzip-compressed, as 'txt_all.gz'.

  6  dist_nam : str
    Name of the variable to set to the optical duplicate distance, '2500', or to 'not used' and the reason.

Returns
-------
  Returns 0 on success or after printing the commands for a dry run; prints an error and returns 1 otherwise.

Notes
-----
  Runtime requirements:
    - awk
    - bash >= 4.4
    - cut
    - gzip
    - mv
    - samtools

  - The diagnostic records are kept by the read names Step 2 listed beside 'fil_wrk', with '.diag.names.txt' for '.bam'.
  - 'samtools merge' renames clashing '@PG' IDs with random suffixes; '-s 0' seeds them, making the header reproducible.
  - A dry run cannot read the names, so it prints the '-d 2500' form, as '_mark_dup_bam()' does.

Examples
--------
  1. Mark duplicates across a paired-end run's files.
    '''bash
    _mark_dup_all_bam \\
        align_fastqs \\
        false \\
        2 \\
        work/pairs.bam \\
        work/pairs.all.markdup.txt \\
        dist_all
    '''

  2. Display the helper's argument contract.
    '''bash
    _mark_dup_all_bam --help
    '''
EOM
    )

    if [[ "${func}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    if [[ "${dry_run}" == "true" ]]; then
        printf '%s\n' \
            "# Mark duplicates in all files together (Markdup 2)."
    fi

    local -n dist_all_out="${dist_nam}"

    # Tell optical duplicates apart only when the read names allow it.
    dist_opt="2500"
    if [[ "${dry_run}" == "false" ]]; then
        rc_opt=0
        _check_optical_names \
            "${func}" \
            "${threads}" \
            "${fil_wrk}" \
            dist_opt \
            "${bam_diag}" \
            || rc_opt=$?

        if (( rc_opt == 2 )); then
            return 1
        fi
    fi

    # shellcheck disable=SC2034  # Read by name in the caller.
    dist_all_out="${dist_opt}"
    arr_opt=()
    if [[ "${dist_opt}" == "2500" ]]; then arr_opt=( -d 2500 ); fi

    # Merge, mark, and keep the diagnostic records in one stream, without
    # writing the merged file.
    arr_mrg=(
        samtools merge
            -@ "${threads}"
            -s 0
            -u
            -o -
            "${fil_wrk}"
            "${bam_diag}"
    )
    arr_mrk=(
        samtools markdup
            -@ "${threads}"
            -t
            "${arr_opt[@]}"
            -u
            -f "${txt_all}"
            -
            -
    )
    arr_view=(
        samtools view
            -@ "${threads}"
            -O "bam,level=1"
            -N "${txt_diag}"
            -o "${bam_mrk}"
            -
    )

    if [[ "${dry_run}" == "true" ]]; then
        _print_dry_run_cmd arr_mrg arr_mrk arr_view
    elif ! (
        set -o pipefail
        "${arr_mrg[@]}" | "${arr_mrk[@]}" | "${arr_view[@]}"
    ); then
        echo_err_func "${func}" \
            "Step #4: failed to mark duplicates in all files together."
        return 1
    fi

    # As the statistics' last line, record whether optical duplicates were told
    # apart and, if not, why.
    if [[ "${dry_run}" == "true" ]]; then
        printf '%s %q\n\n' \
            "printf 'OPTICAL DISTANCE: %s\\n' 2500 >>" \
            "${txt_all}"
    else
        printf 'OPTICAL DISTANCE: %s\n' "${dist_opt}" >> "${txt_all}"
    fi

    arr_cmd=( gzip -f "${txt_all}" )

    if ! \
        _run_or_print_cmd "${dry_run}" arr_cmd
    then
        echo_err_func "${func}" \
            "Step #4: failed to compress the duplicate statistics of all" \
            "files."
        return 1
    fi

    arr_cmd=( mv -f "${bam_mrk}" "${bam_diag}" )

    if ! \
        _run_or_print_cmd "${dry_run}" arr_cmd
    then
        echo_err_func "${func}" \
            "Step #4: failed to rename the duplicate-marked diagnostic work" \
            "file."
        return 1
    fi
}


# Step 4 with '--unal_mode separate': split the duplicate-marked diagnostic
# work file into the unaligned-read, improper, and low-MAPQ files.
function _split_diag_bam() {
    local func="${1:-}"
    local dry_run="${2:-}"
    local threads="${3:-}"
    local fq_2="${4:-}"
    local fil_wrk="${5:-}"
    local fil_unal="${6:-}"
    local fil_improp="${7:-}"
    local fil_mapq="${8:-}"
    local fin="${9:-true}"
    local stm="${fil_wrk%.bam}"
    local bam_diag="${stm}.diag.bam"
    local bam_rest="${stm}.rest.bam"
    local txt_unal="${stm}.unal.names.txt"
    local txt_imp="${stm}.setI.names.txt"
    local txt_mapq="${stm}.setQ.names.txt"
    local txt_diag="${stm}.diag.names.txt"
    local -a arr_cmd arr_lvl
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _split_diag_bam
    [--help] func dry_run threads fq_2 fil_wrk fil_unal fil_improp fil_mapq [fin]

  Split the duplicate-marked diagnostic work file beside 'fil_wrk' ('.diag.bam' for '.bam') by the read names Step 2 listed (Step 4 of 'align_fastqs()' with '--unal_mode separate').

  The unaligned reads and their mates go to 'fil_unal', the improper pairs to 'fil_improp', and the rest (the low-MAPQ templates) to 'fil_mapq'.

  All stay coordinate-sorted, and Step 5 indexes them.

  Single-end data has no improper file, so 'fil_improp' is not written.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  func : str
    Name of the calling function for diagnostics.

  2  dry_run : bool
    Whether to print the commands to stdout instead of running them.

  3  threads : int
    Number of threads to use.

  4  fq_2 : file
    Second FASTQ input file, or empty for single-end data.

  5  fil_wrk : file
    BAM work file of the main output; names the diagnostic work file and the name lists.

  6  fil_unal : file
    BAM file to write the unaligned reads and their mates to.

  7  fil_improp : file
    BAM file to write the improper pairs to.

  8  fil_mapq : file
    BAM file to write the low-MAPQ templates to.

  9  fin : bool
    Whether the files ship as they are, and so are written at the default compression level. 'false' writes them at level 1 for CRAM output, which Step 6 converts them to (default: '${fin}').

Returns
-------
  Returns 0 on success or after printing the commands for a dry run; prints an error and returns 1 otherwise.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - mv
    - rm
    - samtools

  - Step 2 moves templates whole and lists each work file's names, so no read name is in two lists and the split is exact.
  - Deletes the diagnostic work file and the name lists.

Examples
--------
  1. Split a paired-end run's diagnostic work file.
    '''bash
    _split_diag_bam \\
        align_fastqs \\
        false \\
        2 \\
        tests/fixtures/align_fastqs/fastq/pe/pairs_pe_R2.atria.fastq.gz \\
        work/pairs.bam \\
        work/pairs.unal.bam \\
        work/pairs.improp.bam \\
        work/pairs.mapq.bam
    '''

  2. Display the helper's argument contract.
    '''bash
    _split_diag_bam --help
    '''
EOM
    )

    if [[ "${func}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    if [[ "${dry_run}" == "true" ]]; then
        printf '%s\n' \
            "# Split the diagnostic work file into the diagnostic files."
    fi

    arr_lvl=( -b )
    if [[ "${fin}" == "false" ]]; then arr_lvl=( -O "bam,level=1" ); fi

    # For single-end data, the reads not in the side file are the low-MAPQ
    # reads; for paired-end data, the improper pairs are split from them next.
    arr_cmd=(
        samtools view
            -@ "${threads}"
            "${arr_lvl[@]}"
            -N "${txt_unal}"
            -o "${fil_unal}"
    )
    if [[ -n "${fq_2}" ]]; then
        arr_cmd+=( -U "${bam_rest}" "${bam_diag}" )
    else
        arr_cmd+=( -U "${fil_mapq}" "${bam_diag}" )
    fi

    if ! \
        _run_or_print_cmd "${dry_run}" arr_cmd
    then
        echo_err_func "${func}" \
            "Step #4: failed to split the unaligned reads and their mates" \
            "into the side file."
        return 1
    fi

    if [[ -n "${fq_2}" ]]; then
        arr_cmd=(
            samtools view
                -@ "${threads}"
                "${arr_lvl[@]}"
                -N "${txt_imp}"
                -o "${fil_improp}"
                -U "${fil_mapq}"
                "${bam_rest}"
        )

        if ! \
            _run_or_print_cmd "${dry_run}" arr_cmd
        then
            echo_err_func "${func}" \
                "Step #4: failed to split the improper pairs from the" \
                "low-MAPQ templates."
            return 1
        fi
    fi

    arr_cmd=( rm -f "${bam_diag}" )
    if [[ -n "${fq_2}" ]]; then arr_cmd+=( "${bam_rest}" ); fi
    arr_cmd+=( "${txt_unal}" )
    if [[ -n "${fq_2}" ]]; then arr_cmd+=( "${txt_imp}" ); fi
    arr_cmd+=( "${txt_mapq}" "${txt_diag}" )

    if ! \
        _run_or_print_cmd "${dry_run}" arr_cmd
    then
        echo_err_func "${func}" \
            "Step #4: failed to delete the diagnostic work file and the" \
            "read-name lists."
        return 1
    fi
}


# For '.diag.tsv.gz': count the templates and records of one category in a BAM
# file.
function _count_diag_cat() {
    local threads="${1:-}"
    local mode="${2:-}"
    local fil="${3:-}"
    local txt_tmp="${4:-}"
    local n_tpl n_rec
    local -a arr_flt
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _count_diag_cat
    [--help] threads mode fil txt_tmp [flt ...]

  Print the number of templates and of records in one '.diag.tsv.gz' category of a BAM file, tab-separated (for Step 5 of 'align_fastqs()' with '--unal_mode separate').

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  1  threads : int
    Number of threads to use for 'samtools view'.

  2  mode : {'rec', 'tpl'}
    'rec' counts the records that 'flt' selects and the templates holding them; 'tpl' counts the templates with a record that 'flt' selects and all their records.

  3  fil : file
    BAM file to count.

  4  txt_tmp : file
    Temporary file for the read names that mode 'tpl' selects; deleted before returning.

  5+  flt : str
    'samtools view' filter options, such as '-f 12' or '--expr', or none for every record.

Returns
-------
  Prints '<templates><tab><records>' and returns 0; returns 1 if a count fails.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - cut
    - rm
    - samtools
    - sort
    - wc

  Templates are counted by read name, so the file need not be grouped by name.

Examples
--------
  1. Count the pairs with both mates unaligned in a side file.
    '''bash
    _count_diag_cat 2 rec work/pairs.unal.bam work/pairs.tmp.txt -f 12
    '''

  2. Count the improper pairs whose mates map to different references.
    '''bash
    _count_diag_cat \\
        2 tpl work/pairs.improp.bam work/pairs.tmp.txt \\
        -F 2308 --expr 'refid != mrefid'
    '''
EOM
    )

    if [[ "${threads}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    arr_flt=( "${@:5}" )

    if [[ "${mode}" == "rec" ]]; then
        n_rec="$(
            samtools view -c -@ "${threads}" "${arr_flt[@]}" "${fil}"
        )" || return 1
        n_tpl="$(
            set -o pipefail
            samtools view -@ "${threads}" "${arr_flt[@]}" "${fil}" \
                | cut -f 1 \
                | LC_ALL=C sort -u \
                | wc -l
        )" || return 1
    else
        if ! (
            set -o pipefail
            samtools view -@ "${threads}" "${arr_flt[@]}" "${fil}" \
                | cut -f 1 \
                | LC_ALL=C sort -u \
                > "${txt_tmp}"
        ); then
            rm -f "${txt_tmp}"
            return 1
        fi

        n_tpl="$(wc -l < "${txt_tmp}")" || return 1
        n_rec="$(
            samtools view -c -@ "${threads}" -N "${txt_tmp}" "${fil}"
        )" || return 1
        rm -f "${txt_tmp}"
    fi

    printf '%s\t%s\n' "${n_tpl// /}" "${n_rec}"
}


# Step 5 with '--unal_mode separate': write '<stem>.diag.tsv.gz', the templates
# and records of each file and category.
function _write_diag_tsv() {
    local func="${1:-}"
    local dry_run="${2:-}"
    local threads="${3:-}"
    local aligner="${4:-}"
    local bt2_X="${5:-}"
    local fq_2="${6:-}"
    local mapq="${7:-}"
    local fil_wrk="${8:-}"
    local fil_unal="${9:-}"
    local fil_improp="${10:-}"
    local fil_mapq="${11:-}"
    local txt_all="${12:-}"
    local fil_tsv="${13:-}"
    local txt_tmp="${fil_tsv%.tsv.gz}.names.tmp"
    local tsv_tmp="${fil_tsv%.gz}.tmp"
    local rol fil cat cnt n_tpl n_rec n_dup n_pri tot_dup tot_pri
    local all_tpl all_rec ok_dup
    local -a arr_rol arr_row
    local -A arr_fil_rol
    local show_help

    show_help=$(cat << EOM
Usage
-----
  _write_diag_tsv
    [--help] func dry_run threads aligner bt2_X fq_2 mapq fil_wrk fil_unal fil_improp fil_mapq txt_all fil_tsv

  Write '.diag.tsv.gz' for '--unal_mode separate' (Step 5 of 'align_fastqs()').

  It holds a header line, then one row per alignment file and category with its templates and records, in the order that 'docs/design/align_fastqs_separate.md' gives.

  Single-end data has no pair rows and no improper file.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  01  func : str
    Name of the calling function for diagnostics.

  02  dry_run : bool
    Whether to print what the step writes instead of counting; the counts read files that a dry run does not write.

  03  threads : int
    Number of threads to use.

  04  aligner : {'bowtie2', 'bwa', 'bwa-mem2'}
    Aligner of the run, which decides which improper-pair categories can be measured.

  05  bt2_X : int
    Bowtie 2's maximum fragment length ('-X'), for the improper pairs too far apart.

  06  fq_2 : file
    Second FASTQ input file, or empty for single-end data.

  07  mapq : int
    MAPQ threshold, for the improper pairs that would also have failed it.

  08  fil_wrk : file
    Main BAM work file.

  09  fil_unal : file
    BAM work file of the unaligned reads and their mates.

  10  fil_improp : file
    BAM work file of the improper pairs; not read for single-end data.

  11  fil_mapq : file
    BAM work file of the low-MAPQ templates.

  12  txt_all : file
    Gzip-compressed duplicate statistics of Markdup 2.

  13  fil_tsv : file
    Gzip-compressed table to write.

Returns
-------
  Returns 0 after writing the table or printing what it writes for a dry run; prints an error and returns 1 otherwise.

Notes
-----
  Runtime requirements:
    - awk
    - bash >= 4.4
    - cut
    - gawk
    - gzip
    - mv
    - rm
    - samtools
    - sort
    - wc

  - The improper file's rows come from 'lib/awk/classify_improp_pairs.awk', in one pass.
  - The main file's 'duplicate_all' row comes from Markdup 2's statistics:
    + Its records are the duplicate total less the other files' duplicate records.
    + Its templates are the duplicate primary records less theirs, halved for paired-end data, since 'markdup' marks both mates of a duplicate pair.
    + If a count comes out negative, or the paired-end primary count is odd, the row is NA, with a warning.

Examples
--------
  1. Write the table for a paired-end run.
    '''bash
    _write_diag_tsv \\
        align_fastqs \\
        false \\
        2 \\
        bowtie2 \\
        1000 \\
        tests/fixtures/align_fastqs/fastq/pe/pairs_pe_R2.atria.fastq.gz \\
        1 \\
        work/pairs.bam \\
        work/pairs.unal.bam \\
        work/pairs.improp.bam \\
        work/pairs.mapq.bam \\
        work/pairs.all.markdup.txt.gz \\
        work/pairs.diag.tsv.gz
    '''

  2. Display the helper's argument contract.
    '''bash
    _write_diag_tsv --help
    '''
EOM
    )

    if [[ "${func}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    if [[ "${dry_run}" == "true" ]]; then
        printf '%s %q\n\n' \
            "# Count each file's templates and records by category into" \
            "${fil_tsv}"
        return 0
    fi

    arr_fil_rol=( [unal]="${fil_unal}" [mapq]="${fil_mapq}" )
    arr_rol=( unal )
    if [[ -n "${fq_2}" ]]; then
        arr_fil_rol[improp]="${fil_improp}"
        arr_rol+=( improp )
    fi
    arr_rol+=( mapq )

    # The main file holds whole templates only, so its templates are its
    # primary first reads (paired-end) or primary reads (single-end).
    if [[ -n "${fq_2}" ]]; then
        arr_row=( -f 64 -F 2304 )
    else
        arr_row=( -F 2304 )
    fi

    if ! {
        all_rec="$(samtools view -c -@ "${threads}" "${fil_wrk}")" \
            && all_tpl="$(
                samtools view -c -@ "${threads}" "${arr_row[@]}" "${fil_wrk}"
            )"
    }; then
        echo_err_func "${func}" \
            "Step #5: failed to count the main file's records."
        return 1
    fi

    # Markdup 2's duplicates in the main file: its totals less the diagnostic
    # files' duplicates.
    if ! \
        cnt="$(
            set -o pipefail
            gzip -cd "${txt_all}" \
                | awk -F ': ' '
                    $1 == "DUPLICATE TOTAL" { tot = $2 }
                    $1 == "DUPLICATE PRIMARY TOTAL" { pri = $2 }
                    END { print tot + 0, pri + 0 }
                '
        )"
    then
        echo_err_func "${func}" \
            "Step #5: failed to read the duplicate statistics of all files."
        return 1
    fi
    read -r tot_dup tot_pri <<< "${cnt}"

    for rol in "${arr_rol[@]}"; do
        fil="${arr_fil_rol[${rol}]}"
        if ! {
            n_dup="$(samtools view -c -@ "${threads}" -f 1024 "${fil}")" \
                && n_pri="$(
                    samtools view -c -@ "${threads}" -f 1024 -F 2304 "${fil}"
                )"
        }; then
            echo_err_func "${func}" \
                "Step #5: failed to count the duplicates in '${fil}'."
            return 1
        fi

        tot_dup=$(( tot_dup - n_dup ))
        tot_pri=$(( tot_pri - n_pri ))
    done

    # The subtraction can only go wrong through a counting defect or a change
    # in 'samtools markdup', never through the data. Its row is then NA, and
    # every other count stands.
    ok_dup=true
    if (( tot_dup < 0 || tot_pri < 0 )); then
        ok_dup=false
    elif [[ -n "${fq_2}" ]]; then
        if (( tot_pri % 2 != 0 )); then ok_dup=false; fi
        tot_pri=$(( tot_pri / 2 ))
    fi

    if [[ "${ok_dup}" == "false" ]]; then
        echo_warn_func "${func}" \
            "Step #5: the main file's Markdup 2 duplicates do not follow from" \
            "the statistics in '${txt_all}' and the other files' duplicates;" \
            "the main file's 'duplicate_all' row is NA."
        tot_dup="NA"
        tot_pri="NA"
    fi

    # Write the rows in the documented order: each file's 'all' row first, its
    # 'unjudged' and 'duplicate_all' rows last. Only the improper rows have the
    # low-MAPQ columns. The subshell stops at the first failed count.
    (
        printf '%s\t' file category templates records templates_mapq_low
        printf '%s\n' records_mapq_low

        cnt="$(
            _count_diag_cat \
                "${threads}" rec "${fil_wrk}" "${txt_tmp}" --rf 2308
        )" || exit 1
        read -r n_tpl n_rec <<< "${cnt}"
        printf 'main\tall\t%s\t%s\tNA\tNA\n' "${all_tpl}" "${all_rec}"
        printf 'main\tunjudged\t%s\t%s\tNA\tNA\n' "${n_tpl}" "${n_rec}"
        printf 'main\tduplicate_all\t%s\t%s\tNA\tNA\n' \
            "${tot_pri}" "${tot_dup}"

        for rol in "${arr_rol[@]}"; do
            fil="${arr_fil_rol[${rol}]}"

            # One 'gawk' pass counts every improper row.
            if [[ "${rol}" == "improp" ]]; then
                (
                    set -o pipefail
                    samtools view -@ "${threads}" "${fil}" \
                        | gawk -f "${_fil_awk_improp}" \
                            -v aligner="${aligner}" \
                            -v max_frag="${bt2_X}" \
                            -v min_frag=0 \
                            -v mapq="${mapq}"
                ) || exit 1
                continue
            fi

            cnt="$(_count_diag_cat "${threads}" rec "${fil}" "${txt_tmp}")" \
                || exit 1
            read -r all_tpl all_rec <<< "${cnt}"
            printf '%s\tall\t%s\t%s\tNA\tNA\n' \
                "${rol}" "${all_tpl}" "${all_rec}"

            if [[ "${rol}" == "unal" && -n "${fq_2}" ]]; then
                cnt="$(
                    _count_diag_cat \
                        "${threads}" rec "${fil}" "${txt_tmp}" -f 12
                )" || exit 1
                read -r n_tpl n_rec <<< "${cnt}"
                printf 'unal\tboth_unaligned\t%s\t%s\tNA\tNA\n' \
                    "${n_tpl}" "${n_rec}"
                printf 'unal\tone_unaligned\t%s\t%s\tNA\tNA\n' \
                    "$(( all_tpl - n_tpl ))" "$(( all_rec - n_rec ))"
            fi

            for cat in unjudged duplicate_all; do
                if [[ "${cat}" == "unjudged" ]]; then
                    arr_row=( --rf 2308 )
                else
                    arr_row=( -f 1024 )
                fi

                cnt="$(
                    _count_diag_cat \
                        "${threads}" rec "${fil}" "${txt_tmp}" "${arr_row[@]}"
                )" || exit 1
                printf '%s\t%s\t%s\tNA\tNA\n' "${rol}" "${cat}" "${cnt}"
            done
        done
    ) > "${tsv_tmp}" || {
        rm -f "${tsv_tmp}" "${txt_tmp}"
        echo_err_func "${func}" \
            "Step #5: failed to count the diagnostic files' templates and" \
            "records."
        return 1
    }

    if ! \
        gzip -n -c "${tsv_tmp}" > "${fil_tsv}"
    then
        rm -f "${tsv_tmp}" "${fil_tsv}"
        echo_err_func "${func}" \
            "Step #5: failed to compress the diagnostic table."
        return 1
    fi

    rm -f "${tsv_tmp}"
}
