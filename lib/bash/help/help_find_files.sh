#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: help_find_files.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
# GPT-5.6) were used in design, development, and documentation, with all output
# reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


# TODO FIXME: under Usage, do a sensical breakdown of option types per line.
# Maybe like this? (Entails reordering options.)
#   [--help]
#    --dir_fnd <dir> --pattern <str> [--include <csv>] [--exclude <csv>]
#    [--depth <int>] [--follow] [--fastqs]
#    [--chk_con] [--chk_exc]
# shellcheck disable=SC1111
function help_find_files() {
    cat >&2 << EOM
Usage
-----
  find_files.sh
    [--help] --dir_fnd <dir> --pattern <str> [--depth <int>] [--follow] [--fastqs] [--include <csv>] [--exclude <csv>] [--chk_con] [--chk_exc]

  Search for files in a specified directory using the *nix 'find' command. By default, returns the results as a single comma-separated string to stdout.

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  -df, --dir_fnd : dir
    Directory to search with the program 'find'.

  -pa, --pttrn, --pattern : str
    Primary file pattern, including shell wildcard characters, used in the construction of the underlying 'find' call.

  -de, --dpth, --depth : int
    Maximum depth to search within the directory.

  -fl, -sy, --fllw, --follow, --symlink : flag
    Follow symbolic links during the search.

  -fq, --fqs, --fastqs : flag
    Find FASTQ files, returning them in a semicolon-separated string with paired-end files grouped into comma-separated substrings.

  -in, --incld, --include : csv
    Comma-separated list of patterns to include, including shell wildcards. '--include' is subordinate to '--pattern'.

  -ex, --excld, --exclude : csv
    Comma-separated list of patterns to exclude, including shell wildcards. '--exclude' is subordinate to '--pattern'.

  -cn, --chk_con : flag
    Check the construction of the 'find' command and exit.

  -ce, -cu, --chk_exc, --chk_exu : flag
    Check the construction and execution of the 'find' command, output raw list of found files, and then exit.

  -ss, --shw_str, --show_str : flag
    With '--chk_exc', show a one-line “stringified” summary in addition to the raw list; without '--chk_exc', ignored.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - find
    - paste
    - sed
    - sort
    - tr

  - This script does not handle logical 'OR' operations, just 'AND' and 'AND NOT'.
  - 'find_files.sh' will exit with an error message if it is run from the target directory being searched or from one of its subdirectories. This is a conservative safety policy.
  - If the string assigned to '--pattern' matches anything in the current working directory, the script may emit a warning as a reminder that it may be running from an inconvenient location.

Examples
--------
  1. Search one level for BAM files that satisfy include and exclude selectors.
    '''bash
    bash bin/find_files.sh \\
        --dir_fnd "\${HOME}/path/to/directory" \\
        --pattern '*.bam' \\
        --depth 1 \\
        --include '*Hho1*,*Q*' \\
        --exclude '*Hmo1*,*G2M*,*G1*'
    '''

  2. Follow symlinks and serialize paired FASTQ files by sample.
    '''bash
    bash bin/find_files.sh \\
        --dir_fnd "\${HOME}/path/to/another/directory" \\
        --pattern '*.fastq.gz' \\
        --follow \\
        --fastqs
    '''
EOM
}
