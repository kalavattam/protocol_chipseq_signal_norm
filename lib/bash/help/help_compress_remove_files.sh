#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: help_compress_remove_files.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-4- and GPT-5-series models; most recent:
# GPT-5.6) were used in design, development, and documentation, with all output
# reviewed, edited, and approved by the author.
#
# Distributed under the MIT license.


# TODO: Support size thresholds below 1 KB and non-integer thresholds. Once
# supported, document '--size' as type 'size' rather than 'int'.
# shellcheck disable=SC1111,SC2154
function help_compress_remove_files() {
    cat << EOM
Usage
-----
  compress_remove_files.sh
    [--help] [--dry_run]
    [--threads <int>]
    --dir_fnd <dir> --pattern <str> --size <int>
    [--depth <int>] [--include <csv>] [--exclude <csv>]
    [--chk_con]


  Find files in a specified directory that match the given pattern, compress files larger than the specified size, and delete files that are 0 in size.

  Files can be further included or excluded based on comma-separated name patterns, and the search can be limited to a specified directory depth.


Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  -dr, --dry_run : flag
    Run script in dry-run mode. Show the constructed commands and list their matching files, then exit without compressing or deleting files.

  -t, --thr, --threads : int
    Number of threads to use (default: ${threads}).

  -df, --dir_fnd : dir
    Directory to search with the program 'find'.

  -pa, --pttrn, --pattern : str
    File pattern, including shell wildcard characters, used in the construction of the underlying 'find' call (default: '${pattern}').

  -sz, --size : int
    Minimum size in kilobytes for compression (default: ${size}).

  -de, --dpth, --depth : int
    Maximum depth to search within the directory.

  -in, --incld, --include : list of str
    Comma-separated list of patterns to include with respect to '--pattern', including shell wildcards. '--include' is subordinate to '--pattern'.

  -ex, --excld, --exclude : list of str
    Comma-separated list of patterns to exclude with respect to '--pattern', including shell wildcards (default: '${exclude}'). '--exclude' is subordinate to '--pattern'.

  -cn, --chk_con : flag
    Check the construction of the find command and exit.

Notes
-----
  Runtime requirements:
    - bash >= 4.4
    - find
    - gzip
    - parallel (when '--threads' is greater than 1)
    - realpath
    - rm (when '--threads' is greater than 1)
    - sort

  - This script does not handle logical OR operations, just AND and AND NOT.
  - Script will exit with an error message if it is run from the target directory being searched or from one of its subdirectories. This is a conservative safety policy.
  - If the string assigned to '--pattern' matches anything in the current working directory, the script emits a warning as a reminder that it is running from an inconvenient location.
  - If '--threads' is assigned a positive integer greater than 1, then the script will use GNU Parallel to parallelize file handling and processing.


Examples
--------
  1. Compress matching log files larger than 2 KB.
    '''bash
    bash compress_remove_files.sh \\
        --dir_fnd "\${HOME}/path/to/dir" \\
        --pattern "*" \\
        --size 2 \\
        --depth 2 \\
        --include "*.log"
    '''

  2. Check the constructed 'find' command without executing it.
    '''bash
    bash compress_remove_files.sh \\
        --dir_fnd "\${HOME}/path/to/dir_eo" \\
        --pattern "*" \\
        --size 2 \\
        --chk_con
    '''

  3. Run with two threads.
    '''bash
    bash compress_remove_files.sh \\
        --threads 2 \\
        --dir_fnd "\${HOME}/path/to/dir_eo" \\
        --pattern "*" \\
        --size 2
    '''

  4. Preview matching files without compressing or deleting them.
    '''bash
    bash compress_remove_files.sh \\
        --dir_fnd "\${HOME}/path/to/dir_eo" \\
        --pattern "*" \\
        --size 2 \\
        --dry_run
    '''

EOM
}
