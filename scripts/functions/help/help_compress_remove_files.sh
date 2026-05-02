#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: help_compress_remove_files.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
#
# Distributed under the MIT license.


# shellcheck disable=SC1111,SC2154
function help_compress_remove_files() {
    cat << EOM
Usage:
  compress_remove_files.sh
    [--help]
    [--threads <int>]
    --dir_fnd <dir> --pattern <str> --size <int>
    [--depth <int>] [--include <csv:str>] [--exclude <csv:str>]
    [--chk_con] [--chk_exc]


Description:
  Find files in a specified directory that match the given pattern, compress files larger than the specified size, and delete files that are 0 in size.

  Files can be further included or excluded based on comma-separated name patterns, and the search can be limited to a specified directory depth.


Arguments:
  -h, --hlp, --help  <flag>
    Display this help message and exit.

  -t, --thr, --threads  <int>
    Number of threads to use (default: ${threads}).

  -df, --dir_fnd  <dir>
    Directory to search with the program 'find'.

  -pa, --pttrn, --pattern  <str>
    File pattern, including shell wildcard characters, used in the construction of the underlying 'find' call (default: '${pattern}').

  -sz, --size  <int>
    Minimum size in kilobytes for compression (default: ${size}).

  -de, --dpth, --depth  <int>
    Maximum depth to search within the directory.

  -in, --incld, --include  <csv:str>
    Comma-separated list of patterns to include with respect to '--pattern', including shell wildcards. '--include' is subordinate to '--pattern'.

  -ex, --excld, --exclude  <csv:str>
    Comma-separated list of patterns to exclude with respect to '--pattern', including shell wildcards (default: '${exclude}'). '--exclude' is subordinate to '--pattern'.

  -cn, --chk_con  <flag>
    Check the construction of the find command and exit.

  -ce, -cu, --chk_exc, --chk_exu  <flag>
    Check the construction and execution of the find command and exit.


Dependencies:
  External programs:
    - Bash >= 4.4
    - find
    - GNU Parallel, when '--threads > 1'
    - gzip
    - realpath
    - rm, when '--threads > 1'
    - sort

  Sourced function scripts:
    - source_helpers.sh
      + source_helpers
    - check_args.sh
      + check_flags_mut_excl
      + require_optarg
    - check_env.sh
      + check_env_installed
      + check_pgrm_path
    - check_inputs.sh
      + validate_var
      + validate_var_dir
    - check_numbers.sh
      + check_int_pos
    - construct_find.sh
      + construct_find
    - format_outputs.sh
      + echo_err
      + echo_warn
      + print_cmd_array
    - handle_env.sh
      + handle_env
    - help/help_compress_remove_files.sh
      + help_compress_remove_files


Notes:
  - This script does not handle logical OR operations, just AND and AND NOT.
  - Script will exit with an error message if it is run from the target directory being searched or from one of its subdirectories. This is a conservative safety policy.
  - If the string assigned to '--pattern' matches anything in the current working directory, the script may emit a warning as a reminder that it may be running from an inconvenient location.
  - If '--threads' is assigned a positive integer greater than 1, then the script will use GNU Parallel to parallelize file handling and processing.


Examples:
  1. Compress matching log files larger than 2 KB.
    '''bash
    bash compress_remove_files.sh
      --dir_fnd "\${HOME}/path/to/dir"
      --pattern "*"
      --size 2
      --depth 2
      --include "*.log"
    '''

  2. Check the constructed 'find' command without executing it.
    '''bash
    bash compress_remove_files.sh
      --dir_fnd "./err_out"
      --pattern "*"
      --size 2
      --chk_con
    '''

  3. Run with two threads.
    '''bash
    bash compress_remove_files.sh
      --threads 2
      --dir_fnd "./err_out"
      --pattern "*"
      --size 2
    '''


#TODO:
  - Support size thresholds below 1 KB and non-integer size thresholds, e.g., 0.5 KB.
    + Once completed, change from <int> to <size>.
EOM
}
