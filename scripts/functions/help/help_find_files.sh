#!/bin/bash
# -*- coding: utf-8 -*-
#
# Script: help_find_files.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
#
# Distributed under the MIT license.


# shellcheck disable=SC1111
function help_find_files() {
    cat << EOM
Usage:
  find_files.sh
    [--help] --dir_fnd <str> --pattern <str> [--depth <int>] [--follow] [--fastqs] [--include <str>] [--exclude <str>] [--chk_con] [--chk_exc]

Description:
  Search for files in a specified directory using the *nix 'find' command. By default, returns the results as a single comma-separated string to stdout.

Arguments:
  -h, --help
    Display this help message and exit.

  -df, --dir_fnd
    Directory to search with the program 'find'.

  -pa, --pattern
    Primary file pattern, including shell wildcard characters, used in the construction of the underlying 'find' call.

  -de, --depth
    Maximum depth to search within the directory.

  -fl, -sy, --follow, --symlink
    Follow symbolic links during the search.

  -fq, --fastqs
    Find FASTQ files, returning them in a semicolon-separated string with paired-end files grouped into comma-separated substrings.

  -in, --include
    Comma-separated list of patterns to include, including shell wildcards. '--include' is subordinate to '--pattern'.

  -ex, --exclude
    Comma-separated list of patterns to exclude, including shell wildcards. '--exclude' is subordinate to '--pattern'.

  -cn, --chk_con
    Check the construction of the 'find' command and exit.

  -ce, --chk_exc
    Check the construction and execution of the 'find' command, output raw list of found files, and then exit.

  -ss, --shw_str
    With '--chk_exc', show a one-line “stringified” summary in addition to the raw list; without '--chk_exc', ignored.

Dependencies:
  - Programs
    + Bash >= 4.4
    + find
    + paste
    + sed
    + sort
    + tr
  - Functions  #TODO: update, bringing in line with, e.g., 'execute_calculate_scaling_factor.sh' docs
    + check_arg_supplied     ## check_args ##
    + check_file_dir_exists  ## check_inputs ##
    + check_flags_mut_excl   ## check_args ##
    + check_int_pos          ## check_numbers ##
    + check_pgrm_path        ## check_env ##
    + construct_find         ## construct_find ##
    + echo_error             ## format_outputs ##
    + echo_warning           ## format_outputs ##
    + help/help_find_files
    + pair_fastqs            ## process_sequences ##
    + print_banner_pretty    ## format_outputs ##
    + print_cmd_array        ## format_outputs ##

Notes:
  - This script does not handle logical 'OR' operations, just 'AND' and 'AND NOT'.
  - 'find_files.sh' will exit with an error message if it is run from the target directory being searched or from one of its subdirectories. This is a conservative safety policy.
  - If the string assigned to '--pattern' matches anything in the current working directory, the script may emit a warning as a reminder that it may be running from an inconvenient location.

Examples:
  1. Search for specific BAM files
  '''
  bash \${HOME}/path/to/find_files.sh
      --dir_fnd "\${HOME}/path/to/directory"
      --pattern "*.bam"
      --depth 1
      --include "*Hho1*,*Q*"
      --exclude "*Hmo1*,*G2M*,*G1*"
  '''

  2. Search for symlinked FASTQ files
  '''
  bash \${HOME}/path/to/find_files.sh
      --dir_fnd "\${HOME}/path/to/another/directory"
      --pattern "*.fastq.gz"
      --follow
      --fastqs
  '''
EOM
}
