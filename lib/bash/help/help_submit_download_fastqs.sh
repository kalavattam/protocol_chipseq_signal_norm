#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: help_submit_download_fastqs.sh
#
# Copyright 2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT and Codex (GPT-5.5, GPT-5.6) were used in design, development,
# and documentation, with all output reviewed, edited, and approved by the
# author.
#
# Distributed under the MIT license.


function help_submit_download_fastqs() {
    cat << EOM
Usage
-----
  submit_download_fastqs.sh
    [--help]
    --dir_scr <dir>
    srr url_1 url_2 dir_out dir_sym nam_cus dir_eo nam_job

  Download one single-end or paired-end FASTQ entry and create custom symlink(s).

Parameters
----------
  -h, --help : flag
    Display this help message and exit.

  -ds, --dir_scr : dir
    Maintained entrypoint directory. Pass the repository 'bin' directory; shared libraries are resolved from adjacent 'lib/bash'.

  1  srr : str
    NCBI SRA database run accession code.

  2  url_1 : str
    URL (FTP or HTTPS) for FASTQ file.

  3  url_2 : str
    Second FASTQ URL for paired-end data ("NA" for single-end data).

  4  dir_out : dir
    Output directory. Directory to save FASTQ file(s).

  5  dir_sym : dir
    Directory for symlink(s) to FASTQ file(s).

  6  nam_cus : str
    Custom name for symlink(s).

  7  dir_eo : dir
    Directory for stderr and stdout log files.

  8  nam_job : str
    Job name.

Notes
-----
  Runtime requirements:
    - basename
    - bash >= 4.4
    - cat
    - dirname
    - ln
    - Network access
    - wget

  - Use 'NA' for 'url_2' when processing single-end data.

Examples
--------
  1. Download one single-end FASTQ and create a custom symlink.
    '''bash
    bash "\${dir_scr}/submit_download_fastqs.sh" \\
        --dir_scr "\${dir_scr}" \\
        SRR_SINGLE "\${url_single}" NA \\
        "\${dir_out}" "\${dir_sym}" sample_single "\${dir_eo}" download_fastqs
    '''

  2. Download one paired-end FASTQ pair and create custom symlinks.
    '''bash
    bash "\${dir_scr}/submit_download_fastqs.sh" \\
        --dir_scr "\${dir_scr}" \\
        SRR_PAIRED "\${url_r1}" "\${url_r2}" \\
        "\${dir_out}" "\${dir_sym}" sample_paired "\${dir_eo}" download_fastqs
    '''
EOM
}
