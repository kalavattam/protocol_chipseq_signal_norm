#!/bin/bash
# -*- coding: utf-8 -*-
# 
# Script: reset_max_job.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
# 
# Distributed under the MIT license.

#  Reset 'max_job', the maximum number of jobs to be run by SLURM at one time,
#+ if it exceeds the number of input files
function reset_max_job() {
    local max_job="${1:-}"
    local num_fil="${2:-}"

    if ! [[ "${max_job}" =~ ^[0-9]+$ ]]; then
        echo "Error: 'max_job' must be a non-negative integer." >&2
        return 1
    fi

    if ! [[ "${num_fil}" =~ ^[0-9]+$ ]]; then
        echo "Error: 'num_fil' must be a non-negative integer." >&2
        return 1
    fi

    if [[ "${max_job}" -gt "${num_fil}" ]]; then
        #  Cap concurrent jobs at the number of input files
        echo "${num_fil}"
    else
        #  Preserve user-supplied limit when it does not exceed file count
        echo "${max_job}"
    fi
}
