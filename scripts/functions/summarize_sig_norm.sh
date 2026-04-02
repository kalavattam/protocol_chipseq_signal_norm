#!/bin/bash
# -*- coding: utf-8 -*-
# 
# Script: summarize_sig_norm.sh
#
# Copyright 2025-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5-series models) was used in development.
# 
# Distributed under the MIT license.

function summarize_sig_norm() {
    local typ_sig="${1:-}"  # Type of signal computation
    local scl_fct="${2:-}"  # Scaling factor
    local typ_sig_lc        # Lowercase-converted signal type for case matching
    local mth_nrm           # Normalization method message
    local src_scl           # Scaling-factor message
    local show_help         # Help message

    show_help=$(cat << EOM
Usage:
  summarize_sig_norm [-h|--hlp|--help] typ_sig [scl_fct]

Description:
  Summarize resolved signal-normalization and scaling states for the signal-computation workflow.

  Prints a short human-readable summary indicating (i) what normalization mode is implied by 'typ_sig' and (ii) whether multiplicative scaling factors were supplied.

Positional arguments:
  1  typ_sig  <str>  Type of signal computation; e.g., 'unadj', 'frag', or 'norm' (aliases accepted).
  2  scl_fct  <str>  Scaling factor string (optional). If empty, assumes no explicit '--scl_fct' was supplied.

Returns:
  0 after printing the summary to stdout; 1 if required positional argument 1, 'typ_sig', is missing.

Examples:
  1. Summarize normalized coverage with no scaling factors
  '''bash
  summarize_sig_norm "norm"
  '''

  2. Summarize fragment-length normalization with explicit scaling
  '''bash
  summarize_sig_norm "frag" "1.25"
  '''

  3. Show function help
  '''bash
  summarize_sig_norm --help
  '''
EOM
    )

    #  Parse and check function arguments
    if [[ "${typ_sig}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${typ_sig}" ]]; then
        echo "Error: Positional argument 1, 'typ_sig', is missing." >&2
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  Lowercase-convert signal-type input so case matching is case-insensitive
    typ_sig_lc=$(printf '%s' "${typ_sig}" | tr '[:upper:]' '[:lower:]')

    #  Determine normalization method message
    case "${typ_sig_lc}" in
        u|unadj|unadjusted|s|smp|simple|r|raw)
            mth_nrm="No normalization; returning unadjusted signal:"
            mth_nrm+=" '--method ${typ_sig}'."
            ;;

        f|frg|frag|frg_len|frag_len|l|len|len_frg|len_frag)
            mth_nrm="Performing fragment-length normalization:"
            mth_nrm+=" '--method ${typ_sig}'."
            ;;

        n|nrm|norm|normalized)
            mth_nrm="Generating normalized coverage (Dickson et al., Sci Rep"
            mth_nrm+=" 2023): '--method ${typ_sig}'."
            ;;

        *)
            mth_nrm="Unknown normalization method: '--method ${typ_sig}'."
            ;;
    esac

    #  Determine scaling factor message
    if [[ -n "${scl_fct}" ]]; then
        src_scl="Custom multiplicative scaling factor(s):"
        src_scl+=" '--scl_fct ${scl_fct}'."
    else
        src_scl="No multiplicative scaling factor(s)."
    fi

    #  Print resolved argument states
    echo "#################################################"
    echo "## Summary of signal normalization and scaling ##"
    echo "#################################################"
    echo
    echo "- Normalization method: ${mth_nrm}"
    echo "- Scaling factor source: ${src_scl}"
    echo
    echo
}
