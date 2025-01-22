#!/bin/bash

function calculate_factor_depth() {
    local n_in="${1}"      # Alignment count for input (not IP) BAM file
    local siz_bin="${2}"   # Bin size (in bp)
    local siz_gen="${3}"   # Effective genome size for model organism (in bp)
    local mode="${4}"      # "frag" or "norm"
    local rnd="${5:-24}"   # Number of decimal points for rounding
    local show_help        # Help message/documentation
    local num den fct_dep  # Variables for calculations

    show_help=$(cat << EOM
----------------------
calculate_factor_depth
----------------------

Description:
  Computes a depth factor based on the number of alignments in an input BAM 
  file, bin size, and effective genome size. The calculation supports both 
  fragment-length normalization ("frag") and "normalized coverage" ("norm";
  see PMID: 37160995 for more details).

Positional parameters:
  1, n_in    (int): Alignment count from the input BAM file.
  2, siz_bin (int): Bin size (in base pairs).
  3, siz_gen (int): Effective genome size (in base pairs).
  4, mode    (str): Mode of calculation; options: "frag" or "norm".
  5, rnd     (int): Number of decimal places for rounding (default: 24).

Returns:
  The computed depth factor.

Usage:
  calculate_factor_depth
      "\${n_in}" "\${siz_bin}" "\${siz_gen}"
      "\${mode}" "\${rnd}"

Examples:
  \`\`\`
  #  Compute depth factor for fragment-length normalized coverage
  calculate_factor_depth 12851824 20 12157105 "frag" 12

  #  Compute depth factor for "normalized coverage" 
  calculate_factor_depth 12851824 30 12157105 "norm" 24
  \`\`\`
EOM
    )

    #  Parse and check function parameters
    if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
        echo "${show_help}"
        return 0
    fi

    #  Validate 'n_in' is not zero
    if [[ ${n_in} -eq 0 ]]; then
        echo \
            "Error: Positional parameter 1, 'n_in', is zero.." >&2
        return 1
    fi

    #  Validate 'n_in' is a positive integer
    if [[ ! "${n_in}" =~ ^[1-9][0-9]*$ ]]; then
        echo \
            "Error: Positional parameter 1, 'n_in', must be a positive" \
            "integer: '${n_in}'." >&2
        return 1
    fi

    #  Validate 'siz_bin' is a positive integer
    if [[ ! "${siz_bin}" =~ ^[1-9][0-9]*$ ]]; then
        echo \
            "Error: Positional parameter 2, 'siz_bin', must be a positive" \
            "integer: '${siz_bin}'." >&2
        return 1
    fi

    #  Validate 'siz_gen' is a positive integer
    if [[ ! "${siz_gen}" =~ ^[1-9][0-9]*$ ]]; then
        echo \
            "Error: Positional parameter 3, 'siz_gen', must be a positive" \
            "integer: '${siz_gen}'." >&2
        return 1
    fi

    #  Validate 'mode' argument
    case "${mode}" in
        frag|norm) : ;;
        *)
            echo \
                "Error: Positional parameter 4, 'mode', must be 'frag' or" \
                "'norm': '${mode}'." >&2
            return 1
        ;;
    esac

    #  Validate 'rnd' is a positive integer (or zero for no decimal places)
    if [[ ! "${rnd}" =~ ^[0-9]+$ ]]; then
        echo \
            "Error: Positional parameter 5, 'rnd', must be a non-negative" \
            "integer: '${rnd}'." >&2
        return 1
    fi

    #  Compute numerator and denominator
    num=$(echo "scale=${rnd}; ((${n_in} * ${siz_bin}) / ${siz_gen})" | bc -l)
    den=$(echo "scale=${rnd}; (1 - (${siz_bin} / ${siz_gen}))" | bc -l)

    #  Compute depth factor
    fct_dep=$(echo "scale=${rnd}; ${num} / ${den}" | bc -l)

    #  Compute final output based on mode
    if [[ "${mode}" == "norm" ]]; then
        #  For "normalized coverage"
        echo "scale=${rnd}; ${fct_dep} / ${n_in}" | bc -l
    else
        #  For fragment-length normalized coverage
        echo "${fct_dep}"
    fi
}
