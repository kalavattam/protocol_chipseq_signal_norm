#!/bin/bash

#  Function to validate the format of genomic region/range input
check_region() {
    local region="${1}"

    #  Define common components for regex patterns
    local br="I|II|III|IV|V|VI|VII|VIII|IX|X|XI|XII|XIII|XIV|XV|XVI"
    local pi="[0-9]+"

    #  Define regex patterns for valid formats
    local ro="^(${br})$"
    local pi="^${pi}$"
    local cr="^chr(${br})$"
    local ci="^chr${pi}$"
    local ra="^(${br}|chr(${br})|${pi}|chr${pi}):${pi}-${pi}$"

    #  Check against each pattern
    if [[
           "${region}" =~ ${ro} || "${region}" =~ ${pi} \
        || "${region}" =~ ${cr} || "${region}" =~ ${ci} \
        || "${region}" =~ ${ra}
    ]]; then
        return 0
    else
        cat << EOF >&2
Error: Invalid region format: '${region}'. Expected formats:
  - A single Roman numeral (e.g., I, II, III, ..., XVI).
  - A single positive integer (e.g., 1, 2, 3).
  - The string 'chr' followed by a Roman numeral (e.g., chrI, chrII, chrXVI,
    etc.).
  - The string 'chr' followed by a positive integer (e.g., chr1, chr2, chr20,
    etc.).
  - A genomic range in the format <chromosome>:<start_position>-<end_position>,
    where <start_position> and <end_position> are positive integers, and
    <end_position> is greater than <start_position>.
EOF
        return 1
    fi
}