#!/bin/bash

#  Validate the format of genomic region/range input
function check_region() {
    local region="${1:-}"

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


#  Check that a region is within the bounds of the reference genome based on
#+ information encoded in a BAM file; dependency: Samtools
function check_region_bam() {
    local bam="${1:-}"
    local region="${2:-}"

    #  Check that positional parameters are supplied
    if [[ -z "${bam}" ]]; then
        echo "Error: Positional parameter 1, \${bam}, not supplied." >&2
        return 1
    fi

    if [[ -z "${region}" ]]; then
        echo "Error: Positional parameter 2, \${region}, not supplied." >&2
        return 1
    fi

    #  Check existence of BAM file
    if [[ ! -f "${bam}" ]]; then
        echo "Error: BAM file '${bam}' does not exist." >&2
        return 1
    fi

    #  Check that region format is correct
    check_region "${region}" || return 1

    #  Extract reference genome sizes from BAM index, storing them in arrays
    local arr_chr=()
    local arr_siz=()
    while IFS=$'\t' read -r chr size _; do
        arr_chr+=( "${chr}" )
        arr_siz+=( "${size}" )
    done < <(samtools idxstats "${bam}")  #MAYBE Use AWK, not Samtools

    #  Parse the genomic region
    local chr start end
    if [[ "${region}" =~ ^([^:]+):([0-9]+)-([0-9]+)$ ]]; then
        #  Region includes a range
        chr="${BASH_REMATCH[1]}"
        start="${BASH_REMATCH[2]}"
        end="${BASH_REMATCH[3]}"
    else
        #  Region is just a chromosome name
        chr="${region}"
        start=1
        end=0
    fi

    #  Find the chromosome and its size
    local siz_chr=0
    for i in "${!arr_chr[@]}"; do
        if [[ "${arr_chr[$i]}" == "${chr}" ]]; then
            siz_chr="${arr_siz[i]}"
            break
        fi
    done

    #  Check that chromosome exists in BAM file
    if (( siz_chr == 0 )); then
        echo \
            "Error: Chromosome '${chr}' not found in the BAM file '${bam}'." \
            >&2
        return 1
    fi

    #  Set end to chromosome size if not explicitly specified
    if (( end == 0 )); then
        end="${siz_chr}"
    fi

    #  Check that start and end positions are within bounds
    if (( start < 1 || start > siz_chr || end < start || end > siz_chr )); then
        echo \
            "Error: Region '${region}' is out of bounds for chromosome" \
            "'${chr}' in BAM file '${bam}'." >&2
        return 1
    fi

    return 0
}
