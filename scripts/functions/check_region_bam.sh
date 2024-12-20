#!/bin/bash

#  Function to check that a region is within the bounds of the reference genome
#+ based on information encoded in a BAM file; dependency: Samtools
function check_region_bam() {
    local bam="${1}"
    local region="${2}"

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
