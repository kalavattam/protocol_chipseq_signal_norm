#!/bin/bash

function count_alignments_bam() {
    local threads="${1}"          # Number of threads for parallelization
    local fil_in="${2}"           # Input (not IP) BAM file
    local fil_typ="${3:-paired}"  # "paired" or "single"
    local show_help               # Help message/documentation

    show_help=$(cat << EOM
--------------------
count_alignments_bam
--------------------

Description:
  Counts the number of alignments in a BAM file based on whether the data is 
  paired-end ("paired") or single-end ("single"). Uses 'samtools view' with 
  filtering expressions to count specific alignment flags.

Positional parameters:
  1, threads (int): Number of threads for 'samtools view'.
  2, fil_in  (str): Input (not IP) BAM file for which to count alignments.
  3, fil_typ (str): Alignment type; options: 'paired' or 'single' (default:
                    'paired').

Returns:
  An integer representing the count of alignments matching the given type.

Usage:
  count_alignments_bam "\${threads}" "\${fil_in}" "\${fil_typ}"

Examples:
  \`\`\`
  #  Count alignments in a BAM file of paired-end alignments using 8 threads
  count_alignments_bam 8 sample.bam paired

  #  Count alignments in a BAM file of single-end alignments using 4 threads
  count_alignments_bam 4 sample.bam single
  \`\`\`
EOM
    )

    #  Parse and check function parameter
    if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
        echo "${show_help}"
        return 0
    fi

    #  Validate 'threads' is a positive integer
    if [[ ! "${threads}" =~ ^[1-9][0-9]*$ ]]; then
        echo \
            "Error: Positional parameter 1, 'threads', must be a positive" \
            "integer: '${threads}'." >&2
        return 1
    fi

    #  Validate existence of input (not IP) BAM file
    if [[ ! -f "${fil_in}" ]]; then
        echo \
            "Error: Positional parameter 2, 'fil_in', not found:" \
            "'${fil_in}'." >&2
        return 1
    fi

    case "${fil_typ}" in
        single|se) fil_typ="single" ;;
        paired|pe) fil_typ="paired" ;;
        *)
            echo \
                "Error: Positional parameter 3, 'fil_typ', must be 'paired'" \
                "or 'single': '${fil_typ}'." >&2
            return 1
        ;;
    esac

    #  Count alignments based on alignment type
    if [[ "${fil_typ}" == "paired" ]]; then
        n_in=$(
            samtools view \
                -@ "${threads}" -c --expr "$({
                        echo '(flag == 99) || '   | tr -d '\n'
                        echo '(flag == 1123) || ' | tr -d '\n'
                        echo '(flag == 163) || '  | tr -d '\n'
                        echo '(flag == 1187)'
                })" \
                "${fil_in}"
        )
    elif [[ "${fil_typ}" == "single" ]]; then
        n_in=$(
            samtools view \
                -@ "${threads}" -c --expr "$({
                        echo '(flag == 0) || '    | tr -d '\n'
                        echo '(flag == 1024) || ' | tr -d '\n'
                        echo '(flag == 16) || '   | tr -d '\n'
                        echo '(flag == 1040)'
                })" \
                "${fil_in}"
        )
    fi

    #  Return count
    echo "${n_in}"
}
