#!/bin/bash

function calculate_frag_avg() {
    local threads="${1}"      # Number of threads for parallelization
    local fil="${2}"          # Input BAM file
    local fil_typ="${3:-pe}"  # "paired", "pe", "single", or "se"
    local expr=""             # Samtools filtration expression 
    local show_help           # Help message/documentation

    show_help=$(cat << EOM
------------------
calculate_frag_avg
------------------

Description:
  Computes the average fragment length from a BAM file based on whether the 
  data is paired-end ("paired") or single-end ("single"). Uses 'samtools view'
  with filtering expressions and 'awk' to process fragment lengths.

Positional parameters:
  1, threads (int): Number of threads for 'samtools view'.
  2, fil     (str): Input BAM file for which to compute fragment lengths.
  3, fil_typ (str): Alignment type; options: 'paired' or 'single' (default:
                    'paired').

Returns:
  A floating-point value representing the average fragment length.

Usage:
  calculate_average_fragment_length "\${threads}" "\${fil}" "\${fil_typ}"

Examples:
  \`\`\`
  #  Compute average fragment length for paired-end alignments using 8 threads
  calculate_average_fragment_length 8 sample.bam paired

  #  Compute average fragment length for single-end alignments using 4 threads
  calculate_average_fragment_length 4 sample.bam single
  \`\`\`
EOM
    )

    #  Parse and check function parameters
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

    #  Validate existence of input BAM file
    if [[ ! -f "${fil}" ]]; then
        echo \
            "Error: Positional parameter 2, 'fil', not found:" \
            "'${fil}'." >&2
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

    #  Determine filtering flags based on alignment type
    if [[ "${fil_typ}" == "paired" ]]; then
        expr="$({
            echo '(flag == 99) || '   | tr -d '\n'
            echo '(flag == 1123) || ' | tr -d '\n'
            echo '(flag == 163) || '  | tr -d '\n'
            echo '(flag == 1187)'
        })"
    elif [[ "${fil_typ}" == "single" ]]; then
        expr="$({
            echo '(flag == 0) || '    | tr -d '\n'
            echo '(flag == 1024) || ' | tr -d '\n'
            echo '(flag == 16) || '   | tr -d '\n'
            echo '(flag == 1040)'
        })"
    fi

    #  Compute average fragment length using samtools + awk
    samtools view -@ "${threads}" --expr "${expr}" "${fil}" \
        | awk '{
            if ($9 > 0) { sum += $9; count++ }
        } END {
            if (count > 0) { print sum / count }
            else {
                print "Error: No valid fragment lengths found." > "/dev/stderr"
                exit 1
            }
        }'
}
