#!/bin/bash

function filter_bam_sp() {
    local threads=1
    local infile=""
    local outfile=""
    local mito=false
    local tg=false
    local mtr=false
    local chk_chr=false
    local show_help

    show_help=$(cat << EOM
-------------
filter_bam_sp
-------------

Description:
  Filter and reheader a BAM file for S. pombe chromosomes.

Keyword parameters:
   -t, --threads  (int): Number of threads to use (default: ${threads}).
   -i, --infile   (str): Coordinate-sorted BAM infile.
   -o, --outfile  (str): Filtered BAM outfile.
   -m, --mito    (flag): Retain SP_Mito chromosome (optional).
  -tg, --tg      (flag): Retain SP_II_TG chromosome (optional).
  -mr, --mtr     (flag): Retain SP_MTR chromosome (optional).
  -cc, --chk_chr (flag): Check chromosomes in filtered BAM outfile (optional).

Returns:
  Creates a BAM outfile filtered and reheadered for S. pombe chromosomes at the
  specified path.

Dependencies:
  - Programs
    + awk
    + Bash or Zsh
    + grep
    + mv
    + rm
    + Samtools

Example:
  \`\`\`
  #TODO
  \`\`\`
EOM
    )

    #  Parse keyword parameters
    if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
        echo "${show_help}"
        return 0
    fi

    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
             -t|--threads) threads="${2}"; shift 2 ;;
             -i|--infile)  infile="${2}";  shift 2 ;;
             -o|--outfile) outfile="${2}"; shift 2 ;;
             -m|--mito)    mito=true;      shift 1 ;;
            -tg|--tg)      tg=true;        shift 1 ;;
            -mr|--mtr)     mtr=true;       shift 1 ;;
            -cc|--chk_chr) chk_chr=true;   shift 2 ;;
            *)
                echo "## Unknown parameter passed: ${1} ##" >&2
                echo "" >&2
                show_help >&2
                exit 1
                ;;
        esac
    done

    #  Validate keyword parameters
    if [[ -z "${threads}" ]]; then
        echo "Error: --threads is required." >&2
        return 1
    fi

    if [[ ! "${threads}" =~ ^[1-9][0-9]*$ ]]; then
        echo \
            "Error: --threads was assigned '${threads}' but must be a" \
            "positive integer greater than or equal to 1." >&2
        return 1
    fi

    if [[ -z "${infile}" ]]; then
        echo "Error: --infile is a required parameter." >&2
        return 1
    fi

    if [[ ! -f "${infile}" ]]; then
        echo \
            "Error: File associated with --infile does not exist: ${infile}." \
            >&2
        return 1
    fi

    if [[ -z "${outfile}" ]]; then
        echo "Error: --outfile is a required parameter." >&2
        return 1
    fi

    if [[ ! -d "$(dirname "${outfile}")" ]]; then
        echo \
            "Error: Directory associated with --outfile does not exist:" \
            "$(dirname "${outfile}")" >&2
        return 1
    fi

    #  Do the main work
    #  Set the chromosomes to retain
    chromosomes="SP_I SP_II SP_III"
    if ${tg};   then chromosomes="SP_II_TG ${chromosomes}"; fi
    if ${mtr};  then chromosomes="${chromosomes} SP_MTR";   fi
    if ${mito}; then chromosomes="${chromosomes} SP_Mito";  fi

    # #  Convert the space-separated string of chromosomes to a pipe-separated
    # #+ string for use in grep pattern matching for retained chromosomes and
    # #+ program (@PG) lines
    # # shellcheck disable=SC2001,SC2028
    # pattern="^@HD|^@SQ.*SN:($(echo "${chromosomes}" | tr ' ' '|'))|^@PG|^[^@]"

    #  Define variables for subsequent filtering and reheading
    outdir="$(dirname "${outfile}")"
    insam="$(basename "${infile/.bam/.sam}")"
    outsam="$(basename "${outfile/.bam/.sam}")"
    path_insam="${outdir}/${insam}"
    path_outsam="${outdir}/${outsam}"

    #  Filter BAM infile, writing to BAM outfile
    #  First, set a trap to remove intermediate files when and however the
    #+ function exits
    trap 'rm -f "${path_insam}" "${path_outsam}"' EXIT

    #  Begin by converting BAM to SAM while including the current header in the
    #+ SAM
    # shellcheck disable=SC2086
    if ! \
        samtools view \
            -@ "${threads}" \
            -h \
            -o "${path_insam}" \
            "${infile}"
    then
        echo "Error: Failed to generate ${path_insam}." >&2
        return 1
    fi

    #  Filter pertinent lines from the SAM header, and perform filtering for
    #+ lines that meet the following conditions: (1) They do not begin with @
    #+ (i.e., are not in the header) and (2) they contain strings in
    #+ ${chromosomes} within field 3
    # shellcheck disable=SC2002
    if ! \
        cat "${path_insam}" \
            | awk \
                -v chromosomes="${chromosomes}" \
                'BEGIN {
                    split(chromosomes, chrom_arr, " ")
                    for (i in chrom_arr) {
                        chrom_map[chrom_arr[i]] = 1
                    }
                }

                function print_header_line() {
                    if ($0 ~ /^@SQ/ && $0 ~ /SN:([^ \t]+)/) {
                        split($0, a, /SN:/)
                        split(a[2], b, /[ \t]/)
                        sn = b[1]
                        if (!(sn in chrom_map)) {
                            next
                        }
                    }
                    print $0;
                }

                {
                    if ($0 ~ /^@/) {
                        print_header_line()
                    } else {
                        if ($3 in chrom_map) {
                            print $0
                        }
                    }
                }' \
                    > "${path_outsam}"
    then
        echo "Error: Failed to generate ${path_outsam}." >&2
        return 1
    fi

    # if ! \
    #     cat "${path_insam}" \
    #         | grep -E "${pattern}" \
    #         | awk \
    #             -v chromosomes="${chromosomes}" \
    #             'BEGIN {
    #                 split(chromosomes, chrom_arr, " ");
    #                 for (i in chrom_arr) {
    #                     chrom_map[chrom_arr[i]] = 1;
    #                 }
    #             }
    #
    #             {
    #                 if ($0 ~ /^@/) {
    #                     print $0
    #                 } else {
    #                     if ($3 in chrom_map) {
    #                         print $0
    #                     }
    #                 }
    #             }' \
    #                 > "${path_outsam}"
    # then
    #     echo "Error: Failed to generate ${path_outsam}." >&2
    #     return 1
    # fi

    #  Perform SAM-to-BAM conversion for filtered, reheadered alignment
    #+ information
    if ! \
        samtools view -b "${path_outsam}" > "${outfile}"
    then
        echo "Error: Failed to generate ${outfile}." >&2
        return 1
    fi

    #  Index the BAM outfile
    samtools index -@ "${threads}" "${outfile}"

    #  Remove intermediate SAM file if present
    if [[ -f "${outfile}" ]]; then
        if [[ -f "${path_insam}" ]]; then
            rm "${path_insam}"
        fi
    
        if [[ -f "${path_outsam}" ]]; then
            rm "${path_outsam}"
        fi
    fi
    
    #  Check unique entries in field 3 for lines not beginning with @
    if ${chk_chr}; then
        # shellcheck disable=SC2002
        samtools view -h "${outfile}" \
            | awk '!/^@/ { print $3 }' \
            | sort \
            | uniq
    fi
}
