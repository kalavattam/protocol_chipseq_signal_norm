#!/bin/bash

function filter_bam_sc() {
    local threads=1
    local infile=""
    local outfile=""
    local mito=false
    local chk_chr=false
    local show_help

    show_help=$(cat << EOM
Usage:
  filter_bam_sc
    [--help] [--threads <int>] --infile <str> --outfile <str> [--mito] [--chk_chr]

Description:
  Filter and reheader a BAM file for S. cerevisiae chromosomes.

Keyword arguments:
   -t, --threads  <int>  Number of threads to use (default: ${threads}).
   -i, --infile   <str>  Coordinate-sorted BAM infile.
   -o, --outfile  <str>  Filtered BAM outfile.
   -m, --mito     <flg>  Retain mitochondrial chromosome (optional).
  -cc, --chk_chr  <flg>  Check chromosomes in filtered BAM outfile (optional).

Returns:
  Creates a BAM outfile filtered and reheadered for S. cerevisiae chromosomes at the specified path.

Dependencies:
  - Programs
    + Bash or Zsh
    + grep
    + mv
    + rm
    + Samtools

Examples:
  '''bash
  #TODO
  '''

#TODO:
  - Somewhere and somehow, need to handle more than S. cerevisiae as "main" organism.
  - The filter activity needs to be recorded in the BAM header.
  - Support for SAM and CRAM too.
EOM
    )

    #  Parse keyword arguments
    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
             -t|--threads) threads="${2}"; shift 2 ;;
             -i|--infile)  infile="${2}";  shift 2 ;;
             -o|--outfile) outfile="${2}"; shift 2 ;;
             -m|--mito)    mito=true;      shift 1 ;;
            -cc|--chk_chr) chk_chr=true;   shift 1 ;;
            *)
                echo "## Unknown argument passed: '${1}' ##" >&2
                echo "" >&2
                echo "${show_help}" >&2
                exit 1
                ;;
        esac
    done

    #  Validate keyword arguments
    if [[ -z "${threads}" ]]; then
        echo "Error: '--threads' is required." >&2
        return 1
    fi

    if [[ ! "${threads}" =~ ^[1-9][0-9]*$ ]]; then
        echo \
            "Error: '--threads' was assigned '${threads}' but must be a" \
            "positive integer greater than or equal to 1." >&2
        return 1
    fi

    if [[ -z "${infile}" ]]; then
        echo "Error: '--infile' is a required argument." >&2
        return 1
    fi

    if [[ ! -f "${infile}" ]]; then
        echo \
            "Error: File associated with '--infile' does not exist:" \
            "${infile}." >&2
        return 1
    fi

    if [[ -z "${outfile}" ]]; then
        echo "Error: '--outfile' is a required argument." >&2
        return 1
    fi

    if [[ ! -d "$(dirname "${outfile}")" ]]; then
        echo \
            "Error: Directory associated with '--outfile' does not exist:" \
            "$(dirname "${outfile}")" >&2
        return 1
    fi

    #  Do the main work
    #TODO: initialize local variables first
    #  Set the chromosomes to retain
    chromosomes="I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI"
    if ${mito}; then chromosomes="${chromosomes} Mito"; fi

    #  Convert the space-separated string of chromosomes to a pipe-separated
    #+ string for use in grep pattern matching for retained chromosomes and
    #+ program (@PG) lines
    # shellcheck disable=SC2001,SC2028
    pattern="^@HD|^@SQ.*SN:($(echo "${chromosomes}" | tr ' ' '|'))|^@PG"

    #  Define variables for subsequent filtering and reheading
    outdir="$(dirname "${outfile}")"
    outbam="$(basename "${outfile}")"
    bam_rh_init="${outdir}/rehead.${outbam}"
    bam_rh_sort="${outdir}/txt_rh_sort.${outbam}"
    
    #  Filter BAM infile, writing to BAM outfile
    # shellcheck disable=SC2086
    if ! \
        samtools view \
            -@ "${threads}" \
            -b \
            -o "${outfile}" \
            "${infile}" \
            ${chromosomes}
    then
        echo "Error: Failed to filter '${infile}'." >&2
        return 1
    fi

    #  Index the BAM outfile    
    samtools index -@ "${threads}" "${outfile}"
    
    #  Reheader the filtered BAM outfile, including sorting and indexing the
    #+ file
    if \
        samtools reheader \
            -c "grep -E '${pattern}'" \
            "${outfile}" \
                > "${bam_rh_init}"
    then
        samtools sort \
            -@ "${threads}" \
            -o "${bam_rh_sort}" \
            "${bam_rh_init}"

        mv -f \
            "${bam_rh_sort}" \
            "${outfile}"

        rm "${bam_rh_init}"

        samtools index -@ "${threads}" "${outfile}"
    else
        echo "Error: Failed to reheader '${outfile}'." >&2
        return 1
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
Usage:
  filter_bam_sp
    [--help] [--threads <int>] --infile <str> --outfile <str> [--mito] [--tg] [--mtr] [--chk_chr]

Description:
  Filter and reheader a BAM file for S. pombe chromosomes.

Keyword arguments:
   -t, --threads  <int>  Number of threads to use (default: ${threads}).
   -i, --infile   <str>  Coordinate-sorted BAM infile.
   -o, --outfile  <str>  Filtered BAM outfile.
   -m, --mito     <flg>  Retain SP_Mito chromosome (optional).
  -tg, --tg       <flg>  Retain SP_II_TG chromosome (optional).
  -mr, --mtr      <flg>  Retain SP_MTR chromosome (optional).
  -cc, --chk_chr  <flg>  Check chromosomes in filtered BAM outfile (optional).

Returns:
  Creates a BAM outfile filtered and reheadered for S. pombe chromosomes at the specified path.

Dependencies:
  - Programs
    + awk
    + Bash or Zsh
    + grep
    + mv
    + rm
    + Samtools

Examples:
  '''bash
  #TODO
  '''

#TODO:
  - Somewhere and somehow, need to handle more than S. pombe as "spike-in" organism.
  - The filter activity needs to be recorded in the BAM header.
  - Support for SAM and CRAM too.
EOM
    )

    #  Parse keyword arguments
    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
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
            -cc|--chk_chr) chk_chr=true;   shift 1 ;;
            *)
                echo "## Unknown argument passed: '${1}' ##" >&2
                echo "" >&2
                echo "${show_help}" >&2
                exit 1
                ;;
        esac
    done

    #  Validate keyword arguments
    if [[ -z "${threads}" ]]; then
        echo "Error: '--threads' is required." >&2
        return 1
    fi

    if [[ ! "${threads}" =~ ^[1-9][0-9]*$ ]]; then
        echo \
            "Error: '--threads' was assigned '${threads}' but must be a" \
            "positive integer greater than or equal to 1." >&2
        return 1
    fi

    if [[ -z "${infile}" ]]; then
        echo "Error: '--infile' is a required argument." >&2
        return 1
    fi

    if [[ ! -f "${infile}" ]]; then
        echo \
            "Error: File associated with '--infile' does not exist:" \
            "${infile}." >&2
        return 1
    fi

    if [[ -z "${outfile}" ]]; then
        echo "Error: '--outfile' is a required argument." >&2
        return 1
    fi

    if [[ ! -d "$(dirname "${outfile}")" ]]; then
        echo \
            "Error: Directory associated with '--outfile' does not exist:" \
            "$(dirname "${outfile}")" >&2
        return 1
    fi

    #  Do the main work
    #TODO: initialize local variables first
    #  Set the chromosomes to retain
    chromosomes="SP_I SP_II SP_III"
    if ${tg};   then chromosomes="SP_II_TG ${chromosomes}"; fi
    if ${mtr};  then chromosomes="${chromosomes} SP_MTR";   fi
    if ${mito}; then chromosomes="${chromosomes} SP_Mito";  fi

    #  Define variables for subsequent filtering and reheading
    outdir="$(dirname "${outfile}")"
    insam="$(basename "${infile/.bam/.sam}")"
    outsam="$(basename "${outfile/.bam/.sam}")"
    pth_in="${outdir}/${insam}"
    pth_out="${outdir}/${outsam}"

    #  Filter BAM infile, writing to BAM outfile
    #  First, set a trap to remove intermediate files when and however the
    #+ function exits
    trap 'rm -f "${pth_in}" "${pth_out}"' EXIT

    #  Begin by converting BAM to SAM while including the current header in the
    #+ SAM
    # shellcheck disable=SC2086
    if ! \
        samtools view \
            -@ "${threads}" \
            -h \
            -o "${pth_in}" \
            "${infile}"
    then
        echo "Error: Failed to generate ${pth_in}." >&2
        return 1
    fi

    # #  Filter pertinent lines from the SAM header, and perform filtering for
    # #+ lines that meet the following conditions: (1) They do not begin with @
    # #+ (i.e., are not in the header) and (2) they contain strings assigned to
    # #+ variable 'chromosomes' within field 3
    # shellcheck disable=SC2002
    # if ! \
    #     cat "${pth_in}" \
    #         | awk \
    #             -v chromosomes="${chromosomes}" \
    #             'BEGIN {
    #                 split(chromosomes, chrom_arr, " ")
    #                 for (i in chrom_arr) {
    #                     chrom_map[chrom_arr[i]] = 1
    #                 }
    #             }
    #
    #             function print_header_line() {
    #                 if ($0 ~ /^@SQ/ && $0 ~ /SN:([^ \t]+)/) {
    #                     split($0, a, /SN:/)
    #                     split(a[2], b, /[ \t]/)
    #                     sn = b[1]
    #                     if (!(sn in chrom_map)) {
    #                         next
    #                     }
    #                 }
    #                 print $0;
    #             }
    #
    #             {
    #                 if ($0 ~ /^@/) {
    #                     print_header_line()
    #                 } else {
    #                     if ($3 in chrom_map) {
    #                         print $0
    #                     }
    #                 }
    #             }' \
    #                 > "${pth_out}"
    # then
    #     echo "Error: Failed to generate ${pth_out}." >&2
    #     return 1
    # fi
    #TODO: The above only works with GNU Bash but is (relatively) fast; need fast for BSD

    #  Filter pertinent lines from the SAM header and body following these
    #+ conditions:
    #+ (1) Headers ('@') are retained but '@SQ' lines are filtered based on
    #+     'SN:'
    #+ (2) Alignment lines are retained only if their reference sequence (col
    #+     3) is in the chromosome list
    if ! \
        awk -v chromosomes="${chromosomes}" '
            BEGIN {
                split(chromosomes, chrom_arr, " ")
                for (i in chrom_arr) {
                    chrom_map[chrom_arr[i]] = 1
                }
            }

            #  Process header lines (@SQ must match chrom_map)
            /^@SQ/ {
                if ($0 ~ /SN:/) {
                    split($0, a, /SN:/)
                    split(a[2], b, /[ \t]/)
                    sn = b[1]
                    if (sn in chrom_map) { print }
                }
                next
            }

            #  Print all other header lines (@HD, @PG, etc.)
            /^@/ { print; next }

            #  Print alignments only if column 3 (reference sequence) matches
            ($3 in chrom_map) { print }
        ' "${pth_in}" \
            > "${pth_out}"
    then
        echo "Error: Failed to generate ${pth_out}." >&2
        return 1
    fi
    #TODO: The above is portable (both GNU and BSD) but (relatively) slow

    #  Perform SAM-to-BAM conversion for filtered, reheadered alignment
    #+ information
    if ! \
        samtools view -b "${pth_out}" > "${outfile}"
    then
        echo "Error: Failed to generate ${outfile}." >&2
        return 1
    fi

    #  Index the BAM outfile
    samtools index -@ "${threads}" "${outfile}"

    #  Remove intermediate SAM file if present
    if [[ -f "${outfile}" ]]; then
        if [[ -f "${pth_in}" ]];  then rm "${pth_in}";  fi
        if [[ -f "${pth_out}" ]]; then rm "${pth_out}"; fi
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
