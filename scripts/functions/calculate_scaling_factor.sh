#!/bin/bash

#  calculate_scaling_factor.sh
#  KA

#  Support functions for 'submit_calculate_scaling_factor.sh'


#  Parse siQ-ChIP metadata using Python script
function parse_metadata() {
    local scr_met="${1}"  # Python script for parsing metadata
    local tbl_met="${2}"  # siQ-ChIP metadata table
    local eqn="${3}"      # Equation to use for parsing
    local fil_bam="${4}"  # BAM file to process

    if ! \
        eval "$(
            python "${scr_met}" \
                --tbl_met "${tbl_met}" \
                --eqn "${eqn}" \
                --bam "${fil_bam}" \
                --shell
        )"
    then
        echo \
            "Error: Failed to use script '${scr_met}' to parse siQ-ChIP" \
            "metadata in '${tbl_met}' for file '${fil_bam}'." >&2
        return 1
    fi
}


#  Get expression for filtering flags based on alignment type
function get_expr_filter() {
    local fil_typ="${1:-pe}"  # Alignment type for file: "paired", "pe" for
                              # paired-end; "single", "se" for single-end

    case "${fil_typ}" in
        paired|pe) \
            echo "(flag == 99) || (flag == 1123) || (flag == 163) || (flag == 1187)"
        ;;
        single|se) \
            echo "(flag == 0) || (flag == 1024) || (flag == 16) || (flag == 1040)"
        ;;
    esac
}


#  Count the number of alignments in a BAM file
function count_align_bam() {
    local threads="${1}"      # Number of threads for parallelization
    local fil_in="${2}"       # Input (not IP) BAM file
    local fil_typ="${3:-pe}"  # "paired", "pe", "single", "se" (default: "pe")
    local expr                # Alignment filtering expression
    local show_help           # Help message/documentation

    show_help=$(cat << EOM
---------------
count_align_bam
---------------

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
  count_align_bam "\${threads}" "\${fil_in}" "\${fil_typ}"

Examples:
  \`\`\`
  #  Count alignments in a BAM file of paired-end alignments using 8 threads
  count_align_bam 8 sample.bam paired

  #  Count alignments in a BAM file of single-end alignments using 4 threads
  count_align_bam 4 sample.bam single
  \`\`\`
EOM
    )

    if [[ -z "${1:-}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
        echo "${show_help}"
        echo ""
        if [[ -z "${1:-}" ]]; then return 1; else return 0; fi
    fi

    if [[ ! "${threads}" =~ ^[1-9][0-9]*$ ]]; then
        echo \
            "Error: Positional parameter 1, 'threads', must be a positive" \
            "integer: '${threads}'." >&2
        return 1
    fi

    if [[ ! -f "${fil_in}" ]]; then
        echo \
            "Error: Positional parameter 2, 'fil_in', not found:" \
            "'${fil_in}'." >&2
        return 1
    fi

    case "${fil_typ}" in
        single|se|paired|pe) : ;;
        *)
            echo \
                "Error: Positional parameter 3, 'fil_typ', must be 'paired'" \
                "or 'single': '${fil_typ}'." >&2
            return 1
        ;;
    esac

    #  Determine filtering flags based on alignment type
    expr=$(get_expr_filter "${fil_typ}")

    #  Count alignments based on alignment type
    samtools view -@ "${threads}" -c --expr "${expr}" "${fil_in}"
}


#  Computes the average fragment length for a BAM file comprised of paired- or
#+ single-end alignments
function calculate_frag_avg() {
    local threads="${1}"      # Number of threads for parallelization
    local fil="${2}"          # Input BAM file
    local fil_typ="${3:-pe}"  # "paired", "pe", "single", "se" (default: "pe")
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

    if [[ -z "${1:-}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
        echo "${show_help}"
        echo ""
        if [[ -z "${1:-}" ]]; then return 1; else return 0; fi
    fi

    if [[ ! "${threads}" =~ ^[1-9][0-9]*$ ]]; then
        echo \
            "Error: Positional parameter 1, 'threads', must be a positive" \
            "integer: '${threads}'." >&2
        return 1
    fi

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
    expr=$(get_expr_filter "${fil_typ}")

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


#  Compute scaling factor (alpha or spike-in) using the appropriate Python
#+ script
function compute_scl_fct() {
    local mode="${1}"     # Scaling factor mode: "alpha" or "spike"
    local scr_alf="${2}"  # Script to compute alpha scaling factor
    local scr_spk="${3}"  # Script to compute spike-in scaling factor
    shift 3
    local params=( "$@" )
    
    case "${mode}" in
        alpha) python "${scr_alf}" "${params[@]}" ;;
        spike) python "${scr_spk}" "${params[@]}" ;;
        *)
            echo \
                "Error: Positional parameter 1, 'mode', is '${mode}' but" \
                "must be either 'alpha' or 'spike'." >&2
            return 1
            ;;
    esac
}


#  Compute minimum input depth factors based on the number of alignments in an
#+ input BAM file, bin size, and effective genome size
function calculate_dep_fct() {
    local n_in="${1}"      # Alignment count for input (not IP) BAM file
    local siz_bin="${2}"   # Bin size (in bp)
    local siz_gen="${3}"   # Effective genome size for model organism (in bp)
    local mode="${4}"      # "frag" or "norm"
    local rnd="${5:-24}"   # Number of decimal points for rounding
    local fct_dep          # Variable for calculations
    local show_help        # Help message/documentation

    show_help=$(cat << EOM
-----------------
calculate_dep_fct
-----------------

Description:
  Computes a minimum input depth factor based on the number of alignments in an
  input BAM file, bin size, and effective genome size. The calculation supports
  both fragment-length normalization ("frag") and "normalized coverage"
  ("norm"; see PMID: 37160995 for more details).

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
      "\${n_in}" "\${siz_bin}" "\${siz_gen}" "\${mode}" "\${rnd}"

Examples:
  \`\`\`
  #  Compute depth factor for fragment-length normalized coverage
  calculate_factor_depth 12851824 20 12157105 "frag" 12

  #  Compute depth factor for "normalized coverage" 
  calculate_factor_depth 12851824 30 12157105 "norm" 24
  \`\`\`
EOM
    )

    if [[ -z "${1:-}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
        echo "${show_help}"
        echo ""
        if [[ -z "${1:-}" ]]; then return 1; else return 0; fi
    fi

    for var in n_in siz_bin siz_gen rnd; do
        if [[ ! "${!var}" =~ ^[1-9][0-9]*$ ]]; then
            echo "Error: '${var}' must be a positive integer: '${!var}'." >&2
            return 1
        fi
    done

    case "${mode}" in
        frag|norm) : ;;
        *)
            echo \
                "Error: Positional parameter 4, 'mode', must be 'frag' or" \
                "'norm': '${mode}'." >&2
            return 1
        ;;
    esac

    #  Compute depth factor via 'bc' operation based on 'mode'
    if [[ "${mode}" == "norm" ]]; then
        #  For "normalized coverage"
        fct_dep=$(bc -l <<< "
            scale=${rnd}; 
            (${siz_bin}) / (${siz_gen} * (1 - (${siz_bin} / ${siz_gen})))
        ")
    else
        #  For fragment-length adjusted coverage
        fct_dep=$(bc -l <<< "
            scale=${rnd}; 
            (${n_in} * ${siz_bin}) / (${siz_gen} * (1 - (${siz_bin} / ${siz_gen})))
        ")
    fi

    echo "${fct_dep}"
}


#  Compute fragment depths for fragment length-normalized signal ("frag") and
#+ "normalized coverage" per Dickson et al., Sci Rep 2023 ("norm")
function calculate_dep_arr() {
    local dep="${1}"   # Number of mapped reads in the sample
    local rnd="${2}"   # Rounding precision for output
    local mode="${3}"  # Data transformation mode: "frag" or "norm"
    local -a arr_dep   # Initialize empty array for depth values
    
    for var in dep rnd; do
        if ! [[ "${!var}" =~ ^[1-9][0-9]*$ ]]; then
            echo "Error: '${var}' is '${!var}' but must be a positive integer." >&2
            return 1
        fi
    done

    case "${mode}" in
        frag|norm) : ;;
        *)
            echo \
                "Error: 'mode' is '${mode}' but must be either 'frag' or" \
                "'norm'." >&2
            return 1
            ;;
    esac

    for bin in 1 5 10 20 30 40 50; do
        arr_dep+=( "$(
            calculate_dep_fct \
                "${dep}" "${bin}" 12157105 "${mode}" "${rnd}"
        )" )
    done

    IFS=',' printf "%s\n" "${arr_dep[*]}"
}


#  Compute minimum input depth values for fragment-adjusted signal and
#+ normalized coverage
function compute_dep_all() {
    local dep="${1}"  # Number of alignments in sample BAM
    local rnd="${2}"  # Number of decimals to round to
    local -a arr_dm_fr arr_dm_nm

    # shellcheck disable=SC2207
    {
        arr_dm_fr=(
            $(calculate_dep_arr "${dep}" "${rnd}" "frag")
        ) || return 1
        arr_dm_nm=(
            $(calculate_dep_arr "${dep}" "${rnd}" "norm")
        ) || return 1
    }

    output=$(IFS=','; echo "${arr_dm_fr[*]},${arr_dm_nm[*]}")
    echo "${output}"
}


#  Dynamically construct 'printf' format string
function generate_fmt_str() {
    local num_fld="${1}"  # Number of fields in the output row
    local fmt_str=""      # Variable for format string
    
    for ((i = 1; i <= num_fld; i++)); do fmt_str+="%s\t"; done
    printf "%s\n" "${fmt_str%$'\t'}"
}


#  Compute siQ-ChIP alpha scaling factor and related values for a sample
#+ 
#+ Workflow function that processes a sample using global variables; extracts
#+ siQ-ChIP metadata from a TSV table, computes alignment counts and average
#+ fragment lengths, and minimum input depth values
#+ #TODO: Delineate variables expected to be assigned and available globally
# shellcheck disable=SC2154
function process_samp_alpha() {
    #TODO: Write, return 'show_help' when called without argument or with '-h'
    local idx="${1}"  # Array sample index

    #  Declare local variables
    local fil_ip fil_in alpha mass_ip mass_in vol_all vol_in
    local dep_ip dep_in len_ip len_in arr_dm fmt_str
    
    #  Assign BAM files based on sample index
    fil_ip="${arr_mip[idx]}"
    fil_in="${arr_min[idx]}"

    #  Check that BAM files exist
    validate_var_file "fil_ip" "${fil_ip}" "${idx}" || return 1
    validate_var_file "fil_in" "${fil_in}" "${idx}" || return 1

    if ${debug:-false}; then
        debug_var "idx=${idx}" "fil_ip=${fil_ip}" "fil_in=${fil_in}"
    fi

    #  Parse siQ-ChIP metadata, assigning global variables
    parse_metadata "${scr_met}" "${tbl_met}" "${eqn}" "${fil_ip}" || return 1

    #TODO: Add and use function to determine alignment type of BAM ('fil_typ')
    
    #  Count alignments in BAM files
    dep_ip=$(count_align_bam "${threads}" "${fil_ip}") || return 1
    dep_in=$(count_align_bam "${threads}" "${fil_in}") || return 1

    #  Compute average fragment lengths for BAM files
    len_ip=$(calculate_frag_avg "${threads}" "${fil_ip}") || return 1
    len_in=$(calculate_frag_avg "${threads}" "${fil_in}") || return 1

    if ${debug:-false}; then
        debug_var \
            "eqn=${eqn}"         "rnd=${rnd}" \
            "mass_ip=${mass_ip}" "mass_in=${mass_in}" \
            "vol_all=${vol_all}" "vol_in=${vol_in}" \
            "dep_ip=${dep_ip}"   "dep_in=${dep_in}" \
            "len_ip=${len_ip}"   "len_in=${len_in}"
    fi

    #  Compute siQ-ChIP alpha scaling factor
    alpha=$(
        compute_scl_fct \
             "alpha" "${scr_alf}" "${scr_spk}" \
            --eqn     "${eqn}"     --rnd     "${rnd}" \
            --mass_ip "${mass_ip}" --mass_in "${mass_in}" \
            --vol_all "${vol_all}" --vol_in  "${vol_in}" \
            --dep_ip  "${dep_ip}"  --dep_in  "${dep_in}" \
            --len_ip  "${len_ip}"  --len_in  "${len_in}"
    ) || return 1

    #  Compute minimum inputh depth values
    IFS="," read -r -a arr_dm < <(
        compute_dep_all "${dep_in}" "${rnd}"
    ) || return 1

    #  Dynamically construct string for output formatting
    fmt_str=$(generate_fmt_str 13)

    #  Print results using dynamically generated format string
    # shellcheck disable=SC2059
    {
        IFS=$' \t\n'  # Ensure deafult setting for IFS
        printf "${fmt_str}\n" \
            "${fil_ip}" "${fil_in}" "${alpha}" "${eqn}" \
            "${mass_ip}" "${mass_in}" "${vol_all}" "${vol_in}" \
            "${dep_ip}" "${dep_in}" "${len_ip}" "${len_in}" \
            "${arr_dm[*]}" \
        | sed 's:\t$::' \
            >> "${fil_out}"
    }
}


#  Compute spike-in scaling factor and related values for a sample
#+ 
#+ Workflow function that processes a sample using global variables; computes
#+ computes alignment counts and minimum input depth values
#+ #TODO: Delineate variables expected to be assigned and available globally
# shellcheck disable=SC2154
function process_samp_spike() {
    #TODO: Write, return 'show_help' when called without argument or with '-h'
    local idx="${1}"  # Array sample index

    #  Declare local variables
    local mp mn sp sn num_mp num_sp num_mn num_sn sf arr_dm fmt_str

    #  Assign BAM files based on sample index
    mp="${arr_mip[idx]}"
    mn="${arr_min[idx]}"
    sp="${arr_sip[idx]}"
    sn="${arr_sin[idx]}"

    #  Check that BAM files exist
    validate_var_file "mp" "${mp}" "${idx}" || return 1
    validate_var_file "sp" "${sp}" "${idx}" || return 1
    validate_var_file "mn" "${mn}" "${idx}" || return 1
    validate_var_file "sn" "${sn}" "${idx}" || return 1

    if ${debug:-false}; then
        debug_var "idx=${idx}" "mp=${mp}" "mn=${mn}" "sp=${sp}" "sn=${sn}"
    fi

    #TODO: Add and use function to determine alignment type of BAM ('fil_typ')

    #  Count alignments in BAM files
    num_mp=$(count_align_bam "${threads}" "${mp}") || return 1
    num_sp=$(count_align_bam "${threads}" "${sp}") || return 1
    num_mn=$(count_align_bam "${threads}" "${mn}") || return 1
    num_sn=$(count_align_bam "${threads}" "${sn}") || return 1

    if ${debug:-false}; then
        debug_var \
            "num_mp=${num_mp}" "num_sp=${num_sp}" \
            "num_mn=${num_mn}" "num_sn=${num_sn}" \
            "rnd=${rnd}"
    fi

    #  Compute spike-in scaling factor
    sf=$(
        compute_scl_fct \
            "spike" "${scr_alf}" "${scr_spk}" \
            --main_ip "${num_mp}" --spike_ip "${num_sp}" \
            --main_in "${num_mn}" --spike_in "${num_sn}" \
            --rnd     "${rnd}"
    ) || return 1

    #  Compute minimum input depth values
    IFS="," read -r -a arr_dm < <(
        compute_dep_all "${num_mn}" "${rnd}"
    ) || return 1

    #  Dynamically construct string for output formatting
    fmt_str=$(generate_fmt_str 10)

    #  Print results using dynamically generated format string
    # shellcheck disable=SC2059
    {
        IFS=$' \t\n'  # Ensure default setting for IFS
        printf "${fmt_str}\n" \
            "${mp}" "${sp}" "${mn}" "${sn}" "${sf}" \
            "${num_mp}" "${num_sp}" "${num_mn}" "${num_sn}" \
            "${arr_dm[*]}" \
        | sed 's:\t$::' \
            >> "${fil_out}"
    }
}
