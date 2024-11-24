#!/bin/bash


#  Helper function to validate individual parse_table function parameters
function validate_param() {
    local value="${1}"
    local valid="${2}"
    local name="${3}"

    if ! [[ " ${valid} " =~ ${value} ]]; then
        echo \
            "Error: Invalid value '${value}' for '${name}'. Valid options" \
            "are: ${valid}." >&2
        return 1
    fi
}


#  Helper function to validate all parse_table function parameters
function validate_table_params() {
    if [[ ! -f "${table}" ]]; then
        echo "Error: Table file '${table}' does not exist." >&2
        return 1
    fi

    validate_param "${tbl_col}" "alpha scaled sf" "--tbl_col" || return 1
    validate_param "${norm}"    "true false"      "--norm"    || return 1
    validate_param "${raw}"     "true false"      "--raw"     || return 1
}


#  Function to parse a tab-delimited table and dynamically extract values based
#+ on columns.
#+
#+ Positional parameters:
#+   $1: Path to a TSV table file (str).
#+   $2: Table column name for scaling factors (str: 'alpha', 'scaled', 'sf').
#+   $3: Normalized coverage flag (bol: 'true' or 'false').
#+   $4: Raw/unadjusted coverage flag (bol: 'true' or 'false').
#+
#+ Returns:
#+   Typeset declarations of arrays populated (or not) with parsed values from
#+   the table.
function parse_table() {
    local table="${1}"    # str: path to a TSV table file
    local tbl_col="${2}"  # str: column name to use for scaling factors
    local norm=${3}       # bol: normalized coverage flag
    local raw=${4}        # bol: raw/unadjusted coverage flag

    # Validate parse_table parameters
    validate_table_params || return 1

    #  Read the table's header to determine available columns
    local header
    header=$(awk 'NR == 1' "${table}")
    IFS=$'\t' read -r -a arr_header <<< "${header}"

    #  Initialize column indices
    local idx_smp=-1
    local idx_sf=-1  idx_scl=-1
    local idx_map=-1 idx_man=-1 idx_sp=-1  idx_sn=-1
    local idx_alf=-1
    local idx_mp=-1  idx_mn=-1  idx_vp=-1  idx_vn=-1
    local idx_dp=-1  idx_dn=-1  idx_lp=-1  idx_ln=-1

    #  Map column names to indices
    for i in "${!arr_header[@]}"; do
        case "${arr_header[i]}" in
            sample)    idx_smp=${i} ;;
            sf)        idx_sf=${i}  ;;
            scaled)    idx_scl=${i} ;;
            main_ip)   idx_map=${i} ;;
            main_in)   idx_man=${i} ;;
            spike_ip)  idx_sp=${i}  ;;
            spike_in)  idx_sn=${i}  ;;
            alpha)     idx_alf=${i} ;;
            mass_ip)   idx_mp=${i}  ;;
            mass_in)   idx_mn=${i}  ;;
            volume_ip) idx_vp=${i}  ;;
            volume_in) idx_vn=${i}  ;;
            depth_ip)  idx_dp=${i}  ;;
            depth_in)  idx_dn=${i}  ;;
            length_ip) idx_lp=${i}  ;;
            length_in) idx_ln=${i}  ;;
        esac
    done

    #  Initialize arrays to store parsed values
    local -a arr_infiles  arr_stm_out
    local -a arr_sf       arr_scaled
    local -a arr_main_ip  arr_main_in  arr_spike_ip  arr_spike_in
    local -a arr_alpha
    local -a arr_mass_ip  arr_mass_in  arr_volume_ip arr_volume_in
    local -a arr_depth_ip arr_depth_in arr_length_ip arr_length_in


    #  Parse the table rows, skipping the header
    while IFS=$'\t' read -r row; do
        IFS=$'\t' read -r -a fields <<< "${row}"

        #  Extract values if indices are valid, appending them to arrays
        [[ ${idx_sf}   -ne -1 ]] && arr_infiles+=( "${fields[idx_smp]}" )
        [[ ${idx_smp}  -ne -1 ]] \
            && arr_stm_out+=( "$(basename "${fields[idx_smp]}" .bam)" )
        [[ ${idx_sf}   -ne -1 ]] && arr_sf+=( "${fields[idx_sf]}" )
        [[ ${idx_scl}  -ne -1 ]] && arr_scaled+=( "${fields[idx_scl]}" )
        [[ ${idx_map}  -ne -1 ]] && arr_main_ip+=( "${fields[idx_map]}" )
        [[ ${idx_man}  -ne -1 ]] && arr_main_in+=( "${fields[idx_man]}" )
        [[ ${idx_sp}   -ne -1 ]] && arr_spike_ip+=( "${fields[idx_sp]}" )
        [[ ${idx_sn}   -ne -1 ]] && arr_spike_in+=( "${fields[idx_sn]}" )
        [[ ${idx_alf}  -ne -1 ]] && arr_alpha+=( "${fields[idx_alf]}" )
        [[ ${idx_mp}   -ne -1 ]] && arr_mass_ip+=( "${fields[idx_mp]}" )
        [[ ${idx_mn}   -ne -1 ]] && arr_mass_in+=( "${fields[idx_mn]}" )
        [[ ${idx_vp}   -ne -1 ]] && arr_volume_ip+=( "${fields[idx_vp]}" )
        [[ ${idx_vn}   -ne -1 ]] && arr_volume_in+=( "${fields[idx_vn]}" )
        [[ ${idx_dp}   -ne -1 ]] && arr_depth_ip+=( "${fields[idx_dp]}" )
        [[ ${idx_dn}   -ne -1 ]] && arr_depth_in+=( "${fields[idx_dn]}" )
        [[ ${idx_lp}   -ne -1 ]] && arr_length_ip+=( "${fields[idx_lp]}" )
        [[ ${idx_ln}   -ne -1 ]] && arr_length_in+=( "${fields[idx_ln]}" )
    done < <(awk 'NR > 1' "${table}")

    #  Dynamically assign scaling factor array
    unset arr_scl_fct && declare -a arr_scl_fct
    if [[ -z "${scl_fct}" ]] && ! ${norm} && ! ${raw}; then
        case "${tbl_col}" in
            sf)     arr_scl_fct=( "${arr_sf[@]}" )     ;;
            scaled) arr_scl_fct=( "${arr_scaled[@]}" ) ;;
            alpha)  arr_scl_fct=( "${arr_alpha[@]}" )  ;;
        esac
    fi

    #  Return the arrays, populated or not
    typeset -p arr_infiles  arr_stm_out
    typeset -p arr_sf       arr_scaled
    typeset -p arr_main_ip  arr_main_in  arr_spike_ip  arr_spike_in
    typeset -p arr_alpha
    typeset -p arr_mass_ip  arr_mass_in  arr_volume_ip arr_volume_in
    typeset -p arr_depth_ip arr_depth_in arr_length_ip arr_length_in
    typeset -p arr_scl_fct
}


#  Function to parse a TSV table to extract specific columns: sample, scaling
#+ factors.
#+
#+ Positional parameters:
#+   $1: Path to a TSV table file (str).
#+   $2: Table column name for scaling factors (str: 'alpha', 'scaled', 'sf').
#+   $3: Normalized coverage flag (bol: 'true' or 'false').
#+   $4: Raw/unadjusted coverage flag (bol: 'true' or 'false').
#+
#+ Returns:
#+   - arr_infiles: Array of sample values.
#+   - arr_scl_fct: Array of scaling factors for the specified tbl_col.
function parse_table_simple() {
    local table="${1}"    # str: path to a TSV table file
    local tbl_col="${2}"  # str: column name to use for scaling factors
    local norm=${3}       # bol: normalized coverage flag
    local raw=${4}        # bol: raw/unadjusted coverage flag

    # Validate parse_table parameters
    validate_table_params || return 1

    #  Read the table's header to determine available columns
    local header
    header=$(awk 'NR == 1' "${table}")
    IFS=$'\t' read -r -a arr_header <<< "${header}"

    #  Initialize column indices
    local idx_smp=-1 idx_sf=-1 idx_scl=-1 idx_alf=-1

    #  Map column names to indices
    for i in "${!arr_header[@]}"; do
        case "${arr_header[i]}" in
            sample) idx_smp=${i} ;;
            sf)     idx_sf=${i}  ;;
            scaled) idx_scl=${i} ;;
            alpha)  idx_alf=${i} ;;
        esac
    done

    #  Initialize arrays to store values
    local -a arr_infiles arr_sf arr_scaled arr_alpha

    #  Parse the table rows, skipping the header
    while IFS=$'\t' read -r row; do
        IFS=$'\t' read -r -a fields <<< "${row}"

        #  Extract the 'sample' column (always present)
        arr_infiles+=( "${fields[idx_smp]}" )

        #  Extract scaling factor columns if present
        [[ ${idx_sf}  -ne -1 ]] && arr_sf+=( "${fields[idx_sf]}" )
        [[ ${idx_scl} -ne -1 ]] && arr_scaled+=( "${fields[idx_scl]}" )
        [[ ${idx_alf} -ne -1 ]] && arr_alpha+=( "${fields[idx_alf]}" )
    done < <(awk 'NR > 1' "${table}")

    #  Dynamically assign scaling factor array
    unset arr_scl_fct && declare -a arr_scl_fct
    if [[ -z "${scl_fct}" ]] && ! ${norm} && ! ${raw}; then
        case "${tbl_col}" in
            sf)     arr_scl_fct=( "${arr_sf[@]}" )     ;;
            scaled) arr_scl_fct=( "${arr_scaled[@]}" ) ;;
            alpha)  arr_scl_fct=( "${arr_alpha[@]}" )  ;;
        esac
    fi

    #  Return the arrays
    typeset -p arr_infiles arr_scl_fct
}

