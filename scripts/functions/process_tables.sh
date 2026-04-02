#!/bin/bash
# -*- coding: utf-8 -*-
# 
# Script: process_tables.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
# 
# Distributed under the MIT license.


#  Validate TSV table infile is not empty
function check_table() {
    local table="${1}"

    if [[ $(wc -l < "${table}") -le 1 ]]; then
        echo \
            "Error: Table file '${table}' is empty or contains only a" \
            "header." >&2
        return 1
    fi
}


#  Validate existence of column in table
function check_table_column() {
    local table="${1}"
    local column="${2}"

    if \
        ! awk -F '\t' -v col="${column}" '
            NR == 1 {
                for (i = 1; i <= NF; i++) {
                    if ($i == col) exit 0
                } exit 1
            }
        ' "${table}"
    then
        echo \
            "Error: Column '${column}' not found in table header:" \
            "$(awk 'NR == 1' "${table}" | tr '\t' ' ')"
        return 1
    fi
}


#  Note that scaling factor(s) or normalization will be multiplied with those
#+ in the table
function check_table_scaling_factor() {
    local type="${1,,}"   # Expected values: 'string' or 'boolean', etc.
    local table="${2}"    # Path to table file
    local scl_fct="${3}"  # Scaling factor(s) (string) or norm. flag (boolean)
    local name="${4}"     # Option name being checked (e.g., 'scl_fct', 'typ_cvg')
    local msg

    #  Validate positional parameters
    if [[ -z "${type}" ]]; then
        echo "Error: Positional parameter 1, 'type', is required." >&2
        return 1
    fi

    if [[ ! "${type}" =~ ^(str|string|bol|bool|boolean|flg|flag)$ ]]; then
        echo \
            "Error: Invalid positional parameter 1, 'type': '${type}'." \
            "Expected 'str', 'string', 'bol', 'bool' or 'boolean'." >&2
        return 1
    fi

    if [[ -z "${table}" ]]; then
        echo "Error: Positional parameter 2, 'table', is required." >&2
        return 1
    fi

    if [[ ! -f "${table}" ]]; then
        echo \
            "Error: Positional parameter 2, 'table', does not exist:" \
            "'${table}'." >&2
        return 1
    fi

    if [[ -z "${scl_fct}" ]]; then
        echo "Error: Positional parameter 3, 'scl_fct', is required." >&2
        return 1
    fi

    if [[ "${type}" == "string" && "${scl_fct}" == "NA" ]]; then
        echo \
            "Error: Positional parameter 3, 'scl_fct', is assigned 'NA'," \
            "which is invalid for type (positional parameter 1) 'string'." >&2
        return 1
    fi  #TODO: Not sure if this is needed...

    if [[
        "${type}" == "boolean" && "${scl_fct}" != true && "${scl_fct}" != false
    ]]; then
        echo \
            "Error: Positional parameter 3, 'scl_fct', must be 'true' or" \
            "'false' for type (positional parameter 1) 'boolean'." >&2
        return 1
    fi

    if [[ -z "${name}" ]]; then
        echo "Error: Positional parameter 4, 'name', is required." >&2
        return 1
    fi

    #  Determine appropriate message
    if [[ "${name}" == "typ_cvg" ]]; then
        msg="Note: --${name} scaling factors will be multiplied with those"
        msg+=" from --table (--tbl_col) if present; otherwise, --${name}"
        msg+=" scaling factors will be applied directly to raw coverage."
    elif [[ "${name}" == "scl_fct" && -n "${scl_fct}" ]]; then
        msg="Note: --${name} will override scaling factors from --table"
        msg+=" (--tbl_col)."
    fi
    
    #  Output the message if applicable
    if [[ -n "${msg}" ]]; then echo "${msg}"; fi
}


function extract_field_str() {
    local tbl="${1:-}"       # Path to the TSV file
    local fld="${2:-}"       # 1-based index of the column to extract
    local hdr="${3:-false}"  # Skip header (true/false)
    local hdr_lc             # Lowercase-converted skip-header Boolean
    local num_fld            # No. tab-delimited fields in first data line
    local show_help          # Help message

    show_help=$(cat << EOM
Usage:
  extract_field_str [-h|--hlp|--help] tbl fld [hdr]

Description:
  Extract a specific field (column) from a tab-separated value (TSV) file (e.g., a 'table' or 'tbl') and return it as a single, comma-separated string.

  Validates inputs, ensuring file exists, is readable, has more than just a header, and is properly formatted as a TSV file.

  Checks that the specified field index is within the valid range of columns.

Positional arguments:
  1  tbl  <str>  Path to the TSV file.
  2  fld  <int>  The 1-based index of the column to extract.
  3  hdr  <bol>  'true' to skip header, 'false' to include it (default: ${hdr}).

Returns:
  0 and a comma-separated string containing the values from the specified column; otherwise, 1 with an error message.

Examples:
  1. Extract the second column from a valid table
  '''bash
  tbl="example_table.tsv"
  extract_field_str "\${tbl}" 2
  '''

  '''txt
  "value1,value2,value3"
  '''

  2. Confirm invalid field index errors
  '''bash
  extract_field_str "\${tbl}" 0
  '''

  '''txt
  Error: Field index must be a positive integer.
  '''

  3. Confirm nonexistent file errors
  '''bash
  extract_field_str "nonexistent_file.tsv" 1
  '''

  '''txt
  Error: Table file does not exist: 'nonexistent_file.tsv'.
  '''

  4. Confirm out-of-range field index errors
  '''bash
  tbl="example_table.tsv"
  extract_field_str "\${tbl}" 10
  '''

  '''txt
  Error: Field index is out of range: '10'. The table has 3 fields.
  '''
EOM
    )

    #  Parse and check function arguments
    if [[ "${tbl}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${tbl}" ]]; then
        echo "Error: Positional argument 1, 'tbl', is missing." >&2
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    if [[ -z "${fld}" ]]; then
        echo "Error: Positional argument 2, 'fld', is missing." >&2
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  Check that 'tbl' exists, is readable, and is not empty
    if [[ ! -e "${tbl}" ]]; then
        echo "Error: Table file does not exist: '${tbl}'." >&2
        return 1
    elif [[ ! -r "${tbl}" ]]; then
        echo "Error: Table file is not readable: '${tbl}'." >&2
        return 1
    elif [[ ! -s "${tbl}" ]]; then
        echo "Error: Table file is empty: '${tbl}'." >&2
        return 1
    fi

    #  Validate 'fld' is a positive integer
    if [[ ! "${fld}" =~ ^[1-9][0-9]*$ ]]; then
        echo "Error: Field index must be a positive integer." >&2
        return 1
    fi

    #  Convert 'hdr' to lowercase letters so Boolean matching accepts values
    #+ like 'T' / 'tRuE' / 'False' / 'FALSE' / etc.
    hdr_lc=$(printf '%s' "${hdr}" | tr '[:upper:]' '[:lower:]')

    #  Check that 'hdr' is a properly formatted Boolean
    case "${hdr_lc}" in
        t|true)  hdr=true  ;;
        f|false) hdr=false ;;
        *)
            echo \
                "Error: Input for 'hdr', positional argument 3, is" \
                "'${hdr}'. Expected 't', 'true', 'f', or 'false'." >&2
            return 1
            ;;
    esac

    #  Check that the table has enough data rows: 1 or more lines, or header
    #+ and one or more lines
    if ! \
        awk \
            -v hdr="${hdr}" \
            'END {
                min_rows = (hdr == "true" ? 2 : 1)
                if (NR < min_rows) exit 1
            }' \
                "${tbl}"
    then
        echo \
            "Error: Table file does not contain enough data rows:" \
            "'${tbl}'." >&2
        return 1
    fi

    #  Check that 'tbl' is tab-separated by checking the first data line
    if ! \
        awk \
            -F '\t' \
            -v hdr="${hdr}" \
            'NR == (hdr == "true" ? 2 : 1) { exit (NF > 1 ? 0 : 1) }' \
                "${tbl}"
    then
        echo \
            "Error: Table file does not appear to be tab-separated:" \
            "'${tbl}'." >&2
        return 1
    fi

    #  Determine the number of fields in 'tbl' from the first data line
    num_fld="$(
        awk \
            -F '\t' \
            -v hdr="${hdr}" \
            'NR == (hdr == "true" ? 2 : 1) { print NF; exit }' \
                "${tbl}"
    )"

    #  Check that 'fld' is within the valid range of fields
    if (( fld > num_fld )); then
        echo \
            "Error: Field index is out of range: '${fld}'. The table has" \
            "${num_fld} fields." >&2
        return 1
    fi

    #  Extract the specified field and return it as a comma-separated string
    awk \
        -F '\t' \
        -v fld="${fld}" \
        -v hdr="${hdr}" \
        '(hdr == "true" && NR == 1) {
            next
        } {
            val = val ? val "," $fld : $fld
        } END {
            print val
        }' \
            "${tbl}"
}


#  Parse a tab-delimited table and dynamically extract values based on columns.
#+
#+ Positional parameters:
#+   1  table    Path to a TSV table file               <str>
#+   2  tbl_col  Table column name for scaling factors  <str:alpha,scaled,sf>
#+   3  norm     Normalized coverage flag               <bol:true,false>
#+   4  raw      Raw/unadjusted coverage flag           <bol:true,false>
#+
#+ Returns:
#+   Typeset declarations of arrays populated (or not) with parsed values from
#+   the table.
function parse_table() {
    local table="${1}"    # str: path to a TSV table file
    local tbl_col="${2}"  # str: column name to use for scaling factors
    local norm=${3}       # bol: normalized coverage flag
    local raw=${4}        # bol: raw/unadjusted coverage flag

    #  Validate 'parse_table' parameters
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


#  Parse a TSV table to extract specific columns: sample, scaling factors.
#+
#+ Positional parameters:
#+   1  table    Path to a TSV table file               <str>.
#+   2  tbl_col  Table column name for scaling factors  <str:alpha,scaled,sf>
#+   3  norm     Normalized coverage flag               <bol:true,false>
#+   4  raw      Raw/unadjusted coverage flag           <bol:true,false>
#+
#+ Returns:
#+   - 'arr_infiles': Array of sample values.
#+   - 'arr_scl_fct': Array of scaling factors for the specified 'tbl_col'.
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
        [[ ${idx_sf}  -ne -1 ]] && arr_sf+=(     "${fields[idx_sf]}"  )
        [[ ${idx_scl} -ne -1 ]] && arr_scaled+=( "${fields[idx_scl]}" )
        [[ ${idx_alf} -ne -1 ]] && arr_alpha+=(  "${fields[idx_alf]}" )
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


#  Helper function: validate individual 'parse_table' function parameters
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


#  Helper function: validate all 'parse_table' function parameters
function validate_table_params() {
    if [[ ! -f "${table}" ]]; then
        echo "Error: Table file '${table}' does not exist." >&2
        return 1
    fi

    validate_param "${tbl_col}" "alpha scaled sf" "--tbl_col" || return 1
    validate_param "${norm}"    "true false"      "--norm"    || return 1
    validate_param "${raw}"     "true false"      "--raw"     || return 1
}
