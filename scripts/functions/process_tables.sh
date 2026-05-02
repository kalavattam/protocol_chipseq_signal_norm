#!/usr/bin/env bash
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


# check_table
# check_table_column
# check_table_scaling_factor
# extract_field_str
# _validate_arg_csl
# _validate_args_table
# _parse_table_core
# parse_table
# parse_table_simple


#  Require Bash >= 4.4 before defining functions
if [[ -z "${BASH_VERSION:-}" ]]; then
    echo "error(shell):" \
        "this script must be sourced or run under Bash >= 4.4." >&2

    if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
        return 1
    else
        exit 1
    fi
elif ((
    BASH_VERSINFO[0] < 4 || ( BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] < 4 )
)); then
    echo "error($(basename "${BASH_SOURCE[0]}")):" \
        "this script requires Bash >= 4.4; current version is" \
        "'${BASH_VERSION}'." >&2

    if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
        return 1
    else
        exit 1
    fi
fi

#  Source required helper functions if needed
# shellcheck disable=SC1091
{
    _dir_src_tbls="$(
        cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
    )"

    source "${_dir_src_tbls}/source_helpers.sh" || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source '${_dir_src_tbls}/source_helpers.sh'." >&2

        if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
            return 1
        else
            exit 1
        fi
    }

    source_helpers "${_dir_src_tbls}" \
        check_inputs check_numbers check_source format_outputs || {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to source required helper dependencies." >&2

            if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
                return 1
            else
                exit 1
            fi
        }

    unset _dir_src_tbls
}


#TODO: audit current usage; keep this helper even if unused
function check_table() {
    local table="${1:-}"
    local show_help

    show_help=$(cat << EOM
Usage:
  check_table [-h|--hlp|--help] table

Description:
  Check that a table file exists, is readable, non-empty, and contains at least one non-blank line.

Positional argument:
  1  table  <str>  Path to the table file.

Returns:
  0 if the table is usable; otherwise 1.
EOM
    )

    if [[ "${table}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${table}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'table', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    validate_var_file "table" "${table}" || return 1

    if [[ ! -s "${table}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "table file is empty: '${table}'."
        return 1
    fi

    if ! \
        awk 'NF { found = 1; exit } END { exit(found ? 0 : 1) }' "${table}"
    then
        echo_err_func "${FUNCNAME[0]}" \
            "table file '${table}' contains no non-blank lines."
        return 1
    fi
}


#QUESTION: should this remain public?
function check_table_column() {
    local table="${1:-}"
    local column="${2:-}"
    local header
    local show_help

    show_help=$(cat << EOM
Usage:
  check_table_column [-h|--hlp|--help] table column

Description:
  Check that a specified column name is present in the header row of a tab-delimited table.

Positional arguments:
  1  table   <str>  Path to the table file.
  2  column  <str>  Column name to look for in the header row.

Returns:
  0 if the column is present; otherwise 1.
EOM
    )

    if [[ "${table}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${table}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'table', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${column}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'column', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    validate_var_file "table" "${table}" || return 1

    header="$(awk 'NR == 1 { print; exit }' "${table}")"

    if ! \
        awk -F '\t' -v col="${column}" '
            NR == 1 {
                for (i = 1; i <= NF; i++) {
                    if ($i == col) exit 0
                }
                exit 1
            }
        ' "${table}"
    then
        echo_err_func "${FUNCNAME[0]}" \
            "column '${column}' was not found in table header:" \
            "'$(printf '%s' "${header}" | tr '\t' ' ')'."
        return 1
    fi
}


#QUESTION: should this remain public?
function check_table_scaling_factor() {
    local type="${1:-}"
    local table="${2:-}"
    local scl_fct="${3:-}"
    local name="${4:-}"
    local type_lc scl_fct_lc
    local msg=""
    local show_help

    show_help=$(cat << EOM
Usage:
  check_table_scaling_factor [-h|--hlp|--help] type table scl_fct name

Description:
  Check table-related scaling-factor or Boolean-like input and, when applicable, print a note describing how command-line values interact with table-derived values.

Positional arguments:
  1  type     <str>  Expected value type: 'str', 'string', 'bol', 'bool', 'boolean', 'flg', or 'flag'.
  2  table    <str>  Path to the table file.
  3  scl_fct  <str>  Scaling-factor string or Boolean-like value ('true', 't', 'false', 'f') to check.
  4  name     <str>  Option name being checked, such as 'scl_fct' or 'typ_cvg'.

Returns:
  0 if the input is valid; otherwise 1.

Output:
  Prints a note to stdout when applicable.
EOM
    )

    if [[ "${type}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${type}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'type', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${table}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'table', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${scl_fct}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 3, 'scl_fct', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${name}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 4, 'name', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    type_lc="${type,,}"
    scl_fct_lc="${scl_fct,,}"

    case "${type_lc}" in
        str|string)
            type_lc="str"
            ;;
        bol|bool|boolean|flg|flag)
            type_lc="bol"
            ;;
        *)
            echo_err_func "${FUNCNAME[0]}" \
                "positional argument 1, 'type', is '${type}' but must be" \
                "'str', 'string', 'bol', 'bool', 'boolean', 'flg', or 'flag'."
            return 1
            ;;
    esac

    validate_var_file "table" "${table}" || return 1

    if [[ "${type_lc}" == "str" && "${scl_fct_lc}" == "na" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 3, 'scl_fct', must not be 'NA' when" \
            "positional argument 1, 'type', resolves to 'string'."
        return 1
    fi

    if [[ "${type_lc}" == "bol" ]]; then
        case "${scl_fct_lc}" in
            true|t|false|f) : ;;
            *)
                echo_err_func "${FUNCNAME[0]}" \
                    "positional argument 3, 'scl_fct', must be 'true', 't'," \
                    "'false', or 'f' when positional argument 1, 'type'" \
                    "resolves to 'boolean': '${scl_fct}'."
                return 1
                ;;
        esac
    fi

    #NOTE: string, instead of array, accumulation is fine here
    if [[ "${name}" == "typ_cvg" ]]; then
        msg="Note: '--${name}' scaling factors will be multiplied with those"
        msg+=" from '--table' ('--tbl_col') if present; otherwise, '--${name}'"
        msg+=" scaling factors will be applied directly to raw coverage."
    elif [[ "${name}" == "scl_fct" ]]; then
        msg="Note: '--${name}' will override scaling factors from '--table'"
        msg+=" ('--tbl_col')."
    fi

    if [[ -n "${msg}" ]]; then
        echo "${msg}"
    fi
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

  This function validates inputs, ensuring file exists, is readable, has more than just a header, and is properly formatted as a TSV file.

  The function also checks that the specified field index is within the valid range of columns.

Positional arguments:
  1  tbl  <str>  Path to the TSV file.
  2  fld  <int>  The 1-based index of the column to extract.
  3  hdr  <bol>  'true' to skip header, 'false' to include it (default: ${hdr}).

Returns:
  Prints a comma-separated list (<str>) containing the values from the specified column and returns 0; otherwise returns 1.

Examples:
  '''bash
  tbl="example_table.tsv"
  extract_field_str "\${tbl}" 2
  extract_field_str "\${tbl}" 2 true
  extract_field_str "\${tbl}" 10  # Returns 1 and prints an error message to stderr because column 10 does not exist
  '''
EOM
    )

    #  Parse and check function arguments
    if [[ "${tbl}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${tbl}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'tbl', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${fld}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'fld', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  Check that 'tbl' exists, is readable, and is not empty
    validate_var_file "tbl" "${tbl}" || return 1

    if [[ ! -s "${tbl}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "table file is empty: '${tbl}'."
        return 1
    fi

    #  Validate 'fld' is a positive integer
    check_int_pos "${fld}" "fld" || return 1

    #  Convert 'hdr' to lowercase letters so Boolean matching accepts values
    #+ like 'T' / 'tRuE' / 'False' / 'FALSE' / etc
    hdr_lc="${hdr,,}"

    #  Check that 'hdr' is a properly formatted Boolean
    case "${hdr_lc}" in
        t|true)  hdr=true  ;;
        f|false) hdr=false ;;
        *)
            echo_err_func "${FUNCNAME[0]}" \
                "positional argument 3, 'hdr', is '${hdr}' but must be" \
                "'t', 'true', 'f', or 'false'."
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
        echo_err_func "${FUNCNAME[0]}" \
            "table file does not contain enough data rows: '${tbl}'."
        return 1
    fi

    #  Check that 'tbl' is tab-separated by checking the first data line
    if ! \
        awk \
            -F '\t' -v hdr="${hdr}" \
            'NR == (hdr == "true" ? 2 : 1) { exit (NF > 1 ? 0 : 1) }' \
            "${tbl}"
    then
        echo_err_func "${FUNCNAME[0]}" \
            "table file does not appear to be tab-separated: '${tbl}'."
        return 1
    fi

    #  Determine the number of fields in 'tbl' from the first data line
    num_fld="$(
        awk \
            -F '\t' -v hdr="${hdr}" \
            'NR == (hdr == "true" ? 2 : 1) { print NF; exit }' \
            "${tbl}"
    )"

    #  Check that 'fld' is within the valid range of fields
    if (( fld > num_fld )); then
        echo_err_func "${FUNCNAME[0]}" \
            "field index '${fld}' is out of range; table '${tbl}' has" \
            "${num_fld} fields."
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


#MAYBE: move into 'check_inputs.sh'?
#MAYBE: 'validate' should be 'check'?
function _validate_arg_csl() {
    local value="${1:-}"
    local valid="${2:-}"
    local name="${3:-}"
    local value_lc
    local opt opt_lc
    local found=false
    local -a arr_valid
    local show_help

    show_help=$(cat << EOM
Usage:
  _validate_arg_csl [-h|--hlp|--help] value valid name

Description:
  Validate that a value matches one member of a comma-delimited set of allowed values.

Positional arguments:
  1  value  <str>  Value to check.
  2  valid  <str>  Comma-delimited string of allowed values.
  3  name   <str>  Argument/option name associated with 'value'.

Returns:
  0 if 'value' matches one of the allowed values; otherwise 1.

Note:
  - 'csl' stands for "comma-separated list".

Examples:
  '''bash
  _validate_arg_csl sf "sf,scl,scaled,alpha,alf,spike,spk,spike-in,si,siq,siq_chip,siq-chip" "tbl_col"
  _validate_arg_csl TRUE "true,false,t,f" "norm"
  '''
EOM
    )

    if [[ "${value}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${value}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'value', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${valid}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'valid', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${name}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 3, 'name', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    value_lc="${value,,}"
    IFS=',' read -r -a arr_valid <<< "${valid}"

    for opt in "${arr_valid[@]}"; do
        opt_lc="${opt,,}"

        if [[ "${value_lc}" == "${opt_lc}" ]]; then
            found=true
            break
        fi
    done

    if [[ "${found}" != "true" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "invalid value '${value}' for '${name}'. Valid options are:" \
            "${valid}."
        return 1
    fi
}


function _validate_args_table() {
    local table="${1:-}"
    local tbl_col="${2:-}"
    local norm="${3:-}"
    local raw="${4:-}"
    local show_help

    show_help=$(cat << EOM
Usage:
  _validate_args_table [-h|--hlp|--help] table tbl_col norm raw

Description:
  Validate the common positional arguments used by the 'parse_table' family of functions.

Positional arguments:
  1  table    <str>  Path to a TSV table file.
  2  tbl_col  <str>  Table column name for scaling factors; must belong to one of these families:
                       - generic ('sf', 'scl', 'scaled', 'alpha', 'alf'),
                       - spike-in ('spike', 'spk', 'spike-in', 'si'), or
                       - siQ-ChIP ('siq', 'siq_chip', 'siq-chip').
  3  norm     <str>  Normalized-coverage flag; accepts Boolean-like 'true', 'false', 't', or 'f' in any letter case.
  4  raw      <str>  Raw/unadjusted-coverage flag; accepts Boolean-like 'true', 'false', 't', or 'f' in any letter case.

Returns:
  0 if all parameters are valid; otherwise 1.

Notes:
  - This helper checks only the shared parser inputs.
  - Column-presence checks for specific table headers are handled downstream in 'parse_table' functions.
  - 'spike_in' is reserved for the spike-in-organism input-count column and is not accepted as a scaling-factor column name.
  - Additional aliases such as 'mip', 'min', 'sip', 'sin', 'vip', 'vin', 'msp', 'msn', 'dip', 'din', 'lip', and 'lin' are handled downstream during table-header parsing.
EOM
    )

    if [[ "${table}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${table}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'table', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${tbl_col}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'tbl_col', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${norm}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 3, 'norm', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${raw}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 4, 'raw', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    validate_var_file "table" "${table}" || return 1
    check_table "${table}" || return 1

    _validate_arg_csl \
        "${tbl_col}" \
        "sf,scl,scaled,alpha,alf,spike,spk,spike-in,si,siq,siq_chip,siq-chip" \
        "tbl_col" \
        || return 1
    _validate_arg_csl "${norm}" "true,false,t,f" "norm" || return 1
    _validate_arg_csl "${raw}"  "true,false,t,f" "raw"  || return 1
}


function _parse_table_core() {
    local mode="${1:-}"
    local table="${2:-}"
    local tbl_col="${3:-}"
    local norm="${4:-}"
    local raw="${5:-}"
    local norm_lc raw_lc header
    local i
    local show_help

    local idx_smp=-1
    local idx_scl=-1
    local idx_spike=-1
    local idx_siq=-1
    local idx_mp=-1   idx_mn=-1   idx_sp=-1   idx_sn=-1
    local idx_msp=-1  idx_msn=-1  idx_vp=-1   idx_vn=-1
    local idx_dp=-1   idx_dn=-1   idx_lp=-1   idx_ln=-1

    local n_scl=0
    local n_spike=0
    local n_siq=0
    local col_lc

    local -a arr_header
    local -a arr_infiles  arr_stm_out
    local -a arr_scl      arr_spike     arr_siq
    local -a arr_main_ip  arr_main_in   arr_spike_ip arr_spike_in
    local -a arr_mass_ip  arr_mass_in   arr_volume_ip arr_volume_in
    local -a arr_depth_ip arr_depth_in  arr_length_ip arr_length_in

    local n_col fld
    local -a fields

    show_help=$(cat << EOM
Usage:
  _parse_table_core [-h|--hlp|--help] mode table tbl_col norm raw

Description:
  Internal helper for parsing tab-delimited tables in either 'simple' or 'complex' mode.

Positional arguments:
  1  mode     <str>  Parsing mode: 'simple' or 'complex'.
  2  table    <str>  Path to a TSV table file.
  3  tbl_col  <str>  Table column name for scaling factors; must belong to one of these families:
                       - generic ('sf', 'scl', 'scaled', 'alpha', 'alf'),
                       - spike-in ('spike', 'spk', 'spike-in', 'si'), or
                       - siQ-ChIP ('siq', 'siq_chip', 'siq-chip').
  4  norm     <bol>  Normalized-coverage flag; accepts Boolean-like 'true', 'false', 't', or 'f' in any letter case.
  5  raw      <bol>  Raw/unadjusted-coverage flag; accepts Boolean-like 'true', 'false', 't', or 'f' in any letter case.

Returns:
  Prints 'declare -p' declarations for parsed arrays to stdout; otherwise 1.
EOM
    )

    if [[ "${mode}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${mode}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'mode', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${table}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'table', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${tbl_col}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 3, 'tbl_col', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${norm}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 4, 'norm', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${raw}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 5, 'raw', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    case "${mode}" in
        simple|complex) : ;;
        *)
            echo_err_func "${FUNCNAME[0]}" \
                "positional argument 1, 'mode', must be 'simple' or 'complex':" \
                "'${mode}'."
            return 1
            ;;
    esac

    _validate_args_table \
        "${table}" "${tbl_col}" "${norm}" "${raw}" \
        || return 1

    if ! command -v gawk > /dev/null 2>&1; then
        echo_err_func "${FUNCNAME[0]}" \
            "'gawk' is required for robust TSV parsing in '${table}'."
        return 1
    fi

    tbl_col="${tbl_col,,}"
    norm_lc="${norm,,}"
    raw_lc="${raw,,}"

    case "${norm_lc}" in
        t|true)  norm=true  ;;
        f|false) norm=false ;;
    esac

    case "${raw_lc}" in
        t|true)  raw=true  ;;
        f|false) raw=false ;;
    esac

    header="$(awk 'NR == 1 { print; exit }' "${table}")"
    IFS=$'\t' read -r -a arr_header <<< "${header}"
    #NOTE: header parsing still uses Bash splitting; move this to 'gawk' later
    #+     only if trailing-empty header fields ever need to be preserved

    n_col="${#arr_header[@]}"

    for i in "${!arr_header[@]}"; do
        col_lc="${arr_header[i],,}"

        case "${col_lc}" in
            sample)               idx_smp=${i} ;;
            sf|scl|scaled|alpha|alf)
                idx_scl=${i}
                n_scl=$(( n_scl + 1 ))
                ;;
            spike|spk|spike-in|si)
                idx_spike=${i}
                n_spike=$(( n_spike + 1 ))
                ;;
            siq|siq_chip|siq-chip)
                idx_siq=${i}
                n_siq=$(( n_siq + 1 ))
                ;;
            main_ip|mip)          idx_mp=${i}  ;;
            main_in|min)          idx_mn=${i}  ;;
            spike_ip|spk_ip|sip)  idx_sp=${i}  ;;
            spike_in|spk_in|sin)  idx_sn=${i}  ;;
            mass_ip|msp)          idx_msp=${i} ;;
            mass_in|msn)          idx_msn=${i} ;;
            volume_ip|vol_ip|vip) idx_vp=${i}  ;;
            volume_in|vol_in|vin) idx_vn=${i}  ;;
            depth_ip|dep_ip|dip)  idx_dp=${i}  ;;
            depth_in|dep_in|din)  idx_dn=${i}  ;;
            length_ip|len_ip|lip) idx_lp=${i}  ;;
            length_in|len_in|lin) idx_ln=${i}  ;;
        esac
    done

    if (( n_scl > 1 )); then
        echo_err_func "${FUNCNAME[0]}" \
            "table '${table}' contains more than one generic scaling-factor" \
            "column synonym (e.g., 'sf', 'scl', 'scaled', 'alpha', 'alf')."
        return 1
    fi

    if (( n_spike > 1 )); then
        echo_err_func "${FUNCNAME[0]}" \
            "table '${table}' contains more than one spike-in scaling-factor" \
            "column synonym (e.g., 'spike', 'spk', 'spike-in', 'si')."
        return 1
    fi

    if (( n_siq > 1 )); then
        echo_err_func "${FUNCNAME[0]}" \
            "table '${table}' contains more than one siQ-ChIP scaling-factor" \
            "column synonym (e.g., 'siq', 'siq_chip', 'siq-chip')."
        return 1
    fi

    if (( idx_smp == -1 )); then
        echo_err_func "${FUNCNAME[0]}" \
            "required column 'sample' was not found in table '${table}'."
        return 1
    fi

    case "${tbl_col}" in
        sf|scl|scaled|alpha|alf) tbl_col="scl"   ;;
        spike|spk|spike-in|si)   tbl_col="spike" ;;
        siq|siq_chip|siq-chip)   tbl_col="siq"   ;;
    esac

    case "${tbl_col}" in
        scl)
            if (( idx_scl == -1 )); then
                echo_err_func "${FUNCNAME[0]}" \
                    "requested generic scaling-factor column was not found in" \
                    "'${table}'."
                return 1
            fi
            ;;
        spike)
            if (( idx_spike == -1 )); then
                echo_err_func "${FUNCNAME[0]}" \
                    "requested spike-in scaling-factor column was not found" \
                    "in '${table}'."
                return 1
            fi
            ;;
        siq)
            if (( idx_siq == -1 )); then
                echo_err_func "${FUNCNAME[0]}" \
                    "requested siQ-ChIP scaling-factor column was not found" \
                    "in '${table}'."
                return 1
            fi
            ;;
    esac

    while :; do
        fields=()

        for ((i = 0; i < n_col; i++)); do
            if ! IFS= read -r -d '' fld; then
                if (( i == 0 )); then
                    break 2
                fi

                echo_err_func "${FUNCNAME[0]}" \
                    "encountered a truncated parsed row while reading table" \
                    "'${table}'."
                return 1
            fi

            fields+=( "${fld}" )
        done

        arr_infiles+=( "${fields[idx_smp]}" )

        if [[ "${mode}" == "complex" ]]; then
            arr_stm_out+=( "$(basename "${fields[idx_smp]}" .bam)" )
        fi

        [[ ${idx_scl}   -ne -1 ]] && arr_scl+=( "${fields[idx_scl]}" )
        [[ ${idx_spike} -ne -1 ]] && arr_spike+=( "${fields[idx_spike]}" )
        [[ ${idx_siq}   -ne -1 ]] && arr_siq+=( "${fields[idx_siq]}" )

        if [[ "${mode}" == "complex" ]]; then
            [[ ${idx_mp}   -ne -1 ]] && arr_main_ip+=( "${fields[idx_mp]}" )
            [[ ${idx_mn}   -ne -1 ]] && arr_main_in+=( "${fields[idx_mn]}" )
            [[ ${idx_sp}   -ne -1 ]] && arr_spike_ip+=( "${fields[idx_sp]}" )
            [[ ${idx_sn}   -ne -1 ]] && arr_spike_in+=( "${fields[idx_sn]}" )
            [[ ${idx_msp}  -ne -1 ]] && arr_mass_ip+=( "${fields[idx_msp]}" )
            [[ ${idx_msn}  -ne -1 ]] && arr_mass_in+=( "${fields[idx_msn]}" )
            [[ ${idx_vp}   -ne -1 ]] && arr_volume_ip+=( "${fields[idx_vp]}" )
            [[ ${idx_vn}   -ne -1 ]] && arr_volume_in+=( "${fields[idx_vn]}" )
            [[ ${idx_dp}   -ne -1 ]] && arr_depth_ip+=( "${fields[idx_dp]}" )
            [[ ${idx_dn}   -ne -1 ]] && arr_depth_in+=( "${fields[idx_dn]}" )
            [[ ${idx_lp}   -ne -1 ]] && arr_length_ip+=( "${fields[idx_lp]}" )
            [[ ${idx_ln}   -ne -1 ]] && arr_length_in+=( "${fields[idx_ln]}" )
        fi
    done < <(
        #  Use 'gawk' explicitly so TSV rows are parsed faithfully, including
        #+ any trailing empty fields, which can then be detected or handled
        #+ explicitly ('awk' is 'gawk' in 'env_protocol')
        gawk -v n_col="${n_col}" '
            NR == 1 { next }
            {
                n = split($0, fld, /\t/)

                for (i = 1; i <= n_col; i++) {
                    if (i <= n) {
                        printf "%s%c", fld[i], 0
                    } else {
                        printf "%c", 0
                    }
                }
            }
        ' "${table}"
    )

    unset arr_scl_fct && declare -a arr_scl_fct
    if [[ -z "${scl_fct:-}" ]] && ! ${norm} && ! ${raw}; then
        case "${tbl_col}" in
            scl)   arr_scl_fct=( "${arr_scl[@]}" )   ;;
            spike) arr_scl_fct=( "${arr_spike[@]}" ) ;;
            siq)   arr_scl_fct=( "${arr_siq[@]}" )   ;;
        esac
    fi

    if [[ "${mode}" == "complex" ]]; then
        declare -p arr_infiles  arr_stm_out
        declare -p arr_scl      arr_spike    arr_siq
        declare -p arr_main_ip  arr_main_in  arr_spike_ip  arr_spike_in
        declare -p arr_mass_ip  arr_mass_in  arr_volume_ip arr_volume_in
        declare -p arr_depth_ip arr_depth_in arr_length_ip arr_length_in
        declare -p arr_scl_fct
    else
        declare -p arr_infiles arr_scl_fct
    fi
}


function parse_table() {
    local table="${1:-}"
    local tbl_col="${2:-}"
    local norm="${3:-}"
    local raw="${4:-}"
    local show_help

    show_help=$(cat << EOM
Usage:
  parse_table [-h|--hlp|--help] table tbl_col norm raw

Description:
  Parse a tab-delimited table and extract supported columns into shell arrays.

Positional arguments:
  1  table    <str>  Path to a TSV table file.
  2  tbl_col  <str>  Table column name for scaling factors; must belong to one of these families:
                       - generic ('sf', 'scl', 'scaled', 'alpha', 'alf'),
                       - spike-in ('spike', 'spk', 'spike-in', 'si'), or
                       - siQ-ChIP ('siq', 'siq_chip', 'siq-chip').
  3  norm     <bol>  Normalized-coverage flag; accepts Boolean-like 'true', 'false', 't', or 'f' in any letter case.
  4  raw      <bol>  Raw/unadjusted-coverage flag; accepts Boolean-like 'true', 'false', 't', or 'f' in any letter case.

Returns:
  Prints 'declare -p' declarations for parsed arrays to stdout; otherwise 1.

Output arrays:
  arr_infiles, arr_stm_out
  arr_scl, arr_spike, arr_siq
  arr_main_ip, arr_main_in, arr_spike_ip, arr_spike_in
  arr_mass_ip, arr_mass_in, arr_volume_ip, arr_volume_in
  arr_depth_ip, arr_depth_in, arr_length_ip, arr_length_in
  arr_scl_fct

Recognized scaling-factor column names:
  - Generic scaling-factor family:
      'sf', 'scl', 'scaled', 'alpha', 'alf'  ->  arr_scl
  - Spike-in scaling-factor family:
      'spike', 'spk', 'spike-in', 'si'       ->  arr_spike
  - siQ-ChIP scaling-factor family:
      'siq', 'siq_chip', 'siq-chip'          ->  arr_siq

Recognized non-scaling-factor column aliases:
  - 'main_ip',   'mip'            ->  arr_main_ip
  - 'main_in',   'min'            ->  arr_main_in
  - 'spike_ip',  'spk_ip', 'sip'  ->  arr_spike_ip
  - 'spike_in',  'spk_in', 'sin'  ->  arr_spike_in
  - 'mass_ip',   'msp'            ->  arr_mass_ip
  - 'mass_in',   'msn'            ->  arr_mass_in
  - 'volume_ip', 'vol_ip', 'vip'  ->  arr_volume_ip
  - 'volume_in', 'vol_in', 'vin'  ->  arr_volume_in
  - 'depth_ip',  'dep_ip', 'dip'  ->  arr_depth_ip
  - 'depth_in',  'dep_in', 'din'  ->  arr_depth_in
  - 'length_ip', 'len_ip', 'lip'  ->  arr_length_ip
  - 'length_in', 'len_in', 'lin'  ->  arr_length_in

Notes:
  - The table must contain a 'sample' column.
  - The requested 'tbl_col' must also be present in the header.
  - This function reads global variable 'scl_fct' to decide whether table-derived scaling factors should populate 'arr_scl_fct'.
  - 'spike_in' is reserved for the spike-in-organism input-count column and is not accepted as a scaling-factor column name.

Examples:
  '''bash
  parse_table "samples.tsv" sf false false
  parse_table "samples.tsv" alpha true false
  '''
EOM
    )

    if [[ "${table}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${table}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'table', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${tbl_col}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'tbl_col', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${norm}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 3, 'norm', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${raw}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 4, 'raw', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    _parse_table_core "complex" "${table}" "${tbl_col}" "${norm}" "${raw}"
}


function parse_table_simple() {
    local table="${1:-}"
    local tbl_col="${2:-}"
    local norm="${3:-}"
    local raw="${4:-}"
    local show_help

    show_help=$(cat << EOM
Usage:
  parse_table_simple [-h|--hlp|--help] table tbl_col norm raw

Description:
  Parse a tab-delimited table and extract only the sample column together with one table-derived scaling-factor array.

Positional arguments:
  1  table    <str>  Path to a TSV table file.
  2  tbl_col  <str>  Table column name for scaling factors; must belong to one of these families:
                       - generic ('sf', 'scl', 'scaled', 'alpha', 'alf'),
                       - spike-in ('spike', 'spk', 'spike-in', 'si'), or
                       - siQ-ChIP ('siq', 'siq_chip', 'siq-chip').
  3  norm     <bol>  Normalized-coverage flag; accepts Boolean-like 'true', 'false', 't', or 'f' in any letter case.
  4  raw      <bol>  Raw/unadjusted-coverage flag; accepts Boolean-like 'true', 'false', 't', or 'f' in any letter case.

Returns:
  Prints 'declare -p' declarations for 'arr_infiles' and 'arr_scl_fct' to stdout; otherwise 1.

Output arrays:
  arr_infiles
  arr_scl_fct

Notes:
  - The table must contain a 'sample' column.
  - The requested 'tbl_col' must be present in the header.
  - This function reads global variable 'scl_fct' to decide whether table-derived scaling factors should populate 'arr_scl_fct'.
  - 'spike_in' is reserved for the spike-in-organism input-count column and is not accepted as a scaling-factor column name.

Examples:
  '''bash
  parse_table_simple "samples.tsv" sf false false
  parse_table_simple "samples.tsv" alpha true false
  '''
EOM
    )

    if [[ "${table}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${table}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'table', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${tbl_col}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 2, 'tbl_col', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${norm}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 3, 'norm', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    elif [[ -z "${raw}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 4, 'raw', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    _parse_table_core "simple" "${table}" "${tbl_col}" "${norm}" "${raw}"
}


#  Print an error message when function script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    err_source_only "${BASH_SOURCE[0]}"
fi
