#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Script: format_outputs.sh
#
# Copyright 2024-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-4- and GPT-5-series models) was used in development.
#
# Distributed under the MIT license.


# echo_err
# echo_err_func
# echo_warn
# echo_warn_func
# format_print_cmd
# print_banner_pretty
# print_cmd_array
# print_cmd_pretty
# summarize_sig_norm


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
    _dir_src_fmt="$(
        cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
    )"

    source "${_dir_src_fmt}/source_helpers.sh" || {
        echo "error($(basename "${BASH_SOURCE[0]}")):" \
            "failed to source '${_dir_src_fmt}/source_helpers.sh'." >&2

        if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
            return 1
        else
            exit 1
        fi
    }

    source_helpers "${_dir_src_fmt}" \
        check_args check_source || {
            echo "error($(basename "${BASH_SOURCE[0]}")):" \
                "failed to source required helper dependencies." >&2

            if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
                return 1
            else
                exit 1
            fi
        }

    unset _dir_src_fmt
}


#TODO: remove the following functions after they are refactored away in scripts
# #  Write an error message to stderr and return code 1
# function echo_error() { echo "Error: $*" >&2; return 1; }
#
#
# #  Write a warning message to stderr and return code 0
# function echo_warning() { echo "Warning: $*" >&2; }


#  Write a script-level error message to stderr
function echo_err() {
    echo "error($(basename "${0}")): $*" >&2
}


#  Write a function-level error message to stderr
function echo_err_func() {
    local func="${1:-${FUNCNAME[1]:-main}}"
    shift || true
    echo "error($(basename "${0}")::${func}): $*" >&2
}


#  Write a script-level warning message to stderr
function echo_warn() {
    echo "warning($(basename "${0}")): $*" >&2
}


#  Write a function-level warning message to stderr
function echo_warn_func() {
    local func="${1:-${FUNCNAME[1]:-main}}"
    shift || true
    echo "warning($(basename "${0}")::${func}): $*" >&2
}


function format_print_cmd() {
    local slurm="${1:-}"  # Boolean-like: Slurm 'sbatch' command or not
    local scr="${2:-}"    # Script path for Slurm 'sbatch' command
    local slurm_lc
    local show_help       # Help message

    show_help=$(cat << EOM
Usage:
  format_print_cmd [-h|--hlp|--help] slurm [scr]

Description:
  Format and pretty-print a single command line for display, with indentation and line breaks.

  This helper is intended to be used with function 'build_cmd' by piping its output into this function.

Positional arguments:
  1  slurm  <bol>  Boolean-like string: "true", "t", "false", or "f":
                     - "true" or "t": treat the input as an 'sbatch' command and format sbatch flags and script arguments separately.
                     - "false" or "f": treat the input as a plain script call.
  2  scr    <str>  Script path as it appears in the 'sbatch' command. Required when 'slurm=true', otherwise ignored.

Examples:
  '''bash
  #  Plain command
  build_cmd ... | format_print_cmd false

  #  Slurm sbatch command invoking 'execute_compute_signal.sh'
  build_cmd ...
      | format_print_cmd true "\${dir_rep}/scripts/execute_compute_signal.sh"
  '''

Notes:
  - This function is display-oriented only. It does not perform robust shell parsing and should not be used to reconstruct or re-execute commands.
  - For non-Slurm commands, the function does the following:
    + Reads a single-line command string from stdin.
    + Breaks before each '--' and indents the resulting lines.
  - For Slurm-wrapped commands ('sbatch ... <script> ...'), the function:
    + Splits the command into three visual sections:
      - The leading 'sbatch \' line.
      - Indented sbatch flags (e.g., '--cpus-per-task', '--time', etc.).
      - The script path plus its arguments, further indented.
    + Prints a backslash at the end of each non-final line.
EOM
    )

    #  Parse and check function arguments
    if [[ "${slurm}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${slurm}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'slurm', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    slurm_lc="${slurm,,}"

    if [[ "${slurm_lc}" =~ ^(false|f)$ ]]; then
        #  For non-Slurm commands, break before each '--' and indent
        sed -E 's| --| \\\n    --|g'  # Only as robust as incoming command str
    elif [[ "${slurm_lc}" =~ ^(true|t)$ ]]; then
        #  For Slurm commands, require a script path
        if [[ -z "${scr}" ]]; then
            echo_err_func "${FUNCNAME[0]}" \
                "positional argument 2, 'scr', is missing when 'slurm=true'" \
                "or 'slurm=t'."
            echo >&2
            echo "${show_help}" >&2
            return 1
        fi

        #  Perform display-oriented formatting; this logic expects the current
        #+ single-line command layout emitted by 'build_cmd'
        #+
        #+ Standard AWK features are in use here, so 'awk' (rather than 'gawk')
        #+ is used intentionally
        awk \
            -v scr="${scr}" '
            BEGIN {
                indent_1 = "    "      # Set indent for sbatch flags
                indent_2 = "        "  # Set indent for script arguments
                in_scr = 0             # Set to 1 after encountering script
            }

            {
                #  Split the input string into tokens based on " --"
                n = split($0, a, / --/)

                for (i = 1; i <= n; i++) {
                    token = a[i]

                    #  Handle first token: strip and print "sbatch" if present
                    if (i == 1 && token ~ /^sbatch[[:space:]]*/) {
                        sub(/^sbatch[[:space:]]*/, "", token)
                        print "sbatch \\"
                        if (length(token) > 0) {
                            print indent_1 "--" token " \\"
                        }
                        continue
                    }

                    #  If token matches script path, switch indentation context
                    # if (token ~ scr) {  # (avoid accidental regex semantics)
                    if (index(token, scr) > 0) {
                        print indent_1 scr " \\"
                        in_scr = 1
                        continue
                    }

                    #  Print remaining tokens with appropriate indentation
                    if (length(token) > 0) {
                        #  Use script indentation after scr is printed
                        indent = (in_scr ? indent_2 : indent_1)

                        #  Add trailing backslash unless last token
                        print indent "--" token (i == n ? "" : " \\")
                    }
                }
            }'
    else
        echo_err_func "${FUNCNAME[0]}" \
            "unknown 'slurm' value '${slurm}'; expected 'true' / 't' or" \
            "'false' / 'f'."
        return 1
    fi
}


function print_banner_pretty() {
    local text=""    # Text to wrap
    local wrap="#"   # Wrapping character
    local pad=1      # Spaces between side markers and text
    local cols=77    # Max text width per line; 0 means no wrapping
    local show_help  # Help message

    show_help=$(cat << EOM
Usage:
  print_banner_pretty
    [-h|--hlp|--help] -tx|--text <str> [-w|--wrap <str>] [--pad <int>] [-cw|--cols <int>]

Description:
  Pretty-print a single-line banner around user-provided text.

Keyword arguments:
  -tx, --text  <str>  Text to wrap. If omitted, trailing positional arguments are joined with single spaces and used as text.
   -w, --wrap  <str>  Wrapping character; only the first character is used (default: '${wrap}').
  -pd, --pad   <int>  Number of spaces between the side markers and the text on each side (default: '${pad}').
  -cw, --cols  <int>  Maximum text width per line before wrapping. If 0, the text is not wrapped and is printed on a single line (default: '${cols}').

Notes:
  - '--wrap' uses only the first character of its value.
  - For a single line, banner width is as follows: len(text) + (2 * pad) + 4 (for '##' + spaces + text + spaces + '##').

Examples:
  Example 1:
  '''bash
  print_banner_pretty
      --text "This is the text that will be wrapped."
      --wrap "#"
  '''

  '''txt
  ############################################
  ## This is the text that will be wrapped. ##
  ############################################
  '''

  Example 2:
  '''bash
  print_banner_pretty
      --wrap "%"
      --pad 2
      --cols 40
      --text
          "This is a somewhat longer sentence that will be wrapped across"
          "multiple lines."
  '''

  '''txt
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%  This is a somewhat longer sentence that  %%
  %%  will be wrapped across multiple lines.   %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  '''
EOM
    )

    #  Parse and check function arguments
    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    while (( $# > 0 )); do
        case "${1}" in
            -tx|--txt|--text)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                text="${2}"
                shift 2
                ;;

            -w|--wrp|--wrap)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                wrap="${2}"
                shift 2
                ;;

            -pd|--pad)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                pad="${2}"
                shift 2
                ;;

            -cw|--col|--cols)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                cols="${2}"
                shift 2
                ;;

            --)
                shift
                break
                ;;

            *)
                #  Treat any leftover args as part of the text
                if [[ -z "${text}" ]]; then
                    text="${1}"
                else
                    text+=" ${1}"
                fi
                shift
                ;;
        esac
    done

    #  Require text
    if [[ -z "${text}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'--text' is required (or supply text as trailing args). For" \
            "more information, run 'print_banner_pretty -h'."
        return 1
    fi

    #  Require wrap character
    if [[ -z "${wrap}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'--wrap' must not be empty."
        return 1
    fi

    #  Check 'pad'
    if ! [[ ${pad} =~ ^[0-9]+$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'--pad' must be a non-negative integer, but got '${pad}'."
        return 1
    fi

    #  Check 'cols' (if set)
    if ! [[ ${cols} =~ ^[0-9]+$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'--cols' must be a non-negative integer, but got '${cols}'."
        return 1
    fi

    #  Use only the first character of '--wrap'
    local ch="${wrap:0:1}"
    local side="${ch}${ch}"

    #  Split text into wrapped lines
    local -a lines

    if (( cols > 0 )); then
        #  Word-wrapping to 'cols' characters per line (for text only)
        local -a words=()
        local w

        #  Basic whitespace splitting; newlines are treated as spaces
        for w in ${text}; do
            words+=( "${w}" )
        done

        local cur=""
        for w in "${words[@]}"; do
            if [[ -z "${cur}" ]]; then
                cur="${w}"
            elif (( ${#cur} + 1 + ${#w} <= cols )); then
                cur+=" ${w}"
            else
                lines+=( "${cur}" )
                cur="${w}"
            fi
        done
        [[ -n "${cur}" ]] && lines+=( "${cur}" )
    else
        lines=( "${text}" )
    fi

    #  Determine maximum line length
    local len_max=0 line
    for line in "${lines[@]}"; do
        if (( ${#line} > len_max )); then
            len_max=${#line}
        fi
    done

    #  Compute full banner width:
    #+     2 (side) + pad + len_max + pad + 2 (side)
    local width=$(( len_max + 2 * pad + 4 ))

    #  Build border line
    local border="" i
    for (( i = 0; i < width; i++ )); do
        border+="${ch}"
    done

    #  Print banner
    echo "${border}"

    #  For each line, pad out to 'len_max' so the block is rectangular
    local spc_r pad_spc extra

    #  Precompute padding spaces (pad spaces)
    pad_spc=""
    for (( i = 0; i < pad; i++ )); do pad_spc+=" "; done

    for line in "${lines[@]}"; do
        local len_line=${#line}
        extra=$(( len_max - len_line ))
        spc_r=""
        for (( i = 0; i < extra; i++ )); do
            spc_r+=" "
        done
        #  side + pad + text + (extra) + pad + side
        printf '%s%s%s%s%s%s\n' \
            "${side}" \
            "${pad_spc}" \
            "${line}" \
            "${spc_r}" \
            "${pad_spc}" \
            "${side}"
    done

    echo "${border}"
}


function print_cmd_array() {
    local arr_nam="${1:-}"  # Name of command array to print
    local decl              # Output from 'declare -p' for array validation
    local show_help         # Help message

    show_help=$(cat << EOM
Usage:
  print_cmd_array [-h|--hlp|--help] arr_nam

Description:
  Pretty-print a command array as a shell-escaped single-line command.

Positional arguments:
  1  arr_nam  <str>  Name of the indexed array variable to print.

Returns:
  0 after printing the shell-escaped array contents to stdout; 1 if 'arr_nam' is missing, invalid, unset, or not an indexed array.

Examples:
  '''bash
  print_cmd_array "arr_cmd"
  '''
EOM
    )

    if [[ "${arr_nam}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${arr_nam}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'arr_nam', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    if ! [[ "${arr_nam}" =~ ^[a-zA-Z_][a-zA-Z0-9_]*$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "invalid array name '${arr_nam}'."
        return 1
    fi

    if ! decl="$(declare -p "${arr_nam}" 2> /dev/null)"; then
        echo_err_func "${FUNCNAME[0]}" \
            "array '${arr_nam}' is unset."
        return 1
    elif [[ "${decl}" != declare\ -a* ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'${arr_nam}' is not an indexed array."
        return 1
    fi

    local -n arr_cmd_ref="${arr_nam}"

    if (( ${#arr_cmd_ref[@]} == 0 )); then
        echo_err_func "${FUNCNAME[0]}" \
            "array '${arr_nam}' is empty."
        return 1
    fi

    printf '%q ' "${arr_cmd_ref[@]}"
    printf '\n'

    return 0
}


function print_cmd_pretty() {
    local first_n=2       # Number of tokens on first line
    local pair_flg=false  # Boolean-like: print flag token plus next token
    local indent=""       # Extra indentation
    local cont="    "     # Extra indentation: continuation lines
    local show_help       # Help message

    show_help=$(cat << EOM
Usage:
  print_cmd_pretty
    [-h|--hlp|--help] [-hd|--head <int>] [-pf|--pair_flg|--pair-flg] [-it|--indent]
    cmd arg1 arg2 ...

Description:
  Pretty-print a command in a readable, multi-line form.

Options:
  -hd, --head  <int>
    Number of tokens to keep on the first line (default: ${first_n}).

  -pf, --pair_flg, --pair-flg
    On continuation lines, if the first token starts with '-' or '--', print that token plus the next token on the same line.

  -it, --indent
    Use extra indentation suitable for, e.g., embedding inside a "Command:" block.

Examples:
  '''bash
  print_cmd_pretty --pair_flg bash ./foo.sh --a 1 --b 2
  '''

  '''bash
  print_cmd_pretty --pair-flg --indent --head=1
      sbatch
          --account tsukiyama
          --time 12:00:00
          --cpus-per-task 4
              ./submit_compute_signal.sh
              --mode signal
              --method norm
  '''
EOM
    )

    #  Parse and check function arguments
    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    while (( $# > 0 )); do
        case "${1}" in
            -hd|--hd|--head)
                require_optarg "${1}" "${2:-}" "${FUNCNAME[0]}" || {
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                }
                first_n="${2}"
                shift 2
                ;;

            -pf|--pair[_-]flg)
                pair_flg=true
                shift
                ;;

            -it|--indnt|--indent)
                indent="  "
                cont="      "
                shift
                ;;

            --)
                shift
                break
                ;;

            *)
                if [[ "${1}" == -* ]]; then
                    echo "## Unknown argument passed: '${1}' ##" >&2
                    echo >&2
                    echo "${show_help}" >&2
                    return 1
                fi

                #  First non-option is the command
                break
                ;;
        esac
    done

    #  Check that 'first_n' is a non-negative integer
    if ! [[ ${first_n} =~ ^[0-9]+$ ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "'first_n' ('--head <int>') must be a non-negative integer, but" \
            "got '${first_n}'."
        return 1
    elif (( first_n < 1 )); then
        #  Handle 'first_n=0'
        first_n=1
    fi

    #  Handle remaining args, which are the actual command
    local cmd=( "$@" )
    local n=${#cmd[@]}
    if (( n == 0 )); then
        echo_err_func "${FUNCNAME[0]}" \
            "a command must be supplied after any formatting options."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    # Clamp 'first_n' so that it never exceeds the number of tokens
    if (( first_n > n )); then first_n=${n}; fi

    #  First line: print up to 'first_n' tokens
    local i=0
    printf '%s%q' "${indent}" "${cmd[0]}"
    i=1
    while (( i < n && i < first_n )); do
        printf ' %q' "${cmd[i]}"
        (( i++ ))
    done

    #  Handle subsequent lines, including final line
    if (( i < n )); then
        printf ' \\\n'

        #  Handle continuation lines
        local start_i end_i k
        while (( i < n )); do
            start_i=${i}
            end_i=${i}

            #  Pair “flag and value” if requested and possible:
            #+     - current token starts with '-' or '--'
            #+     - there is a next token
            #+     - the next token does not start with '-'
            if (
                   [[ ${pair_flg} == "true" ]] \
                && [[ ${cmd[start_i]} == -* ]] \
                && (( start_i + 1 < n )) \
                && [[ ${cmd[start_i+1]} != -* ]]
            ); then
                end_i=$(( start_i + 1 ))
            fi

            #  Print tokens 'cmd[start_i..end_i]'
            printf '%s%q' "${cont}" "${cmd[start_i]}"
            for (( k = start_i + 1; k <= end_i; k++ )); do
                printf ' %q' "${cmd[k]}"
            done

            if (( end_i == n - 1 )); then
                #  Handle last line: no trailing backslash
                printf '\n'
            else
                printf ' \\\n'
            fi

            i=$(( end_i + 1 ))
        done
    else
        printf '\n'
    fi
}


#MAYBE: 'format_outputs.sh' may not be the best place for this function
#MAYBE: make this function "private", i.e., '_summarize_sig_norm'
function summarize_sig_norm() {
    local typ_sig="${1:-}"  # Type of signal computation
    local scl_fct="${2:-}"  # Scaling factor
    local typ_sig_lc        # Lowercase-converted signal type for case matching
    local mth_nrm           # Normalization method message
    local src_scl           # Scaling-factor message
    local show_help         # Help message

    show_help=$(cat << EOM
Usage:
  summarize_sig_norm [-h|--hlp|--help] typ_sig [scl_fct]

Description:
  Summarize resolved signal-normalization and scaling states for the signal-computation workflow.

  Prints a short human-readable summary indicating (i) what normalization mode is implied by 'typ_sig' and (ii) whether multiplicative scaling factors were supplied.

Positional arguments:
  1  typ_sig  <str>  Type of signal computation; e.g., 'unadj', 'frag', or 'norm' (aliases accepted).
  2  scl_fct  <str>  Scaling factor string (optional). If empty, assumes no explicit '--scl_fct' was supplied.

Returns:
  0 after printing the summary to stdout; otherwise 1 if positional argument 1, 'typ_sig', is missing.

Examples:
  1. Summarize normalized coverage with no scaling factors
  '''bash
  summarize_sig_norm "norm"
  '''

  2. Summarize fragment-length normalization with explicit scaling
  '''bash
  summarize_sig_norm "frag" "1.25"
  '''
EOM
    )

    #  Parse and check function arguments
    if [[ "${typ_sig}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    elif [[ -z "${typ_sig}" ]]; then
        echo_err_func "${FUNCNAME[0]}" \
            "positional argument 1, 'typ_sig', is missing."
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  Lowercase-convert signal-type input so case matching is case-insensitive
    typ_sig_lc="${typ_sig,,}"

    #  Determine normalization method message
    case "${typ_sig_lc}" in
        u|unadj|unadjusted|s|smp|simple|r|raw)
            mth_nrm="No normalization; returning unadjusted signal:"
            mth_nrm+=" '--method ${typ_sig}'."
            ;;

        f|frg|frag|frg[_-]len|frag[_-]len|l|len|len[_-]frg|len[_-]frag)
            mth_nrm="Performing fragment-length normalization:"
            mth_nrm+=" '--method ${typ_sig}'."
            ;;

        n|nrm|norm|normalized)
            mth_nrm="Generating normalized coverage (Dickson et al., Sci Rep"
            mth_nrm+=" 2023): '--method ${typ_sig}'."
            ;;

        *)
            mth_nrm="Unknown normalization method: '--method ${typ_sig}'."
            ;;
    esac

    #  Determine scaling factor message
    if [[ -n "${scl_fct}" ]]; then
        src_scl="Custom multiplicative scaling factor(s):"
        src_scl+=" '--scl_fct ${scl_fct}'."
    else
        src_scl="No multiplicative scaling factor(s)."
    fi

    #  Print resolved argument states
    echo "#################################################"
    echo "## Summary of signal normalization and scaling ##"
    echo "#################################################"
    echo
    echo "- Normalization method: ${mth_nrm}"
    echo "- Scaling factor source: ${src_scl}"
    echo
    echo
}


#  Print an error message when function script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    err_source_only "${BASH_SOURCE[0]}"
fi
