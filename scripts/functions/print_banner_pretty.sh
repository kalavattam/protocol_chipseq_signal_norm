#!/bin/bash

function print_banner_pretty() {
    local text=""    # Text to wrap
    local wrap="#"   # Wrapping character
    local pad=1      # Spaces between side markers and text
    local cols=77    # Max text width per line; 0 means no wrapping
    local show_help  # Help message

    show_help=$(cat << EOM
Usage:
  print_banner_pretty [<empty>|-h|--hlp|--help] -t|--text <str> [-w|--wrap <str>] [-p|--pad <int>] [-c|--cols <int>]

Description:
  Pretty-print a single-line banner around user-provided text.

Keyword arguments:
  -t, --text  <str>  Text to wrap (required). If omitted, any trailing positional arguments will be joined with single spaces and used as text.
  -w, --wrap  <str>  Wrapping character; only the first character is used (default: '${wrap}').
  -p, --pad   <int>  Number of spaces between the side markers and the text on each side (default: '${pad}').
  -c, --cols  <int>  Maximum text width per line before wrapping. If 0, the text is not wrapped and is printed on a single line (default: '${cols}').

Notes:
  - '--wrap' uses only the first character of its value.
  - For a single line, banner width is as follows: len(text) + (2 * pad) + 4 (for '##' + spaces + text + spaces + '##').

Examples:
  Example 1:
  '''bash
  print_banner_pretty
      --text "This is the text that will be 'wrapped'."
      --wrap "#"
  '''

  '''txt
  ##############################################
  ## This is the text that will be 'wrapped'. ##
  ##############################################
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

#TODO:
  - The parser does not check for missing values after options '--text', '--wrap', '--pad', and '--cols'.
    + For example, 'print_banner_pretty --pad' will silently assign an empty value and later fail.
    + This needs to be fixed by adding (i) empty-call help behavior and (ii) stricter missing-argument checks.
EOM
    )

    #  Parse and check function arguments
    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    while (( $# )); do
        case "$1" in
            -t|--text) text="${2:-}";  shift 2      ;;
            -w|--wrap) wrap="${2:-#}"; shift 2      ;;
            -p|--pad)  pad="${2:-1}";  shift 2      ;;
            -c|--cols) cols="${2:-0}"; shift 2      ;;
            --)                        shift; break ;;
            *)
                #  Treat any leftover args as part of the text
                if [[ -z "${text}" ]]; then
                    text="${1}"
                else
                    text="${text} ${1}"
                fi
                shift
                ;;
        esac
    done

    #  Require text
    if [[ -z "${text}" ]]; then
        echo \
            "Error: '--text' is required (or supply text as trailing" \
            "args). For more information, run 'print_banner_pretty -h'." >&2
        return 1
    fi

    #  Check 'pad'
    if ! [[ ${pad} =~ ^[0-9]+$ ]]; then
        echo "Error: '--pad' must be a non-negative integer, got '${pad}'." >&2
        return 1
    fi

    #  Check 'cols' (if set)
    if ! [[ ${cols} =~ ^[0-9]+$ ]]; then
        echo \
            "Error: '--cols' must be a non-negative integer, got" \
            "'${cols}'." >&2
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
    #      2 (side) + pad + len_max + pad + 2 (side)
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
