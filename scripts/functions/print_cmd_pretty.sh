#!/bin/bash

function print_cmd_pretty() {
    local first_n=2       # Number of tokens on first line
    local pair_flg=false  # Boolean: print flag token plus next token
    local indent=""       # Extra indenation
    local cont="    "     # Extra indentation: continuation lines
    local show_help       # Help message

    show_help=$(cat << EOM
Usage:
  print_cmd_pretty [<empty>|-h|--hlp|--help] [-hd|--head <int>|--head=<int>] [-pf|--pair_flg|--pair-flg] [-it|--indent] cmd arg1 arg2 ...

Description:
  Pretty-print a command in a readable, multi-line form.

Options:
  -hd <int> / -hd=<int>, --head <int> / --head=<int>
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
  print_cmd_pretty --pair-flg --indent --head=1 \\
      sbatch \\
          --account tsukiyama \\
          --time 12:00:00 \\
          --cpus-per-task 4 \\
              ./submit_compute_signal.sh \\
              --mode signal \\
              --method norm
  '''
EOM
    )

    #  Parse and check function arguments
    if [[ -z "${1:-}" || "${1}" =~ ^(-h|--h[e]?lp)$ ]]; then
        echo "${show_help}" >&2
        return 0
    fi

    while (( $# )); do
        case "${1}" in
            -hd|--head|--head=*)
                local val
                if [[ ${1} =~ ^(-hd|--head)$ ]]; then
                    if (( $# < 2 )); then
                        echo "Error: '--head' requires an argument." >&2
                        return 1
                    fi
                    val=${2}
                    shift 2
                else
                    if [[ "${1}" == -hd=* ]]; then
                        val=${1#-hd=}
                    else
                        val=${1#--head=}
                    fi
                    shift
                fi
                first_n="${val}"
                ;;

            -pf|--pair[_-]flg)
                pair_flg=true
                shift
                ;;

            -it|--indent)
                indent="  "
                cont="      "
                shift
                ;;

            --)
                shift
                break
                ;;

            *)
                #  First non-option is the command
                break
                ;;
        esac
    done

    #  Check that 'first_n' is a non-negative integer
    if ! [[ ${first_n} =~ ^[0-9]+$ ]]; then
        echo \
            "Error: 'first_n' ('--head <int>') must be a non-negative integer," \
            "but got '${first_n}'." >&2
        return 1
    elif (( first_n < 1 )); then
        #  Handle 'first_n=0'
        first_n=1
    fi

    #  Handle remaining args, which are the actual command
    local cmd=( "$@" )
    local n=${#cmd[@]}
    if (( n == 0 )); then return 0; fi

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
            if \
                   [[ ${pair_flg} == "true" ]] \
                && [[ ${cmd[start_i]} == -* ]] \
                && (( start_i + 1 < n )) \
                && [[ ${cmd[start_i+1]} != -* ]]
            then
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
