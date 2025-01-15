#!/bin/bash

function check_scl_fct {
    local scl_fct="${1}"  # Comma-separated scaling factors to validate
    local entries=()      # Array for individual comma-separated components
    local entry           # Individual scaling factor: 'num:den' or 'num'
    local num den         # Values for validation and formatting
    local val_fmt=()      # Array for validated and formatted entries
    local out_str         # Final formatted, comma-separated string output
    local show_help       # Help message/function documentation

    show_help=$(cat << EOM
-------------
check_scl_fct
-------------

Description:
  Function (subroutine) to validate and format positive floats for precomputed
  scaling factors used with deepTools 'bamCompare' or 'bigwigCompare'. This
  subroutine ensures the following:
    - Scaling factors are provided as a single comma-separated string.
    - Each factor pair (numerator and denominator) is in the format 'num:den'.
    - Single values (numerators) are automatically paired with a denominator
      of 1.
    - All values are validated as positive floats.
    - Leading decimal points (e.g., ".5") are zero-padded (e.g., "0.5").

Positional argument:
  1, scl_fct (str): Comma-separated scaling factors to validate.

Returns:
  0 and a validated, comma-separated string of formatted scaling factors;
  example: '0.4:1,1.2:1,0.65:0.8'. Otherwise, returns 1 with an error message.

Usage:
  check_scl_fct "\${scl_fct}"

Examples:
  \`\`\`
  scl_fct=0.5,1.2:3,.65
  check_scl_fct \${scl_fct}
  0.5:1,1.2:3,0.65:1

  scl_fct=0.00456,0.00789
  check_scl_fct \${scl_fct}
  0.00456:1,0.00789:1

  scl_fct="num:den"
  check_scl_fct \${scl_fct}
  Error: Invalid 'num' in 'num:den' (must be a positive float).

  scl_fct="0.33:den"
  check_scl_fct \${scl_fct}
  Error: Invalid 'den' in '0.33:den' (must be a positive float).

  scl_fct=0.5,1.2a,.65
  check_scl_fct \${scl_fct}
  Error: Invalid numerator-only scaling factor '1.2a' (must be a positive float).

  scl_fct=hotdog
  check_scl_fct \${scl_fct}
  Error: Invalid numerator-only scaling factor 'hotdog' (must be a positive float)
  \`\`\`
EOM
    )

    #  Parse and check function parameter
    if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
        echo "${show_help}"
        return 0
    fi

    #  Split the input string by commas
    IFS=',' read -ra entries <<< "${scl_fct}"

    for entry in "${entries[@]}"; do
        #  Check that entry is a 'num:den' pair (containing a colon)
        if [[ "${entry}" == *:* ]]; then
            #  Split into num and den
            IFS=':' read -r num den <<< "${entry}"

            #  Validate that both 'num' and 'den' are positive floats
            if [[ ! "${num}" =~ ^[0-9]*\.?[0-9]+$ ]]; then
                echo \
                    "Error: Invalid 'num' in '${entry}' (must be a positive" \
                    "float)." >&2
                return 1
            fi

            if [[ ! "${den}" =~ ^[0-9]*\.?[0-9]+$ ]]; then
                echo \
                    "Error: Invalid 'den' in '${entry}' (must be a positive" \
                    "float)." >&2
                return 1
            fi

            #  0-pad values starting with a decimal point
            if [[ "${num}" == .* ]]; then num="0${num}"; fi
            if [[ "${den}" == .* ]]; then den="0${den}"; fi

            #  Add back to the formatted values as 'num:den'
            val_fmt+=( "${num}:${den}" )
        else
            #  Validate single numerator-only entry as a positive float
            if [[ ! "${entry}" =~ ^[0-9]*\.?[0-9]+$ ]]; then
                echo \
                    "Error: Invalid numerator-only scaling factor '${entry}'" \
                    "(must be a positive float)." >&2
                return 1
            fi

            #  Zero-pad values starting with a decimal point
            if [[ "${entry}" == .* ]]; then entry="0${entry}"; fi

            #  Append as 'entry:1' (numerator-only case)
            val_fmt+=( "${entry}:1" )
        fi
    done

    #  Join the validated, formatted values back into a comma-separated string
    out_str=$(IFS=','; echo "${val_fmt[*]}")

    #  Output the final formatted string
    echo "${out_str}"  
}
