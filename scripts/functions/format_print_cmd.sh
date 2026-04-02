#!/bin/bash

function format_print_cmd() {
    local slurm="${1:-}"  # Boolean: whether or not 'sbatch' SLURM command
    local scr="${2:-}"    # Script path for 'sbatch' SLURM command
    local show_help       # Help message

    show_help=$(cat << EOM
Usage:
  format_print_cmd [-h|--hlp|--help] slurm [scr]

Description:
  Format and pretty-print a single command line with indentation and line breaks. Intended to be used with function 'build_cmd' by piping its output into this function.

Positional arguments:
  1  slurm  <bol>  "true" / "t" or "false" / "f":
                     - "true" / "t": treat the input as an 'sbatch' command and format sbatch flags and script arguments separately.
                     - "false" / "f": treat the input as a plain script call.
  2  scr    <str>  Script path as it appears in the 'sbatch' command. Required when 'slurm=true', otherwise ignored.

Examples:
  '''bash
  #  Plain command
  build_cmd ... | format_print_cmd false

  #  SLURM sbatch command invoking 'execute_compute_signal.sh'
  build_cmd ...
      | format_print_cmd true "\${dir_rep}/scripts/execute_compute_signal.sh"
  '''

Notes:
- For non-SLURM commands, the function does the following:
    + Reads a single-line command string from stdin.
    + Breaks before each '--' and indents the resulting lines.
- For SLURM-wrapped commands ('sbatch ... <script> ...'), the function:
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
        echo "Error: Positional argument 1, 'slurm', is missing." >&2
        echo >&2
        echo "${show_help}" >&2
        return 1
    fi

    #  Lower-case flag 'slurm' in a Bash-3.2-safe way
    local slurm_lc
    slurm_lc=$(printf '%s' "${slurm}" | tr '[:upper:]' '[:lower:]')

    if [[ "${slurm_lc}" =~ ^(false|f)$ ]]; then
        #  For non-SLURM commands, break before each '--' and indent
        sed -E 's| --| \\\n    --|g'  # Only as robust as incoming command str
    elif [[ "${slurm_lc}" =~ ^(true|t)$ ]]; then
        #  For SLURM commands, require a script path
        if [[ -z "${scr}" ]]; then
            echo \
                "Error in 'format_print_cmd': script path 'scr' (positional" \
                "argument 2) must be supplied when 'slurm=true' or" \
                "'slurm=t'." >&2
            return 1
        fi

        #  ...and perform more complicated logic; formatter expects current
        #+ string layout emitted by 'build_cmd'
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
        echo \
            "Error in 'format_print_cmd': unknown 'slurm' value '${slurm}';" \
            "expected 'true' / 't' or 'false' / 'f'." >&2
        return 1
    fi
}
