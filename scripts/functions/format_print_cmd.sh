#!/bin/bash

#  Format and print command string from 'build_cmd' with indentation and line
#+ breaks, handling both plain script calls and SLURM 'sbatch' wrappers
#+ 
#+ Positional parameters:
#+   1, slurm (bol): Command invokes 'sbatch': "true" or "false"
#+   2, scr   (str): Script, including path (only required if 'slurm=true')
function format_print_cmd() {
    local slurm="${1}"
    local scr="${2}"

    if [[ "${slurm}" == "false" ]]; then
        #  For non-SLURM commands, break before each '--' and indent
        sed -E 's| --| \\\n    --|g'
    else
        #  Perform more complicated logic for SLURM commands
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
                    if (token ~ scr) {
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
    fi
}
