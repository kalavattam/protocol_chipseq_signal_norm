#!/bin/bash

#  write_header.sh
#  KA

#  Run script in interactive mode (true) or command-line mode (false)
interactive=false

#  Exit on errors, unset variables, or pipe failures if not in "interactive
#+ mode"
if ! ${interactive}; then set -euo pipefail; fi

#  Set the path to the "scripts" directory
if ${interactive:-false}; then
    ## WARNING: If 'interactive=true', change path as needed ##
    dir_scr="${HOME}/repos/protocol_chipseq_signal_norm/scripts"
else
    dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi


#  Source and define functions ================================================
#  Set the path to the "functions" directory
dir_fnc="${dir_scr}/functions"

# shellcheck disable=SC1090
for fnc in check_exists_file_dir echo_error echo_warning exit_0 exit_1; do
    source "${dir_fnc}/${fnc}.sh"
done
unset fnc


#  Set paths, values, arguments, etc. for interactive mode
function set_interactive() {
    #  Set base repository paths
    dir_bas="${HOME}/repos"  ## WARNING: Change as needed ##
    dir_rep="${dir_bas}/protocol_chipseq_signal_norm"
    dir_scr="${dir_rep}/scripts"

    #  Set alignment parameters
    aligner="bowtie2"
    a_type="global"
    req_flg=true
    flg="$(if ${req_flg}; then echo "2"; else echo "NA"; fi)"
    mapq=1
    str_det="${aligner}_${a_type}_flag-${flg}_mapq-${mapq}"

    #  Set output directories
    dir_dat="${dir_rep}/data"
    dir_pro="${dir_dat}/processed"
    dir_out="${dir_pro}/compute_signal/${str_det}/tables"

    #  Set hardcoded argument assignments
    verbose=true
    dry_run=true
    mode="spike"
    fil_out="${dir_out}/ChIP_WT_G1-G2M-Q_Hho1-Hmo1_${mode}_6nd.tsv"
}


#  Initialize argument variables, check and parse arguments, etc. =============
#  Initialize argument variables, assigning default values where applicable
verbose=false
dry_run=false
mode="alpha"
fil_out=""

#  Assign variable for help message
show_help=$(cat << EOM
Usage:
  write_header_alpha.sh [--verbose] [--dry-run] --mode <str> --fil_out <str>

Description:
  Writes a predefined tab-delimited header to the specified output file.
  Implemented to pre-write the header before running SLURM jobs, preventing
  race conditions that occurred when header-writing logic was in 'submit'
  scripts.

Options:
   -h, --help     Display this help message and exit.
   -v, --verbose  Print the header before writing.
  -dr, --dry-run  Print the header but do not write to a file.
   -m, --mode     Type of header to write: 'alpha' or 'spike' (default:
                  '${mode}').
  -fo, --fil_out  Output file where the header should be written.

Dependencies:
  - Bash
  - echo
  - printf

Example:
  \`\`\`
  write_header.sh -v -t alpha -o results/ChIP_samples_alpha_6nd.tsv
  \`\`\`
EOM
)

#  Parse arguments
if [[ -z "${1:-}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    echo "${show_help}"
    exit_0
fi

if ${interactive:-false}; then
    set_interactive
else
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
             -v|--verbose) verbose=true;   shift 1 ;;
            -dr|--dry-run) dry_run=true;   shift 1 ;;
             -m|--mode)    mode="${2}";    shift 2 ;;
            -fo|--fil_out) fil_out="${2}"; shift 2 ;;
            *)
                echo "## Unknown parameter passed: '${1}' ##" >&2
                echo "${show_help}" >&2
                exit_1
                ;;
        esac
    done
fi

#  Check arguments
case "${mode}" in
    alpha|spike) : ;;
    *)
        echo_error \
            "Header mode ('--mode') was assigned '${mode}' but must be" \
            "'alpha' or 'spike'."
        exit_1
        ;;
esac

check_exists_file_dir "f" "${fil_out}" "fil_out"


#  Do the main work ===========================================================
#  Define the header column names as an array
case "${type}" in
    alpha)
        nam_col=(
            "fil_ip" "fil_in" "alpha" "eqn"
            "mass_ip" "mass_in" "vol_all" "vol_in" "dep_ip" "dep_in"
            "len_ip" "len_in"
            "dm_fr_1" "dm_fr_5" "dm_fr_10" "dm_fr_20" "dm_fr_30" "dm_fr_40"
            "dm_fr_50"
            "dm_nm_1" "dm_nm_5" "dm_nm_10" "dm_nm_20" "dm_nm_30" "dm_nm_40"
            "dm_nm_50"
        )
        ;;
        spike) nam_col=(
            "main_ip" "spike_ip" "main_in" "spike_in" "sf"
            "num_mp" "num_sp" "num_mn" "num_sn"
            "dm_fr_1" "dm_fr_5" "dm_fr_10" "dm_fr_20" "dm_fr_30" "dm_fr_40"
            "dm_fr_50"
            "dm_nm_1" "dm_nm_5" "dm_nm_10" "dm_nm_20" "dm_nm_30" "dm_nm_40"
            "dm_nm_50"
        )
        ;;
esac

#  Generate 'printf' format string dynamically
fmt_str=$(printf "%s\t" "${nam_col[@]}")
fmt_str="${fmt_str%$'\t'}\n"  # Remove trailing tab and add newline

#  Print the formatted header line
# shellcheck disable=SC2059
header=$(printf "${fmt_str}" "${nam_col[@]}")

#  Print the header (if in dry-run or verbose modes)
if ${dry_run} || ${verbose}; then
    echo "##################"
    echo "## Table header ##"
    echo "##################"
    echo ""
    echo "${header}"
    echo ""
    echo ""
fi

#  Write the header to the file (if not in dry-run mode)
if ! ${dry_run}; then echo "${header}" > "${fil_out}"; fi
