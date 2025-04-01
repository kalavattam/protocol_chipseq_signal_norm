#!/bin/bash

#  install_envs.sh
#  KA


#  Run script in interactive mode (true) or command-line mode (false)
interactive=false

#  Exit on errors, unset variables, or pipe failures if not in "interactive
#+ mode"
if ! ${interactive:-false}; then set -euo pipefail; fi

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
for script in \
    check_installed_env.sh \
    check_installed_mamba.sh \
    echo_error.sh \
    echo_warning.sh \
    handle_env.sh
do
    source "${dir_fnc}/${script}"
done


#  Set values for interactive mode
function set_interactive() {
    #  Hardcoded argument assignments
    env_nam="env_align"
    yes=false
}


#  Parse arguments ============================================================
#  Initialize variables along with default assignments
env_nam=""
yes=false

show_help=$(cat << EOM
install_envs.sh
  --env_nam <str> [--yes]

Description:
  Sets up one of multiple Mamba environments with necessary programs and
  dependencies used in this project. For more details, see "Notes" below.

Arguments:
  -h, --help     Display this help message and exit.
 -en, --env_nam  Mamba environment to create: 'env_align', 'env_analyze',
                 'env_protocol', or 'env_siqchip'.
  -y, --yes      Automatically answer yes to mamba prompts (optional).

Dependencies:
  - Programs
    + Bash or Zsh
    + Conda or Mamba
  - Functions
    + check_installed_env
    + check_installed_mamba
    + echo_error
    + echo_warning
    + handle_env_deactivate

Notes:
  - The call to 'mamba' will be adapted to allow Rosetta translation if the
    script detects that the system has "arm64" architecture.
  - The following packages are installed via a call to 'mamba' create:
    + env_align
      - bamtools
      - bbmap
      - bedtools
      - bowtie2
      - bwa
      - fastqc
      - macs3
      - minimap
      - mosdepth
      - parallel
      - picard
      - preseq
      - rename
      - samtools
      - subread
      - ucsc-bedgraphtobigwig
      - ucsc-bedsort
      - ucsc-facount
      - wget
    + env_analyze
      - bioconductor-annotationdbi
      - bioconductor-chipqc
      - bioconductor-chipseeker
      - bioconductor-clusterprofiler
      - bioconductor-deseq2
      - bioconductor-diffbind
      - bioconductor-edger
      - bioconductor-enhancedvolcano
      - bioconductor-genomicfeatures
      - bioconductor-genomicranges
      - bioconductor-ihw
      - bioconductor-iranges
      - bioconductor-pcatools
      - bioconductor-sva
      - deeptools
      - ipython
      - pandas
      - parallel
      - pbzip2
      - phantompeakqualtools
      - pigz
      - r-argparse
      - r-dendextend
      - r-devtools
      - r-ggalt
      - r-ggpubr
      - r-ggrepel
      - r-ggsci
      - r-pheatmap
      - r-plotly
      - r-readxl
      - r-rjson
      - r-tidyverse
      - r-upsetr
      - r-venneuler
      - r-writexl
      - r-xml2
      - rename
      - tree
    + env_protocol
      - bc
      - bowtie2
      - fastqc
      - gawk
      - ipython
      - parallel
      - pbzip2
      - pigz
      - pysam
      - python=3.11  # Restrict to version 3.11 for 'sequali' installation
      - r-argparse
      - r-plotly
      - r-ggsci
      - rename
      - samtools
      - sequali
      - wget
    + env_siqchip
      - bc
      - bedtools
      - gfortran
      - gnuplot
      - parallel
      - samtools
      - ucsc-bedgraphtobigwig
      - ucsc-bedclip
  - Given that >100 packages and dependencies may be installed depending on
    the environment specified, the execution and completion of the 'mamba
    create' operation will take >10 minutes if making use of cached packages
    from previous installations/separate environments. If this is a fresh
    installation that does not make use of cached packages, then the creation 
    of the environment may take even longer: e.g., more than 20 or 30 minutes.

Example:
  \`\`\`
  bash install_envs.sh
      --env_nam "env_protocol"
      --yes
  \`\`\`
EOM
)

#  Parse arguments
if [[ -z "${1:-}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    echo "${show_help}"
    if ! ${interactive}; then exit 0; fi
fi

if ${interactive:-false}; then
    set_interactive
else
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
           -en|--env_nam) env_nam="${2}"; shift 2 ;;
            -y|--yes)      yes=true;        shift 1 ;;
            *)
                echo "## Unknown parameter passed: '${1}' ##" >&2
                echo "" >&2
                echo "${show_help}" >&2
                exit 1
                ;;
        esac
    done
fi


#  Check that required arguments are provided, appropriate, formatted, etc.
if [[ -z "${env_nam}" ]]; then
    echo_error "Environment name was not specified."
fi

case "${env_nam}" in
    env_align|env_analyze|env_protocol|env_siqchip) : ;;
    *) 
        echo_error \
            "Invalid environment name specified. Must be 'env_align', " \
            "'env_analyze', 'env_protocol', or 'env_siqchip'."
        ;;
esac

if check_installed_env "${env_nam}"; then
    echo_error \
        "An environment with the name '${env_nam}' is already installed."
fi

#  Check that dependencies are in PATH
check_installed_mamba


#  Do the main work ===========================================================
echo "Creating environment '${env_nam}'."

case "${env_nam}" in
    env_align|env_analyze)
        echo_warning \
            "Creating '${env_nam}' will take some time given the >100" \
            "packages to install. Don't worry if this script is still" \
            "running without any apparent progress after 10, 20, or even 30" \
            "minutes: It will eventually complete."
    ;;
esac

#  If not in base environment, then deactivate current environment
handle_env_deactivate

#  Construct the mamba command
cmd_mamba="mamba create -n ${env_nam}"

if ${yes}; then cmd_mamba+=" --yes"; fi

#  Assign an array of packages to install
if [[ "${env_nam}" == "env_align" ]]; then
    packages=(
        bamtools
        bbmap
        bedtools
        bowtie2
        bwa
        fastqc
        macs3
        minimap
        mosdepth
        parallel
        picard
        preseq
        rename
        samtools
        subread
        ucsc-bedgraphtobigwig
        ucsc-bedsort
        ucsc-facount
        wget
    )
elif [[ "${env_nam}" == "env_analyze" ]]; then
    packages=(
        bioconductor-annotationdbi
        bioconductor-chipqc
        bioconductor-chipseeker
        bioconductor-clusterprofiler
        bioconductor-deseq2
        bioconductor-diffbind
        bioconductor-edger
        bioconductor-enhancedvolcano
        bioconductor-genomicfeatures
        bioconductor-genomicranges
        bioconductor-ihw
        bioconductor-iranges
        bioconductor-pcatools
        bioconductor-sva
        ipython
        deeptools
        pandas
        parallel
        pbzip2
        phantompeakqualtools
        pigz
        r-argparse
        r-dendextend
        r-devtools
        r-ggalt
        r-ggpubr
        r-ggrepel
        r-ggsci
        r-pheatmap
        r-plotly
        r-readxl
        r-rjson
        r-tidyverse
        r-upsetr
        r-venneuler
        r-writexl
        r-xml2
        rename
        tree
    )
elif [[ "${env_nam}" == "env_protocol" ]]; then
    packages=(
        bc
        bowtie2
        fastqc
        gawk
        ipython
        multiqc
        parallel
        pbzip2
        pigz
        pysam
        python=3.11  # Restrict to version 3.11 for 'sequali' installation
        r-argparse
        r-plotly
        r-ggsci
        rename
        samtools
        sequali
        wget
    )
elif [[ "${env_nam}" == "env_siqchip" ]]; then
    packages=(
        bc
        bedtools
        gfortran
        gnuplot
        parallel
        samtools
        ucsc-bedgraphtobigwig
        ucsc-bedclip
    )
fi

#  Run the mamba environment installation
if ! eval "${cmd_mamba} ${packages[*]}"; then
    echo_error \
        "Failed to create environment '${env_nam}'. Please check the error" \
        "message(s) above."
    exit 1
fi
