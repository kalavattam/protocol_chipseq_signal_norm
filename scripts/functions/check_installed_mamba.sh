#!/bin/bash

#  Function to check that either Mamba (priority) or Conda (fallback) is
#+ installed
function check_installed_mamba() {
    #  Check for Mamba
    if command -v mamba &> /dev/null; then
        :
    #  If Mamba is not found, check for Conda
    elif command -v conda &> /dev/null; then
        :
    #  If neither Mamba nor Conda is found, provide an error message
    else
        echo "Error: Neither Mamba nor Conda is installed on your system." >&2
        echo "" >&2
        echo \
            "Mamba is a package manager that makes package installations" \
            "faster and more reliable, e.g., in comparison to Conda." >&2
        echo "" >&2
        echo \
            "For Mamba installation instructions, please check the following" \
            "link: https://github.com/mamba-org/mamba#installation" >&2
        echo "" >&2
        return 1
    fi
}
