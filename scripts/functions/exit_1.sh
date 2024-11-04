#!/bin/bash

#  Function to exit script with code 1 if not in "interactive mode"
# shellcheck disable=SC2154
exit_1() { if ! ${interactive}; then exit 1; fi }
