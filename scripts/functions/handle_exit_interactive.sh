#!/bin/bash

#  Function to exit script with code 0 if not in "interactive mode"
# shellcheck disable=SC2154
exit_0() { if ! ${interactive}; then exit 0; fi }
