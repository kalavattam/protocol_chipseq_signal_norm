#!/bin/bash

#  Exit script with code 1 if not in “interactive mode”
# shellcheck disable=SC2154
exit_1() { if [[ "${interactive}" != "true" ]]; then exit 1; fi; }
