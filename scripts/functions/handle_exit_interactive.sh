#!/bin/bash

#  Exit script with code 0 if not in “interactive mode”
# shellcheck disable=SC2154
function exit_0() { if [[ "${interactive}" != "true" ]]; then exit 0; fi; }

#  Exit script with code 1 if not in “interactive mode”
function exit_1() { if [[ "${interactive}" != "true" ]]; then exit 1; fi; }
