#!/bin/bash

#  Write a warning message to stderr and return code 0
function echo_warning() { echo "Warning: $*" >&2; }
