#!/bin/bash

#  Function to write a warning message to stderr and return code 0
function echo_warning() {
    echo "Warning: $*" >&2
}
