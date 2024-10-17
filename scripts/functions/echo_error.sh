#!/bin/bash

#  Function to return a message to stderr and return exit code 1
function echo_error() {
    echo "Error: $*" >&2
    return 1
}
