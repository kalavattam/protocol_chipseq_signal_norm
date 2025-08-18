#!/bin/bash

#  Convert letters in a user-supplied string to lowercase format
function make_lower() { echo "${1}" | tr '[:upper:]' '[:lower:]'; }
