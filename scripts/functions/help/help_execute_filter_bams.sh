#!/bin/bash
# -*- coding: utf-8 -*-
#
# Compatibility alias for the renamed filter_alignments help function.

# shellcheck disable=SC1091
source "$(
    cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd
)/help_execute_filter_alignments.sh"

function help_execute_filter_bams() {
    help_execute_filter_alignments "$@"
}
