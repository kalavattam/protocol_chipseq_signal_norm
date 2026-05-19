#!/bin/bash
# -*- coding: utf-8 -*-
#
# Compatibility wrapper for the renamed filter_alignments chain.

set -euo pipefail

dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"

exec "${BASH}" "${dir_scr}/submit_filter_alignments.sh" "$@"
