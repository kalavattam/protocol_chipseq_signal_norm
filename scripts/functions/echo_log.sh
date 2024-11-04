#!/bin/bash

echo_log() {
    # shellcheck disable=SC2154
    if ${dry_run} || ${verbose}; then
        echo "$@"
    fi
}
