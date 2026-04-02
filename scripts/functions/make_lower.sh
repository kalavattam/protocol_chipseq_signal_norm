#!/bin/bash
# -*- coding: utf-8 -*-
# 
# Script: make_lower.sh
#
# Copyright 2025-2026 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# OpenAI ChatGPT (GPT-5-series models) was used in development.
# 
# Distributed under the MIT license.

#  Convert a user-supplied string to lowercase
function make_lower() { printf '%s\n' "${1:-}" | tr '[:upper:]' '[:lower:]'; }
