#!/usr/bin/env bash
# macOS: double-click this file in Finder (the first time: right-click, Open). It runs start_calculator.sh.
cd "$(dirname "$0")"
./start_calculator.sh "$@"
status=$?
if [ $status -ne 0 ]; then echo; echo "The calculator stopped with an error (code $status). Press Return to close."; read -r _; fi
exit $status
