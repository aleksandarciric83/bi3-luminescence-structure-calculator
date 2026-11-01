#!/usr/bin/env bash
# macOS / Linux: build the self-contained program dist/Bi3_Calculator (all libraries inside).
# Needs internet once on this computer. Run it on the system the program is meant for.
set -e
cd "$(dirname "$0")"
[ -x .venv/bin/python ] || ./start_calculator.sh check
.venv/bin/python build_app.py
