#!/usr/bin/env bash
# One command on macOS and Linux:  ./start_calculator.sh          (the first start sets everything up)
#                                  ./start_calculator.sh check    (set up, self-test, stop)
# Needs Python 3.9 or newer and, once, an internet connection to install four Python libraries.
set -e
cd "$(dirname "$0")"
VENV=".venv"

make_venv() {
  if "$1" -m venv "$VENV" 2>/dev/null; then return 0; fi
  # Debian / Ubuntu without the package python3-venv: make the environment without pip and add pip from the web
  rm -rf "$VENV"
  if "$1" -m venv --without-pip "$VENV" 2>/dev/null; then
    echo "The Python venv module has no pip here; fetching pip (needs internet) ..."
    if command -v curl >/dev/null 2>&1; then curl -fsSL https://bootstrap.pypa.io/get-pip.py -o "$VENV/get-pip.py"
    elif command -v wget >/dev/null 2>&1; then wget -q https://bootstrap.pypa.io/get-pip.py -O "$VENV/get-pip.py"
    else echo "Neither curl nor wget found."; fi
    if [ -f "$VENV/get-pip.py" ] && "$VENV/bin/python" "$VENV/get-pip.py" --quiet; then rm -f "$VENV/get-pip.py"; return 0; fi
  fi
  rm -rf "$VENV"
  echo "Could not create the Python environment."
  echo "On Debian / Ubuntu install it with:  sudo apt install python3-venv python3-pip"
  exit 1
}

if [ ! -x "$VENV/bin/python" ]; then
  PY="${BI3_PYTHON:-}"
  if [ -z "$PY" ]; then
    for c in python3.12 python3.13 python3.14 python3.11 python3.10 python3.9 python3; do
      if command -v "$c" >/dev/null 2>&1 && "$c" -c 'import sys; sys.exit(0 if sys.version_info[:2]>=(3,9) else 1)'; then PY="$c"; break; fi
    done
  fi
  if [ -z "$PY" ]; then echo "Python 3.9 or newer was not found. Install it from https://www.python.org/downloads/ (macOS) or with your package manager (Linux)."; exit 1; fi
  echo "First start: setting up the calculator (about 2 to 5 minutes, one time only) using $PY"
  make_venv "$PY"
fi
if ! cmp -s requirements.txt "$VENV/requirements.installed.txt"; then
  "$VENV/bin/python" -m pip install --disable-pip-version-check -r requirements.txt
  cp requirements.txt "$VENV/requirements.installed.txt"
  "$VENV/bin/python" bi3_launcher.py --check
fi
if [ "${1:-}" = "check" ]; then "$VENV/bin/python" bi3_launcher.py --check; exit $?; fi
exec "$VENV/bin/python" bi3_launcher.py
