#!/usr/bin/env bash
# Launch the interactive supplement (Mo-W Thermodynamic model + Two-stage fractionation model).
# Creates a local virtual environment on first run, installs the pinned requirements,
# then starts Streamlit and opens the browser. Re-run the script to start the app again.
set -euo pipefail

cd "$(dirname "$0")"

PY="${PYTHON:-python3}"
if ! command -v "$PY" >/dev/null 2>&1; then
    echo "python3 not found. Install Python 3.9 or newer and re-run." >&2
    exit 1
fi

if [ ! -x .venv/bin/python ]; then
    echo "Creating virtual environment (.venv) ..."
    "$PY" -m venv .venv
fi

echo "Installing / checking requirements ..."
.venv/bin/python -m pip install --quiet --upgrade pip
.venv/bin/python -m pip install --quiet -r requirements.txt

echo "Starting the dashboard at http://localhost:8501  (Ctrl+C to stop)"
exec .venv/bin/python -m streamlit run dashboard.py --server.port "${PORT:-8501}" --server.address localhost
