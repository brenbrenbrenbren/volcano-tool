#!/bin/bash
set -e
cd "$(dirname "$0")"

LOG_FILE="$HOME/Library/Logs/volcano_tool.log"

# Activate venv if present
if [ -f "venv/bin/activate" ]; then
  source "venv/bin/activate"
fi

# Choose Streamlit launcher
if command -v streamlit >/dev/null 2>&1; then
  LAUNCH_CMD=(streamlit)
else
  LAUNCH_CMD=(python3 -m streamlit)
fi

# Launch Streamlit in the background and log output
nohup "${LAUNCH_CMD[@]}" run app.py --server.headless true --browser.serverAddress localhost --server.port 8501 \
  > "$LOG_FILE" 2>&1 &

# Wait briefly for the server to start, then open the browser
if command -v nc >/dev/null 2>&1; then
  for _ in {1..20}; do
    if nc -z localhost 8501; then
      open "http://localhost:8501"
      exit 0
    fi
    sleep 0.5
  done
fi

sleep 1
open "http://localhost:8501"
