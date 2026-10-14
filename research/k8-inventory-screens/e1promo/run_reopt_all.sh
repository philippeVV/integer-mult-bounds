#!/bin/bash
# Exact frame re-optimisation (PR1148's frame_model LP, integral here) for each single E1->E2 promotion
# and for all 32 jointly.  usage (from the repo root): run_reopt_all.sh FINAL-WORD.json.gz OUT_DIR
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); W=$1; OUT=$2; mkdir -p "$OUT"
PY="uv run -q --no-project --python 3.12 --with scipy --with numpy python -I $HERE/reopt.py"
P=$(cat "$HERE/evidence/PROMOS-32.txt")
$PY "$W" "$OUT/all.json.gz" --promote $P > "$OUT/all.log" 2>&1
for p in $P; do $PY "$W" "$OUT/single-${p%%:*}.json.gz" --promote "$p" >> "$OUT/singles.log" 2>&1; done
