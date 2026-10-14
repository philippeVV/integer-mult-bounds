#!/bin/sh
# Run from the repository root. Writes JSON outputs to ${OUT:-/tmp/bit-helper-stop-obstruction}.
set -e
D=research/bit-helper-stop-obstruction; OUT=${OUT:-/tmp/bit-helper-stop-obstruction}; mkdir -p "$OUT"
python3 -I -B $D/census.py . "$OUT/paths.json"
python3 -I -B $D/stopenum.py . "$OUT/enum.json"
python3 -I -B $D/relocate.py . "$OUT/relocate.json"
python3 -I -B $D/relocate.py . "$OUT/relocate_dirty.json" --dirty
python3 -I -B $D/witness.py . "$OUT/witness_rows.json"
python3 -I -B $D/replay_merges.py . | tail -3
